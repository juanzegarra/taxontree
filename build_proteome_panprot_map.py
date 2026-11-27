# Author: Juan Zegarra
# Date 19/11/26
#
# Builds a SQLite index of UniProt complete proteomes and pan-proteomes
# (byte offsets into the FASTA files), reads which reference proteomes make up
# each pan-proteome from the FTP *_stats.json files, and produces a combined
# table: complete-proteome sequences that have NO pan-proteome + every
# pan-proteome sequence.
#
# Usage:
#   python build_proteome_panprot_map.py PANPROT.fasta COMPLETE.fasta [--db pp_prot.db] [--overwrite]

###################################################################

import argparse
import gzip
import io
import json
import multiprocessing as mp
import re
from ftplib import FTP, error_perm
from itertools import islice
from pathlib import Path

from sqlalchemy import BigInteger, Integer, String, create_engine, func, insert, literal, select, text, union_all
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column

PP_FTP = "ftp.uniprot.org"
PP_DIR = "pub/databases/uniprot/knowledgebase/pan_proteomes/"
CORES = 2
DB_NAME = "pp_prot.db"
FASTA_PP = "panprot_db"
FASTA_COMPLETE = "complete_prot_db"
FTP_PP = "ftp_pp_db"
COMBINED = "combined_pp_complete_prot_db"
CHUNK = 50_000  # rows per bulk insert

# complete:    >tr|A0A0A7KU74|A0A0A7KU74_9VIRU VP1 OS=... OX=1581151 PE=3 SV=1 UP=UP000097907
# panproteome: >pp100177|UPI003AB4B090 67% UP001301963 XPK47238 OX=100177 OS=... ; tr|A0ACD5GM08|...
OX_RE = re.compile(r"\bOX=(\d+)")
UP_COMPLETE_RE = re.compile(r"\bUP=(UP\d+)")
PP_ID_RE = re.compile(r"^(pp\d+)\|")
UP_PP_RE = re.compile(r"\b(UP\d+)\b")


class Base(DeclarativeBase):
    pass


class CompleteProteomes(Base):
    __tablename__ = FASTA_COMPLETE

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    proteome_id: Mapped[str] = mapped_column(String(50))
    byte_offset: Mapped[int] = mapped_column(BigInteger)  # start of the record in the FASTA
    length: Mapped[int] = mapped_column(BigInteger)       # byte length of the whole record
    taxid: Mapped[int] = mapped_column(Integer)


class PanProteomes(Base):
    __tablename__ = FASTA_PP

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    pp_id: Mapped[str] = mapped_column(String(50))
    reference_proteome: Mapped[str] = mapped_column(String(50))
    byte_offset: Mapped[int] = mapped_column(BigInteger)
    length: Mapped[int] = mapped_column(BigInteger)
    taxid: Mapped[int] = mapped_column(Integer)


class PPftp(Base):
    """One row per (pan-proteome, reference proteome) pair from the FTP JSONs."""
    __tablename__ = FTP_PP

    pp_id: Mapped[str] = mapped_column(String(50), primary_key=True)
    complete_prots: Mapped[str] = mapped_column(String(100), primary_key=True)


class CombinedProteins(Base):
    __tablename__ = COMBINED

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    source: Mapped[str] = mapped_column(String(10), primary_key=True)  # "complete" | "pan": which FASTA the offset refers to
    byte_offset: Mapped[int] = mapped_column(BigInteger)
    length: Mapped[int] = mapped_column(BigInteger)


# Created after the bulk load (much faster than maintaining them during inserts)
INDEXES = [
    f"CREATE INDEX IF NOT EXISTS ix_{FASTA_COMPLETE}_proteome_id ON {FASTA_COMPLETE}(proteome_id)",
    f"CREATE INDEX IF NOT EXISTS ix_{FASTA_COMPLETE}_taxid ON {FASTA_COMPLETE}(taxid)",
    f"CREATE INDEX IF NOT EXISTS ix_{FASTA_PP}_pp_id ON {FASTA_PP}(pp_id)",
    f"CREATE INDEX IF NOT EXISTS ix_{FASTA_PP}_taxid ON {FASTA_PP}(taxid)",
    f"CREATE INDEX IF NOT EXISTS ix_{FTP_PP}_complete_prots ON {FTP_PP}(complete_prots)",
]


def create_database(base, db_name=DB_NAME):
    engine = create_engine(f"sqlite:///{db_name}")
    base.metadata.create_all(engine)
    return engine


###################################################################
# FTP: which reference proteomes belong to each pan-proteome
###################################################################

_ftp = None  # one connection per worker process


def init_worker():
    global _ftp
    _ftp = FTP(PP_FTP, timeout=60)
    _ftp.login()


def _retrieve(path):
    for attempt in range(3):
        try:
            buf = io.BytesIO()
            _ftp.retrbinary(f"RETR {path}", buf.write)
            return buf.getvalue()
        except (EOFError, OSError):  # connection dropped: reconnect and retry
            if attempt == 2:
                raise
            init_worker()

def list_pp_dirs():
    with FTP(PP_FTP) as ftp:
        ftp.login()
        ftp.cwd(PP_DIR)
        try:
            names = ftp.nlst()
        except error_perm:  # server without NLST: parse a plain LIST instead
            lines = []
            ftp.retrlines("LIST", lines.append)
            names = [l.split()[-1] for l in lines if l.startswith("d")]

    # NLST may return bare names or full paths; keep only ppNNN entries
    names = (n.rstrip("/").rsplit("/", 1)[-1] for n in names)
    pp_dirs = sorted(n for n in names if re.fullmatch(r"pp\d+", n))
    if not pp_dirs:
        raise RuntimeError(f"No ppNNN directories found under {PP_DIR}; check the listing manually")
    return pp_dirs

def fetch_pp(pp_dir):
    """Download <pp_dir>/<pp_dir>_stats.json -> (pp_id, [reference proteome ids]).

    Returns (pp_id, None) if the file could not be found.
    """
    raw = None
    for suffix in ("_stats.json", "_stats.json.gz"):
        try:
            raw = _retrieve(f"{PP_DIR}{pp_dir}/{pp_dir}{suffix}")
            break
        except error_perm:  # file not there under this name
            continue
    if raw is None:
        return pp_dir, None

    if raw[:2] == b"\x1f\x8b":  # gzip magic number
        raw = gzip.decompress(raw)
    data = json.loads(raw.decode("utf-8"))

    # NB: rp_composition lists only the *reference* proteomes of the pan-proteome
    # (e.g. 13 of the 129 proteomes clustered for pp9).
    composition = (data.get("pan_proteome") or {}).get("rp_composition") or []
    return pp_dir, [rp["rp_id"] for rp in composition if rp.get("rp_id")]


def traverse_ftp_dir(cores=CORES):
    pp_dirs = list_pp_dirs()
    print(f"{len(pp_dirs)} pan-proteome directories found")

    with mp.Pool(processes=cores, initializer=init_worker) as pool:
        results = pool.map(fetch_pp, pp_dirs, chunksize=50)

    missing = [pp for pp, rps in results if rps is None]
    if missing:
        print(f"WARNING: no stats JSON found for {len(missing)} pan-proteomes (e.g. {missing[:5]})")
    return {pp: rps for pp, rps in results if rps}


###################################################################
# FASTA indexing
###################################################################

def scan_fasta(path):
    """Yield (header, byte_offset, byte_length) for every record of a FASTA file.

    The record covers the header line plus all its sequence lines, so a sequence
    can later be read with f.seek(byte_offset); f.read(byte_length).
    """
    with open(path, "rb") as f:
        offset = start = 0
        header = None
        for line in f:
            if line.startswith(b">"):
                if header is not None:
                    yield header, start, offset - start
                header = line[1:].decode("utf-8", errors="replace").rstrip()
                start = offset
            offset += len(line)
        if header is not None:
            yield header, start, offset - start


def parse_complete(header, offset, length):
    ox, up = OX_RE.search(header), UP_COMPLETE_RE.search(header)
    if not (ox and up):
        return None
    return dict(id=header.split(None, 1)[0], proteome_id=up.group(1),
                byte_offset=offset, length=length, taxid=int(ox.group(1)))


def parse_pan(header, offset, length):
    pp, up, ox = PP_ID_RE.match(header), UP_PP_RE.search(header), OX_RE.search(header)
    if not (pp and up and ox):
        return None
    return dict(id=header.split(None, 1)[0], pp_id=pp.group(1), reference_proteome=up.group(1),
                byte_offset=offset, length=length, taxid=int(ox.group(1)))


def bulk_insert(session, model, rows):
    """Insert an iterable of dicts in chunks (duplicate keys are skipped)."""
    stmt = insert(model.__table__).prefix_with("OR IGNORE")
    it = iter(rows)
    while chunk := list(islice(it, CHUNK)):
        session.execute(stmt, chunk)


def count_rows(session, model):
    return session.scalar(select(func.count()).select_from(model.__table__))


def build_fasta_db(model, fasta_file, session, parser):
    stats = {"read": 0, "skipped": 0}

    def rows():
        for header, offset, length in scan_fasta(fasta_file):
            stats["read"] += 1
            row = parser(header, offset, length)
            if row is None:
                stats["skipped"] += 1
            else:
                yield row

    bulk_insert(session, model, rows())
    session.commit()

    stored = count_rows(session, model)
    duplicates = stats["read"] - stats["skipped"] - stored
    print(f"{model.__tablename__}: {stats['read']} records read, {stats['skipped']} skipped "
          f"(header did not match), {duplicates} duplicate ids, {stored} stored")


def populate_pp_json_db(json_dict, session):
    rows = ({"pp_id": pp, "complete_prots": rp} for pp, rps in json_dict.items() for rp in set(rps))
    bulk_insert(session, PPftp, rows)
    session.commit()
    print(f"{FTP_PP}: {count_rows(session, PPftp)} (pan-proteome, reference proteome) pairs stored")


###################################################################
# Combine: complete proteomes WITHOUT a pan-proteome + all pan-proteome sequences
###################################################################

def build_combined(session):
    no_pp = select(
        CompleteProteomes.id,
        literal("complete").label("source"),
        CompleteProteomes.byte_offset,
        CompleteProteomes.length,
    ).where(CompleteProteomes.proteome_id.not_in(select(PPftp.complete_prots)))

    pan = select(
        PanProteomes.id,
        literal("pan").label("source"),
        PanProteomes.byte_offset,
        PanProteomes.length,
    )

    session.execute(
        insert(CombinedProteins).from_select(
            ["id", "source", "byte_offset", "length"], union_all(no_pp, pan)
        )
    )
    session.commit()
    print(f"{COMBINED}: {count_rows(session, CombinedProteins)} sequences stored")


def fetch_record(combined_row, fasta_pp_path, fasta_complete_path):
    """Read one FASTA record back, given a row (source, byte_offset, length) of the combined table."""
    path = fasta_pp_path if combined_row.source == "pan" else fasta_complete_path
    with open(path, "rb") as f:
        f.seek(combined_row.byte_offset)
        return f.read(combined_row.length).decode("utf-8")

def export_combined_fasta(engine, fasta_pp_path, fasta_complete_path, out_path):
    q = select(
        CombinedProteins.source, CombinedProteins.byte_offset, CombinedProteins.length
    ).order_by(CombinedProteins.source, CombinedProteins.byte_offset)

    with open(fasta_pp_path, "rb") as pan, open(fasta_complete_path, "rb") as comp, \
         open(out_path, "wb") as out, Session(engine) as session:
        files = {"pan": pan, "complete": comp}
        n = 0
        for source, offset, length in session.execute(q).yield_per(10_000):
            f = files[source]
            f.seek(offset)
            data = f.read(length)
            if not data.endswith(b"\n"):   # last record of a file may lack a newline
                data += b"\n"
            out.write(data)
            n += 1
    print(f"{n} sequences written to {out_path}")


###################################################################

def main():
    parser = argparse.ArgumentParser(description=__doc__ or "Build UniProt proteome / pan-proteome SQLite index")
    parser.add_argument("fasta_pp", help="pan-proteomes FASTA file")
    parser.add_argument("fasta_complete", help="complete proteomes FASTA file")
    parser.add_argument("--db", default=DB_NAME, help=f"SQLite file (default: {DB_NAME})")
    parser.add_argument("--cores", type=int, default=CORES, help="parallel FTP connections")
    parser.add_argument("--overwrite", action="store_true", help="replace the database if it already exists")
    parser.add_argument("--export-fasta", help="write the combined sequences to this FASTA file")
    args = parser.parse_args()

    if not Path(args.fasta_pp).is_file() or not Path(args.fasta_complete).is_file():
        print("Could not find Panproteomes or complete proteomes file. exiting...")
        return

    db_path = Path(args.db)
    if db_path.exists():
        if not args.overwrite:
            print(f"{db_path} already exists; use --overwrite to rebuild it. exiting...")
            return
        db_path.unlink()

    pp_composition = traverse_ftp_dir(args.cores)
    engine = create_database(Base, args.db)

    with Session(engine) as session:
        populate_pp_json_db(pp_composition, session)
        build_fasta_db(CompleteProteomes, args.fasta_complete, session, parse_complete)
        build_fasta_db(PanProteomes, args.fasta_pp, session, parse_pan)

    with engine.begin() as conn:
        for ddl in INDEXES:
            conn.execute(text(ddl))

    with Session(engine) as session:
        build_combined(session)

    if args.export_fasta:
        export_combined_fasta(engine, args.fasta_pp, args.fasta_complete, args.export_fasta)


if __name__ == "__main__":
    main()
