# Author: Juan Zegarra
# Date 19/11/26

###################################################################

from ftplib import FTP
import os
import re
import json
import gzip
import io
import multiprocessing as mp
import pandas as pd
from sqlalchemy import String, BigInteger, Integer, create_engine, select, unionall, Table, MetaData, Column
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column
from sqlalchemy import JSON
from Bio import  SeqIO

PP_FTP = "ftp.uniprot.org"
CORES = 2
DB_NAME = "pp_prot.db"
FASTA_PP = "panprot_db"
FASTA_COMPLETE  = "complete_prot_db"
FTP_PP = "ftp_pp_db"
proteome_id_regex = re.compile(r"\bOX=([0-9+])\b.*UP=\b(UP[0-9]+)\b")
panproteome_id_regex = re.compile(r"(pp[0-9]+).*\b(UP[0-9]+)\b.*\bOX=([0-9]+)\b")



class Base(DeclarativeBase):
    pass

class CompleteProteomes(Base):
    __tablename__ = FASTA_COMPLETE

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    proteome_id: Mapped[str] = mapped_column(String(50), index=True)
    byte_offset: Mapped[int] = mapped_column(BigInteger)  # start byte offset
    length: Mapped[int] = mapped_column(BigInteger)      # byte length of record
    taxid: Mapped[int] = mapped_column(Integer, index=True)


class PanProteomes(Base):
    __tablename__ = FASTA_PP

    id: Mapped[str] = mapped_column(String(100), primary_key=True)
    pp_id: Mapped[str] = mapped_column(String(50), index=True)
    reference_proteome : Mapped[str] = mapped_column(String(50))
    byte_offset: Mapped[int] = mapped_column(BigInteger)      # byte offset
    length: Mapped[int] = mapped_column(BigInteger)
    taxid: Mapped[int] = mapped_column(Integer, index=True)

class PPftp(Base):
    __tablename__ = FTP_PP
    pp_id: Mapped[str] = mapped_column(String(50), primary_key=True,index=True)
    complete_prots:  Mapped[str] = mapped_column(String(100))

def create_database(Base):
    engine = create_engine(f"sqlite:///{DB_NAME}")
    Base.metadata.create_all(engine)
    return engine

# Trasverse the url to find directories for each pan_proteome
def trasverse_ftp_dir(url):
    with FTP(url) as ftp:
        subdirs = [
            name for name, facts in ftp.mlsd("pub/databases/uniprot/knowledgebase/pan_proteomes/")
            if facts.get("type") == "dir"
        ]

        def build_json_path(dir):
            json = dir.replace("/","_stats.json")
            return f"{dir}{json}"

        json_files = list(map(build_json_path, subdirs))

        batch_size = (len(json_files)/CORES)
        batches = []
        i = len(json_files)
        start = 0

        while i > 0:
            if i >= batch_size:
                end = start + batch_size
            else:
                end = start + i
            batch = (start, end)
            batches.append(batch)
            i = i - batch_size
            start = end

        with mp.Pool(processes = CORES, initializer=init_worker) as pool:
            results = pool.starmap_async(worker_execute_batches, batches)

        result_dict = {panproteome : ref_prots for (panproteome, ref_prots) in list(results)}
        return result_dict


def worker_reader_json(file, panproteome, ftp):
    refprots = []
    buffer = io.BytesIO()
    ftp.retrbinary(f"RETR {file}", buffer.write)
    buffer.seek(0)
    with gzip.open(buffer, mode="rt", encoding="utf-8") as gz_file:
        data = json.load(gz_file)
        for refprot in data.get("pan_proteome", []).get("rp_composition", []):
            proteome_id = refprot.get("rp_id", None)
            if proteome_id:
                refprots.append(proteome_id)

    return (panproteome, refprots)

def worker_execute_batches(batch:tuple, json_files):
    start, end = batch
    for i in range(start,end):
        panproteome = re.search(r"/([a-z0-9]+_stats\.json)",json_files[i])
        if panproteome:
            worker_reader_json(json_files[i], panproteome[1], ftp=ftp)

def init_worker(ftp, json_files):
    global ftp_global
    global jsons
    ftp_global = ftp
    jsons = json_files

def build_fasta_db(DB_class: type[PanProteomes] | type[CompleteProteomes], fasta_file : str, session : type[Session], db : str) -> None:
    seq_dict = None
    try:
        seq_dict = SeqIO.index(fasta_file, "fasta")
        proxy = seq_dict._proxy
        rows = []

        for key in seq_dict.keys():
            record = seq_dict[key]
            header = record.description
            offset, seq_length = proxy.get(key)
            if db == FASTA_COMPLETE:
                match = proteome_id_regex.search(header)
                if match:
                    rows.append(DB_class(id= record.id,
                    proteome_id = match[2],
                    byte_offset = offset,
                    length = seq_length,
                    taxid = match[1] ))
                else:
                   continue
            elif db == FASTA_PP:
                match = panproteome_id_regex.search(header)
                if match:
                    rows.append(DB_class(id = record.id,
                    pp_id = match[1],
                    reference_proteome = match[2],
                    byte_offset = offset,   # byte offset
                    length = seq_length,
                    taxid = match[3]))
                else:
                   continue

        session.add_all(rows)
        session.commit()
    except Exception as e:
        print(f"Error parsing fasta into database: {e}")

    finally:
        if seq_dict:
            seq_dict.close()

def populate_pp_json_db(json_dict, pp_db : type[PPftp], session : Session)-> None:
    rows = []
    for key in json_dict:
        for ref_prot in json_dict[key]:
            rows.append(pp_db(pp_id = str(key), complete_prots = ref_prot))

    session.add_all(rows)
    session.commit()


def main():

    ## config
    fasta_pp_path = ""

    fasta_complete_path = ""

    result_table = "combined_pp_complete_prot_db"

    if not os.path.isfile(fasta_pp_path) or not os.path.isfile(fasta_complete_path):
        print("Could not find Panproteomes or complete proteomes file. exiting...")
        return

    result_dict = trasverse_ftp_dir(PP_FTP)
    engine = create_database(Base)

    with Session(engine) as session:
        metadata = MetaData()
        populate_pp_json_db(result_dict, PPftp, session)
        build_fasta_db(CompleteProteomes, fasta_complete_path, session, FASTA_COMPLETE)
        build_fasta_db(PanProteomes, fasta_pp_path, session, FASTA_PP)

        # Use complete_proteomes -> panproteome information to exclude complete proteomes sequences that have a panproteome
        merge_complete_no_pp = select(CompleteProteomes, PPftp).join(CompleteProteomes, CompleteProteomes.proteome_id != PPftp.complete_prots)

        results = session.execute(merge_complete_no_pp).all()
        # Stack sequences from complete proteomes (no pp) and panproteomes (substituting sequences that have a panproteome)
        for filtered, _ in results:
            part1=select(filtered.id, filtered.byte_offset, filtered.length)
            part2=select(PanProteomes.id, PanProteomes.byte_offset, PanProteomes.length)

            combined = unionall(part1,part2)

            combined_table = Table(result_table, metadata,  Column("id", String(100), primary_key=True, index = True),
               Column("byte_offset", BigInteger),
               Column("length", BigInteger))

            metadata.create_all(engine)

            insert_stmt = insert(combined_table).from_select(
                ["id", "byte_offset", "length"],
                combined
            )
            session.execute(insert_stmt)
            session.commit()

if __name__ == "__main__":
    main()
