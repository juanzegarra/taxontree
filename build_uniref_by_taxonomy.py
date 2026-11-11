import sqlite3
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from pathlib import Path
import gzip
import subprocess
import sys
from itertools import islice
import csv


def init_biosql(db_dir):
    with open(db_dir, "r") as f:
        sql_script = f.read()

    conn = sqlite3.connect(path_to_biosql_db)

    conn.executescript(sql_script)
    conn.commit()
    conn.close()
    return

def load_fasta_to_db(fasta_file, fasta_db, server,batch_size=5000):

    with open(fasta_file, "r") as handle:
            records = SeqIO.parse(handle, "fasta")

            while True:
                batch = list(islice(records, batch_size))
                if not batch:
                    break

                # Ensure molecule_type annotation exists for BioSQL compatibility
                for rec in batch:
                    rec.annotations["molecule_type"] = "protein"

                fasta_db.load(batch)
                server.commit()

def load_idmap_to_db(idmapping_file, idmapping_db, server,batch_size=5000):

     with open(idmapping_file, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                record = SeqRecord(
                    Seq(row["sequence"]),
                    id=row["id"],
                    name=row["id"],
                    description=row.get("description", ""),
                    annotations={
                        "molecule_type": "DNA",
                        "ncbi_taxid": [row["taxon_id"]] if "taxon_id" in row else []
                    }
                )
                yield record

fasta_file = Path("")
idmapping_file = Path("")
path_to_biosql_db = ""
idmapping_db_name = "idmapping_db"
fasta_db_name = "fasta_db"

conn = sqlite3.connect(path_to_biosql_db)
repo_biosql = "https://github.com/biosql/biosql.git"
idmapping_columns = ["UniProtKB-AC","UniProtKB-ID","GeneID","RefSeq","GI","PDB","GO","UniRef100",
    "UniRef90","1UniRef50","1UniParc","1PIR","1NCBI-taxon",
    "1MIM","1UniGene","1PubMed","1EMBL","1EMBL-CDS",
    "1Ensembl","2Ensembl_TRS","2Ensembl_PRO","2Additional_PubMed"]

try:
    DB = BioSeqDatabase.open_database(driver="sqlite3", database= path_to_biosql_db)

except Exception as e:
    print(f"No BioSQL database set up, trying to set up, error: {e}")
    dir = Path("/home")
    db_dir = next(dir.rglob("biosql/sql/biosqldb-sqlite.sql"), None)
    if db_dir:
        init_biosql(db_dir)
    else:
        subprocess.run(["git", "clone", repo_biosql], check = True)
        db_dir = "biosql/sql/biosqldb-sqlite.sql"
        init_biosql(db_dir)

DB = BioSeqDatabase.open_database(driver="sqlite3", database= path_to_biosql_db)

if not fasta_file.exists() or not idmapping_file.exists():
    print("Fasta or idmapping files not found, update script to match real paths for files")
    sys.exit()

if not fasta_db_name in DB:

    fasta_db = DB.new_database(fasta_db_name)
    load_fasta_to_db(fasta_file, fasta_db, DB)


elif not idmapping_db_name in DB:

    idmapping_db = DB.new_database(idmapping_db_name)
