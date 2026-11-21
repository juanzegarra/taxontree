# Author: Juan Zegarra
# Date 19/11/26

###################################################################

from ftplib import FTP
from multiprocessing import process
import re
import json
import gzip
import io
import multiprocessing as mp
import sqlite3
import pandas as pd


PP_FTP = "ftp.uniprot.org"
CORES = 2
DB_NAME = "pp_prot.db"
FASTA_PP = ""
FASTA_WHOLE  = ""
proteome_id_regex = re.compile(r"\bUP[0-9]+\b")

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


def worker_reader(file, panproteome, ftp):
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

def worker_execute_batches(batch:tuple, json_files = jsons):
    start, end = batch
    for i in range(start,end):
        panproteome = re.search(r"/([a-z0-9]+_stats\.json)",json_files[i])
        if panproteome:
            worker_reader(json_files[i], panproteome[1], ftp=ftp)

def init_worker(ftp, json_files):
    global ftp_global
    global jsons
    ftp_global = ftp
    jsons = json_files

def build_fasta_db(db, fasta_file, cursor):
    cursor.execute( f"""CREATE TABLE IF NOT EXISTS {db} (
                   key TEXT NOT NULL,
                   item TEXT NOT NULL,
                   PRIMARY KEY (key, item)
               ) WITHOUT ROWID
           """)
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                match = proteome_id_regex.search(line)
                if match:
                    match[1]

def uniprot_worker():

def fasta_to_proteome(fasta_file, fasta_db, conn):
    records = SeqIO.parse()

def main():

    result_dict = trasverse_ftp_dir(PP_FTP)

    with sqlite3.connect(DB_NAME) as conn:
        cursor = conn.cursor()
        cursor.execute( """CREATE TABLE IF NOT EXISTS pp_prot (
                       key TEXT NOT NULL,
                       item TEXT NOT NULL,
                       PRIMARY KEY (key, item)
                   ) WITHOUT ROWID
               """)
        rows = [(k, item) for k, items in result_dict.items() for item in items]
        cursor.executemany("INSERT OR IGNORE INTO key_sets (key, item) VALUES (?, ?)", rows)
