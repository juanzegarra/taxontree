import re
from Bio import SeqIO
import pandas as pd

names_totaxid_file = "names.dmp"
fasta_file = "all_uniprot_conc.fasta"
pattern = r"OS=([^=]+?)(?=\s+[A-Z]{2}=|$)"
fasta_to_name = []
id_to_name = []

with open(fasta_file, "r") as f:
    fastas = list(SeqIO.parse(f, "fasta"))
    for fasta in fastas:
        tax_name = re.search(pattern, str(fasta.description))
        if tax_name:
            entry = (fasta.id, tax_name.group(1))
            fasta_to_name.append(entry)

with open(names_totaxid_file, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        taxid, name = parts[0], parts[1]
        id_to_name.append((taxid, name))

df_fasta_to_name = pd.DataFrame(fasta_to_name, columns=["id", "name"])
df_id_to_name = pd.DataFrame(id_to_name, columns=["taxid", "name"])

full_df = pd.merge(df_fasta_to_name, df_id_to_name, left_on="name", right_on="name", how="inner")
