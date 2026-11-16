import sqlite3
import re
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import sys

def setup_database(db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    # Create simple, flat tables
    c.executescript("""
        CREATE TABLE IF NOT EXISTS taxonomy (
            uniprot_ac TEXT PRIMARY KEY,
            genus TEXT,
            species TEXT
        );
        CREATE TABLE IF NOT EXISTS uniref_map (
            uniprot_ac TEXT PRIMARY KEY,
            uniref100 TEXT,
            uniref90 TEXT,
            uniref50 TEXT
        );
    """)
    conn.commit()
    return conn

def load_fasta_to_db(fasta_file, conn, batch_size=10000):
    print("Loading FASTA taxonomy into SQLite...")
    pattern = re.compile(r"OS=(.*?)(?=\s+(?:OX|GN|PE|SV)=|\s*;|$)")
    c = conn.cursor()

    data_batch = []

    with open(fasta_file, "r") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            # Extract actual UniProt Accession from header (e.g., >sp|P12345|NAME -> P12345)
            parts = rec.id.split('|')
            uniprot_ac = parts[1] if len(parts) > 1 else rec.id

            sp_match = pattern.search(rec.description)
            if sp_match:
                species = sp_match.group(1).strip()
                genus = species.split(" ")[0]
                data_batch.append((uniprot_ac, genus, species))

            if len(data_batch) >= batch_size:
                c.executemany("INSERT OR IGNORE INTO taxonomy (uniprot_ac, genus, species) VALUES (?, ?, ?)", data_batch)
                data_batch = []

        # Insert any remaining records
        if data_batch:
            c.executemany("INSERT OR IGNORE INTO taxonomy (uniprot_ac, genus, species) VALUES (?, ?, ?)", data_batch)

    conn.commit()

def load_idmap_to_db(idmapping_file, conn):
    print("Loading ID mapping file in chunks...")
    # Removed the trailing numbers (1, 2) from your columns for cleaner DB headers
    idmapping_columns = [
        "UniProtKB-AC", "UniProtKB-ID", "GeneID", "RefSeq", "GI", "PDB", "GO", "UniRef100",
        "UniRef90", "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM", "UniGene", "PubMed",
        "EMBL", "EMBL-CDS", "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional_PubMed"
    ]

    # Only keep the columns we actually need for the join to save massive amounts of space
    cols_to_use = ["UniProtKB-AC", "UniRef100", "UniRef90", "UniRef50"]

    # Process in chunks of 100,000 rows to prevent RAM crashes
    chunk_size = 100000
    for chunk in pd.read_csv(idmapping_file, sep="\t", header=None, names=idmapping_columns, usecols=cols_to_use, chunksize=chunk_size):
        # Rename columns to match SQLite schema
        chunk.columns = ["uniprot_ac", "uniref100", "uniref90", "uniref50"]
        chunk.to_sql("uniref_map", conn, if_exists="append", index=False)

    # Create indexes to make the JOIN extremely fast
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tax_ac ON taxonomy(uniprot_ac);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_map_ac ON uniref_map(uniprot_ac);")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_genus ON taxonomy(genus);")
    conn.commit()

# --- Configuration ---
fasta_file = Path("your_proteome.fasta")
idmapping_file = Path("idmapping.dat")
db_path = "proteome_mapping.db"

if not fasta_file.exists() or not idmapping_file.exists():
    print("Fasta or idmapping files not found. Update script paths.")
    sys.exit(1)

# --- Execution ---
conn = setup_database(db_path)

# Check if tables are empty before loading
c = conn.cursor()
c.execute("SELECT COUNT(*) FROM taxonomy")
if c.fetchone()[0] == 0:
    load_fasta_to_db(fasta_file, conn)

c.execute("SELECT COUNT(*) FROM uniref_map")
if c.fetchone()[0] == 0:
    load_idmap_to_db(idmapping_file, conn)

print("Database built successfully.")

# --- Calculate reduction at the SPECIES level ---
print("Calculating unique UniRef counts per species...")
query_species_summary = """
    SELECT
        t.genus,
        t.species,
        COUNT(DISTINCT t.uniprot_ac) as total_original_proteins,
        COUNT(DISTINCT u.uniref100) as unique_uniref100,
        COUNT(DISTINCT u.uniref90) as unique_uniref90,
        COUNT(DISTINCT u.uniref50) as unique_uniref50
    FROM taxonomy t
    LEFT JOIN uniref_map u ON t.uniprot_ac = u.uniprot_ac
    GROUP BY t.genus, t.species
"""
species_summary_df = pd.read_sql_query(query_species_summary, conn)

# Export to CSV so you can review the sequence reduction for every species
species_summary_df.to_csv("uniref_reduction_by_species.csv", index=False)
print("Saved species-level reduction to 'uniref_reduction_by_species.csv'")


# --- Calculate reduction at the GENUS level ---
print("Calculating unique UniRef counts per genus...")
query_genus_summary = """
    SELECT
        t.genus,
        COUNT(DISTINCT t.uniprot_ac) as total_original_proteins,
        COUNT(DISTINCT u.uniref100) as unique_uniref100,
        COUNT(DISTINCT u.uniref90) as unique_uniref90,
        COUNT(DISTINCT u.uniref50) as unique_uniref50
    FROM taxonomy t
    LEFT JOIN uniref_map u ON t.uniprot_ac = u.uniprot_ac
    GROUP BY t.genus
"""
genus_summary_df = pd.read_sql_query(query_genus_summary, conn)

# Export to CSV so you can review the sequence reduction for every genus
genus_summary_df.to_csv("uniref_reduction_by_genus.csv", index=False)
print("Saved genus-level reduction to 'uniref_reduction_by_genus.csv'")

# Close the database connection
conn.close()
