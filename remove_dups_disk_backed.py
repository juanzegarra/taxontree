# remove_dups_disk_backed.py
import sqlite3
import os

FASTA_FILE = "all_uniprot_concat.fasta"
DB_FILE = "uniprot_cache.db"
OUTPUT_FILE = f"deduplicated_renamed_{FASTA_FILE}"

def iter_fasta(filepath):
    """Generator that reads a FASTA file sequence-by-sequence (Streaming)."""
    with open(filepath, "r") as f:
        id_ = None
        ann = ""
        seq_lines = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith(">"):
                # Yield the previously accumulated sequence
                if id_ is not None:
                    yield id_, ann, "".join(seq_lines)
                
                parts = line.split(maxsplit=1)
                id_ = parts[0][1:] # Strip '>'
                ann = parts[1] if len(parts) > 1 else ""
                seq_lines = []
            else:
                seq_lines.append(line)
                
        # Yield the very last sequence in the file
        if id_ is not None:
            yield id_, ann, "".join(seq_lines)


def main():
    print("1. Initializing disk-backed SQLite database...")
    # Connect to a file-based database (creates it if it doesn't exist)
    conn = sqlite3.connect(DB_FILE)
    cursor = conn.cursor()
    
    # Create a table with a UNIQUE constraint. 
    # This automatically rejects exact duplicates of (id, annotation, sequence)
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS seqs (
            id TEXT,
            annotation TEXT,
            sequence TEXT,
            UNIQUE(id, annotation, sequence)
        )
    """)
    
    print("2. Streaming FASTA file and inserting unique records to disk...")
    batch = []
    batch_size = 100_000  # Write to disk in chunks for massive speed gains
    
    for id_, ann, seq in iter_fasta(FASTA_FILE):
        batch.append((id_, ann, seq))
        
        if len(batch) >= batch_size:
            # INSERT OR IGNORE silently drops exact duplicates
            cursor.executemany("INSERT OR IGNORE INTO seqs VALUES (?, ?, ?)", batch)
            conn.commit()
            batch.clear()
            
    # Insert any remaining records in the final batch
    if batch:
        cursor.executemany("INSERT OR IGNORE INTO seqs VALUES (?, ?, ?)", batch)
        conn.commit()

    print("3. Querying disk and writing renamed sequences...")
    # We use SQLite Window Functions to count duplicates and assign an index (rn) automatically
    query = """
        SELECT 
            id, 
            annotation, 
            sequence, 
            COUNT(*) OVER (PARTITION BY id) as total_count,
            ROW_NUMBER() OVER (PARTITION BY id ORDER BY sequence) - 1 as index_num
        FROM seqs
        ORDER BY id
    """
    
    with open(OUTPUT_FILE, "w") as out_f:
        # cursor.execute is an iterator; it streams results from disk without loading all into RAM
        for row in cursor.execute(query):
            id_, ann, seq, total_count, index_num = row
            
            # If the ID appears more than once with different sequences/annotations, append _0, _1, etc.
            final_id = f"{id_}_{index_num}" if total_count > 1 else id_
            
            out_f.write(f">{final_id} {ann}\n{seq}\n")

    # Cleanup
    conn.close()
    if os.path.exists(DB_FILE):
        os.remove(DB_FILE) # Delete the temporary database file
        
    print("Done! Output saved to:", OUTPUT_FILE)

if __name__ == "__main__":
    main()