import sqlite3
import pandas as pd
import os

DB_PATH     = "immune_cells.db"
OUTPUT_DIR  = "outputs"


def get_connection():
    return sqlite3.connect(DB_PATH)


def ensure_output_dir():
    os.makedirs(OUTPUT_DIR, exist_ok=True)


# Frequency table 

def compute_frequency_table(conn):
    """
    Query the pre-built view to get relative frequency of each cell population
    per sample. Returns a DataFrame with columns:
      sample | total_count | population | count | percentage
    """
    query = "SELECT sample, total_count, population, count, percentage FROM sample_population_frequencies"
    df = pd.read_sql_query(query, conn)
    return df


def run_freq_table(conn):
    print("\nPart 2: Frequency table")
    df = compute_frequency_table(conn)
    out_path = os.path.join(OUTPUT_DIR, "frequency_table.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} rows to {out_path}")
    print(df.head(10).to_string(index=False))
    return df


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ensure_output_dir()
    conn = get_connection()
    try:
        run_freq_table(conn)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
