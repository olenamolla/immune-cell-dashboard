import sqlite3
import pandas as pd
import os


CSV_PATH = "data/cell-count.csv"
DB_PATH  = "immune_cells.db"

# The five cell populations 
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]



SCHEMA = """
CREATE TABLE IF NOT EXISTS projects (
    project_id  TEXT PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT NOT NULL,
    condition   TEXT,
    age         INTEGER,
    sex         TEXT,
    treatment   TEXT,
    response    TEXT,
    FOREIGN KEY (project_id) REFERENCES projects (project_id)
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                   TEXT PRIMARY KEY,
    subject_id                  TEXT NOT NULL,
    sample_type                 TEXT,
    time_from_treatment_start   INTEGER,
    FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
);

CREATE TABLE IF NOT EXISTS cell_counts (
    sample_id   TEXT    NOT NULL,
    population  TEXT    NOT NULL,
    count       INTEGER NOT NULL,
    PRIMARY KEY (sample_id, population),
    FOREIGN KEY (sample_id) REFERENCES samples (sample_id)
);

CREATE VIEW IF NOT EXISTS sample_population_frequencies AS
SELECT
    cc.sample_id                                                                   AS sample,
    SUM(cc.count) OVER (PARTITION BY cc.sample_id)                                AS total_count,
    cc.population,
    cc.count,
    ROUND(cc.count * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample_id), 2)  AS percentage
FROM cell_counts cc
ORDER BY cc.sample_id, cc.population;
"""


def create_schema(conn):
    """Create all four tables if they don't already exist."""
    conn.executescript(SCHEMA)
    conn.commit()
    print("Schema created.")


def load_csv():
    """Read the CSV into a pandas DataFrame."""
    df = pd.read_csv(CSV_PATH)
    print(f"Loaded {len(df)} rows from {CSV_PATH}")
    return df


def insert_projects(conn, df):
    """Insert unique project IDs."""
    projects = df[["project"]].drop_duplicates().rename(columns={"project": "project_id"})
    projects.to_sql("projects", conn, if_exists="append", index=False)
    print(f"Inserted {len(projects)} projects.")


def insert_subjects(conn, df):
    """
    Insert one row per unique subject.
    """
    subject_cols = ["subject", "project", "condition", "age", "sex", "treatment", "response"]
    subjects = (
        df[subject_cols]
        .drop_duplicates(subset=["subject"])
        .rename(columns={"subject": "subject_id", "project": "project_id"})
    )
    subjects.to_sql("subjects", conn, if_exists="append", index=False)
    print(f"Inserted {len(subjects)} subjects.")


def insert_samples(conn, df):
    """Insert one row per sample (the collection event at a given timepoint)."""
    sample_cols = ["sample", "subject", "sample_type", "time_from_treatment_start"]
    samples = df[sample_cols].rename(columns={"sample": "sample_id", "subject": "subject_id"})
    samples.to_sql("samples", conn, if_exists="append", index=False)
    print(f"Inserted {len(samples)} samples.")


def insert_cell_counts(conn, df):
    """
    Convert the wide-format CSV (one column per population) to long format
    (one row per sample + population), then insert into cell_counts.

    pandas melt() does the wide -> long conversion:
      id_vars   = the columns we keep as-is (sample_id)
      value_vars = the columns we want to unpivot (the 5 populations)
      var_name  = name for the new 'which column was this' column
      value_name= name for the actual numeric value
    """
    cell_counts = df[["sample"] + POPULATIONS].melt(
        id_vars="sample",
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count"
    ).rename(columns={"sample": "sample_id"})

    cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)
    print(f"Inserted {len(cell_counts)} cell count rows.")


def main():
    # Remove old DB if it exists
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print(f"Removed existing {DB_PATH}")

    conn = sqlite3.connect(DB_PATH)

    try:
        create_schema(conn)
        df = load_csv()
        insert_projects(conn, df)
        insert_subjects(conn, df)
        insert_samples(conn, df)
        insert_cell_counts(conn, df)
        conn.commit()
        print(f"\nDone. Database saved to {DB_PATH}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
