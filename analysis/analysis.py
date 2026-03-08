import sqlite3
import pandas as pd
import os
from scipy import stats

DB_PATH     = "immune_cells.db"
OUTPUT_DIR  = "outputs"

CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

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

def get_part3_data(conn):
    """
    Query filtered data for Part 3 statistical analysis.
    Returns percentages per population for melanoma PBMC miraclib patients at baseline.
    """
    query = """
        SELECT
            spf.sample,
            spf.population,
            spf.percentage,
            sub.response
        FROM sample_population_frequencies spf
        JOIN samples sa   ON sa.sample_id    = spf.sample
        JOIN subjects sub ON sub.subject_id  = sa.subject_id
        WHERE sub.condition                = 'melanoma'
          AND sub.treatment               = 'miraclib'
          AND sa.sample_type              = 'PBMC'
          AND sa.time_from_treatment_start = 0
          AND sub.response IN ('yes', 'no')
    """
    df = pd.read_sql_query(query, conn)
    return df

def check_normality(df):
    """
    Run Shapiro-Wilk normality test for each population and response group.
    Returns a dict: { population: {"yes": p_value, "no": p_value} }
    Answers the question: "Is the percentage data for each cell population normally distributed?"
    """
    results = {}
    for pop in CELL_POPULATIONS:
        pop_df = df[df["population"] == pop]
        yes_vals = pop_df[pop_df["response"] == "yes"]["percentage"]
        no_vals  = pop_df[pop_df["response"] == "no"]["percentage"]
        _, p_yes = stats.shapiro(yes_vals)
        _, p_no  = stats.shapiro(no_vals)
        results[pop] = {"yes": p_yes, "no": p_no}
    return results
    
def run_statistical_tests(df):
    """
    Run Mann-Whitney U test for each cell population comparing
    responders vs non-responders. Applies Bonferroni correction.
    Returns a list of dicts with results per population.
    """
    n_tests = len(CELL_POPULATIONS)
    rows = []

    for pop in CELL_POPULATIONS:
        pop_df   = df[df["population"] == pop]
        yes_vals = pop_df[pop_df["response"] == "yes"]["percentage"]
        no_vals  = pop_df[pop_df["response"] == "no"]["percentage"]

        _, p_raw = stats.mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
        p_corr   = min(p_raw * n_tests, 1.0)

        rows.append({
            "population":          pop,
            "responders_median":   round(yes_vals.median(), 4),
            "nonresponders_median": round(no_vals.median(), 4),
            "p_value":             round(p_raw, 4),
            "p_corrected":         round(p_corr, 4),
            "significant":         "yes" if p_corr < 0.05 else "no"
        })

    return rows

def run_statistical_comparison(conn):
    print("\nPart 3: Statistical analysis")
    df      = get_part3_data(conn)
    results = run_statistical_tests(df)
    stats_df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "statistical_comparison.csv")
    stats_df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")
    print(stats_df.to_string(index=False))
    return df, stats_df

# Main

def main():
    ensure_output_dir()
    conn = get_connection()
    try:
        run_freq_table(conn)
        run_statistical_comparison(conn)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
