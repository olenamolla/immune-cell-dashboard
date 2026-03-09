import sqlite3
import pandas as pd
import os
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    df.to_csv(os.path.join(OUTPUT_DIR, "melanoma_pbmc_response_frequencies.csv"), index=False)
    results = run_statistical_tests(df)
    stats_df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "statistical_comparison.csv")
    stats_df.to_csv(out_path, index=False)
    print(f"Saved to {out_path}")
    print(stats_df.to_string(index=False))
    save_boxplot(df, stats_df)
    return df, stats_df

def save_boxplot(df, stats_df):
    """
    Save boxplot PNG comparing responders vs non-responders per cell population.
    """
    fig, axes = plt.subplots(1, 5, figsize=(18, 6), sharey=False)
    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-responders\n"
        "Melanoma · miraclib · PBMC · Baseline (day 0)",
        fontsize=13, y=1.02
    )

    for i, pop in enumerate(CELL_POPULATIONS):
        ax      = axes[i]
        pop_df  = df[df["population"] == pop]
        yes_pct = pop_df[pop_df["response"] == "yes"]["percentage"].values
        no_pct  = pop_df[pop_df["response"] == "no"]["percentage"].values

        bp = ax.boxplot([yes_pct, no_pct], patch_artist=True, widths=0.5)
        bp["boxes"][0].set(facecolor="#4C72B0", alpha=0.7)
        bp["boxes"][1].set(facecolor="#DD8452", alpha=0.7)

        ax.scatter([1] * len(yes_pct), yes_pct, color="#4C72B0", alpha=0.3, s=10, zorder=3)
        ax.scatter([2] * len(no_pct),  no_pct,  color="#DD8452", alpha=0.3, s=10, zorder=3)

        p_corr = stats_df[stats_df["population"] == pop]["p_corrected"].values[0]
        ax.set_title(f"{pop}\np={p_corr}", fontsize=10)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Responders", "Non-responders"], fontsize=8)
        ax.set_ylabel("Relative frequency (%)" if i == 0 else "")

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "statistical_comparison_boxplot.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved boxplot to {out_path}")

# Part 4: Cohort description

def describe_melanoma_cohort(conn):
    """
    Describe the melanoma/miraclib/PBMC/day-0 cohort from four angles:
      Q1 - sample count per project
      Q2 - responders vs non-responders count
      Q3 - males vs females count
      Q4 - average b_cell absolute count for melanoma male responders at day 0 (all treatments/sample types)
    Saves results to outputs/melanoma_miraclib_subset.csv in long format:
      category | group | value
    """
    base_filter = """
        FROM samples sa
        JOIN subjects sub ON sub.subject_id = sa.subject_id
        WHERE sub.condition                 = 'melanoma'
          AND sub.treatment                 = 'miraclib'
          AND sa.sample_type                = 'PBMC'
          AND sa.time_from_treatment_start  = 0
    """

    q1 = pd.read_sql_query(f"""
        SELECT sub.project_id AS grp, COUNT(DISTINCT sa.sample_id) AS value
        {base_filter}
        GROUP BY sub.project_id
    """, conn)
    q1.insert(0, "category", "by_project")

    q2 = pd.read_sql_query(f"""
        SELECT sub.response AS grp, COUNT(DISTINCT sa.sample_id) AS value
        {base_filter}
          AND sub.response IN ('yes', 'no')
        GROUP BY sub.response
    """, conn)
    q2.insert(0, "category", "by_response")

    q3 = pd.read_sql_query(f"""
        SELECT sub.sex AS grp, COUNT(DISTINCT sa.sample_id) AS value
        {base_filter}
        GROUP BY sub.sex
    """, conn)
    q3.insert(0, "category", "by_sex")

    q4 = pd.read_sql_query("""
        SELECT ROUND(AVG(cc.count), 2) AS value
        FROM cell_counts cc
        JOIN samples sa   ON sa.sample_id   = cc.sample_id
        JOIN subjects sub ON sub.subject_id = sa.subject_id
        WHERE sub.condition                = 'melanoma'
          AND sa.time_from_treatment_start = 0
          AND sub.response                = 'yes'
          AND sub.sex                     = 'M'
          AND cc.population               = 'b_cell'
    """, conn)
    q4.insert(0, "category", "avg_b_cell_male_responders")
    q4.insert(1, "grp", "overall")

    result = pd.concat([q1, q2, q3, q4], ignore_index=True)
    result.columns = ["category", "group", "value"]

    out_path = os.path.join(OUTPUT_DIR, "melanoma_miraclib_subset.csv")
    result.to_csv(out_path, index=False)

    print("\nPart 4: Melanoma cohort breakdown")
    print(result.to_string(index=False))
    print(f"\nSaved to {out_path}")
    return result


# Main

def main():
    ensure_output_dir()
    conn = get_connection()
    try:
        run_freq_table(conn)
        run_statistical_comparison(conn)
        describe_melanoma_cohort(conn)
    finally:
        conn.close()


if __name__ == "__main__":
    main()
