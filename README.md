# Immune Cell Dashboard

Interactive dashboard for exploring immune cell population data from a clinical trial. Built with Python, SQLite, and Streamlit.

**Live dashboard:** _[link to be added after deployment]_

---

## Quickstart (GitHub Codespaces)

```bash
make setup      # install dependencies
make pipeline   # build database + generate all outputs
make dashboard  # launch Streamlit on port 8501
```

All outputs are written to `outputs/`. The SQLite database (`immune_cells.db`) is regenerated from scratch on every `make pipeline` run.


---

## Database Schema

The relational schema uses **four tables** and one **SQL view**:

```
projects   (project_id PK)
subjects   (subject_id PK, project_id FK, condition, age, sex, treatment, response)
samples    (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
cell_counts (sample_id FK, population, count  — composite PK on both columns)

VIEW: sample_population_frequencies
      (sample, total_count, population, count, percentage)
```

### Design rationale

**Normalization.** Patient-level attributes (condition, treatment, response, sex) live in `subjects`, not repeated on every sample row. Timepoint and sample-type metadata live in `samples`.

**Long-format `cell_counts`.** The raw CSV stores one column per cell population (wide format). The database stores one row per `(sample, population)` pair instead. This is the scalability decision: adding a new cell population requires zero schema changes — it is just a new value in the `population` column. A wide-format table would require an ALTER TABLE for every new biomarker. At thousands of samples and dozens of populations, long format also makes aggregate queries (filter by population, compare across populations) simpler and avoids wide table scans.

**SQL VIEW for frequencies.** The relative frequency calculation (count / sum-over-sample × 100) is defined once in the `sample_population_frequencies` view using a window function. Every downstream query: Part 2 summary table, Part 3 statistical analysis, the dashboard — reads from this view. There is a single source of truth for the percentage calculation.


### Scalability

With hundreds of projects, thousands of subjects, and varied analytics:

- **Indexes** on `cell_counts.population`, `samples.subject_id`, `subjects.project_id`, and `samples.time_from_treatment_start` would keep analytical queries fast as row counts grow. These are omitted here because SQLite on this dataset size does not need them, but they are the first addition at scale.
- **Long-format `cell_counts`** means new biomarkers (e.g., dendritic cells, mast cells) slot in without schema migration.
- The `projects` → `subjects` → `samples` → `cell_counts` hierarchy maps cleanly to a partitioned data warehouse for parallel analytical queries at scale.
- For multi-site trials with hundreds of projects, `subjects` would gain a `site_id` foreign key and `projects` would store trial metadata. The hierarchy deepens without breaking existing queries.

---

## Code Structure

```
immune-cell-dashboard/
├── load_data.py          # Part 1: ETL — builds immune_cells.db from data/cell-count.csv
├── analysis/
│   └── analysis.py       # Parts 2–4: frequency table, statistics, cohort breakdown
├── dashboard.py          # Interactive Streamlit dashboard
├── data/
│   └── cell-count.csv    # Raw input data
├── outputs/              # Generated CSVs and plots (committed for reference)
├── requirements.txt      # Python dependencies
└── Makefile              # setup / pipeline / dashboard targets
```

**`load_data.py`** was a project requirement — it initializes and populates the database in one standalone script. It drops and recreates the database on every run (idempotent), reads the CSV with pandas, and inserts into all four tables using `to_sql`.

**`analysis/analysis.py`** handles all analytical work in three sequential steps driven by a `main()` function:
- *Part 2* - queries the frequency view and saves `frequency_table.csv` (52,500 rows).
- *Part 3* - filters to melanoma / miraclib / PBMC / day 0, checks normality with Shapiro-Wilk across all 10 groups (5 populations × 2 response groups), then runs Mann-Whitney U tests with Bonferroni correction. The statistical test choice was data-driven: Shapiro-Wilk showed 4 of 5 populations violate normality in at least one response group, ruling out a t-test. Mann-Whitney U is applied uniformly across all populations for consistency. Bonferroni correction (multiply p by number of tests, cap at 1.0) is appropriate here because there are only 5 tests - a  conservative correction is the right call at this scale.
- *Part 4* - four SQL queries against the filtered cohort produce counts by project, response, sex, and the average B cell count for melanoma male responders at baseline.

**`dashboard.py`** is a Streamlit app with three tabs mirroring Parts 2–4. It reads directly from the pre-generated CSVs in `outputs/` and queries the live database only for the interactive frequency filter in Tab 1.

---

**Average B cell count — melanoma male responders at baseline: 10,206.15**
