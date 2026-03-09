import os
import sqlite3


import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

DB_PATH     = "immune_cells.db"
OUTPUT_DIR  = "outputs"

CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
POP_COLORS = {
    "b_cell":      "#4C72B0",
    "cd8_t_cell":  "#DD8452",
    "cd4_t_cell":  "#55A868",
    "nk_cell":     "#C44E52",
    "monocyte":    "#8172B2",
}

# Data loaders 

@st.cache_data
def load_stats():
    return pd.read_csv(os.path.join(OUTPUT_DIR, "statistical_comparison.csv"))


@st.cache_data
def load_response_frequencies():
    return pd.read_csv(os.path.join(OUTPUT_DIR, "melanoma_pbmc_response_frequencies.csv"))


@st.cache_data
def load_cohort():
    return pd.read_csv(os.path.join(OUTPUT_DIR, "melanoma_miraclib_subset.csv"))


@st.cache_data
def total_sample_count():
    conn = sqlite3.connect(DB_PATH)
    result = pd.read_sql_query("SELECT COUNT(DISTINCT sample_id) AS n FROM samples", conn)
    conn.close()
    return int(result["n"].iloc[0])


def query_frequency(projects, conditions, treatments, sample_types, populations):
    placeholders = lambda lst: ",".join(f"'{v}'" for v in lst)
    sql = f"""
        SELECT
            spf.sample,
            spf.total_count,
            spf.population,
            spf.count,
            spf.percentage
        FROM sample_population_frequencies spf
        JOIN samples sa   ON sa.sample_id   = spf.sample
        JOIN subjects sub ON sub.subject_id = sa.subject_id
        WHERE sub.project_id  IN ({placeholders(projects)})
          AND sub.condition   IN ({placeholders(conditions)})
          AND sub.treatment   IN ({placeholders(treatments)})
          AND sa.sample_type  IN ({placeholders(sample_types)})
          AND spf.population  IN ({placeholders(populations)})
    """
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(sql, conn)
    conn.close()
    return df


# Page config 

st.set_page_config(
    page_title="Immune Cell Dashboard",
    layout="wide",
)
st.title("Immune Cell Dashboard")

tab1, tab2, tab3 = st.tabs([
    "Frequency Overview",
    "Statistical Comparison",
    "Cohort Breakdown",
])

# Tab 1: Frequency Overview (Part 2) 

with tab1:
    st.header("Cell Population Frequencies by Sample")

    with st.sidebar:
        st.header("Filters")
        sel_projects   = st.multiselect("Project",     ["prj1", "prj2", "prj3"],                        default=["prj1", "prj2", "prj3"])
        sel_conditions = st.multiselect("Condition",   ["melanoma", "carcinoma", "healthy"],             default=["melanoma", "carcinoma", "healthy"])
        sel_treatments = st.multiselect("Treatment",   ["miraclib", "phauximab", "none"],                default=["miraclib", "phauximab", "none"])
        sel_types      = st.multiselect("Sample type", ["PBMC", "WB"],                                   default=["PBMC", "WB"])
        sel_pops       = st.multiselect("Population",  CELL_POPULATIONS,                                 default=CELL_POPULATIONS)

    if not all([sel_projects, sel_conditions, sel_treatments, sel_types, sel_pops]):
        st.warning("Select at least one value in each filter.")
    else:
        df = query_frequency(sel_projects, sel_conditions, sel_treatments, sel_types, sel_pops)

        st.subheader("Relative Frequency per Sample")

        view_mode = st.radio(
            "view", ["Summary Table", "Frequency Chart"],
            horizontal=True, label_visibility="collapsed",
        )

        if view_mode == "Summary Table":
            n_filtered = df["sample"].nunique()
            n_total    = total_sample_count()
            st.info(f"Showing **{n_filtered:,}** of **{n_total:,}** samples. Use sidebar filters to narrow results.")
            st.dataframe(
                df[["sample", "total_count", "population", "count", "percentage"]],
                use_container_width=True,
            )
        else:
            # Pivot to wide: rows = samples, cols = populations
            pivot = df.pivot_table(index="sample", columns="population", values="percentage", aggfunc="first").reset_index()

            # Attach total_count (one value per sample)
            totals = df.drop_duplicates("sample")[["sample", "total_count"]]
            pivot = pivot.merge(totals, on="sample")

            # Cap at 500 samples for performance
            MAX_SAMPLES = 500
            if len(pivot) > MAX_SAMPLES:
                st.info(f"Showing first {MAX_SAMPLES} of {len(pivot)} samples. Apply filters to narrow the view.")
                pivot = pivot.head(MAX_SAMPLES)

            ROW_HEIGHT = 22   # px per row - controls bar height and figure height
            X_SAMPLE   = -30  # x position for sample ID text column
            X_TOTAL    = -8   # x position for total count text column

            pops_present = [p for p in CELL_POPULATIONS if p in pivot.columns]
            n_rows = len(pivot)

            # Use integer y positions for reliable row alignment and height control
            # Categorical y-axis in Plotly doesn't respect figure height predictably
            y_pos = list(range(n_rows))

            fig = go.Figure()

            # Text column: Sample ID
            fig.add_trace(go.Scatter(
                x=[X_SAMPLE] * n_rows,
                y=y_pos,
                mode="text",
                text=pivot["sample"].tolist(),
                textposition="middle center",
                textfont=dict(size=9, color="#333"),
                showlegend=False,
                hoverinfo="skip",
            ))

            # Text column: Total count
            fig.add_trace(go.Scatter(
                x=[X_TOTAL] * n_rows,
                y=y_pos,
                mode="text",
                text=pivot["total_count"].apply(lambda v: f"{int(v):,}").tolist(),
                textposition="middle center",
                textfont=dict(size=9, color="#333"),
                showlegend=False,
                hoverinfo="skip",
            ))

            # Stacked bars (0–100)
            for pop in pops_present:
                fig.add_trace(go.Bar(
                    name=pop,
                    y=y_pos,
                    x=pivot[pop].tolist(),
                    orientation="h",
                    marker_color=POP_COLORS[pop],
                    text=pivot[pop].round(1).astype(str) + "%",
                    textposition="inside",
                    insidetextanchor="middle",
                    textfont=dict(size=9, color="white"),
                    hovertemplate=f"<b>{pop}</b>: %{{x:.1f}}%<extra></extra>",
                ))

            fig.update_layout(
                barmode="stack",
                xaxis=dict(
                    range=[-40, 106],
                    tickvals=[0, 25, 50, 75, 100],
                    ticktext=["0%", "25%", "50%", "75%", "100%"],
                    title="Relative frequency (%)",
                    side="top",
                ),
                yaxis=dict(
                    range=[n_rows - 0.5, -0.5],  # explicit reversed range: row 0 at top
                    showticklabels=False,
                    showgrid=False,
                ),
                height=max(400, n_rows * ROW_HEIGHT + 80),
                legend_title="Population",
                margin=dict(l=10, r=20, t=60, b=10),
                annotations=[
                    dict(x=X_SAMPLE, y=1.05, xref="x", yref="paper",
                         text="<b>Sample</b>", showarrow=False,
                         font=dict(size=11, color="#222"), xanchor="center"),
                    dict(x=X_TOTAL, y=1.05, xref="x", yref="paper",
                         text="<b>Total count</b>", showarrow=False,
                         font=dict(size=11, color="#222"), xanchor="center"),
                    dict(x=50, y=1.05, xref="x", yref="paper",
                         text="<b>Cell population frequencies</b>", showarrow=False,
                         font=dict(size=11, color="#222"), xanchor="center"),
                ],
            )

            st.plotly_chart(fig, use_container_width=True)

# Tab 2: Statistical Comparison (Part 3) 

with tab2:
    st.header("Statistical Comparison: Responders vs Non-responders")
    st.caption("Cohort: Melanoma · miraclib · PBMC · Baseline (day 0) | Test: Mann-Whitney U · Bonferroni correction (×5)")

    stats_df = load_stats()
    freq_df  = load_response_frequencies()

    # Findings summary 
    st.subheader("Findings")
    sig_pops = stats_df[stats_df["significant"] == "yes"]["population"].tolist()
    if sig_pops:
        st.success(
            f"**Significant difference found** in: {', '.join(sig_pops)}. "
            "These populations may have predictive value for miraclib response."
        )
    else:
        st.info(
            "**No statistically significant differences detected** across all five cell populations "
            "(all p_corrected = 1.00 after Bonferroni correction). "
            "Baseline immune cell composition does not predict miraclib response in this cohort."
        )

    st.subheader("Statistical Results")
    st.dataframe(
        stats_df,
        column_config={
            "responders_median":    st.column_config.NumberColumn(format="%.2f"),
            "nonresponders_median": st.column_config.NumberColumn(format="%.2f"),
            "p_value":              st.column_config.NumberColumn(format="%.4f"),
            "p_corrected":          st.column_config.NumberColumn(format="%.2f"),
        },
        use_container_width=True,
    )

    st.subheader("Boxplot: Cell Population Frequencies")

    RESPONDER_COLOR    = "#4C72B0"
    NONRESPONDER_COLOR = "#DD8452"

    # Build titles with p-values upfront and pass to make_subplots directly
    subplot_titles = [
        f"{pop}<br>p={stats_df[stats_df['population'] == pop]['p_corrected'].values[0]:.2f}"
        for pop in CELL_POPULATIONS
    ]

    fig2 = make_subplots(
        rows=1, cols=5,
        subplot_titles=subplot_titles,
        shared_yaxes=False,
    )

    for i, pop in enumerate(CELL_POPULATIONS, start=1):
        pop_df   = freq_df[freq_df["population"] == pop]
        yes_vals = pop_df[pop_df["response"] == "yes"]["percentage"]
        no_vals  = pop_df[pop_df["response"] == "no"]["percentage"]

        fig2.add_trace(go.Box(
            y=yes_vals,
            name="Responders",
            marker_color=RESPONDER_COLOR,
            boxpoints="outliers",  # only show outlier points — box stays visible
            pointpos=-1.5,         # offset points to left of box
            showlegend=(i == 1),
            legendgroup="Responders",
            hovertemplate="Responders<br>%{y:.2f}%<extra></extra>",
        ), row=1, col=i)

        fig2.add_trace(go.Box(
            y=no_vals,
            name="Non-responders",
            marker_color=NONRESPONDER_COLOR,
            boxpoints="outliers",
            pointpos=-1.5,
            showlegend=(i == 1),
            legendgroup="Non-responders",
            hovertemplate="Non-responders<br>%{y:.2f}%<extra></extra>",
        ), row=1, col=i)

    fig2.update_layout(
        height=500,
        boxmode="group",
        legend=dict(orientation="h", y=-0.15),
        margin=dict(t=60, b=60),
    )
    fig2.update_yaxes(title_text="Relative frequency (%)", col=1)

    st.plotly_chart(fig2, use_container_width=True)


# Tab 3: Cohort Breakdown (Part 4)

with tab3:
    st.header("Cohort Breakdown")
    st.caption("Melanoma · miraclib · PBMC · Baseline (day 0)")

    cohort_df = load_cohort()

    by_project  = cohort_df[cohort_df["category"] == "by_project"]
    by_response = cohort_df[cohort_df["category"] == "by_response"]
    by_sex      = cohort_df[cohort_df["category"] == "by_sex"]
    avg_b_cell  = cohort_df[cohort_df["category"] == "avg_b_cell_male_responders"]["value"].iloc[0]

    col1, col2, col3 = st.columns(3)

    def small_bar(df_slice, title, color):
        fig = go.Figure(go.Bar(
            x=df_slice["group"],
            y=df_slice["value"],
            marker_color=color,
            text=df_slice["value"].astype(int),
            textposition="outside",
            cliponaxis=False,
        ))
        fig.update_layout(
            title=title,
            height=300,
            margin=dict(t=50, b=20, l=20, r=20),
            yaxis=dict(title="Sample count"),
        )
        return fig

    with col1:
        st.plotly_chart(small_bar(by_project,  "Samples per Project",            "#4C72B0"), use_container_width=True)
    with col2:
        st.plotly_chart(small_bar(by_response, "Responders vs Non-responders",   "#55A868"), use_container_width=True)
    with col3:
        st.plotly_chart(small_bar(by_sex,      "Males vs Females",               "#C44E52"), use_container_width=True)

    st.metric(
        label="Avg B cells — melanoma male responders at baseline (day 0)",
        value=f"{avg_b_cell:,.2f}",
    )
