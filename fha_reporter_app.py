"""
FHA Reporter Explorer - interactive web app
============================================

A Streamlit app for exploring the candidate table produced by the notebook
(FHA_reporter_pipeline.ipynb). It does NOT re-run the 4-minute atlas stream;
instead it loads a pre-computed candidate table and lets you filter/plot live.

------------------------------------------------------------------------------
SETUP (one time)
------------------------------------------------------------------------------
1. Install dependencies:
       pip install streamlit 
       pip install plotly 
       pip install pandas 
       pip install pyarrow

2. Export the candidate table from the notebook. After running notebook
   section 8 (which builds `df_all`), run this in a new cell:

       df_all.to_csv("candidate_table.csv", index=False)
       # also save metadata so the app knows which kinases to offer and
       # which are HEK-expressed:
       import json
       with open("app_metadata.json", "w") as f:
           json.dump({
               "output_tag": OUTPUT_TAG,
               "pathway_kinases": sorted(pathway_set),
               "all_kinases": all_atlas_kinases,
               "hek_expressed_kinases": [k for k in all_atlas_kinases
                                         if kinase_is_hek_expressed[k]],
           }, f)

   (CSV is ~3 MB and needs no extra dependency. For faster loading you can
    use df_all.to_parquet("candidate_table.parquet") instead -- requires
    `pip install pyarrow` -- and set DATA_FILE below to the .parquet name.)

------------------------------------------------------------------------------
RUN
------------------------------------------------------------------------------
       streamlit run fha_reporter_app.py

   Then open the URL it prints (usually http://localhost:8501).

------------------------------------------------------------------------------
WHAT IT DOES
------------------------------------------------------------------------------
- Pick a pathway kinase from a dropdown.
- Slide the percentile cutoff, PSSM cutoff, and max promiscuity.
- See an interactive scatter (hover for peptide details).
- See the filtered candidate table; download it as CSV.
- Toggle "HEK-expressed kinases only" for the kinase-ranking columns.
"""

import json
from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.express as px

# ------------------------------------------------------------------
# Config / data loading
# ------------------------------------------------------------------

st.set_page_config(page_title="FHA Reporter Explorer", layout="wide")

DATA_FILE = "candidate_table.csv"   # change to .parquet if you exported parquet (needs pyarrow)
META_FILE = "app_metadata.json"


@st.cache_data
def load_data():
    """Load the pre-computed candidate table and metadata (cached across reruns)."""
    if DATA_FILE.endswith(".parquet"):
        df = pd.read_parquet(DATA_FILE)
    else:
        df = pd.read_csv(DATA_FILE)
    with open(META_FILE) as f:
        meta = json.load(f)
    return df, meta


try:
    df_all, meta = load_data()
except FileNotFoundError:
    st.error(
        f"Could not find `{DATA_FILE}` or `{META_FILE}`.\n\n"
        "Export them from the notebook first (see the setup instructions at the "
        "top of this file)."
    )
    st.stop()

OUTPUT_TAG = meta["output_tag"]
PATHWAY_KINASES = meta["pathway_kinases"]
HEK_EXPRESSED = set(meta["hek_expressed_kinases"])

# ------------------------------------------------------------------
# Sidebar controls
# ------------------------------------------------------------------

st.sidebar.title(f"{OUTPUT_TAG} Reporter Explorer")
st.sidebar.markdown(
    "Explore phospho-Thr peptides predicted to bind the "
    f"**{OUTPUT_TAG}** domain and be phosphorylated by pathway kinases."
)

selected_kinase = st.sidebar.selectbox(
    "Pathway kinase", PATHWAY_KINASES,
    help="Peptides are ranked by this kinase's percentile.",
)

min_percentile = st.sidebar.slider(
    "Minimum percentile (selected kinase)", 0.0, 100.0, 90.0, 1.0,
    help="Keep peptides where the selected kinase scores at or above this percentile.",
)

min_pssm = st.sidebar.slider(
    "Minimum binding PSSM score", 0.0, 4.0, 1.0, 0.1,
    help="Keep peptides whose domain binding PSSM is at or above this value.",
)

max_promiscuity = st.sidebar.slider(
    "Maximum promiscuity index", 0, int(df_all["promiscuity_index"].max()),
    int(df_all["promiscuity_index"].max()), 1,
    help="Lower = more selective sites (fewer kinases hit them). "
         "Drag down to keep only clean reporter substrates.",
)

color_axis = st.sidebar.selectbox(
    "Color points by", ["promiscuity_index", "binding_PSSM_score",
                         f"{selected_kinase}_percentile"],
)

# ------------------------------------------------------------------
# Filter the data
# ------------------------------------------------------------------

pct_col = f"{selected_kinase}_percentile"
rank_col = f"{selected_kinase}_rank"

work = df_all.dropna(subset=[pct_col]).copy()
work = work[
    (work[pct_col] >= min_percentile)
    & (work["binding_PSSM_score"] >= min_pssm)
    & (work["promiscuity_index"] <= max_promiscuity)
]

st.title(f"{selected_kinase} candidates")
st.markdown(
    f"**{len(work)}** peptides match: {selected_kinase} percentile >= {min_percentile}, "
    f"PSSM >= {min_pssm}, promiscuity <= {max_promiscuity}."
)

if len(work) == 0:
    st.warning("No peptides match these filters. Loosen the sliders.")
    st.stop()

# ------------------------------------------------------------------
# Interactive scatter
# ------------------------------------------------------------------

hover_cols = ["Gene", "Phosphosite", "SITE_+/-7_AA", pct_col, "binding_PSSM_score",
              "promiscuity_index"]
if rank_col in work.columns:
    hover_cols.append(rank_col)

fig = px.scatter(
    work, x="binding_PSSM_score", y=pct_col, color=color_axis,
    color_continuous_scale="spring_r" if color_axis == "promiscuity_index" else "spring",
    hover_data={c: True for c in hover_cols},
    labels={"binding_PSSM_score": "Binding PSSM score", pct_col: f"{selected_kinase} percentile"},
    height=550,
)
fig.update_traces(marker=dict(size=11, line=dict(width=0.5, color="black")))
fig.update_layout(template="plotly_white")
st.plotly_chart(fig, use_container_width=True)

# ------------------------------------------------------------------
# Results table + download
# ------------------------------------------------------------------

st.subheader("Filtered peptides")

# Show pathway kinase percentile + rank columns by default
display_cols = ["Gene", "Phosphosite", "phos_res", "SITE_+/-7_AA", "Protein",
                "binding_PSSM_score", "promiscuity_index", "median_percentile"]
for pk in PATHWAY_KINASES:
    for suffix in ("_percentile", "_rank"):
        col = f"{pk}{suffix}"
        if col in work.columns:
            display_cols.append(col)

show_hek_only = st.checkbox(
    "Show only HEK293-expressed pathway kinases in the table", value=False)
if show_hek_only:
    display_cols = [c for c in display_cols
                    if not any(c.startswith(f"{pk}_") for pk in PATHWAY_KINASES
                               if pk not in HEK_EXPRESSED)]

table = work[[c for c in display_cols if c in work.columns]].sort_values(
    pct_col, ascending=False).reset_index(drop=True)
st.dataframe(table, use_container_width=True, height=400)

st.download_button(
    "Download filtered table (CSV)",
    table.to_csv(index=False).encode("utf-8"),
    file_name=f"{OUTPUT_TAG}_{selected_kinase}_filtered.csv",
    mime="text/csv",
)

# ------------------------------------------------------------------
# Top-10 quick view
# ------------------------------------------------------------------

st.subheader(f"Top 10 peptides by {selected_kinase} percentile")
top10 = work.sort_values(pct_col, ascending=False).head(10)
st.dataframe(
    top10[["Gene", "Phosphosite", "SITE_+/-7_AA", "binding_PSSM_score",
           pct_col, "promiscuity_index"]].reset_index(drop=True),
    use_container_width=True,
)
