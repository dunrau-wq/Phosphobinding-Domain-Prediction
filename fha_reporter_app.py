"""
Phospho-Binding-Domain Reporter Explorer - interactive web app (V5)
===================================================================

Explore candidate phosphopeptides that (1) bind a phospho-binding domain and
(2) are phosphorylated by kinases of interest. Multi-domain and fully generic:
the app auto-discovers whatever domains you've exported and works for any
kinase set.

It does NOT re-run the atlas stream - it loads pre-computed candidate tables.

------------------------------------------------------------------------------
SETUP
------------------------------------------------------------------------------
1. Dependencies (also list these in requirements.txt for Streamlit Cloud):
       pip install streamlit plotly pandas

2. From the notebook, export one pair of files PER DOMAIN using the V5 export
   cell (export_cell_v5.py). You get, e.g.:
       candidate_table_FHA1.csv   app_metadata_FHA1.json
       candidate_table_FHA2.csv   app_metadata_FHA2.json
   Put all of them in the SAME folder as this script (the repo root on
   Streamlit Cloud).

3. Run:
       streamlit run fha_reporter_app_v5.py

------------------------------------------------------------------------------
WHAT IT DOES
------------------------------------------------------------------------------
- Auto-discovers every candidate_table_<DOMAIN>.csv next to this script.
- "Single domain" view: pick one domain + one kinase; filter by percentile,
  PSSM, promiscuity; interactive scatter (cool colorscale) + table + download.
- "Overlay domains" view: plot multiple domains on one graph for the same
  kinase, colored by domain, so you can compare binding landscapes.
"""

import json
import glob
from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go

# ------------------------------------------------------------------
# Config
# ------------------------------------------------------------------

st.set_page_config(page_title="Reporter Explorer", layout="wide")

# Resolve files relative to THIS script so it works locally and on Streamlit Cloud.
SCRIPT_DIR = Path(__file__).resolve().parent

# Colorscale for single-domain continuous coloring.
# matplotlib-style "cool" (cyan -> magenta). Change to any Plotly named scale
# ("Viridis", "Plasma", "Cividis", ...) or a custom [[pos,color],...] list.
COOL_COLORSCALE = [[0.0, "#00FFFF"], [1.0, "#FF00FF"]]


def reverse_colorscale(scale):
    """
    Reverse a [[pos,color],...] colorscale correctly: positions must stay
    ascending 0->1, so we flip the colors but keep the position order.
    (A plain list[::-1] would put positions in descending order, which Plotly
    rejects.)
    """
    positions = [p for p, _ in scale]
    colors = [c for _, c in scale]
    return [[p, c] for p, c in zip(positions, colors[::-1])]

# Palette for overlaying multiple domains (categorical, visually distinct).
DOMAIN_PALETTE = px.colors.qualitative.Plotly


# ------------------------------------------------------------------
# Discover available domains
# ------------------------------------------------------------------

@st.cache_data
def discover_domains():
    """
    Find all candidate_table_<DOMAIN>.csv files next to this script and pair
    them with their metadata JSON. Returns {domain_tag: {"csv":path,"meta":path}}.
    """
    found = {}
    for csv_path in sorted(glob.glob(str(SCRIPT_DIR / "candidate_table_*.csv"))):
        tag = Path(csv_path).stem.replace("candidate_table_", "")
        meta_path = SCRIPT_DIR / f"app_metadata_{tag}.json"
        if meta_path.exists():
            found[tag] = {"csv": csv_path, "meta": str(meta_path)}
    return found


@st.cache_data
def load_domain(tag, paths):
    """Load one domain's candidate table + metadata (cached)."""
    df = pd.read_csv(paths["csv"])
    with open(paths["meta"]) as f:
        meta = json.load(f)
    return df, meta


domains = discover_domains()

if not domains:
    st.error(
        "No candidate tables found.\n\n"
        "Expected files named `candidate_table_<DOMAIN>.csv` (with matching "
        "`app_metadata_<DOMAIN>.json`) in the same folder as this script.\n\n"
        "Export them from the notebook using the V5 export cell, once per domain."
    )
    st.stop()

# ------------------------------------------------------------------
# Sidebar: view mode + shared controls
# ------------------------------------------------------------------

st.sidebar.title("Reporter Explorer")
st.sidebar.caption(f"Discovered domains: {', '.join(domains.keys())}")

view_mode = st.sidebar.radio(
    "View mode",
    ["Single domain", "Overlay domains"],
    help="Single: explore one domain in depth. Overlay: compare domains on one graph.",
)

# Shared filter sliders
min_percentile = st.sidebar.slider(
    "Minimum percentile (selected kinase)", 0.0, 100.0, 90.0, 1.0,
    help="Keep peptides where the selected kinase scores at or above this percentile.",
)
min_pssm = st.sidebar.slider(
    "Minimum binding PSSM score", 0.0, 4.0, 1.0, 0.1,
    help="Keep peptides whose domain binding PSSM is at or above this value.",
)


def filter_frame(df, kinase, min_pct, min_ps, max_prom):
    """Apply the standard percentile / PSSM / promiscuity filters for one kinase."""
    pct_col = f"{kinase}_percentile"
    if pct_col not in df.columns:
        return df.iloc[0:0]  # empty - kinase not in this domain's table
    out = df.dropna(subset=[pct_col]).copy()
    out = out[
        (out[pct_col] >= min_pct)
        & (out["binding_PSSM_score"] >= min_ps)
        & (out["promiscuity_index"] <= max_prom)
    ]
    return out


# ==================================================================
# SINGLE DOMAIN VIEW
# ==================================================================

if view_mode == "Single domain":
    domain_tag = st.sidebar.selectbox("Domain", list(domains.keys()))
    df_all, meta = load_domain(domain_tag, domains[domain_tag])

    pathway_kinases = meta["pathway_kinases"]
    hek_expressed = set(meta["hek_expressed_kinases"])

    selected_kinase = st.sidebar.selectbox(
        "Kinase", pathway_kinases,
        help="Peptides are ranked by this kinase's percentile.",
    )

    max_prom_val = int(df_all["promiscuity_index"].max())
    max_promiscuity = st.sidebar.slider(
        "Maximum promiscuity index", 0, max_prom_val, max_prom_val, 1,
        help="Lower = more selective sites. Drag down to keep clean reporter substrates.",
    )
    color_axis = st.sidebar.selectbox(
        "Color points by",
        ["promiscuity_index", "binding_PSSM_score", f"{selected_kinase}_percentile"],
    )

    pct_col = f"{selected_kinase}_percentile"
    rank_col = f"{selected_kinase}_rank"
    work = filter_frame(df_all, selected_kinase, min_percentile, min_pssm, max_promiscuity)

    st.title(f"{domain_tag}: {selected_kinase} candidates")
    st.caption(meta.get("domain_description", ""))
    st.markdown(
        f"**{len(work)}** peptides match: {selected_kinase} percentile >= {min_percentile}, "
        f"PSSM >= {min_pssm}, promiscuity <= {max_promiscuity}."
    )

    if len(work) == 0:
        st.warning("No peptides match these filters. Loosen the sliders.")
        st.stop()

    # --- Interactive scatter (cool colorscale) ---
    hover_cols = ["Gene", "Phosphosite", "SITE_+/-7_AA", pct_col, "binding_PSSM_score",
                  "promiscuity_index"]
    if rank_col in work.columns:
        hover_cols.append(rank_col)

    # Reverse the scale for promiscuity so low (good) reads as the bright end.
    scale = COOL_COLORSCALE if color_axis != "promiscuity_index" else reverse_colorscale(COOL_COLORSCALE)
    fig = px.scatter(
        work, x="binding_PSSM_score", y=pct_col, color=color_axis,
        color_continuous_scale=scale,
        hover_data={c: True for c in hover_cols},
        labels={"binding_PSSM_score": "Binding PSSM score",
                pct_col: f"{selected_kinase} percentile"},
        height=550,
    )
    fig.update_traces(marker=dict(size=11, line=dict(width=0.5, color="black")))
    fig.update_layout(template="plotly_white")
    st.plotly_chart(fig, use_container_width=True)

    # --- Table + download ---
    st.subheader("Filtered peptides")
    display_cols = ["Gene", "Phosphosite", "phos_res", "SITE_+/-7_AA", "Protein",
                    "binding_PSSM_score", "promiscuity_index", "median_percentile"]
    for pk in pathway_kinases:
        for suffix in ("_percentile", "_rank"):
            col = f"{pk}{suffix}"
            if col in work.columns:
                display_cols.append(col)

    show_hek_only = st.checkbox(
        "Show only HEK293-expressed pathway kinases in the table", value=False)
    if show_hek_only:
        display_cols = [c for c in display_cols
                        if not any(c.startswith(f"{pk}_") for pk in pathway_kinases
                                   if pk not in hek_expressed)]

    table = work[[c for c in display_cols if c in work.columns]].sort_values(
        pct_col, ascending=False).reset_index(drop=True)
    st.dataframe(table, use_container_width=True, height=400)
    st.download_button(
        "Download filtered table (CSV)",
        table.to_csv(index=False).encode("utf-8"),
        file_name=f"{domain_tag}_{selected_kinase}_filtered.csv",
        mime="text/csv",
    )

    # --- Top-10 quick view ---
    st.subheader(f"Top 10 peptides by {selected_kinase} percentile")
    top10 = work.sort_values(pct_col, ascending=False).head(10)
    st.dataframe(
        top10[["Gene", "Phosphosite", "SITE_+/-7_AA", "binding_PSSM_score",
               pct_col, "promiscuity_index"]].reset_index(drop=True),
        use_container_width=True,
    )


# ==================================================================
# OVERLAY DOMAINS VIEW
# ==================================================================

else:
    st.sidebar.markdown("---")
    chosen_domains = st.sidebar.multiselect(
        "Domains to overlay", list(domains.keys()), default=list(domains.keys()),
    )
    if not chosen_domains:
        st.warning("Select at least one domain in the sidebar.")
        st.stop()

    # Kinase choices = union of pathway kinases across chosen domains
    kinase_union = sorted(set().union(
        *[set(load_domain(t, domains[t])[1]["pathway_kinases"]) for t in chosen_domains]
    ))
    selected_kinase = st.sidebar.selectbox("Kinase", kinase_union)

    # Promiscuity max across chosen domains
    global_max_prom = int(max(
        load_domain(t, domains[t])[0]["promiscuity_index"].max() for t in chosen_domains
    ))
    max_promiscuity = st.sidebar.slider(
        "Maximum promiscuity index", 0, global_max_prom, global_max_prom, 1)

    pct_col = f"{selected_kinase}_percentile"

    st.title(f"Overlay: {selected_kinase} across {', '.join(chosen_domains)}")
    st.markdown(
        f"Filters: {selected_kinase} percentile >= {min_percentile}, "
        f"PSSM >= {min_pssm}, promiscuity <= {max_promiscuity}."
    )

    # Build one combined frame tagged by domain
    frames = []
    for t in chosen_domains:
        df_d, _ = load_domain(t, domains[t])
        sub = filter_frame(df_d, selected_kinase, min_percentile, min_pssm, max_promiscuity)
        if len(sub) > 0:
            sub = sub.copy()
            sub["domain"] = t
            frames.append(sub)

    if not frames:
        st.warning("No peptides match these filters in any chosen domain.")
        st.stop()

    combined = pd.concat(frames, ignore_index=True)
    st.markdown(f"**{len(combined)}** peptides total - "
                + ", ".join(f"{t}: {(combined['domain']==t).sum()}" for t in chosen_domains))

    # Scatter colored by domain (categorical palette)
    hover_cols = ["Gene", "Phosphosite", "SITE_+/-7_AA", pct_col,
                  "binding_PSSM_score", "promiscuity_index"]
    fig = px.scatter(
        combined, x="binding_PSSM_score", y=pct_col, color="domain",
        color_discrete_sequence=DOMAIN_PALETTE,
        hover_data={c: True for c in hover_cols if c in combined.columns},
        labels={"binding_PSSM_score": "Binding PSSM score",
                pct_col: f"{selected_kinase} percentile"},
        height=600,
    )
    fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color="black")))
    fig.update_layout(template="plotly_white", legend_title_text="Domain")
    st.plotly_chart(fig, use_container_width=True)

    # Combined table + download
    st.subheader("Combined filtered peptides")
    show_cols = ["domain", "Gene", "Phosphosite", "phos_res", "SITE_+/-7_AA", "Protein",
                 "binding_PSSM_score", "promiscuity_index", pct_col, f"{selected_kinase}_rank"]
    table = combined[[c for c in show_cols if c in combined.columns]].sort_values(
        ["domain", pct_col], ascending=[True, False]).reset_index(drop=True)
    st.dataframe(table, use_container_width=True, height=400)
    st.download_button(
        "Download combined table (CSV)",
        table.to_csv(index=False).encode("utf-8"),
        file_name=f"overlay_{selected_kinase}_{'_'.join(chosen_domains)}.csv",
        mime="text/csv",
    )
