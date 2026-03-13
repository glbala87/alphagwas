#!/usr/bin/env python3
"""
AlphaGWAS Streamlit Dashboard

Interactive web application for:
- Uploading and configuring GWAS studies
- Running the prioritization pipeline
- Exploring and visualizing results
- Downloading outputs

Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import json
import tempfile
import time
from typing import Optional, Dict, Any
import io

# Page configuration
st.set_page_config(
    page_title="AlphaGWAS",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Try to import visualization libraries
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 1rem;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        text-align: center;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 2rem;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 0.5rem 1rem;
    }
</style>
""", unsafe_allow_html=True)


def load_sample_data() -> pd.DataFrame:
    """Generate sample GWAS data for demonstration."""
    np.random.seed(42)
    n = 100

    data = {
        'chromosome': np.random.choice([str(i) for i in range(1, 23)], n),
        'position': np.random.randint(1000000, 200000000, n),
        'rsid': [f'rs{np.random.randint(1000, 9999999)}' for _ in range(n)],
        'effect_allele': np.random.choice(['A', 'C', 'G', 'T'], n),
        'other_allele': np.random.choice(['A', 'C', 'G', 'T'], n),
        'beta': np.random.normal(0, 0.5, n),
        'se': np.abs(np.random.normal(0.1, 0.05, n)),
        'pvalue': np.concatenate([
            np.random.uniform(1e-15, 1e-8, 20),
            np.random.uniform(1e-5, 1, n - 20)
        ])
    }

    return pd.DataFrame(data)


def load_sample_results() -> tuple:
    """Generate sample results for demonstration."""
    np.random.seed(42)

    # Ranked variants
    n_variants = 50
    ranked = pd.DataFrame({
        'final_rank': range(1, n_variants + 1),
        'variant_id': [f'chr{np.random.randint(1,23)}:{np.random.randint(1e6, 2e8)}' for _ in range(n_variants)],
        'rsid': [f'rs{np.random.randint(1000, 9999999)}' for _ in range(n_variants)],
        'chromosome': [str(np.random.randint(1, 23)) for _ in range(n_variants)],
        'position': [np.random.randint(1000000, 200000000) for _ in range(n_variants)],
        'consensus_score': np.sort(np.random.uniform(0.1, 1.0, n_variants))[::-1],
        'max_effect': np.sort(np.random.uniform(0.05, 0.8, n_variants))[::-1],
        'n_significant_tissues': np.random.randint(0, 10, n_variants),
        'top_tissues': ['Heart_Left_Ventricle,Whole_Blood'] * n_variants,
        'nearby_genes': ['GENE1,GENE2'] * n_variants
    })

    # Tissue scores
    tissues = ['Whole_Blood', 'Heart_Left_Ventricle', 'Liver', 'Brain_Cortex', 'Adipose_Subcutaneous']
    tissue_data = []
    for var_id in ranked['variant_id'].head(20):
        for tissue in tissues:
            tissue_data.append({
                'variant_id': var_id,
                'tissue': tissue,
                'tissue_score': np.random.uniform(0, 1),
                'max_effect': np.random.uniform(0, 0.5),
                'mean_confidence': np.random.uniform(0.5, 1)
            })
    tissue_scores = pd.DataFrame(tissue_data)

    # Predictions
    predictions = pd.DataFrame({
        'variant_id': np.repeat(ranked['variant_id'].head(20).values, 5),
        'tissue': tissues * 20,
        'effect_size': np.random.uniform(0, 0.5, 100),
        'confidence': np.random.uniform(0.5, 1, 100),
        'prediction': np.random.normal(0, 0.3, 100)
    })

    return ranked, tissue_scores, predictions


def create_manhattan_plot(df: pd.DataFrame, score_col: str = 'consensus_score') -> go.Figure:
    """Create an interactive Manhattan plot."""
    if not PLOTLY_AVAILABLE:
        st.warning("Plotly not available for interactive plots")
        return None

    # Prepare data
    df = df.copy()
    if 'chromosome' not in df.columns and 'variant_id' in df.columns:
        df[['chromosome', 'position']] = df['variant_id'].str.extract(r'chr(\w+):(\d+)')
        df['position'] = pd.to_numeric(df['position'])

    chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    df['chrom_num'] = df['chromosome'].astype(str).apply(
        lambda x: chrom_order.index(x) if x in chrom_order else 23
    )
    df = df.sort_values(['chrom_num', 'position'])

    # Calculate cumulative position
    cumpos = []
    offset = 0
    prev_chrom = None
    chrom_centers = {}

    for _, row in df.iterrows():
        if prev_chrom is not None and row['chrom_num'] != prev_chrom:
            offset += 5e7
        cumpos.append(row['position'] + offset)
        prev_chrom = row['chrom_num']

    df['cumpos'] = cumpos

    for chrom in df['chrom_num'].unique():
        chrom_centers[chrom] = df[df['chrom_num'] == chrom]['cumpos'].median()

    # Create plot
    fig = px.scatter(
        df,
        x='cumpos',
        y=score_col,
        color='chrom_num',
        hover_data=['rsid', 'chromosome', 'position', score_col],
        color_continuous_scale='Viridis'
    )

    fig.update_layout(
        title='Variant Prioritization Manhattan Plot',
        xaxis_title='Chromosome',
        yaxis_title='Prioritization Score',
        showlegend=False,
        xaxis=dict(
            tickmode='array',
            tickvals=[chrom_centers[c] for c in sorted(chrom_centers.keys())],
            ticktext=[chrom_order[int(c)] if c < len(chrom_order) else '?' for c in sorted(chrom_centers.keys())]
        ),
        height=500
    )

    return fig


def create_tissue_heatmap(tissue_scores: pd.DataFrame) -> go.Figure:
    """Create a tissue score heatmap."""
    if not PLOTLY_AVAILABLE:
        return None

    score_col = 'tissue_score' if 'tissue_score' in tissue_scores.columns else 'max_effect'

    # Pivot data
    pivot = tissue_scores.pivot_table(
        index='variant_id',
        columns='tissue',
        values=score_col,
        aggfunc='max'
    ).fillna(0)

    # Limit to top variants
    top_variants = pivot.mean(axis=1).nlargest(20).index
    pivot = pivot.loc[top_variants]

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale='YlOrRd',
        hovertemplate='Variant: %{y}<br>Tissue: %{x}<br>Score: %{z:.3f}<extra></extra>'
    ))

    fig.update_layout(
        title='Tissue-Specific Effect Scores',
        xaxis_title='Tissue',
        yaxis_title='Variant',
        height=600
    )

    return fig


def create_effect_distribution(predictions: pd.DataFrame) -> go.Figure:
    """Create effect size distribution plot."""
    if not PLOTLY_AVAILABLE:
        return None

    effect_col = 'effect_size' if 'effect_size' in predictions.columns else 'max_effect'

    fig = make_subplots(rows=1, cols=2, subplot_titles=('Effect Size Distribution', 'By Tissue'))

    # Histogram
    fig.add_trace(
        go.Histogram(x=predictions[effect_col], nbinsx=30, name='Effect Size'),
        row=1, col=1
    )

    # Box plot by tissue
    if 'tissue' in predictions.columns:
        for tissue in predictions['tissue'].unique():
            fig.add_trace(
                go.Box(y=predictions[predictions['tissue'] == tissue][effect_col], name=tissue),
                row=1, col=2
            )

    fig.update_layout(height=400, showlegend=False)
    return fig


# =============================================================================
# Main App
# =============================================================================

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 AlphaGWAS Dashboard</h1>', unsafe_allow_html=True)
    st.markdown("**Variant Prioritization using AlphaGenome**")

    # Sidebar
    with st.sidebar:
        st.header("Navigation")
        page = st.radio(
            "Select Page",
            ["🏠 Home", "📤 Upload Data", "⚙️ Configure", "🔬 Run Pipeline", "📊 Results", "📥 Download"],
            label_visibility="collapsed"
        )

        st.divider()
        st.header("Quick Actions")

        if st.button("Load Demo Data", use_container_width=True):
            st.session_state['demo_mode'] = True
            st.session_state['gwas_data'] = load_sample_data()
            ranked, tissue_scores, predictions = load_sample_results()
            st.session_state['ranked_variants'] = ranked
            st.session_state['tissue_scores'] = tissue_scores
            st.session_state['predictions'] = predictions
            st.success("Demo data loaded!")

        if st.button("Clear All Data", use_container_width=True):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()

    # ==========================================================================
    # Home Page
    # ==========================================================================
    if page == "🏠 Home":
        col1, col2 = st.columns([2, 1])

        with col1:
            st.markdown("""
            ### Welcome to AlphaGWAS

            AlphaGWAS is a pipeline for prioritizing genetic variants at GWAS loci
            using Google DeepMind's **AlphaGenome** AI model.

            #### Pipeline Steps:
            1. **Upload** your GWAS summary statistics
            2. **Configure** tissues, modalities, and parameters
            3. **Run** the prioritization pipeline
            4. **Explore** interactive visualizations
            5. **Download** ranked variants and reports

            #### Features:
            - 🚀 Parallel processing for fast predictions
            - 📊 Interactive Plotly visualizations
            - 🧬 ClinVar, gnomAD, GTEx annotations
            - 📄 Downloadable HTML reports
            """)

        with col2:
            st.markdown("#### Quick Stats")

            if 'gwas_data' in st.session_state:
                data = st.session_state['gwas_data']
                st.metric("Variants Loaded", f"{len(data):,}")

                if 'pvalue' in data.columns:
                    sig = (data['pvalue'] < 5e-8).sum()
                    st.metric("Significant (p<5e-8)", f"{sig:,}")

            if 'ranked_variants' in st.session_state:
                ranked = st.session_state['ranked_variants']
                st.metric("Variants Ranked", f"{len(ranked):,}")

            if 'demo_mode' in st.session_state:
                st.info("Demo mode active")

    # ==========================================================================
    # Upload Data Page
    # ==========================================================================
    elif page == "📤 Upload Data":
        st.header("Upload GWAS Data")

        upload_tab, paste_tab = st.tabs(["📁 Upload File", "📋 Paste Data"])

        with upload_tab:
            uploaded_file = st.file_uploader(
                "Upload GWAS summary statistics",
                type=['tsv', 'csv', 'txt', 'gz'],
                help="Supported formats: TSV, CSV, TXT (tab or comma separated)"
            )

            if uploaded_file:
                try:
                    # Detect separator
                    sep = '\t' if uploaded_file.name.endswith('.tsv') else ','

                    if uploaded_file.name.endswith('.gz'):
                        df = pd.read_csv(uploaded_file, sep=sep, compression='gzip')
                    else:
                        df = pd.read_csv(uploaded_file, sep=sep)

                    st.success(f"Loaded {len(df):,} variants")
                    st.session_state['gwas_data'] = df

                    # Preview
                    st.subheader("Data Preview")
                    st.dataframe(df.head(10), use_container_width=True)

                    # Column info
                    st.subheader("Columns Detected")
                    st.write(list(df.columns))

                except Exception as e:
                    st.error(f"Error loading file: {e}")

        with paste_tab:
            pasted_data = st.text_area(
                "Paste tab-separated data",
                height=200,
                placeholder="CHR\tBP\tSNP\tA1\tA2\tP\n1\t100000\trs123\tA\tG\t1e-10"
            )

            if st.button("Parse Pasted Data"):
                if pasted_data:
                    try:
                        df = pd.read_csv(io.StringIO(pasted_data), sep='\t')
                        st.session_state['gwas_data'] = df
                        st.success(f"Parsed {len(df):,} variants")
                        st.dataframe(df.head(), use_container_width=True)
                    except Exception as e:
                        st.error(f"Error parsing data: {e}")

    # ==========================================================================
    # Configure Page
    # ==========================================================================
    elif page == "⚙️ Configure":
        st.header("Pipeline Configuration")

        if 'gwas_data' not in st.session_state:
            st.warning("Please upload GWAS data first")
            st.stop()

        df = st.session_state['gwas_data']

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Column Mapping")

            config = {}
            config['chrom_col'] = st.selectbox("Chromosome column", df.columns, index=0)
            config['pos_col'] = st.selectbox("Position column", df.columns, index=min(1, len(df.columns)-1))
            config['rsid_col'] = st.selectbox("rsID column", df.columns, index=min(2, len(df.columns)-1))
            config['pval_col'] = st.selectbox("P-value column", df.columns, index=min(5, len(df.columns)-1))

        with col2:
            st.subheader("Analysis Settings")

            config['pval_threshold'] = st.number_input(
                "P-value threshold",
                value=5e-8,
                format="%.0e"
            )

            config['tissues'] = st.multiselect(
                "Tissues",
                ['Whole_Blood', 'Heart_Left_Ventricle', 'Heart_Atrial_Appendage',
                 'Liver', 'Brain_Cortex', 'Adipose_Subcutaneous', 'Muscle_Skeletal'],
                default=['Whole_Blood', 'Heart_Left_Ventricle']
            )

            config['modalities'] = st.multiselect(
                "Modalities",
                ['expression', 'chromatin_accessibility', 'histone_marks', 'splicing'],
                default=['expression']
            )

        st.session_state['config'] = config

        if st.button("Save Configuration", type="primary"):
            st.success("Configuration saved!")

    # ==========================================================================
    # Run Pipeline Page
    # ==========================================================================
    elif page == "🔬 Run Pipeline":
        st.header("Run Pipeline")

        if 'gwas_data' not in st.session_state:
            st.warning("Please upload GWAS data first")
            st.stop()

        st.info("""
        **Note:** In this demo, the pipeline simulation runs locally.
        For actual AlphaGenome predictions, you need:
        1. AlphaGenome API key
        2. Run via command line: `python run_pipeline.py`
        """)

        if st.button("Run Demo Pipeline", type="primary"):
            progress_bar = st.progress(0)
            status = st.empty()

            steps = [
                ("Extracting significant variants...", 0.2),
                ("Finding LD proxies...", 0.4),
                ("Running liftover...", 0.5),
                ("Generating mock predictions...", 0.8),
                ("Scoring variants...", 0.9),
                ("Complete!", 1.0)
            ]

            for msg, progress in steps:
                status.text(msg)
                progress_bar.progress(progress)
                time.sleep(0.5)

            # Load demo results
            ranked, tissue_scores, predictions = load_sample_results()
            st.session_state['ranked_variants'] = ranked
            st.session_state['tissue_scores'] = tissue_scores
            st.session_state['predictions'] = predictions

            st.success("Pipeline completed! Go to Results page to explore.")

    # ==========================================================================
    # Results Page
    # ==========================================================================
    elif page == "📊 Results":
        st.header("Results Explorer")

        if 'ranked_variants' not in st.session_state:
            st.warning("No results available. Run the pipeline or load demo data first.")
            st.stop()

        ranked = st.session_state['ranked_variants']
        tissue_scores = st.session_state.get('tissue_scores', pd.DataFrame())
        predictions = st.session_state.get('predictions', pd.DataFrame())

        # Metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Variants", len(ranked))
        with col2:
            st.metric("Top Score", f"{ranked['consensus_score'].max():.3f}")
        with col3:
            if 'n_significant_tissues' in ranked.columns:
                st.metric("Avg Tissues", f"{ranked['n_significant_tissues'].mean():.1f}")
        with col4:
            if not tissue_scores.empty:
                st.metric("Tissues Analyzed", tissue_scores['tissue'].nunique())

        # Tabs for different views
        tab1, tab2, tab3, tab4 = st.tabs(["📋 Table", "🗺️ Manhattan", "🔥 Heatmap", "📈 Distributions"])

        with tab1:
            st.subheader("Ranked Variants")

            # Filters
            col1, col2 = st.columns(2)
            with col1:
                min_score = st.slider("Min Score", 0.0, 1.0, 0.0)
            with col2:
                max_rows = st.slider("Max Rows", 10, 500, 50)

            filtered = ranked[ranked['consensus_score'] >= min_score].head(max_rows)
            st.dataframe(filtered, use_container_width=True, height=400)

        with tab2:
            st.subheader("Manhattan Plot")
            if PLOTLY_AVAILABLE:
                fig = create_manhattan_plot(ranked)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Install plotly for interactive plots: `pip install plotly`")

        with tab3:
            st.subheader("Tissue Heatmap")
            if not tissue_scores.empty and PLOTLY_AVAILABLE:
                fig = create_tissue_heatmap(tissue_scores)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Tissue scores not available")

        with tab4:
            st.subheader("Effect Distributions")
            if not predictions.empty and PLOTLY_AVAILABLE:
                fig = create_effect_distribution(predictions)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Predictions not available")

    # ==========================================================================
    # Download Page
    # ==========================================================================
    elif page == "📥 Download":
        st.header("Download Results")

        if 'ranked_variants' not in st.session_state:
            st.warning("No results available to download.")
            st.stop()

        ranked = st.session_state['ranked_variants']
        tissue_scores = st.session_state.get('tissue_scores', pd.DataFrame())

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Ranked Variants")
            csv = ranked.to_csv(index=False)
            st.download_button(
                "📥 Download TSV",
                csv,
                "ranked_variants.tsv",
                "text/tab-separated-values",
                use_container_width=True
            )

            st.download_button(
                "📥 Download Top 20",
                ranked.head(20).to_csv(index=False),
                "top20_variants.tsv",
                "text/tab-separated-values",
                use_container_width=True
            )

        with col2:
            st.subheader("Tissue Scores")
            if not tissue_scores.empty:
                st.download_button(
                    "📥 Download Tissue Scores",
                    tissue_scores.to_csv(index=False),
                    "tissue_scores.tsv",
                    "text/tab-separated-values",
                    use_container_width=True
                )

        st.subheader("Summary Report")
        if st.button("Generate HTML Report"):
            # Create simple HTML report
            html = f"""
            <html>
            <head><title>AlphaGWAS Report</title></head>
            <body>
            <h1>AlphaGWAS Results Summary</h1>
            <p>Total variants: {len(ranked)}</p>
            <p>Top score: {ranked['consensus_score'].max():.4f}</p>
            <h2>Top 20 Variants</h2>
            {ranked.head(20).to_html()}
            </body>
            </html>
            """
            st.download_button(
                "📥 Download HTML Report",
                html,
                "alphagwas_report.html",
                "text/html",
                use_container_width=True
            )


if __name__ == "__main__":
    main()
