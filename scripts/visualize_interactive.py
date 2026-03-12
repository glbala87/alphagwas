#!/usr/bin/env python3
"""
Interactive visualization module for AlphaGWAS pipeline.

Creates interactive Plotly visualizations:
- Interactive Manhattan plots with hover info
- Zoomable tissue heatmaps
- Dynamic effect size distributions
- Interactive HTML dashboards
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Optional, Any
import yaml
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Try to import Plotly
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not installed. Install with: pip install plotly")


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def check_plotly():
    """Check if Plotly is available."""
    if not PLOTLY_AVAILABLE:
        raise ImportError(
            "Plotly is required for interactive visualizations. "
            "Install with: pip install plotly"
        )


def create_interactive_manhattan(
    variants: pd.DataFrame,
    score_col: str = 'consensus_score',
    output_path: Optional[str] = None,
    title: str = "Interactive Variant Prioritization",
    highlight_top_n: int = 20
) -> Optional[go.Figure]:
    """
    Create an interactive Manhattan plot with Plotly.

    Features:
    - Hover to see variant details
    - Click to select variants
    - Zoom and pan
    - Highlight top variants

    Args:
        variants: DataFrame with chromosome, position, and score columns
        score_col: Column name for the score to plot
        output_path: Path to save the HTML file
        title: Plot title
        highlight_top_n: Number of top variants to highlight

    Returns:
        Plotly Figure object
    """
    check_plotly()

    df = variants.copy()

    # Handle variant_id if chromosome/position not present
    if 'chromosome' not in df.columns and 'variant_id' in df.columns:
        df[['chromosome', 'position']] = df['variant_id'].str.extract(r'chr(\w+):(\d+)')
        df['position'] = pd.to_numeric(df['position'])

    # Convert chromosome to numeric for ordering
    chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    df['chrom_num'] = df['chromosome'].astype(str).apply(
        lambda x: chrom_order.index(x) if x in chrom_order else 23
    )

    # Sort and calculate cumulative position
    df = df.sort_values(['chrom_num', 'position'])

    cumulative_pos = []
    offset = 0
    chrom_centers = {}
    prev_chrom = None

    for idx, row in df.iterrows():
        chrom = row['chrom_num']
        if prev_chrom is not None and chrom != prev_chrom:
            offset += df[df['chrom_num'] == prev_chrom]['position'].max() + 1e7
        cumulative_pos.append(row['position'] + offset)
        prev_chrom = chrom

    df['cumulative_pos'] = cumulative_pos

    # Calculate chromosome centers for x-axis labels
    for chrom in df['chrom_num'].unique():
        chrom_data = df[df['chrom_num'] == chrom]
        chrom_centers[chrom] = chrom_data['cumulative_pos'].median()

    # Identify top variants
    top_mask = df[score_col].rank(ascending=False) <= highlight_top_n

    # Create hover text
    df['hover_text'] = df.apply(
        lambda row: (
            f"<b>{row.get('rsid', row.get('variant_id', 'Unknown'))}</b><br>"
            f"Chr{row['chromosome']}:{row['position']:,}<br>"
            f"Score: {row[score_col]:.4f}<br>"
            f"Rank: {int(row.get('final_rank', 0))}" +
            (f"<br>Genes: {row.get('nearby_genes', 'N/A')}" if 'nearby_genes' in row else "")
        ),
        axis=1
    )

    # Create figure
    fig = go.Figure()

    # Color palette for chromosomes
    colors = ['#1f77b4', '#aec7e8'] * 12

    # Plot each chromosome
    for i, chrom in enumerate(sorted(df['chrom_num'].unique())):
        chrom_data = df[(df['chrom_num'] == chrom) & ~top_mask]

        fig.add_trace(go.Scatter(
            x=chrom_data['cumulative_pos'],
            y=chrom_data[score_col],
            mode='markers',
            marker=dict(
                size=6,
                color=colors[i % 2],
                opacity=0.6
            ),
            text=chrom_data['hover_text'],
            hoverinfo='text',
            name=f'Chr {chrom_order[int(chrom)] if chrom < len(chrom_order) else "?"}',
            showlegend=False
        ))

    # Add top variants with different styling
    top_variants = df[top_mask]
    fig.add_trace(go.Scatter(
        x=top_variants['cumulative_pos'],
        y=top_variants[score_col],
        mode='markers+text',
        marker=dict(
            size=12,
            color='red',
            symbol='star',
            line=dict(color='darkred', width=1)
        ),
        text=top_variants.apply(
            lambda row: row.get('rsid', '')[:10] if pd.notna(row.get('rsid')) else '',
            axis=1
        ),
        textposition='top center',
        textfont=dict(size=9),
        hovertext=top_variants['hover_text'],
        hoverinfo='text',
        name=f'Top {highlight_top_n} Variants'
    ))

    # Update layout
    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis=dict(
            title='Chromosome',
            tickmode='array',
            tickvals=[chrom_centers[c] for c in sorted(chrom_centers.keys())],
            ticktext=[chrom_order[int(c)] if c < len(chrom_order) else '?' for c in sorted(chrom_centers.keys())],
            showgrid=False
        ),
        yaxis=dict(
            title=f'Prioritization Score ({score_col})',
            showgrid=True,
            gridcolor='lightgray'
        ),
        hovermode='closest',
        plot_bgcolor='white',
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        ),
        height=600
    )

    # Add threshold line
    if len(df) > highlight_top_n:
        threshold = df.nlargest(highlight_top_n, score_col)[score_col].min()
        fig.add_hline(
            y=threshold,
            line_dash="dash",
            line_color="red",
            opacity=0.5,
            annotation_text=f"Top {highlight_top_n} threshold"
        )

    if output_path:
        fig.write_html(output_path)
        logger.info(f"Saved interactive Manhattan plot to {output_path}")

    return fig


def create_interactive_heatmap(
    tissue_scores: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Tissue-Specific Effects",
    top_n_variants: int = 30
) -> Optional[go.Figure]:
    """
    Create an interactive heatmap of variant scores across tissues.

    Args:
        tissue_scores: DataFrame with variant_id, tissue, and score columns
        output_path: Path to save the HTML file
        title: Plot title
        top_n_variants: Number of top variants to show

    Returns:
        Plotly Figure object
    """
    check_plotly()

    df = tissue_scores.copy()
    score_col = 'tissue_score' if 'tissue_score' in df.columns else 'max_effect'

    # Get top variants by mean score
    variant_means = df.groupby('variant_id')[score_col].mean().nlargest(top_n_variants)
    top_variants = variant_means.index.tolist()

    # Filter and pivot
    df_top = df[df['variant_id'].isin(top_variants)]
    pivot = df_top.pivot_table(
        index='variant_id',
        columns='tissue',
        values=score_col,
        aggfunc='max'
    ).fillna(0)

    # Reorder by mean score
    pivot = pivot.loc[top_variants]

    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(),
        y=pivot.index.tolist(),
        colorscale='YlOrRd',
        hovertemplate=(
            "Variant: %{y}<br>"
            "Tissue: %{x}<br>"
            "Score: %{z:.4f}<extra></extra>"
        ),
        colorbar=dict(title="Effect Score")
    ))

    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis=dict(
            title="Tissue",
            tickangle=45
        ),
        yaxis=dict(
            title="Variant",
            autorange="reversed"
        ),
        height=max(500, top_n_variants * 20),
        margin=dict(l=150, b=150)
    )

    if output_path:
        fig.write_html(output_path)
        logger.info(f"Saved interactive heatmap to {output_path}")

    return fig


def create_interactive_distribution(
    predictions: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Effect Size Distribution"
) -> Optional[go.Figure]:
    """
    Create interactive distribution plots.

    Args:
        predictions: DataFrame with effect_size column
        output_path: Path to save the HTML file
        title: Plot title

    Returns:
        Plotly Figure object
    """
    check_plotly()

    df = predictions.copy()
    effect_col = 'effect_size' if 'effect_size' in df.columns else 'max_effect'

    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Effect Size Distribution", "Effect by Tissue"),
        specs=[[{"type": "histogram"}, {"type": "box"}]]
    )

    # Histogram
    fig.add_trace(
        go.Histogram(
            x=df[effect_col],
            nbinsx=50,
            name="Effect Size",
            marker_color='steelblue',
            opacity=0.7
        ),
        row=1, col=1
    )

    # Add mean line
    mean_val = df[effect_col].mean()
    median_val = df[effect_col].median()

    fig.add_vline(
        x=mean_val, line_dash="dash", line_color="red",
        annotation_text=f"Mean: {mean_val:.3f}",
        row=1, col=1
    )

    # Box plots by tissue
    if 'tissue' in df.columns:
        for tissue in df['tissue'].unique():
            tissue_data = df[df['tissue'] == tissue][effect_col]
            fig.add_trace(
                go.Box(
                    y=tissue_data,
                    name=tissue,
                    boxmean=True
                ),
                row=1, col=2
            )

    fig.update_layout(
        title=dict(text=title, x=0.5),
        showlegend=False,
        height=500
    )

    if output_path:
        fig.write_html(output_path)
        logger.info(f"Saved interactive distribution to {output_path}")

    return fig


def create_interactive_scatter(
    predictions: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Confidence vs Effect Size"
) -> Optional[go.Figure]:
    """
    Create an interactive scatter plot of confidence vs effect size.

    Args:
        predictions: DataFrame with confidence and effect_size columns
        output_path: Path to save the HTML file
        title: Plot title

    Returns:
        Plotly Figure object
    """
    check_plotly()

    df = predictions.copy()

    if 'confidence' not in df.columns or 'effect_size' not in df.columns:
        logger.warning("Required columns not found for scatter plot")
        return None

    # Create figure with color by tissue if available
    if 'tissue' in df.columns:
        fig = px.scatter(
            df,
            x='effect_size',
            y='confidence',
            color='tissue',
            hover_data=['variant_id'] if 'variant_id' in df.columns else None,
            opacity=0.6,
            title=title
        )
    else:
        fig = px.scatter(
            df,
            x='effect_size',
            y='confidence',
            opacity=0.6,
            title=title
        )

    # Add quadrant lines
    fig.add_hline(y=df['confidence'].median(), line_dash="dash", line_color="gray", opacity=0.5)
    fig.add_vline(x=df['effect_size'].median(), line_dash="dash", line_color="gray", opacity=0.5)

    fig.update_layout(
        xaxis_title="Effect Size",
        yaxis_title="Confidence",
        height=600
    )

    if output_path:
        fig.write_html(output_path)
        logger.info(f"Saved interactive scatter plot to {output_path}")

    return fig


def create_interactive_dashboard(
    ranked_variants: pd.DataFrame,
    tissue_scores: pd.DataFrame,
    predictions: pd.DataFrame,
    output_path: str,
    study_name: str = "",
    phenotype: str = ""
) -> str:
    """
    Create a comprehensive interactive HTML dashboard.

    Args:
        ranked_variants: DataFrame with ranked variants
        tissue_scores: DataFrame with tissue scores
        predictions: DataFrame with predictions
        output_path: Path to save the dashboard HTML
        study_name: Study name for title
        phenotype: Phenotype being studied

    Returns:
        Path to generated HTML file
    """
    check_plotly()

    from datetime import datetime

    # Create all figures
    manhattan_fig = create_interactive_manhattan(
        ranked_variants,
        title=f"Variant Prioritization - {phenotype}" if phenotype else "Variant Prioritization"
    )

    heatmap_fig = create_interactive_heatmap(
        tissue_scores,
        title="Tissue-Specific Effects"
    )

    dist_fig = create_interactive_distribution(
        predictions,
        title="Effect Size Distribution"
    )

    scatter_fig = create_interactive_scatter(
        predictions,
        title="Confidence vs Effect"
    )

    # Get top variants table
    top_variants = ranked_variants.head(20)
    display_cols = ['final_rank', 'rsid', 'variant_id', 'consensus_score', 'max_effect', 'top_tissues']
    display_cols = [c for c in display_cols if c in top_variants.columns]

    # Create summary stats
    n_variants = len(ranked_variants)
    n_tissues = tissue_scores['tissue'].nunique() if 'tissue' in tissue_scores.columns else 0
    mean_score = ranked_variants['consensus_score'].mean() if 'consensus_score' in ranked_variants.columns else 0

    # Build HTML
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AlphaGWAS Dashboard - {study_name}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 20px;
            text-align: center;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .stat-value {{ font-size: 2em; font-weight: bold; color: #667eea; }}
        .stat-label {{ color: #666; margin-top: 5px; }}
        .plot-container {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #eee; }}
        th {{ background: #f8f9fa; font-weight: 600; }}
        tr:hover {{ background: #f8f9fa; }}
        .footer {{ text-align: center; color: #666; margin-top: 30px; padding: 20px; }}
        h2 {{ color: #333; border-bottom: 2px solid #667eea; padding-bottom: 10px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>AlphaGWAS Interactive Dashboard</h1>
        <p>Study: {study_name or 'Unnamed'} | Phenotype: {phenotype or 'Unknown'}</p>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>

    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-value">{n_variants:,}</div>
            <div class="stat-label">Total Variants</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{n_tissues}</div>
            <div class="stat-label">Tissues Analyzed</div>
        </div>
        <div class="stat-card">
            <div class="stat-value">{mean_score:.3f}</div>
            <div class="stat-label">Mean Score</div>
        </div>
    </div>

    <div class="plot-container">
        <h2>Manhattan Plot</h2>
        <div id="manhattan"></div>
    </div>

    <div class="plot-grid">
        <div class="plot-container">
            <h2>Tissue Heatmap</h2>
            <div id="heatmap"></div>
        </div>
        <div class="plot-container">
            <h2>Effect Distribution</h2>
            <div id="distribution"></div>
        </div>
    </div>

    <div class="plot-container">
        <h2>Confidence vs Effect</h2>
        <div id="scatter"></div>
    </div>

    <div class="plot-container">
        <h2>Top 20 Prioritized Variants</h2>
        <table>
            <tr>{''.join(f'<th>{c}</th>' for c in display_cols)}</tr>
"""

    for _, row in top_variants.iterrows():
        html += "<tr>"
        for col in display_cols:
            val = row.get(col, '')
            if isinstance(val, float):
                val = f"{val:.4f}"
            html += f"<td>{val}</td>"
        html += "</tr>\n"

    html += f"""
        </table>
    </div>

    <div class="footer">
        Generated by AlphaGWAS | Powered by AlphaGenome & Plotly
    </div>

    <script>
        var manhattan = {manhattan_fig.to_json()};
        var heatmap = {heatmap_fig.to_json()};
        var distribution = {dist_fig.to_json()};
        var scatter = {scatter_fig.to_json() if scatter_fig else '{}'};

        Plotly.newPlot('manhattan', manhattan.data, manhattan.layout, {{responsive: true}});
        Plotly.newPlot('heatmap', heatmap.data, heatmap.layout, {{responsive: true}});
        Plotly.newPlot('distribution', distribution.data, distribution.layout, {{responsive: true}});
        if (scatter.data) {{
            Plotly.newPlot('scatter', scatter.data, scatter.layout, {{responsive: true}});
        }}
    </script>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)

    logger.info(f"Generated interactive dashboard: {output_path}")
    return output_path


def main(config_path: str = "config/config.yaml"):
    """Generate all interactive visualizations for a completed pipeline run."""
    config = load_config(config_path)

    prefix = config['output']['prefix']
    output_dir = Path(config['output']['dir'])
    intermediate_dir = Path("data/intermediate")

    study_name = config.get('study', {}).get('name', '')
    phenotype = config.get('study', {}).get('phenotype', '')

    # Load data
    ranked_file = output_dir / f"{prefix}_ranked_variants.tsv"
    tissue_file = output_dir / f"{prefix}_tissue_scores.tsv"
    pred_file = intermediate_dir / f"{prefix}_alphagenome_predictions.parquet"

    if not ranked_file.exists():
        logger.error(f"Ranked variants file not found: {ranked_file}")
        return

    ranked_variants = pd.read_csv(ranked_file, sep='\t')

    tissue_scores = pd.DataFrame()
    if tissue_file.exists():
        tissue_scores = pd.read_csv(tissue_file, sep='\t')

    predictions = pd.DataFrame()
    if pred_file.exists():
        predictions = pd.read_parquet(pred_file)
    else:
        tsv_file = intermediate_dir / f"{prefix}_alphagenome_predictions.tsv"
        if tsv_file.exists():
            predictions = pd.read_csv(tsv_file, sep='\t')

    # Generate interactive plots
    output_dir.mkdir(parents=True, exist_ok=True)

    # Individual plots
    create_interactive_manhattan(
        ranked_variants,
        output_path=str(output_dir / f"{prefix}_manhattan_interactive.html")
    )

    if not tissue_scores.empty:
        create_interactive_heatmap(
            tissue_scores,
            output_path=str(output_dir / f"{prefix}_heatmap_interactive.html")
        )

    if not predictions.empty:
        create_interactive_distribution(
            predictions,
            output_path=str(output_dir / f"{prefix}_distribution_interactive.html")
        )

        create_interactive_scatter(
            predictions,
            output_path=str(output_dir / f"{prefix}_scatter_interactive.html")
        )

    # Full dashboard
    if not tissue_scores.empty and not predictions.empty:
        create_interactive_dashboard(
            ranked_variants,
            tissue_scores,
            predictions,
            str(output_dir / f"{prefix}_dashboard.html"),
            study_name,
            phenotype
        )

    logger.info(f"All interactive visualizations saved to {output_dir}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate interactive AlphaGWAS visualizations")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
