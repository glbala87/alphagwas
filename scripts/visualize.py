#!/usr/bin/env python3
"""
Visualization module for AlphaGWAS pipeline.

Creates publication-ready plots:
- Manhattan plots for variant prioritization
- Tissue score heatmaps
- Effect size distributions
- Summary reports
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Optional, Tuple, Any
import yaml

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Try to import visualization libraries
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("matplotlib not installed. Visualization features disabled.")

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def check_dependencies():
    """Check if visualization dependencies are available."""
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError(
            "matplotlib is required for visualization. "
            "Install with: pip install matplotlib seaborn"
        )


def setup_plot_style():
    """Configure matplotlib style for publication-quality plots."""
    if not MATPLOTLIB_AVAILABLE:
        return

    plt.style.use('seaborn-v0_8-whitegrid' if hasattr(plt.style, 'use') else 'ggplot')

    plt.rcParams.update({
        'figure.figsize': (12, 8),
        'figure.dpi': 150,
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16
    })


def create_manhattan_plot(
    variants: pd.DataFrame,
    score_col: str = 'consensus_score',
    output_path: Optional[str] = None,
    title: str = "Variant Prioritization Manhattan Plot",
    highlight_top_n: int = 10,
    figsize: Tuple[int, int] = (14, 6)
) -> Optional[plt.Figure]:
    """
    Create a Manhattan-style plot for variant prioritization scores.

    Args:
        variants: DataFrame with chromosome, position, and score columns
        score_col: Column name for the score to plot
        output_path: Path to save the plot (optional)
        title: Plot title
        highlight_top_n: Number of top variants to highlight
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    check_dependencies()
    setup_plot_style()

    df = variants.copy()

    # Ensure required columns exist
    required = ['chromosome', 'position', score_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        # Try alternative column names
        if 'variant_id' in df.columns and 'chromosome' not in df.columns:
            df[['chromosome', 'position']] = df['variant_id'].str.extract(r'chr(\w+):(\d+)')
            df['position'] = pd.to_numeric(df['position'])

    if score_col not in df.columns:
        logger.error(f"Score column '{score_col}' not found")
        return None

    # Convert chromosome to numeric for plotting
    chrom_order = [str(i) for i in range(1, 23)] + ['X', 'Y']
    df['chrom_num'] = df['chromosome'].astype(str).apply(
        lambda x: chrom_order.index(x) if x in chrom_order else 23
    )

    # Sort by chromosome and position
    df = df.sort_values(['chrom_num', 'position'])

    # Calculate cumulative position for x-axis
    df['cumulative_pos'] = 0
    offset = 0
    chrom_centers = {}

    for chrom in df['chrom_num'].unique():
        mask = df['chrom_num'] == chrom
        df.loc[mask, 'cumulative_pos'] = df.loc[mask, 'position'] + offset
        chrom_centers[chrom] = df.loc[mask, 'cumulative_pos'].median()
        offset = df.loc[mask, 'cumulative_pos'].max() + 1e7

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Color scheme for chromosomes
    colors = ['#1f77b4', '#ff7f0e'] * 12

    # Plot each chromosome
    for i, chrom in enumerate(sorted(df['chrom_num'].unique())):
        mask = df['chrom_num'] == chrom
        ax.scatter(
            df.loc[mask, 'cumulative_pos'],
            df.loc[mask, score_col],
            c=colors[i % 2],
            alpha=0.6,
            s=20,
            edgecolors='none'
        )

    # Highlight top variants
    if highlight_top_n > 0:
        top_variants = df.nlargest(highlight_top_n, score_col)
        ax.scatter(
            top_variants['cumulative_pos'],
            top_variants[score_col],
            c='red',
            s=100,
            alpha=0.9,
            edgecolors='black',
            linewidths=1,
            zorder=5,
            label=f'Top {highlight_top_n}'
        )

        # Add labels for top variants
        for _, row in top_variants.head(5).iterrows():
            label = row.get('rsid', row.get('variant_id', ''))
            if pd.notna(label) and label:
                ax.annotate(
                    label,
                    (row['cumulative_pos'], row[score_col]),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=8,
                    alpha=0.8
                )

    # Set axis labels
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(f'Prioritization Score ({score_col})')
    ax.set_title(title)

    # Set x-axis ticks to chromosome centers
    chrom_labels = [chrom_order[int(c)] if c < len(chrom_order) else '?' for c in sorted(chrom_centers.keys())]
    ax.set_xticks([chrom_centers[c] for c in sorted(chrom_centers.keys())])
    ax.set_xticklabels(chrom_labels)

    # Add legend
    ax.legend(loc='upper right')

    # Add threshold line for top variants
    if len(df) > highlight_top_n:
        threshold = df.nlargest(highlight_top_n, score_col)[score_col].min()
        ax.axhline(y=threshold, color='red', linestyle='--', alpha=0.5, label='Top threshold')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved Manhattan plot to {output_path}")

    return fig


def create_tissue_heatmap(
    tissue_scores: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Tissue-Specific Variant Effect Scores",
    top_n_variants: int = 30,
    figsize: Tuple[int, int] = (12, 10)
) -> Optional[plt.Figure]:
    """
    Create a heatmap of variant scores across tissues.

    Args:
        tissue_scores: DataFrame with variant_id, tissue, and score columns
        output_path: Path to save the plot
        title: Plot title
        top_n_variants: Number of top variants to show
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    check_dependencies()
    setup_plot_style()

    df = tissue_scores.copy()

    # Create pivot table
    score_col = 'tissue_score' if 'tissue_score' in df.columns else 'max_effect'

    # Get top variants by mean score
    variant_means = df.groupby('variant_id')[score_col].mean().nlargest(top_n_variants)
    top_variants = variant_means.index.tolist()

    # Filter to top variants
    df_top = df[df['variant_id'].isin(top_variants)]

    # Create pivot table
    pivot = df_top.pivot_table(
        index='variant_id',
        columns='tissue',
        values=score_col,
        aggfunc='max'
    ).fillna(0)

    # Reorder by mean score
    pivot = pivot.loc[top_variants]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create heatmap
    if SEABORN_AVAILABLE:
        sns.heatmap(
            pivot,
            cmap='YlOrRd',
            ax=ax,
            linewidths=0.5,
            cbar_kws={'label': 'Effect Score'}
        )
    else:
        im = ax.imshow(pivot.values, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        plt.colorbar(im, label='Effect Score')

    ax.set_title(title)
    ax.set_xlabel('Tissue')
    ax.set_ylabel('Variant')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved tissue heatmap to {output_path}")

    return fig


def create_effect_distribution(
    predictions: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Distribution of Predicted Effect Sizes",
    figsize: Tuple[int, int] = (12, 5)
) -> Optional[plt.Figure]:
    """
    Create distribution plots for predicted effect sizes.

    Args:
        predictions: DataFrame with effect_size column
        output_path: Path to save the plot
        title: Plot title
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    check_dependencies()
    setup_plot_style()

    df = predictions.copy()

    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Overall effect size distribution
    ax1 = axes[0]
    effect_col = 'effect_size' if 'effect_size' in df.columns else 'max_effect'

    if SEABORN_AVAILABLE:
        sns.histplot(df[effect_col], bins=50, ax=ax1, kde=True)
    else:
        ax1.hist(df[effect_col].dropna(), bins=50, edgecolor='black', alpha=0.7)

    ax1.set_xlabel('Effect Size')
    ax1.set_ylabel('Count')
    ax1.set_title('Effect Size Distribution')
    ax1.axvline(df[effect_col].median(), color='red', linestyle='--', label=f'Median: {df[effect_col].median():.3f}')
    ax1.legend()

    # Effect by tissue (if available)
    ax2 = axes[1]
    if 'tissue' in df.columns:
        tissue_means = df.groupby('tissue')[effect_col].mean().sort_values(ascending=True)

        if SEABORN_AVAILABLE:
            sns.barplot(x=tissue_means.values, y=tissue_means.index, ax=ax2, palette='viridis')
        else:
            ax2.barh(tissue_means.index, tissue_means.values)

        ax2.set_xlabel('Mean Effect Size')
        ax2.set_ylabel('Tissue')
        ax2.set_title('Mean Effect by Tissue')
    else:
        # Plot by modality if available
        if 'modality' in df.columns:
            modality_means = df.groupby('modality')[effect_col].mean().sort_values(ascending=True)
            ax2.barh(modality_means.index, modality_means.values)
            ax2.set_xlabel('Mean Effect Size')
            ax2.set_ylabel('Modality')
            ax2.set_title('Mean Effect by Modality')

    plt.suptitle(title)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved effect distribution plot to {output_path}")

    return fig


def create_confidence_vs_effect_plot(
    predictions: pd.DataFrame,
    output_path: Optional[str] = None,
    title: str = "Confidence vs Effect Size",
    figsize: Tuple[int, int] = (10, 8)
) -> Optional[plt.Figure]:
    """
    Create a scatter plot of confidence vs effect size.

    Args:
        predictions: DataFrame with confidence and effect_size columns
        output_path: Path to save the plot
        title: Plot title
        figsize: Figure size

    Returns:
        matplotlib Figure object
    """
    check_dependencies()
    setup_plot_style()

    df = predictions.copy()

    if 'confidence' not in df.columns or 'effect_size' not in df.columns:
        logger.warning("Required columns 'confidence' and 'effect_size' not found")
        return None

    fig, ax = plt.subplots(figsize=figsize)

    # Color by tissue if available
    if 'tissue' in df.columns and SEABORN_AVAILABLE:
        sns.scatterplot(
            data=df,
            x='effect_size',
            y='confidence',
            hue='tissue',
            alpha=0.6,
            ax=ax
        )
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.scatter(df['effect_size'], df['confidence'], alpha=0.5, s=20)

    ax.set_xlabel('Effect Size')
    ax.set_ylabel('Confidence')
    ax.set_title(title)

    # Add quadrant lines
    ax.axhline(y=df['confidence'].median(), color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=df['effect_size'].median(), color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved confidence vs effect plot to {output_path}")

    return fig


def create_summary_report(
    ranked_variants: pd.DataFrame,
    tissue_scores: pd.DataFrame,
    predictions: pd.DataFrame,
    output_dir: str,
    prefix: str = "study",
    study_name: str = "",
    phenotype: str = ""
) -> Dict[str, str]:
    """
    Generate all visualizations and create a summary report.

    Args:
        ranked_variants: DataFrame with ranked variants
        tissue_scores: DataFrame with tissue-specific scores
        predictions: DataFrame with AlphaGenome predictions
        output_dir: Directory to save outputs
        prefix: File prefix
        study_name: Study name for report
        phenotype: Phenotype being studied

    Returns:
        Dictionary mapping plot names to file paths
    """
    check_dependencies()

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plots = {}

    # Generate Manhattan plot
    try:
        manhattan_path = output_path / f"{prefix}_manhattan.png"
        create_manhattan_plot(
            ranked_variants,
            output_path=str(manhattan_path),
            title=f"Variant Prioritization - {phenotype}" if phenotype else "Variant Prioritization"
        )
        plots['manhattan'] = str(manhattan_path)
    except Exception as e:
        logger.error(f"Failed to create Manhattan plot: {e}")

    # Generate tissue heatmap
    try:
        heatmap_path = output_path / f"{prefix}_tissue_heatmap.png"
        create_tissue_heatmap(
            tissue_scores,
            output_path=str(heatmap_path),
            title=f"Tissue Effects - {phenotype}" if phenotype else "Tissue Effects"
        )
        plots['tissue_heatmap'] = str(heatmap_path)
    except Exception as e:
        logger.error(f"Failed to create tissue heatmap: {e}")

    # Generate effect distribution
    try:
        dist_path = output_path / f"{prefix}_effect_distribution.png"
        create_effect_distribution(
            predictions,
            output_path=str(dist_path),
            title=f"Effect Size Distribution - {phenotype}" if phenotype else "Effect Size Distribution"
        )
        plots['effect_distribution'] = str(dist_path)
    except Exception as e:
        logger.error(f"Failed to create effect distribution: {e}")

    # Generate confidence vs effect plot
    try:
        conf_path = output_path / f"{prefix}_confidence_effect.png"
        create_confidence_vs_effect_plot(
            predictions,
            output_path=str(conf_path)
        )
        plots['confidence_effect'] = str(conf_path)
    except Exception as e:
        logger.error(f"Failed to create confidence plot: {e}")

    # Close all figures to free memory
    plt.close('all')

    logger.info(f"Generated {len(plots)} visualization plots")
    return plots


def generate_html_report(
    ranked_variants: pd.DataFrame,
    tissue_scores: pd.DataFrame,
    plots: Dict[str, str],
    output_path: str,
    study_name: str = "",
    phenotype: str = "",
    config: Optional[Dict] = None
) -> str:
    """
    Generate an HTML summary report.

    Args:
        ranked_variants: DataFrame with ranked variants
        tissue_scores: DataFrame with tissue scores
        plots: Dictionary of plot paths
        output_path: Path for output HTML file
        study_name: Study name
        phenotype: Phenotype
        config: Pipeline configuration

    Returns:
        Path to generated HTML file
    """
    import base64
    from datetime import datetime

    def embed_image(path: str) -> str:
        """Embed image as base64."""
        try:
            with open(path, 'rb') as f:
                data = base64.b64encode(f.read()).decode()
            return f'data:image/png;base64,{data}'
        except Exception:
            return ''

    # Get top variants
    top_variants = ranked_variants.head(20)

    # Calculate summary statistics
    n_variants = len(ranked_variants)
    n_tissues = tissue_scores['tissue'].nunique() if 'tissue' in tissue_scores.columns else 0
    mean_score = ranked_variants['consensus_score'].mean() if 'consensus_score' in ranked_variants.columns else 0
    max_score = ranked_variants['consensus_score'].max() if 'consensus_score' in ranked_variants.columns else 0

    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>AlphaGWAS Report - {study_name}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 20px;
        }}
        h1 {{ margin: 0; }}
        .subtitle {{ opacity: 0.9; margin-top: 10px; }}
        .card {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 20px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 15px;
            margin-bottom: 20px;
        }}
        .stat-box {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .stat-label {{ color: #666; margin-top: 5px; }}
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{ background: #f8f9fa; font-weight: 600; }}
        tr:hover {{ background: #f8f9fa; }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .section-title {{
            color: #333;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .footer {{
            text-align: center;
            color: #666;
            margin-top: 30px;
            padding: 20px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>AlphaGWAS Analysis Report</h1>
        <div class="subtitle">
            Study: {study_name or 'Unnamed'} | Phenotype: {phenotype or 'Unknown'}<br>
            Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>

    <div class="stats-grid">
        <div class="stat-box">
            <div class="stat-value">{n_variants:,}</div>
            <div class="stat-label">Total Variants</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{n_tissues}</div>
            <div class="stat-label">Tissues Analyzed</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{mean_score:.3f}</div>
            <div class="stat-label">Mean Score</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{max_score:.3f}</div>
            <div class="stat-label">Max Score</div>
        </div>
    </div>

    <div class="card">
        <h2 class="section-title">Top 20 Prioritized Variants</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>Variant</th>
                <th>rsID</th>
                <th>Score</th>
                <th>Max Effect</th>
                <th>Top Tissues</th>
            </tr>
"""

    for _, row in top_variants.iterrows():
        rank = row.get('final_rank', '')
        variant = row.get('variant_id', '')
        rsid = row.get('rsid', '')
        score = row.get('consensus_score', 0)
        max_effect = row.get('max_effect', 0)
        tissues = row.get('top_tissues', '')

        html += f"""
            <tr>
                <td>{rank}</td>
                <td>{variant}</td>
                <td>{rsid}</td>
                <td>{score:.4f}</td>
                <td>{max_effect:.4f}</td>
                <td>{tissues}</td>
            </tr>
"""

    html += """
        </table>
    </div>
"""

    # Add plots
    for plot_name, plot_path in plots.items():
        img_data = embed_image(plot_path)
        if img_data:
            html += f"""
    <div class="card">
        <h2 class="section-title">{plot_name.replace('_', ' ').title()}</h2>
        <div class="plot-container">
            <img src="{img_data}" alt="{plot_name}">
        </div>
    </div>
"""

    html += """
    <div class="footer">
        Generated by AlphaGWAS | Powered by AlphaGenome
    </div>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)

    logger.info(f"Generated HTML report: {output_path}")
    return output_path


def main(config_path: str = "config/config.yaml"):
    """Generate all visualizations for a completed pipeline run."""
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
        logger.error("Run the full pipeline first.")
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

    # Generate plots
    plots = create_summary_report(
        ranked_variants,
        tissue_scores,
        predictions,
        str(output_dir),
        prefix,
        study_name,
        phenotype
    )

    # Generate HTML report
    html_path = output_dir / f"{prefix}_report.html"
    generate_html_report(
        ranked_variants,
        tissue_scores,
        plots,
        str(html_path),
        study_name,
        phenotype,
        config
    )

    logger.info(f"All visualizations saved to {output_dir}")
    return plots


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate AlphaGWAS visualizations")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
