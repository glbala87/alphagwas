#!/usr/bin/env python3
"""
LocusZoom-style regional association plots for AlphaGWAS.

Creates publication-ready regional plots showing:
- Association p-values
- LD structure (r² coloring)
- Gene annotations
- AlphaGenome prediction overlay
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
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


class LocusZoomPlot:
    """
    Generate LocusZoom-style regional association plots.
    """

    # LD color scheme (similar to LocusZoom)
    LD_COLORS = {
        (0.8, 1.0): '#FF0000',    # High LD (red)
        (0.6, 0.8): '#FFA500',    # Medium-high (orange)
        (0.4, 0.6): '#00FF00',    # Medium (green)
        (0.2, 0.4): '#87CEEB',    # Low-medium (light blue)
        (0.0, 0.2): '#0000FF',    # Low LD (blue)
    }

    def __init__(self, config: Optional[Dict] = None):
        """Initialize LocusZoom plotter."""
        self.config = config or {}
        self.flank_kb = self.config.get('flank_kb', 500)

    def get_ld_color(self, r2: float) -> str:
        """Get color based on LD r² value."""
        for (low, high), color in self.LD_COLORS.items():
            if low <= r2 < high:
                return color
        return '#808080'  # Gray for missing

    def create_static_plot(
        self,
        variants: pd.DataFrame,
        lead_variant: str,
        chrom: str,
        start: int,
        end: int,
        output_path: Optional[str] = None,
        title: str = "",
        show_genes: bool = True,
        show_alphagwas: bool = True,
        figsize: Tuple[int, int] = (12, 8)
    ) -> Optional[plt.Figure]:
        """
        Create a static LocusZoom-style plot using matplotlib.

        Args:
            variants: DataFrame with position, pvalue, r2 columns
            lead_variant: Lead variant ID for LD reference
            chrom: Chromosome
            start: Start position
            end: End position
            output_path: Path to save plot
            title: Plot title
            show_genes: Show gene track
            show_alphagwas: Show AlphaGWAS score track
            figsize: Figure size

        Returns:
            matplotlib Figure object
        """
        if not MATPLOTLIB_AVAILABLE:
            logger.error("matplotlib not available")
            return None

        # Filter variants to region
        df = variants.copy()
        if 'position' in df.columns:
            df = df[(df['position'] >= start) & (df['position'] <= end)]

        if df.empty:
            logger.warning(f"No variants in region {chrom}:{start}-{end}")
            return None

        # Calculate -log10(p)
        df['neglog10p'] = -np.log10(df['pvalue'].clip(lower=1e-300))

        # Set up figure with subplots
        n_panels = 1 + int(show_alphagwas) + int(show_genes)
        height_ratios = [3] + ([1] if show_alphagwas else []) + ([0.5] if show_genes else [])

        fig, axes = plt.subplots(
            n_panels, 1,
            figsize=figsize,
            height_ratios=height_ratios,
            sharex=True,
            gridspec_kw={'hspace': 0.05}
        )

        if n_panels == 1:
            axes = [axes]

        ax_idx = 0

        # Main association plot
        ax_assoc = axes[ax_idx]
        ax_idx += 1

        # Color by LD if available
        if 'r2' in df.columns:
            colors = df['r2'].apply(self.get_ld_color)
            sizes = 50 + 100 * df['r2']
        else:
            colors = '#808080'
            sizes = 50

        # Plot variants
        scatter = ax_assoc.scatter(
            df['position'],
            df['neglog10p'],
            c=colors,
            s=sizes,
            alpha=0.8,
            edgecolors='black',
            linewidths=0.5,
            zorder=2
        )

        # Highlight lead variant
        if lead_variant in df['variant_id'].values if 'variant_id' in df.columns else False:
            lead_data = df[df['variant_id'] == lead_variant].iloc[0]
            ax_assoc.scatter(
                lead_data['position'],
                lead_data['neglog10p'],
                c='purple',
                s=200,
                marker='D',
                edgecolors='black',
                linewidths=2,
                zorder=3,
                label='Lead variant'
            )

        # Genome-wide significance line
        ax_assoc.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.7, label='p = 5e-8')

        # Suggestive significance line
        ax_assoc.axhline(y=-np.log10(1e-5), color='blue', linestyle=':', alpha=0.5, label='p = 1e-5')

        ax_assoc.set_ylabel('-log₁₀(p-value)')
        ax_assoc.set_ylim(bottom=0)
        ax_assoc.legend(loc='upper right', fontsize=8)

        # Add LD legend
        ld_patches = [
            mpatches.Patch(color=color, label=f'r² {low:.1f}-{high:.1f}')
            for (low, high), color in sorted(self.LD_COLORS.items(), reverse=True)
        ]
        ax_assoc.legend(handles=ld_patches, loc='upper left', fontsize=7, title='LD (r²)')

        # AlphaGWAS score track
        if show_alphagwas and 'consensus_score' in df.columns:
            ax_alpha = axes[ax_idx]
            ax_idx += 1

            ax_alpha.bar(
                df['position'],
                df['consensus_score'],
                width=(end - start) / len(df) * 0.8,
                color='purple',
                alpha=0.7
            )
            ax_alpha.set_ylabel('AlphaGWAS\nScore')
            ax_alpha.set_ylim(bottom=0)

        # Gene track (simplified)
        if show_genes:
            ax_genes = axes[ax_idx]
            self._draw_gene_track(ax_genes, chrom, start, end)
            ax_genes.set_ylabel('Genes')

        # X-axis formatting
        axes[-1].set_xlabel(f'Position on chromosome {chrom} (Mb)')
        axes[-1].set_xlim(start, end)

        # Format x-axis as Mb
        def format_mb(x, pos):
            return f'{x/1e6:.2f}'
        axes[-1].xaxis.set_major_formatter(plt.FuncFormatter(format_mb))

        # Title
        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold')
        else:
            fig.suptitle(f'Regional Association Plot: {chrom}:{start:,}-{end:,}', fontsize=14)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            logger.info(f"Saved LocusZoom plot to {output_path}")

        return fig

    def _draw_gene_track(self, ax, chrom: str, start: int, end: int):
        """Draw a simplified gene track."""
        # This is a placeholder - in production, load from GTF/GFF
        ax.set_xlim(start, end)
        ax.set_ylim(-0.5, 0.5)
        ax.axhline(y=0, color='gray', linewidth=0.5)
        ax.set_yticks([])
        ax.text(
            (start + end) / 2, 0,
            'Gene annotations (load from GTF)',
            ha='center', va='center',
            fontsize=8, style='italic', color='gray'
        )

    def create_interactive_plot(
        self,
        variants: pd.DataFrame,
        lead_variant: str,
        chrom: str,
        start: int,
        end: int,
        output_path: Optional[str] = None,
        title: str = ""
    ) -> Optional[go.Figure]:
        """
        Create an interactive LocusZoom-style plot using Plotly.

        Args:
            variants: DataFrame with position, pvalue, r2 columns
            lead_variant: Lead variant ID
            chrom: Chromosome
            start: Start position
            end: End position
            output_path: Path to save HTML
            title: Plot title

        Returns:
            Plotly Figure object
        """
        if not PLOTLY_AVAILABLE:
            logger.error("plotly not available")
            return None

        df = variants.copy()
        if 'position' in df.columns:
            df = df[(df['position'] >= start) & (df['position'] <= end)]

        if df.empty:
            logger.warning(f"No variants in region {chrom}:{start}-{end}")
            return None

        df['neglog10p'] = -np.log10(df['pvalue'].clip(lower=1e-300))

        # Create figure with secondary y-axis for AlphaGWAS score
        fig = make_subplots(
            rows=2, cols=1,
            row_heights=[0.7, 0.3],
            shared_xaxes=True,
            vertical_spacing=0.05,
            subplot_titles=('Association', 'AlphaGWAS Score')
        )

        # Color scale for LD
        if 'r2' in df.columns:
            colors = df['r2']
            colorscale = [
                [0.0, '#0000FF'],
                [0.2, '#87CEEB'],
                [0.4, '#00FF00'],
                [0.6, '#FFA500'],
                [0.8, '#FF0000'],
                [1.0, '#FF0000']
            ]
        else:
            colors = '#808080'
            colorscale = None

        # Hover text
        hover_text = df.apply(
            lambda row: (
                f"<b>{row.get('rsid', row.get('variant_id', 'Unknown'))}</b><br>"
                f"Position: {int(row['position']):,}<br>"
                f"P-value: {row['pvalue']:.2e}<br>"
                f"r²: {row.get('r2', 'N/A'):.2f}<br>"
                f"Score: {row.get('consensus_score', 'N/A')}"
            ),
            axis=1
        )

        # Association scatter plot
        fig.add_trace(
            go.Scatter(
                x=df['position'],
                y=df['neglog10p'],
                mode='markers',
                marker=dict(
                    size=8,
                    color=colors,
                    colorscale=colorscale,
                    showscale=True,
                    colorbar=dict(title='r²', x=1.02),
                    line=dict(width=0.5, color='black')
                ),
                text=hover_text,
                hoverinfo='text',
                name='Variants'
            ),
            row=1, col=1
        )

        # Significance lines
        fig.add_hline(y=-np.log10(5e-8), line_dash='dash', line_color='red',
                     annotation_text='p=5e-8', row=1, col=1)
        fig.add_hline(y=-np.log10(1e-5), line_dash='dot', line_color='blue',
                     annotation_text='p=1e-5', row=1, col=1)

        # AlphaGWAS score bar plot
        if 'consensus_score' in df.columns:
            fig.add_trace(
                go.Bar(
                    x=df['position'],
                    y=df['consensus_score'],
                    marker_color='purple',
                    opacity=0.7,
                    name='AlphaGWAS Score',
                    hovertext=hover_text,
                    hoverinfo='text'
                ),
                row=2, col=1
            )

        # Layout
        fig.update_layout(
            title=title or f'Regional Association: chr{chrom}:{start:,}-{end:,}',
            height=600,
            showlegend=False,
            hovermode='closest'
        )

        fig.update_xaxes(title_text=f'Position on chr{chrom} (bp)', row=2, col=1)
        fig.update_yaxes(title_text='-log₁₀(p)', row=1, col=1)
        fig.update_yaxes(title_text='Score', row=2, col=1)

        if output_path:
            fig.write_html(output_path)
            logger.info(f"Saved interactive LocusZoom to {output_path}")

        return fig


def generate_locus_plots(
    ranked_variants: pd.DataFrame,
    config: Dict[str, Any],
    output_dir: str,
    interactive: bool = True
) -> List[str]:
    """
    Generate LocusZoom plots for all defined loci.

    Args:
        ranked_variants: DataFrame with ranked variants
        config: Pipeline configuration
        output_dir: Output directory
        interactive: Generate interactive (Plotly) plots

    Returns:
        List of generated plot paths
    """
    loci = config.get('loci', [])
    if not loci:
        logger.info("No loci defined in config, skipping LocusZoom plots")
        return []

    plotter = LocusZoomPlot(config)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    generated = []

    for locus in loci:
        name = locus['name']
        chrom = str(locus['chromosome'])
        start = locus['start']
        end = locus['end']
        lead_snp = locus.get('lead_snp', '')

        logger.info(f"Generating LocusZoom for {name} ({chrom}:{start}-{end})")

        # Filter variants to this locus
        locus_variants = ranked_variants[
            (ranked_variants['chromosome'].astype(str) == chrom) &
            (ranked_variants['position'] >= start) &
            (ranked_variants['position'] <= end)
        ]

        if locus_variants.empty:
            logger.warning(f"No variants found in locus {name}")
            continue

        # Generate static plot
        static_path = output_path / f"locuszoom_{name}.png"
        plotter.create_static_plot(
            locus_variants,
            lead_variant=lead_snp,
            chrom=chrom,
            start=start,
            end=end,
            output_path=str(static_path),
            title=f"{name} Locus"
        )
        generated.append(str(static_path))

        # Generate interactive plot
        if interactive and PLOTLY_AVAILABLE:
            interactive_path = output_path / f"locuszoom_{name}_interactive.html"
            plotter.create_interactive_plot(
                locus_variants,
                lead_variant=lead_snp,
                chrom=chrom,
                start=start,
                end=end,
                output_path=str(interactive_path),
                title=f"{name} Locus"
            )
            generated.append(str(interactive_path))

    logger.info(f"Generated {len(generated)} LocusZoom plots")
    return generated


def main(config_path: str = "config/config.yaml"):
    """Generate LocusZoom plots for all loci."""
    config = load_config(config_path)

    prefix = config['output']['prefix']
    output_dir = Path(config['output']['dir'])

    # Load ranked variants
    ranked_file = output_dir / f"{prefix}_ranked_variants.tsv"
    if not ranked_file.exists():
        logger.error(f"Ranked variants not found: {ranked_file}")
        return

    ranked_variants = pd.read_csv(ranked_file, sep='\t')

    # Generate plots
    plots = generate_locus_plots(
        ranked_variants,
        config,
        str(output_dir),
        interactive=True
    )

    print(f"\nGenerated {len(plots)} LocusZoom plots")
    for p in plots:
        print(f"  - {p}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate LocusZoom plots")
    parser.add_argument("--config", default="config/config.yaml", help="Config file")
    args = parser.parse_args()

    main(args.config)
