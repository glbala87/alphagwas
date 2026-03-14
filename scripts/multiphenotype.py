"""
Multi-Phenotype Comparison Module for AlphaGWAS.

Compare and visualize variant effects across multiple GWAS traits:
- Cross-trait effect correlations
- Shared vs trait-specific variants
- Pleiotropic variant identification
- Multi-trait visualization
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster import hierarchy

logger = logging.getLogger(__name__)


@dataclass
class PhenotypeData:
    """Container for phenotype GWAS data."""

    name: str
    gwas_df: pd.DataFrame
    ranked_variants: Optional[pd.DataFrame] = None
    tissue_scores: Optional[pd.DataFrame] = None
    color: str = "#2E86AB"

    @property
    def n_significant(self) -> int:
        """Count significant variants."""
        if "pvalue" in self.gwas_df.columns:
            return (self.gwas_df["pvalue"] < 5e-8).sum()
        return 0

    @property
    def variant_ids(self) -> set:
        """Get set of variant IDs."""
        if "variant_id" in self.gwas_df.columns:
            return set(self.gwas_df["variant_id"])
        return set()


class MultiPhenotypeAnalyzer:
    """
    Analyze variants across multiple GWAS phenotypes.

    Identifies:
    - Pleiotropic variants (affecting multiple traits)
    - Trait-specific variants
    - Cross-trait correlations
    - Shared genetic architecture
    """

    def __init__(self, phenotypes: list[PhenotypeData]):
        """
        Initialize analyzer with phenotype data.

        Args:
            phenotypes: List of PhenotypeData objects
        """
        self.phenotypes = phenotypes
        self.n_phenotypes = len(phenotypes)

        # Assign default colors if not set
        colors = ["#2E86AB", "#E94F37", "#44AF69", "#F18F01", "#C73E1D", "#9B59B6"]
        for i, pheno in enumerate(self.phenotypes):
            if pheno.color == "#2E86AB" and i > 0:
                pheno.color = colors[i % len(colors)]

    def identify_pleiotropic_variants(
        self,
        pval_threshold: float = 5e-8,
        min_phenotypes: int = 2,
    ) -> pd.DataFrame:
        """
        Identify variants significant in multiple phenotypes.

        Args:
            pval_threshold: P-value threshold for significance
            min_phenotypes: Minimum number of phenotypes

        Returns:
            DataFrame with pleiotropic variants
        """
        # Collect significant variants from each phenotype
        variant_pheno_map: dict[str, list[str]] = {}
        variant_effects: dict[str, dict[str, dict]] = {}

        for pheno in self.phenotypes:
            df = pheno.gwas_df

            # Get significant variants
            if "pvalue" in df.columns:
                sig_df = df[df["pvalue"] < pval_threshold].copy()
            else:
                continue

            # Create variant ID if needed
            if "variant_id" not in sig_df.columns:
                sig_df["variant_id"] = (
                    "chr" + sig_df["chromosome"].astype(str) +
                    ":" + sig_df["position"].astype(str)
                )

            for _, row in sig_df.iterrows():
                vid = row["variant_id"]
                if vid not in variant_pheno_map:
                    variant_pheno_map[vid] = []
                    variant_effects[vid] = {}

                variant_pheno_map[vid].append(pheno.name)
                variant_effects[vid][pheno.name] = {
                    "beta": row.get("beta", np.nan),
                    "se": row.get("se", np.nan),
                    "pvalue": row.get("pvalue", np.nan),
                }

        # Filter to pleiotropic variants
        pleiotropic = []
        for vid, phenotypes_list in variant_pheno_map.items():
            if len(phenotypes_list) >= min_phenotypes:
                effects = variant_effects[vid]

                # Calculate effect consistency
                betas = [effects[p]["beta"] for p in phenotypes_list if not np.isnan(effects[p]["beta"])]
                effect_direction = "consistent" if len(set(np.sign(betas))) == 1 else "opposing"

                pleiotropic.append({
                    "variant_id": vid,
                    "n_phenotypes": len(phenotypes_list),
                    "phenotypes": ",".join(phenotypes_list),
                    "effect_direction": effect_direction,
                    "min_pvalue": min(effects[p]["pvalue"] for p in phenotypes_list),
                    **{f"{p}_beta": effects.get(p, {}).get("beta", np.nan) for p in [ph.name for ph in self.phenotypes]},
                    **{f"{p}_pval": effects.get(p, {}).get("pvalue", np.nan) for p in [ph.name for ph in self.phenotypes]},
                })

        result = pd.DataFrame(pleiotropic)
        if len(result) > 0:
            result = result.sort_values("n_phenotypes", ascending=False)

        logger.info(f"Found {len(result)} pleiotropic variants affecting {min_phenotypes}+ phenotypes")
        return result

    def calculate_genetic_correlation(
        self,
        method: str = "zscore",
    ) -> pd.DataFrame:
        """
        Calculate genetic correlations between phenotypes.

        Args:
            method: Correlation method ('zscore', 'beta', 'rank')

        Returns:
            DataFrame with pairwise correlations
        """
        # Merge all phenotypes on variant ID
        merged = None

        for pheno in self.phenotypes:
            df = pheno.gwas_df.copy()

            # Create variant ID
            if "variant_id" not in df.columns:
                df["variant_id"] = (
                    "chr" + df["chromosome"].astype(str) +
                    ":" + df["position"].astype(str)
                )

            # Calculate z-score
            if method == "zscore" and "beta" in df.columns and "se" in df.columns:
                df[f"zscore_{pheno.name}"] = df["beta"] / df["se"]
            elif method == "beta" and "beta" in df.columns:
                df[f"beta_{pheno.name}"] = df["beta"]
            elif method == "rank" and "pvalue" in df.columns:
                df[f"rank_{pheno.name}"] = df["pvalue"].rank()

            keep_cols = ["variant_id"] + [c for c in df.columns if c.startswith(("zscore_", "beta_", "rank_"))]
            df = df[keep_cols]

            if merged is None:
                merged = df
            else:
                merged = pd.merge(merged, df, on="variant_id", how="inner")

        if merged is None or len(merged) < 10:
            logger.warning("Insufficient overlapping variants for correlation")
            return pd.DataFrame()

        # Calculate pairwise correlations
        value_cols = [c for c in merged.columns if c != "variant_id"]
        corr_matrix = merged[value_cols].corr()

        # Clean up column names
        corr_matrix.columns = [c.split("_", 1)[1] for c in corr_matrix.columns]
        corr_matrix.index = [c.split("_", 1)[1] for c in corr_matrix.index]

        logger.info(f"Calculated {method} correlations for {len(merged)} shared variants")
        return corr_matrix

    def compare_effect_sizes(
        self,
        pheno1_name: str,
        pheno2_name: str,
        variants: Optional[list[str]] = None,
    ) -> pd.DataFrame:
        """
        Compare effect sizes between two phenotypes.

        Args:
            pheno1_name: Name of first phenotype
            pheno2_name: Name of second phenotype
            variants: Optional list of variants to compare

        Returns:
            DataFrame with effect size comparison
        """
        pheno1 = next((p for p in self.phenotypes if p.name == pheno1_name), None)
        pheno2 = next((p for p in self.phenotypes if p.name == pheno2_name), None)

        if not pheno1 or not pheno2:
            raise ValueError(f"Phenotype not found: {pheno1_name} or {pheno2_name}")

        # Prepare data
        df1 = pheno1.gwas_df.copy()
        df2 = pheno2.gwas_df.copy()

        for df in [df1, df2]:
            if "variant_id" not in df.columns:
                df["variant_id"] = (
                    "chr" + df["chromosome"].astype(str) +
                    ":" + df["position"].astype(str)
                )

        # Merge
        merged = pd.merge(
            df1[["variant_id", "beta", "se", "pvalue"]],
            df2[["variant_id", "beta", "se", "pvalue"]],
            on="variant_id",
            suffixes=(f"_{pheno1_name}", f"_{pheno2_name}"),
        )

        if variants:
            merged = merged[merged["variant_id"].isin(variants)]

        # Calculate comparison metrics
        merged["beta_ratio"] = merged[f"beta_{pheno1_name}"] / merged[f"beta_{pheno2_name}"]
        merged["same_direction"] = np.sign(merged[f"beta_{pheno1_name}"]) == np.sign(merged[f"beta_{pheno2_name}"])

        # Z-score for heterogeneity
        beta_diff = merged[f"beta_{pheno1_name}"] - merged[f"beta_{pheno2_name}"]
        se_combined = np.sqrt(merged[f"se_{pheno1_name}"]**2 + merged[f"se_{pheno2_name}"]**2)
        merged["heterogeneity_z"] = beta_diff / se_combined
        merged["heterogeneity_p"] = 2 * stats.norm.sf(np.abs(merged["heterogeneity_z"]))

        return merged

    def create_upset_data(
        self,
        pval_threshold: float = 5e-8,
    ) -> pd.DataFrame:
        """
        Create data for UpSet plot of variant overlap.

        Args:
            pval_threshold: P-value threshold

        Returns:
            DataFrame with set membership for each variant
        """
        all_variants: set[str] = set()
        pheno_variants: dict[str, set[str]] = {}

        for pheno in self.phenotypes:
            df = pheno.gwas_df
            if "pvalue" in df.columns:
                sig = df[df["pvalue"] < pval_threshold].copy()
            else:
                sig = df.copy()

            if "variant_id" not in sig.columns:
                sig["variant_id"] = (
                    "chr" + sig["chromosome"].astype(str) +
                    ":" + sig["position"].astype(str)
                )

            pheno_variants[pheno.name] = set(sig["variant_id"])
            all_variants.update(pheno_variants[pheno.name])

        # Create membership matrix
        upset_data = []
        for vid in all_variants:
            row = {"variant_id": vid}
            for pheno in self.phenotypes:
                row[pheno.name] = vid in pheno_variants[pheno.name]
            upset_data.append(row)

        return pd.DataFrame(upset_data)

    def cluster_phenotypes(
        self,
        method: str = "average",
        metric: str = "correlation",
    ) -> dict:
        """
        Hierarchically cluster phenotypes by genetic similarity.

        Args:
            method: Clustering method
            metric: Distance metric

        Returns:
            Dictionary with clustering results
        """
        corr_matrix = self.calculate_genetic_correlation()

        if corr_matrix.empty:
            return {"error": "Could not calculate correlations"}

        # Convert correlation to distance
        distance_matrix = 1 - corr_matrix.abs()

        # Hierarchical clustering
        linkage = hierarchy.linkage(
            distance_matrix.values[np.triu_indices(len(distance_matrix), k=1)],
            method=method,
        )

        # Get dendrogram order
        dendro = hierarchy.dendrogram(linkage, no_plot=True, labels=list(corr_matrix.columns))

        return {
            "linkage": linkage,
            "order": dendro["leaves"],
            "phenotype_order": [corr_matrix.columns[i] for i in dendro["leaves"]],
            "correlation_matrix": corr_matrix,
        }

    def generate_summary(self) -> dict:
        """Generate summary statistics across phenotypes."""
        summary = {
            "n_phenotypes": self.n_phenotypes,
            "phenotypes": [],
            "total_unique_variants": 0,
            "shared_variants": 0,
        }

        all_variants: set[str] = set()
        shared_variants: Optional[set[str]] = None

        for pheno in self.phenotypes:
            pheno_summary = {
                "name": pheno.name,
                "n_variants": len(pheno.gwas_df),
                "n_significant": pheno.n_significant,
            }
            summary["phenotypes"].append(pheno_summary)

            if "variant_id" not in pheno.gwas_df.columns:
                continue

            variants = set(pheno.gwas_df["variant_id"])
            all_variants.update(variants)

            if shared_variants is None:
                shared_variants = variants
            else:
                shared_variants &= variants

        summary["total_unique_variants"] = len(all_variants)
        summary["shared_variants"] = len(shared_variants) if shared_variants else 0

        return summary


class MultiPhenotypeVisualizer:
    """Visualization tools for multi-phenotype analysis."""

    def __init__(self, analyzer: MultiPhenotypeAnalyzer):
        """Initialize visualizer."""
        self.analyzer = analyzer

    def plot_effect_comparison(
        self,
        pheno1_name: str,
        pheno2_name: str,
        highlight_pleiotropic: bool = True,
        figsize: tuple = (8, 8),
    ):
        """
        Create scatter plot comparing effects between two phenotypes.

        Args:
            pheno1_name: Name of first phenotype
            pheno2_name: Name of second phenotype
            highlight_pleiotropic: Highlight pleiotropic variants
            figsize: Figure size

        Returns:
            matplotlib Figure
        """
        import matplotlib.pyplot as plt

        comparison = self.analyzer.compare_effect_sizes(pheno1_name, pheno2_name)

        fig, ax = plt.subplots(figsize=figsize)

        # Plot all variants
        ax.scatter(
            comparison[f"beta_{pheno1_name}"],
            comparison[f"beta_{pheno2_name}"],
            alpha=0.5,
            s=20,
            c="#CCCCCC",
            label="All variants",
        )

        # Highlight significant in both
        sig_both = comparison[
            (comparison[f"pvalue_{pheno1_name}"] < 5e-8) &
            (comparison[f"pvalue_{pheno2_name}"] < 5e-8)
        ]

        if len(sig_both) > 0:
            ax.scatter(
                sig_both[f"beta_{pheno1_name}"],
                sig_both[f"beta_{pheno2_name}"],
                alpha=0.8,
                s=50,
                c="#E94F37",
                label=f"Significant in both (n={len(sig_both)})",
                edgecolors="black",
                linewidth=0.5,
            )

        # Add diagonal line
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, "--", color="gray", alpha=0.5)
        ax.axhline(0, color="gray", linewidth=0.5)
        ax.axvline(0, color="gray", linewidth=0.5)

        # Correlation annotation
        corr, pval = stats.pearsonr(
            comparison[f"beta_{pheno1_name}"],
            comparison[f"beta_{pheno2_name}"],
        )
        ax.annotate(
            f"r = {corr:.3f}\np = {pval:.2e}",
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        ax.set_xlabel(f"Effect size ({pheno1_name})", fontsize=12)
        ax.set_ylabel(f"Effect size ({pheno2_name})", fontsize=12)
        ax.set_title(f"Effect Size Comparison: {pheno1_name} vs {pheno2_name}", fontsize=14)
        ax.legend(loc="lower right")

        plt.tight_layout()
        return fig

    def plot_correlation_heatmap(self, figsize: tuple = (10, 8)):
        """
        Create correlation heatmap across phenotypes.

        Args:
            figsize: Figure size

        Returns:
            matplotlib Figure
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        corr_matrix = self.analyzer.calculate_genetic_correlation()

        if corr_matrix.empty:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Insufficient data for correlation", ha="center")
            return fig

        fig, ax = plt.subplots(figsize=figsize)

        sns.heatmap(
            corr_matrix,
            annot=True,
            fmt=".2f",
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            ax=ax,
            cbar_kws={"label": "Genetic Correlation"},
        )

        ax.set_title("Cross-Phenotype Genetic Correlations", fontsize=14)
        plt.tight_layout()

        return fig

    def plot_pleiotropic_summary(self, pleiotropic_df: pd.DataFrame, figsize: tuple = (12, 5)):
        """
        Create summary visualization of pleiotropic variants.

        Args:
            pleiotropic_df: DataFrame from identify_pleiotropic_variants
            figsize: Figure size

        Returns:
            matplotlib Figure
        """
        import matplotlib.pyplot as plt

        if len(pleiotropic_df) == 0:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "No pleiotropic variants found", ha="center")
            return fig

        fig, axes = plt.subplots(1, 2, figsize=figsize)

        # Bar chart of variants by number of phenotypes
        pheno_counts = pleiotropic_df["n_phenotypes"].value_counts().sort_index()
        axes[0].bar(pheno_counts.index, pheno_counts.values, color="#2E86AB", edgecolor="black")
        axes[0].set_xlabel("Number of Phenotypes", fontsize=11)
        axes[0].set_ylabel("Number of Variants", fontsize=11)
        axes[0].set_title("Pleiotropic Variant Distribution", fontsize=12)

        # Effect direction pie chart
        direction_counts = pleiotropic_df["effect_direction"].value_counts()
        colors = ["#44AF69", "#E94F37"]
        axes[1].pie(
            direction_counts.values,
            labels=direction_counts.index,
            autopct="%1.1f%%",
            colors=colors,
            startangle=90,
        )
        axes[1].set_title("Effect Direction Consistency", fontsize=12)

        plt.tight_layout()
        return fig

    def plot_miami(
        self,
        pheno1_name: str,
        pheno2_name: str,
        figsize: tuple = (14, 8),
    ):
        """
        Create Miami plot (mirrored Manhattan plots) for two phenotypes.

        Args:
            pheno1_name: Name of first phenotype (top)
            pheno2_name: Name of second phenotype (bottom)
            figsize: Figure size

        Returns:
            matplotlib Figure
        """
        import matplotlib.pyplot as plt

        pheno1 = next((p for p in self.analyzer.phenotypes if p.name == pheno1_name), None)
        pheno2 = next((p for p in self.analyzer.phenotypes if p.name == pheno2_name), None)

        if not pheno1 or not pheno2:
            raise ValueError("Phenotype not found")

        fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)

        for ax, pheno, direction in zip(axes, [pheno1, pheno2], [1, -1]):
            df = pheno.gwas_df.copy()

            if "pvalue" not in df.columns:
                continue

            # Calculate -log10(p)
            df["-log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

            # Create cumulative position
            if "chromosome" in df.columns and "position" in df.columns:
                df["chromosome"] = df["chromosome"].astype(str)
                chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
                df["chrom_num"] = df["chromosome"].map(
                    {c: i for i, c in enumerate(chrom_order)}
                )
                df = df.sort_values(["chrom_num", "position"])

                # Add offset for each chromosome
                chrom_sizes = df.groupby("chromosome")["position"].max()
                offsets = {}
                cumulative = 0
                for chrom in chrom_order:
                    if chrom in chrom_sizes.index:
                        offsets[chrom] = cumulative
                        cumulative += chrom_sizes[chrom]

                df["cumulative_pos"] = df.apply(
                    lambda r: offsets.get(r["chromosome"], 0) + r["position"], axis=1
                )

                # Plot
                y_values = df["-log10p"] * direction
                colors = df["chrom_num"] % 2

                ax.scatter(
                    df["cumulative_pos"],
                    y_values,
                    c=colors,
                    cmap="Set1",
                    s=10,
                    alpha=0.6,
                )

                # Significance line
                ax.axhline(y=-np.log10(5e-8) * direction, color="red", linestyle="--", alpha=0.7)

            ax.set_ylabel(f"-log10(p)\n{pheno.name}", fontsize=10)
            ax.spines["top" if direction == 1 else "bottom"].set_visible(False)

        # Format
        axes[0].set_title("Miami Plot: Multi-Phenotype Comparison", fontsize=14)
        axes[1].set_xlabel("Genomic Position", fontsize=11)
        axes[1].invert_yaxis()

        plt.tight_layout()
        return fig


def run_multiphenotype_analysis(
    gwas_files: list[Path],
    phenotype_names: list[str],
    output_dir: Path,
    pval_threshold: float = 5e-8,
) -> dict[str, pd.DataFrame]:
    """
    Run full multi-phenotype analysis.

    Args:
        gwas_files: List of paths to GWAS files
        phenotype_names: Names for each phenotype
        output_dir: Output directory
        pval_threshold: P-value threshold

    Returns:
        Dictionary of result DataFrames
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load phenotypes
    phenotypes = []
    for path, name in zip(gwas_files, phenotype_names):
        logger.info(f"Loading {name} from {path}")
        df = pd.read_csv(path, sep="\t")
        phenotypes.append(PhenotypeData(name=name, gwas_df=df))

    # Initialize analyzer
    analyzer = MultiPhenotypeAnalyzer(phenotypes)

    # Run analyses
    results = {}

    # Pleiotropic variants
    pleiotropic = analyzer.identify_pleiotropic_variants(pval_threshold=pval_threshold)
    pleiotropic.to_csv(output_dir / "pleiotropic_variants.tsv", sep="\t", index=False)
    results["pleiotropic"] = pleiotropic

    # Genetic correlations
    correlations = analyzer.calculate_genetic_correlation()
    correlations.to_csv(output_dir / "genetic_correlations.tsv", sep="\t")
    results["correlations"] = correlations

    # UpSet data
    upset = analyzer.create_upset_data(pval_threshold=pval_threshold)
    upset.to_csv(output_dir / "variant_overlap.tsv", sep="\t", index=False)
    results["upset"] = upset

    # Summary
    summary = analyzer.generate_summary()
    results["summary"] = summary

    # Visualizations
    visualizer = MultiPhenotypeVisualizer(analyzer)

    try:
        fig = visualizer.plot_correlation_heatmap()
        fig.savefig(output_dir / "correlation_heatmap.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        logger.warning(f"Could not create correlation heatmap: {e}")

    if len(pleiotropic) > 0:
        try:
            fig = visualizer.plot_pleiotropic_summary(pleiotropic)
            fig.savefig(output_dir / "pleiotropic_summary.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
        except Exception as e:
            logger.warning(f"Could not create pleiotropic summary: {e}")

    logger.info(f"Multi-phenotype analysis complete. Results in {output_dir}")
    return results
