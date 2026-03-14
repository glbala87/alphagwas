"""
Colocalization Analysis Module for AlphaGWAS.

Implements statistical colocalization methods to test whether two traits
share a common causal variant at a genomic locus.

Methods implemented:
- COLOC: Bayesian colocalization analysis (Giambartolomei et al., 2014)
- eCAVIAR: Probabilistic colocalization (Hormozdiari et al., 2016)
- Simple correlation-based colocalization
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


@dataclass
class ColocResult:
    """Colocalization analysis result."""

    locus_id: str
    trait1: str
    trait2: str
    n_variants: int
    pp_h0: float  # No association
    pp_h1: float  # Association with trait 1 only
    pp_h2: float  # Association with trait 2 only
    pp_h3: float  # Association with both, different variants
    pp_h4: float  # Association with both, shared variant (colocalization)
    top_coloc_variant: Optional[str] = None
    top_coloc_pp: float = 0.0
    method: str = "coloc"

    @property
    def is_colocalized(self) -> bool:
        """Check if locus shows evidence of colocalization (PP.H4 > 0.8)."""
        return self.pp_h4 > 0.8

    @property
    def summary(self) -> str:
        """Generate human-readable summary."""
        if self.pp_h4 > 0.8:
            return f"Strong colocalization (PP.H4={self.pp_h4:.3f})"
        elif self.pp_h4 > 0.5:
            return f"Moderate colocalization (PP.H4={self.pp_h4:.3f})"
        elif self.pp_h3 > 0.5:
            return f"Distinct causal variants (PP.H3={self.pp_h3:.3f})"
        elif self.pp_h1 > 0.5:
            return f"Association with {self.trait1} only (PP.H1={self.pp_h1:.3f})"
        elif self.pp_h2 > 0.5:
            return f"Association with {self.trait2} only (PP.H2={self.pp_h2:.3f})"
        else:
            return f"No clear signal (PP.H0={self.pp_h0:.3f})"


class ColocAnalyzer:
    """
    Colocalization analysis using COLOC method.

    The COLOC method tests five hypotheses:
    - H0: No association with either trait
    - H1: Association with trait 1 only
    - H2: Association with trait 2 only
    - H3: Association with both traits, different causal variants
    - H4: Association with both traits, shared causal variant
    """

    def __init__(
        self,
        p1: float = 1e-4,
        p2: float = 1e-4,
        p12: float = 1e-5,
    ):
        """
        Initialize COLOC analyzer.

        Args:
            p1: Prior probability of association with trait 1
            p2: Prior probability of association with trait 2
            p12: Prior probability of shared causal variant
        """
        self.p1 = p1
        self.p2 = p2
        self.p12 = p12

    def compute_abf(
        self, beta: np.ndarray, se: np.ndarray, prior_var: float = 0.04
    ) -> np.ndarray:
        """
        Compute Approximate Bayes Factors.

        Args:
            beta: Effect sizes
            se: Standard errors
            prior_var: Prior variance on effect size (W)

        Returns:
            Log ABF values for each variant
        """
        # Variance of beta
        var_beta = se**2

        # Bayes factor calculation
        # ABF = sqrt(V / (V + W)) * exp(z^2 * W / (2 * (V + W)))
        z = beta / se
        r = prior_var / (var_beta + prior_var)

        log_abf = 0.5 * np.log(1 - r) + 0.5 * r * z**2

        return log_abf

    def coloc_analysis(
        self,
        trait1_df: pd.DataFrame,
        trait2_df: pd.DataFrame,
        trait1_name: str = "trait1",
        trait2_name: str = "trait2",
        locus_id: str = "locus",
    ) -> ColocResult:
        """
        Run COLOC analysis for a single locus.

        Args:
            trait1_df: Summary stats for trait 1 (columns: variant_id, beta, se)
            trait2_df: Summary stats for trait 2 (columns: variant_id, beta, se)
            trait1_name: Name of trait 1
            trait2_name: Name of trait 2
            locus_id: Identifier for this locus

        Returns:
            ColocResult with posterior probabilities
        """
        # Merge on variant ID
        merged = pd.merge(
            trait1_df[["variant_id", "beta", "se"]],
            trait2_df[["variant_id", "beta", "se"]],
            on="variant_id",
            suffixes=("_1", "_2"),
        )

        if len(merged) == 0:
            logger.warning(f"No overlapping variants for locus {locus_id}")
            return ColocResult(
                locus_id=locus_id,
                trait1=trait1_name,
                trait2=trait2_name,
                n_variants=0,
                pp_h0=1.0,
                pp_h1=0.0,
                pp_h2=0.0,
                pp_h3=0.0,
                pp_h4=0.0,
            )

        n_variants = len(merged)

        # Compute ABFs for each trait
        log_abf1 = self.compute_abf(merged["beta_1"].values, merged["se_1"].values)
        log_abf2 = self.compute_abf(merged["beta_2"].values, merged["se_2"].values)

        # Compute likelihoods for each hypothesis
        # Use log-sum-exp for numerical stability
        def logsumexp(x):
            max_x = np.max(x)
            return max_x + np.log(np.sum(np.exp(x - max_x)))

        # Prior odds
        log_p1 = np.log(self.p1)
        log_p2 = np.log(self.p2)
        log_p12 = np.log(self.p12)

        # H0: No association
        log_lh0 = 0.0

        # H1: Association with trait 1 only
        log_lh1 = logsumexp(log_abf1) + log_p1

        # H2: Association with trait 2 only
        log_lh2 = logsumexp(log_abf2) + log_p2

        # H3: Both traits, different variants
        log_lh3 = logsumexp(log_abf1) + logsumexp(log_abf2) + log_p1 + log_p2

        # H4: Both traits, same variant
        log_lh4 = logsumexp(log_abf1 + log_abf2) + log_p12

        # Normalize to get posterior probabilities
        log_lhs = np.array([log_lh0, log_lh1, log_lh2, log_lh3, log_lh4])
        log_sum = logsumexp(log_lhs)
        pps = np.exp(log_lhs - log_sum)

        # Find top colocalized variant
        coloc_scores = log_abf1 + log_abf2
        top_idx = np.argmax(coloc_scores)
        top_variant = merged.iloc[top_idx]["variant_id"]

        # Compute per-variant PP for H4
        log_pp_snp = log_abf1 + log_abf2 - logsumexp(log_abf1 + log_abf2)
        pp_snp = np.exp(log_pp_snp)
        top_pp = pp_snp[top_idx] * pps[4]  # Scale by overall H4 probability

        return ColocResult(
            locus_id=locus_id,
            trait1=trait1_name,
            trait2=trait2_name,
            n_variants=n_variants,
            pp_h0=float(pps[0]),
            pp_h1=float(pps[1]),
            pp_h2=float(pps[2]),
            pp_h3=float(pps[3]),
            pp_h4=float(pps[4]),
            top_coloc_variant=top_variant,
            top_coloc_pp=float(top_pp),
            method="coloc",
        )

    def run_multi_locus(
        self,
        trait1_df: pd.DataFrame,
        trait2_df: pd.DataFrame,
        loci: pd.DataFrame,
        trait1_name: str = "trait1",
        trait2_name: str = "trait2",
    ) -> list[ColocResult]:
        """
        Run colocalization analysis across multiple loci.

        Args:
            trait1_df: Full summary stats for trait 1
            trait2_df: Full summary stats for trait 2
            loci: DataFrame with locus definitions (locus_id, chromosome, start, end)
            trait1_name: Name of trait 1
            trait2_name: Name of trait 2

        Returns:
            List of ColocResult for each locus
        """
        results = []

        for _, locus in loci.iterrows():
            locus_id = locus["locus_id"]
            chrom = str(locus["chromosome"])
            start = locus["start"]
            end = locus["end"]

            # Extract variants in locus
            t1_locus = trait1_df[
                (trait1_df["chromosome"].astype(str) == chrom)
                & (trait1_df["position"] >= start)
                & (trait1_df["position"] <= end)
            ].copy()

            t2_locus = trait2_df[
                (trait2_df["chromosome"].astype(str) == chrom)
                & (trait2_df["position"] >= start)
                & (trait2_df["position"] <= end)
            ].copy()

            if len(t1_locus) == 0 or len(t2_locus) == 0:
                logger.warning(f"Skipping locus {locus_id}: insufficient variants")
                continue

            # Create variant IDs if not present
            if "variant_id" not in t1_locus.columns:
                t1_locus["variant_id"] = (
                    "chr" + t1_locus["chromosome"].astype(str) + ":" + t1_locus["position"].astype(str)
                )
            if "variant_id" not in t2_locus.columns:
                t2_locus["variant_id"] = (
                    "chr" + t2_locus["chromosome"].astype(str) + ":" + t2_locus["position"].astype(str)
                )

            result = self.coloc_analysis(
                t1_locus, t2_locus, trait1_name, trait2_name, locus_id
            )
            results.append(result)

            logger.info(f"Locus {locus_id}: {result.summary}")

        return results


class EcaviarAnalyzer:
    """
    eCAVIAR-style colocalization analysis.

    Uses fine-mapping posterior probabilities to compute
    colocalization posterior probability (CLPP).
    """

    def __init__(self, clpp_threshold: float = 0.01):
        """
        Initialize eCAVIAR analyzer.

        Args:
            clpp_threshold: Threshold for significant colocalization
        """
        self.clpp_threshold = clpp_threshold

    def compute_clpp(
        self,
        pip1: np.ndarray,
        pip2: np.ndarray,
    ) -> tuple[float, np.ndarray]:
        """
        Compute colocalization posterior probability.

        Args:
            pip1: Posterior inclusion probabilities for trait 1
            pip2: Posterior inclusion probabilities for trait 2

        Returns:
            Tuple of (total CLPP, per-variant CLPP)
        """
        # Per-variant CLPP is product of PIPs
        clpp_per_variant = pip1 * pip2

        # Total CLPP is sum across variants
        total_clpp = np.sum(clpp_per_variant)

        return total_clpp, clpp_per_variant

    def ecaviar_analysis(
        self,
        trait1_df: pd.DataFrame,
        trait2_df: pd.DataFrame,
        trait1_name: str = "trait1",
        trait2_name: str = "trait2",
        locus_id: str = "locus",
    ) -> ColocResult:
        """
        Run eCAVIAR-style analysis.

        Expects DataFrames with 'variant_id' and 'pip' columns.
        """
        # Merge on variant ID
        merged = pd.merge(
            trait1_df[["variant_id", "pip"]],
            trait2_df[["variant_id", "pip"]],
            on="variant_id",
            suffixes=("_1", "_2"),
        )

        if len(merged) == 0:
            return ColocResult(
                locus_id=locus_id,
                trait1=trait1_name,
                trait2=trait2_name,
                n_variants=0,
                pp_h0=1.0,
                pp_h1=0.0,
                pp_h2=0.0,
                pp_h3=0.0,
                pp_h4=0.0,
                method="ecaviar",
            )

        pip1 = merged["pip_1"].values
        pip2 = merged["pip_2"].values

        total_clpp, clpp_per_variant = self.compute_clpp(pip1, pip2)

        # Find top colocalized variant
        top_idx = np.argmax(clpp_per_variant)
        top_variant = merged.iloc[top_idx]["variant_id"]

        # Convert CLPP to COLOC-style posteriors (approximate)
        # This is a simplification for consistent output format
        pp_h4 = min(total_clpp, 1.0)
        pp_h3 = min((1 - pp_h4) * 0.3, 1.0)
        pp_h1 = min(np.sum(pip1) * (1 - pp_h4 - pp_h3) * 0.5, 1.0)
        pp_h2 = min(np.sum(pip2) * (1 - pp_h4 - pp_h3) * 0.5, 1.0)
        pp_h0 = max(0, 1 - pp_h1 - pp_h2 - pp_h3 - pp_h4)

        # Normalize
        total = pp_h0 + pp_h1 + pp_h2 + pp_h3 + pp_h4
        if total > 0:
            pp_h0, pp_h1, pp_h2, pp_h3, pp_h4 = [p / total for p in [pp_h0, pp_h1, pp_h2, pp_h3, pp_h4]]

        return ColocResult(
            locus_id=locus_id,
            trait1=trait1_name,
            trait2=trait2_name,
            n_variants=len(merged),
            pp_h0=pp_h0,
            pp_h1=pp_h1,
            pp_h2=pp_h2,
            pp_h3=pp_h3,
            pp_h4=pp_h4,
            top_coloc_variant=top_variant,
            top_coloc_pp=float(clpp_per_variant[top_idx]),
            method="ecaviar",
        )


class SimpleColocAnalyzer:
    """
    Simple correlation-based colocalization test.

    Uses correlation of z-scores as a quick colocalization check.
    """

    def analyze(
        self,
        trait1_df: pd.DataFrame,
        trait2_df: pd.DataFrame,
        trait1_name: str = "trait1",
        trait2_name: str = "trait2",
        locus_id: str = "locus",
    ) -> dict:
        """
        Run simple correlation analysis.

        Args:
            trait1_df: Summary stats with beta, se columns
            trait2_df: Summary stats with beta, se columns

        Returns:
            Dictionary with correlation results
        """
        # Merge variants
        merged = pd.merge(
            trait1_df[["variant_id", "beta", "se"]],
            trait2_df[["variant_id", "beta", "se"]],
            on="variant_id",
            suffixes=("_1", "_2"),
        )

        if len(merged) < 3:
            return {
                "locus_id": locus_id,
                "n_variants": len(merged),
                "correlation": np.nan,
                "p_value": np.nan,
                "colocalized": False,
            }

        # Compute z-scores
        z1 = merged["beta_1"] / merged["se_1"]
        z2 = merged["beta_2"] / merged["se_2"]

        # Correlation test
        corr, pval = stats.pearsonr(z1, z2)

        return {
            "locus_id": locus_id,
            "trait1": trait1_name,
            "trait2": trait2_name,
            "n_variants": len(merged),
            "correlation": corr,
            "p_value": pval,
            "colocalized": corr > 0.8 and pval < 0.05,
        }


def run_colocalization(
    gwas1_path: Path,
    gwas2_path: Path,
    loci_path: Optional[Path] = None,
    output_path: Optional[Path] = None,
    method: str = "coloc",
    trait1_name: str = "GWAS1",
    trait2_name: str = "GWAS2",
) -> pd.DataFrame:
    """
    Run colocalization analysis between two GWAS studies.

    Args:
        gwas1_path: Path to GWAS 1 summary statistics
        gwas2_path: Path to GWAS 2 summary statistics
        loci_path: Path to locus definitions (optional)
        output_path: Path to save results (optional)
        method: Analysis method ('coloc', 'ecaviar', 'simple')
        trait1_name: Name for GWAS 1
        trait2_name: Name for GWAS 2

    Returns:
        DataFrame with colocalization results
    """
    logger.info(f"Loading GWAS data from {gwas1_path} and {gwas2_path}")

    # Load data
    gwas1 = pd.read_csv(gwas1_path, sep="\t")
    gwas2 = pd.read_csv(gwas2_path, sep="\t")

    # Standardize column names
    for df in [gwas1, gwas2]:
        if "BETA" in df.columns:
            df.rename(columns={"BETA": "beta"}, inplace=True)
        if "SE" in df.columns:
            df.rename(columns={"SE": "se"}, inplace=True)

    # Define loci if not provided
    if loci_path and loci_path.exists():
        loci = pd.read_csv(loci_path, sep="\t")
    else:
        # Auto-detect significant loci from GWAS1
        logger.info("Auto-detecting loci from significant variants")
        loci = _auto_detect_loci(gwas1)

    logger.info(f"Analyzing {len(loci)} loci with method: {method}")

    # Run analysis
    if method == "coloc":
        analyzer = ColocAnalyzer()
        results = analyzer.run_multi_locus(gwas1, gwas2, loci, trait1_name, trait2_name)
    elif method == "ecaviar":
        analyzer = EcaviarAnalyzer()
        results = []
        for _, locus in loci.iterrows():
            # Would need PIP columns for real eCAVIAR
            result = analyzer.ecaviar_analysis(
                gwas1, gwas2, trait1_name, trait2_name, locus["locus_id"]
            )
            results.append(result)
    else:
        analyzer = SimpleColocAnalyzer()
        results = []
        for _, locus in loci.iterrows():
            result = analyzer.analyze(
                gwas1, gwas2, trait1_name, trait2_name, locus["locus_id"]
            )
            results.append(result)

    # Convert to DataFrame
    if method in ["coloc", "ecaviar"]:
        results_df = pd.DataFrame(
            [
                {
                    "locus_id": r.locus_id,
                    "trait1": r.trait1,
                    "trait2": r.trait2,
                    "n_variants": r.n_variants,
                    "PP.H0": r.pp_h0,
                    "PP.H1": r.pp_h1,
                    "PP.H2": r.pp_h2,
                    "PP.H3": r.pp_h3,
                    "PP.H4": r.pp_h4,
                    "top_variant": r.top_coloc_variant,
                    "top_variant_pp": r.top_coloc_pp,
                    "colocalized": r.is_colocalized,
                    "summary": r.summary,
                    "method": r.method,
                }
                for r in results
            ]
        )
    else:
        results_df = pd.DataFrame(results)

    # Save results
    if output_path:
        results_df.to_csv(output_path, sep="\t", index=False)
        logger.info(f"Results saved to {output_path}")

    # Summary
    if "colocalized" in results_df.columns:
        n_coloc = results_df["colocalized"].sum()
        logger.info(f"Found {n_coloc}/{len(results_df)} colocalized loci")

    return results_df


def _auto_detect_loci(
    gwas_df: pd.DataFrame,
    pval_threshold: float = 5e-8,
    window_kb: int = 500,
) -> pd.DataFrame:
    """Auto-detect loci from significant variants."""
    # Get significant variants
    if "pvalue" in gwas_df.columns:
        sig = gwas_df[gwas_df["pvalue"] < pval_threshold].copy()
    elif "P" in gwas_df.columns:
        sig = gwas_df[gwas_df["P"] < pval_threshold].copy()
    else:
        logger.warning("No p-value column found, using all variants")
        sig = gwas_df.copy()

    if len(sig) == 0:
        logger.warning("No significant variants found")
        return pd.DataFrame(columns=["locus_id", "chromosome", "start", "end"])

    # Create loci around significant variants
    loci = []
    window = window_kb * 1000

    for chrom in sig["chromosome"].unique():
        chrom_sig = sig[sig["chromosome"] == chrom].sort_values("position")

        # Merge overlapping windows
        current_start = None
        current_end = None
        locus_num = 1

        for _, row in chrom_sig.iterrows():
            pos = row["position"]

            if current_start is None:
                current_start = max(0, pos - window)
                current_end = pos + window
            elif pos - window <= current_end:
                # Extend current locus
                current_end = pos + window
            else:
                # Save current locus and start new one
                loci.append(
                    {
                        "locus_id": f"chr{chrom}_locus{locus_num}",
                        "chromosome": chrom,
                        "start": current_start,
                        "end": current_end,
                    }
                )
                locus_num += 1
                current_start = max(0, pos - window)
                current_end = pos + window

        # Save last locus
        if current_start is not None:
            loci.append(
                {
                    "locus_id": f"chr{chrom}_locus{locus_num}",
                    "chromosome": chrom,
                    "start": current_start,
                    "end": current_end,
                }
            )

    return pd.DataFrame(loci)
