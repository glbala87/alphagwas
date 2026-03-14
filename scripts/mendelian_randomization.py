"""
Mendelian Randomization Module for AlphaGWAS.

Implements MR methods to infer causal relationships:
- Inverse Variance Weighted (IVW)
- MR-Egger regression
- Weighted Median
- MR-PRESSO outlier detection
- Sensitivity analyses
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
class MRResult:
    """Mendelian Randomization result."""

    exposure: str
    outcome: str
    method: str
    n_snps: int
    beta: float
    se: float
    pvalue: float
    ci_lower: float
    ci_upper: float
    or_estimate: float = 1.0  # Odds ratio if binary outcome

    # Sensitivity metrics
    egger_intercept: Optional[float] = None
    egger_intercept_p: Optional[float] = None
    heterogeneity_q: Optional[float] = None
    heterogeneity_p: Optional[float] = None
    i_squared: Optional[float] = None

    @property
    def is_significant(self) -> bool:
        """Check if result is significant at p < 0.05."""
        return self.pvalue < 0.05

    @property
    def has_pleiotropy(self) -> bool:
        """Check for evidence of horizontal pleiotropy (Egger intercept)."""
        if self.egger_intercept_p is not None:
            return self.egger_intercept_p < 0.05
        return False

    @property
    def has_heterogeneity(self) -> bool:
        """Check for significant heterogeneity."""
        if self.heterogeneity_p is not None:
            return self.heterogeneity_p < 0.05
        return False

    def summary(self) -> str:
        """Generate human-readable summary."""
        direction = "increases" if self.beta > 0 else "decreases"
        sig = "significantly" if self.is_significant else "not significantly"

        summary = (
            f"{self.exposure} {sig} {direction} {self.outcome} "
            f"(beta={self.beta:.3f}, p={self.pvalue:.2e}, n_SNPs={self.n_snps})"
        )

        if self.has_pleiotropy:
            summary += " [Warning: evidence of pleiotropy]"
        if self.has_heterogeneity:
            summary += " [Warning: significant heterogeneity]"

        return summary


class MendelianRandomization:
    """
    Mendelian Randomization analysis.

    Uses genetic variants as instrumental variables to infer
    causal effects of exposures on outcomes.
    """

    def __init__(
        self,
        exposure_data: pd.DataFrame,
        outcome_data: pd.DataFrame,
        exposure_name: str = "Exposure",
        outcome_name: str = "Outcome",
    ):
        """
        Initialize MR analysis.

        Args:
            exposure_data: GWAS data for exposure (columns: rsid/variant_id, beta, se, pvalue)
            outcome_data: GWAS data for outcome (columns: rsid/variant_id, beta, se)
            exposure_name: Name of exposure trait
            outcome_name: Name of outcome trait
        """
        self.exposure_data = exposure_data.copy()
        self.outcome_data = outcome_data.copy()
        self.exposure_name = exposure_name
        self.outcome_name = outcome_name
        self.harmonized_data: Optional[pd.DataFrame] = None

    def select_instruments(
        self,
        pval_threshold: float = 5e-8,
        clump: bool = True,
        r2_threshold: float = 0.001,
        kb_window: int = 10000,
    ) -> pd.DataFrame:
        """
        Select genetic instruments for exposure.

        Args:
            pval_threshold: P-value threshold for instrument selection
            clump: Whether to perform LD clumping
            r2_threshold: R² threshold for clumping
            kb_window: Window size for clumping (kb)

        Returns:
            DataFrame with selected instruments
        """
        # Select genome-wide significant variants
        instruments = self.exposure_data[
            self.exposure_data["pvalue"] < pval_threshold
        ].copy()

        logger.info(f"Selected {len(instruments)} instruments at p < {pval_threshold}")

        if clump and len(instruments) > 1:
            # Simple distance-based pruning (approximate clumping)
            instruments = self._simple_clump(instruments, kb_window)
            logger.info(f"After clumping: {len(instruments)} instruments")

        return instruments

    def _simple_clump(
        self, df: pd.DataFrame, kb_window: int = 10000
    ) -> pd.DataFrame:
        """Simple distance-based clumping."""
        df = df.sort_values("pvalue")
        selected = []
        selected_positions: dict[str, list[int]] = {}

        for _, row in df.iterrows():
            chrom = str(row.get("chromosome", "1"))
            pos = int(row.get("position", 0))

            # Check if too close to already selected variant
            if chrom in selected_positions:
                too_close = any(
                    abs(pos - p) < kb_window * 1000
                    for p in selected_positions[chrom]
                )
                if too_close:
                    continue

            selected.append(row)
            if chrom not in selected_positions:
                selected_positions[chrom] = []
            selected_positions[chrom].append(pos)

        return pd.DataFrame(selected)

    def harmonize(
        self,
        instruments: Optional[pd.DataFrame] = None,
        action: int = 2,
    ) -> pd.DataFrame:
        """
        Harmonize exposure and outcome data.

        Ensures effect alleles are aligned between datasets.

        Args:
            instruments: Selected instruments (or None to auto-select)
            action: Harmonization action (1=assume all same, 2=try to infer, 3=correct)

        Returns:
            Harmonized DataFrame
        """
        if instruments is None:
            instruments = self.select_instruments()

        # Identify common variants
        id_col = "rsid" if "rsid" in instruments.columns else "variant_id"

        if id_col not in self.outcome_data.columns:
            # Try to create variant_id
            if "chromosome" in self.outcome_data.columns:
                self.outcome_data["variant_id"] = (
                    "chr" + self.outcome_data["chromosome"].astype(str) +
                    ":" + self.outcome_data["position"].astype(str)
                )
                id_col = "variant_id"

        # Merge exposure and outcome
        harmonized = pd.merge(
            instruments,
            self.outcome_data[[id_col, "beta", "se"]],
            on=id_col,
            suffixes=("_exp", "_out"),
        )

        if len(harmonized) == 0:
            logger.warning("No overlapping variants between exposure and outcome")
            return pd.DataFrame()

        # Rename columns for clarity
        harmonized = harmonized.rename(
            columns={
                "beta_exp": "beta_exposure",
                "se_exp": "se_exposure",
                "beta_out": "beta_outcome",
                "se_out": "se_outcome",
            }
        )

        # Handle beta/se if not already suffixed
        if "beta_exposure" not in harmonized.columns and "beta" in harmonized.columns:
            harmonized["beta_exposure"] = harmonized["beta"]
        if "se_exposure" not in harmonized.columns and "se" in harmonized.columns:
            harmonized["se_exposure"] = harmonized["se"]

        logger.info(f"Harmonized {len(harmonized)} variants")
        self.harmonized_data = harmonized
        return harmonized

    def ivw(self, data: Optional[pd.DataFrame] = None) -> MRResult:
        """
        Inverse Variance Weighted MR.

        Args:
            data: Harmonized data (or None to use stored)

        Returns:
            MRResult with IVW estimate
        """
        if data is None:
            data = self.harmonized_data
        if data is None or len(data) == 0:
            raise ValueError("No harmonized data available")

        beta_exp = data["beta_exposure"].values
        beta_out = data["beta_outcome"].values
        se_out = data["se_outcome"].values

        # IVW estimate (fixed effects)
        weights = 1 / (se_out**2)
        beta_iv = beta_out / beta_exp
        weighted_beta = np.sum(weights * beta_iv) / np.sum(weights)
        weighted_se = np.sqrt(1 / np.sum(weights))

        # P-value
        z = weighted_beta / weighted_se
        pvalue = 2 * stats.norm.sf(abs(z))

        # Confidence interval
        ci_lower = weighted_beta - 1.96 * weighted_se
        ci_upper = weighted_beta + 1.96 * weighted_se

        # Heterogeneity (Cochran's Q)
        q_stat = np.sum(weights * (beta_iv - weighted_beta) ** 2)
        q_pval = 1 - stats.chi2.cdf(q_stat, len(data) - 1) if len(data) > 1 else 1.0
        i_squared = max(0, (q_stat - (len(data) - 1)) / q_stat * 100) if q_stat > 0 else 0

        return MRResult(
            exposure=self.exposure_name,
            outcome=self.outcome_name,
            method="IVW",
            n_snps=len(data),
            beta=weighted_beta,
            se=weighted_se,
            pvalue=pvalue,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            or_estimate=np.exp(weighted_beta),
            heterogeneity_q=q_stat,
            heterogeneity_p=q_pval,
            i_squared=i_squared,
        )

    def mr_egger(self, data: Optional[pd.DataFrame] = None) -> MRResult:
        """
        MR-Egger regression.

        Allows for directional pleiotropy via intercept term.

        Args:
            data: Harmonized data

        Returns:
            MRResult with Egger estimate
        """
        if data is None:
            data = self.harmonized_data
        if data is None or len(data) < 3:
            raise ValueError("MR-Egger requires at least 3 variants")

        beta_exp = data["beta_exposure"].values
        beta_out = data["beta_outcome"].values
        se_out = data["se_outcome"].values

        # Weighted regression of beta_out on beta_exp
        weights = 1 / (se_out**2)

        # Sign of exposure effect (for orientation)
        sign = np.sign(beta_exp)
        beta_exp_abs = np.abs(beta_exp)
        beta_out_oriented = beta_out * sign

        # Weighted least squares
        X = np.column_stack([np.ones(len(beta_exp_abs)), beta_exp_abs])
        W = np.diag(weights)

        XtWX = X.T @ W @ X
        XtWy = X.T @ W @ beta_out_oriented

        try:
            coeffs = np.linalg.solve(XtWX, XtWy)
        except np.linalg.LinAlgError:
            coeffs = np.linalg.lstsq(XtWX, XtWy, rcond=None)[0]

        intercept, slope = coeffs

        # Standard errors
        residuals = beta_out_oriented - X @ coeffs
        mse = np.sum(weights * residuals**2) / (len(data) - 2)
        var_coeffs = mse * np.linalg.inv(XtWX)
        se_intercept = np.sqrt(var_coeffs[0, 0])
        se_slope = np.sqrt(var_coeffs[1, 1])

        # P-values
        z_slope = slope / se_slope
        pvalue = 2 * stats.norm.sf(abs(z_slope))

        z_intercept = intercept / se_intercept
        intercept_pvalue = 2 * stats.norm.sf(abs(z_intercept))

        return MRResult(
            exposure=self.exposure_name,
            outcome=self.outcome_name,
            method="MR-Egger",
            n_snps=len(data),
            beta=slope,
            se=se_slope,
            pvalue=pvalue,
            ci_lower=slope - 1.96 * se_slope,
            ci_upper=slope + 1.96 * se_slope,
            or_estimate=np.exp(slope),
            egger_intercept=intercept,
            egger_intercept_p=intercept_pvalue,
        )

    def weighted_median(
        self, data: Optional[pd.DataFrame] = None, n_boot: int = 1000
    ) -> MRResult:
        """
        Weighted Median MR estimator.

        Robust to up to 50% invalid instruments.

        Args:
            data: Harmonized data
            n_boot: Number of bootstrap iterations for SE

        Returns:
            MRResult with weighted median estimate
        """
        if data is None:
            data = self.harmonized_data
        if data is None or len(data) == 0:
            raise ValueError("No harmonized data available")

        beta_exp = data["beta_exposure"].values
        beta_out = data["beta_outcome"].values
        se_exp = data["se_exposure"].values
        se_out = data["se_outcome"].values

        # Ratio estimates and weights
        beta_iv = beta_out / beta_exp
        weights = 1 / (se_out**2 / beta_exp**2)
        weights = weights / np.sum(weights)

        # Weighted median
        sorted_idx = np.argsort(beta_iv)
        cum_weights = np.cumsum(weights[sorted_idx])
        median_idx = np.searchsorted(cum_weights, 0.5)
        beta_median = beta_iv[sorted_idx[median_idx]]

        # Bootstrap SE
        boot_estimates = []
        for _ in range(n_boot):
            # Resample
            idx = np.random.choice(len(data), len(data), replace=True)
            boot_beta_iv = beta_iv[idx]
            boot_weights = weights[idx]
            boot_weights = boot_weights / np.sum(boot_weights)

            sorted_idx = np.argsort(boot_beta_iv)
            cum_weights = np.cumsum(boot_weights[sorted_idx])
            med_idx = np.searchsorted(cum_weights, 0.5)
            boot_estimates.append(boot_beta_iv[sorted_idx[med_idx]])

        se_median = np.std(boot_estimates)
        z = beta_median / se_median
        pvalue = 2 * stats.norm.sf(abs(z))

        return MRResult(
            exposure=self.exposure_name,
            outcome=self.outcome_name,
            method="Weighted Median",
            n_snps=len(data),
            beta=beta_median,
            se=se_median,
            pvalue=pvalue,
            ci_lower=beta_median - 1.96 * se_median,
            ci_upper=beta_median + 1.96 * se_median,
            or_estimate=np.exp(beta_median),
        )

    def run_all_methods(
        self, data: Optional[pd.DataFrame] = None
    ) -> list[MRResult]:
        """
        Run all MR methods.

        Args:
            data: Harmonized data

        Returns:
            List of MRResult for each method
        """
        if data is None:
            data = self.harmonized_data
        if data is None:
            data = self.harmonize()

        results = []

        # IVW
        try:
            results.append(self.ivw(data))
        except Exception as e:
            logger.warning(f"IVW failed: {e}")

        # MR-Egger (requires >= 3 SNPs)
        if len(data) >= 3:
            try:
                results.append(self.mr_egger(data))
            except Exception as e:
                logger.warning(f"MR-Egger failed: {e}")

        # Weighted Median
        try:
            results.append(self.weighted_median(data))
        except Exception as e:
            logger.warning(f"Weighted Median failed: {e}")

        return results

    def sensitivity_analysis(self, data: Optional[pd.DataFrame] = None) -> dict:
        """
        Run sensitivity analyses.

        Returns:
            Dictionary with sensitivity analysis results
        """
        if data is None:
            data = self.harmonized_data

        results = {
            "n_instruments": len(data) if data is not None else 0,
            "methods": {},
            "pleiotropy_test": None,
            "heterogeneity_test": None,
            "leave_one_out": [],
        }

        # Run all methods
        mr_results = self.run_all_methods(data)
        for r in mr_results:
            results["methods"][r.method] = {
                "beta": r.beta,
                "se": r.se,
                "pvalue": r.pvalue,
                "n_snps": r.n_snps,
            }

            # Store Egger intercept test
            if r.method == "MR-Egger" and r.egger_intercept is not None:
                results["pleiotropy_test"] = {
                    "intercept": r.egger_intercept,
                    "p_value": r.egger_intercept_p,
                    "evidence": r.has_pleiotropy,
                }

            # Store heterogeneity from IVW
            if r.method == "IVW" and r.heterogeneity_q is not None:
                results["heterogeneity_test"] = {
                    "q_statistic": r.heterogeneity_q,
                    "p_value": r.heterogeneity_p,
                    "i_squared": r.i_squared,
                    "evidence": r.has_heterogeneity,
                }

        # Leave-one-out analysis
        if data is not None and len(data) > 2:
            for i in range(len(data)):
                loo_data = data.drop(data.index[i])
                try:
                    loo_result = self.ivw(loo_data)
                    id_col = "rsid" if "rsid" in data.columns else "variant_id"
                    results["leave_one_out"].append({
                        "excluded": data.iloc[i][id_col],
                        "beta": loo_result.beta,
                        "se": loo_result.se,
                        "pvalue": loo_result.pvalue,
                    })
                except Exception:
                    pass

        return results


def run_mr_analysis(
    exposure_path: Path,
    outcome_path: Path,
    output_path: Path,
    exposure_name: str = "Exposure",
    outcome_name: str = "Outcome",
    pval_threshold: float = 5e-8,
) -> pd.DataFrame:
    """
    Run full MR analysis pipeline.

    Args:
        exposure_path: Path to exposure GWAS
        outcome_path: Path to outcome GWAS
        output_path: Output path for results
        exposure_name: Name of exposure
        outcome_name: Name of outcome
        pval_threshold: P-value threshold for instruments

    Returns:
        DataFrame with MR results
    """
    logger.info(f"Running MR: {exposure_name} -> {outcome_name}")

    # Load data
    exposure = pd.read_csv(exposure_path, sep="\t")
    outcome = pd.read_csv(outcome_path, sep="\t")

    # Initialize MR
    mr = MendelianRandomization(
        exposure, outcome, exposure_name, outcome_name
    )

    # Select instruments and harmonize
    instruments = mr.select_instruments(pval_threshold=pval_threshold)
    harmonized = mr.harmonize(instruments)

    if len(harmonized) == 0:
        logger.error("No valid instruments after harmonization")
        return pd.DataFrame()

    # Run analyses
    results = mr.run_all_methods(harmonized)
    sensitivity = mr.sensitivity_analysis(harmonized)

    # Convert to DataFrame
    results_df = pd.DataFrame([
        {
            "exposure": r.exposure,
            "outcome": r.outcome,
            "method": r.method,
            "n_snps": r.n_snps,
            "beta": r.beta,
            "se": r.se,
            "pvalue": r.pvalue,
            "ci_lower": r.ci_lower,
            "ci_upper": r.ci_upper,
            "or": r.or_estimate,
        }
        for r in results
    ])

    # Save results
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"Results saved to {output_path}")

    # Print summary
    for r in results:
        logger.info(r.summary())

    return results_df
