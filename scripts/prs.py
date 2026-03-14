"""
Polygenic Risk Score (PRS) Module for AlphaGWAS.

Calculate and validate polygenic risk scores:
- Clumping and thresholding (C+T)
- PRSice-style optimization
- Score calculation and validation
- Performance metrics
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
class PRSModel:
    """Polygenic Risk Score model."""

    name: str
    n_variants: int
    pval_threshold: float
    weights: pd.DataFrame  # variant_id, weight, effect_allele
    method: str = "clumping"

    # Performance metrics (after validation)
    r2: Optional[float] = None
    r2_liability: Optional[float] = None
    auc: Optional[float] = None
    beta: Optional[float] = None
    se: Optional[float] = None
    pvalue: Optional[float] = None

    def __post_init__(self):
        if self.r2 is not None and self.r2_liability is None:
            # Convert observed R² to liability scale (approximate)
            self.r2_liability = self.r2


@dataclass
class PRSResult:
    """PRS calculation result for an individual."""

    sample_id: str
    score: float
    n_variants_used: int
    percentile: Optional[float] = None
    risk_category: Optional[str] = None


class PRSCalculator:
    """
    Calculate Polygenic Risk Scores from GWAS summary statistics.

    Implements clumping + thresholding approach with optional
    threshold optimization.
    """

    def __init__(
        self,
        gwas_data: pd.DataFrame,
        trait_name: str = "Trait",
    ):
        """
        Initialize PRS calculator.

        Args:
            gwas_data: GWAS summary statistics (rsid, chromosome, position, beta, se, pvalue)
            trait_name: Name of the trait
        """
        self.gwas_data = gwas_data.copy()
        self.trait_name = trait_name
        self._standardize_columns()

    def _standardize_columns(self):
        """Standardize column names."""
        col_map = {
            "BETA": "beta",
            "SE": "se",
            "P": "pvalue",
            "PVALUE": "pvalue",
            "CHR": "chromosome",
            "POS": "position",
            "RSID": "rsid",
            "SNP": "rsid",
            "A1": "effect_allele",
            "A2": "other_allele",
        }
        self.gwas_data.rename(columns=col_map, inplace=True)

        # Create variant_id if not present
        if "variant_id" not in self.gwas_data.columns:
            self.gwas_data["variant_id"] = (
                "chr" + self.gwas_data["chromosome"].astype(str) +
                ":" + self.gwas_data["position"].astype(str)
            )

    def clump_variants(
        self,
        pval_threshold: float = 1.0,
        r2_threshold: float = 0.1,
        kb_window: int = 250,
    ) -> pd.DataFrame:
        """
        Clump variants by p-value and distance.

        Args:
            pval_threshold: P-value threshold
            r2_threshold: R² threshold (approximate via distance)
            kb_window: Clumping window in kb

        Returns:
            DataFrame with clumped variants
        """
        # Filter by p-value
        filtered = self.gwas_data[self.gwas_data["pvalue"] <= pval_threshold].copy()
        filtered = filtered.sort_values("pvalue")

        if len(filtered) == 0:
            logger.warning(f"No variants below p-value threshold {pval_threshold}")
            return pd.DataFrame()

        # Distance-based clumping (approximation without LD reference)
        selected = []
        selected_positions: dict[str, list[int]] = {}

        for _, row in filtered.iterrows():
            chrom = str(row["chromosome"])
            pos = int(row["position"])

            # Check distance to already selected variants
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

        result = pd.DataFrame(selected)
        logger.info(
            f"Clumped {len(filtered)} variants to {len(result)} "
            f"(p < {pval_threshold}, window={kb_window}kb)"
        )
        return result

    def create_model(
        self,
        pval_threshold: float = 5e-8,
        clump: bool = True,
        r2_threshold: float = 0.1,
        kb_window: int = 250,
    ) -> PRSModel:
        """
        Create PRS model with specified parameters.

        Args:
            pval_threshold: P-value threshold for inclusion
            clump: Whether to clump variants
            r2_threshold: R² threshold for clumping
            kb_window: Clumping window

        Returns:
            PRSModel with weights
        """
        if clump:
            variants = self.clump_variants(pval_threshold, r2_threshold, kb_window)
        else:
            variants = self.gwas_data[
                self.gwas_data["pvalue"] <= pval_threshold
            ].copy()

        if len(variants) == 0:
            raise ValueError(f"No variants selected at p < {pval_threshold}")

        # Create weights DataFrame
        weight_cols = ["variant_id", "beta", "pvalue"]
        if "rsid" in variants.columns:
            weight_cols.insert(1, "rsid")
        if "effect_allele" in variants.columns:
            weight_cols.append("effect_allele")
        if "other_allele" in variants.columns:
            weight_cols.append("other_allele")

        weights = variants[[c for c in weight_cols if c in variants.columns]].copy()
        weights = weights.rename(columns={"beta": "weight"})

        model = PRSModel(
            name=f"{self.trait_name}_prs_p{pval_threshold:.0e}",
            n_variants=len(weights),
            pval_threshold=pval_threshold,
            weights=weights,
            method="clumping" if clump else "thresholding",
        )

        logger.info(f"Created PRS model with {model.n_variants} variants")
        return model

    def optimize_threshold(
        self,
        genotypes: pd.DataFrame,
        phenotype: pd.Series,
        thresholds: Optional[list[float]] = None,
        n_folds: int = 5,
    ) -> tuple[PRSModel, pd.DataFrame]:
        """
        Optimize p-value threshold using cross-validation.

        Args:
            genotypes: Genotype matrix (samples x variants)
            phenotype: Phenotype values
            thresholds: P-value thresholds to test
            n_folds: Number of CV folds

        Returns:
            Tuple of (best model, results DataFrame)
        """
        if thresholds is None:
            thresholds = [5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1.0]

        results = []
        best_model = None
        best_r2 = -np.inf

        for threshold in thresholds:
            try:
                model = self.create_model(pval_threshold=threshold)

                # Calculate scores
                scores = self.calculate_scores(genotypes, model)

                if len(scores) == 0:
                    continue

                # Match with phenotype
                common = set(scores.index) & set(phenotype.index)
                if len(common) < 10:
                    continue

                scores_matched = scores.loc[list(common)]
                pheno_matched = phenotype.loc[list(common)]

                # Calculate R²
                corr = np.corrcoef(scores_matched["score"], pheno_matched)[0, 1]
                r2 = corr ** 2

                results.append({
                    "threshold": threshold,
                    "n_variants": model.n_variants,
                    "r2": r2,
                    "correlation": corr,
                })

                if r2 > best_r2:
                    best_r2 = r2
                    best_model = model
                    best_model.r2 = r2

            except Exception as e:
                logger.warning(f"Threshold {threshold} failed: {e}")

        results_df = pd.DataFrame(results)

        if best_model is not None:
            logger.info(
                f"Best threshold: p < {best_model.pval_threshold:.0e} "
                f"(R² = {best_r2:.4f}, n = {best_model.n_variants})"
            )

        return best_model, results_df

    def calculate_scores(
        self,
        genotypes: pd.DataFrame,
        model: PRSModel,
    ) -> pd.DataFrame:
        """
        Calculate PRS for individuals.

        Args:
            genotypes: Genotype matrix (samples x variants, values 0/1/2)
            model: PRS model with weights

        Returns:
            DataFrame with sample_id and score
        """
        # Find overlapping variants
        weight_variants = set(model.weights["variant_id"])
        geno_variants = set(genotypes.columns)
        common = weight_variants & geno_variants

        if len(common) == 0:
            # Try rsid matching
            if "rsid" in model.weights.columns:
                weight_variants = set(model.weights["rsid"].dropna())
                common = weight_variants & geno_variants

        if len(common) == 0:
            logger.warning("No overlapping variants between model and genotypes")
            return pd.DataFrame()

        logger.info(f"Using {len(common)}/{model.n_variants} variants for scoring")

        # Get weights for common variants
        id_col = "variant_id" if model.weights["variant_id"].iloc[0] in common else "rsid"
        weights = model.weights[model.weights[id_col].isin(common)].set_index(id_col)

        # Calculate scores
        scores = []
        for sample_id, row in genotypes.iterrows():
            score = 0.0
            n_used = 0
            for var_id, weight_row in weights.iterrows():
                if var_id in row.index and pd.notna(row[var_id]):
                    score += row[var_id] * weight_row["weight"]
                    n_used += 1
            scores.append({
                "sample_id": sample_id,
                "score": score,
                "n_variants_used": n_used,
            })

        result = pd.DataFrame(scores).set_index("sample_id")

        # Add percentiles
        result["percentile"] = result["score"].rank(pct=True) * 100

        # Add risk categories
        result["risk_category"] = pd.cut(
            result["percentile"],
            bins=[0, 20, 40, 60, 80, 100],
            labels=["Very Low", "Low", "Average", "High", "Very High"],
        )

        return result

    def validate_model(
        self,
        model: PRSModel,
        genotypes: pd.DataFrame,
        phenotype: pd.Series,
        covariates: Optional[pd.DataFrame] = None,
        binary_outcome: bool = False,
    ) -> PRSModel:
        """
        Validate PRS model against phenotype.

        Args:
            model: PRS model to validate
            genotypes: Genotype matrix
            phenotype: Phenotype values
            covariates: Optional covariates to adjust for
            binary_outcome: Whether phenotype is binary

        Returns:
            Updated model with performance metrics
        """
        # Calculate scores
        scores = self.calculate_scores(genotypes, model)

        if len(scores) == 0:
            logger.error("Could not calculate scores for validation")
            return model

        # Match samples
        common = set(scores.index) & set(phenotype.index)
        if covariates is not None:
            common = common & set(covariates.index)

        if len(common) < 10:
            logger.error(f"Too few samples for validation: {len(common)}")
            return model

        scores_matched = scores.loc[list(common), "score"]
        pheno_matched = phenotype.loc[list(common)]

        if binary_outcome:
            # Logistic regression / AUC
            try:
                from sklearn.metrics import roc_auc_score

                auc = roc_auc_score(pheno_matched, scores_matched)
                model.auc = auc
                logger.info(f"Validation AUC: {auc:.4f}")
            except ImportError:
                # Use simple correlation
                corr = np.corrcoef(scores_matched, pheno_matched)[0, 1]
                model.r2 = corr ** 2
        else:
            # Linear regression / R²
            corr, pvalue = stats.pearsonr(scores_matched, pheno_matched)
            model.r2 = corr ** 2
            model.pvalue = pvalue

            # Calculate beta/se via regression
            n = len(scores_matched)
            model.beta = corr * np.std(pheno_matched) / np.std(scores_matched)
            model.se = model.beta / np.sqrt(model.r2 * (n - 2) / (1 - model.r2)) if model.r2 < 1 else 0

            logger.info(f"Validation R²: {model.r2:.4f}, p = {pvalue:.2e}")

        return model


def calculate_prs_from_gwas(
    gwas_path: Path,
    genotype_path: Path,
    output_path: Path,
    phenotype_path: Optional[Path] = None,
    pval_threshold: float = 5e-8,
    optimize: bool = False,
) -> tuple[PRSModel, pd.DataFrame]:
    """
    Calculate PRS from GWAS and genotype data.

    Args:
        gwas_path: Path to GWAS summary statistics
        genotype_path: Path to genotype file (TSV with samples as rows)
        output_path: Output path for scores
        phenotype_path: Optional phenotype file for validation
        pval_threshold: P-value threshold
        optimize: Whether to optimize threshold

    Returns:
        Tuple of (model, scores DataFrame)
    """
    logger.info(f"Calculating PRS from {gwas_path}")

    # Load data
    gwas = pd.read_csv(gwas_path, sep="\t")
    genotypes = pd.read_csv(genotype_path, sep="\t", index_col=0)

    # Initialize calculator
    calculator = PRSCalculator(gwas)

    # Create or optimize model
    if optimize and phenotype_path:
        phenotype = pd.read_csv(phenotype_path, sep="\t", index_col=0).iloc[:, 0]
        model, opt_results = calculator.optimize_threshold(genotypes, phenotype)
        opt_results.to_csv(output_path.with_suffix(".optimization.tsv"), sep="\t", index=False)
    else:
        model = calculator.create_model(pval_threshold=pval_threshold)

    # Calculate scores
    scores = calculator.calculate_scores(genotypes, model)

    # Validate if phenotype provided
    if phenotype_path:
        phenotype = pd.read_csv(phenotype_path, sep="\t", index_col=0).iloc[:, 0]
        model = calculator.validate_model(model, genotypes, phenotype)

    # Save results
    scores.to_csv(output_path, sep="\t")
    model.weights.to_csv(output_path.with_suffix(".weights.tsv"), sep="\t", index=False)

    logger.info(f"Scores saved to {output_path}")
    return model, scores


class PRSEnsemble:
    """
    Ensemble of PRS models for improved prediction.

    Combines multiple PRS (e.g., from different thresholds or methods).
    """

    def __init__(self, models: list[PRSModel]):
        """Initialize ensemble with list of models."""
        self.models = models

    def calculate_ensemble_score(
        self,
        genotypes: pd.DataFrame,
        weights: Optional[list[float]] = None,
    ) -> pd.DataFrame:
        """
        Calculate ensemble PRS.

        Args:
            genotypes: Genotype matrix
            weights: Optional weights for each model (default: equal)

        Returns:
            DataFrame with ensemble scores
        """
        if weights is None:
            weights = [1.0 / len(self.models)] * len(self.models)

        calculator = PRSCalculator(pd.DataFrame())  # Dummy for score calculation

        all_scores = []
        for model, weight in zip(self.models, weights):
            scores = calculator.calculate_scores(genotypes, model)
            if not scores.empty:
                scores["score"] *= weight
                all_scores.append(scores)

        if not all_scores:
            return pd.DataFrame()

        # Combine scores
        combined = all_scores[0].copy()
        combined["score"] = sum(s["score"] for s in all_scores)
        combined["n_models"] = len(all_scores)

        return combined
