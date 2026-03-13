#!/usr/bin/env python3
"""
Fine-mapping integration module for AlphaGWAS pipeline.

Integrates statistical fine-mapping methods:
- SuSiE (Sum of Single Effects)
- FINEMAP
- CARMA (if available)

Combines fine-mapping posterior probabilities with AlphaGenome predictions
for enhanced variant prioritization.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import subprocess
import tempfile
import shutil
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import yaml

# Import utilities
try:
    from .utils import retry_with_backoff, progress_iterator, ProgressTracker, print_summary_table
except ImportError:
    from utils import retry_with_backoff, progress_iterator, ProgressTracker, print_summary_table

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass
class FinemapResult:
    """Container for fine-mapping results."""
    variant_id: str
    pip: float  # Posterior Inclusion Probability
    credible_set: int  # Credible set membership (0 if not in any)
    method: str
    log10bf: Optional[float] = None  # Log10 Bayes Factor


class SuSiERunner:
    """
    Interface for running SuSiE fine-mapping.

    SuSiE (Sum of Single Effects) performs Bayesian fine-mapping
    to identify likely causal variants at GWAS loci.

    Requires: R with susieR package installed
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize SuSiE runner.

        Args:
            config: Configuration dictionary
        """
        self.config = config.get('finemapping', {}).get('susie', {})
        self.L = self.config.get('L', 10)  # Max number of causal variants
        self.coverage = self.config.get('coverage', 0.95)  # Credible set coverage
        self.min_abs_corr = self.config.get('min_abs_corr', 0.5)

        self._check_r_available()

    def _check_r_available(self) -> bool:
        """Check if R and susieR are available."""
        try:
            result = subprocess.run(
                ['Rscript', '-e', 'library(susieR); cat("OK")'],
                capture_output=True,
                text=True,
                timeout=30
            )
            if 'OK' in result.stdout:
                logger.info("SuSiE (R) is available")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        logger.warning("R with susieR not found. SuSiE fine-mapping will be simulated.")
        return False

    def run(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n_samples: int,
        variant_ids: List[str]
    ) -> List[FinemapResult]:
        """
        Run SuSiE fine-mapping.

        Args:
            z_scores: Array of Z-scores for each variant
            ld_matrix: LD correlation matrix (R, not R^2)
            n_samples: Sample size from GWAS
            variant_ids: List of variant identifiers

        Returns:
            List of FinemapResult objects
        """
        if len(z_scores) != len(variant_ids):
            raise ValueError("z_scores and variant_ids must have same length")

        if ld_matrix.shape != (len(z_scores), len(z_scores)):
            raise ValueError("LD matrix dimensions must match number of variants")

        # Try to run actual SuSiE
        try:
            return self._run_susie_r(z_scores, ld_matrix, n_samples, variant_ids)
        except Exception as e:
            logger.warning(f"SuSiE R execution failed: {e}")
            logger.info("Using simulated fine-mapping results")
            return self._simulate_finemapping(z_scores, variant_ids)

    def _run_susie_r(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n_samples: int,
        variant_ids: List[str]
    ) -> List[FinemapResult]:
        """Run SuSiE via R."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save inputs
            z_file = Path(tmpdir) / 'z_scores.txt'
            ld_file = Path(tmpdir) / 'ld_matrix.txt'
            out_file = Path(tmpdir) / 'susie_results.txt'

            np.savetxt(z_file, z_scores)
            np.savetxt(ld_file, ld_matrix)

            # R script
            r_script = f"""
            library(susieR)

            z <- scan("{z_file}")
            R <- as.matrix(read.table("{ld_file}"))
            n <- {n_samples}

            # Run SuSiE
            fit <- susie_rss(
                z = z,
                R = R,
                n = n,
                L = {self.L},
                coverage = {self.coverage},
                min_abs_corr = {self.min_abs_corr}
            )

            # Extract PIPs
            pip <- fit$pip

            # Get credible sets
            cs <- rep(0, length(pip))
            if (!is.null(fit$sets$cs)) {{
                for (i in seq_along(fit$sets$cs)) {{
                    cs[fit$sets$cs[[i]]] <- i
                }}
            }}

            # Save results
            results <- data.frame(pip = pip, cs = cs)
            write.table(results, "{out_file}", row.names = FALSE, quote = FALSE)
            """

            script_file = Path(tmpdir) / 'run_susie.R'
            with open(script_file, 'w') as f:
                f.write(r_script)

            # Execute R script
            result = subprocess.run(
                ['Rscript', str(script_file)],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                raise RuntimeError(f"SuSiE failed: {result.stderr}")

            # Read results
            results_df = pd.read_csv(out_file, sep=' ')

            return [
                FinemapResult(
                    variant_id=var_id,
                    pip=float(results_df.loc[i, 'pip']),
                    credible_set=int(results_df.loc[i, 'cs']),
                    method='SuSiE'
                )
                for i, var_id in enumerate(variant_ids)
            ]

    def _simulate_finemapping(
        self,
        z_scores: np.ndarray,
        variant_ids: List[str]
    ) -> List[FinemapResult]:
        """
        Simulate fine-mapping results based on Z-scores.

        This is a fallback when R/SuSiE is not available.
        Uses approximate Bayes factors to estimate PIPs.
        """
        # Calculate approximate Bayes factors (Wakefield approximation)
        # ABF = sqrt(1 + W/V) * exp(z^2 * W / (2 * (W + V)))
        # where W is prior variance on effect size, V = 1 (for z-scores)
        W = 0.04  # Prior variance (assumes effect size ~ N(0, 0.2^2))

        log_abf = 0.5 * np.log(1 / (1 + W)) + (z_scores ** 2) * W / (2 * (W + 1))
        abf = np.exp(log_abf - np.max(log_abf))  # Normalize for numerical stability

        # Calculate PIPs
        pip = abf / np.sum(abf)

        # Determine credible set
        sorted_idx = np.argsort(pip)[::-1]
        cumsum = 0
        credible_set = np.zeros(len(pip), dtype=int)

        for rank, idx in enumerate(sorted_idx):
            cumsum += pip[idx]
            credible_set[idx] = 1  # In 95% credible set
            if cumsum >= 0.95:
                break

        return [
            FinemapResult(
                variant_id=var_id,
                pip=float(pip[i]),
                credible_set=int(credible_set[i]),
                method='SuSiE_approx',
                log10bf=float(log_abf[i] / np.log(10))
            )
            for i, var_id in enumerate(variant_ids)
        ]


class FINEMAPRunner:
    """
    Interface for running FINEMAP.

    FINEMAP performs shotgun stochastic search for fine-mapping.

    Requires: FINEMAP executable in PATH
    """

    def __init__(self, config: Dict[str, Any]):
        """Initialize FINEMAP runner."""
        self.config = config.get('finemapping', {}).get('finemap', {})
        self.n_causal = self.config.get('n_causal', 5)
        self.n_iter = self.config.get('n_iter', 100000)

        self.finemap_path = self._find_finemap()

    def _find_finemap(self) -> Optional[str]:
        """Find FINEMAP executable."""
        result = shutil.which('finemap')
        if result:
            logger.info(f"FINEMAP found at: {result}")
            return result

        logger.warning("FINEMAP executable not found in PATH")
        return None

    def run(
        self,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n_samples: int,
        variant_ids: List[str],
        mafs: Optional[np.ndarray] = None
    ) -> List[FinemapResult]:
        """
        Run FINEMAP.

        Args:
            z_scores: Array of Z-scores
            ld_matrix: LD correlation matrix
            n_samples: Sample size
            variant_ids: Variant identifiers
            mafs: Minor allele frequencies (optional)

        Returns:
            List of FinemapResult objects
        """
        if self.finemap_path is None:
            logger.info("Using simulated FINEMAP results")
            return self._simulate_finemap(z_scores, variant_ids)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Prepare input files
            self._write_finemap_inputs(
                tmpdir, z_scores, ld_matrix, n_samples, variant_ids, mafs
            )

            # Run FINEMAP
            master_file = tmpdir / 'master'
            result = subprocess.run(
                [self.finemap_path, '--sss', '--in-files', str(master_file)],
                capture_output=True,
                text=True,
                timeout=600
            )

            if result.returncode != 0:
                logger.warning(f"FINEMAP failed: {result.stderr}")
                return self._simulate_finemap(z_scores, variant_ids)

            # Parse results
            return self._parse_finemap_results(tmpdir, variant_ids)

    def _write_finemap_inputs(
        self,
        tmpdir: Path,
        z_scores: np.ndarray,
        ld_matrix: np.ndarray,
        n_samples: int,
        variant_ids: List[str],
        mafs: Optional[np.ndarray]
    ):
        """Write FINEMAP input files."""
        # Z-scores file
        z_df = pd.DataFrame({
            'rsid': variant_ids,
            'chromosome': 1,
            'position': range(len(variant_ids)),
            'allele1': 'A',
            'allele2': 'G',
            'maf': mafs if mafs is not None else 0.1,
            'beta': z_scores / np.sqrt(n_samples),
            'se': 1 / np.sqrt(n_samples)
        })
        z_df.to_csv(tmpdir / 'data.z', sep=' ', index=False)

        # LD matrix
        np.savetxt(tmpdir / 'data.ld', ld_matrix, fmt='%.6f')

        # Master file
        with open(tmpdir / 'master', 'w') as f:
            f.write('z;ld;snp;config;cred;log;n_samples\n')
            f.write(f'data.z;data.ld;data.snp;data.config;data.cred;data.log;{n_samples}\n')

    def _parse_finemap_results(
        self,
        tmpdir: Path,
        variant_ids: List[str]
    ) -> List[FinemapResult]:
        """Parse FINEMAP output."""
        snp_file = tmpdir / 'data.snp'

        if not snp_file.exists():
            return self._simulate_finemap(np.zeros(len(variant_ids)), variant_ids)

        results_df = pd.read_csv(snp_file, sep=' ')

        # Map PIPs back to original variant order
        pip_map = dict(zip(results_df['rsid'], results_df['prob']))

        return [
            FinemapResult(
                variant_id=var_id,
                pip=pip_map.get(var_id, 0.0),
                credible_set=1 if pip_map.get(var_id, 0) > 0.01 else 0,
                method='FINEMAP',
                log10bf=None
            )
            for var_id in variant_ids
        ]

    def _simulate_finemap(
        self,
        z_scores: np.ndarray,
        variant_ids: List[str]
    ) -> List[FinemapResult]:
        """Simulate FINEMAP results."""
        # Use same approximation as SuSiE simulation
        W = 0.04
        log_abf = 0.5 * np.log(1 / (1 + W)) + (z_scores ** 2) * W / (2 * (W + 1))
        abf = np.exp(log_abf - np.max(log_abf))
        pip = abf / np.sum(abf)

        return [
            FinemapResult(
                variant_id=var_id,
                pip=float(pip[i]),
                credible_set=1 if pip[i] > 0.01 else 0,
                method='FINEMAP_approx'
            )
            for i, var_id in enumerate(variant_ids)
        ]


class FinemappingIntegrator:
    """
    Integrates fine-mapping results with AlphaGenome predictions.

    Combines statistical fine-mapping PIPs with functional predictions
    for enhanced variant prioritization.
    """

    def __init__(self, config: Dict[str, Any]):
        """Initialize integrator."""
        self.config = config.get('finemapping', {})
        self.method = self.config.get('method', 'susie')
        self.pip_weight = self.config.get('pip_weight', 0.5)

        if self.method == 'susie':
            self.runner = SuSiERunner(config)
        elif self.method == 'finemap':
            self.runner = FINEMAPRunner(config)
        else:
            self.runner = SuSiERunner(config)  # Default

    def calculate_z_scores(
        self,
        variants: pd.DataFrame,
        beta_col: str = 'beta',
        se_col: str = 'se',
        pval_col: str = 'pvalue'
    ) -> np.ndarray:
        """
        Calculate Z-scores from GWAS summary statistics.

        Args:
            variants: DataFrame with effect estimates
            beta_col: Column name for effect size
            se_col: Column name for standard error
            pval_col: Column name for p-value

        Returns:
            Array of Z-scores
        """
        if beta_col in variants.columns and se_col in variants.columns:
            return (variants[beta_col] / variants[se_col]).values

        elif pval_col in variants.columns:
            # Convert p-values to Z-scores
            from scipy import stats
            pvals = variants[pval_col].values
            # Handle edge cases
            pvals = np.clip(pvals, 1e-300, 1 - 1e-10)
            z = stats.norm.ppf(1 - pvals / 2)
            return z

        else:
            raise ValueError("Need beta/se or pvalue columns to calculate Z-scores")

    def estimate_ld_matrix(
        self,
        variants: pd.DataFrame,
        ld_file: Optional[str] = None
    ) -> np.ndarray:
        """
        Estimate or load LD matrix.

        Args:
            variants: DataFrame with variant info
            ld_file: Optional path to pre-computed LD matrix

        Returns:
            LD correlation matrix
        """
        n = len(variants)

        if ld_file and Path(ld_file).exists():
            ld_matrix = np.loadtxt(ld_file)
            if ld_matrix.shape == (n, n):
                return ld_matrix

        # Generate approximate LD based on distance
        # This is a rough approximation - real LD should be computed from genotypes
        logger.warning("Using distance-based LD approximation. "
                      "For accurate results, provide pre-computed LD matrix.")

        positions = variants['position'].values
        ld_matrix = np.eye(n)

        for i in range(n):
            for j in range(i + 1, n):
                dist = abs(positions[i] - positions[j])
                # Approximate LD decay: r^2 decays with distance
                # Using exponential decay with ~50kb half-life
                r = np.exp(-dist / 50000)
                ld_matrix[i, j] = r
                ld_matrix[j, i] = r

        return ld_matrix

    def run_finemapping(
        self,
        variants: pd.DataFrame,
        n_samples: int,
        ld_matrix: Optional[np.ndarray] = None
    ) -> pd.DataFrame:
        """
        Run fine-mapping and return results.

        Args:
            variants: DataFrame with variant data
            n_samples: GWAS sample size
            ld_matrix: Optional pre-computed LD matrix

        Returns:
            DataFrame with fine-mapping results
        """
        logger.info(f"Running fine-mapping with {self.method} on {len(variants)} variants")

        # Calculate Z-scores
        z_scores = self.calculate_z_scores(variants)

        # Get or estimate LD matrix
        if ld_matrix is None:
            ld_matrix = self.estimate_ld_matrix(variants)

        # Get variant IDs
        if 'variant_id' in variants.columns:
            variant_ids = variants['variant_id'].tolist()
        elif 'rsid' in variants.columns:
            variant_ids = variants['rsid'].tolist()
        else:
            variant_ids = [f"var_{i}" for i in range(len(variants))]

        # Run fine-mapping
        results = self.runner.run(z_scores, ld_matrix, n_samples, variant_ids)

        # Convert to DataFrame
        results_df = pd.DataFrame([
            {
                'variant_id': r.variant_id,
                'pip': r.pip,
                'credible_set': r.credible_set,
                'finemapping_method': r.method,
                'log10bf': r.log10bf
            }
            for r in results
        ])

        return results_df

    def integrate_with_predictions(
        self,
        ranked_variants: pd.DataFrame,
        finemapping_results: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Integrate fine-mapping PIPs with AlphaGenome predictions.

        Args:
            ranked_variants: DataFrame with AlphaGenome-based rankings
            finemapping_results: DataFrame with fine-mapping PIPs

        Returns:
            DataFrame with integrated scores
        """
        # Merge results
        merged = ranked_variants.merge(
            finemapping_results,
            on='variant_id',
            how='left'
        )

        # Fill missing PIPs with low value
        merged['pip'] = merged['pip'].fillna(0.001)

        # Calculate integrated score
        # Combines AlphaGenome consensus score with fine-mapping PIP
        alphagwas_weight = 1 - self.pip_weight

        if 'consensus_score' in merged.columns:
            # Normalize scores to 0-1 range
            max_score = merged['consensus_score'].max()
            if max_score > 0:
                norm_score = merged['consensus_score'] / max_score
            else:
                norm_score = 0

            merged['integrated_score'] = (
                alphagwas_weight * norm_score +
                self.pip_weight * merged['pip']
            )
        else:
            merged['integrated_score'] = merged['pip']

        # Re-rank by integrated score
        merged['integrated_rank'] = merged['integrated_score'].rank(
            ascending=False, method='min'
        ).astype(int)

        # Sort by integrated rank
        merged = merged.sort_values('integrated_rank')

        return merged


def main(config_path: str = "config/config.yaml"):
    """Main fine-mapping workflow."""
    config = load_config(config_path)

    prefix = config['output']['prefix']
    output_dir = Path(config['output']['dir'])
    intermediate_dir = Path("data/intermediate")

    # Load ranked variants
    ranked_file = output_dir / f"{prefix}_ranked_variants.tsv"
    if not ranked_file.exists():
        logger.error(f"Ranked variants file not found: {ranked_file}")
        logger.error("Run the main pipeline first.")
        return

    ranked_variants = pd.read_csv(ranked_file, sep='\t')
    logger.info(f"Loaded {len(ranked_variants)} ranked variants")

    # Load original GWAS data for Z-scores
    gwas_files = [
        intermediate_dir / f"{prefix}_significant_variants.tsv",
        intermediate_dir / f"{prefix}_variants_with_ld.tsv"
    ]

    gwas_data = None
    for f in gwas_files:
        if f.exists():
            gwas_data = pd.read_csv(f, sep='\t')
            break

    if gwas_data is None:
        logger.warning("Original GWAS data not found. Using ranked variants for Z-score estimation.")
        gwas_data = ranked_variants

    # Get sample size from config or estimate
    n_samples = config.get('gwas', {}).get('n_samples', 100000)

    # Initialize integrator
    integrator = FinemappingIntegrator(config)

    # Run fine-mapping
    finemapping_results = integrator.run_finemapping(gwas_data, n_samples)

    # Save fine-mapping results
    fm_file = output_dir / f"{prefix}_finemapping_results.tsv"
    finemapping_results.to_csv(fm_file, sep='\t', index=False)
    logger.info(f"Saved fine-mapping results to {fm_file}")

    # Integrate with AlphaGenome predictions
    integrated = integrator.integrate_with_predictions(ranked_variants, finemapping_results)

    # Save integrated results
    int_file = output_dir / f"{prefix}_integrated_ranking.tsv"
    integrated.to_csv(int_file, sep='\t', index=False)
    logger.info(f"Saved integrated ranking to {int_file}")

    # Print summary
    print_summary_table({
        'Variants analyzed': len(finemapping_results),
        'In 95% credible set': (finemapping_results['credible_set'] > 0).sum(),
        'PIP > 0.5': (finemapping_results['pip'] > 0.5).sum(),
        'PIP > 0.1': (finemapping_results['pip'] > 0.1).sum(),
        'Method': finemapping_results['finemapping_method'].iloc[0]
    }, title="Fine-mapping Summary")

    return integrated


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run statistical fine-mapping")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    parser.add_argument("--method", choices=['susie', 'finemap'], default='susie',
                       help="Fine-mapping method to use")
    args = parser.parse_args()

    main(args.config)
