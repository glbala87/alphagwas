#!/usr/bin/env python3
"""
Truthset Validation for AlphaGWAS Pipeline.

Compare pipeline results against known causal variants (truthset) to evaluate
discovery performance. Supports matching by rsid, genomic position, or
locus window, and reports precision, recall, F1, sensitivity, and specificity.

Truthset sources:
- NHGRI-EBI GWAS Catalog exports
- Simulated datasets with known causal variants
- Published replication sets
- OpenTargets or other curated databases

Usage:
    python validate_truthset.py --config config.yaml
    python validate_truthset.py --results ranked.tsv --truthset truth.tsv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import yaml
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class MatchResult:
    """Result of matching a single truthset variant against pipeline results."""
    truth_id: str
    matched: bool
    match_type: Optional[str] = None  # exact_rsid, exact_position, window
    matched_variant_id: Optional[str] = None
    distance_bp: Optional[int] = None
    pipeline_rank: Optional[int] = None
    pipeline_score: Optional[float] = None


@dataclass
class ValidationMetrics:
    """Aggregate validation metrics."""
    n_truth: int = 0
    n_results: int = 0
    true_positives: int = 0
    false_positives: int = 0
    false_negatives: int = 0
    true_negatives: int = 0
    precision: float = 0.0
    recall: float = 0.0
    f1_score: float = 0.0
    sensitivity: float = 0.0
    specificity: float = 0.0
    mean_rank_of_true: float = 0.0
    median_rank_of_true: float = 0.0
    fraction_in_top10: float = 0.0
    fraction_in_top20: float = 0.0
    fraction_in_top50: float = 0.0
    match_details: List[MatchResult] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Truthset loading
# ---------------------------------------------------------------------------

def load_truthset(
    path: str,
    columns: Optional[Dict[str, str]] = None,
) -> pd.DataFrame:
    """
    Load a truthset file.

    Supports TSV/CSV with flexible column names. Returns a DataFrame with
    normalised columns: chromosome, position, rsid, gene, trait, source.

    Args:
        path: Path to truthset file.
        columns: Optional mapping of file column names to standard names.
                 e.g. {"CHR": "chromosome", "BP": "position"}

    Returns:
        DataFrame with standardised truthset columns.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Truthset file not found: {path}")

    sep = '\t' if path.suffix in ('.tsv', '.txt') else ','
    df = pd.read_csv(path, sep=sep)
    logger.info(f"Loaded truthset with {len(df)} variants from {path}")

    # Apply column mapping
    if columns:
        df = df.rename(columns=columns)

    # Auto-detect common column name variants
    col_aliases = {
        'chromosome': ['chr', 'chrom', '#chrom', 'CHR', 'CHROM', 'Chr'],
        'position': ['pos', 'bp', 'BP', 'POS', 'base_pair_location'],
        'rsid': ['snp', 'SNP', 'rsID', 'RSID', 'variant_id', 'SNPS'],
        'gene': ['GENE', 'mapped_gene', 'MAPPED_GENE', 'gene_name'],
        'trait': ['TRAIT', 'phenotype', 'DISEASE/TRAIT', 'disease_trait'],
        'pvalue': ['p_value', 'P', 'PVALUE', 'P-VALUE', 'p-value'],
        'source': ['SOURCE', 'study', 'STUDY'],
    }

    for standard, aliases in col_aliases.items():
        if standard not in df.columns:
            for alias in aliases:
                if alias in df.columns:
                    df = df.rename(columns={alias: standard})
                    break

    # Normalise chromosome values
    if 'chromosome' in df.columns:
        df['chromosome'] = (
            df['chromosome'].astype(str)
            .str.replace('^chr', '', regex=True)
            .str.strip()
        )

    # Ensure position is integer where present
    if 'position' in df.columns:
        df['position'] = pd.to_numeric(df['position'], errors='coerce')
        df = df.dropna(subset=['position'])
        df['position'] = df['position'].astype(int)

    return df


def load_results(
    path: str,
    rank_col: str = 'final_rank',
    score_col: str = 'consensus_score',
) -> pd.DataFrame:
    """
    Load pipeline results (ranked variants).

    Args:
        path: Path to ranked variants TSV.
        rank_col: Column containing variant rank.
        score_col: Column containing variant score.

    Returns:
        DataFrame with variant_id, chromosome, position, rsid, rank, score.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Results file not found: {path}")

    sep = '\t' if path.suffix in ('.tsv', '.txt') else ','
    df = pd.read_csv(path, sep=sep)
    logger.info(f"Loaded {len(df)} ranked variants from {path}")

    # Parse chromosome and position from variant_id if not already present
    if 'chromosome' not in df.columns and 'variant_id' in df.columns:
        parsed = df['variant_id'].str.extract(r'chr(\w+):(\d+)')
        df['chromosome'] = parsed[0]
        df['position'] = parsed[1].astype(int)

    # Normalise chromosome
    if 'chromosome' in df.columns:
        df['chromosome'] = df['chromosome'].astype(str).str.replace('^chr', '', regex=True)

    # Ensure rank and score columns exist
    if rank_col not in df.columns:
        df[rank_col] = range(1, len(df) + 1)
    if score_col not in df.columns:
        df[score_col] = np.nan

    return df


# ---------------------------------------------------------------------------
# Matching strategies
# ---------------------------------------------------------------------------

def match_by_rsid(
    truth: pd.DataFrame,
    results: pd.DataFrame,
) -> List[MatchResult]:
    """Match truthset to results by rsid."""
    matches = []
    if 'rsid' not in truth.columns or 'rsid' not in results.columns:
        return matches

    result_rsids = set(results['rsid'].dropna().astype(str))

    for _, row in truth.iterrows():
        rsid = str(row.get('rsid', ''))
        if not rsid or rsid == 'nan':
            continue

        if rsid in result_rsids:
            hit = results[results['rsid'] == rsid].iloc[0]
            matches.append(MatchResult(
                truth_id=rsid,
                matched=True,
                match_type='exact_rsid',
                matched_variant_id=hit.get('variant_id', rsid),
                distance_bp=0,
                pipeline_rank=int(hit.get('final_rank', 0)),
                pipeline_score=float(hit.get('consensus_score', 0)),
            ))
        else:
            matches.append(MatchResult(
                truth_id=rsid,
                matched=False,
            ))

    return matches


def match_by_position(
    truth: pd.DataFrame,
    results: pd.DataFrame,
    window_bp: int = 0,
) -> List[MatchResult]:
    """
    Match truthset to results by genomic position.

    Args:
        truth: Truthset DataFrame (must have chromosome, position).
        results: Results DataFrame (must have chromosome, position).
        window_bp: Matching window in base pairs (0 = exact match).

    Returns:
        List of MatchResult objects.
    """
    matches = []
    required = {'chromosome', 'position'}
    if not required.issubset(truth.columns) or not required.issubset(results.columns):
        logger.warning("Position matching requires chromosome and position columns")
        return matches

    # Index results by chromosome for fast lookup
    results_by_chr: Dict[str, pd.DataFrame] = {
        chrom: group for chrom, group in results.groupby('chromosome')
    }

    for _, row in truth.iterrows():
        chrom = str(row['chromosome'])
        pos = int(row['position'])
        truth_id = str(row.get('rsid', f'{chrom}:{pos}'))

        if chrom not in results_by_chr:
            matches.append(MatchResult(truth_id=truth_id, matched=False))
            continue

        chr_results = results_by_chr[chrom]
        distances = (chr_results['position'] - pos).abs()
        min_dist = distances.min()

        if min_dist <= window_bp:
            idx = distances.idxmin()
            hit = chr_results.loc[idx]
            match_type = 'exact_position' if min_dist == 0 else 'window'
            matches.append(MatchResult(
                truth_id=truth_id,
                matched=True,
                match_type=match_type,
                matched_variant_id=hit.get('variant_id', f'chr{chrom}:{int(hit["position"])}'),
                distance_bp=int(min_dist),
                pipeline_rank=int(hit.get('final_rank', 0)),
                pipeline_score=float(hit.get('consensus_score', 0)),
            ))
        else:
            matches.append(MatchResult(truth_id=truth_id, matched=False))

    return matches


# ---------------------------------------------------------------------------
# Validation engine
# ---------------------------------------------------------------------------

class TruthsetValidator:
    """
    Validate pipeline results against a known truthset.

    Supports multiple matching strategies, calculates classification metrics
    and ranking quality metrics.
    """

    def __init__(
        self,
        match_strategy: str = 'auto',
        window_bp: int = 500_000,
        rank_col: str = 'final_rank',
        score_col: str = 'consensus_score',
        n_total_variants: Optional[int] = None,
    ):
        """
        Args:
            match_strategy: 'rsid', 'position', 'window', or 'auto'.
            window_bp: Window size for window-based matching.
            rank_col: Column name for variant rank in results.
            score_col: Column name for variant score in results.
            n_total_variants: Total variants tested (for specificity). If None,
                              estimated from results size.
        """
        self.match_strategy = match_strategy
        self.window_bp = window_bp
        self.rank_col = rank_col
        self.score_col = score_col
        self.n_total_variants = n_total_variants

    def validate(
        self,
        truth: pd.DataFrame,
        results: pd.DataFrame,
    ) -> ValidationMetrics:
        """
        Run truthset validation.

        Args:
            truth: Truthset DataFrame.
            results: Ranked results DataFrame.

        Returns:
            ValidationMetrics with all computed metrics.
        """
        # Determine matching strategy
        strategy = self._resolve_strategy(truth, results)
        logger.info(f"Using match strategy: {strategy}")

        # Perform matching
        if strategy == 'rsid':
            matches = match_by_rsid(truth, results)
        elif strategy == 'position':
            matches = match_by_position(truth, results, window_bp=0)
        elif strategy == 'window':
            matches = match_by_position(truth, results, window_bp=self.window_bp)
        elif strategy == 'auto':
            # Try rsid first, fall back to window
            matches = match_by_rsid(truth, results)
            rsid_matched = sum(1 for m in matches if m.matched)

            pos_matches = match_by_position(truth, results, window_bp=self.window_bp)
            pos_matched = sum(1 for m in pos_matches if m.matched)

            if pos_matched > rsid_matched:
                matches = pos_matches
                logger.info(f"Auto-selected position/window matching ({pos_matched} vs {rsid_matched} rsid matches)")
            else:
                logger.info(f"Auto-selected rsid matching ({rsid_matched} vs {pos_matched} position matches)")
        else:
            raise ValueError(f"Unknown match strategy: {strategy}")

        # Calculate metrics
        metrics = self._compute_metrics(matches, truth, results)
        return metrics

    def _resolve_strategy(self, truth: pd.DataFrame, results: pd.DataFrame) -> str:
        """Resolve 'auto' strategy to a concrete one."""
        if self.match_strategy != 'auto':
            return self.match_strategy

        has_rsid = 'rsid' in truth.columns and 'rsid' in results.columns
        has_pos = 'chromosome' in truth.columns and 'position' in truth.columns
        has_pos = has_pos and 'chromosome' in results.columns and 'position' in results.columns

        if has_rsid and has_pos:
            return 'auto'  # will try both
        elif has_rsid:
            return 'rsid'
        elif has_pos:
            return 'window'
        else:
            raise ValueError("Truthset must have either rsid or chromosome+position columns")

    def _compute_metrics(
        self,
        matches: List[MatchResult],
        truth: pd.DataFrame,
        results: pd.DataFrame,
    ) -> ValidationMetrics:
        """Compute classification and ranking metrics from match results."""
        n_truth = len(truth)
        n_results = len(results)
        n_total = self.n_total_variants or n_results

        tp = sum(1 for m in matches if m.matched)
        fn = sum(1 for m in matches if not m.matched)
        # Variants in results that don't match any truth variant
        matched_result_ids = {m.matched_variant_id for m in matches if m.matched}
        fp = n_results - len(matched_result_ids)
        tn = max(0, n_total - tp - fp - fn)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        sensitivity = recall  # same as recall / TPR
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        # Ranking quality for matched truth variants
        ranks = [m.pipeline_rank for m in matches if m.matched and m.pipeline_rank]
        mean_rank = float(np.mean(ranks)) if ranks else 0.0
        median_rank = float(np.median(ranks)) if ranks else 0.0
        frac_top10 = sum(1 for r in ranks if r <= 10) / n_truth if n_truth > 0 else 0.0
        frac_top20 = sum(1 for r in ranks if r <= 20) / n_truth if n_truth > 0 else 0.0
        frac_top50 = sum(1 for r in ranks if r <= 50) / n_truth if n_truth > 0 else 0.0

        return ValidationMetrics(
            n_truth=n_truth,
            n_results=n_results,
            true_positives=tp,
            false_positives=fp,
            false_negatives=fn,
            true_negatives=tn,
            precision=precision,
            recall=recall,
            f1_score=f1,
            sensitivity=sensitivity,
            specificity=specificity,
            mean_rank_of_true=mean_rank,
            median_rank_of_true=median_rank,
            fraction_in_top10=frac_top10,
            fraction_in_top20=frac_top20,
            fraction_in_top50=frac_top50,
            match_details=matches,
        )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def format_validation_report(metrics: ValidationMetrics) -> str:
    """Format validation metrics as a human-readable report."""
    lines = [
        "",
        "=" * 64,
        "  TRUTHSET VALIDATION REPORT",
        "=" * 64,
        "",
        "  Dataset Summary",
        "  " + "-" * 40,
        f"  Truth variants:          {metrics.n_truth}",
        f"  Pipeline results:        {metrics.n_results}",
        "",
        "  Classification Metrics",
        "  " + "-" * 40,
        f"  True positives (TP):     {metrics.true_positives}",
        f"  False positives (FP):    {metrics.false_positives}",
        f"  False negatives (FN):    {metrics.false_negatives}",
        f"  True negatives (TN):     {metrics.true_negatives}",
        "",
        f"  Precision:               {metrics.precision:.4f}",
        f"  Recall / Sensitivity:    {metrics.recall:.4f}",
        f"  F1 Score:                {metrics.f1_score:.4f}",
        f"  Specificity:             {metrics.specificity:.4f}",
        "",
        "  Ranking Quality",
        "  " + "-" * 40,
        f"  Mean rank of true hits:  {metrics.mean_rank_of_true:.1f}",
        f"  Median rank of true hits:{metrics.median_rank_of_true:.1f}",
        f"  Fraction in top 10:      {metrics.fraction_in_top10:.2%}",
        f"  Fraction in top 20:      {metrics.fraction_in_top20:.2%}",
        f"  Fraction in top 50:      {metrics.fraction_in_top50:.2%}",
        "",
    ]

    # Per-variant details
    matched = [m for m in metrics.match_details if m.matched]
    missed = [m for m in metrics.match_details if not m.matched]

    if matched:
        lines.append("  Matched Variants")
        lines.append("  " + "-" * 40)
        lines.append(f"  {'Variant':<20} {'Match Type':<16} {'Rank':>6} {'Score':>8} {'Dist(bp)':>10}")
        for m in sorted(matched, key=lambda x: x.pipeline_rank or 9999):
            lines.append(
                f"  {m.truth_id:<20} {m.match_type or '':<16} "
                f"{m.pipeline_rank or 0:>6} {m.pipeline_score or 0:>8.4f} "
                f"{m.distance_bp or 0:>10}"
            )
        lines.append("")

    if missed:
        lines.append("  Missed Variants (in truthset but not found)")
        lines.append("  " + "-" * 40)
        for m in missed:
            lines.append(f"  {m.truth_id}")
        lines.append("")

    lines.append("=" * 64)
    return "\n".join(lines)


def save_validation_results(
    metrics: ValidationMetrics,
    output_dir: str,
    prefix: str = "study",
) -> Dict[str, str]:
    """
    Save validation results to files.

    Returns:
        Dict mapping file type to file path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    saved = {}

    # 1. Metrics summary TSV
    metrics_path = output_dir / f"{prefix}_truthset_metrics.tsv"
    metrics_data = {
        'metric': [
            'n_truth', 'n_results', 'true_positives', 'false_positives',
            'false_negatives', 'true_negatives', 'precision', 'recall',
            'f1_score', 'sensitivity', 'specificity',
            'mean_rank_of_true', 'median_rank_of_true',
            'fraction_in_top10', 'fraction_in_top20', 'fraction_in_top50',
        ],
        'value': [
            metrics.n_truth, metrics.n_results, metrics.true_positives,
            metrics.false_positives, metrics.false_negatives,
            metrics.true_negatives, metrics.precision, metrics.recall,
            metrics.f1_score, metrics.sensitivity, metrics.specificity,
            metrics.mean_rank_of_true, metrics.median_rank_of_true,
            metrics.fraction_in_top10, metrics.fraction_in_top20,
            metrics.fraction_in_top50,
        ],
    }
    pd.DataFrame(metrics_data).to_csv(metrics_path, sep='\t', index=False)
    saved['metrics'] = str(metrics_path)
    logger.info(f"Saved metrics to {metrics_path}")

    # 2. Per-variant match details TSV
    details_path = output_dir / f"{prefix}_truthset_details.tsv"
    details_data = [{
        'truth_id': m.truth_id,
        'matched': m.matched,
        'match_type': m.match_type or '',
        'matched_variant_id': m.matched_variant_id or '',
        'distance_bp': m.distance_bp if m.distance_bp is not None else '',
        'pipeline_rank': m.pipeline_rank or '',
        'pipeline_score': m.pipeline_score or '',
    } for m in metrics.match_details]
    pd.DataFrame(details_data).to_csv(details_path, sep='\t', index=False)
    saved['details'] = str(details_path)
    logger.info(f"Saved match details to {details_path}")

    # 3. Human-readable report
    report_path = output_dir / f"{prefix}_truthset_report.txt"
    report_path.write_text(format_validation_report(metrics))
    saved['report'] = str(report_path)
    logger.info(f"Saved report to {report_path}")

    return saved


# ---------------------------------------------------------------------------
# Power analysis
# ---------------------------------------------------------------------------

def power_by_effect_size(
    truth: pd.DataFrame,
    results: pd.DataFrame,
    validator: TruthsetValidator,
    effect_col: str = 'beta',
    bins: int = 5,
) -> pd.DataFrame:
    """
    Calculate recall stratified by effect size bins.

    Args:
        truth: Truthset with an effect size column.
        results: Pipeline results.
        validator: Configured TruthsetValidator.
        effect_col: Column in truthset containing effect sizes.
        bins: Number of effect size bins.

    Returns:
        DataFrame with columns: bin_label, n_variants, n_found, recall.
    """
    if effect_col not in truth.columns:
        logger.warning(f"Effect column '{effect_col}' not in truthset; skipping power analysis")
        return pd.DataFrame()

    truth = truth.copy()
    truth['abs_effect'] = truth[effect_col].abs()
    truth['effect_bin'] = pd.qcut(truth['abs_effect'], q=bins, duplicates='drop')

    rows = []
    for bin_label, group in truth.groupby('effect_bin', observed=True):
        sub_metrics = validator.validate(group, results)
        rows.append({
            'effect_bin': str(bin_label),
            'n_variants': len(group),
            'n_found': sub_metrics.true_positives,
            'recall': sub_metrics.recall,
        })

    return pd.DataFrame(rows)


def power_by_maf(
    truth: pd.DataFrame,
    results: pd.DataFrame,
    validator: TruthsetValidator,
    maf_col: str = 'maf',
    bins: int = 5,
) -> pd.DataFrame:
    """
    Calculate recall stratified by minor allele frequency bins.

    Args:
        truth: Truthset with a MAF column.
        results: Pipeline results.
        validator: Configured TruthsetValidator.
        maf_col: Column in truthset containing MAF.
        bins: Number of MAF bins.

    Returns:
        DataFrame with columns: maf_bin, n_variants, n_found, recall.
    """
    if maf_col not in truth.columns:
        logger.warning(f"MAF column '{maf_col}' not in truthset; skipping MAF analysis")
        return pd.DataFrame()

    truth = truth.copy()
    truth['maf_bin'] = pd.qcut(truth[maf_col], q=bins, duplicates='drop')

    rows = []
    for bin_label, group in truth.groupby('maf_bin', observed=True):
        sub_metrics = validator.validate(group, results)
        rows.append({
            'maf_bin': str(bin_label),
            'n_variants': len(group),
            'n_found': sub_metrics.true_positives,
            'recall': sub_metrics.recall,
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Pipeline integration
# ---------------------------------------------------------------------------

def main(config_path: str = "config/config.yaml"):
    """Run truthset validation as a pipeline step."""
    with open(config_path) as f:
        config = yaml.safe_load(f)

    validation_config = config.get('validation', {})
    truthset_path = validation_config.get('truthset_file')

    if not truthset_path:
        logger.info("No truthset_file specified in config validation section; skipping validation")
        return None

    # Load truthset
    truth_columns = validation_config.get('columns', {})
    truth = load_truthset(truthset_path, columns=truth_columns)

    # Load results
    output_dir = config.get('output', {}).get('dir', 'data/output')
    prefix = config.get('output', {}).get('prefix', 'study')
    results_path = Path(output_dir) / f"{prefix}_ranked_variants.tsv"
    results = load_results(str(results_path))

    # Configure validator
    match_strategy = validation_config.get('match_strategy', 'auto')
    window_bp = validation_config.get('window_bp', 500_000)
    n_total = validation_config.get('n_total_variants')

    validator = TruthsetValidator(
        match_strategy=match_strategy,
        window_bp=window_bp,
        n_total_variants=n_total,
    )

    # Validate
    metrics = validator.validate(truth, results)

    # Print report
    report = format_validation_report(metrics)
    print(report)

    # Save results
    saved = save_validation_results(metrics, output_dir, prefix)

    # Power analysis (if effect sizes available)
    effect_col = validation_config.get('effect_column', 'beta')
    power_df = power_by_effect_size(truth, results, validator, effect_col=effect_col)
    if not power_df.empty:
        power_path = Path(output_dir) / f"{prefix}_truthset_power.tsv"
        power_df.to_csv(power_path, sep='\t', index=False)
        logger.info(f"Saved power analysis to {power_path}")
        print("\n  Power by Effect Size")
        print("  " + "-" * 40)
        for _, row in power_df.iterrows():
            print(f"  {row['effect_bin']:<20} n={int(row['n_variants']):>4}  "
                  f"found={int(row['n_found']):>4}  recall={row['recall']:.2%}")

    # MAF analysis
    maf_col = validation_config.get('maf_column', 'maf')
    maf_df = power_by_maf(truth, results, validator, maf_col=maf_col)
    if not maf_df.empty:
        maf_path = Path(output_dir) / f"{prefix}_truthset_maf.tsv"
        maf_df.to_csv(maf_path, sep='\t', index=False)
        logger.info(f"Saved MAF analysis to {maf_path}")

    return metrics


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate pipeline results against truthset")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    parser.add_argument("--results", help="Path to ranked variants (overrides config)")
    parser.add_argument("--truthset", help="Path to truthset file (overrides config)")
    parser.add_argument("--strategy", choices=['rsid', 'position', 'window', 'auto'],
                        default='auto', help="Matching strategy")
    parser.add_argument("--window", type=int, default=500_000,
                        help="Window size in bp for window matching (default: 500000)")
    args = parser.parse_args()

    if args.results and args.truthset:
        # Standalone mode: skip config
        truth = load_truthset(args.truthset)
        results = load_results(args.results)
        validator = TruthsetValidator(
            match_strategy=args.strategy,
            window_bp=args.window,
        )
        metrics = validator.validate(truth, results)
        print(format_validation_report(metrics))
        save_validation_results(metrics, "data/output", "standalone")
    else:
        main(args.config)
