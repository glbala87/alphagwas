#!/usr/bin/env python3
"""
Step 5: Score and prioritize variants based on AlphaGenome predictions.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import yaml
from typing import Dict, List, Optional
from scipy import stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


class VariantScorer:
    """Score and rank variants based on predicted functional effects."""

    def __init__(self, config: dict):
        self.config = config
        self.scoring_config = config.get('scoring', {})
        self.consensus_method = self.scoring_config.get('consensus_method', 'mean')
        self.min_effect = self.scoring_config.get('min_effect_threshold', 0.1)

    def calculate_tissue_scores(self, predictions: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate aggregate scores per variant per tissue.

        Combines predictions across modalities within each tissue.
        """
        # Group by variant and tissue
        grouped = predictions.groupby(['variant_id', 'tissue'])

        tissue_scores = grouped.agg({
            'effect_size': ['mean', 'max', 'std'],
            'confidence': 'mean',
            'prediction': lambda x: np.abs(x).max()
        }).reset_index()

        # Flatten column names
        tissue_scores.columns = [
            'variant_id', 'tissue',
            'mean_effect', 'max_effect', 'effect_std',
            'mean_confidence', 'max_abs_prediction'
        ]

        # Calculate composite tissue score
        tissue_scores['tissue_score'] = (
            tissue_scores['max_effect'] *
            tissue_scores['mean_confidence'] *
            (1 + tissue_scores['max_abs_prediction'])
        )

        return tissue_scores

    def calculate_consensus_scores(self, tissue_scores: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate consensus score across all tissues.
        """
        grouped = tissue_scores.groupby('variant_id')

        if self.consensus_method == 'mean':
            agg_func = 'mean'
        elif self.consensus_method == 'median':
            agg_func = 'median'
        elif self.consensus_method == 'max':
            agg_func = 'max'
        else:
            agg_func = 'mean'

        consensus = grouped.agg({
            'tissue_score': [agg_func, 'max', 'std'],
            'mean_effect': 'mean',
            'max_effect': 'max',
            'mean_confidence': 'mean'
        }).reset_index()

        consensus.columns = [
            'variant_id',
            'consensus_score', 'best_tissue_score', 'tissue_score_std',
            'avg_effect', 'max_effect', 'avg_confidence'
        ]

        # Add number of tissues with significant effects
        significant_tissues = tissue_scores[tissue_scores['max_effect'] >= self.min_effect]
        tissue_counts = significant_tissues.groupby('variant_id').size().reset_index(name='n_significant_tissues')
        consensus = consensus.merge(tissue_counts, on='variant_id', how='left')
        consensus['n_significant_tissues'] = consensus['n_significant_tissues'].fillna(0).astype(int)

        return consensus

    def rank_variants(self, consensus_scores: pd.DataFrame) -> pd.DataFrame:
        """
        Rank variants by consensus score and other metrics.
        """
        df = consensus_scores.copy()

        # Primary rank by consensus score
        df['rank_consensus'] = df['consensus_score'].rank(ascending=False, method='min')

        # Secondary rank by max effect across tissues
        df['rank_max_effect'] = df['max_effect'].rank(ascending=False, method='min')

        # Combined rank (average of ranks)
        df['rank_combined'] = (df['rank_consensus'] + df['rank_max_effect']) / 2

        # Final rank
        df['final_rank'] = df['rank_combined'].rank(method='min').astype(int)

        # Sort by final rank
        df = df.sort_values('final_rank')

        return df

    def identify_top_tissues(
        self,
        tissue_scores: pd.DataFrame,
        top_n: int = 3
    ) -> pd.DataFrame:
        """Identify top tissues for each variant."""
        top_tissues = []

        for variant_id in tissue_scores['variant_id'].unique():
            variant_data = tissue_scores[tissue_scores['variant_id'] == variant_id]
            top = variant_data.nlargest(top_n, 'tissue_score')

            tissues = top['tissue'].tolist()
            scores = top['tissue_score'].tolist()

            top_tissues.append({
                'variant_id': variant_id,
                'top_tissues': ','.join(tissues),
                'top_tissue_scores': ','.join([f"{s:.3f}" for s in scores])
            })

        return pd.DataFrame(top_tissues)


def annotate_with_genes(variants: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    """Add gene annotations from predictions."""
    if 'nearby_genes' not in predictions.columns:
        return variants

    def extract_genes(gene_series):
        """Extract unique genes from a series of gene lists/strings."""
        all_genes = set()
        for genes in gene_series:
            if genes is None:
                continue
            if isinstance(genes, list):
                all_genes.update(genes)
            elif isinstance(genes, str):
                all_genes.add(genes)
        return ','.join(sorted(all_genes)) if all_genes else ''

    # Get unique genes per variant
    gene_info = predictions.groupby('variant_id')['nearby_genes'].apply(extract_genes).reset_index()
    gene_info.columns = ['variant_id', 'nearby_genes']

    return variants.merge(gene_info, on='variant_id', how='left')


def create_summary_table(
    ranked_variants: pd.DataFrame,
    tissue_scores: pd.DataFrame,
    original_variants: pd.DataFrame
) -> pd.DataFrame:
    """Create final summary table with all relevant information."""

    # Get top tissues per variant
    scorer = VariantScorer({})
    top_tissues = scorer.identify_top_tissues(tissue_scores)

    # Merge all information
    summary = ranked_variants.merge(top_tissues, on='variant_id', how='left')

    # Add original variant info (rsid, position, etc.)
    if 'rsid' in original_variants.columns:
        variant_info = original_variants[['rsid', 'chromosome', 'position']].drop_duplicates()
        variant_info['variant_id'] = 'chr' + variant_info['chromosome'].astype(str) + ':' + variant_info['position'].astype(str)
        summary = summary.merge(variant_info, on='variant_id', how='left')

    # Reorder columns
    priority_cols = ['final_rank', 'rsid', 'variant_id', 'chromosome', 'position',
                     'consensus_score', 'max_effect', 'n_significant_tissues',
                     'top_tissues', 'nearby_genes']
    other_cols = [c for c in summary.columns if c not in priority_cols]
    available_cols = [c for c in priority_cols if c in summary.columns] + other_cols
    summary = summary[available_cols]

    return summary


def main(config_path: str = "config/config.yaml"):
    """Main scoring and prioritization workflow."""
    config = load_config(config_path)

    # Load predictions from previous step
    prefix = config['output']['prefix']
    intermediate_dir = Path("data/intermediate")
    output_dir = Path(config['output']['dir'])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load AlphaGenome predictions
    pred_file = intermediate_dir / f"{prefix}_alphagenome_predictions.parquet"
    if pred_file.exists():
        predictions = pd.read_parquet(pred_file)
    else:
        pred_file = intermediate_dir / f"{prefix}_alphagenome_predictions.tsv"
        if not pred_file.exists():
            logger.error(f"Predictions file not found: {pred_file}")
            logger.error("Run 04_alphagenome_predict.py first")
            return
        predictions = pd.read_csv(pred_file, sep='\t')

    logger.info(f"Loaded {len(predictions)} predictions")

    # Load original variants for annotation
    variant_files = [
        intermediate_dir / f"{prefix}_variants_hg38.tsv",
        intermediate_dir / f"{prefix}_variants_with_ld.tsv",
        intermediate_dir / f"{prefix}_significant_variants.tsv"
    ]
    original_variants = None
    for vf in variant_files:
        if vf.exists():
            original_variants = pd.read_csv(vf, sep='\t')
            break

    # Initialize scorer
    scorer = VariantScorer(config)

    # Calculate tissue-specific scores
    logger.info("Calculating tissue-specific scores...")
    tissue_scores = scorer.calculate_tissue_scores(predictions)

    # Calculate consensus scores
    logger.info("Calculating consensus scores...")
    consensus_scores = scorer.calculate_consensus_scores(tissue_scores)

    # Rank variants
    logger.info("Ranking variants...")
    ranked_variants = scorer.rank_variants(consensus_scores)

    # Add gene annotations
    ranked_variants = annotate_with_genes(ranked_variants, predictions)

    # Create summary table
    if original_variants is not None:
        summary = create_summary_table(ranked_variants, tissue_scores, original_variants)
    else:
        summary = ranked_variants

    # Save results
    # 1. Full ranked list
    ranked_file = output_dir / f"{prefix}_ranked_variants.tsv"
    summary.to_csv(ranked_file, sep='\t', index=False)
    logger.info(f"Saved ranked variants to {ranked_file}")

    # 2. Tissue-specific scores
    tissue_file = output_dir / f"{prefix}_tissue_scores.tsv"
    tissue_scores.to_csv(tissue_file, sep='\t', index=False)
    logger.info(f"Saved tissue scores to {tissue_file}")

    # 3. Top variants (top 20)
    top_file = output_dir / f"{prefix}_top20_variants.tsv"
    summary.head(20).to_csv(top_file, sep='\t', index=False)
    logger.info(f"Saved top 20 variants to {top_file}")

    # Print summary
    print("\n" + "="*60)
    print("VARIANT PRIORITIZATION SUMMARY")
    print("="*60)
    print(f"Total variants scored: {len(summary)}")
    print(f"Variants with significant effects: {len(summary[summary['max_effect'] >= scorer.min_effect])}")
    print(f"\nTop 10 prioritized variants:")
    print("-"*60)

    display_cols = ['final_rank', 'rsid', 'consensus_score', 'max_effect', 'top_tissues']
    display_cols = [c for c in display_cols if c in summary.columns]
    if display_cols:
        print(summary[display_cols].head(10).to_string(index=False))
    print("="*60 + "\n")

    return summary


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Score and prioritize variants")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
