"""
Tests for score_variants module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.score_variants import (
    VariantScorer,
    annotate_with_genes,
    create_summary_table
)


class TestVariantScorer:
    """Tests for VariantScorer class."""

    def test_init_default(self):
        """Test default initialization."""
        scorer = VariantScorer({})

        assert scorer.consensus_method == 'mean'
        assert scorer.min_effect == 0.1

    def test_init_custom(self):
        """Test initialization with custom config."""
        config = {
            'scoring': {
                'consensus_method': 'max',
                'min_effect_threshold': 0.2
            }
        }
        scorer = VariantScorer(config)

        assert scorer.consensus_method == 'max'
        assert scorer.min_effect == 0.2

    def test_calculate_tissue_scores(self, sample_predictions):
        """Test tissue score calculation."""
        scorer = VariantScorer({})
        result = scorer.calculate_tissue_scores(sample_predictions)

        assert 'variant_id' in result.columns
        assert 'tissue' in result.columns
        assert 'tissue_score' in result.columns
        assert 'mean_effect' in result.columns
        assert 'max_effect' in result.columns

        # Check that we have scores for each variant-tissue combination
        unique_combos = result[['variant_id', 'tissue']].drop_duplicates()
        assert len(unique_combos) > 0

    def test_calculate_consensus_scores(self, sample_tissue_scores):
        """Test consensus score calculation."""
        scorer = VariantScorer({})
        result = scorer.calculate_consensus_scores(sample_tissue_scores)

        assert 'variant_id' in result.columns
        assert 'consensus_score' in result.columns
        assert 'best_tissue_score' in result.columns
        assert 'n_significant_tissues' in result.columns

        # Each variant should appear once
        assert len(result) == len(sample_tissue_scores['variant_id'].unique())

    def test_consensus_methods(self, sample_tissue_scores):
        """Test different consensus methods."""
        for method in ['mean', 'median', 'max']:
            config = {'scoring': {'consensus_method': method}}
            scorer = VariantScorer(config)
            result = scorer.calculate_consensus_scores(sample_tissue_scores)

            assert 'consensus_score' in result.columns
            assert not result['consensus_score'].isna().all()

    def test_rank_variants(self, sample_tissue_scores):
        """Test variant ranking."""
        scorer = VariantScorer({})
        consensus = scorer.calculate_consensus_scores(sample_tissue_scores)
        result = scorer.rank_variants(consensus)

        assert 'rank_consensus' in result.columns
        assert 'rank_max_effect' in result.columns
        assert 'rank_combined' in result.columns
        assert 'final_rank' in result.columns

        # Ranks should be sequential
        assert result['final_rank'].min() == 1
        assert result['final_rank'].max() == len(result)

        # Should be sorted by final rank
        assert list(result['final_rank']) == sorted(result['final_rank'])

    def test_identify_top_tissues(self, sample_tissue_scores):
        """Test identification of top tissues per variant."""
        scorer = VariantScorer({})
        result = scorer.identify_top_tissues(sample_tissue_scores, top_n=2)

        assert 'variant_id' in result.columns
        assert 'top_tissues' in result.columns
        assert 'top_tissue_scores' in result.columns

        # Each variant should have top tissues listed
        for _, row in result.iterrows():
            tissues = row['top_tissues'].split(',')
            assert len(tissues) <= 2


class TestAnnotateWithGenes:
    """Tests for gene annotation function."""

    def test_annotate_with_genes(self):
        """Test adding gene annotations."""
        variants = pd.DataFrame({
            'variant_id': ['chr1:1000000', 'chr1:1500000'],
            'score': [0.5, 0.3]
        })

        predictions = pd.DataFrame({
            'variant_id': ['chr1:1000000', 'chr1:1000000', 'chr1:1500000'],
            'nearby_genes': [['GENE1', 'GENE2'], ['GENE1'], 'GENE3']
        })

        result = annotate_with_genes(variants, predictions)

        assert 'nearby_genes' in result.columns
        # Check that genes are merged
        assert 'GENE1' in result[result['variant_id'] == 'chr1:1000000']['nearby_genes'].values[0]

    def test_annotate_no_genes_column(self):
        """Test when predictions don't have genes column."""
        variants = pd.DataFrame({
            'variant_id': ['chr1:1000000'],
            'score': [0.5]
        })

        predictions = pd.DataFrame({
            'variant_id': ['chr1:1000000'],
            'effect': [0.3]
        })

        result = annotate_with_genes(variants, predictions)

        # Should return variants unchanged
        assert len(result) == len(variants)
        assert 'nearby_genes' not in result.columns


class TestCreateSummaryTable:
    """Tests for summary table creation."""

    def test_create_summary(self, sample_tissue_scores):
        """Test summary table creation."""
        scorer = VariantScorer({})
        consensus = scorer.calculate_consensus_scores(sample_tissue_scores)
        ranked = scorer.rank_variants(consensus)

        original_variants = pd.DataFrame({
            'rsid': ['rs1', 'rs2'],
            'chromosome': ['1', '1'],
            'position': [1000000, 1500000]
        })

        result = create_summary_table(ranked, sample_tissue_scores, original_variants)

        assert 'final_rank' in result.columns
        # Result should have combined info
        assert len(result) > 0
