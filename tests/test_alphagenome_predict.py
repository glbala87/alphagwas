"""
Tests for alphagenome_predict module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.alphagenome_predict import (
    AlphaGenomePredictor,
    TISSUE_ONTOLOGY,
    VariantPrediction
)


class TestTissueOntology:
    """Tests for tissue ontology mapping."""

    def test_all_tissues_have_uberon(self):
        """Test that all tissues have UBERON IDs."""
        for tissue, uberon in TISSUE_ONTOLOGY.items():
            assert uberon.startswith('UBERON:')
            assert len(uberon) > 7  # UBERON: + ID

    def test_key_tissues_present(self):
        """Test that key tissues are in the mapping."""
        required = ['Whole_Blood', 'Heart_Left_Ventricle', 'Liver', 'Brain_Cortex']
        for tissue in required:
            assert tissue in TISSUE_ONTOLOGY


class TestAlphaGenomePredictor:
    """Tests for AlphaGenomePredictor class."""

    def test_init_without_api(self, sample_config):
        """Test initialization without AlphaGenome API."""
        predictor = AlphaGenomePredictor(sample_config)

        assert predictor.tissues == ['Whole_Blood', 'Heart_Left_Ventricle']
        assert predictor.modalities == ['expression']
        assert predictor.max_workers == 2
        assert predictor.model is None  # No API key

    def test_mock_predict(self, sample_config):
        """Test mock prediction (without API)."""
        predictor = AlphaGenomePredictor(sample_config)

        result = predictor.predict_variant(
            chrom='1',
            pos=1000000,
            ref='A',
            alt='G',
            tissue='Whole_Blood',
            modality='expression'
        )

        assert result is not None
        assert 'variant_id' in result
        assert 'tissue' in result
        assert 'modality' in result
        assert 'prediction' in result
        assert 'effect_size' in result
        assert 'confidence' in result
        assert result['tissue'] == 'Whole_Blood'
        assert result['modality'] == 'expression'

    def test_mock_predict_deterministic(self, sample_config):
        """Test that mock predictions are deterministic for same input."""
        predictor = AlphaGenomePredictor(sample_config)

        result1 = predictor._mock_predict('1', 1000000, 'A', 'G', 'Whole_Blood', 'expression')
        result2 = predictor._mock_predict('1', 1000000, 'A', 'G', 'Whole_Blood', 'expression')

        assert result1['prediction'] == result2['prediction']
        assert result1['effect_size'] == result2['effect_size']

    def test_mock_predict_varies_by_input(self, sample_config):
        """Test that mock predictions vary for different inputs."""
        predictor = AlphaGenomePredictor(sample_config)

        result1 = predictor._mock_predict('1', 1000000, 'A', 'G', 'Whole_Blood', 'expression')
        result2 = predictor._mock_predict('1', 2000000, 'A', 'G', 'Whole_Blood', 'expression')

        # Predictions should be different for different positions
        assert result1['prediction'] != result2['prediction']

    def test_get_ontology_term(self, sample_config):
        """Test ontology term retrieval."""
        predictor = AlphaGenomePredictor(sample_config)

        term = predictor._get_ontology_term('Whole_Blood')
        assert term == 'UBERON:0000178'

        # Unknown tissue should default to blood
        unknown_term = predictor._get_ontology_term('Unknown_Tissue')
        assert unknown_term == 'UBERON:0000178'

    def test_predict_batch_sequential(self, sample_config, sample_variants, temp_dir):
        """Test sequential batch prediction."""
        sample_config['alphagenome']['max_workers'] = 1
        predictor = AlphaGenomePredictor(sample_config)

        results = predictor.predict_batch(
            sample_variants,
            checkpoint_path=None,
            use_parallel=False
        )

        expected_count = len(sample_variants) * len(predictor.tissues) * len(predictor.modalities)
        assert len(results) == expected_count

    def test_predict_batch_parallel(self, sample_config, sample_variants, temp_dir):
        """Test parallel batch prediction."""
        predictor = AlphaGenomePredictor(sample_config)

        results = predictor.predict_batch(
            sample_variants,
            checkpoint_path=None,
            use_parallel=True
        )

        expected_count = len(sample_variants) * len(predictor.tissues) * len(predictor.modalities)
        assert len(results) == expected_count

    def test_checkpoint_save_and_resume(self, sample_config, sample_variants, temp_dir):
        """Test checkpoint saving and resuming."""
        checkpoint_path = Path(temp_dir) / 'checkpoint.json'
        predictor = AlphaGenomePredictor(sample_config)

        # First run - process half
        variants_half = sample_variants.head(2)
        results1 = predictor.predict_batch(
            variants_half,
            checkpoint_path=str(checkpoint_path),
            use_parallel=False
        )

        assert checkpoint_path.exists()

        # Second run - should resume
        results2 = predictor.predict_batch(
            sample_variants,
            checkpoint_path=str(checkpoint_path),
            use_parallel=False
        )

        # Should have more results than first run
        assert len(results2) >= len(results1)


class TestVariantPrediction:
    """Tests for VariantPrediction dataclass."""

    def test_create_prediction(self):
        """Test creating a VariantPrediction."""
        pred = VariantPrediction(
            variant_id='chr1:1000000',
            chromosome='1',
            position=1000000,
            ref='A',
            alt='G',
            tissue='Whole_Blood',
            modality='expression',
            prediction=0.5,
            effect_size=0.3,
            confidence=0.9
        )

        assert pred.variant_id == 'chr1:1000000'
        assert pred.effect_size == 0.3


class TestModalityMapping:
    """Tests for modality to output type mapping."""

    def test_modality_mapping(self, sample_config):
        """Test that modalities map correctly."""
        predictor = AlphaGenomePredictor(sample_config)

        # Check mapping exists
        assert 'expression' in predictor.MODALITY_MAPPING
        assert predictor.MODALITY_MAPPING['expression'] == 'RNA_SEQ'
        assert predictor.MODALITY_MAPPING['chromatin_accessibility'] == 'ATAC_SEQ'
