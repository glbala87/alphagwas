"""
Pytest configuration and shared fixtures for AlphaGWAS tests.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import shutil
import yaml


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_config(temp_dir):
    """Create a sample configuration dictionary."""
    config = {
        'study': {
            'name': 'test_study',
            'phenotype': 'test_phenotype'
        },
        'gwas': {
            'input_file': f'{temp_dir}/test_gwas.tsv',
            'columns': {
                'chromosome': 'CHR',
                'position': 'BP',
                'rsid': 'SNP',
                'effect_allele': 'A1',
                'other_allele': 'A2',
                'pvalue': 'P'
            },
            'pvalue_threshold': 5e-8
        },
        'reference': {
            'input_build': 'hg19',
            'output_build': 'hg38'
        },
        'ld': {
            'vcf_path': f'{temp_dir}/vcf',
            'population': 'EUR',
            'r2_threshold': 0.8,
            'window_kb': 500
        },
        'alphagenome': {
            'tissues': ['Whole_Blood', 'Heart_Left_Ventricle'],
            'modalities': ['expression'],
            'max_workers': 2,
            'max_retries': 1,
            'batch_size': 10,
            'checkpoint_every': 5
        },
        'output': {
            'dir': f'{temp_dir}/output',
            'prefix': 'test'
        }
    }
    return config


@pytest.fixture
def config_file(temp_dir, sample_config):
    """Create a temporary config file."""
    config_path = Path(temp_dir) / 'config.yaml'
    with open(config_path, 'w') as f:
        yaml.dump(sample_config, f)
    return str(config_path)


@pytest.fixture
def sample_gwas_data():
    """Create sample GWAS summary statistics."""
    np.random.seed(42)
    n_variants = 100

    data = {
        'CHR': np.random.choice(['1', '2', '3'], n_variants),
        'BP': np.random.randint(1000000, 50000000, n_variants),
        'SNP': [f'rs{i}' for i in range(n_variants)],
        'A1': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'A2': np.random.choice(['A', 'C', 'G', 'T'], n_variants),
        'P': np.concatenate([
            np.random.uniform(1e-10, 1e-8, 10),  # Significant
            np.random.uniform(1e-5, 1, n_variants - 10)  # Not significant
        ])
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_gwas_file(temp_dir, sample_gwas_data):
    """Create a sample GWAS file."""
    gwas_path = Path(temp_dir) / 'test_gwas.tsv'
    sample_gwas_data.to_csv(gwas_path, sep='\t', index=False)
    return str(gwas_path)


@pytest.fixture
def sample_variants():
    """Create sample variant DataFrame."""
    return pd.DataFrame({
        'chromosome': ['1', '1', '2', '2', '3'],
        'position': [1000000, 1500000, 2000000, 2500000, 3000000],
        'rsid': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
        'ref': ['A', 'C', 'G', 'T', 'A'],
        'alt': ['G', 'T', 'A', 'C', 'G'],
        'pvalue': [1e-10, 2e-9, 5e-9, 1e-8, 3e-8],
        'effect_allele': ['G', 'T', 'A', 'C', 'G'],
        'other_allele': ['A', 'C', 'G', 'T', 'A']
    })


@pytest.fixture
def sample_predictions():
    """Create sample AlphaGenome predictions."""
    np.random.seed(42)

    variants = ['chr1:1000000', 'chr1:1500000', 'chr2:2000000']
    tissues = ['Whole_Blood', 'Heart_Left_Ventricle']
    modalities = ['expression', 'chromatin_accessibility']

    predictions = []
    for var in variants:
        for tissue in tissues:
            for modality in modalities:
                predictions.append({
                    'variant_id': var,
                    'tissue': tissue,
                    'modality': modality,
                    'prediction': np.random.normal(0, 0.5),
                    'effect_size': abs(np.random.normal(0, 0.3)),
                    'confidence': np.random.uniform(0.5, 1.0),
                    'nearby_genes': ['GENE1', 'GENE2']
                })

    return pd.DataFrame(predictions)


@pytest.fixture
def sample_tissue_scores():
    """Create sample tissue scores."""
    return pd.DataFrame({
        'variant_id': ['chr1:1000000', 'chr1:1000000', 'chr1:1500000', 'chr1:1500000'],
        'tissue': ['Whole_Blood', 'Heart_Left_Ventricle', 'Whole_Blood', 'Heart_Left_Ventricle'],
        'mean_effect': [0.3, 0.5, 0.2, 0.4],
        'max_effect': [0.5, 0.7, 0.3, 0.6],
        'effect_std': [0.1, 0.15, 0.08, 0.12],
        'mean_confidence': [0.8, 0.85, 0.75, 0.9],
        'max_abs_prediction': [0.4, 0.6, 0.25, 0.5],
        'tissue_score': [0.5, 0.8, 0.3, 0.7]
    })


@pytest.fixture
def mock_alphagenome_response():
    """Create mock AlphaGenome API response."""
    return {
        'reference': {'prediction': 0.5},
        'alternate': {'prediction': 0.8},
        'confidence': 0.9,
        'genes': ['GENE1', 'GENE2'],
        'elements': ['enhancer']
    }


class MockAlphaGenomeModel:
    """Mock AlphaGenome model for testing."""

    def __init__(self):
        self.call_count = 0

    def predict_variant(self, interval, variant, ontology_terms, requested_outputs):
        """Return mock predictions."""
        self.call_count += 1
        np.random.seed(self.call_count)

        return {
            'reference': {'prediction': np.random.uniform(0, 1)},
            'alternate': {'prediction': np.random.uniform(0, 1)},
            'confidence': np.random.uniform(0.5, 1.0),
            'genes': ['MOCK_GENE'],
            'elements': ['mock_element']
        }


@pytest.fixture
def mock_model():
    """Create a mock AlphaGenome model."""
    return MockAlphaGenomeModel()
