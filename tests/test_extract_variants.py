"""
Tests for extract_variants module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.extract_variants import (
    load_gwas_sumstats,
    extract_significant_variants,
    extract_locus_variants,
    identify_lead_snps
)


class TestLoadGwasSumstats:
    """Tests for loading GWAS summary statistics."""

    def test_load_basic_file(self, temp_dir, sample_gwas_data, sample_config):
        """Test loading a basic TSV file."""
        # Save sample data
        gwas_path = Path(temp_dir) / 'test.tsv'
        sample_gwas_data.to_csv(gwas_path, sep='\t', index=False)
        sample_config['gwas']['input_file'] = str(gwas_path)

        df = load_gwas_sumstats(sample_config)

        assert len(df) == len(sample_gwas_data)
        assert 'chromosome' in df.columns
        assert 'position' in df.columns
        assert 'pvalue' in df.columns

    def test_column_renaming(self, temp_dir, sample_gwas_data, sample_config):
        """Test that columns are properly renamed."""
        gwas_path = Path(temp_dir) / 'test.tsv'
        sample_gwas_data.to_csv(gwas_path, sep='\t', index=False)
        sample_config['gwas']['input_file'] = str(gwas_path)

        df = load_gwas_sumstats(sample_config)

        # Original column names should be renamed
        assert 'CHR' not in df.columns
        assert 'BP' not in df.columns
        assert 'chromosome' in df.columns
        assert 'position' in df.columns


class TestExtractSignificantVariants:
    """Tests for extracting significant variants."""

    def test_extract_significant(self, sample_gwas_data, sample_config):
        """Test extraction of genome-wide significant variants."""
        # Rename columns to match expected format
        df = sample_gwas_data.rename(columns={
            'CHR': 'chromosome',
            'BP': 'position',
            'SNP': 'rsid',
            'P': 'pvalue'
        })

        result = extract_significant_variants(df, sample_config)

        # Should have 10 significant variants (from fixture)
        assert len(result) == 10
        assert all(result['pvalue'] < 5e-8)

    def test_no_significant_variants(self, sample_config):
        """Test when no variants pass threshold."""
        df = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, 200],
            'pvalue': [0.5, 0.8]
        })

        result = extract_significant_variants(df, sample_config)
        assert len(result) == 0

    def test_custom_threshold(self, sample_config):
        """Test with custom p-value threshold."""
        df = pd.DataFrame({
            'chromosome': ['1', '1', '2'],
            'position': [100, 200, 300],
            'pvalue': [1e-6, 1e-5, 1e-4]
        })

        sample_config['gwas']['pvalue_threshold'] = 1e-5
        result = extract_significant_variants(df, sample_config)

        assert len(result) == 1


class TestExtractLocusVariants:
    """Tests for extracting variants within loci."""

    def test_extract_locus(self, sample_config):
        """Test extraction of variants within defined loci."""
        df = pd.DataFrame({
            'chromosome': ['1', '1', '1', '2'],
            'position': [1000000, 1500000, 2000000, 1000000],
            'pvalue': [1e-10, 1e-9, 1e-8, 1e-10]
        })

        sample_config['loci'] = [{
            'name': 'test_locus',
            'chromosome': '1',
            'start': 900000,
            'end': 1600000
        }]

        result = extract_locus_variants(df, sample_config)

        assert len(result) == 2
        assert all(result['chromosome'].astype(str) == '1')
        assert all(result['position'] <= 1600000)

    def test_no_loci_defined(self, sample_config):
        """Test behavior when no loci are defined."""
        df = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, 200],
            'pvalue': [1e-10, 1e-9]
        })

        sample_config['loci'] = []
        result = extract_locus_variants(df, sample_config)

        # Should return original dataframe
        assert len(result) == 2


class TestIdentifyLeadSnps:
    """Tests for identifying lead SNPs."""

    def test_identify_leads_by_locus(self):
        """Test lead SNP identification when loci are defined."""
        df = pd.DataFrame({
            'chromosome': ['1', '1', '2', '2'],
            'position': [1000, 2000, 3000, 4000],
            'pvalue': [1e-10, 1e-8, 1e-9, 1e-7],
            'locus': ['locus1', 'locus1', 'locus2', 'locus2']
        })

        result = identify_lead_snps(df)

        assert 'is_lead' in result.columns
        leads = result[result['is_lead']]
        assert len(leads) == 2  # One per locus

        # Check that most significant per locus was selected
        locus1_lead = result[(result['locus'] == 'locus1') & (result['is_lead'])]
        assert locus1_lead['pvalue'].values[0] == 1e-10

    def test_identify_leads_by_clumping(self):
        """Test lead SNP identification by distance clumping."""
        df = pd.DataFrame({
            'chromosome': ['1', '1', '1'],
            'position': [1000000, 1100000, 2000000],  # First two close, third far
            'pvalue': [1e-10, 1e-9, 1e-8]
        })

        result = identify_lead_snps(df, window_kb=500)

        leads = result[result['is_lead']]
        # Should have 2 leads: first (most sig) and third (far from first)
        assert len(leads) == 2
