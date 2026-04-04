"""
Tests for truthset validation module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import shutil

from scripts.validate_truthset import (
    load_truthset,
    load_results,
    match_by_rsid,
    match_by_position,
    TruthsetValidator,
    ValidationMetrics,
    format_validation_report,
    save_validation_results,
    power_by_effect_size,
    power_by_maf,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def sample_truthset():
    """Create a sample truthset DataFrame."""
    return pd.DataFrame({
        'chromosome': ['1', '2', '3', '5', '7'],
        'position': [1000000, 2000000, 3000000, 5000000, 7000000],
        'rsid': ['rs1', 'rs2', 'rs3', 'rs5', 'rs7'],
        'gene': ['GENE_A', 'GENE_B', 'GENE_C', 'GENE_E', 'GENE_G'],
        'trait': ['T1', 'T1', 'T1', 'T1', 'T1'],
        'beta': [0.15, -0.22, 0.08, 0.35, 0.05],
        'maf': [0.30, 0.18, 0.22, 0.07, 0.45],
        'source': ['GWAS_Catalog'] * 5,
    })


@pytest.fixture
def sample_results():
    """Create sample pipeline results DataFrame."""
    return pd.DataFrame({
        'variant_id': [
            'chr1:1000000:A:G', 'chr2:2000050:C:T', 'chr3:3000000:G:A',
            'chr4:4000000:T:C', 'chr5:5500000:A:G', 'chr6:6000000:C:G',
        ],
        'rsid': ['rs1', 'rs99', 'rs3', 'rs4', 'rs55', 'rs66'],
        'chromosome': ['1', '2', '3', '4', '5', '6'],
        'position': [1000000, 2000050, 3000000, 4000000, 5500000, 6000000],
        'consensus_score': [0.93, 0.64, 0.45, 0.40, 0.35, 0.30],
        'final_rank': [1, 2, 3, 4, 5, 6],
    })


@pytest.fixture
def truthset_file(temp_dir, sample_truthset):
    """Write truthset to a TSV file."""
    path = Path(temp_dir) / 'truthset.tsv'
    sample_truthset.to_csv(path, sep='\t', index=False)
    return str(path)


@pytest.fixture
def results_file(temp_dir, sample_results):
    """Write results to a TSV file."""
    path = Path(temp_dir) / 'ranked_variants.tsv'
    sample_results.to_csv(path, sep='\t', index=False)
    return str(path)


# ---------------------------------------------------------------------------
# Tests: Loading
# ---------------------------------------------------------------------------

class TestLoadTruthset:
    def test_load_tsv(self, truthset_file):
        df = load_truthset(truthset_file)
        assert len(df) == 5
        assert 'chromosome' in df.columns
        assert 'position' in df.columns
        assert 'rsid' in df.columns

    def test_load_with_column_mapping(self, temp_dir):
        """Test loading with custom column names."""
        df = pd.DataFrame({
            'CHR': ['1', '2'], 'BP': [100, 200], 'SNP': ['rs1', 'rs2']
        })
        path = Path(temp_dir) / 'custom.tsv'
        df.to_csv(path, sep='\t', index=False)

        result = load_truthset(str(path))
        assert 'chromosome' in result.columns
        assert 'position' in result.columns
        assert 'rsid' in result.columns

    def test_chromosome_normalisation(self, temp_dir):
        """Test chr prefix stripping."""
        df = pd.DataFrame({
            'chromosome': ['chr1', 'chr2', 'chrX'],
            'position': [100, 200, 300],
            'rsid': ['rs1', 'rs2', 'rs3'],
        })
        path = Path(temp_dir) / 'chrprefix.tsv'
        df.to_csv(path, sep='\t', index=False)

        result = load_truthset(str(path))
        assert list(result['chromosome']) == ['1', '2', 'X']

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_truthset('/nonexistent/file.tsv')


class TestLoadResults:
    def test_load_results(self, results_file):
        df = load_results(results_file)
        assert len(df) == 6
        assert 'final_rank' in df.columns
        assert 'consensus_score' in df.columns

    def test_parse_variant_id(self, temp_dir):
        """Test chromosome/position parsing from variant_id."""
        df = pd.DataFrame({
            'variant_id': ['chr1:1000:A:G', 'chr2:2000:C:T'],
            'final_rank': [1, 2],
            'consensus_score': [0.9, 0.5],
        })
        path = Path(temp_dir) / 'results.tsv'
        df.to_csv(path, sep='\t', index=False)

        result = load_results(str(path))
        assert list(result['chromosome']) == ['1', '2']
        assert list(result['position']) == [1000, 2000]


# ---------------------------------------------------------------------------
# Tests: Matching
# ---------------------------------------------------------------------------

class TestMatchByRsid:
    def test_exact_matches(self, sample_truthset, sample_results):
        matches = match_by_rsid(sample_truthset, sample_results)
        matched = [m for m in matches if m.matched]
        # rs1 and rs3 should match
        assert len(matched) == 2
        matched_ids = {m.truth_id for m in matched}
        assert 'rs1' in matched_ids
        assert 'rs3' in matched_ids

    def test_match_details(self, sample_truthset, sample_results):
        matches = match_by_rsid(sample_truthset, sample_results)
        rs1_match = next(m for m in matches if m.truth_id == 'rs1')
        assert rs1_match.matched
        assert rs1_match.match_type == 'exact_rsid'
        assert rs1_match.pipeline_rank == 1
        assert rs1_match.pipeline_score == pytest.approx(0.93)

    def test_unmatched(self, sample_truthset, sample_results):
        matches = match_by_rsid(sample_truthset, sample_results)
        unmatched = [m for m in matches if not m.matched]
        # rs2, rs5, rs7 should not match
        assert len(unmatched) == 3

    def test_no_rsid_column(self):
        truth = pd.DataFrame({'chromosome': ['1'], 'position': [100]})
        results = pd.DataFrame({'variant_id': ['chr1:100:A:G']})
        matches = match_by_rsid(truth, results)
        assert len(matches) == 0


class TestMatchByPosition:
    def test_exact_match(self, sample_truthset, sample_results):
        matches = match_by_position(sample_truthset, sample_results, window_bp=0)
        matched = [m for m in matches if m.matched]
        # chr1:1000000 and chr3:3000000 are exact
        assert len(matched) == 2

    def test_window_match(self, sample_truthset, sample_results):
        matches = match_by_position(sample_truthset, sample_results, window_bp=100)
        matched = [m for m in matches if m.matched]
        # chr1:1000000 exact, chr2:2000000 within 50bp of 2000050, chr3:3000000 exact
        assert len(matched) == 3
        rs2_match = next(m for m in matches if m.truth_id == 'rs2')
        assert rs2_match.matched
        assert rs2_match.distance_bp == 50
        assert rs2_match.match_type == 'window'

    def test_large_window(self, sample_truthset, sample_results):
        matches = match_by_position(sample_truthset, sample_results, window_bp=1_000_000)
        matched = [m for m in matches if m.matched]
        # More variants should match with a large window
        assert len(matched) >= 3

    def test_different_chromosome(self, sample_truthset, sample_results):
        """Variants on chr7 should not match anything."""
        matches = match_by_position(sample_truthset, sample_results, window_bp=100)
        rs7_match = next(m for m in matches if m.truth_id == 'rs7')
        assert not rs7_match.matched


# ---------------------------------------------------------------------------
# Tests: Validator
# ---------------------------------------------------------------------------

class TestTruthsetValidator:
    def test_rsid_strategy(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        assert metrics.true_positives == 2
        assert metrics.false_negatives == 3
        assert metrics.recall == pytest.approx(2 / 5)

    def test_window_strategy(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='window', window_bp=100)
        metrics = validator.validate(sample_truthset, sample_results)
        assert metrics.true_positives == 3
        assert metrics.recall == pytest.approx(3 / 5)

    def test_auto_strategy(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='auto', window_bp=100)
        metrics = validator.validate(sample_truthset, sample_results)
        # Auto should pick whichever gives more matches (window: 3 > rsid: 2)
        assert metrics.true_positives == 3

    def test_precision(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        # TP=2, total results=6, matched results=2, FP=6-2=4
        assert metrics.precision == pytest.approx(2 / 6)

    def test_f1_score(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        p = metrics.precision
        r = metrics.recall
        expected_f1 = 2 * p * r / (p + r)
        assert metrics.f1_score == pytest.approx(expected_f1)

    def test_ranking_quality(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        # rs1 has rank 1, rs3 has rank 3
        assert metrics.mean_rank_of_true == pytest.approx(2.0)
        assert metrics.median_rank_of_true == pytest.approx(2.0)

    def test_perfect_recall(self):
        """All truth variants found."""
        truth = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, 200],
            'rsid': ['rs1', 'rs2'],
        })
        results = pd.DataFrame({
            'variant_id': ['chr1:100:A:G', 'chr2:200:C:T'],
            'rsid': ['rs1', 'rs2'],
            'chromosome': ['1', '2'],
            'position': [100, 200],
            'final_rank': [1, 2],
            'consensus_score': [0.9, 0.8],
        })
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(truth, results)
        assert metrics.recall == 1.0
        assert metrics.true_positives == 2
        assert metrics.false_negatives == 0

    def test_zero_recall(self):
        """No truth variants found."""
        truth = pd.DataFrame({
            'chromosome': ['1'],
            'position': [100],
            'rsid': ['rs1'],
        })
        results = pd.DataFrame({
            'variant_id': ['chr2:200:C:T'],
            'rsid': ['rs99'],
            'chromosome': ['2'],
            'position': [200],
            'final_rank': [1],
            'consensus_score': [0.9],
        })
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(truth, results)
        assert metrics.recall == 0.0
        assert metrics.f1_score == 0.0

    def test_n_total_variants(self, sample_truthset, sample_results):
        """Test specificity with explicit total variant count."""
        validator = TruthsetValidator(
            match_strategy='rsid', n_total_variants=1000
        )
        metrics = validator.validate(sample_truthset, sample_results)
        assert metrics.true_negatives > 0
        assert metrics.specificity > 0


# ---------------------------------------------------------------------------
# Tests: Reporting
# ---------------------------------------------------------------------------

class TestReporting:
    def test_format_report(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        report = format_validation_report(metrics)
        assert 'TRUTHSET VALIDATION REPORT' in report
        assert 'Precision' in report
        assert 'Recall' in report
        assert 'Matched Variants' in report
        assert 'Missed Variants' in report

    def test_save_results(self, temp_dir, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        saved = save_validation_results(metrics, temp_dir, prefix='test')

        assert Path(saved['metrics']).exists()
        assert Path(saved['details']).exists()
        assert Path(saved['report']).exists()

        # Verify metrics file content
        metrics_df = pd.read_csv(saved['metrics'], sep='\t')
        assert 'metric' in metrics_df.columns
        assert 'value' in metrics_df.columns
        assert 'precision' in metrics_df['metric'].values

    def test_details_file(self, temp_dir, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, sample_results)
        saved = save_validation_results(metrics, temp_dir, prefix='test')

        details = pd.read_csv(saved['details'], sep='\t')
        assert len(details) == 5  # all truth variants
        assert 'truth_id' in details.columns
        assert 'matched' in details.columns


# ---------------------------------------------------------------------------
# Tests: Power analysis
# ---------------------------------------------------------------------------

class TestPowerAnalysis:
    def test_power_by_effect_size(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        result = power_by_effect_size(
            sample_truthset, sample_results, validator,
            effect_col='beta', bins=2
        )
        assert len(result) > 0
        assert 'effect_bin' in result.columns
        assert 'recall' in result.columns
        assert all(0 <= r <= 1 for r in result['recall'])

    def test_power_missing_column(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        result = power_by_effect_size(
            sample_truthset, sample_results, validator,
            effect_col='nonexistent'
        )
        assert result.empty

    def test_power_by_maf(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        result = power_by_maf(
            sample_truthset, sample_results, validator,
            maf_col='maf', bins=2
        )
        assert len(result) > 0
        assert 'maf_bin' in result.columns
        assert 'recall' in result.columns

    def test_maf_missing_column(self, sample_truthset, sample_results):
        validator = TruthsetValidator(match_strategy='rsid')
        result = power_by_maf(
            sample_truthset, sample_results, validator,
            maf_col='nonexistent'
        )
        assert result.empty


# ---------------------------------------------------------------------------
# Tests: Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_empty_truthset(self, sample_results):
        truth = pd.DataFrame(columns=['chromosome', 'position', 'rsid'])
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(truth, sample_results)
        assert metrics.n_truth == 0
        assert metrics.recall == 0.0

    def test_empty_results(self, sample_truthset):
        results = pd.DataFrame(columns=[
            'variant_id', 'rsid', 'chromosome', 'position',
            'final_rank', 'consensus_score'
        ])
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(sample_truthset, results)
        assert metrics.n_results == 0
        assert metrics.true_positives == 0
        assert metrics.false_negatives == 5

    def test_duplicate_rsids_in_results(self):
        """Handle results with duplicate rsids gracefully."""
        truth = pd.DataFrame({
            'chromosome': ['1'], 'position': [100], 'rsid': ['rs1']
        })
        results = pd.DataFrame({
            'variant_id': ['chr1:100:A:G', 'chr1:100:A:T'],
            'rsid': ['rs1', 'rs1'],
            'chromosome': ['1', '1'],
            'position': [100, 100],
            'final_rank': [1, 2],
            'consensus_score': [0.9, 0.8],
        })
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(truth, results)
        assert metrics.true_positives == 1

    def test_top_fraction_metrics(self):
        """Test top-N fraction calculations."""
        truth = pd.DataFrame({
            'rsid': [f'rs{i}' for i in range(20)],
            'chromosome': ['1'] * 20,
            'position': list(range(100, 2100, 100)),
        })
        results = pd.DataFrame({
            'variant_id': [f'chr1:{p}:A:G' for p in range(100, 5100, 100)],
            'rsid': [f'rs{i}' for i in range(50)],
            'chromosome': ['1'] * 50,
            'position': list(range(100, 5100, 100)),
            'final_rank': list(range(1, 51)),
            'consensus_score': [1.0 - i * 0.02 for i in range(50)],
        })
        validator = TruthsetValidator(match_strategy='rsid')
        metrics = validator.validate(truth, results)
        assert metrics.true_positives == 20
        assert metrics.recall == 1.0
        assert metrics.fraction_in_top10 == pytest.approx(10 / 20)
        assert metrics.fraction_in_top20 == pytest.approx(20 / 20)
        assert metrics.fraction_in_top50 == pytest.approx(20 / 20)
