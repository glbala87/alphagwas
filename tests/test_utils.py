"""
Tests for utils module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import time
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils import (
    retry_with_backoff,
    progress_iterator,
    validate_gwas_input,
    validate_config,
    format_runtime,
    estimate_memory_usage
)


class TestRetryWithBackoff:
    """Tests for retry decorator."""

    def test_success_first_try(self):
        """Test that function succeeds on first try."""
        call_count = 0

        @retry_with_backoff(max_retries=3)
        def succeed():
            nonlocal call_count
            call_count += 1
            return "success"

        result = succeed()

        assert result == "success"
        assert call_count == 1

    def test_success_after_retry(self):
        """Test that function succeeds after retry."""
        call_count = 0

        @retry_with_backoff(max_retries=3, base_delay=0.01)
        def fail_then_succeed():
            nonlocal call_count
            call_count += 1
            if call_count < 2:
                raise ValueError("Temporary error")
            return "success"

        result = fail_then_succeed()

        assert result == "success"
        assert call_count == 2

    def test_max_retries_exceeded(self):
        """Test that exception is raised after max retries."""
        call_count = 0

        @retry_with_backoff(max_retries=2, base_delay=0.01)
        def always_fail():
            nonlocal call_count
            call_count += 1
            raise ValueError("Persistent error")

        with pytest.raises(ValueError):
            always_fail()

        assert call_count == 3  # Initial + 2 retries

    def test_specific_exception(self):
        """Test that only specified exceptions are caught."""
        @retry_with_backoff(max_retries=3, exceptions=(ValueError,), base_delay=0.01)
        def raise_type_error():
            raise TypeError("Wrong type")

        with pytest.raises(TypeError):
            raise_type_error()

    def test_on_retry_callback(self):
        """Test that on_retry callback is called."""
        retry_calls = []

        def on_retry(exc, attempt):
            retry_calls.append((str(exc), attempt))

        @retry_with_backoff(max_retries=2, base_delay=0.01, on_retry=on_retry)
        def fail_twice():
            if len(retry_calls) < 2:
                raise ValueError("Error")
            return "success"

        result = fail_twice()

        assert result == "success"
        assert len(retry_calls) == 2


class TestProgressIterator:
    """Tests for progress iterator."""

    def test_iterate_list(self):
        """Test iterating over a list."""
        items = [1, 2, 3, 4, 5]
        result = list(progress_iterator(items, desc="Test", disable=True))

        assert result == items

    def test_iterate_with_total(self):
        """Test iterating with explicit total."""
        items = range(5)
        result = list(progress_iterator(items, total=5, disable=True))

        assert result == [0, 1, 2, 3, 4]


class TestValidateGwasInput:
    """Tests for GWAS input validation."""

    def test_valid_input(self):
        """Test validation of valid input."""
        df = pd.DataFrame({
            'chromosome': ['1', '2', '3'],
            'position': [100, 200, 300],
            'pvalue': [0.01, 0.02, 0.03]
        })

        result = validate_gwas_input(df, ['chromosome', 'position'])

        assert result['valid'] is True
        assert len(result['issues']) == 0

    def test_missing_columns(self):
        """Test validation with missing columns."""
        df = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, 200]
        })

        with pytest.raises(ValueError) as exc_info:
            validate_gwas_input(df, ['chromosome', 'position', 'pvalue'])

        assert 'Missing required columns' in str(exc_info.value)

    def test_invalid_positions(self):
        """Test validation with invalid positions."""
        df = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, -5],  # Negative position
            'pvalue': [0.01, 0.02]
        })

        with pytest.raises(ValueError) as exc_info:
            validate_gwas_input(df, ['chromosome', 'position'])

        assert 'non-positive' in str(exc_info.value)

    def test_invalid_pvalues(self):
        """Test validation with invalid p-values."""
        df = pd.DataFrame({
            'chromosome': ['1', '2'],
            'position': [100, 200],
            'pvalue': [0.5, 1.5]  # >1 p-value
        })

        with pytest.raises(ValueError) as exc_info:
            validate_gwas_input(df, ['chromosome', 'position', 'pvalue'])

        assert 'outside [0, 1]' in str(exc_info.value)

    def test_empty_dataframe(self):
        """Test validation of empty DataFrame."""
        df = pd.DataFrame()

        with pytest.raises(ValueError) as exc_info:
            validate_gwas_input(df, ['chromosome'])

        assert 'empty' in str(exc_info.value).lower()

    def test_warnings_only(self):
        """Test validation that produces warnings but not errors."""
        df = pd.DataFrame({
            'chromosome': ['1', '25'],  # 25 is non-standard
            'position': [100, 200]
        })

        result = validate_gwas_input(df, ['chromosome', 'position'], raise_on_error=False)

        assert result['valid'] is True
        assert len(result['warnings']) > 0


class TestValidateConfig:
    """Tests for config validation."""

    def test_valid_config(self, temp_dir):
        """Test validation of valid config."""
        config = {
            'study': {'name': 'test'},
            'gwas': {'input_file': __file__},  # Use this file as dummy
            'output': {'dir': temp_dir}
        }

        result = validate_config(config)

        assert result['valid'] is True

    def test_missing_sections(self):
        """Test validation with missing sections."""
        config = {'study': {'name': 'test'}}

        with pytest.raises(ValueError) as exc_info:
            validate_config(config)

        assert 'missing required sections' in str(exc_info.value).lower()

    def test_default_output_values(self, temp_dir):
        """Test that defaults are applied."""
        config = {
            'study': {'name': 'my_study'},
            'gwas': {'input_file': __file__},
            'output': {}
        }

        result = validate_config(config)

        assert result['config']['output']['dir'] == 'data/output'
        assert result['config']['output']['prefix'] == 'my_study'


class TestFormatRuntime:
    """Tests for runtime formatting."""

    def test_format_seconds(self):
        """Test formatting seconds."""
        assert format_runtime(30.5) == "30.5s"
        assert format_runtime(0.5) == "0.5s"

    def test_format_minutes(self):
        """Test formatting minutes."""
        assert format_runtime(120) == "2.0m"
        assert format_runtime(90) == "1.5m"

    def test_format_hours(self):
        """Test formatting hours."""
        assert format_runtime(3600) == "1.0h"
        assert format_runtime(5400) == "1.5h"


class TestEstimateMemoryUsage:
    """Tests for memory estimation."""

    def test_estimate_small_df(self):
        """Test memory estimation for small DataFrame."""
        df = pd.DataFrame({
            'a': [1, 2, 3],
            'b': ['x', 'y', 'z']
        })

        result = estimate_memory_usage(df)

        assert 'B' in result or 'KB' in result

    def test_estimate_larger_df(self):
        """Test memory estimation for larger DataFrame."""
        df = pd.DataFrame({
            'a': range(10000),
            'b': ['x' * 100] * 10000
        })

        result = estimate_memory_usage(df)

        assert any(unit in result for unit in ['KB', 'MB'])
