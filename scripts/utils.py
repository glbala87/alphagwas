#!/usr/bin/env python3
"""
Utility functions for AlphaGWAS pipeline.

Includes:
- Retry logic with exponential backoff
- Progress bar wrappers
- Logging configuration
- Input validation
- Performance utilities
"""

import functools
import logging
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, TypeVar, Union
import sys

import pandas as pd
import numpy as np

# Type variable for generic retry decorator
T = TypeVar('T')

# Try to import optional dependencies
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
    from rich.table import Table
    from rich.panel import Panel
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False


# Global console for rich output
console = Console() if RICH_AVAILABLE else None


def setup_logging(
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    use_rich: bool = True
) -> logging.Logger:
    """
    Configure logging with optional rich formatting and file output.

    Args:
        level: Logging level (default: INFO)
        log_file: Optional file path to write logs
        use_rich: Use rich formatting if available

    Returns:
        Configured logger instance
    """
    handlers = []

    if use_rich and RICH_AVAILABLE:
        handlers.append(RichHandler(
            console=console,
            show_time=True,
            show_path=False,
            rich_tracebacks=True
        ))
    else:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        )
        handlers.append(stream_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        )
        handlers.append(file_handler)

    logging.basicConfig(
        level=level,
        handlers=handlers,
        force=True
    )

    return logging.getLogger('alphagwas')


def retry_with_backoff(
    max_retries: int = 3,
    base_delay: float = 1.0,
    max_delay: float = 60.0,
    exponential_base: float = 2.0,
    exceptions: tuple = (Exception,),
    on_retry: Optional[Callable[[Exception, int], None]] = None
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator for retrying functions with exponential backoff.

    Args:
        max_retries: Maximum number of retry attempts
        base_delay: Initial delay between retries (seconds)
        max_delay: Maximum delay between retries (seconds)
        exponential_base: Base for exponential backoff
        exceptions: Tuple of exceptions to catch and retry
        on_retry: Optional callback function called on each retry

    Returns:
        Decorated function with retry logic

    Example:
        @retry_with_backoff(max_retries=3, exceptions=(ConnectionError,))
        def fetch_data(url):
            return requests.get(url)
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> T:
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e

                    if attempt < max_retries:
                        delay = min(base_delay * (exponential_base ** attempt), max_delay)
                        # Add jitter to prevent thundering herd
                        delay = delay * (0.5 + np.random.random())

                        if on_retry:
                            on_retry(e, attempt + 1)

                        logging.getLogger(__name__).warning(
                            f"Attempt {attempt + 1}/{max_retries + 1} failed: {e}. "
                            f"Retrying in {delay:.1f}s..."
                        )
                        time.sleep(delay)

            raise last_exception

        return wrapper
    return decorator


def progress_iterator(
    iterable,
    desc: str = "",
    total: Optional[int] = None,
    disable: bool = False,
    unit: str = "it"
):
    """
    Wrap an iterable with a progress bar.

    Uses tqdm if available, otherwise falls back to simple logging.

    Args:
        iterable: The iterable to wrap
        desc: Description prefix for the progress bar
        total: Total number of iterations (auto-detected if possible)
        disable: Disable progress bar
        unit: Unit name for iterations

    Yields:
        Items from the iterable
    """
    if disable:
        yield from iterable
        return

    if total is None:
        try:
            total = len(iterable)
        except TypeError:
            total = None

    if TQDM_AVAILABLE:
        yield from tqdm(iterable, desc=desc, total=total, unit=unit)
    else:
        # Fallback to simple logging
        logger = logging.getLogger(__name__)
        for i, item in enumerate(iterable):
            if total and (i + 1) % max(1, total // 10) == 0:
                logger.info(f"{desc}: {i + 1}/{total} ({100 * (i + 1) / total:.0f}%)")
            yield item


class ProgressTracker:
    """
    Context manager for tracking progress with rich or tqdm.

    Example:
        with ProgressTracker(total=100, desc="Processing") as tracker:
            for item in items:
                process(item)
                tracker.advance()
    """

    def __init__(
        self,
        total: int,
        desc: str = "Processing",
        disable: bool = False
    ):
        self.total = total
        self.desc = desc
        self.disable = disable
        self.progress = None
        self.task_id = None
        self.count = 0

    def __enter__(self):
        if self.disable:
            return self

        if RICH_AVAILABLE:
            self.progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                TimeElapsedColumn(),
                console=console
            )
            self.progress.__enter__()
            self.task_id = self.progress.add_task(self.desc, total=self.total)
        elif TQDM_AVAILABLE:
            self.progress = tqdm(total=self.total, desc=self.desc)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.progress:
            if RICH_AVAILABLE:
                self.progress.__exit__(exc_type, exc_val, exc_tb)
            elif TQDM_AVAILABLE:
                self.progress.close()

    def advance(self, amount: int = 1):
        """Advance the progress bar by the specified amount."""
        self.count += amount
        if self.progress:
            if RICH_AVAILABLE:
                self.progress.update(self.task_id, advance=amount)
            elif TQDM_AVAILABLE:
                self.progress.update(amount)


def validate_gwas_input(
    df: pd.DataFrame,
    required_columns: List[str],
    raise_on_error: bool = True
) -> Dict[str, Any]:
    """
    Validate GWAS input data.

    Args:
        df: Input DataFrame
        required_columns: List of required column names
        raise_on_error: Raise ValueError on validation failure

    Returns:
        Dictionary with validation results

    Raises:
        ValueError: If validation fails and raise_on_error is True
    """
    logger = logging.getLogger(__name__)
    issues = []
    warnings = []

    # Check for required columns
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        issues.append(f"Missing required columns: {missing_cols}")

    # Check for empty DataFrame
    if df.empty:
        issues.append("DataFrame is empty")

    # Check chromosome values
    if 'chromosome' in df.columns:
        valid_chroms = set(str(i) for i in range(1, 23)) | {'X', 'Y', 'MT', 'M'}
        invalid_chroms = set(df['chromosome'].astype(str).unique()) - valid_chroms
        if invalid_chroms:
            warnings.append(f"Non-standard chromosome values: {invalid_chroms}")

    # Check position values
    if 'position' in df.columns:
        if df['position'].isna().any():
            issues.append("Position column contains NA values")
        if (df['position'] <= 0).any():
            issues.append("Position column contains non-positive values")

    # Check p-value column
    if 'pvalue' in df.columns:
        if df['pvalue'].isna().any():
            warnings.append("P-value column contains NA values")
        if (df['pvalue'] < 0).any() or (df['pvalue'] > 1).any():
            issues.append("P-value column contains values outside [0, 1]")

    result = {
        'valid': len(issues) == 0,
        'issues': issues,
        'warnings': warnings,
        'n_rows': len(df),
        'n_columns': len(df.columns)
    }

    # Log validation results
    if issues:
        for issue in issues:
            logger.error(f"Validation error: {issue}")
    for warning in warnings:
        logger.warning(f"Validation warning: {warning}")

    if issues and raise_on_error:
        raise ValueError(f"Input validation failed: {'; '.join(issues)}")

    return result


def validate_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate pipeline configuration.

    Args:
        config: Configuration dictionary

    Returns:
        Dictionary with validation results

    Raises:
        ValueError: If required sections are missing
    """
    logger = logging.getLogger(__name__)

    required_sections = ['study', 'gwas', 'output']
    missing = [s for s in required_sections if s not in config]

    if missing:
        raise ValueError(f"Config missing required sections: {missing}")

    # Validate GWAS section
    if 'gwas' in config:
        gwas = config['gwas']
        if 'input_file' not in gwas:
            raise ValueError("Config 'gwas' section missing 'input_file'")
        if not Path(gwas['input_file']).exists():
            logger.warning(f"GWAS input file not found: {gwas['input_file']}")

    # Validate output section
    if 'output' in config:
        output = config['output']
        if 'dir' not in output:
            config['output']['dir'] = 'data/output'
        if 'prefix' not in output:
            config['output']['prefix'] = config.get('study', {}).get('name', 'study')

    return {'valid': True, 'config': config}


def print_summary_table(
    data: Dict[str, Any],
    title: str = "Summary"
):
    """
    Print a formatted summary table.

    Args:
        data: Dictionary of key-value pairs to display
        title: Table title
    """
    if RICH_AVAILABLE:
        table = Table(title=title)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green")

        for key, value in data.items():
            table.add_row(str(key), str(value))

        console.print(table)
    else:
        print(f"\n{title}")
        print("-" * 40)
        for key, value in data.items():
            print(f"  {key}: {value}")
        print()


def print_panel(message: str, title: str = "", style: str = "blue"):
    """
    Print a message in a styled panel.

    Args:
        message: Message to display
        title: Panel title
        style: Panel border style/color
    """
    if RICH_AVAILABLE:
        console.print(Panel(message, title=title, border_style=style))
    else:
        print(f"\n{'=' * 60}")
        if title:
            print(f" {title}")
        print(f"{'=' * 60}")
        print(message)
        print(f"{'=' * 60}\n")


def format_runtime(seconds: float) -> str:
    """Format runtime in human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def estimate_memory_usage(df: pd.DataFrame) -> str:
    """Estimate memory usage of a DataFrame."""
    bytes_used = df.memory_usage(deep=True).sum()

    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_used < 1024:
            return f"{bytes_used:.1f} {unit}"
        bytes_used /= 1024

    return f"{bytes_used:.1f} TB"
