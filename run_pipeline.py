#!/usr/bin/env python3
"""
AlphaGWAS - GWAS Variant Prioritization Pipeline using AlphaGenome

Main entry point for running the complete pipeline or individual steps.

Features:
- Parallel AlphaGenome predictions
- Progress tracking with tqdm/rich
- Automatic retry with exponential backoff
- Visualization generation
- HTML summary reports

Usage:
    python run_pipeline.py                    # Run complete pipeline
    python run_pipeline.py --step 1           # Run specific step
    python run_pipeline.py --step 1,2,3       # Run multiple steps
    python run_pipeline.py --config my.yaml   # Use custom config
    python run_pipeline.py --visualize        # Generate visualizations only
"""

import argparse
import sys
from pathlib import Path
import logging
import yaml
import time

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / "scripts"))

from scripts import (
    extract_variants as step1,
    get_ld_proxies as step2,
    liftover as step3,
    alphagenome_predict as step4,
    score_variants as step5
)

# Try to import truthset validation module
try:
    from scripts import validate_truthset as step7
    VALIDATION_AVAILABLE = True
except ImportError:
    VALIDATION_AVAILABLE = False

# Try to import visualization module
try:
    from scripts import visualize as step6
    VISUALIZE_AVAILABLE = True
except ImportError:
    VISUALIZE_AVAILABLE = False

# Try to import utilities for better output
try:
    from scripts.utils import setup_logging, print_panel, print_summary_table, format_runtime
    UTILS_AVAILABLE = True
except ImportError:
    UTILS_AVAILABLE = False

# Try to import rich for better console output
try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None

if UTILS_AVAILABLE:
    logger = setup_logging()
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
    )
    logger = logging.getLogger(__name__)


STEPS = {
    1: ("Extract GWAS variants", step1.main, "Extracting GWAS-significant variants"),
    2: ("Get LD proxies", step2.main, "Finding LD proxy variants"),
    3: ("Liftover coordinates", step3.main, "Converting coordinates to hg38"),
    4: ("AlphaGenome predictions", step4.main, "Running AlphaGenome predictions"),
    5: ("Score and prioritize", step5.main, "Scoring and ranking variants"),
}

# Add visualization step if available
if VISUALIZE_AVAILABLE:
    STEPS[6] = ("Generate visualizations", step6.main, "Creating plots and reports")

# Add truthset validation step if available
if VALIDATION_AVAILABLE:
    STEPS[7] = ("Truthset validation", step7.main, "Validating results against known variants")


def validate_config(config_path: str) -> dict:
    """Validate configuration file."""
    if not Path(config_path).exists():
        logger.error(f"Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Check required sections
    required = ['study', 'gwas', 'reference', 'alphagenome', 'output']
    missing = [s for s in required if s not in config]
    if missing:
        logger.warning(f"Config missing sections: {missing}")

    return config


def check_dependencies():
    """Check for required Python packages."""
    required_packages = [
        ('pandas', 'pandas'),
        ('numpy', 'numpy'),
        ('yaml', 'pyyaml'),
    ]

    optional_packages = [
        ('pyliftover', 'pyliftover'),
        ('cyvcf2', 'cyvcf2'),
        ('pyarrow', 'pyarrow'),
    ]

    missing_required = []
    missing_optional = []

    for import_name, package_name in required_packages:
        try:
            __import__(import_name)
        except ImportError:
            missing_required.append(package_name)

    for import_name, package_name in optional_packages:
        try:
            __import__(import_name)
        except ImportError:
            missing_optional.append(package_name)

    if missing_required:
        logger.error(f"Missing required packages: {missing_required}")
        logger.error(f"Install with: pip install {' '.join(missing_required)}")
        sys.exit(1)

    if missing_optional:
        logger.warning(f"Missing optional packages: {missing_optional}")
        logger.warning("Some features may be limited. Install with:")
        logger.warning(f"  pip install {' '.join(missing_optional)}")


def run_step(step_num: int, config_path: str, verbose: bool = False) -> tuple:
    """
    Run a single pipeline step.

    Args:
        step_num: Step number to run
        config_path: Path to configuration file
        verbose: Enable verbose output

    Returns:
        Tuple of (success: bool, elapsed_time: float, result: any)
    """
    if step_num not in STEPS:
        logger.error(f"Invalid step number: {step_num}")
        return False, 0, None

    step_info = STEPS[step_num]
    step_name = step_info[0]
    step_func = step_info[1]
    step_description = step_info[2] if len(step_info) > 2 else step_name

    # Print step header
    if RICH_AVAILABLE:
        console.print(Panel(
            f"[bold]{step_description}[/bold]",
            title=f"Step {step_num}: {step_name}",
            border_style="blue"
        ))
    else:
        logger.info(f"\n{'='*60}")
        logger.info(f"Step {step_num}: {step_name}")
        logger.info(f"{'='*60}")

    start_time = time.time()

    try:
        result = step_func(config_path)
        elapsed = time.time() - start_time

        # Log success
        if RICH_AVAILABLE:
            console.print(f"[green]Step {step_num} completed in {elapsed:.1f}s[/green]")
        else:
            logger.info(f"Step {step_num} completed in {elapsed:.1f}s")

        return True, elapsed, result

    except Exception as e:
        elapsed = time.time() - start_time
        logger.error(f"Step {step_num} failed: {e}")

        if verbose:
            import traceback
            traceback.print_exc()

        return False, elapsed, None


def run_pipeline(
    config_path: str,
    steps: list = None,
    verbose: bool = False,
    skip_visualize: bool = False
):
    """
    Run the complete pipeline or specified steps.

    Args:
        config_path: Path to configuration file
        steps: List of step numbers to run (None = all)
        verbose: Enable verbose output
        skip_visualize: Skip visualization step
    """
    # Check dependencies
    check_dependencies()

    # Validate config
    config = validate_config(config_path)

    study_name = config.get('study', {}).get('name', 'unnamed')
    phenotype = config.get('study', {}).get('phenotype', 'unknown')

    # Print header
    if RICH_AVAILABLE:
        console.print(Panel(
            f"[bold cyan]Study:[/bold cyan] {study_name}\n"
            f"[bold cyan]Phenotype:[/bold cyan] {phenotype}\n"
            f"[bold cyan]Config:[/bold cyan] {config_path}",
            title="[bold white]AlphaGWAS Pipeline[/bold white]",
            border_style="cyan"
        ))
    else:
        logger.info(f"\n{'#'*60}")
        logger.info(f"# AlphaGWAS Pipeline")
        logger.info(f"# Study: {study_name}")
        logger.info(f"# Phenotype: {phenotype}")
        logger.info(f"{'#'*60}\n")

    # Determine which steps to run
    if steps is None:
        steps = [s for s in STEPS.keys() if s <= 5]  # Default: steps 1-5
        if VISUALIZE_AVAILABLE and not skip_visualize:
            steps.append(6)

    # Remove visualization step if requested
    if skip_visualize and 6 in steps:
        steps.remove(6)

    total_start = time.time()
    step_times = {}
    failed_steps = []
    completed_steps = []

    # Run steps
    for step_num in steps:
        success, elapsed, result = run_step(step_num, config_path, verbose)
        step_times[step_num] = elapsed

        if success:
            completed_steps.append(step_num)
        else:
            failed_steps.append(step_num)
            logger.error(f"Pipeline stopped at step {step_num}")
            break

    total_elapsed = time.time() - total_start

    # Print summary
    if RICH_AVAILABLE:
        # Create summary table
        table = Table(title="Pipeline Summary")
        table.add_column("Step", style="cyan")
        table.add_column("Name", style="white")
        table.add_column("Status", style="white")
        table.add_column("Time", style="white")

        for step_num in steps:
            name = STEPS[step_num][0]
            if step_num in completed_steps:
                status = "[green]Completed[/green]"
            elif step_num in failed_steps:
                status = "[red]Failed[/red]"
            else:
                status = "[yellow]Skipped[/yellow]"
            elapsed = step_times.get(step_num, 0)
            table.add_row(str(step_num), name, status, f"{elapsed:.1f}s")

        console.print(table)
        console.print(f"\n[bold]Total runtime:[/bold] {total_elapsed:.1f}s ({total_elapsed/60:.1f} minutes)")

    else:
        logger.info(f"\n{'='*60}")
        logger.info("PIPELINE SUMMARY")
        logger.info(f"{'='*60}")
        for step_num in steps:
            name = STEPS[step_num][0]
            status = "Completed" if step_num in completed_steps else "Failed" if step_num in failed_steps else "Skipped"
            elapsed = step_times.get(step_num, 0)
            logger.info(f"  Step {step_num}: {name} - {status} ({elapsed:.1f}s)")
        logger.info(f"\nTotal runtime: {total_elapsed:.1f}s")

    if failed_steps:
        logger.error(f"Failed steps: {failed_steps}")
        logger.info("Fix the errors and re-run with: --step " + ",".join(map(str, failed_steps)))
        return False
    else:
        if RICH_AVAILABLE:
            console.print("\n[bold green]All steps completed successfully![/bold green]")
        else:
            logger.info("All steps completed successfully!")

        output_dir = config.get('output', {}).get('dir', 'data/output')
        prefix = config.get('output', {}).get('prefix', 'study')
        logger.info(f"Results saved to: {output_dir}/{prefix}_*")

        # List output files
        output_path = Path(output_dir)
        if output_path.exists():
            output_files = list(output_path.glob(f"{prefix}_*"))
            if output_files and RICH_AVAILABLE:
                console.print("\n[bold]Output files:[/bold]")
                for f in sorted(output_files)[:10]:
                    console.print(f"  {f.name}")

        return True


def main():
    parser = argparse.ArgumentParser(
        description="AlphaGWAS - GWAS Variant Prioritization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_pipeline.py                      # Run complete pipeline
  python run_pipeline.py --step 4             # Run only AlphaGenome predictions
  python run_pipeline.py --step 4,5           # Run steps 4 and 5
  python run_pipeline.py --config my.yaml     # Use custom config
  python run_pipeline.py --visualize          # Generate visualizations only
  python run_pipeline.py --step 7             # Run truthset validation only

Steps:
  1. Extract GWAS-significant variants
  2. Get LD proxies from 1000 Genomes
  3. Liftover coordinates (hg19 -> hg38)
  4. Run AlphaGenome predictions
  5. Score and prioritize variants
  6. Generate visualizations (optional)
  7. Truthset validation (optional, requires validation config)
        """
    )

    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="Path to configuration file (default: config/config.yaml)"
    )

    parser.add_argument(
        "--step", "-s",
        default=None,
        help="Run specific step(s), comma-separated (e.g., '1,2,3')"
    )

    parser.add_argument(
        "--check",
        action="store_true",
        help="Check dependencies and config without running"
    )

    parser.add_argument(
        "--visualize", "-V",
        action="store_true",
        help="Generate visualizations only (requires completed pipeline)"
    )

    parser.add_argument(
        "--skip-visualize",
        action="store_true",
        help="Skip visualization step in full pipeline run"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output with full error tracebacks"
    )

    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress non-essential output"
    )

    args = parser.parse_args()

    # Set logging level
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    elif args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.check:
        check_dependencies()
        config = validate_config(args.config)

        if RICH_AVAILABLE:
            console.print("[green]All checks passed![/green]")

            # Show config summary
            table = Table(title="Configuration Summary")
            table.add_column("Setting", style="cyan")
            table.add_column("Value", style="white")

            table.add_row("Study", config.get('study', {}).get('name', 'N/A'))
            table.add_row("Phenotype", config.get('study', {}).get('phenotype', 'N/A'))
            table.add_row("Input file", config.get('gwas', {}).get('input_file', 'N/A'))
            table.add_row("P-value threshold", str(config.get('gwas', {}).get('pvalue_threshold', 'N/A')))
            table.add_row("Tissues", str(len(config.get('alphagenome', {}).get('tissues', []))))
            table.add_row("Output dir", config.get('output', {}).get('dir', 'N/A'))

            console.print(table)
        else:
            logger.info("All checks passed!")
        return

    # Visualize only mode
    if args.visualize:
        if not VISUALIZE_AVAILABLE:
            logger.error("Visualization module not available. Install matplotlib: pip install matplotlib seaborn")
            sys.exit(1)

        logger.info("Generating visualizations...")
        try:
            step6.main(args.config)
            logger.info("Visualizations generated successfully!")
            sys.exit(0)
        except Exception as e:
            logger.error(f"Visualization failed: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)

    # Parse steps
    steps = None
    if args.step:
        try:
            steps = [int(s.strip()) for s in args.step.split(",")]
        except ValueError:
            logger.error(f"Invalid step specification: {args.step}")
            sys.exit(1)

    # Run pipeline
    success = run_pipeline(
        args.config,
        steps,
        verbose=args.verbose,
        skip_visualize=args.skip_visualize
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
