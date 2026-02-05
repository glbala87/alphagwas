#!/usr/bin/env python3
"""
AlphaGWAS - GWAS Variant Prioritization Pipeline using AlphaGenome

Main entry point for running the complete pipeline or individual steps.

Usage:
    python run_pipeline.py                    # Run complete pipeline
    python run_pipeline.py --step 1           # Run specific step
    python run_pipeline.py --step 1,2,3       # Run multiple steps
    python run_pipeline.py --config my.yaml   # Use custom config
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

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
)
logger = logging.getLogger(__name__)


STEPS = {
    1: ("Extract GWAS variants", step1.main),
    2: ("Get LD proxies", step2.main),
    3: ("Liftover coordinates", step3.main),
    4: ("AlphaGenome predictions", step4.main),
    5: ("Score and prioritize", step5.main),
}


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


def run_step(step_num: int, config_path: str) -> bool:
    """Run a single pipeline step."""
    if step_num not in STEPS:
        logger.error(f"Invalid step number: {step_num}")
        return False

    step_name, step_func = STEPS[step_num]

    logger.info(f"\n{'='*60}")
    logger.info(f"Step {step_num}: {step_name}")
    logger.info(f"{'='*60}")

    start_time = time.time()

    try:
        result = step_func(config_path)
        elapsed = time.time() - start_time
        logger.info(f"Step {step_num} completed in {elapsed:.1f}s")
        return True
    except Exception as e:
        logger.error(f"Step {step_num} failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def run_pipeline(config_path: str, steps: list = None):
    """Run the complete pipeline or specified steps."""

    # Check dependencies
    check_dependencies()

    # Validate config
    config = validate_config(config_path)

    study_name = config.get('study', {}).get('name', 'unnamed')
    phenotype = config.get('study', {}).get('phenotype', 'unknown')

    logger.info(f"\n{'#'*60}")
    logger.info(f"# AlphaGWAS Pipeline")
    logger.info(f"# Study: {study_name}")
    logger.info(f"# Phenotype: {phenotype}")
    logger.info(f"{'#'*60}\n")

    # Determine which steps to run
    if steps is None:
        steps = list(STEPS.keys())

    total_start = time.time()
    failed_steps = []

    for step_num in steps:
        success = run_step(step_num, config_path)
        if not success:
            failed_steps.append(step_num)
            logger.error(f"Pipeline stopped at step {step_num}")
            break

    total_elapsed = time.time() - total_start

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("PIPELINE SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Total runtime: {total_elapsed:.1f}s")

    if failed_steps:
        logger.error(f"Failed steps: {failed_steps}")
        logger.info("Fix the errors and re-run with: --step " + ",".join(map(str, failed_steps)))
        return False
    else:
        logger.info("All steps completed successfully!")
        output_dir = config.get('output', {}).get('dir', 'data/output')
        prefix = config.get('output', {}).get('prefix', 'study')
        logger.info(f"\nResults saved to: {output_dir}/{prefix}_*")
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

Steps:
  1. Extract GWAS-significant variants
  2. Get LD proxies from 1000 Genomes
  3. Liftover coordinates (hg19 -> hg38)
  4. Run AlphaGenome predictions
  5. Score and prioritize variants
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

    args = parser.parse_args()

    if args.check:
        check_dependencies()
        validate_config(args.config)
        logger.info("All checks passed!")
        return

    # Parse steps
    steps = None
    if args.step:
        try:
            steps = [int(s.strip()) for s in args.step.split(",")]
        except ValueError:
            logger.error(f"Invalid step specification: {args.step}")
            sys.exit(1)

    # Run pipeline
    success = run_pipeline(args.config, steps)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
