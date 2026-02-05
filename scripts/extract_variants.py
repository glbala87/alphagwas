#!/usr/bin/env python3
"""
Step 1: Extract GWAS-significant variants from summary statistics.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import yaml

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    """Load pipeline configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_gwas_sumstats(config: dict) -> pd.DataFrame:
    """Load GWAS summary statistics file."""
    gwas_config = config['gwas']
    input_file = gwas_config['input_file']

    logger.info(f"Loading GWAS summary statistics from {input_file}")

    # Detect file format
    if input_file.endswith('.gz'):
        df = pd.read_csv(input_file, sep='\t', compression='gzip')
    else:
        df = pd.read_csv(input_file, sep='\t')

    # Rename columns to standard names
    col_map = gwas_config['columns']
    rename_map = {v: k for k, v in col_map.items()}
    df = df.rename(columns=rename_map)

    logger.info(f"Loaded {len(df):,} variants")
    return df


def extract_significant_variants(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Extract genome-wide significant variants."""
    pval_threshold = config['gwas']['pvalue_threshold']

    significant = df[df['pvalue'] < pval_threshold].copy()
    logger.info(f"Found {len(significant):,} genome-wide significant variants (p < {pval_threshold})")

    return significant


def extract_locus_variants(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Extract variants within defined loci."""
    loci = config.get('loci', [])
    if not loci:
        logger.warning("No loci defined in config, returning all significant variants")
        return df

    locus_variants = []

    for locus in loci:
        chrom = str(locus['chromosome'])
        start = locus['start']
        end = locus['end']
        name = locus['name']

        # Filter variants in this locus
        mask = (
            (df['chromosome'].astype(str) == chrom) &
            (df['position'] >= start) &
            (df['position'] <= end)
        )
        locus_df = df[mask].copy()
        locus_df['locus'] = name

        logger.info(f"Locus {name} (chr{chrom}:{start}-{end}): {len(locus_df):,} variants")
        locus_variants.append(locus_df)

    if locus_variants:
        return pd.concat(locus_variants, ignore_index=True)
    return pd.DataFrame()


def identify_lead_snps(df: pd.DataFrame, window_kb: int = 500) -> pd.DataFrame:
    """Identify lead SNPs by taking the most significant variant per locus."""
    if 'locus' not in df.columns:
        # Define loci by clumping
        df = df.sort_values('pvalue')
        df['is_lead'] = False

        lead_positions = []
        for idx, row in df.iterrows():
            chrom = row['chromosome']
            pos = row['position']

            # Check if within window of existing lead
            is_independent = True
            for lead_chrom, lead_pos in lead_positions:
                if chrom == lead_chrom and abs(pos - lead_pos) < window_kb * 1000:
                    is_independent = False
                    break

            if is_independent:
                df.loc[idx, 'is_lead'] = True
                lead_positions.append((chrom, pos))
    else:
        # Take most significant per locus
        df['is_lead'] = False
        for locus in df['locus'].unique():
            locus_mask = df['locus'] == locus
            lead_idx = df.loc[locus_mask, 'pvalue'].idxmin()
            df.loc[lead_idx, 'is_lead'] = True

    lead_count = df['is_lead'].sum()
    logger.info(f"Identified {lead_count} lead SNPs")

    return df


def main(config_path: str = "config/config.yaml"):
    """Main extraction workflow."""
    config = load_config(config_path)

    # Create output directory
    output_dir = Path(config['output']['dir'])
    intermediate_dir = Path("data/intermediate")
    output_dir.mkdir(parents=True, exist_ok=True)
    intermediate_dir.mkdir(parents=True, exist_ok=True)

    # Load and process GWAS data
    df = load_gwas_sumstats(config)

    # Extract significant variants
    significant = extract_significant_variants(df, config)

    # Filter to loci of interest
    locus_variants = extract_locus_variants(significant, config)

    # If no loci defined, use all significant variants
    if locus_variants.empty:
        locus_variants = significant

    # Identify lead SNPs
    result = identify_lead_snps(locus_variants)

    # Save results
    prefix = config['output']['prefix']
    output_file = intermediate_dir / f"{prefix}_significant_variants.tsv"
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(result):,} variants to {output_file}")

    # Save lead SNPs separately
    lead_snps = result[result['is_lead']]
    lead_file = intermediate_dir / f"{prefix}_lead_snps.tsv"
    lead_snps.to_csv(lead_file, sep='\t', index=False)
    logger.info(f"Saved {len(lead_snps)} lead SNPs to {lead_file}")

    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Extract GWAS-significant variants")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
