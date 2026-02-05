#!/usr/bin/env python3
"""
Step 2: Identify LD proxies from 1000 Genomes reference panel.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import logging
import yaml
from typing import List, Tuple, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def calculate_ld_plink(
    snp_list: List[str],
    chrom: str,
    vcf_path: str,
    population: str,
    r2_threshold: float,
    window_kb: int,
    output_prefix: str
) -> pd.DataFrame:
    """
    Calculate LD using PLINK2 (if available) or fall back to manual calculation.
    """
    # Try using plink2 for LD calculation
    vcf_file = Path(vcf_path) / f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

    if not vcf_file.exists():
        logger.warning(f"VCF file not found: {vcf_file}")
        logger.info("Falling back to LDlink API or returning lead SNPs only")
        return pd.DataFrame()

    # Write SNP list to file
    snp_file = f"{output_prefix}_snps.txt"
    with open(snp_file, 'w') as f:
        for snp in snp_list:
            f.write(f"{snp}\n")

    # Run PLINK2 LD calculation
    cmd = [
        "plink2",
        "--vcf", str(vcf_file),
        "--extract", snp_file,
        "--r2",
        "--ld-window-kb", str(window_kb),
        "--ld-window-r2", str(r2_threshold),
        "--out", output_prefix
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        ld_file = f"{output_prefix}.ld"
        if Path(ld_file).exists():
            ld_df = pd.read_csv(ld_file, sep=r'\s+')
            return ld_df
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"PLINK2 not available or failed: {e}")

    return pd.DataFrame()


def query_ldlink_api(
    snp: str,
    population: str,
    r2_threshold: float,
    window: int = 500000
) -> pd.DataFrame:
    """
    Query LDlink API for LD proxies (requires internet connection).
    Note: Rate limited, use sparingly.
    """
    import requests
    import time

    base_url = "https://ldlink.nih.gov/LDlinkRest/ldproxy"
    params = {
        "var": snp,
        "pop": population,
        "r2_d": "r2",
        "window": window,
        "genome_build": "grch37",
        "token": "your_ldlink_token"  # Get from https://ldlink.nih.gov/?tab=apiaccess
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                from io import StringIO
                df = pd.read_csv(StringIO(response.text), sep='\t')
                df = df[df['R2'] >= r2_threshold]
                time.sleep(1)  # Rate limiting
                return df
    except Exception as e:
        logger.warning(f"LDlink API query failed for {snp}: {e}")

    return pd.DataFrame()


def get_ld_proxies_cyvcf2(
    lead_snps: pd.DataFrame,
    vcf_path: str,
    population: str,
    r2_threshold: float,
    window_kb: int
) -> pd.DataFrame:
    """
    Calculate LD using cyvcf2 for direct VCF parsing.
    More memory efficient for large files.
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        logger.warning("cyvcf2 not installed. Install with: pip install cyvcf2")
        return pd.DataFrame()

    all_proxies = []

    for _, lead in lead_snps.iterrows():
        chrom = str(lead['chromosome'])
        pos = int(lead['position'])
        rsid = lead.get('rsid', f"chr{chrom}:{pos}")

        vcf_file = Path(vcf_path) / f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

        if not vcf_file.exists():
            logger.warning(f"VCF not found for chr{chrom}, skipping")
            continue

        try:
            vcf = VCF(str(vcf_file))

            # Define region
            start = max(1, pos - window_kb * 1000)
            end = pos + window_kb * 1000
            region = f"{chrom}:{start}-{end}"

            lead_geno = None
            proxies = []

            for variant in vcf(region):
                geno = np.array([sum(g[:2]) for g in variant.genotypes])

                if variant.POS == pos:
                    lead_geno = geno
                    continue

                proxies.append({
                    'proxy_rsid': variant.ID if variant.ID else f"chr{chrom}:{variant.POS}",
                    'proxy_chrom': chrom,
                    'proxy_pos': variant.POS,
                    'proxy_ref': variant.REF,
                    'proxy_alt': variant.ALT[0] if variant.ALT else '',
                    'genotypes': geno
                })

            # Calculate R2 with lead SNP
            if lead_geno is not None:
                for proxy in proxies:
                    proxy_geno = proxy.pop('genotypes')
                    r2 = np.corrcoef(lead_geno, proxy_geno)[0, 1] ** 2
                    if r2 >= r2_threshold:
                        proxy['r2'] = r2
                        proxy['lead_rsid'] = rsid
                        proxy['lead_pos'] = pos
                        all_proxies.append(proxy)

            vcf.close()

        except Exception as e:
            logger.error(f"Error processing chr{chrom}: {e}")

    if all_proxies:
        return pd.DataFrame(all_proxies)
    return pd.DataFrame()


def expand_with_proxies(lead_snps: pd.DataFrame, ld_proxies: pd.DataFrame) -> pd.DataFrame:
    """Combine lead SNPs with their LD proxies."""
    # Mark lead SNPs
    leads = lead_snps.copy()
    leads['is_proxy'] = False
    leads['r2'] = 1.0
    leads['lead_rsid'] = leads['rsid']

    if ld_proxies.empty:
        return leads

    # Process proxies
    proxies = ld_proxies.copy()
    proxies['is_proxy'] = True
    proxies['chromosome'] = proxies['proxy_chrom']
    proxies['position'] = proxies['proxy_pos']
    proxies['rsid'] = proxies['proxy_rsid']

    # Combine
    combined = pd.concat([leads, proxies], ignore_index=True)
    combined = combined.drop_duplicates(subset=['chromosome', 'position'])

    return combined


def main(config_path: str = "config/config.yaml"):
    """Main LD proxy identification workflow."""
    config = load_config(config_path)

    # Load lead SNPs from previous step
    prefix = config['output']['prefix']
    intermediate_dir = Path("data/intermediate")

    lead_file = intermediate_dir / f"{prefix}_lead_snps.tsv"
    if not lead_file.exists():
        logger.error(f"Lead SNPs file not found: {lead_file}")
        logger.error("Run 01_extract_variants.py first")
        return

    lead_snps = pd.read_csv(lead_file, sep='\t')
    logger.info(f"Loaded {len(lead_snps)} lead SNPs")

    # LD configuration
    ld_config = config['ld']
    vcf_path = ld_config['vcf_path']
    population = ld_config['population']
    r2_threshold = ld_config['r2_threshold']
    window_kb = ld_config['window_kb']

    # Get LD proxies
    logger.info(f"Finding LD proxies (r2 >= {r2_threshold}, window = {window_kb}kb)")

    ld_proxies = get_ld_proxies_cyvcf2(
        lead_snps,
        vcf_path,
        population,
        r2_threshold,
        window_kb
    )

    if ld_proxies.empty:
        logger.warning("No LD proxies found via cyvcf2, trying LDlink API")
        # Fallback to LDlink (rate limited)
        proxy_list = []
        for _, snp in lead_snps.iterrows():
            rsid = snp.get('rsid', '')
            if rsid.startswith('rs'):
                proxies = query_ldlink_api(rsid, population, r2_threshold)
                if not proxies.empty:
                    proxies['lead_rsid'] = rsid
                    proxy_list.append(proxies)

        if proxy_list:
            ld_proxies = pd.concat(proxy_list, ignore_index=True)

    # Combine lead SNPs with proxies
    result = expand_with_proxies(lead_snps, ld_proxies)

    # Save results
    output_file = intermediate_dir / f"{prefix}_variants_with_ld.tsv"
    result.to_csv(output_file, sep='\t', index=False)

    n_leads = len(lead_snps)
    n_proxies = len(result) - n_leads
    logger.info(f"Saved {len(result)} variants ({n_leads} leads + {n_proxies} proxies) to {output_file}")

    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Identify LD proxies for lead SNPs")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
