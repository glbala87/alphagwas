#!/usr/bin/env python3
"""
Setup helper for AlphaGWAS pipeline.
Helps you configure the pipeline for your GWAS study.
"""

import pandas as pd
from pathlib import Path
import yaml
import sys


def inspect_gwas_file(file_path: str):
    """Inspect a GWAS summary statistics file and show column info."""
    print(f"\n{'='*60}")
    print(f"Inspecting: {file_path}")
    print(f"{'='*60}\n")

    # Detect file format
    if file_path.endswith('.gz'):
        df = pd.read_csv(file_path, sep=None, engine='python', nrows=5, compression='gzip')
    else:
        df = pd.read_csv(file_path, sep=None, engine='python', nrows=5)

    print("COLUMNS FOUND:")
    print("-" * 40)
    for i, col in enumerate(df.columns):
        sample_vals = df[col].head(3).tolist()
        print(f"  {i+1}. {col:20} | Examples: {sample_vals}")

    print(f"\n{'='*60}")
    print("COLUMN MAPPING GUIDE")
    print("="*60)
    print("""
Update your config.yaml 'columns' section to match your file:

gwas:
  columns:
    chromosome: "YOUR_CHR_COLUMN"    # e.g., CHR, #CHROM, chromosome
    position: "YOUR_POS_COLUMN"      # e.g., BP, POS, position
    rsid: "YOUR_SNP_COLUMN"          # e.g., SNP, rsid, ID
    effect_allele: "YOUR_A1_COLUMN"  # e.g., A1, ALT, effect_allele
    other_allele: "YOUR_A2_COLUMN"   # e.g., A2, REF, other_allele
    pvalue: "YOUR_P_COLUMN"          # e.g., P, PVALUE, pval
    beta: "YOUR_BETA_COLUMN"         # e.g., BETA, beta, effect (optional)
    se: "YOUR_SE_COLUMN"             # e.g., SE, stderr (optional)
""")

    # Try to auto-detect columns
    print("\nAUTO-DETECTED MAPPINGS (verify these!):")
    print("-" * 40)

    col_lower = {c.lower(): c for c in df.columns}

    mappings = {
        'chromosome': ['chr', 'chrom', 'chromosome', '#chrom'],
        'position': ['bp', 'pos', 'position', 'base_pair_location'],
        'rsid': ['snp', 'rsid', 'id', 'variant_id', 'snpid'],
        'effect_allele': ['a1', 'alt', 'effect_allele', 'ea', 'allele1'],
        'other_allele': ['a2', 'ref', 'other_allele', 'oa', 'allele2', 'nea'],
        'pvalue': ['p', 'pvalue', 'pval', 'p_value', 'p-value'],
        'beta': ['beta', 'effect', 'b', 'effect_size'],
        'se': ['se', 'stderr', 'standard_error'],
    }

    detected = {}
    for field, patterns in mappings.items():
        for pattern in patterns:
            if pattern in col_lower:
                detected[field] = col_lower[pattern]
                break

    for field, col in detected.items():
        print(f"  {field:15} -> {col}")

    missing = set(mappings.keys()) - set(detected.keys())
    if missing:
        print(f"\n  Could not auto-detect: {missing}")
        print("  Please set these manually in your config file.")

    return detected


def download_liftover_chain():
    """Download UCSC liftover chain file."""
    import urllib.request

    chain_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    output_path = Path("data/input/hg19ToHg38.over.chain.gz")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        print(f"Chain file already exists: {output_path}")
        return

    print(f"Downloading liftover chain from UCSC...")
    try:
        urllib.request.urlretrieve(chain_url, output_path)
        print(f"Downloaded to: {output_path}")
    except Exception as e:
        print(f"Download failed: {e}")
        print(f"Manual download: {chain_url}")


def create_config_from_template(
    gwas_file: str,
    phenotype: str,
    output_name: str = None
):
    """Create a customized config file."""
    # Load template
    template_path = Path("config/cardiovascular_config.yaml")
    if not template_path.exists():
        template_path = Path("config/config.yaml")

    with open(template_path) as f:
        config = yaml.safe_load(f)

    # Update with user values
    config['study']['phenotype'] = phenotype
    config['study']['name'] = f"{phenotype}_gwas"
    config['gwas']['input_file'] = gwas_file

    # Detect columns
    detected = inspect_gwas_file(gwas_file)
    if detected:
        config['gwas']['columns'] = detected

    # Save new config
    output_name = output_name or f"{phenotype}_config.yaml"
    output_path = Path("config") / output_name

    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"\nConfig saved to: {output_path}")
    print("Please review and adjust the column mappings if needed.")

    return output_path


def check_setup():
    """Check if everything is set up correctly."""
    print("\n" + "="*60)
    print("SETUP CHECKLIST")
    print("="*60)

    checks = []

    # Check directories
    dirs = ["data/input", "data/intermediate", "data/output", "config"]
    for d in dirs:
        exists = Path(d).exists()
        status = "[OK]" if exists else "[MISSING]"
        print(f"  {status} Directory: {d}")
        checks.append(exists)

    # Check for chain file
    chain = Path("data/input/hg19ToHg38.over.chain.gz")
    status = "[OK]" if chain.exists() else "[MISSING]"
    print(f"  {status} Liftover chain: {chain}")
    if not chain.exists():
        print("        Run: python setup_study.py --download-chain")
    checks.append(chain.exists())

    # Check Python packages
    packages = ['pandas', 'numpy', 'yaml', 'scipy']
    for pkg in packages:
        try:
            __import__(pkg)
            print(f"  [OK] Package: {pkg}")
            checks.append(True)
        except ImportError:
            print(f"  [MISSING] Package: {pkg}")
            checks.append(False)

    # Optional packages
    optional = ['pyliftover', 'cyvcf2', 'pyarrow']
    for pkg in optional:
        try:
            __import__(pkg)
            print(f"  [OK] Optional: {pkg}")
        except ImportError:
            print(f"  [OPTIONAL] Package: {pkg} (some features limited)")

    print("\n" + "="*60)
    if all(checks):
        print("All required components are ready!")
    else:
        print("Some components are missing. Please install/configure them.")
    print("="*60)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Setup helper for AlphaGWAS pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check if everything is set up
  python setup_study.py --check

  # Inspect your GWAS file to see columns
  python setup_study.py --inspect my_gwas.tsv

  # Create a new config from your GWAS file
  python setup_study.py --gwas my_gwas.tsv --phenotype blood_pressure

  # Download the liftover chain file
  python setup_study.py --download-chain
        """
    )

    parser.add_argument("--check", action="store_true", help="Check setup status")
    parser.add_argument("--inspect", metavar="FILE", help="Inspect a GWAS sumstats file")
    parser.add_argument("--gwas", metavar="FILE", help="Path to GWAS sumstats file")
    parser.add_argument("--phenotype", default="trait", help="Phenotype name")
    parser.add_argument("--download-chain", action="store_true", help="Download liftover chain")

    args = parser.parse_args()

    if args.check:
        check_setup()
    elif args.download_chain:
        download_liftover_chain()
    elif args.inspect:
        inspect_gwas_file(args.inspect)
    elif args.gwas:
        create_config_from_template(args.gwas, args.phenotype)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
