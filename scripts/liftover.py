#!/usr/bin/env python3
"""
Step 3: Liftover genomic coordinates between reference builds (hg19 <-> hg38).
"""

import pandas as pd
from pathlib import Path
import logging
import yaml
from typing import Optional, Tuple

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


class LiftoverConverter:
    """Handle coordinate conversion between genome builds."""

    def __init__(self, chain_file: str):
        """
        Initialize with a UCSC chain file.

        Chain files can be downloaded from:
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
        """
        self.chain_file = chain_file
        self.converter = None
        self._load_chain()

    def _load_chain(self):
        """Load the chain file using pyliftover."""
        try:
            from pyliftover import LiftOver
            self.converter = LiftOver(self.chain_file)
            logger.info(f"Loaded liftover chain: {self.chain_file}")
        except ImportError:
            logger.error("pyliftover not installed. Install with: pip install pyliftover")
            raise
        except Exception as e:
            logger.error(f"Failed to load chain file: {e}")
            raise

    def convert_coordinate(
        self,
        chrom: str,
        pos: int
    ) -> Optional[Tuple[str, int]]:
        """
        Convert a single coordinate.

        Returns (new_chrom, new_pos) or None if conversion fails.
        """
        if self.converter is None:
            return None

        # Ensure chromosome format matches chain file (usually 'chr1' format)
        if not chrom.startswith('chr'):
            chrom = f'chr{chrom}'

        try:
            result = self.converter.convert_coordinate(chrom, pos)
            if result and len(result) > 0:
                new_chrom, new_pos, strand, _ = result[0]
                # Remove 'chr' prefix if needed
                new_chrom = new_chrom.replace('chr', '')
                return (new_chrom, new_pos)
        except Exception as e:
            logger.debug(f"Conversion failed for {chrom}:{pos}: {e}")

        return None


def download_chain_file(output_path: str, source_build: str = "hg19", target_build: str = "hg38"):
    """Download UCSC chain file if not present."""
    import urllib.request
    import gzip

    url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{source_build}/liftOver/{source_build}To{target_build.capitalize()}.over.chain.gz"

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        logger.info(f"Chain file already exists: {output_path}")
        return

    logger.info(f"Downloading chain file from {url}")

    try:
        urllib.request.urlretrieve(url, output_path)
        logger.info(f"Downloaded chain file to {output_path}")
    except Exception as e:
        logger.error(f"Failed to download chain file: {e}")
        raise


def liftover_variants(
    df: pd.DataFrame,
    converter: LiftoverConverter,
    chrom_col: str = 'chromosome',
    pos_col: str = 'position'
) -> pd.DataFrame:
    """
    Liftover all variants in a dataframe.

    Adds new columns: lifted_chrom, lifted_pos, liftover_success
    """
    df = df.copy()

    lifted_chroms = []
    lifted_positions = []
    success = []

    total = len(df)
    for i, (_, row) in enumerate(df.iterrows()):
        if (i + 1) % 1000 == 0:
            logger.info(f"Progress: {i + 1}/{total} variants")

        chrom = str(row[chrom_col])
        pos = int(row[pos_col])

        result = converter.convert_coordinate(chrom, pos)

        if result:
            lifted_chroms.append(result[0])
            lifted_positions.append(result[1])
            success.append(True)
        else:
            lifted_chroms.append(None)
            lifted_positions.append(None)
            success.append(False)

    df['lifted_chrom'] = lifted_chroms
    df['lifted_pos'] = lifted_positions
    df['liftover_success'] = success

    success_count = sum(success)
    fail_count = total - success_count
    logger.info(f"Liftover complete: {success_count}/{total} successful ({fail_count} failed)")

    return df


def main(config_path: str = "config/config.yaml"):
    """Main liftover workflow."""
    config = load_config(config_path)

    # Check if liftover is needed
    ref_config = config['reference']
    input_build = ref_config['input_build']
    output_build = ref_config['output_build']

    if input_build == output_build:
        logger.info(f"Input and output builds are the same ({input_build}), skipping liftover")
        return

    logger.info(f"Converting coordinates from {input_build} to {output_build}")

    # Get chain file
    chain_file = ref_config.get('liftover_chain')
    if not chain_file or not Path(chain_file).exists():
        # Download chain file
        chain_file = f"data/input/{input_build}To{output_build.capitalize()}.over.chain.gz"
        download_chain_file(chain_file, input_build, output_build)

    # Initialize converter
    converter = LiftoverConverter(chain_file)

    # Load variants from previous step
    prefix = config['output']['prefix']
    intermediate_dir = Path("data/intermediate")

    input_file = intermediate_dir / f"{prefix}_variants_with_ld.tsv"
    if not input_file.exists():
        # Fall back to just significant variants
        input_file = intermediate_dir / f"{prefix}_significant_variants.tsv"

    if not input_file.exists():
        logger.error(f"No input file found. Run previous steps first.")
        return

    variants = pd.read_csv(input_file, sep='\t')
    logger.info(f"Loaded {len(variants)} variants for liftover")

    # Perform liftover
    result = liftover_variants(variants, converter)

    # Update coordinate columns for successful liftovers
    result.loc[result['liftover_success'], 'original_pos'] = result.loc[result['liftover_success'], 'position']
    result.loc[result['liftover_success'], 'position'] = result.loc[result['liftover_success'], 'lifted_pos']
    result['genome_build'] = output_build

    # Filter failed liftovers (optional - keep for troubleshooting)
    failed = result[~result['liftover_success']]
    if len(failed) > 0:
        failed_file = intermediate_dir / f"{prefix}_liftover_failed.tsv"
        failed.to_csv(failed_file, sep='\t', index=False)
        logger.warning(f"Saved {len(failed)} failed liftovers to {failed_file}")

    # Keep only successful liftovers for downstream analysis
    result = result[result['liftover_success']]

    # Save results
    output_file = intermediate_dir / f"{prefix}_variants_{output_build}.tsv"
    result.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(result)} variants ({output_build}) to {output_file}")

    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Liftover variant coordinates")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
