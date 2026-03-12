#!/usr/bin/env python3
"""
Variant annotation module for AlphaGWAS pipeline.

Integrates external annotations:
- ClinVar pathogenicity
- gnomAD allele frequencies
- GTEx eQTL associations
- CADD scores

All annotations are fetched via public APIs or local databases.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import requests
import time
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass
import yaml
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import utilities
try:
    from .utils import retry_with_backoff, progress_iterator, ProgressTracker
except ImportError:
    from utils import retry_with_backoff, progress_iterator, ProgressTracker

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass
class VariantAnnotation:
    """Container for variant annotations."""
    variant_id: str
    clinvar_significance: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    gnomad_af: Optional[float] = None
    gnomad_af_popmax: Optional[float] = None
    gtex_egenes: Optional[List[str]] = None
    gtex_tissues: Optional[List[str]] = None
    cadd_phred: Optional[float] = None
    cadd_raw: Optional[float] = None


class ClinVarAnnotator:
    """
    Fetch ClinVar pathogenicity annotations.

    Uses NCBI E-utilities API.
    """

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Pathogenicity levels for scoring
    PATHOGENICITY_SCORES = {
        'Pathogenic': 1.0,
        'Likely pathogenic': 0.8,
        'Pathogenic/Likely pathogenic': 0.9,
        'Uncertain significance': 0.3,
        'Likely benign': 0.1,
        'Benign': 0.0,
        'Benign/Likely benign': 0.05
    }

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize ClinVar annotator.

        Args:
            api_key: NCBI API key (increases rate limit)
        """
        self.api_key = api_key
        self.session = requests.Session()
        self.cache = {}

    @retry_with_backoff(max_retries=3, base_delay=0.5, exceptions=(requests.RequestException,))
    def fetch_by_rsid(self, rsid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch ClinVar annotation for an rsID.

        Args:
            rsid: dbSNP rsID (e.g., 'rs12345')

        Returns:
            Dictionary with ClinVar annotation or None
        """
        if rsid in self.cache:
            return self.cache[rsid]

        if not rsid.startswith('rs'):
            return None

        # Search ClinVar for the rsID
        search_url = f"{self.BASE_URL}/esearch.fcgi"
        params = {
            'db': 'clinvar',
            'term': f'{rsid}[variant name]',
            'retmode': 'json'
        }
        if self.api_key:
            params['api_key'] = self.api_key

        try:
            response = self.session.get(search_url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()

            id_list = data.get('esearchresult', {}).get('idlist', [])
            if not id_list:
                self.cache[rsid] = None
                return None

            # Fetch details for first result
            fetch_url = f"{self.BASE_URL}/esummary.fcgi"
            params = {
                'db': 'clinvar',
                'id': id_list[0],
                'retmode': 'json'
            }
            if self.api_key:
                params['api_key'] = self.api_key

            response = self.session.get(fetch_url, params=params, timeout=10)
            response.raise_for_status()
            summary = response.json()

            result = summary.get('result', {}).get(id_list[0], {})

            annotation = {
                'clinvar_significance': result.get('clinical_significance', {}).get('description', 'Unknown'),
                'clinvar_review_status': result.get('clinical_significance', {}).get('review_status', 'Unknown'),
                'clinvar_id': id_list[0]
            }

            self.cache[rsid] = annotation
            time.sleep(0.34)  # Rate limit: 3 requests/second without API key

            return annotation

        except Exception as e:
            logger.debug(f"ClinVar lookup failed for {rsid}: {e}")
            return None

    def get_pathogenicity_score(self, significance: str) -> float:
        """Convert ClinVar significance to numeric score."""
        if not significance:
            return 0.0
        for key, score in self.PATHOGENICITY_SCORES.items():
            if key.lower() in significance.lower():
                return score
        return 0.0


class GnomADAnnotator:
    """
    Fetch gnomAD allele frequency annotations.

    Uses gnomAD API (GraphQL).
    """

    API_URL = "https://gnomad.broadinstitute.org/api"

    def __init__(self, dataset: str = "gnomad_r4"):
        """
        Initialize gnomAD annotator.

        Args:
            dataset: gnomAD dataset version
        """
        self.dataset = dataset
        self.session = requests.Session()
        self.cache = {}

    @retry_with_backoff(max_retries=3, base_delay=1.0, exceptions=(requests.RequestException,))
    def fetch_by_position(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch gnomAD annotation by genomic position.

        Args:
            chrom: Chromosome (1-22, X, Y)
            pos: Position (1-based)
            ref: Reference allele
            alt: Alternate allele

        Returns:
            Dictionary with gnomAD annotation or None
        """
        variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        if variant_id in self.cache:
            return self.cache[variant_id]

        query = """
        query VariantQuery($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                variant_id
                genome {
                    ac
                    an
                    af
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
                exome {
                    ac
                    an
                    af
                }
            }
        }
        """

        variables = {
            "variantId": variant_id,
            "dataset": self.dataset
        }

        try:
            response = self.session.post(
                self.API_URL,
                json={"query": query, "variables": variables},
                timeout=15
            )
            response.raise_for_status()
            data = response.json()

            variant_data = data.get('data', {}).get('variant')
            if not variant_data:
                self.cache[variant_id] = None
                return None

            # Extract allele frequencies
            genome = variant_data.get('genome', {}) or {}
            exome = variant_data.get('exome', {}) or {}

            # Calculate combined AF
            genome_af = genome.get('af', 0) or 0
            exome_af = exome.get('af', 0) or 0
            combined_af = max(genome_af, exome_af)

            # Get population max AF
            popmax_af = 0
            populations = genome.get('populations', []) or []
            for pop in populations:
                pop_af = pop.get('af', 0) or 0
                if pop_af > popmax_af:
                    popmax_af = pop_af

            annotation = {
                'gnomad_af': combined_af,
                'gnomad_af_genome': genome_af,
                'gnomad_af_exome': exome_af,
                'gnomad_af_popmax': popmax_af,
                'gnomad_ac': (genome.get('ac', 0) or 0) + (exome.get('ac', 0) or 0),
                'gnomad_an': (genome.get('an', 0) or 0) + (exome.get('an', 0) or 0)
            }

            self.cache[variant_id] = annotation
            return annotation

        except Exception as e:
            logger.debug(f"gnomAD lookup failed for {variant_id}: {e}")
            return None


class GTExAnnotator:
    """
    Fetch GTEx eQTL annotations.

    Uses GTEx Portal API.
    """

    API_URL = "https://gtexportal.org/api/v2"

    def __init__(self):
        """Initialize GTEx annotator."""
        self.session = requests.Session()
        self.cache = {}

    @retry_with_backoff(max_retries=3, base_delay=1.0, exceptions=(requests.RequestException,))
    def fetch_eqtls(self, rsid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch GTEx eQTL associations for a variant.

        Args:
            rsid: dbSNP rsID

        Returns:
            Dictionary with GTEx eQTL annotation or None
        """
        if rsid in self.cache:
            return self.cache[rsid]

        if not rsid.startswith('rs'):
            return None

        url = f"{self.API_URL}/association/singleTissueEqtl"
        params = {
            'snpId': rsid,
            'datasetId': 'gtex_v8'
        }

        try:
            response = self.session.get(url, params=params, timeout=15)
            response.raise_for_status()
            data = response.json()

            associations = data.get('data', [])
            if not associations:
                self.cache[rsid] = None
                return None

            # Extract significant eQTLs (p < 0.05)
            significant = [a for a in associations if a.get('pValue', 1) < 0.05]

            if not significant:
                self.cache[rsid] = None
                return None

            # Get unique genes and tissues
            egenes = list(set(a.get('geneSymbol', '') for a in significant if a.get('geneSymbol')))
            tissues = list(set(a.get('tissueSiteDetailId', '') for a in significant if a.get('tissueSiteDetailId')))

            # Get best eQTL
            best_eqtl = min(significant, key=lambda x: x.get('pValue', 1))

            annotation = {
                'gtex_egenes': egenes,
                'gtex_tissues': tissues,
                'gtex_n_eqtls': len(significant),
                'gtex_best_pvalue': best_eqtl.get('pValue'),
                'gtex_best_gene': best_eqtl.get('geneSymbol'),
                'gtex_best_tissue': best_eqtl.get('tissueSiteDetailId')
            }

            self.cache[rsid] = annotation
            return annotation

        except Exception as e:
            logger.debug(f"GTEx lookup failed for {rsid}: {e}")
            return None


class VariantAnnotator:
    """
    Main class for variant annotation.

    Combines annotations from multiple sources.
    """

    def __init__(
        self,
        config: Dict[str, Any],
        enable_clinvar: bool = True,
        enable_gnomad: bool = True,
        enable_gtex: bool = True
    ):
        """
        Initialize variant annotator.

        Args:
            config: Configuration dictionary
            enable_clinvar: Enable ClinVar annotations
            enable_gnomad: Enable gnomAD annotations
            enable_gtex: Enable GTEx annotations
        """
        self.config = config.get('annotations', {})
        self.max_workers = self.config.get('max_workers', 2)

        self.clinvar = ClinVarAnnotator() if enable_clinvar else None
        self.gnomad = GnomADAnnotator() if enable_gnomad else None
        self.gtex = GTExAnnotator() if enable_gtex else None

    def annotate_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        rsid: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Annotate a single variant.

        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            rsid: Optional rsID

        Returns:
            Dictionary with all annotations
        """
        annotations = {
            'variant_id': f"chr{chrom}:{pos}:{ref}:{alt}",
            'rsid': rsid
        }

        # ClinVar (requires rsID)
        if self.clinvar and rsid:
            clinvar_data = self.clinvar.fetch_by_rsid(rsid)
            if clinvar_data:
                annotations.update(clinvar_data)
                annotations['clinvar_score'] = self.clinvar.get_pathogenicity_score(
                    clinvar_data.get('clinvar_significance', '')
                )

        # gnomAD
        if self.gnomad:
            gnomad_data = self.gnomad.fetch_by_position(chrom, pos, ref, alt)
            if gnomad_data:
                annotations.update(gnomad_data)

        # GTEx (requires rsID)
        if self.gtex and rsid:
            gtex_data = self.gtex.fetch_eqtls(rsid)
            if gtex_data:
                annotations.update(gtex_data)

        return annotations

    def annotate_dataframe(
        self,
        df: pd.DataFrame,
        chrom_col: str = 'chromosome',
        pos_col: str = 'position',
        ref_col: str = 'ref',
        alt_col: str = 'alt',
        rsid_col: str = 'rsid'
    ) -> pd.DataFrame:
        """
        Annotate all variants in a DataFrame.

        Args:
            df: Input DataFrame
            chrom_col: Chromosome column name
            pos_col: Position column name
            ref_col: Reference allele column name
            alt_col: Alternate allele column name
            rsid_col: rsID column name

        Returns:
            DataFrame with annotation columns added
        """
        logger.info(f"Annotating {len(df)} variants...")

        annotations = []

        with ProgressTracker(total=len(df), desc="Annotating variants") as tracker:
            for idx, row in df.iterrows():
                chrom = str(row.get(chrom_col, ''))
                pos = int(row.get(pos_col, 0))
                ref = str(row.get(ref_col, row.get('other_allele', 'N')))
                alt = str(row.get(alt_col, row.get('effect_allele', 'N')))
                rsid = row.get(rsid_col) if rsid_col in df.columns else None

                annotation = self.annotate_variant(chrom, pos, ref, alt, rsid)
                annotations.append(annotation)
                tracker.advance()

        # Convert to DataFrame and merge
        ann_df = pd.DataFrame(annotations)

        # Merge on variant_id
        df = df.copy()
        df['_variant_id'] = 'chr' + df[chrom_col].astype(str) + ':' + df[pos_col].astype(str) + ':' + \
                           df.get(ref_col, df.get('other_allele', 'N')).astype(str) + ':' + \
                           df.get(alt_col, df.get('effect_allele', 'N')).astype(str)

        # Add annotation columns
        ann_cols = [c for c in ann_df.columns if c not in ['variant_id', 'rsid']]
        for col in ann_cols:
            if col in ann_df.columns:
                df[col] = ann_df[col].values

        df = df.drop(columns=['_variant_id'], errors='ignore')

        return df

    def calculate_annotation_score(self, row: pd.Series) -> float:
        """
        Calculate a combined annotation score for variant prioritization.

        Higher scores indicate more likely functional/pathogenic variants.

        Args:
            row: Series with annotation columns

        Returns:
            Combined annotation score (0-1)
        """
        score = 0.0
        weights = {
            'clinvar': 0.4,
            'gnomad': 0.3,
            'gtex': 0.3
        }

        # ClinVar contribution
        clinvar_score = row.get('clinvar_score', 0)
        score += weights['clinvar'] * clinvar_score

        # gnomAD contribution (rare variants score higher)
        gnomad_af = row.get('gnomad_af', 0)
        if gnomad_af is not None and gnomad_af > 0:
            # Log-scale rarity score
            rarity_score = max(0, 1 - np.log10(gnomad_af + 1e-10) / np.log10(1e-6))
            rarity_score = min(1, rarity_score)
        else:
            rarity_score = 1.0  # Novel variants get high score
        score += weights['gnomad'] * rarity_score

        # GTEx contribution
        gtex_n = row.get('gtex_n_eqtls', 0) or 0
        gtex_score = min(1.0, gtex_n / 10)  # Cap at 10 eQTLs
        score += weights['gtex'] * gtex_score

        return score


def main(config_path: str = "config/config.yaml"):
    """Main annotation workflow."""
    config = load_config(config_path)

    prefix = config['output']['prefix']
    intermediate_dir = Path("data/intermediate")
    output_dir = Path(config['output']['dir'])

    # Load variants
    input_files = [
        output_dir / f"{prefix}_ranked_variants.tsv",
        intermediate_dir / f"{prefix}_variants_hg38.tsv",
        intermediate_dir / f"{prefix}_variants_with_ld.tsv"
    ]

    variants = None
    for input_file in input_files:
        if input_file.exists():
            variants = pd.read_csv(input_file, sep='\t')
            logger.info(f"Loaded {len(variants)} variants from {input_file}")
            break

    if variants is None:
        logger.error("No variant file found")
        return

    # Initialize annotator
    annotator = VariantAnnotator(config)

    # Annotate variants
    annotated = annotator.annotate_dataframe(variants)

    # Calculate annotation scores
    annotated['annotation_score'] = annotated.apply(
        annotator.calculate_annotation_score, axis=1
    )

    # Save results
    output_file = output_dir / f"{prefix}_annotated_variants.tsv"
    annotated.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved annotated variants to {output_file}")

    # Print summary
    n_clinvar = annotated['clinvar_significance'].notna().sum()
    n_gnomad = annotated['gnomad_af'].notna().sum()
    n_gtex = annotated['gtex_n_eqtls'].notna().sum()

    print("\n" + "=" * 60)
    print("ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"Total variants: {len(annotated)}")
    print(f"ClinVar annotations: {n_clinvar}")
    print(f"gnomAD annotations: {n_gnomad}")
    print(f"GTEx eQTL annotations: {n_gtex}")
    print("=" * 60 + "\n")

    return annotated


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Annotate variants with external databases")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    parser.add_argument("--no-clinvar", action="store_true", help="Disable ClinVar annotations")
    parser.add_argument("--no-gnomad", action="store_true", help="Disable gnomAD annotations")
    parser.add_argument("--no-gtex", action="store_true", help="Disable GTEx annotations")
    args = parser.parse_args()

    main(args.config)
