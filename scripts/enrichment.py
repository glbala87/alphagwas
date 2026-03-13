#!/usr/bin/env python3
"""
Pathway enrichment analysis module for AlphaGWAS pipeline.

Performs gene set enrichment analysis on prioritized variants:
- GO (Gene Ontology) terms
- KEGG pathways
- Reactome pathways
- Custom gene sets

Uses Enrichr API for enrichment calculations.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import requests
import json
import time
from typing import Dict, List, Optional, Tuple, Any, Set
from dataclasses import dataclass
import yaml
from collections import defaultdict

# Import utilities
try:
    from .utils import retry_with_backoff, progress_iterator, ProgressTracker, print_summary_table
except ImportError:
    from utils import retry_with_backoff, progress_iterator, ProgressTracker, print_summary_table

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


@dataclass
class EnrichmentResult:
    """Container for enrichment analysis result."""
    term: str
    term_id: str
    database: str
    p_value: float
    adjusted_p_value: float
    odds_ratio: float
    combined_score: float
    genes: List[str]
    gene_count: int
    background_count: int


class EnrichrClient:
    """
    Client for Enrichr gene set enrichment API.

    Enrichr is a web-based tool for gene set enrichment analysis.
    https://maayanlab.cloud/Enrichr/
    """

    BASE_URL = "https://maayanlab.cloud/Enrichr"

    # Available libraries
    LIBRARIES = {
        # Gene Ontology
        'GO_Biological_Process_2023': 'GO:BP',
        'GO_Molecular_Function_2023': 'GO:MF',
        'GO_Cellular_Component_2023': 'GO:CC',

        # Pathways
        'KEGG_2021_Human': 'KEGG',
        'Reactome_2022': 'Reactome',
        'WikiPathways_2023_Human': 'WikiPathways',
        'BioPlanet_2019': 'BioPlanet',

        # Disease/Phenotype
        'DisGeNET': 'DisGeNET',
        'GWAS_Catalog_2023': 'GWAS_Catalog',
        'OMIM_Disease': 'OMIM',

        # Transcription factors
        'ENCODE_TF_ChIP-seq_2015': 'ENCODE_TF',
        'ChEA_2022': 'ChEA',

        # Tissue expression
        'GTEx_Tissue_Expression_Up': 'GTEx_Up',
        'GTEx_Tissue_Expression_Down': 'GTEx_Down',

        # Drug/compound
        'DrugMatrix': 'DrugMatrix',
        'Drug_Perturbations_from_GEO_2014': 'Drug_GEO',
    }

    def __init__(self, timeout: int = 30):
        """Initialize Enrichr client."""
        self.timeout = timeout
        self.session = requests.Session()
        self._list_id = None

    @retry_with_backoff(max_retries=3, base_delay=1.0, exceptions=(requests.RequestException,))
    def submit_gene_list(self, genes: List[str], description: str = "") -> int:
        """
        Submit a gene list to Enrichr.

        Args:
            genes: List of gene symbols
            description: Optional description

        Returns:
            List ID for querying results
        """
        url = f"{self.BASE_URL}/addList"

        payload = {
            'list': '\n'.join(genes),
            'description': description or 'AlphaGWAS analysis'
        }

        response = self.session.post(url, data=payload, timeout=self.timeout)
        response.raise_for_status()

        data = response.json()
        self._list_id = data.get('userListId')

        if not self._list_id:
            raise ValueError("Failed to get list ID from Enrichr")

        logger.info(f"Submitted {len(genes)} genes to Enrichr (ID: {self._list_id})")
        return self._list_id

    @retry_with_backoff(max_retries=3, base_delay=1.0, exceptions=(requests.RequestException,))
    def get_enrichment(
        self,
        library: str,
        list_id: Optional[int] = None
    ) -> List[EnrichmentResult]:
        """
        Get enrichment results for a library.

        Args:
            library: Name of the gene set library
            list_id: Optional list ID (uses last submitted if not provided)

        Returns:
            List of EnrichmentResult objects
        """
        list_id = list_id or self._list_id

        if not list_id:
            raise ValueError("No list ID available. Submit genes first.")

        url = f"{self.BASE_URL}/enrich"
        params = {
            'userListId': list_id,
            'backgroundType': library
        }

        response = self.session.get(url, params=params, timeout=self.timeout)
        response.raise_for_status()

        data = response.json()
        results = []

        for item in data.get(library, []):
            # Enrichr returns: [rank, term, p-value, z-score, combined_score,
            #                   overlapping_genes, adjusted_p-value, old_p-value,
            #                   old_adjusted_p-value]
            if len(item) >= 7:
                term = item[1]
                p_value = item[2]
                z_score = item[3]
                combined_score = item[4]
                genes = item[5] if isinstance(item[5], list) else item[5].split(';')
                adj_p_value = item[6]

                # Extract term ID if present
                term_id = ""
                if '(' in term and ')' in term:
                    term_id = term[term.rfind('(')+1:term.rfind(')')]

                results.append(EnrichmentResult(
                    term=term,
                    term_id=term_id,
                    database=library,
                    p_value=p_value,
                    adjusted_p_value=adj_p_value,
                    odds_ratio=np.exp(z_score) if z_score else 1.0,
                    combined_score=combined_score,
                    genes=genes,
                    gene_count=len(genes),
                    background_count=0  # Not provided by Enrichr
                ))

        return results

    def run_enrichment(
        self,
        genes: List[str],
        libraries: Optional[List[str]] = None,
        description: str = ""
    ) -> Dict[str, List[EnrichmentResult]]:
        """
        Run enrichment analysis on multiple libraries.

        Args:
            genes: List of gene symbols
            libraries: List of library names (default: all main libraries)
            description: Analysis description

        Returns:
            Dictionary mapping library names to results
        """
        if not genes:
            logger.warning("No genes provided for enrichment")
            return {}

        # Default libraries
        if libraries is None:
            libraries = [
                'GO_Biological_Process_2023',
                'GO_Molecular_Function_2023',
                'KEGG_2021_Human',
                'Reactome_2022'
            ]

        # Submit gene list
        self.submit_gene_list(genes, description)

        # Get results for each library
        all_results = {}

        for library in libraries:
            try:
                results = self.get_enrichment(library)
                all_results[library] = results
                logger.info(f"{library}: {len(results)} enriched terms")
                time.sleep(0.5)  # Rate limiting
            except Exception as e:
                logger.warning(f"Failed to get results for {library}: {e}")
                all_results[library] = []

        return all_results


class LocalEnrichment:
    """
    Local gene set enrichment using Fisher's exact test.

    Use when Enrichr API is not available or for custom gene sets.
    """

    def __init__(self):
        """Initialize local enrichment."""
        self.gene_sets: Dict[str, Dict[str, Set[str]]] = {}
        self._load_builtin_genesets()

    def _load_builtin_genesets(self):
        """Load built-in minimal gene sets."""
        # Minimal example gene sets (in practice, load from files)
        self.gene_sets['cardiovascular'] = {
            'Heart_development': {'GATA4', 'NKX2-5', 'TBX5', 'MYH6', 'MYH7', 'TNNT2', 'TNNI3'},
            'Cardiac_conduction': {'SCN5A', 'KCNQ1', 'KCNH2', 'KCNE1', 'KCNJ2', 'HCN4'},
            'Blood_pressure': {'ACE', 'AGT', 'AGTR1', 'NOS3', 'EDN1', 'NPPA', 'NPPB'},
            'Lipid_metabolism': {'LDLR', 'APOB', 'PCSK9', 'APOE', 'CETP', 'LIPC', 'LPL'},
        }

        self.gene_sets['metabolic'] = {
            'Glucose_homeostasis': {'INS', 'GCK', 'GLUT4', 'IRS1', 'IRS2', 'PPARG', 'TCF7L2'},
            'Fatty_acid_metabolism': {'FASN', 'ACACA', 'CPT1A', 'ACADM', 'HADH'},
            'Cholesterol_synthesis': {'HMGCR', 'SQLE', 'DHCR7', 'SREBF2'},
        }

    def load_gmt_file(self, gmt_path: str, name: str):
        """
        Load gene sets from GMT file.

        GMT format: term_name<tab>description<tab>gene1<tab>gene2<tab>...
        """
        gene_sets = {}

        with open(gmt_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    term = parts[0]
                    genes = set(parts[2:])
                    gene_sets[term] = genes

        self.gene_sets[name] = gene_sets
        logger.info(f"Loaded {len(gene_sets)} gene sets from {gmt_path}")

    def fisher_test(
        self,
        query_genes: Set[str],
        gene_set: Set[str],
        background_size: int = 20000
    ) -> Tuple[float, float]:
        """
        Perform Fisher's exact test for enrichment.

        Args:
            query_genes: Set of query genes
            gene_set: Set of genes in the pathway
            background_size: Total number of genes in background

        Returns:
            Tuple of (odds_ratio, p_value)
        """
        from scipy import stats

        # Contingency table
        a = len(query_genes & gene_set)  # In query AND in pathway
        b = len(gene_set - query_genes)   # Not in query, in pathway
        c = len(query_genes - gene_set)   # In query, not in pathway
        d = background_size - a - b - c    # Not in query, not in pathway

        # Fisher's exact test
        contingency = [[a, b], [c, d]]
        odds_ratio, p_value = stats.fisher_exact(contingency, alternative='greater')

        return odds_ratio, p_value

    def run_enrichment(
        self,
        genes: List[str],
        gene_set_names: Optional[List[str]] = None,
        min_genes: int = 2,
        max_genes: int = 500
    ) -> List[EnrichmentResult]:
        """
        Run local enrichment analysis.

        Args:
            genes: List of query genes
            gene_set_names: Which gene set collections to use
            min_genes: Minimum overlap for testing
            max_genes: Maximum pathway size to test

        Returns:
            List of EnrichmentResult objects
        """
        query_genes = set(genes)
        results = []

        if gene_set_names is None:
            gene_set_names = list(self.gene_sets.keys())

        for gs_name in gene_set_names:
            if gs_name not in self.gene_sets:
                continue

            for term, pathway_genes in self.gene_sets[gs_name].items():
                if len(pathway_genes) < min_genes or len(pathway_genes) > max_genes:
                    continue

                overlap = query_genes & pathway_genes
                if len(overlap) < min_genes:
                    continue

                odds_ratio, p_value = self.fisher_test(query_genes, pathway_genes)

                results.append(EnrichmentResult(
                    term=term,
                    term_id="",
                    database=gs_name,
                    p_value=p_value,
                    adjusted_p_value=p_value,  # Will be adjusted later
                    odds_ratio=odds_ratio,
                    combined_score=odds_ratio * (-np.log10(p_value + 1e-300)),
                    genes=list(overlap),
                    gene_count=len(overlap),
                    background_count=len(pathway_genes)
                ))

        # Multiple testing correction (Benjamini-Hochberg)
        if results:
            results = self._adjust_pvalues(results)

        # Sort by adjusted p-value
        results.sort(key=lambda x: x.adjusted_p_value)

        return results

    def _adjust_pvalues(self, results: List[EnrichmentResult]) -> List[EnrichmentResult]:
        """Apply Benjamini-Hochberg FDR correction."""
        n = len(results)
        p_values = [r.p_value for r in results]

        # Sort indices by p-value
        sorted_indices = np.argsort(p_values)

        # Calculate adjusted p-values
        adjusted = [0.0] * n
        prev_adj = 1.0

        for i in range(n - 1, -1, -1):
            idx = sorted_indices[i]
            adj = min(prev_adj, p_values[idx] * n / (i + 1))
            adjusted[idx] = adj
            prev_adj = adj

        # Update results
        for i, r in enumerate(results):
            r.adjusted_p_value = adjusted[i]

        return results


class EnrichmentAnalyzer:
    """
    Main class for running enrichment analysis on AlphaGWAS results.
    """

    def __init__(self, config: Dict[str, Any]):
        """Initialize analyzer."""
        self.config = config.get('enrichment', {})
        self.use_enrichr = self.config.get('use_enrichr', True)
        self.p_threshold = self.config.get('p_threshold', 0.05)
        self.min_genes = self.config.get('min_genes', 3)

        if self.use_enrichr:
            self.enrichr = EnrichrClient()
        else:
            self.enrichr = None

        self.local = LocalEnrichment()

    def extract_genes(
        self,
        variants: pd.DataFrame,
        gene_col: str = 'nearby_genes',
        top_n: Optional[int] = None
    ) -> List[str]:
        """
        Extract unique genes from variant data.

        Args:
            variants: DataFrame with variant and gene information
            gene_col: Column containing gene names
            top_n: Only use top N variants

        Returns:
            List of unique gene symbols
        """
        if gene_col not in variants.columns:
            logger.warning(f"Gene column '{gene_col}' not found")
            return []

        if top_n:
            variants = variants.head(top_n)

        genes = set()

        for gene_str in variants[gene_col].dropna():
            if isinstance(gene_str, str):
                # Handle comma-separated or semicolon-separated genes
                for g in gene_str.replace(';', ',').split(','):
                    g = g.strip()
                    if g and g != 'nan':
                        genes.add(g)
            elif isinstance(gene_str, list):
                genes.update(gene_str)

        return list(genes)

    def run_analysis(
        self,
        variants: pd.DataFrame,
        gene_col: str = 'nearby_genes',
        top_n: int = 100,
        libraries: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Run full enrichment analysis.

        Args:
            variants: DataFrame with ranked variants
            gene_col: Column containing gene names
            top_n: Number of top variants to analyze
            libraries: Enrichr libraries to use

        Returns:
            DataFrame with enrichment results
        """
        # Extract genes
        genes = self.extract_genes(variants, gene_col, top_n)
        logger.info(f"Extracted {len(genes)} unique genes from top {top_n} variants")

        if len(genes) < self.min_genes:
            logger.warning(f"Too few genes ({len(genes)}) for enrichment analysis")
            return pd.DataFrame()

        all_results = []

        # Run Enrichr if available
        if self.enrichr:
            try:
                enrichr_results = self.enrichr.run_enrichment(
                    genes,
                    libraries=libraries,
                    description=f"AlphaGWAS top {top_n} variants"
                )

                for library, results in enrichr_results.items():
                    for r in results:
                        if r.adjusted_p_value <= self.p_threshold:
                            all_results.append({
                                'term': r.term,
                                'term_id': r.term_id,
                                'database': r.database,
                                'p_value': r.p_value,
                                'adjusted_p_value': r.adjusted_p_value,
                                'odds_ratio': r.odds_ratio,
                                'combined_score': r.combined_score,
                                'genes': ','.join(r.genes),
                                'gene_count': r.gene_count,
                                'source': 'Enrichr'
                            })

            except Exception as e:
                logger.warning(f"Enrichr analysis failed: {e}")

        # Run local enrichment
        local_results = self.local.run_enrichment(genes)
        for r in local_results:
            if r.adjusted_p_value <= self.p_threshold:
                all_results.append({
                    'term': r.term,
                    'term_id': r.term_id,
                    'database': r.database,
                    'p_value': r.p_value,
                    'adjusted_p_value': r.adjusted_p_value,
                    'odds_ratio': r.odds_ratio,
                    'combined_score': r.combined_score,
                    'genes': ','.join(r.genes),
                    'gene_count': r.gene_count,
                    'source': 'Local'
                })

        # Convert to DataFrame
        results_df = pd.DataFrame(all_results)

        if not results_df.empty:
            results_df = results_df.sort_values('adjusted_p_value')

        logger.info(f"Found {len(results_df)} significantly enriched terms")

        return results_df


def main(config_path: str = "config/config.yaml"):
    """Main enrichment analysis workflow."""
    config = load_config(config_path)

    prefix = config['output']['prefix']
    output_dir = Path(config['output']['dir'])

    # Load ranked variants
    ranked_file = output_dir / f"{prefix}_ranked_variants.tsv"

    if not ranked_file.exists():
        logger.error(f"Ranked variants file not found: {ranked_file}")
        return

    ranked_variants = pd.read_csv(ranked_file, sep='\t')
    logger.info(f"Loaded {len(ranked_variants)} ranked variants")

    # Initialize analyzer
    analyzer = EnrichmentAnalyzer(config)

    # Run enrichment on top variants
    top_n = config.get('enrichment', {}).get('top_n', 100)
    results = analyzer.run_analysis(ranked_variants, top_n=top_n)

    if results.empty:
        logger.info("No significant enrichments found")
        return

    # Save results
    output_file = output_dir / f"{prefix}_enrichment_results.tsv"
    results.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved enrichment results to {output_file}")

    # Print summary
    print("\n" + "=" * 70)
    print("PATHWAY ENRICHMENT SUMMARY")
    print("=" * 70)
    print(f"Input genes: from top {top_n} variants")
    print(f"Significant pathways (adj. p < 0.05): {len(results)}")
    print("\nTop 10 enriched terms:")
    print("-" * 70)

    display_cols = ['term', 'database', 'adjusted_p_value', 'gene_count']
    display_cols = [c for c in display_cols if c in results.columns]

    if display_cols:
        print(results[display_cols].head(10).to_string(index=False))

    print("=" * 70 + "\n")

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run pathway enrichment analysis")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    parser.add_argument("--top-n", type=int, default=100, help="Number of top variants to analyze")
    parser.add_argument("--no-enrichr", action="store_true", help="Skip Enrichr API (use local only)")
    args = parser.parse_args()

    main(args.config)
