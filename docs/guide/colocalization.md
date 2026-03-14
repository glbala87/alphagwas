# Colocalization Analysis

Test whether two GWAS traits share a common causal variant at a locus.

## Overview

AlphaGWAS implements several colocalization methods:

- **COLOC**: Bayesian colocalization (Giambartolomei et al., 2014)
- **eCAVIAR**: Uses fine-mapping PIPs (Hormozdiari et al., 2016)
- **Simple**: Correlation-based quick test

## COLOC Method

COLOC tests five hypotheses:

| Hypothesis | Description |
|------------|-------------|
| H0 | No association with either trait |
| H1 | Association with trait 1 only |
| H2 | Association with trait 2 only |
| H3 | Association with both, different causal variants |
| H4 | Association with both, shared causal variant |

**PP.H4 > 0.8** indicates strong evidence for colocalization.

## Usage

### Command Line

```bash
python -c "
from scripts.colocalization import run_colocalization
from pathlib import Path

run_colocalization(
    gwas1_path=Path('trait1_sumstats.tsv'),
    gwas2_path=Path('trait2_sumstats.tsv'),
    output_path=Path('coloc_results.tsv'),
    method='coloc',
    trait1_name='LDL',
    trait2_name='CAD'
)
"
```

### Python API

```python
from scripts.colocalization import ColocAnalyzer, run_colocalization
import pandas as pd

# Load GWAS data
gwas1 = pd.read_csv("ldl_gwas.tsv", sep="\t")
gwas2 = pd.read_csv("cad_gwas.tsv", sep="\t")

# Initialize analyzer
analyzer = ColocAnalyzer(
    p1=1e-4,   # Prior for trait 1 association
    p2=1e-4,   # Prior for trait 2 association
    p12=1e-5   # Prior for shared causal variant
)

# Run for single locus
result = analyzer.coloc_analysis(
    trait1_df=gwas1_locus,
    trait2_df=gwas2_locus,
    trait1_name="LDL",
    trait2_name="CAD",
    locus_id="chr19_45000000"
)

print(f"PP.H4 (colocalization): {result.pp_h4:.3f}")
print(f"Top variant: {result.top_coloc_variant}")
print(f"Summary: {result.summary}")
```

### Multi-locus Analysis

```python
# Define loci
loci = pd.DataFrame({
    'locus_id': ['APOE', 'LPA', 'PCSK9'],
    'chromosome': ['19', '6', '1'],
    'start': [44900000, 160500000, 55000000],
    'end': [45500000, 161500000, 56000000]
})

# Run across all loci
results = analyzer.run_multi_locus(
    gwas1, gwas2, loci,
    trait1_name="LDL",
    trait2_name="CAD"
)

# Filter to colocalized
colocalized = [r for r in results if r.is_colocalized]
print(f"Found {len(colocalized)} colocalized loci")
```

## Input Format

GWAS files should contain:

| Column | Description |
|--------|-------------|
| chromosome | Chromosome (1-22, X, Y) |
| position | Base pair position |
| beta | Effect size |
| se | Standard error |
| variant_id | Unique identifier (optional) |

## Output

| Column | Description |
|--------|-------------|
| locus_id | Locus identifier |
| PP.H0-H4 | Posterior probabilities for each hypothesis |
| top_variant | Variant with highest colocalization evidence |
| colocalized | Boolean: PP.H4 > 0.8 |

## Integration with eQTL

```python
# GWAS vs eQTL colocalization
results = run_colocalization(
    gwas1_path=Path('disease_gwas.tsv'),
    gwas2_path=Path('gene_eqtl.tsv'),
    method='coloc',
    trait1_name='Disease',
    trait2_name='GENE_expression'
)

# Filter to colocalized genes
coloc_genes = results[results['colocalized']]
```

## References

- Giambartolomei C et al. (2014) Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics. PLOS Genetics.
- Hormozdiari F et al. (2016) Colocalization of GWAS and eQTL Signals Detects Target Genes. AJHG.
