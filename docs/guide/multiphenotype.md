# Multi-Phenotype Comparison

Compare variant effects across multiple GWAS traits to identify pleiotropic variants and shared genetic architecture.

## Features

- **Pleiotropic variant identification**: Find variants affecting multiple traits
- **Cross-trait correlations**: Genetic correlation matrix
- **Effect comparison**: Compare effect sizes between traits
- **Miami plots**: Mirrored Manhattan plots
- **Clustering**: Hierarchically cluster related phenotypes

## Usage

### Python API

```python
from scripts.multiphenotype import (
    MultiPhenotypeAnalyzer,
    MultiPhenotypeVisualizer,
    PhenotypeData
)
import pandas as pd

# Load GWAS data
phenotypes = [
    PhenotypeData(name="LDL", gwas_df=pd.read_csv("ldl.tsv", sep="\t")),
    PhenotypeData(name="HDL", gwas_df=pd.read_csv("hdl.tsv", sep="\t")),
    PhenotypeData(name="TG", gwas_df=pd.read_csv("tg.tsv", sep="\t")),
    PhenotypeData(name="CAD", gwas_df=pd.read_csv("cad.tsv", sep="\t")),
]

# Initialize analyzer
analyzer = MultiPhenotypeAnalyzer(phenotypes)
```

### Identify Pleiotropic Variants

```python
# Find variants significant in 2+ traits
pleiotropic = analyzer.identify_pleiotropic_variants(
    pval_threshold=5e-8,
    min_phenotypes=2
)

print(f"Found {len(pleiotropic)} pleiotropic variants")
print(pleiotropic[['variant_id', 'n_phenotypes', 'phenotypes', 'effect_direction']])
```

Output:
```
variant_id          n_phenotypes  phenotypes     effect_direction
chr19:44908684      4             LDL,HDL,TG,CAD consistent
chr6:161010118      3             LDL,TG,CAD     opposing
chr1:55505647       2             LDL,CAD        consistent
```

### Genetic Correlations

```python
# Calculate cross-trait correlations
correlations = analyzer.calculate_genetic_correlation(method='zscore')
print(correlations)
```

Output:
```
        LDL    HDL     TG    CAD
LDL   1.00  -0.45   0.32   0.55
HDL  -0.45   1.00  -0.38  -0.42
TG    0.32  -0.38   1.00   0.28
CAD   0.55  -0.42   0.28   1.00
```

### Effect Size Comparison

```python
# Compare effects between two traits
comparison = analyzer.compare_effect_sizes("LDL", "CAD")

# Filter to variants with heterogeneous effects
heterogeneous = comparison[comparison['heterogeneity_p'] < 0.05]
```

### Visualizations

```python
visualizer = MultiPhenotypeVisualizer(analyzer)

# Correlation heatmap
fig = visualizer.plot_correlation_heatmap()
fig.savefig("correlation_heatmap.png", dpi=150)

# Effect comparison scatter
fig = visualizer.plot_effect_comparison("LDL", "CAD")
fig.savefig("ldl_vs_cad.png", dpi=150)

# Miami plot
fig = visualizer.plot_miami("LDL", "CAD")
fig.savefig("miami_plot.png", dpi=150)

# Pleiotropic summary
fig = visualizer.plot_pleiotropic_summary(pleiotropic)
fig.savefig("pleiotropic_summary.png", dpi=150)
```

### Full Analysis Pipeline

```python
from scripts.multiphenotype import run_multiphenotype_analysis
from pathlib import Path

results = run_multiphenotype_analysis(
    gwas_files=[
        Path("ldl_gwas.tsv"),
        Path("hdl_gwas.tsv"),
        Path("tg_gwas.tsv"),
        Path("cad_gwas.tsv")
    ],
    phenotype_names=["LDL", "HDL", "TG", "CAD"],
    output_dir=Path("multiphenotype_results"),
    pval_threshold=5e-8
)

# Access results
print(f"Pleiotropic variants: {len(results['pleiotropic'])}")
print(f"Genetic correlations:\n{results['correlations']}")
```

## Output Files

The full pipeline generates:

| File | Description |
|------|-------------|
| pleiotropic_variants.tsv | Variants significant in multiple traits |
| genetic_correlations.tsv | Pairwise genetic correlation matrix |
| variant_overlap.tsv | Set membership for UpSet plots |
| correlation_heatmap.png | Correlation visualization |
| pleiotropic_summary.png | Pleiotropy distribution |

## Interpreting Results

### Pleiotropic Variants

- **Consistent direction**: Variant affects traits in same direction (e.g., increases both LDL and CAD risk)
- **Opposing direction**: Variant has opposite effects (e.g., increases LDL but decreases HDL)

### Genetic Correlations

| Correlation | Interpretation |
|-------------|----------------|
| > 0.7 | Strong positive correlation |
| 0.3 - 0.7 | Moderate correlation |
| -0.3 - 0.3 | Weak/no correlation |
| < -0.3 | Negative correlation |

## Integration with AlphaGWAS Scoring

```python
# Score pleiotropic variants
from scripts import score_variants

scorer = score_variants.VariantScorer({})

# Get AlphaGenome predictions for pleiotropic variants
# ... (run predictions)

# Rank by functional impact
ranked = scorer.rank_variants(predictions)

# Merge with pleiotropy info
combined = pd.merge(
    ranked,
    pleiotropic[['variant_id', 'n_phenotypes', 'phenotypes']],
    on='variant_id'
)
```
