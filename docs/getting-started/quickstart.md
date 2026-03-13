# Quick Start

This guide walks you through running AlphaGWAS on your GWAS data in 5 minutes.

## Prerequisites

- AlphaGWAS installed ([Installation Guide](installation.md))
- GWAS summary statistics file
- AlphaGenome API key

## Step 1: Prepare Your Data

Place your GWAS summary statistics in `data/input/`:

```bash
cp /path/to/your_gwas.tsv data/input/
```

### Supported Formats

- Tab-separated values (`.tsv`, `.txt`)
- Comma-separated values (`.csv`)
- Gzipped files (`.gz`)

### Required Columns

Your file should contain:

| Column | Description | Example |
|--------|-------------|---------|
| Chromosome | Chromosome number | 1, 2, ..., 22, X |
| Position | Genomic position | 12345678 |
| rsID | Variant identifier | rs12345 |
| Effect allele | Effect/risk allele | A |
| Other allele | Reference allele | G |
| P-value | Association p-value | 1.5e-10 |

## Step 2: Create Configuration

Copy and edit the example config:

```bash
cp config/example_config.yaml config/my_study.yaml
```

Edit `config/my_study.yaml`:

```yaml
study:
  name: "my_gwas_study"
  phenotype: "Blood Pressure"

gwas:
  input_file: "data/input/your_gwas.tsv"
  columns:
    chromosome: "CHR"      # Your column name
    position: "BP"         # Your column name
    rsid: "SNP"           # Your column name
    effect_allele: "A1"   # Your column name
    other_allele: "A2"    # Your column name
    pvalue: "P"           # Your column name
  pvalue_threshold: 5.0e-8

alphagenome:
  tissues:
    - "Heart_Left_Ventricle"
    - "Whole_Blood"
  modalities:
    - "expression"
    - "chromatin_accessibility"

output:
  dir: "data/output"
  prefix: "my_study"
```

## Step 3: Set API Key

```bash
export ALPHAGENOME_API_KEY="your_key_here"
```

## Step 4: Run the Pipeline

```bash
python run_pipeline.py --config config/my_study.yaml
```

### Pipeline Output

```
============================================================
# AlphaGWAS Pipeline
# Study: my_gwas_study
# Phenotype: Blood Pressure
============================================================

Step 1: Extract GWAS variants
  Found 150 genome-wide significant variants

Step 2: Get LD proxies
  Added 450 LD proxy variants (r² ≥ 0.8)

Step 3: Liftover coordinates
  Successfully lifted 598/600 variants to hg38

Step 4: AlphaGenome predictions
  Running 2392 predictions (598 variants × 2 tissues × 2 modalities)
  Progress: 100%

Step 5: Score and prioritize
  Ranked 598 variants by functional impact

============================================================
PIPELINE COMPLETE
Results saved to: data/output/my_study_*
============================================================
```

## Step 5: Explore Results

### View Top Variants

```bash
head data/output/my_study_top20_variants.tsv
```

### Generate Visualizations

```bash
python run_pipeline.py --visualize
```

### Launch Dashboard

```bash
streamlit run app.py
```

## Next Steps

- [Configuration Guide](configuration.md) - Customize pipeline settings
- [Visualizations](../guide/visualizations.md) - Create publication-ready plots
- [Fine-mapping](../guide/finemapping.md) - Statistical fine-mapping integration
- [API Reference](../api/alphagenome_predict.md) - Detailed API documentation
