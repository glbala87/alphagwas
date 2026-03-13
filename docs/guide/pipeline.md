# Pipeline Overview

AlphaGWAS processes GWAS summary statistics through a series of steps to prioritize variants by their predicted functional impact.

## Pipeline Steps

### Step 1: Extract Variants

Extracts genome-wide significant variants from GWAS summary statistics.

```bash
python run_pipeline.py --step 1
```

**Input:** GWAS summary statistics (TSV/CSV)
**Output:** `*_significant_variants.tsv`, `*_lead_snps.tsv`

### Step 2: Get LD Proxies

Identifies variants in linkage disequilibrium with lead SNPs.

```bash
python run_pipeline.py --step 2
```

**Input:** Lead SNPs from Step 1
**Output:** `*_variants_with_ld.tsv`

### Step 3: Liftover Coordinates

Converts genomic coordinates from hg19 to hg38 (required for AlphaGenome).

```bash
python run_pipeline.py --step 3
```

**Input:** Variants from Step 2
**Output:** `*_variants_hg38.tsv`

### Step 4: AlphaGenome Predictions

Runs AlphaGenome variant effect predictions across tissues and modalities.

```bash
python run_pipeline.py --step 4
```

**Input:** hg38 variants
**Output:** `*_alphagenome_predictions.parquet`

### Step 5: Score and Prioritize

Calculates consensus scores and ranks variants by functional impact.

```bash
python run_pipeline.py --step 5
```

**Input:** AlphaGenome predictions
**Output:** `*_ranked_variants.tsv`, `*_tissue_scores.tsv`

### Step 6: Visualizations (Optional)

Generates plots and HTML reports.

```bash
python run_pipeline.py --step 6
```

**Output:** PNG plots, HTML dashboard

## Data Flow

```
GWAS Summary Stats
        │
        ▼
┌───────────────┐
│ Step 1: Extract│
└───────────────┘
        │
        ▼
┌───────────────┐
│ Step 2: LD    │
└───────────────┘
        │
        ▼
┌───────────────┐
│ Step 3: Liftover│
└───────────────┘
        │
        ▼
┌───────────────┐
│ Step 4: Predict│
└───────────────┘
        │
        ▼
┌───────────────┐
│ Step 5: Score │
└───────────────┘
        │
        ▼
  Ranked Variants
```
