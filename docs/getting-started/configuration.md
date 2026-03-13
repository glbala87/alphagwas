# Configuration

AlphaGWAS uses YAML configuration files to customize pipeline behavior.

## Configuration File Structure

```yaml
# config/config.yaml

# Study Information
study:
  name: "my_study"
  phenotype: "Type 2 Diabetes"
  description: "T2D GWAS meta-analysis"

# GWAS Input Settings
gwas:
  input_file: "data/input/t2d_gwas.tsv"
  columns:
    chromosome: "CHR"
    position: "BP"
    rsid: "SNP"
    effect_allele: "A1"
    other_allele: "A2"
    pvalue: "P"
    beta: "BETA"           # Optional
    se: "SE"               # Optional
  pvalue_threshold: 5.0e-8
  n_samples: 100000        # For fine-mapping

# Reference Genome Settings
reference:
  input_build: "hg19"
  output_build: "hg38"
  liftover_chain: "data/input/hg19ToHg38.over.chain.gz"

# LD Proxy Settings
ld:
  vcf_path: "/path/to/1000genomes/"
  population: "EUR"
  r2_threshold: 0.8
  window_kb: 500

# Locus Definitions (optional)
loci:
  - name: "TCF7L2"
    chromosome: "10"
    start: 114710000
    end: 114930000
    lead_snp: "rs7903146"
  - name: "PPARG"
    chromosome: "3"
    start: 12300000
    end: 12500000

# AlphaGenome Settings
alphagenome:
  api_key: "${ALPHAGENOME_API_KEY}"  # From environment
  tissues:
    - "Pancreas"
    - "Liver"
    - "Adipose_Subcutaneous"
    - "Muscle_Skeletal"
  modalities:
    - "expression"
    - "chromatin_accessibility"
  max_workers: 4
  max_retries: 3
  batch_size: 100
  checkpoint_every: 50

# Scoring Settings
scoring:
  consensus_method: "mean"  # mean, median, or max
  min_effect_threshold: 0.1

# Fine-mapping Settings
finemapping:
  method: "susie"           # susie or finemap
  pip_weight: 0.5           # Weight for PIP vs AlphaGenome score
  susie:
    L: 10                   # Max causal variants
    coverage: 0.95

# Enrichment Settings
enrichment:
  use_enrichr: true
  top_n: 100
  p_threshold: 0.05

# Annotation Settings
annotations:
  max_workers: 2

# Output Settings
output:
  dir: "data/output"
  prefix: "t2d_analysis"
```

## Environment Variables

Configuration values can reference environment variables:

```yaml
alphagenome:
  api_key: "${ALPHAGENOME_API_KEY}"
```

## Tissue Options

Available GTEx tissues:

| Category | Tissues |
|----------|---------|
| **Cardiovascular** | Heart_Atrial_Appendage, Heart_Left_Ventricle, Artery_Aorta, Artery_Coronary |
| **Metabolic** | Liver, Pancreas, Adipose_Subcutaneous, Adipose_Visceral_Omentum |
| **Brain** | Brain_Cortex, Brain_Cerebellum, Brain_Hippocampus |
| **Other** | Whole_Blood, Muscle_Skeletal, Lung, Kidney_Cortex |

## Modality Options

| Modality | Description |
|----------|-------------|
| `expression` | Gene expression (RNA-seq) |
| `chromatin_accessibility` | Chromatin openness (ATAC-seq) |
| `histone_marks` | Histone modifications (ChIP-seq) |
| `transcription_factor` | TF binding |
| `contact` | 3D chromatin contacts (Hi-C) |
| `splicing` | Splice site effects |

## Validating Configuration

```bash
python run_pipeline.py --check --config config/my_study.yaml
```
