# AlphaGWAS

**Prioritizing Variants at an Estimated Phenotype GWAS Locus with AlphaGenome**

A reproducible bioinformatics pipeline for identifying and ranking genetic variants at GWAS loci using Google DeepMind's [AlphaGenome](https://github.com/google-deepmind/alphagenome) AI model.

## Overview

AlphaGWAS enables researchers to prioritize candidate causal variants from any GWAS by leveraging AlphaGenome's deep learning predictions of variant functional effects across multiple tissues and genomic modalities.

### Pipeline Workflow

1. **Variant Discovery** - Extract genome-wide significant SNPs from GWAS summary statistics
2. **LD Proxy Identification** - Find linked variants using 1000 Genomes reference panel
3. **Coordinate Liftover** - Convert genomic positions from hg19 to hg38
4. **Functional Predictions** - Run AlphaGenome's `predict_variant()` across tissues and modalities
5. **Variant Prioritization** - Score and rank variants by predicted functional impact

## Installation

```bash
# Clone this repository
git clone https://github.com/glbala87/alphagwas.git
cd alphagwas

# Install dependencies
pip install -r requirements.txt

# Install AlphaGenome
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

## Quick Start

```bash
# 1. Set up your API key (get from https://deepmind.google.com/science/alphagenome)
export ALPHAGENOME_API_KEY="your_key_here"

# 2. Copy your GWAS summary statistics
cp /path/to/your_gwas.tsv data/input/

# 3. Inspect your file to see column names
python setup_study.py --inspect data/input/your_gwas.tsv

# 4. Edit the config file with your column mappings and phenotype
nano config/config.yaml

# 5. Run the pipeline
python run_pipeline.py --config config/config.yaml
```

## Configuration

Edit `config/config.yaml` to match your study:

```yaml
study:
  name: "my_gwas"
  phenotype: "my_trait"

gwas:
  input_file: "data/input/my_gwas.tsv"
  columns:
    chromosome: "CHR"
    position: "BP"
    rsid: "SNP"
    effect_allele: "A1"
    other_allele: "A2"
    pvalue: "P"
  pvalue_threshold: 5.0e-8

# Define loci of interest
loci:
  - name: "my_locus"
    chromosome: "1"
    start: 39000000
    end: 40000000
    lead_snp: "rs12345"

alphagenome:
  tissues:
    - "Heart_Left_Ventricle"
    - "Artery_Coronary"
    - "Whole_Blood"
  modalities:
    - "expression"
    - "chromatin_accessibility"
```

## Pipeline Steps

Run individual steps or the complete pipeline:

```bash
# Run complete pipeline
python run_pipeline.py --config config/config.yaml

# Run specific steps
python run_pipeline.py --step 1      # Extract variants
python run_pipeline.py --step 2      # Get LD proxies
python run_pipeline.py --step 3      # Liftover coordinates
python run_pipeline.py --step 4      # AlphaGenome predictions
python run_pipeline.py --step 5      # Score and rank

# Run multiple steps
python run_pipeline.py --step 4,5

# Check setup
python run_pipeline.py --check
```

## Output

Results are saved to `data/output/`:

| File | Description |
|------|-------------|
| `*_ranked_variants.tsv` | All variants ranked by predicted functional impact |
| `*_tissue_scores.tsv` | Tissue-specific scores for each variant |
| `*_top20_variants.tsv` | Top 20 prioritized variants |

## Supported Tissues

The pipeline maps GTEx tissue names to UBERON ontology terms for AlphaGenome:

| Category | Tissues |
|----------|---------|
| Cardiovascular | Heart (Atrial Appendage, Left Ventricle), Arteries (Aorta, Coronary, Tibial) |
| Metabolic | Liver, Adipose (Subcutaneous, Visceral), Pancreas |
| Brain | Cortex, Cerebellum, Hippocampus |
| Other | Whole Blood, Skeletal Muscle, Lung, Kidney, Thyroid |

## Prediction Modalities

AlphaGenome predicts variant effects on:

- **Gene Expression** (RNA-seq)
- **Chromatin Accessibility** (ATAC-seq)
- **Histone Modifications** (ChIP-seq)
- **Transcription Factor Binding**
- **3D Chromatin Contacts** (Hi-C)
- **Splicing**

## Requirements

- Python 3.8+
- AlphaGenome API access ([request here](https://deepmind.google.com/science/alphagenome))
- See `requirements.txt` for Python packages

## Project Structure

```
alphagwas/
├── config/              # Configuration files
│   ├── config.yaml
│   └── cardiovascular_config.yaml
├── data/
│   ├── input/          # GWAS files, reference data
│   ├── intermediate/   # Pipeline checkpoints
│   └── output/         # Final results
├── scripts/            # Pipeline step scripts
├── run_pipeline.py     # Main entry point
├── setup_study.py      # Setup helper
└── requirements.txt
```

## Author

**BalaSubramani Gattu Linga** ([@glbala87](https://github.com/glbala87))

## License

MIT License

## Acknowledgments

- [AlphaGenome](https://github.com/google-deepmind/alphagenome) by Google DeepMind
- Inspired by [alphagenome-test](https://github.com/aa9gj/alphagenome-test)
