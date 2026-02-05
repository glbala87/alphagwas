# AlphaGWAS

A GWAS variant prioritization pipeline using Google DeepMind's [AlphaGenome](https://github.com/google-deepmind/alphagenome) for functional effect prediction.

## Overview

AlphaGWAS takes GWAS summary statistics and prioritizes candidate causal variants by:

1. **Extracting significant variants** from GWAS summary statistics
2. **Finding LD proxies** using 1000 Genomes reference panel
3. **Lifting coordinates** from hg19 to hg38 (if needed)
4. **Predicting functional effects** using AlphaGenome across tissues and modalities
5. **Scoring and ranking** variants by predicted impact

## Installation

```bash
# Clone this repository
git clone https://github.com/YOUR_USERNAME/alphagwas.git
cd alphagwas

# Install dependencies
pip install -r requirements.txt

# Install AlphaGenome
git clone https://github.com/google-deepmind/alphagenome.git
pip install ./alphagenome
```

## Quick Start

```bash
# 1. Set up your API key
export ALPHAGENOME_API_KEY="your_key_here"

# 2. Copy your GWAS summary statistics
cp /path/to/your_gwas.tsv data/input/

# 3. Inspect your file to see column names
python setup_study.py --inspect data/input/your_gwas.tsv

# 4. Edit the config file with your column mappings
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

- Cardiovascular: Heart, Arteries (Aorta, Coronary, Tibial)
- Metabolic: Liver, Adipose, Pancreas
- Brain: Cortex, Cerebellum, Hippocampus
- Other: Whole Blood, Muscle, Lung, Kidney, etc.

## Requirements

- Python 3.8+
- AlphaGenome API access ([request here](https://deepmind.google.com/science/alphagenome))
- See `requirements.txt` for Python packages

## Project Structure

```
alphagwas/
├── config/              # Configuration files
├── data/
│   ├── input/          # GWAS files, reference data
│   ├── intermediate/   # Pipeline checkpoints
│   └── output/         # Final results
├── scripts/            # Pipeline step scripts
├── run_pipeline.py     # Main entry point
├── setup_study.py      # Setup helper
└── requirements.txt
```

## License

MIT License

## Acknowledgments

- [AlphaGenome](https://github.com/google-deepmind/alphagenome) by Google DeepMind
- Inspired by [alphagenome-test](https://github.com/aa9gj/alphagenome-test)
