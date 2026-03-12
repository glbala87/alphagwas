# AlphaGWAS

**Prioritizing Variants at an Estimated Phenotype GWAS Locus with AlphaGenome**

A reproducible bioinformatics pipeline for identifying and ranking genetic variants at GWAS loci using Google DeepMind's [AlphaGenome](https://github.com/google-deepmind/alphagenome) AI model.

## Features

- **Parallel Processing** - Multi-threaded AlphaGenome predictions for faster execution
- **Progress Tracking** - Real-time progress bars with tqdm/rich
- **Automatic Retry** - Exponential backoff for robust API calls
- **Checkpoint/Resume** - Continue interrupted runs from the last checkpoint
- **Visualizations** - Manhattan plots, tissue heatmaps, effect distributions
- **Interactive Plots** - Zoomable Plotly visualizations with hover details
- **HTML Reports** - Publication-ready summary reports and dashboards
- **Variant Annotations** - ClinVar pathogenicity, gnomAD frequencies, GTEx eQTLs
- **Docker Support** - Containerized deployment for reproducibility
- **Unit Tests** - Comprehensive pytest test suite

## Overview

AlphaGWAS enables researchers to prioritize candidate causal variants from any GWAS by leveraging AlphaGenome's deep learning predictions of variant functional effects across multiple tissues and genomic modalities.

### Pipeline Workflow

1. **Variant Discovery** - Extract genome-wide significant SNPs from GWAS summary statistics
2. **LD Proxy Identification** - Find linked variants using 1000 Genomes reference panel
3. **Coordinate Liftover** - Convert genomic positions from hg19 to hg38
4. **Functional Predictions** - Run AlphaGenome's `predict_variant()` across tissues and modalities
5. **Variant Prioritization** - Score and rank variants by predicted functional impact
6. **Visualization** - Generate plots and HTML reports (optional)

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
  max_workers: 4        # Parallel workers
  max_retries: 3        # API retry attempts
  checkpoint_every: 50  # Save checkpoint every N variants
```

## Pipeline Usage

```bash
# Run complete pipeline
python run_pipeline.py --config config/config.yaml

# Run specific steps
python run_pipeline.py --step 1      # Extract variants
python run_pipeline.py --step 2      # Get LD proxies
python run_pipeline.py --step 3      # Liftover coordinates
python run_pipeline.py --step 4      # AlphaGenome predictions
python run_pipeline.py --step 5      # Score and rank
python run_pipeline.py --step 6      # Generate visualizations

# Run multiple steps
python run_pipeline.py --step 4,5,6

# Check setup
python run_pipeline.py --check

# Generate visualizations only (after pipeline completion)
python run_pipeline.py --visualize

# Run with verbose output
python run_pipeline.py --verbose

# Skip visualization step
python run_pipeline.py --skip-visualize
```

## AlphaGenome Prediction Options

```bash
# Run predictions with specific settings
python scripts/alphagenome_predict.py --config config.yaml

# Use 8 parallel workers
python scripts/alphagenome_predict.py --workers 8

# Disable parallel processing
python scripts/alphagenome_predict.py --no-parallel

# Dry run (validate without running)
python scripts/alphagenome_predict.py --dry-run

# Disable checkpoints
python scripts/alphagenome_predict.py --no-checkpoint
```

## Output

Results are saved to `data/output/`:

| File | Description |
|------|-------------|
| `*_ranked_variants.tsv` | All variants ranked by predicted functional impact |
| `*_tissue_scores.tsv` | Tissue-specific scores for each variant |
| `*_top20_variants.tsv` | Top 20 prioritized variants |
| `*_annotated_variants.tsv` | Variants with ClinVar/gnomAD/GTEx annotations |
| `*_manhattan.png` | Manhattan-style prioritization plot |
| `*_tissue_heatmap.png` | Heatmap of tissue-specific effects |
| `*_effect_distribution.png` | Effect size distributions |
| `*_report.html` | Static HTML summary report |
| `*_dashboard.html` | Interactive Plotly dashboard |
| `*_manhattan_interactive.html` | Interactive Manhattan plot |

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

### Core Dependencies
- pandas, numpy, scipy
- pyyaml
- pyliftover, cyvcf2
- pyarrow

### Recommended (for better UX)
- tqdm - Progress bars
- rich - Rich console output

### Visualization
- matplotlib
- seaborn

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
│   ├── extract_variants.py
│   ├── get_ld_proxies.py
│   ├── liftover.py
│   ├── alphagenome_predict.py
│   ├── score_variants.py
│   ├── visualize.py    # Visualization module
│   └── utils.py        # Utility functions
├── run_pipeline.py     # Main entry point
├── setup_study.py      # Setup helper
└── requirements.txt
```

## Docker Usage

```bash
# Build Docker image
docker build -t alphagwas:latest .

# Run pipeline with Docker
docker run -v $(pwd)/data:/app/data -v $(pwd)/config:/app/config \
  -e ALPHAGENOME_API_KEY=$ALPHAGENOME_API_KEY \
  alphagwas:latest --config /app/config/config.yaml

# Using docker-compose
docker-compose run --rm alphagwas --config /app/config/config.yaml

# Run tests in Docker
docker-compose run --rm alphagwas-test

# Interactive development shell
docker-compose run --rm alphagwas-dev /bin/bash
```

## Variant Annotations

Add external annotations to prioritized variants:

```bash
# Run annotation step
python scripts/annotate.py --config config/config.yaml

# Disable specific annotation sources
python scripts/annotate.py --no-clinvar --no-gnomad
```

### Annotation Sources

| Source | Data | API |
|--------|------|-----|
| **ClinVar** | Pathogenicity classifications | NCBI E-utilities |
| **gnomAD** | Population allele frequencies | gnomAD GraphQL API |
| **GTEx** | eQTL associations | GTEx Portal API |

### Annotation Scoring

Variants receive an `annotation_score` (0-1) based on:
- **ClinVar** (40%): Pathogenic/likely pathogenic variants score higher
- **gnomAD** (30%): Rare variants (low AF) score higher
- **GTEx** (30%): Variants with more eQTL associations score higher

## Interactive Visualizations

Generate interactive Plotly visualizations:

```bash
# Generate interactive dashboard
python scripts/visualize_interactive.py --config config/config.yaml
```

Features:
- Hover to see variant details (rsID, position, score, genes)
- Click and drag to zoom
- Double-click to reset view
- Export as PNG/SVG

## Testing

```bash
# Run all tests
make test

# Run tests with coverage
make test-cov

# Run specific test file
python -m pytest tests/test_alphagenome_predict.py -v

# Run fast tests only
make test-fast
```

## Troubleshooting

### AlphaGenome API Issues
- Ensure `ALPHAGENOME_API_KEY` is set
- Check your API quota and rate limits
- Use `--max-retries 5` for unreliable connections

### Memory Issues
- Reduce `--workers` for parallel processing
- Process chromosomes separately with locus configuration
- Use checkpoint/resume for large variant sets

### Missing Dependencies
```bash
# Install all optional dependencies
pip install tqdm rich matplotlib seaborn
```

## Author

**BalaSubramani Gattu Linga** ([@glbala87](https://github.com/glbala87))

## License

MIT License

## Acknowledgments

- [AlphaGenome](https://github.com/google-deepmind/alphagenome) by Google DeepMind
- Inspired by [alphagenome-test](https://github.com/aa9gj/alphagenome-test)
