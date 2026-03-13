# Installation

## Requirements

- Python 3.9 or higher
- AlphaGenome API access ([request here](https://deepmind.google.com/science/alphagenome))

## Installation Methods

### Using pip (Recommended)

```bash
pip install alphagwas
```

### With optional dependencies

```bash
# Full installation with all features
pip install alphagwas[all]

# Specific feature sets
pip install alphagwas[viz]        # Visualization (matplotlib, plotly)
pip install alphagwas[genomics]   # Genomics tools (pyliftover, cyvcf2)
pip install alphagwas[streamlit]  # Web dashboard
pip install alphagwas[dev]        # Development tools
```

### From source

```bash
git clone https://github.com/glbala87/alphagwas.git
cd alphagwas
pip install -e .
```

### Using Docker

```bash
# Pull from GitHub Container Registry
docker pull ghcr.io/glbala87/alphagwas:latest

# Or build locally
docker build -t alphagwas .
```

## AlphaGenome Setup

AlphaGWAS requires the AlphaGenome package from Google DeepMind:

```bash
# Clone AlphaGenome
git clone https://github.com/google-deepmind/alphagenome.git

# Install
pip install ./alphagenome
```

### API Key

1. Request access at [AlphaGenome Portal](https://deepmind.google.com/science/alphagenome)
2. Set your API key:

```bash
export ALPHAGENOME_API_KEY="your_key_here"
```

Or add to your config file:

```yaml
alphagenome:
  api_key: "your_key_here"
```

## Verify Installation

```bash
# Check dependencies
python run_pipeline.py --check

# Run tests
python -m pytest tests/ -v
```

## Optional Dependencies

### For Fine-mapping (SuSiE)

Install R and susieR:

```r
install.packages("susieR")
```

### For Fine-mapping (FINEMAP)

Download from [FINEMAP website](http://www.christianbenner.com/) and add to PATH.

## Troubleshooting

### Common Issues

**ImportError: No module named 'cyvcf2'**

```bash
# On macOS
brew install htslib
pip install cyvcf2
```

**pyliftover chain file errors**

Chain files are downloaded automatically. If issues persist:

```bash
# Manual download
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
mv hg19ToHg38.over.chain data/input/
```

**Memory issues with large GWAS**

Reduce parallel workers:

```bash
python run_pipeline.py --workers 2
```
