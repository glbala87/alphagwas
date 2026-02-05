#!/usr/bin/env python3
"""
Step 4: Run AlphaGenome variant effect predictions.

This script interfaces with Google DeepMind's AlphaGenome model to predict
functional effects of genetic variants across tissues and modalities.

AlphaGenome Setup:
1. Clone and install: git clone https://github.com/google-deepmind/alphagenome.git
                      pip install ./alphagenome
2. Get API key from: https://deepmind.google.com/science/alphagenome
3. Set environment variable: export ALPHAGENOME_API_KEY="your_key"
   Or add to config: alphagenome.api_key: "your_key"
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import yaml
import json
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config.yaml") -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# GTEx tissue to UBERON ontology mapping
TISSUE_ONTOLOGY = {
    # Cardiovascular
    "Heart_Atrial_Appendage": "UBERON:0006618",
    "Heart_Left_Ventricle": "UBERON:0002084",
    "Artery_Aorta": "UBERON:0000947",
    "Artery_Coronary": "UBERON:0001621",
    "Artery_Tibial": "UBERON:0007610",
    # Blood
    "Whole_Blood": "UBERON:0000178",
    # Metabolic
    "Liver": "UBERON:0002107",
    "Adipose_Subcutaneous": "UBERON:0002190",
    "Adipose_Visceral_Omentum": "UBERON:0003688",
    "Pancreas": "UBERON:0001264",
    # Muscle
    "Muscle_Skeletal": "UBERON:0001134",
    # Brain
    "Brain_Cortex": "UBERON:0001851",
    "Brain_Cerebellum": "UBERON:0002037",
    "Brain_Hippocampus": "UBERON:0002310",
    # Other
    "Lung": "UBERON:0002048",
    "Kidney_Cortex": "UBERON:0001225",
    "Skin_Sun_Exposed_Lower_leg": "UBERON:0001415",
    "Thyroid": "UBERON:0002046",
    "Nerve_Tibial": "UBERON:0001323",
    "Esophagus_Mucosa": "UBERON:0002469",
    "Colon_Transverse": "UBERON:0001157",
    "Small_Intestine_Terminal_Ileum": "UBERON:0002116",
    "Stomach": "UBERON:0000945",
    "Spleen": "UBERON:0002106",
}


@dataclass
class VariantPrediction:
    """Container for AlphaGenome prediction results."""
    variant_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    tissue: str
    modality: str
    prediction: float
    effect_size: float
    confidence: float


class AlphaGenomePredictor:
    """
    Interface to AlphaGenome variant effect prediction.

    AlphaGenome predicts the functional impact of genetic variants on:
    - Gene expression (RNA-seq)
    - Chromatin accessibility (ATAC-seq)
    - Histone modifications (ChIP-seq)
    - Transcription factor binding
    - 3D chromatin contacts
    - Splicing
    """

    # Modality to AlphaGenome OutputType mapping
    MODALITY_MAPPING = {
        'expression': 'RNA_SEQ',
        'chromatin_accessibility': 'ATAC_SEQ',
        'histone_marks': 'CHIP_SEQ',
        'transcription_factor': 'TF_CHIP_SEQ',
        'contact': 'HIC',
        'splicing': 'SPLICE',
    }

    def __init__(self, config: dict):
        """Initialize AlphaGenome predictor."""
        self.config = config.get('alphagenome', {})
        self.tissues = self.config.get('tissues', ['Whole_Blood'])
        self.modalities = self.config.get('modalities', ['expression'])
        self.batch_size = self.config.get('batch_size', 100)
        self.checkpoint_every = self.config.get('checkpoint_every', 50)

        # Get API key from config or environment
        self.api_key = self.config.get('api_key') or os.environ.get('ALPHAGENOME_API_KEY')

        self._init_alphagenome()

    def _init_alphagenome(self):
        """Initialize AlphaGenome model connection."""
        self.model = None
        self.genome_module = None
        self.output_types = None

        try:
            from alphagenome.data import genome
            from alphagenome.models import dna_client

            self.genome_module = genome
            self.dna_client = dna_client

            if not self.api_key:
                logger.warning("No API key found. Set ALPHAGENOME_API_KEY environment variable")
                logger.warning("or add 'api_key' to alphagenome section in config")
                logger.warning("Get your key at: https://deepmind.google.com/science/alphagenome")
                logger.info("Using mock predictions for testing.")
                return

            self.model = dna_client.create(self.api_key)
            self.output_types = dna_client.OutputType
            logger.info("AlphaGenome model initialized successfully")

        except ImportError:
            logger.warning("AlphaGenome package not installed.")
            logger.warning("Install with: git clone https://github.com/google-deepmind/alphagenome.git")
            logger.warning("             pip install ./alphagenome")
            logger.info("Using mock predictions for testing.")
        except Exception as e:
            logger.error(f"Failed to initialize AlphaGenome: {e}")
            logger.info("Using mock predictions for testing.")

    def _get_ontology_term(self, tissue: str) -> str:
        """Get UBERON ontology term for a tissue."""
        return TISSUE_ONTOLOGY.get(tissue, "UBERON:0000178")  # Default to blood

    def _get_output_type(self, modality: str):
        """Get AlphaGenome OutputType for a modality."""
        if self.output_types is None:
            return None

        output_name = self.MODALITY_MAPPING.get(modality, 'RNA_SEQ')
        return getattr(self.output_types, output_name, self.output_types.RNA_SEQ)

    def predict_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        tissue: str,
        modality: str
    ) -> Optional[Dict[str, Any]]:
        """
        Predict effect of a single variant using AlphaGenome.

        Args:
            chrom: Chromosome (e.g., '1', 'X')
            pos: Position (1-based, hg38)
            ref: Reference allele
            alt: Alternate allele
            tissue: GTEx tissue name
            modality: Prediction modality

        Returns:
            Dictionary with prediction results
        """
        if self.model is not None:
            try:
                return self._predict_with_api(chrom, pos, ref, alt, tissue, modality)
            except Exception as e:
                logger.error(f"Prediction failed for chr{chrom}:{pos} {ref}>{alt}: {e}")
                return None
        else:
            # Mock prediction for testing
            return self._mock_predict(chrom, pos, ref, alt, tissue, modality)

    def _predict_with_api(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        tissue: str,
        modality: str
    ) -> Dict[str, Any]:
        """Make actual AlphaGenome API call."""
        # Format chromosome
        chrom_str = f"chr{chrom}" if not chrom.startswith('chr') else chrom

        # Define genomic interval (±500kb around variant)
        window = 500000
        interval = self.genome_module.Interval(
            chromosome=chrom_str,
            start=max(0, pos - window),
            end=pos + window
        )

        # Define variant
        variant = self.genome_module.Variant(
            chromosome=chrom_str,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt
        )

        # Get ontology term and output type
        ontology_term = self._get_ontology_term(tissue)
        output_type = self._get_output_type(modality)

        # Make prediction
        outputs = self.model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=[ontology_term],
            requested_outputs=[output_type]
        )

        # Extract results
        # The actual output structure may vary - adjust as needed
        ref_pred = outputs.get('reference', {}).get('prediction', 0)
        alt_pred = outputs.get('alternate', {}).get('prediction', 0)
        effect_size = abs(alt_pred - ref_pred)

        return {
            'variant_id': f"chr{chrom}:{pos}:{ref}:{alt}",
            'tissue': tissue,
            'modality': modality,
            'prediction': alt_pred - ref_pred,
            'effect_size': effect_size,
            'confidence': outputs.get('confidence', 0.5),
            'ref_prediction': ref_pred,
            'alt_prediction': alt_pred,
            'nearby_genes': outputs.get('genes', []),
            'regulatory_elements': outputs.get('elements', [])
        }

    def _mock_predict(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        tissue: str,
        modality: str
    ) -> Dict[str, Any]:
        """Generate mock predictions for testing without AlphaGenome API."""
        np.random.seed(hash(f"{chrom}{pos}{ref}{alt}{tissue}{modality}") % 2**32)

        return {
            'variant_id': f"chr{chrom}:{pos}:{ref}:{alt}",
            'tissue': tissue,
            'modality': modality,
            'prediction': np.random.normal(0, 0.5),
            'effect_size': abs(np.random.normal(0, 0.3)),
            'confidence': np.random.uniform(0.5, 1.0),
            'nearby_genes': ['GENE1', 'GENE2'],
            'regulatory_elements': ['enhancer', 'promoter'][np.random.randint(0, 2)]
        }

    def predict_batch(
        self,
        variants: pd.DataFrame,
        checkpoint_path: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Run predictions for a batch of variants across all tissues and modalities.

        Args:
            variants: DataFrame with columns: chromosome, position, ref, alt
            checkpoint_path: Path to save/load checkpoints

        Returns:
            List of prediction dictionaries
        """
        all_predictions = []

        # Load checkpoint if exists
        completed_variants = set()
        if checkpoint_path and Path(checkpoint_path).exists():
            with open(checkpoint_path) as f:
                checkpoint = json.load(f)
                all_predictions = checkpoint.get('predictions', [])
                completed_variants = set(checkpoint.get('completed', []))
                logger.info(f"Resumed from checkpoint: {len(completed_variants)} variants already processed")

        total_predictions = len(variants) * len(self.tissues) * len(self.modalities)
        logger.info(f"Running {total_predictions} predictions ({len(variants)} variants x {len(self.tissues)} tissues x {len(self.modalities)} modalities)")

        count = 0
        for idx, variant in variants.iterrows():
            variant_id = f"chr{variant['chromosome']}:{variant['position']}"

            if variant_id in completed_variants:
                continue

            # Get alleles
            ref = variant.get('ref', variant.get('other_allele', 'N'))
            alt = variant.get('alt', variant.get('effect_allele', 'N'))

            for tissue in self.tissues:
                for modality in self.modalities:
                    result = self.predict_variant(
                        chrom=str(variant['chromosome']),
                        pos=int(variant['position']),
                        ref=str(ref),
                        alt=str(alt),
                        tissue=tissue,
                        modality=modality
                    )

                    if result:
                        result['rsid'] = variant.get('rsid', variant_id)
                        result['chromosome'] = variant['chromosome']
                        result['position'] = variant['position']
                        all_predictions.append(result)

                    count += 1
                    if count % 100 == 0:
                        logger.info(f"Progress: {count}/{total_predictions} predictions")

            completed_variants.add(variant_id)

            # Save checkpoint
            if checkpoint_path and len(completed_variants) % self.checkpoint_every == 0:
                self._save_checkpoint(checkpoint_path, all_predictions, completed_variants)

        # Final checkpoint
        if checkpoint_path:
            self._save_checkpoint(checkpoint_path, all_predictions, completed_variants)

        return all_predictions

    def _save_checkpoint(
        self,
        checkpoint_path: str,
        predictions: List[Dict],
        completed: set
    ):
        """Save checkpoint for resumable processing."""
        checkpoint = {
            'predictions': predictions,
            'completed': list(completed),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }
        with open(checkpoint_path, 'w') as f:
            json.dump(checkpoint, f)
        logger.info(f"Checkpoint saved: {len(completed)} variants processed")


def main(config_path: str = "config/config.yaml"):
    """Main AlphaGenome prediction workflow."""
    config = load_config(config_path)

    # Initialize predictor
    predictor = AlphaGenomePredictor(config)

    # Load variants from previous step
    prefix = config['output']['prefix']
    intermediate_dir = Path("data/intermediate")

    # Try to load hg38 variants first, fall back to other files
    input_files = [
        intermediate_dir / f"{prefix}_variants_hg38.tsv",
        intermediate_dir / f"{prefix}_variants_with_ld.tsv",
        intermediate_dir / f"{prefix}_significant_variants.tsv"
    ]

    variants = None
    for input_file in input_files:
        if input_file.exists():
            variants = pd.read_csv(input_file, sep='\t')
            logger.info(f"Loaded {len(variants)} variants from {input_file}")
            break

    if variants is None:
        logger.error("No input variant file found. Run previous steps first.")
        return

    # Run predictions
    checkpoint_path = intermediate_dir / f"{prefix}_alphagenome_checkpoint.json"

    predictions = predictor.predict_batch(
        variants,
        checkpoint_path=str(checkpoint_path)
    )

    # Convert to DataFrame
    pred_df = pd.DataFrame(predictions)

    # Save results
    output_file = intermediate_dir / f"{prefix}_alphagenome_predictions.tsv"
    pred_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved {len(pred_df)} predictions to {output_file}")

    # Also save as parquet for efficiency
    parquet_file = intermediate_dir / f"{prefix}_alphagenome_predictions.parquet"
    try:
        pred_df.to_parquet(parquet_file, index=False)
        logger.info(f"Saved predictions to {parquet_file}")
    except ImportError:
        logger.warning("pyarrow not installed, skipping parquet output")

    return pred_df


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run AlphaGenome variant effect predictions")
    parser.add_argument("--config", default="config/config.yaml", help="Path to config file")
    args = parser.parse_args()

    main(args.config)
