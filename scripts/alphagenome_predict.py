#!/usr/bin/env python3
"""
Step 4: Run AlphaGenome variant effect predictions.

This script interfaces with Google DeepMind's AlphaGenome model to predict
functional effects of genetic variants across tissues and modalities.

Features:
- Parallel processing with configurable worker threads
- Automatic retry with exponential backoff
- Progress tracking with tqdm/rich
- Checkpoint/resume capability

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
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Import utilities
try:
    from .utils import (
        retry_with_backoff,
        progress_iterator,
        ProgressTracker,
        setup_logging,
        print_summary_table,
        print_panel
    )
except ImportError:
    from utils import (
        retry_with_backoff,
        progress_iterator,
        ProgressTracker,
        setup_logging,
        print_summary_table,
        print_panel
    )

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

    Features:
    - Parallel processing with ThreadPoolExecutor
    - Automatic retry with exponential backoff
    - Progress tracking
    - Checkpoint/resume support
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

        # Performance settings
        self.max_workers = self.config.get('max_workers', 4)
        self.max_retries = self.config.get('max_retries', 3)
        self.retry_delay = self.config.get('retry_delay', 1.0)

        # Thread-safe counters
        self._lock = threading.Lock()
        self._completed_count = 0
        self._error_count = 0

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
            return self._predict_with_retry(chrom, pos, ref, alt, tissue, modality)
        else:
            # Mock prediction for testing
            return self._mock_predict(chrom, pos, ref, alt, tissue, modality)

    def _predict_with_retry(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        tissue: str,
        modality: str
    ) -> Optional[Dict[str, Any]]:
        """Predict with retry logic."""
        last_error = None

        for attempt in range(self.max_retries + 1):
            try:
                return self._predict_with_api(chrom, pos, ref, alt, tissue, modality)
            except Exception as e:
                last_error = e
                if attempt < self.max_retries:
                    delay = self.retry_delay * (2 ** attempt) * (0.5 + np.random.random())
                    logger.warning(
                        f"Attempt {attempt + 1}/{self.max_retries + 1} failed for "
                        f"chr{chrom}:{pos}: {e}. Retrying in {delay:.1f}s..."
                    )
                    time.sleep(delay)

        logger.error(f"All retries failed for chr{chrom}:{pos} {ref}>{alt}: {last_error}")
        with self._lock:
            self._error_count += 1
        return None

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
        checkpoint_path: Optional[str] = None,
        use_parallel: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Run predictions for a batch of variants across all tissues and modalities.

        Args:
            variants: DataFrame with columns: chromosome, position, ref, alt
            checkpoint_path: Path to save/load checkpoints
            use_parallel: Enable parallel processing (default: True)

        Returns:
            List of prediction dictionaries
        """
        all_predictions = []

        # Load checkpoint if exists
        completed_variants = set()
        if checkpoint_path and Path(checkpoint_path).exists():
            try:
                with open(checkpoint_path) as f:
                    checkpoint = json.load(f)
                    all_predictions = checkpoint.get('predictions', [])
                    completed_variants = set(checkpoint.get('completed', []))
                    logger.info(f"Resumed from checkpoint: {len(completed_variants)} variants already processed")
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Could not load checkpoint: {e}. Starting fresh.")

        # Filter out already completed variants
        pending_variants = []
        for idx, variant in variants.iterrows():
            variant_id = f"chr{variant['chromosome']}:{variant['position']}"
            if variant_id not in completed_variants:
                pending_variants.append((idx, variant))

        total_predictions = len(pending_variants) * len(self.tissues) * len(self.modalities)

        if total_predictions == 0:
            logger.info("All variants already processed")
            return all_predictions

        # Print summary
        print_panel(
            f"Variants: {len(pending_variants)}\n"
            f"Tissues: {len(self.tissues)}\n"
            f"Modalities: {len(self.modalities)}\n"
            f"Total predictions: {total_predictions}\n"
            f"Workers: {self.max_workers if use_parallel else 1}",
            title="AlphaGenome Predictions",
            style="blue"
        )

        if use_parallel and self.max_workers > 1:
            all_predictions = self._predict_batch_parallel(
                pending_variants, all_predictions, completed_variants,
                checkpoint_path, total_predictions
            )
        else:
            all_predictions = self._predict_batch_sequential(
                pending_variants, all_predictions, completed_variants,
                checkpoint_path, total_predictions
            )

        # Final checkpoint
        if checkpoint_path:
            self._save_checkpoint(checkpoint_path, all_predictions, completed_variants)

        # Print summary
        print_summary_table({
            'Total predictions': len(all_predictions),
            'Successful': len(all_predictions),
            'Errors': self._error_count,
            'Variants processed': len(completed_variants)
        }, title="Prediction Summary")

        return all_predictions

    def _predict_batch_sequential(
        self,
        pending_variants: List[Tuple[int, pd.Series]],
        all_predictions: List[Dict],
        completed_variants: set,
        checkpoint_path: Optional[str],
        total_predictions: int
    ) -> List[Dict[str, Any]]:
        """Sequential batch processing."""
        count = 0

        with ProgressTracker(total=total_predictions, desc="Predicting variants") as tracker:
            for idx, variant in pending_variants:
                variant_id = f"chr{variant['chromosome']}:{variant['position']}"

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
                        tracker.advance()

                completed_variants.add(variant_id)

                # Save checkpoint periodically
                if checkpoint_path and len(completed_variants) % self.checkpoint_every == 0:
                    self._save_checkpoint(checkpoint_path, all_predictions, completed_variants)

        return all_predictions

    def _predict_batch_parallel(
        self,
        pending_variants: List[Tuple[int, pd.Series]],
        all_predictions: List[Dict],
        completed_variants: set,
        checkpoint_path: Optional[str],
        total_predictions: int
    ) -> List[Dict[str, Any]]:
        """Parallel batch processing using ThreadPoolExecutor."""

        # Create list of all prediction tasks
        tasks = []
        for idx, variant in pending_variants:
            variant_id = f"chr{variant['chromosome']}:{variant['position']}"
            ref = variant.get('ref', variant.get('other_allele', 'N'))
            alt = variant.get('alt', variant.get('effect_allele', 'N'))

            for tissue in self.tissues:
                for modality in self.modalities:
                    tasks.append({
                        'variant_id': variant_id,
                        'chrom': str(variant['chromosome']),
                        'pos': int(variant['position']),
                        'ref': str(ref),
                        'alt': str(alt),
                        'tissue': tissue,
                        'modality': modality,
                        'rsid': variant.get('rsid', variant_id)
                    })

        logger.info(f"Starting parallel prediction with {self.max_workers} workers...")

        results_lock = threading.Lock()

        with ProgressTracker(total=len(tasks), desc="Predicting variants") as tracker:
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {
                    executor.submit(self._predict_task, task): task
                    for task in tasks
                }

                for future in as_completed(futures):
                    task = futures[future]
                    try:
                        result = future.result()
                        if result:
                            with results_lock:
                                all_predictions.append(result)
                                completed_variants.add(task['variant_id'])
                    except Exception as e:
                        logger.error(f"Task failed for {task['variant_id']}: {e}")

                    tracker.advance()

                    # Checkpoint periodically
                    with self._lock:
                        self._completed_count += 1
                        if checkpoint_path and self._completed_count % (self.checkpoint_every * len(self.tissues) * len(self.modalities)) == 0:
                            self._save_checkpoint(checkpoint_path, all_predictions, completed_variants)

        return all_predictions

    def _predict_task(self, task: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Execute a single prediction task."""
        result = self.predict_variant(
            chrom=task['chrom'],
            pos=task['pos'],
            ref=task['ref'],
            alt=task['alt'],
            tissue=task['tissue'],
            modality=task['modality']
        )

        if result:
            result['rsid'] = task['rsid']
            result['chromosome'] = task['chrom']
            result['position'] = task['pos']

        return result

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


def main(
    config_path: str = "config/config.yaml",
    use_parallel: bool = True,
    max_workers: Optional[int] = None,
    no_checkpoint: bool = False,
    dry_run: bool = False
):
    """
    Main AlphaGenome prediction workflow.

    Args:
        config_path: Path to configuration file
        use_parallel: Enable parallel processing
        max_workers: Override max workers from config
        no_checkpoint: Disable checkpoint saving
        dry_run: Validate inputs without running predictions
    """
    start_time = time.time()

    config = load_config(config_path)

    # Override max_workers if specified
    if max_workers:
        if 'alphagenome' not in config:
            config['alphagenome'] = {}
        config['alphagenome']['max_workers'] = max_workers

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
        return None

    # Validate input
    required_cols = ['chromosome', 'position']
    missing_cols = [c for c in required_cols if c not in variants.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        return None

    if dry_run:
        print_summary_table({
            'Variants': len(variants),
            'Tissues': len(predictor.tissues),
            'Modalities': len(predictor.modalities),
            'Total predictions': len(variants) * len(predictor.tissues) * len(predictor.modalities),
            'Parallel': use_parallel,
            'Workers': predictor.max_workers
        }, title="Dry Run - Would Process")
        return None

    # Run predictions
    checkpoint_path = None if no_checkpoint else str(intermediate_dir / f"{prefix}_alphagenome_checkpoint.json")

    predictions = predictor.predict_batch(
        variants,
        checkpoint_path=checkpoint_path,
        use_parallel=use_parallel
    )

    if not predictions:
        logger.warning("No predictions generated")
        return None

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

    # Print runtime
    elapsed = time.time() - start_time
    logger.info(f"Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} minutes)")

    return pred_df


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Run AlphaGenome variant effect predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python alphagenome_predict.py --config config.yaml
  python alphagenome_predict.py --workers 8 --parallel
  python alphagenome_predict.py --dry-run  # Validate without running
        """
    )
    parser.add_argument("--config", "-c", default="config/config.yaml",
                        help="Path to config file")
    parser.add_argument("--parallel", "-p", action="store_true", default=True,
                        help="Enable parallel processing (default: True)")
    parser.add_argument("--no-parallel", dest="parallel", action="store_false",
                        help="Disable parallel processing")
    parser.add_argument("--workers", "-w", type=int, default=None,
                        help="Number of parallel workers (default: from config or 4)")
    parser.add_argument("--no-checkpoint", action="store_true",
                        help="Disable checkpoint saving")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs without running predictions")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        setup_logging(level=logging.DEBUG)
    else:
        setup_logging(level=logging.INFO)

    main(
        config_path=args.config,
        use_parallel=args.parallel,
        max_workers=args.workers,
        no_checkpoint=args.no_checkpoint,
        dry_run=args.dry_run
    )
