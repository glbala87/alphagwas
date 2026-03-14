"""
Database Backend Module for AlphaGWAS.

Provides persistent storage for:
- GWAS analysis results
- Variant predictions and scores
- Job history and metadata
- Query and export capabilities

Supports SQLite (local) and PostgreSQL (production).
"""

import json
import logging
from contextlib import contextmanager
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Generator, Optional

import pandas as pd

logger = logging.getLogger(__name__)


class DatabaseType(str, Enum):
    """Supported database types."""

    SQLITE = "sqlite"
    POSTGRESQL = "postgresql"


class AlphaGWASDatabase:
    """
    Database interface for AlphaGWAS results storage.

    Supports both SQLite for local development and PostgreSQL for production.
    """

    def __init__(
        self,
        db_url: str = "sqlite:///alphagwas.db",
        echo: bool = False,
    ):
        """
        Initialize database connection.

        Args:
            db_url: Database URL (sqlite:///path.db or postgresql://user:pass@host/db)
            echo: Echo SQL statements for debugging
        """
        self.db_url = db_url
        self.echo = echo
        self.engine = None
        self.Session = None
        self._init_db()

    def _init_db(self):
        """Initialize database engine and create tables."""
        try:
            from sqlalchemy import create_engine
            from sqlalchemy.orm import sessionmaker

            self.engine = create_engine(self.db_url, echo=self.echo)
            self.Session = sessionmaker(bind=self.engine)

            # Create tables
            from .database_models import Base

            Base.metadata.create_all(self.engine)
            logger.info(f"Database initialized: {self.db_url}")

        except ImportError:
            logger.warning(
                "SQLAlchemy not installed. Install with: pip install sqlalchemy"
            )
            self.engine = None

    @contextmanager
    def session_scope(self) -> Generator:
        """Provide transactional scope for database operations."""
        if self.Session is None:
            raise RuntimeError("Database not initialized")

        session = self.Session()
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

    # ==================== Analysis Operations ====================

    def create_analysis(
        self,
        name: str,
        description: str = "",
        config: Optional[dict] = None,
        metadata: Optional[dict] = None,
    ) -> int:
        """
        Create a new analysis record.

        Args:
            name: Analysis name
            description: Analysis description
            config: Pipeline configuration
            metadata: Additional metadata

        Returns:
            Analysis ID
        """
        from .database_models import Analysis

        with self.session_scope() as session:
            analysis = Analysis(
                name=name,
                description=description,
                config=json.dumps(config or {}),
                metadata=json.dumps(metadata or {}),
                status="created",
                created_at=datetime.utcnow(),
            )
            session.add(analysis)
            session.flush()
            analysis_id = analysis.id
            logger.info(f"Created analysis: {name} (ID: {analysis_id})")
            return analysis_id

    def update_analysis_status(
        self,
        analysis_id: int,
        status: str,
        message: str = "",
    ):
        """Update analysis status."""
        from .database_models import Analysis

        with self.session_scope() as session:
            analysis = session.query(Analysis).get(analysis_id)
            if analysis:
                analysis.status = status
                analysis.updated_at = datetime.utcnow()
                if status == "completed":
                    analysis.completed_at = datetime.utcnow()
                logger.info(f"Analysis {analysis_id} status: {status}")

    def get_analysis(self, analysis_id: int) -> Optional[dict]:
        """Get analysis by ID."""
        from .database_models import Analysis

        with self.session_scope() as session:
            analysis = session.query(Analysis).get(analysis_id)
            if analysis:
                return {
                    "id": analysis.id,
                    "name": analysis.name,
                    "description": analysis.description,
                    "status": analysis.status,
                    "config": json.loads(analysis.config),
                    "metadata": json.loads(analysis.metadata),
                    "created_at": analysis.created_at.isoformat(),
                    "completed_at": (
                        analysis.completed_at.isoformat()
                        if analysis.completed_at
                        else None
                    ),
                }
        return None

    def list_analyses(
        self,
        status: Optional[str] = None,
        limit: int = 100,
        offset: int = 0,
    ) -> list[dict]:
        """List analyses with optional filtering."""
        from .database_models import Analysis

        with self.session_scope() as session:
            query = session.query(Analysis)
            if status:
                query = query.filter(Analysis.status == status)
            query = query.order_by(Analysis.created_at.desc())
            query = query.offset(offset).limit(limit)

            return [
                {
                    "id": a.id,
                    "name": a.name,
                    "status": a.status,
                    "created_at": a.created_at.isoformat(),
                }
                for a in query.all()
            ]

    # ==================== Variant Operations ====================

    def store_variants(
        self,
        analysis_id: int,
        variants_df: pd.DataFrame,
        variant_type: str = "ranked",
    ):
        """
        Store variants for an analysis.

        Args:
            analysis_id: Analysis ID
            variants_df: DataFrame with variant data
            variant_type: Type of variants (ranked, significant, all)
        """
        from .database_models import Variant

        with self.session_scope() as session:
            for _, row in variants_df.iterrows():
                variant = Variant(
                    analysis_id=analysis_id,
                    variant_id=row.get("variant_id", ""),
                    chromosome=str(row.get("chromosome", "")),
                    position=int(row.get("position", 0)),
                    rsid=row.get("rsid"),
                    ref_allele=row.get("ref_allele") or row.get("other_allele"),
                    alt_allele=row.get("alt_allele") or row.get("effect_allele"),
                    pvalue=float(row.get("pvalue", 1.0)) if pd.notna(row.get("pvalue")) else None,
                    beta=float(row.get("beta", 0)) if pd.notna(row.get("beta")) else None,
                    consensus_score=float(row.get("consensus_score", 0)) if pd.notna(row.get("consensus_score")) else None,
                    rank=int(row.get("final_rank", 0)) if pd.notna(row.get("final_rank")) else None,
                    variant_type=variant_type,
                    metadata=json.dumps(row.to_dict(), default=str),
                )
                session.add(variant)

            logger.info(f"Stored {len(variants_df)} variants for analysis {analysis_id}")

    def get_variants(
        self,
        analysis_id: int,
        variant_type: Optional[str] = None,
        min_score: Optional[float] = None,
        max_rank: Optional[int] = None,
        chromosome: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Query variants with filtering.

        Args:
            analysis_id: Analysis ID
            variant_type: Filter by type
            min_score: Minimum consensus score
            max_rank: Maximum rank (top N)
            chromosome: Filter by chromosome

        Returns:
            DataFrame with variants
        """
        from .database_models import Variant

        with self.session_scope() as session:
            query = session.query(Variant).filter(Variant.analysis_id == analysis_id)

            if variant_type:
                query = query.filter(Variant.variant_type == variant_type)
            if min_score is not None:
                query = query.filter(Variant.consensus_score >= min_score)
            if max_rank is not None:
                query = query.filter(Variant.rank <= max_rank)
            if chromosome:
                query = query.filter(Variant.chromosome == chromosome)

            query = query.order_by(Variant.rank)

            variants = query.all()

            if not variants:
                return pd.DataFrame()

            return pd.DataFrame(
                [
                    {
                        "variant_id": v.variant_id,
                        "chromosome": v.chromosome,
                        "position": v.position,
                        "rsid": v.rsid,
                        "pvalue": v.pvalue,
                        "beta": v.beta,
                        "consensus_score": v.consensus_score,
                        "rank": v.rank,
                    }
                    for v in variants
                ]
            )

    # ==================== Prediction Operations ====================

    def store_predictions(
        self,
        analysis_id: int,
        predictions_df: pd.DataFrame,
    ):
        """Store AlphaGenome predictions."""
        from .database_models import Prediction

        with self.session_scope() as session:
            for _, row in predictions_df.iterrows():
                prediction = Prediction(
                    analysis_id=analysis_id,
                    variant_id=row.get("variant_id", ""),
                    tissue=row.get("tissue", ""),
                    modality=row.get("modality", ""),
                    effect_size=float(row.get("effect_size", 0)),
                    confidence=float(row.get("confidence", 0)) if pd.notna(row.get("confidence")) else None,
                    prediction_value=float(row.get("prediction", 0)) if pd.notna(row.get("prediction")) else None,
                )
                session.add(prediction)

            logger.info(f"Stored {len(predictions_df)} predictions for analysis {analysis_id}")

    def get_predictions(
        self,
        analysis_id: int,
        variant_id: Optional[str] = None,
        tissue: Optional[str] = None,
    ) -> pd.DataFrame:
        """Query predictions with filtering."""
        from .database_models import Prediction

        with self.session_scope() as session:
            query = session.query(Prediction).filter(
                Prediction.analysis_id == analysis_id
            )

            if variant_id:
                query = query.filter(Prediction.variant_id == variant_id)
            if tissue:
                query = query.filter(Prediction.tissue == tissue)

            predictions = query.all()

            return pd.DataFrame(
                [
                    {
                        "variant_id": p.variant_id,
                        "tissue": p.tissue,
                        "modality": p.modality,
                        "effect_size": p.effect_size,
                        "confidence": p.confidence,
                    }
                    for p in predictions
                ]
            )

    # ==================== Export Operations ====================

    def export_analysis(
        self,
        analysis_id: int,
        output_dir: Path,
        formats: list[str] = ["tsv", "json"],
    ) -> dict[str, Path]:
        """
        Export analysis results to files.

        Args:
            analysis_id: Analysis ID
            output_dir: Output directory
            formats: Export formats (tsv, json, parquet)

        Returns:
            Dictionary of format -> file path
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        exported = {}

        # Get data
        variants = self.get_variants(analysis_id)
        predictions = self.get_predictions(analysis_id)
        analysis = self.get_analysis(analysis_id)

        # Export variants
        if not variants.empty:
            for fmt in formats:
                if fmt == "tsv":
                    path = output_dir / f"analysis_{analysis_id}_variants.tsv"
                    variants.to_csv(path, sep="\t", index=False)
                elif fmt == "json":
                    path = output_dir / f"analysis_{analysis_id}_variants.json"
                    variants.to_json(path, orient="records", indent=2)
                elif fmt == "parquet":
                    path = output_dir / f"analysis_{analysis_id}_variants.parquet"
                    variants.to_parquet(path, index=False)
                exported[f"variants_{fmt}"] = path

        # Export predictions
        if not predictions.empty and "parquet" in formats:
            path = output_dir / f"analysis_{analysis_id}_predictions.parquet"
            predictions.to_parquet(path, index=False)
            exported["predictions_parquet"] = path

        # Export metadata
        if analysis and "json" in formats:
            path = output_dir / f"analysis_{analysis_id}_metadata.json"
            with open(path, "w") as f:
                json.dump(analysis, f, indent=2)
            exported["metadata_json"] = path

        logger.info(f"Exported analysis {analysis_id} to {output_dir}")
        return exported

    # ==================== Statistics ====================

    def get_analysis_stats(self, analysis_id: int) -> dict:
        """Get statistics for an analysis."""
        from .database_models import Prediction, Variant

        with self.session_scope() as session:
            n_variants = (
                session.query(Variant)
                .filter(Variant.analysis_id == analysis_id)
                .count()
            )
            n_predictions = (
                session.query(Prediction)
                .filter(Prediction.analysis_id == analysis_id)
                .count()
            )

            # Get score distribution
            variants = self.get_variants(analysis_id)
            stats = {
                "n_variants": n_variants,
                "n_predictions": n_predictions,
            }

            if not variants.empty and "consensus_score" in variants.columns:
                scores = variants["consensus_score"].dropna()
                stats.update(
                    {
                        "mean_score": float(scores.mean()),
                        "median_score": float(scores.median()),
                        "max_score": float(scores.max()),
                        "min_score": float(scores.min()),
                    }
                )

            return stats

    def get_global_stats(self) -> dict:
        """Get global database statistics."""
        from .database_models import Analysis, Prediction, Variant

        with self.session_scope() as session:
            return {
                "total_analyses": session.query(Analysis).count(),
                "completed_analyses": (
                    session.query(Analysis)
                    .filter(Analysis.status == "completed")
                    .count()
                ),
                "total_variants": session.query(Variant).count(),
                "total_predictions": session.query(Prediction).count(),
            }


# Convenience function for quick database setup
def get_database(db_url: Optional[str] = None) -> AlphaGWASDatabase:
    """
    Get database instance.

    Args:
        db_url: Database URL. If None, uses SQLite in current directory.

    Returns:
        AlphaGWASDatabase instance
    """
    if db_url is None:
        db_url = "sqlite:///alphagwas.db"
    return AlphaGWASDatabase(db_url)
