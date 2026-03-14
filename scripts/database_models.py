"""
SQLAlchemy ORM Models for AlphaGWAS Database.

Defines the database schema for storing:
- Analyses (jobs/runs)
- Variants and their scores
- AlphaGenome predictions
- Tissue scores
"""

from datetime import datetime
from typing import Optional

from sqlalchemy import (
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class Analysis(Base):
    """Analysis/job record."""

    __tablename__ = "analyses"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False)
    description = Column(Text, default="")
    status = Column(String(50), default="created")  # created, running, completed, failed
    config = Column(Text, default="{}")  # JSON configuration
    metadata = Column(Text, default="{}")  # Additional metadata

    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    variants = relationship("Variant", back_populates="analysis", cascade="all, delete-orphan")
    predictions = relationship("Prediction", back_populates="analysis", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Analysis(id={self.id}, name='{self.name}', status='{self.status}')>"


class Variant(Base):
    """Variant record with scores."""

    __tablename__ = "variants"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=False)

    # Variant identification
    variant_id = Column(String(100), nullable=False)  # chr:pos format
    chromosome = Column(String(10), nullable=False)
    position = Column(Integer, nullable=False)
    rsid = Column(String(50), nullable=True)
    ref_allele = Column(String(500), nullable=True)
    alt_allele = Column(String(500), nullable=True)

    # GWAS statistics
    pvalue = Column(Float, nullable=True)
    beta = Column(Float, nullable=True)
    se = Column(Float, nullable=True)
    maf = Column(Float, nullable=True)

    # AlphaGWAS scores
    consensus_score = Column(Float, nullable=True)
    max_effect = Column(Float, nullable=True)
    n_significant_tissues = Column(Integer, nullable=True)
    rank = Column(Integer, nullable=True)

    # Classification
    variant_type = Column(String(50), default="ranked")  # ranked, significant, lead, proxy

    # Additional data as JSON
    metadata = Column(Text, default="{}")

    # Relationships
    analysis = relationship("Analysis", back_populates="variants")

    # Indexes for common queries
    __table_args__ = (
        Index("idx_variant_analysis", "analysis_id"),
        Index("idx_variant_chrom_pos", "chromosome", "position"),
        Index("idx_variant_rsid", "rsid"),
        Index("idx_variant_score", "consensus_score"),
        Index("idx_variant_rank", "rank"),
    )

    def __repr__(self):
        return f"<Variant(id={self.id}, variant_id='{self.variant_id}', score={self.consensus_score})>"


class Prediction(Base):
    """AlphaGenome prediction record."""

    __tablename__ = "predictions"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=False)

    # Variant reference
    variant_id = Column(String(100), nullable=False)

    # Prediction details
    tissue = Column(String(100), nullable=False)
    modality = Column(String(50), nullable=False)  # expression, chromatin_accessibility
    effect_size = Column(Float, nullable=False)
    confidence = Column(Float, nullable=True)
    prediction_value = Column(Float, nullable=True)

    # Relationships
    analysis = relationship("Analysis", back_populates="predictions")

    # Indexes
    __table_args__ = (
        Index("idx_pred_analysis", "analysis_id"),
        Index("idx_pred_variant", "variant_id"),
        Index("idx_pred_tissue", "tissue"),
        Index("idx_pred_analysis_variant", "analysis_id", "variant_id"),
    )

    def __repr__(self):
        return f"<Prediction(variant='{self.variant_id}', tissue='{self.tissue}', effect={self.effect_size})>"


class TissueScore(Base):
    """Tissue-level aggregated scores."""

    __tablename__ = "tissue_scores"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=False)

    variant_id = Column(String(100), nullable=False)
    tissue = Column(String(100), nullable=False)

    # Aggregated scores
    mean_effect = Column(Float, nullable=True)
    max_effect = Column(Float, nullable=True)
    expression_effect = Column(Float, nullable=True)
    chromatin_effect = Column(Float, nullable=True)

    __table_args__ = (
        Index("idx_tissue_analysis_variant", "analysis_id", "variant_id"),
    )


class PRSModel(Base):
    """Polygenic Risk Score model."""

    __tablename__ = "prs_models"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=False)

    name = Column(String(255), nullable=False)
    n_variants = Column(Integer, nullable=False)
    method = Column(String(50), default="clumping")  # clumping, pruning, bayesian

    # Performance metrics
    r2 = Column(Float, nullable=True)
    auc = Column(Float, nullable=True)
    liability_r2 = Column(Float, nullable=True)

    # Model weights stored as JSON
    weights = Column(Text, default="{}")
    metadata = Column(Text, default="{}")

    created_at = Column(DateTime, default=datetime.utcnow)


class MRResult(Base):
    """Mendelian Randomization result."""

    __tablename__ = "mr_results"

    id = Column(Integer, primary_key=True, autoincrement=True)
    analysis_id = Column(Integer, ForeignKey("analyses.id"), nullable=True)

    exposure = Column(String(255), nullable=False)
    outcome = Column(String(255), nullable=False)
    method = Column(String(50), nullable=False)  # IVW, MR-Egger, weighted_median

    n_snps = Column(Integer, nullable=False)
    beta = Column(Float, nullable=False)
    se = Column(Float, nullable=False)
    pvalue = Column(Float, nullable=False)

    # Sensitivity analyses
    egger_intercept = Column(Float, nullable=True)
    egger_intercept_p = Column(Float, nullable=True)
    heterogeneity_q = Column(Float, nullable=True)
    heterogeneity_p = Column(Float, nullable=True)

    metadata = Column(Text, default="{}")
    created_at = Column(DateTime, default=datetime.utcnow)

    __table_args__ = (
        Index("idx_mr_exposure_outcome", "exposure", "outcome"),
    )
