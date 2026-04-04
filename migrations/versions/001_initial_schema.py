"""Initial schema - all AlphaGWAS tables.

Revision ID: 001_initial
Revises: None
Create Date: 2026-04-04
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa

revision: str = "001_initial"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # Analyses table
    op.create_table(
        "analyses",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("description", sa.Text(), server_default=""),
        sa.Column("status", sa.String(50), server_default="created"),
        sa.Column("config", sa.Text(), server_default="{}"),
        sa.Column("metadata", sa.Text(), server_default="{}"),
        sa.Column("created_at", sa.DateTime()),
        sa.Column("updated_at", sa.DateTime()),
        sa.Column("completed_at", sa.DateTime(), nullable=True),
    )

    # Variants table
    op.create_table(
        "variants",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("analysis_id", sa.Integer(), sa.ForeignKey("analyses.id"), nullable=False),
        sa.Column("variant_id", sa.String(100), nullable=False),
        sa.Column("chromosome", sa.String(10), nullable=False),
        sa.Column("position", sa.Integer(), nullable=False),
        sa.Column("rsid", sa.String(50), nullable=True),
        sa.Column("ref_allele", sa.String(500), nullable=True),
        sa.Column("alt_allele", sa.String(500), nullable=True),
        sa.Column("pvalue", sa.Float(), nullable=True),
        sa.Column("beta", sa.Float(), nullable=True),
        sa.Column("se", sa.Float(), nullable=True),
        sa.Column("maf", sa.Float(), nullable=True),
        sa.Column("consensus_score", sa.Float(), nullable=True),
        sa.Column("max_effect", sa.Float(), nullable=True),
        sa.Column("n_significant_tissues", sa.Integer(), nullable=True),
        sa.Column("rank", sa.Integer(), nullable=True),
        sa.Column("variant_type", sa.String(50), server_default="ranked"),
        sa.Column("metadata", sa.Text(), server_default="{}"),
    )
    op.create_index("idx_variant_analysis", "variants", ["analysis_id"])
    op.create_index("idx_variant_chrom_pos", "variants", ["chromosome", "position"])
    op.create_index("idx_variant_rsid", "variants", ["rsid"])
    op.create_index("idx_variant_score", "variants", ["consensus_score"])
    op.create_index("idx_variant_rank", "variants", ["rank"])

    # Predictions table
    op.create_table(
        "predictions",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("analysis_id", sa.Integer(), sa.ForeignKey("analyses.id"), nullable=False),
        sa.Column("variant_id", sa.String(100), nullable=False),
        sa.Column("tissue", sa.String(100), nullable=False),
        sa.Column("modality", sa.String(50), nullable=False),
        sa.Column("effect_size", sa.Float(), nullable=False),
        sa.Column("confidence", sa.Float(), nullable=True),
        sa.Column("prediction_value", sa.Float(), nullable=True),
    )
    op.create_index("idx_pred_analysis", "predictions", ["analysis_id"])
    op.create_index("idx_pred_variant", "predictions", ["variant_id"])
    op.create_index("idx_pred_tissue", "predictions", ["tissue"])
    op.create_index("idx_pred_analysis_variant", "predictions", ["analysis_id", "variant_id"])

    # Tissue scores table
    op.create_table(
        "tissue_scores",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("analysis_id", sa.Integer(), sa.ForeignKey("analyses.id"), nullable=False),
        sa.Column("variant_id", sa.String(100), nullable=False),
        sa.Column("tissue", sa.String(100), nullable=False),
        sa.Column("mean_effect", sa.Float(), nullable=True),
        sa.Column("max_effect", sa.Float(), nullable=True),
        sa.Column("expression_effect", sa.Float(), nullable=True),
        sa.Column("chromatin_effect", sa.Float(), nullable=True),
    )
    op.create_index("idx_tissue_analysis_variant", "tissue_scores", ["analysis_id", "variant_id"])

    # PRS models table
    op.create_table(
        "prs_models",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("analysis_id", sa.Integer(), sa.ForeignKey("analyses.id"), nullable=False),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("n_variants", sa.Integer(), nullable=False),
        sa.Column("method", sa.String(50), server_default="clumping"),
        sa.Column("r2", sa.Float(), nullable=True),
        sa.Column("auc", sa.Float(), nullable=True),
        sa.Column("liability_r2", sa.Float(), nullable=True),
        sa.Column("weights", sa.Text(), server_default="{}"),
        sa.Column("metadata", sa.Text(), server_default="{}"),
        sa.Column("created_at", sa.DateTime()),
    )

    # MR results table
    op.create_table(
        "mr_results",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("analysis_id", sa.Integer(), sa.ForeignKey("analyses.id"), nullable=True),
        sa.Column("exposure", sa.String(255), nullable=False),
        sa.Column("outcome", sa.String(255), nullable=False),
        sa.Column("method", sa.String(50), nullable=False),
        sa.Column("n_snps", sa.Integer(), nullable=False),
        sa.Column("beta", sa.Float(), nullable=False),
        sa.Column("se", sa.Float(), nullable=False),
        sa.Column("pvalue", sa.Float(), nullable=False),
        sa.Column("egger_intercept", sa.Float(), nullable=True),
        sa.Column("egger_intercept_p", sa.Float(), nullable=True),
        sa.Column("heterogeneity_q", sa.Float(), nullable=True),
        sa.Column("heterogeneity_p", sa.Float(), nullable=True),
        sa.Column("metadata", sa.Text(), server_default="{}"),
        sa.Column("created_at", sa.DateTime()),
    )
    op.create_index("idx_mr_exposure_outcome", "mr_results", ["exposure", "outcome"])

    # API jobs table (for persistent job storage)
    op.create_table(
        "api_jobs",
        sa.Column("job_id", sa.String(8), primary_key=True),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("status", sa.String(20), nullable=False, server_default="pending"),
        sa.Column("created_at", sa.String(50), nullable=False),
        sa.Column("updated_at", sa.String(50), nullable=False),
        sa.Column("progress", sa.Float(), server_default="0"),
        sa.Column("current_step", sa.String(50), nullable=True),
        sa.Column("message", sa.Text(), nullable=True),
        sa.Column("results_url", sa.String(255), nullable=True),
        sa.Column("config_json", sa.Text(), server_default="{}"),
        sa.Column("input_file", sa.String(500), nullable=True),
        sa.Column("results_json", sa.Text(), nullable=True),
    )


def downgrade() -> None:
    op.drop_table("api_jobs")
    op.drop_table("mr_results")
    op.drop_table("prs_models")
    op.drop_table("tissue_scores")
    op.drop_table("predictions")
    op.drop_table("variants")
    op.drop_table("analyses")
