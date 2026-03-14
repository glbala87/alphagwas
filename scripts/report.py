"""
PDF Report Generation Module for AlphaGWAS.

Generates professional PDF reports with:
- Executive summary
- Top variant tables
- Publication-ready figures
- Detailed methodology
"""

import io
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # Non-interactive backend

logger = logging.getLogger(__name__)


@dataclass
class ReportConfig:
    """Configuration for report generation."""

    title: str = "AlphaGWAS Variant Prioritization Report"
    author: str = "AlphaGWAS Pipeline"
    study_name: str = "GWAS Analysis"
    date: str = ""
    logo_path: Optional[Path] = None
    include_methods: bool = True
    include_qc: bool = True
    top_n_variants: int = 20
    figure_dpi: int = 150

    def __post_init__(self):
        if not self.date:
            self.date = datetime.now().strftime("%Y-%m-%d")


class ReportGenerator:
    """
    Generate PDF reports for AlphaGWAS results.

    Uses ReportLab for PDF generation with matplotlib figures.
    """

    def __init__(self, config: Optional[ReportConfig] = None):
        """Initialize report generator."""
        self.config = config or ReportConfig()
        self._check_dependencies()

    def _check_dependencies(self):
        """Check if ReportLab is available."""
        try:
            from reportlab.lib import colors
            from reportlab.lib.pagesizes import letter
            from reportlab.platypus import SimpleDocTemplate

            self.has_reportlab = True
        except ImportError:
            logger.warning(
                "ReportLab not installed. Install with: pip install reportlab"
            )
            self.has_reportlab = False

    def generate_report(
        self,
        ranked_variants: pd.DataFrame,
        tissue_scores: pd.DataFrame,
        output_path: Path,
        gwas_summary: Optional[dict] = None,
        predictions_summary: Optional[dict] = None,
    ) -> Path:
        """
        Generate comprehensive PDF report.

        Args:
            ranked_variants: DataFrame with ranked variants
            tissue_scores: DataFrame with tissue-level scores
            output_path: Path for output PDF
            gwas_summary: Optional GWAS QC summary statistics
            predictions_summary: Optional prediction summary

        Returns:
            Path to generated PDF
        """
        if not self.has_reportlab:
            # Fallback to HTML report
            return self._generate_html_report(
                ranked_variants, tissue_scores, output_path, gwas_summary
            )

        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            Image,
            PageBreak,
            Paragraph,
            SimpleDocTemplate,
            Spacer,
            Table,
            TableStyle,
        )

        # Create document
        doc = SimpleDocTemplate(
            str(output_path),
            pagesize=letter,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=0.75 * inch,
            bottomMargin=0.75 * inch,
        )

        # Styles
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            "CustomTitle",
            parent=styles["Heading1"],
            fontSize=24,
            spaceAfter=30,
            alignment=1,  # Center
        )
        heading_style = ParagraphStyle(
            "CustomHeading",
            parent=styles["Heading2"],
            fontSize=16,
            spaceBefore=20,
            spaceAfter=10,
        )
        body_style = styles["Normal"]

        # Build content
        content = []

        # Title page
        content.append(Spacer(1, 2 * inch))
        content.append(Paragraph(self.config.title, title_style))
        content.append(Spacer(1, 0.5 * inch))
        content.append(
            Paragraph(f"<b>Study:</b> {self.config.study_name}", body_style)
        )
        content.append(Paragraph(f"<b>Date:</b> {self.config.date}", body_style))
        content.append(Paragraph(f"<b>Author:</b> {self.config.author}", body_style))
        content.append(PageBreak())

        # Executive Summary
        content.append(Paragraph("Executive Summary", heading_style))
        summary_text = self._generate_executive_summary(
            ranked_variants, tissue_scores, gwas_summary
        )
        content.append(Paragraph(summary_text, body_style))
        content.append(Spacer(1, 0.3 * inch))

        # Key findings box
        key_findings = self._get_key_findings(ranked_variants, tissue_scores)
        content.append(Paragraph("<b>Key Findings:</b>", body_style))
        for finding in key_findings:
            content.append(Paragraph(f"• {finding}", body_style))
        content.append(Spacer(1, 0.5 * inch))

        # Top Variants Table
        content.append(Paragraph("Top Prioritized Variants", heading_style))
        top_variants_table = self._create_variants_table(
            ranked_variants.head(self.config.top_n_variants)
        )
        content.append(top_variants_table)
        content.append(PageBreak())

        # Visualizations
        content.append(Paragraph("Visualizations", heading_style))

        # Score distribution figure
        fig_score = self._create_score_distribution_figure(ranked_variants)
        if fig_score:
            content.append(fig_score)
            content.append(Spacer(1, 0.3 * inch))

        # Tissue heatmap
        fig_heatmap = self._create_tissue_heatmap_figure(
            tissue_scores, ranked_variants
        )
        if fig_heatmap:
            content.append(fig_heatmap)
            content.append(PageBreak())

        # QC Section
        if self.config.include_qc and gwas_summary:
            content.append(Paragraph("Quality Control Summary", heading_style))
            qc_text = self._generate_qc_summary(gwas_summary)
            content.append(Paragraph(qc_text, body_style))
            content.append(Spacer(1, 0.5 * inch))

        # Methods Section
        if self.config.include_methods:
            content.append(Paragraph("Methods", heading_style))
            methods_text = self._generate_methods_section()
            content.append(Paragraph(methods_text, body_style))
            content.append(PageBreak())

        # Appendix: Full variant list
        content.append(Paragraph("Appendix: Complete Variant Rankings", heading_style))
        content.append(
            Paragraph(
                f"Showing top {min(50, len(ranked_variants))} variants. "
                "Full results available in TSV output files.",
                body_style,
            )
        )
        full_table = self._create_variants_table(
            ranked_variants.head(50), compact=True
        )
        content.append(full_table)

        # Build PDF
        doc.build(content)
        logger.info(f"PDF report generated: {output_path}")

        return output_path

    def _generate_executive_summary(
        self,
        ranked_variants: pd.DataFrame,
        tissue_scores: pd.DataFrame,
        gwas_summary: Optional[dict],
    ) -> str:
        """Generate executive summary text."""
        n_variants = len(ranked_variants)
        top_variant = ranked_variants.iloc[0] if n_variants > 0 else None

        # Get unique tissues
        if "tissue" in tissue_scores.columns:
            n_tissues = tissue_scores["tissue"].nunique()
        else:
            n_tissues = len([c for c in tissue_scores.columns if c not in ["variant_id", "rsid"]])

        summary = f"""
        This report presents the results of variant prioritization analysis using the
        AlphaGWAS pipeline, which integrates AlphaGenome functional predictions with
        GWAS summary statistics.

        <b>Analysis Overview:</b>
        A total of <b>{n_variants}</b> variants were analyzed across <b>{n_tissues}</b>
        tissues and multiple molecular modalities.
        """

        if top_variant is not None:
            top_score = top_variant.get("consensus_score", top_variant.get("score", 0))
            top_id = top_variant.get("variant_id", "Unknown")
            summary += f"""

            The top-ranked variant is <b>{top_id}</b> with a consensus score of
            <b>{top_score:.4f}</b>.
            """

        if gwas_summary:
            n_sig = gwas_summary.get("n_significant", 0)
            n_loci = gwas_summary.get("n_loci", 0)
            summary += f"""

            The input GWAS contained <b>{n_sig}</b> genome-wide significant variants
            across <b>{n_loci}</b> independent loci.
            """

        return summary

    def _get_key_findings(
        self, ranked_variants: pd.DataFrame, tissue_scores: pd.DataFrame
    ) -> list[str]:
        """Extract key findings from results."""
        findings = []

        # Top variant
        if len(ranked_variants) > 0:
            top = ranked_variants.iloc[0]
            findings.append(
                f"Top prioritized variant: {top.get('variant_id', 'N/A')} "
                f"(score: {top.get('consensus_score', 0):.3f})"
            )

        # High-confidence variants
        if "consensus_score" in ranked_variants.columns:
            high_conf = (ranked_variants["consensus_score"] > 0.5).sum()
            findings.append(f"{high_conf} variants with consensus score > 0.5")

        # Top tissue
        if "tissue" in tissue_scores.columns and "effect_size" in tissue_scores.columns:
            tissue_means = tissue_scores.groupby("tissue")["effect_size"].mean()
            top_tissue = tissue_means.idxmax()
            findings.append(f"Most impacted tissue: {top_tissue}")

        # Multi-tissue effects
        if "n_significant_tissues" in ranked_variants.columns:
            multi_tissue = (ranked_variants["n_significant_tissues"] >= 3).sum()
            findings.append(f"{multi_tissue} variants affect 3+ tissues")

        return findings

    def _create_variants_table(
        self, df: pd.DataFrame, compact: bool = False
    ) -> Any:
        """Create formatted table for variants."""
        from reportlab.lib import colors
        from reportlab.lib.units import inch
        from reportlab.platypus import Table, TableStyle

        # Select columns to display
        display_cols = ["variant_id", "consensus_score", "final_rank"]

        if "rsid" in df.columns:
            display_cols.insert(1, "rsid")
        if "max_effect" in df.columns and not compact:
            display_cols.append("max_effect")
        if "n_significant_tissues" in df.columns and not compact:
            display_cols.append("n_significant_tissues")
        if "top_tissue" in df.columns and not compact:
            display_cols.append("top_tissue")

        # Filter to available columns
        display_cols = [c for c in display_cols if c in df.columns]

        # Prepare data
        data = [display_cols]  # Header row

        for _, row in df.iterrows():
            row_data = []
            for col in display_cols:
                val = row[col]
                if isinstance(val, float):
                    val = f"{val:.4f}" if col != "final_rank" else f"{int(val)}"
                row_data.append(str(val)[:20])  # Truncate long values
            data.append(row_data)

        # Create table
        col_widths = [1.5 * inch] + [1.0 * inch] * (len(display_cols) - 1)
        if compact:
            col_widths = [1.2 * inch] * len(display_cols)

        table = Table(data, colWidths=col_widths)

        # Style
        style = TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2E86AB")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, 0), 9 if compact else 10),
                ("FONTSIZE", (0, 1), (-1, -1), 8 if compact else 9),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
                ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#F5F5F5")),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F0F0F0")]),
            ]
        )
        table.setStyle(style)

        return table

    def _create_score_distribution_figure(self, ranked_variants: pd.DataFrame) -> Optional[Any]:
        """Create score distribution figure."""
        from reportlab.lib.units import inch
        from reportlab.platypus import Image

        if "consensus_score" not in ranked_variants.columns:
            return None

        fig, ax = plt.subplots(figsize=(6, 4))

        scores = ranked_variants["consensus_score"].dropna()
        ax.hist(scores, bins=30, edgecolor="black", alpha=0.7, color="#2E86AB")
        ax.axvline(
            x=scores.median(), color="red", linestyle="--", label=f"Median: {scores.median():.3f}"
        )
        ax.set_xlabel("Consensus Score", fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.set_title("Distribution of Variant Consensus Scores", fontsize=12)
        ax.legend()

        # Save to buffer
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=self.config.figure_dpi, bbox_inches="tight")
        buf.seek(0)
        plt.close(fig)

        return Image(buf, width=5 * inch, height=3.5 * inch)

    def _create_tissue_heatmap_figure(
        self, tissue_scores: pd.DataFrame, ranked_variants: pd.DataFrame
    ) -> Optional[Any]:
        """Create tissue heatmap figure."""
        from reportlab.lib.units import inch
        from reportlab.platypus import Image

        try:
            # Get top variants
            top_variants = ranked_variants.head(15)["variant_id"].tolist()

            # Pivot tissue scores
            if "variant_id" in tissue_scores.columns and "tissue" in tissue_scores.columns:
                score_col = "effect_size" if "effect_size" in tissue_scores.columns else "score"
                if score_col not in tissue_scores.columns:
                    return None

                pivot = tissue_scores[tissue_scores["variant_id"].isin(top_variants)].pivot_table(
                    index="variant_id", columns="tissue", values=score_col, aggfunc="mean"
                )
            else:
                # Already pivoted format
                pivot = tissue_scores[tissue_scores.index.isin(top_variants)]

            if pivot.empty:
                return None

            fig, ax = plt.subplots(figsize=(8, 6))
            im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd")

            ax.set_xticks(range(len(pivot.columns)))
            ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(len(pivot.index)))
            ax.set_yticklabels(pivot.index, fontsize=8)

            plt.colorbar(im, ax=ax, label="Effect Size")
            ax.set_title("Tissue-Specific Effects (Top Variants)", fontsize=12)

            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=self.config.figure_dpi, bbox_inches="tight")
            buf.seek(0)
            plt.close(fig)

            return Image(buf, width=6 * inch, height=4.5 * inch)

        except Exception as e:
            logger.warning(f"Could not create tissue heatmap: {e}")
            return None

    def _generate_qc_summary(self, gwas_summary: dict) -> str:
        """Generate QC summary text."""
        qc_text = "<b>Input Data Quality:</b><br/>"

        if "n_variants_input" in gwas_summary:
            qc_text += f"• Input variants: {gwas_summary['n_variants_input']:,}<br/>"
        if "n_significant" in gwas_summary:
            qc_text += f"• Genome-wide significant: {gwas_summary['n_significant']:,}<br/>"
        if "n_loci" in gwas_summary:
            qc_text += f"• Independent loci: {gwas_summary['n_loci']}<br/>"
        if "lambda_gc" in gwas_summary:
            qc_text += f"• Genomic inflation (λ): {gwas_summary['lambda_gc']:.3f}<br/>"
        if "n_with_predictions" in gwas_summary:
            qc_text += f"• Variants with predictions: {gwas_summary['n_with_predictions']:,}<br/>"

        return qc_text

    def _generate_methods_section(self) -> str:
        """Generate methods section text."""
        return """
        <b>Variant Prioritization Pipeline:</b><br/><br/>

        1. <b>Variant Extraction:</b> Genome-wide significant variants (p < 5×10⁻⁸)
        were extracted from GWAS summary statistics. Lead SNPs were identified using
        LD-based clumping (r² < 0.1, 500kb window).<br/><br/>

        2. <b>LD Proxy Identification:</b> Variants in high LD (r² > 0.8) with lead
        SNPs were identified using the 1000 Genomes reference panel.<br/><br/>

        3. <b>Functional Predictions:</b> AlphaGenome was used to predict variant
        effects on gene expression and chromatin accessibility across multiple tissues.<br/><br/>

        4. <b>Scoring and Ranking:</b> A consensus score was calculated by aggregating
        predictions across tissues using a weighted mean approach. Variants were ranked
        by their consensus functional impact scores.<br/><br/>

        <b>Software:</b> AlphaGWAS pipeline (https://github.com/glbala87/alphagwas)<br/>
        <b>AlphaGenome:</b> Google DeepMind variant effect prediction model
        """

    def _generate_html_report(
        self,
        ranked_variants: pd.DataFrame,
        tissue_scores: pd.DataFrame,
        output_path: Path,
        gwas_summary: Optional[dict],
    ) -> Path:
        """Fallback HTML report when ReportLab is not available."""
        html_path = output_path.with_suffix(".html")

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>{self.config.title}</title>
    <style>
        body {{ font-family: Arial, sans-serif; max-width: 1000px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #2E86AB; }}
        h2 {{ color: #333; border-bottom: 2px solid #2E86AB; padding-bottom: 5px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #2E86AB; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .summary-box {{ background: #f5f5f5; padding: 15px; border-radius: 5px; margin: 20px 0; }}
        .key-finding {{ margin: 5px 0; }}
    </style>
</head>
<body>
    <h1>{self.config.title}</h1>
    <p><strong>Study:</strong> {self.config.study_name}</p>
    <p><strong>Date:</strong> {self.config.date}</p>

    <h2>Executive Summary</h2>
    <div class="summary-box">
        <p>Analyzed <strong>{len(ranked_variants)}</strong> variants.</p>
        {"<p>Top variant: <strong>" + str(ranked_variants.iloc[0].get('variant_id', 'N/A')) +
         "</strong> (score: " + f"{ranked_variants.iloc[0].get('consensus_score', 0):.4f}" + ")</p>"
         if len(ranked_variants) > 0 else ""}
    </div>

    <h2>Top Variants</h2>
    {ranked_variants.head(20).to_html(index=False, classes='variants-table')}

    <h2>Methods</h2>
    <p>Variant prioritization performed using AlphaGWAS pipeline with AlphaGenome predictions.</p>

    <footer>
        <p><em>Generated by AlphaGWAS on {self.config.date}</em></p>
    </footer>
</body>
</html>
        """

        with open(html_path, "w") as f:
            f.write(html_content)

        logger.info(f"HTML report generated: {html_path}")
        return html_path


def generate_report(
    ranked_variants_path: Path,
    tissue_scores_path: Path,
    output_path: Path,
    config: Optional[dict] = None,
) -> Path:
    """
    Convenience function to generate report from file paths.

    Args:
        ranked_variants_path: Path to ranked variants TSV
        tissue_scores_path: Path to tissue scores TSV
        output_path: Path for output PDF
        config: Optional configuration dict

    Returns:
        Path to generated report
    """
    ranked = pd.read_csv(ranked_variants_path, sep="\t")
    tissue_scores = pd.read_csv(tissue_scores_path, sep="\t")

    report_config = ReportConfig(**(config or {}))
    generator = ReportGenerator(report_config)

    return generator.generate_report(ranked, tissue_scores, output_path)
