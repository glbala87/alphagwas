"""
FastAPI REST Interface for AlphaGWAS Pipeline.

Provides endpoints for:
- Submitting GWAS data for analysis
- Running individual pipeline steps
- Retrieving results and visualizations
- Job status monitoring
"""

import asyncio
import hashlib
import tempfile
import uuid
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from fastapi import BackgroundTasks, FastAPI, File, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field

from . import extract_variants, score_variants, utils

# Initialize FastAPI app
app = FastAPI(
    title="AlphaGWAS API",
    description="REST API for GWAS variant prioritization using AlphaGenome predictions",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job storage (use Redis/database in production)
jobs: dict[str, dict] = {}
results_cache: dict[str, Any] = {}

# Configuration
UPLOAD_DIR = Path(tempfile.gettempdir()) / "alphagwas_uploads"
RESULTS_DIR = Path(tempfile.gettempdir()) / "alphagwas_results"
UPLOAD_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)


class JobStatus(str, Enum):
    """Job status enumeration."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class PipelineStep(str, Enum):
    """Pipeline step enumeration."""

    EXTRACT = "extract"
    LD_PROXIES = "ld_proxies"
    LIFTOVER = "liftover"
    PREDICT = "predict"
    SCORE = "score"
    VISUALIZE = "visualize"
    ALL = "all"


# Request/Response Models
class JobSubmission(BaseModel):
    """Job submission request."""

    name: str = Field(..., description="Name for this analysis job")
    pvalue_threshold: float = Field(5e-8, description="P-value significance threshold")
    ld_r2_threshold: float = Field(0.8, description="LD R² threshold for proxy variants")
    tissues: list[str] = Field(
        default_factory=lambda: ["Whole_Blood", "Liver", "Brain_Cortex"],
        description="Tissues for AlphaGenome predictions",
    )
    consensus_method: str = Field("weighted_mean", description="Scoring method")
    steps: list[PipelineStep] = Field(
        default_factory=lambda: [PipelineStep.ALL], description="Pipeline steps to run"
    )


class JobResponse(BaseModel):
    """Job response model."""

    job_id: str
    name: str
    status: JobStatus
    created_at: str
    updated_at: str
    progress: float = 0.0
    current_step: Optional[str] = None
    message: Optional[str] = None
    results_url: Optional[str] = None


class VariantQuery(BaseModel):
    """Variant query request."""

    chromosome: str
    position: int
    ref: Optional[str] = None
    alt: Optional[str] = None


class VariantScore(BaseModel):
    """Variant score response."""

    variant_id: str
    chromosome: str
    position: int
    rsid: Optional[str] = None
    consensus_score: float
    rank: int
    top_tissues: list[dict]
    effect_sizes: dict[str, float]


class HealthResponse(BaseModel):
    """Health check response."""

    status: str
    version: str
    timestamp: str


# Helper functions
def generate_job_id() -> str:
    """Generate unique job ID."""
    return str(uuid.uuid4())[:8]


def get_file_hash(content: bytes) -> str:
    """Generate hash for file content."""
    return hashlib.md5(content).hexdigest()[:12]


async def run_pipeline_async(
    job_id: str, input_file: Path, config: JobSubmission
) -> dict[str, Any]:
    """Run pipeline steps asynchronously."""
    job = jobs[job_id]
    results = {}

    try:
        steps_to_run = config.steps
        if PipelineStep.ALL in steps_to_run:
            steps_to_run = [
                PipelineStep.EXTRACT,
                PipelineStep.LD_PROXIES,
                PipelineStep.LIFTOVER,
                PipelineStep.PREDICT,
                PipelineStep.SCORE,
            ]

        total_steps = len(steps_to_run)

        for i, step in enumerate(steps_to_run):
            job["current_step"] = step.value
            job["progress"] = (i / total_steps) * 100
            job["updated_at"] = datetime.now().isoformat()

            if step == PipelineStep.EXTRACT:
                # Load and extract significant variants
                gwas_df = pd.read_csv(input_file, sep="\t")
                significant = extract_variants.extract_significant_variants(
                    gwas_df, pval_threshold=config.pvalue_threshold
                )
                lead_snps = extract_variants.identify_lead_snps(significant)
                results["significant_count"] = len(significant)
                results["lead_snp_count"] = len(lead_snps)
                results["significant_variants"] = significant.to_dict(orient="records")

            elif step == PipelineStep.SCORE:
                # Score variants (using mock predictions for API demo)
                if "significant_variants" in results:
                    scorer = score_variants.VariantScorer(
                        {"scoring": {"consensus_method": config.consensus_method}}
                    )
                    # Generate mock predictions for demonstration
                    predictions = _generate_mock_predictions(
                        results["significant_variants"], config.tissues
                    )
                    tissue_scores = scorer.calculate_tissue_scores(predictions)
                    consensus = scorer.calculate_consensus_scores(tissue_scores)
                    ranked = scorer.rank_variants(consensus)
                    results["ranked_variants"] = ranked.head(100).to_dict(orient="records")
                    results["top_variant"] = ranked.iloc[0].to_dict() if len(ranked) > 0 else None

            # Simulate processing time for other steps
            await asyncio.sleep(0.5)

        # Save results
        results_file = RESULTS_DIR / f"{job_id}_results.json"
        pd.DataFrame(results.get("ranked_variants", [])).to_json(results_file, orient="records")

        job["status"] = JobStatus.COMPLETED
        job["progress"] = 100.0
        job["results_url"] = f"/api/v1/jobs/{job_id}/results"
        job["message"] = "Pipeline completed successfully"

    except Exception as e:
        job["status"] = JobStatus.FAILED
        job["message"] = str(e)
        results["error"] = str(e)

    job["updated_at"] = datetime.now().isoformat()
    results_cache[job_id] = results
    return results


def _generate_mock_predictions(variants: list[dict], tissues: list[str]) -> pd.DataFrame:
    """Generate mock AlphaGenome predictions for API demonstration."""
    import numpy as np

    predictions = []
    modalities = ["expression", "chromatin_accessibility"]

    for variant in variants:
        variant_id = f"chr{variant.get('chromosome', '1')}:{variant.get('position', 0)}"
        pval = variant.get("pvalue", 0.05)

        for tissue in tissues:
            for modality in modalities:
                base_effect = -np.log10(max(pval, 1e-300)) / 15
                predictions.append(
                    {
                        "variant_id": variant_id,
                        "rsid": variant.get("rsid", ""),
                        "tissue": tissue,
                        "modality": modality,
                        "effect_size": abs(np.random.normal(base_effect, 0.1)),
                        "confidence": np.random.uniform(0.6, 1.0),
                    }
                )

    return pd.DataFrame(predictions)


# API Endpoints
@app.get("/", response_model=HealthResponse)
async def root():
    """Root endpoint with health status."""
    return HealthResponse(
        status="healthy", version="1.0.0", timestamp=datetime.now().isoformat()
    )


@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint."""
    return HealthResponse(
        status="healthy", version="1.0.0", timestamp=datetime.now().isoformat()
    )


@app.post("/api/v1/jobs", response_model=JobResponse)
async def submit_job(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(..., description="GWAS summary statistics file (TSV/CSV)"),
    name: str = Query("analysis", description="Job name"),
    pvalue_threshold: float = Query(5e-8, description="P-value threshold"),
    consensus_method: str = Query("weighted_mean", description="Scoring method"),
):
    """
    Submit a new GWAS analysis job.

    Upload GWAS summary statistics and configure analysis parameters.
    The job will be processed asynchronously.
    """
    # Validate file
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file provided")

    content = await file.read()
    if len(content) == 0:
        raise HTTPException(status_code=400, detail="Empty file")

    # Save uploaded file
    job_id = generate_job_id()
    file_hash = get_file_hash(content)
    input_file = UPLOAD_DIR / f"{job_id}_{file_hash}.tsv"

    with open(input_file, "wb") as f:
        f.write(content)

    # Create job record
    now = datetime.now().isoformat()
    config = JobSubmission(
        name=name, pvalue_threshold=pvalue_threshold, consensus_method=consensus_method
    )

    job = {
        "job_id": job_id,
        "name": name,
        "status": JobStatus.PENDING,
        "created_at": now,
        "updated_at": now,
        "progress": 0.0,
        "current_step": None,
        "message": "Job submitted",
        "results_url": None,
        "config": config.model_dump(),
        "input_file": str(input_file),
    }
    jobs[job_id] = job

    # Start background processing
    job["status"] = JobStatus.RUNNING
    background_tasks.add_task(run_pipeline_async, job_id, input_file, config)

    return JobResponse(**job)


@app.get("/api/v1/jobs", response_model=list[JobResponse])
async def list_jobs(
    status: Optional[JobStatus] = Query(None, description="Filter by status"),
    limit: int = Query(50, ge=1, le=100, description="Maximum results"),
):
    """List all analysis jobs."""
    job_list = list(jobs.values())

    if status:
        job_list = [j for j in job_list if j["status"] == status]

    # Sort by created_at descending
    job_list.sort(key=lambda x: x["created_at"], reverse=True)

    return [JobResponse(**j) for j in job_list[:limit]]


@app.get("/api/v1/jobs/{job_id}", response_model=JobResponse)
async def get_job(job_id: str):
    """Get job status and details."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return JobResponse(**jobs[job_id])


@app.get("/api/v1/jobs/{job_id}/results")
async def get_job_results(job_id: str):
    """Get job results."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    job = jobs[job_id]
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(
            status_code=400, detail=f"Job not completed. Status: {job['status']}"
        )

    if job_id not in results_cache:
        raise HTTPException(status_code=404, detail="Results not found")

    return JSONResponse(content=results_cache[job_id])


@app.delete("/api/v1/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job and its results."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Clean up files
    job = jobs[job_id]
    input_file = Path(job.get("input_file", ""))
    if input_file.exists():
        input_file.unlink()

    results_file = RESULTS_DIR / f"{job_id}_results.json"
    if results_file.exists():
        results_file.unlink()

    # Remove from memory
    del jobs[job_id]
    if job_id in results_cache:
        del results_cache[job_id]

    return {"message": f"Job {job_id} deleted"}


@app.post("/api/v1/variants/score", response_model=list[VariantScore])
async def score_variants_endpoint(
    variants: list[VariantQuery],
    tissues: list[str] = Query(
        default=["Whole_Blood", "Liver"], description="Tissues for scoring"
    ),
):
    """
    Score a list of variants.

    Provides quick scoring for individual variants without full pipeline execution.
    """
    import numpy as np

    results = []
    for i, v in enumerate(variants):
        variant_id = f"chr{v.chromosome}:{v.position}"

        # Generate mock scores (replace with real AlphaGenome in production)
        tissue_effects = {t: abs(np.random.normal(0.5, 0.2)) for t in tissues}
        consensus = np.mean(list(tissue_effects.values()))

        results.append(
            VariantScore(
                variant_id=variant_id,
                chromosome=v.chromosome,
                position=v.position,
                rsid=None,
                consensus_score=round(consensus, 4),
                rank=i + 1,
                top_tissues=[
                    {"tissue": t, "effect": round(e, 4)} for t, e in tissue_effects.items()
                ],
                effect_sizes=tissue_effects,
            )
        )

    # Sort by score
    results.sort(key=lambda x: x.consensus_score, reverse=True)
    for i, r in enumerate(results):
        r.rank = i + 1

    return results


@app.get("/api/v1/tissues")
async def list_tissues():
    """List available tissues for AlphaGenome predictions."""
    from .alphagenome_predict import TISSUE_ONTOLOGY

    return {
        "tissues": [
            {"name": name, "uberon_id": uberon} for name, uberon in TISSUE_ONTOLOGY.items()
        ],
        "count": len(TISSUE_ONTOLOGY),
    }


@app.get("/api/v1/jobs/{job_id}/download")
async def download_results(job_id: str, format: str = Query("tsv", enum=["tsv", "json"])):
    """Download job results as a file."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job_id not in results_cache:
        raise HTTPException(status_code=404, detail="Results not found")

    results = results_cache[job_id]
    ranked = results.get("ranked_variants", [])

    if not ranked:
        raise HTTPException(status_code=404, detail="No ranked variants in results")

    df = pd.DataFrame(ranked)

    # Create temp file
    suffix = ".tsv" if format == "tsv" else ".json"
    output_file = RESULTS_DIR / f"{job_id}_download{suffix}"

    if format == "tsv":
        df.to_csv(output_file, sep="\t", index=False)
    else:
        df.to_json(output_file, orient="records", indent=2)

    return FileResponse(
        output_file,
        filename=f"alphagwas_{job_id}_results{suffix}",
        media_type="application/octet-stream",
    )


# WebSocket for real-time updates (optional)
@app.websocket("/ws/jobs/{job_id}")
async def websocket_job_updates(websocket, job_id: str):
    """WebSocket endpoint for real-time job updates."""
    from fastapi import WebSocket

    await websocket.accept()

    if job_id not in jobs:
        await websocket.send_json({"error": f"Job {job_id} not found"})
        await websocket.close()
        return

    # Poll for updates
    last_status = None
    while True:
        job = jobs.get(job_id)
        if not job:
            break

        if job["status"] != last_status or job["status"] == JobStatus.RUNNING:
            await websocket.send_json(
                {
                    "job_id": job_id,
                    "status": job["status"],
                    "progress": job["progress"],
                    "current_step": job["current_step"],
                    "message": job["message"],
                }
            )
            last_status = job["status"]

        if job["status"] in [JobStatus.COMPLETED, JobStatus.FAILED]:
            break

        await asyncio.sleep(1)

    await websocket.close()


def create_app() -> FastAPI:
    """Factory function for creating the FastAPI app."""
    return app


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
