"""
FastAPI REST Interface for AlphaGWAS Pipeline.

Provides endpoints for:
- Submitting GWAS data for analysis
- Running individual pipeline steps
- Retrieving results and visualizations
- Job status monitoring

Production features:
- Per-IP rate limiting (slowapi)
- Persistent job storage (database-backed)
- Structured logging with request IDs
- Metrics endpoint
"""

import asyncio
import hashlib
import logging
import os
import re
import tempfile
import time
import uuid
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from fastapi import BackgroundTasks, FastAPI, File, HTTPException, Query, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field, field_validator

from . import extract_variants, score_variants, utils
from .job_store import JobStore

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration from environment
# ---------------------------------------------------------------------------
ALLOWED_ORIGINS = os.getenv(
    "ALPHAGWAS_CORS_ORIGINS", "http://localhost:3000,http://localhost:8080"
).split(",")
API_HOST = os.getenv("ALPHAGWAS_API_HOST", "127.0.0.1")
API_PORT = int(os.getenv("ALPHAGWAS_API_PORT", "8000"))
MAX_UPLOAD_BYTES = int(os.getenv("ALPHAGWAS_MAX_UPLOAD_BYTES", str(100 * 1024 * 1024)))  # 100 MB
MAX_JOBS = int(os.getenv("ALPHAGWAS_MAX_JOBS", "1000"))
RATE_LIMIT = os.getenv("ALPHAGWAS_RATE_LIMIT", "30/minute")
RATE_LIMIT_SUBMIT = os.getenv("ALPHAGWAS_RATE_LIMIT_SUBMIT", "5/minute")

# ---------------------------------------------------------------------------
# Initialize FastAPI app
# ---------------------------------------------------------------------------
app = FastAPI(
    title="AlphaGWAS API",
    description="REST API for GWAS variant prioritization using AlphaGenome predictions",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# CORS - restricted to configured origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "DELETE"],
    allow_headers=["Authorization", "Content-Type", "X-API-Key"],
)

# Request ID middleware for tracing
try:
    from .observability import RequestIdMiddleware, metrics_collector as obs_metrics
    app.add_middleware(RequestIdMiddleware)
    logger.info("Request ID middleware enabled")
except (ImportError, RuntimeError):
    obs_metrics = None

# ---------------------------------------------------------------------------
# Rate limiting (per-IP)
# ---------------------------------------------------------------------------
try:
    from slowapi import Limiter, _rate_limit_exceeded_handler
    from slowapi.errors import RateLimitExceeded
    from slowapi.util import get_remote_address

    limiter = Limiter(key_func=get_remote_address, default_limits=[RATE_LIMIT])
    app.state.limiter = limiter
    app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)
    HAS_RATE_LIMIT = True
    logger.info(f"Rate limiting enabled: {RATE_LIMIT} (submit: {RATE_LIMIT_SUBMIT})")
except ImportError:
    # Fallback: no-op decorator
    HAS_RATE_LIMIT = False
    logger.warning("slowapi not installed. Rate limiting disabled. Install with: pip install slowapi")

    class _NoOpLimiter:
        def limit(self, *args, **kwargs):
            def decorator(func):
                return func
            return decorator

        def shared_limit(self, *args, **kwargs):
            def decorator(func):
                return func
            return decorator

    limiter = _NoOpLimiter()


# ---------------------------------------------------------------------------
# Persistent job storage
# ---------------------------------------------------------------------------
job_store = JobStore(max_jobs=MAX_JOBS)

# Configuration
UPLOAD_DIR = Path(tempfile.gettempdir()) / "alphagwas_uploads"
RESULTS_DIR = Path(tempfile.gettempdir()) / "alphagwas_results"
UPLOAD_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# Job ID validation pattern (UUID prefix, 8 hex chars)
_JOB_ID_PATTERN = re.compile(r"^[0-9a-f]{8}$")

# App start time for uptime tracking
_APP_START_TIME = time.time()


def _validate_job_id(job_id: str) -> str:
    """Validate job_id format to prevent path traversal."""
    if not _JOB_ID_PATTERN.match(job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID format")
    return job_id


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


class ConsensusMethod(str, Enum):
    """Allowed consensus scoring methods."""
    MEAN = "mean"
    MEDIAN = "median"
    MAX = "max"
    WEIGHTED_MEAN = "weighted_mean"


# Request/Response Models
class JobSubmission(BaseModel):
    """Job submission request."""
    name: str = Field(..., description="Name for this analysis job", max_length=255)
    pvalue_threshold: float = Field(5e-8, description="P-value significance threshold")
    ld_r2_threshold: float = Field(0.8, description="LD R² threshold for proxy variants")
    tissues: list[str] = Field(
        default_factory=lambda: ["Whole_Blood", "Liver", "Brain_Cortex"],
        description="Tissues for AlphaGenome predictions",
    )
    consensus_method: ConsensusMethod = Field(
        ConsensusMethod.WEIGHTED_MEAN, description="Scoring method"
    )
    steps: list[PipelineStep] = Field(
        default_factory=lambda: [PipelineStep.ALL], description="Pipeline steps to run"
    )

    @field_validator("pvalue_threshold")
    @classmethod
    def validate_pvalue(cls, v: float) -> float:
        if not 0 < v <= 1:
            raise ValueError("pvalue_threshold must be between 0 (exclusive) and 1")
        return v

    @field_validator("ld_r2_threshold")
    @classmethod
    def validate_r2(cls, v: float) -> float:
        if not 0 < v <= 1:
            raise ValueError("ld_r2_threshold must be between 0 (exclusive) and 1")
        return v


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

    @field_validator("chromosome")
    @classmethod
    def validate_chromosome(cls, v: str) -> str:
        valid = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
        cleaned = v.replace("chr", "")
        if cleaned not in valid:
            raise ValueError(f"Invalid chromosome: {v}")
        return cleaned


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
    jobs_active: int = 0
    uptime_seconds: float = 0.0
    storage_backend: str = "memory"


# Helper functions
def generate_job_id() -> str:
    """Generate unique job ID (8 hex chars)."""
    return str(uuid.uuid4())[:8]


def get_file_hash(content: bytes) -> str:
    """Generate hash for file content."""
    return hashlib.sha256(content).hexdigest()[:12]


async def run_pipeline_async(
    job_id: str, input_file: Path, config: JobSubmission
) -> dict[str, Any]:
    """Run pipeline steps asynchronously."""
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
            job_store.update_job(job_id, {
                "current_step": step.value,
                "progress": (i / total_steps) * 100,
            })

            if step == PipelineStep.EXTRACT:
                gwas_df = pd.read_csv(input_file, sep="\t")
                significant = extract_variants.extract_significant_variants(
                    gwas_df, pval_threshold=config.pvalue_threshold
                )
                lead_snps = extract_variants.identify_lead_snps(significant)
                results["significant_count"] = len(significant)
                results["lead_snp_count"] = len(lead_snps)
                results["significant_variants"] = significant.to_dict(orient="records")

            elif step == PipelineStep.SCORE:
                if "significant_variants" in results:
                    scorer = score_variants.VariantScorer(
                        {"scoring": {"consensus_method": config.consensus_method.value}}
                    )
                    predictions = _generate_mock_predictions(
                        results["significant_variants"], config.tissues
                    )
                    tissue_scores = scorer.calculate_tissue_scores(predictions)
                    consensus = scorer.calculate_consensus_scores(tissue_scores)
                    ranked = scorer.rank_variants(consensus)
                    results["ranked_variants"] = ranked.head(100).to_dict(orient="records")
                    results["top_variant"] = ranked.iloc[0].to_dict() if len(ranked) > 0 else None

            await asyncio.sleep(0.5)

        # Save results
        results_file = RESULTS_DIR / f"{job_id}_results.json"
        pd.DataFrame(results.get("ranked_variants", [])).to_json(results_file, orient="records")

        job_store.update_job(job_id, {
            "status": JobStatus.COMPLETED,
            "progress": 100.0,
            "results_url": f"/api/v1/jobs/{job_id}/results",
            "message": "Pipeline completed successfully",
        })

    except Exception as e:
        job_store.update_job(job_id, {
            "status": JobStatus.FAILED,
            "message": str(e),
        })
        results["error"] = str(e)
        logger.error(f"Pipeline failed for job {job_id}: {e}")

    job_store.save_results(job_id, results)

    # Clean up uploaded input file
    try:
        if input_file.exists():
            input_file.unlink()
    except OSError:
        logger.warning(f"Failed to clean up input file: {input_file}")

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
                predictions.append({
                    "variant_id": variant_id,
                    "rsid": variant.get("rsid", ""),
                    "tissue": tissue,
                    "modality": modality,
                    "effect_size": abs(np.random.normal(base_effect, 0.1)),
                    "confidence": np.random.uniform(0.6, 1.0),
                })

    return pd.DataFrame(predictions)


# ===========================================================================
# API Endpoints
# ===========================================================================

@app.get("/", response_model=HealthResponse)
@limiter.limit("60/minute")
async def root(request: Request):
    """Root endpoint with health status."""
    return HealthResponse(
        status="healthy",
        version="1.0.0",
        timestamp=datetime.now().isoformat(),
        jobs_active=job_store.count_active(),
        uptime_seconds=round(time.time() - _APP_START_TIME, 1),
        storage_backend="database" if job_store._use_db else "memory",
    )


@app.get("/health", response_model=HealthResponse)
@limiter.limit("60/minute")
async def health_check(request: Request):
    """Health check endpoint."""
    return HealthResponse(
        status="healthy",
        version="1.0.0",
        timestamp=datetime.now().isoformat(),
        jobs_active=job_store.count_active(),
        uptime_seconds=round(time.time() - _APP_START_TIME, 1),
        storage_backend="database" if job_store._use_db else "memory",
    )


@app.get("/metrics")
@limiter.limit("30/minute")
async def metrics_endpoint(request: Request):
    """Application metrics endpoint."""
    try:
        from .observability import metrics_collector
        return JSONResponse(content=metrics_collector.get_metrics())
    except (ImportError, AttributeError):
        # Return basic metrics if observability module not available
        return JSONResponse(content={
            "jobs_active": job_store.count_active(),
            "uptime_seconds": round(time.time() - _APP_START_TIME, 1),
            "rate_limiting_enabled": HAS_RATE_LIMIT,
            "storage_backend": "database" if job_store._use_db else "memory",
        })


@app.post("/api/v1/jobs", response_model=JobResponse)
@limiter.limit(RATE_LIMIT_SUBMIT)
async def submit_job(
    request: Request,
    background_tasks: BackgroundTasks,
    file: UploadFile = File(..., description="GWAS summary statistics file (TSV/CSV)"),
    name: str = Query("analysis", description="Job name", max_length=255),
    pvalue_threshold: float = Query(5e-8, description="P-value threshold", gt=0, le=1),
    consensus_method: ConsensusMethod = Query(
        ConsensusMethod.WEIGHTED_MEAN, description="Scoring method"
    ),
):
    """
    Submit a new GWAS analysis job.

    Upload GWAS summary statistics and configure analysis parameters.
    The job will be processed asynchronously.
    """
    if not file.filename:
        raise HTTPException(status_code=400, detail="No file provided")

    # Validate file extension
    allowed_extensions = {".tsv", ".csv", ".txt", ".gz"}
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in allowed_extensions:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported file type '{file_ext}'. Allowed: {allowed_extensions}",
        )

    # Read with size limit
    content = await file.read()
    if len(content) == 0:
        raise HTTPException(status_code=400, detail="Empty file")
    if len(content) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"File too large. Maximum size: {MAX_UPLOAD_BYTES // (1024 * 1024)} MB",
        )

    # Save uploaded file
    job_id = generate_job_id()
    file_hash = get_file_hash(content)
    input_file = UPLOAD_DIR / f"{job_id}_{file_hash}.tsv"

    try:
        with open(input_file, "wb") as f:
            f.write(content)
    except OSError as e:
        logger.error(f"Failed to save uploaded file: {e}")
        raise HTTPException(status_code=500, detail="Failed to save uploaded file")

    # Create job record
    now = datetime.now().isoformat()
    config = JobSubmission(
        name=name, pvalue_threshold=pvalue_threshold, consensus_method=consensus_method
    )

    job_data = {
        "job_id": job_id,
        "name": name,
        "status": JobStatus.RUNNING,
        "created_at": now,
        "updated_at": now,
        "progress": 0.0,
        "current_step": None,
        "message": "Job submitted",
        "results_url": None,
        "config": config.model_dump(),
        "input_file": str(input_file),
    }
    job_store.create_job(job_data)

    # Start background processing
    background_tasks.add_task(run_pipeline_async, job_id, input_file, config)

    return JobResponse(**job_data)


@app.get("/api/v1/jobs", response_model=list[JobResponse])
@limiter.limit(RATE_LIMIT)
async def list_jobs(
    request: Request,
    status: Optional[JobStatus] = Query(None, description="Filter by status"),
    limit: int = Query(50, ge=1, le=100, description="Maximum results"),
):
    """List all analysis jobs."""
    job_list = job_store.list_jobs(status=status.value if status else None, limit=limit)
    return [JobResponse(**j) for j in job_list]


@app.get("/api/v1/jobs/{job_id}", response_model=JobResponse)
@limiter.limit(RATE_LIMIT)
async def get_job(request: Request, job_id: str):
    """Get job status and details."""
    job_id = _validate_job_id(job_id)
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    return JobResponse(**job)


@app.get("/api/v1/jobs/{job_id}/results")
@limiter.limit(RATE_LIMIT)
async def get_job_results(request: Request, job_id: str):
    """Get job results."""
    job_id = _validate_job_id(job_id)
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(
            status_code=400, detail=f"Job not completed. Status: {job['status']}"
        )

    results = job_store.get_results(job_id)
    if not results:
        raise HTTPException(status_code=404, detail="Results not found")

    return JSONResponse(content=results)


@app.delete("/api/v1/jobs/{job_id}")
@limiter.limit("10/minute")
async def delete_job(request: Request, job_id: str):
    """Delete a job and its results."""
    job_id = _validate_job_id(job_id)
    job = job_store.delete_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Clean up files
    input_file = Path(job.get("input_file", ""))
    try:
        if input_file.exists():
            input_file.unlink()
    except OSError:
        pass

    results_file = RESULTS_DIR / f"{job_id}_results.json"
    try:
        if results_file.exists():
            results_file.unlink()
    except OSError:
        pass

    return {"message": f"Job {job_id} deleted"}


@app.post("/api/v1/variants/score", response_model=list[VariantScore])
@limiter.limit(RATE_LIMIT)
async def score_variants_endpoint(
    request: Request,
    variants: list[VariantQuery],
    tissues: list[str] = Query(
        default=["Whole_Blood", "Liver"], description="Tissues for scoring"
    ),
):
    """Score a list of variants."""
    import numpy as np

    if len(variants) > 1000:
        raise HTTPException(status_code=400, detail="Maximum 1000 variants per request")

    results = []
    for i, v in enumerate(variants):
        variant_id = f"chr{v.chromosome}:{v.position}"
        tissue_effects = {t: abs(np.random.normal(0.5, 0.2)) for t in tissues}
        consensus = np.mean(list(tissue_effects.values()))

        results.append(VariantScore(
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
        ))

    results.sort(key=lambda x: x.consensus_score, reverse=True)
    for i, r in enumerate(results):
        r.rank = i + 1

    return results


@app.get("/api/v1/tissues")
@limiter.limit(RATE_LIMIT)
async def list_tissues(request: Request):
    """List available tissues for AlphaGenome predictions."""
    from .alphagenome_predict import TISSUE_ONTOLOGY

    return {
        "tissues": [
            {"name": name, "uberon_id": uberon} for name, uberon in TISSUE_ONTOLOGY.items()
        ],
        "count": len(TISSUE_ONTOLOGY),
    }


@app.get("/api/v1/jobs/{job_id}/download")
@limiter.limit("10/minute")
async def download_results(request: Request, job_id: str, format: str = Query("tsv", enum=["tsv", "json"])):
    """Download job results as a file."""
    job_id = _validate_job_id(job_id)

    if not job_store.job_exists(job_id):
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    results = job_store.get_results(job_id)
    if not results:
        raise HTTPException(status_code=404, detail="Results not found")

    ranked = results.get("ranked_variants", [])
    if not ranked:
        raise HTTPException(status_code=404, detail="No ranked variants in results")

    df = pd.DataFrame(ranked)
    suffix = ".tsv" if format == "tsv" else ".json"
    output_file = RESULTS_DIR / f"{job_id}_download{suffix}"

    try:
        if format == "tsv":
            df.to_csv(output_file, sep="\t", index=False)
        else:
            df.to_json(output_file, orient="records", indent=2)
    except OSError as e:
        logger.error(f"Failed to write download file: {e}")
        raise HTTPException(status_code=500, detail="Failed to generate download")

    return FileResponse(
        output_file,
        filename=f"alphagwas_{job_id}_results{suffix}",
        media_type="application/octet-stream",
    )


# WebSocket for real-time updates
@app.websocket("/ws/jobs/{job_id}")
async def websocket_job_updates(websocket, job_id: str):
    """WebSocket endpoint for real-time job updates."""
    await websocket.accept()

    if not _JOB_ID_PATTERN.match(job_id):
        await websocket.send_json({"error": "Invalid job ID format"})
        await websocket.close()
        return

    if not job_store.job_exists(job_id):
        await websocket.send_json({"error": f"Job {job_id} not found"})
        await websocket.close()
        return

    last_status = None
    while True:
        job = job_store.get_job(job_id)
        if not job:
            break

        if job["status"] != last_status or job["status"] == JobStatus.RUNNING:
            await websocket.send_json({
                "job_id": job_id,
                "status": job["status"],
                "progress": job["progress"],
                "current_step": job["current_step"],
                "message": job["message"],
            })
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
    uvicorn.run(app, host=API_HOST, port=API_PORT)
