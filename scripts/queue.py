"""
Job Queue System for AlphaGWAS.

Provides scalable async processing using Celery with Redis.

Features:
- Task queuing and execution
- Progress tracking
- Result storage
- Retry handling
- Priority queues
"""

import json
import logging
import os
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Optional

from pydantic import BaseModel

logger = logging.getLogger(__name__)

# Try to import Celery
try:
    from celery import Celery, Task
    from celery.result import AsyncResult
    HAS_CELERY = True
except ImportError:
    HAS_CELERY = False
    logger.warning("Celery not installed. Install with: pip install celery[redis]")
    Celery = None
    Task = object
    AsyncResult = None


class JobStatus(str, Enum):
    """Job status enumeration."""
    PENDING = "pending"
    STARTED = "started"
    RUNNING = "running"
    SUCCESS = "success"
    FAILURE = "failure"
    REVOKED = "revoked"


class JobInfo(BaseModel):
    """Job information model."""
    job_id: str
    task_name: str
    status: JobStatus
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    progress: float = 0.0
    message: str = ""
    result: Optional[Any] = None
    error: Optional[str] = None


class QueueConfig(BaseModel):
    """Queue configuration."""
    broker_url: str = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
    result_backend: str = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")
    task_serializer: str = "json"
    result_serializer: str = "json"
    accept_content: list[str] = ["json"]
    timezone: str = "UTC"
    task_track_started: bool = True
    task_time_limit: int = 3600  # 1 hour
    task_soft_time_limit: int = 3300  # 55 minutes
    worker_prefetch_multiplier: int = 1
    worker_concurrency: int = 4


def create_celery_app(config: Optional[QueueConfig] = None) -> "Celery":
    """
    Create and configure Celery application.

    Args:
        config: Queue configuration

    Returns:
        Configured Celery app
    """
    if not HAS_CELERY:
        raise RuntimeError("Celery not installed. Install with: pip install celery[redis]")

    config = config or QueueConfig()

    app = Celery(
        "alphagwas",
        broker=config.broker_url,
        backend=config.result_backend,
    )

    app.conf.update(
        task_serializer=config.task_serializer,
        result_serializer=config.result_serializer,
        accept_content=config.accept_content,
        timezone=config.timezone,
        task_track_started=config.task_track_started,
        task_time_limit=config.task_time_limit,
        task_soft_time_limit=config.task_soft_time_limit,
        worker_prefetch_multiplier=config.worker_prefetch_multiplier,
        worker_concurrency=config.worker_concurrency,
        # Task routes for priority queues
        task_routes={
            "alphagwas.tasks.high_priority.*": {"queue": "high"},
            "alphagwas.tasks.low_priority.*": {"queue": "low"},
        },
        # Default queue
        task_default_queue="default",
    )

    return app


# Create default app instance
celery_app = create_celery_app() if HAS_CELERY else None


class ProgressTask(Task if HAS_CELERY else object):
    """
    Base task with progress tracking support.

    Usage:
        @celery_app.task(base=ProgressTask, bind=True)
        def my_task(self, arg1, arg2):
            self.update_progress(0.5, "Halfway done")
            ...
    """

    def update_progress(self, progress: float, message: str = ""):
        """Update task progress."""
        if hasattr(self, 'request'):
            self.update_state(
                state="PROGRESS",
                meta={
                    "progress": progress,
                    "message": message,
                    "updated_at": datetime.utcnow().isoformat(),
                }
            )


# ==================== Pipeline Tasks ====================

if HAS_CELERY:

    @celery_app.task(base=ProgressTask, bind=True, name="alphagwas.extract_variants")
    def task_extract_variants(
        self,
        input_file: str,
        output_dir: str,
        pval_threshold: float = 5e-8,
    ) -> dict:
        """Extract significant variants from GWAS data."""
        import pandas as pd
        from scripts import extract_variants

        self.update_progress(0.1, "Loading GWAS data")

        gwas_df = pd.read_csv(input_file, sep="\t")
        self.update_progress(0.3, "Extracting significant variants")

        significant = extract_variants.extract_significant_variants(
            gwas_df, pval_threshold=pval_threshold
        )
        self.update_progress(0.6, "Identifying lead SNPs")

        lead_snps = extract_variants.identify_lead_snps(significant)
        self.update_progress(0.8, "Saving results")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        sig_file = output_path / "significant_variants.tsv"
        lead_file = output_path / "lead_snps.tsv"

        significant.to_csv(sig_file, sep="\t", index=False)
        lead_snps.to_csv(lead_file, sep="\t", index=False)

        self.update_progress(1.0, "Complete")

        return {
            "n_significant": len(significant),
            "n_lead_snps": len(lead_snps),
            "output_files": [str(sig_file), str(lead_file)],
        }

    @celery_app.task(base=ProgressTask, bind=True, name="alphagwas.predict_variants")
    def task_predict_variants(
        self,
        variants_file: str,
        output_file: str,
        tissues: list[str] = None,
        batch_size: int = 100,
    ) -> dict:
        """Run AlphaGenome predictions on variants."""
        import pandas as pd
        from scripts import alphagenome_predict

        self.update_progress(0.1, "Loading variants")

        variants_df = pd.read_csv(variants_file, sep="\t")
        n_variants = len(variants_df)

        self.update_progress(0.2, f"Initializing predictor for {n_variants} variants")

        predictor = alphagenome_predict.AlphaGenomePredictor(config={})

        # Predict with progress updates
        predictions = []
        for i in range(0, n_variants, batch_size):
            batch = variants_df.iloc[i:i + batch_size]
            batch_preds = predictor.predict_batch(batch.to_dict("records"))
            predictions.extend(batch_preds)

            progress = 0.2 + 0.7 * (i + len(batch)) / n_variants
            self.update_progress(progress, f"Predicted {i + len(batch)}/{n_variants} variants")

        self.update_progress(0.95, "Saving predictions")

        pred_df = pd.DataFrame([p.__dict__ for p in predictions])
        pred_df.to_parquet(output_file, index=False)

        self.update_progress(1.0, "Complete")

        return {
            "n_predictions": len(predictions),
            "output_file": output_file,
        }

    @celery_app.task(base=ProgressTask, bind=True, name="alphagwas.score_variants")
    def task_score_variants(
        self,
        predictions_file: str,
        output_dir: str,
        consensus_method: str = "weighted_mean",
    ) -> dict:
        """Score and rank variants."""
        import pandas as pd
        from scripts import score_variants

        self.update_progress(0.1, "Loading predictions")

        predictions_df = pd.read_parquet(predictions_file)

        self.update_progress(0.3, "Calculating tissue scores")

        scorer = score_variants.VariantScorer(
            {"scoring": {"consensus_method": consensus_method}}
        )

        tissue_scores = scorer.calculate_tissue_scores(predictions_df)
        self.update_progress(0.5, "Calculating consensus scores")

        consensus = scorer.calculate_consensus_scores(tissue_scores)
        self.update_progress(0.7, "Ranking variants")

        ranked = scorer.rank_variants(consensus)

        self.update_progress(0.9, "Saving results")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        ranked.to_csv(output_path / "ranked_variants.tsv", sep="\t", index=False)
        tissue_scores.to_csv(output_path / "tissue_scores.tsv", sep="\t", index=False)

        self.update_progress(1.0, "Complete")

        return {
            "n_ranked": len(ranked),
            "top_variant": ranked.iloc[0].to_dict() if len(ranked) > 0 else None,
            "output_dir": str(output_path),
        }

    @celery_app.task(base=ProgressTask, bind=True, name="alphagwas.run_pipeline")
    def task_run_full_pipeline(
        self,
        input_file: str,
        output_dir: str,
        config: dict = None,
    ) -> dict:
        """Run full AlphaGWAS pipeline."""
        self.update_progress(0.0, "Starting pipeline")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract variants
        self.update_progress(0.1, "Step 1: Extracting variants")
        extract_result = task_extract_variants(
            input_file,
            str(output_path / "step1_extract"),
            pval_threshold=config.get("pval_threshold", 5e-8) if config else 5e-8,
        )

        # Step 2: Predict (using significant variants)
        self.update_progress(0.4, "Step 2: Running predictions")
        predict_result = task_predict_variants(
            extract_result["output_files"][0],
            str(output_path / "step2_predictions.parquet"),
        )

        # Step 3: Score
        self.update_progress(0.7, "Step 3: Scoring variants")
        score_result = task_score_variants(
            predict_result["output_file"],
            str(output_path / "step3_scores"),
        )

        self.update_progress(1.0, "Pipeline complete")

        return {
            "extract": extract_result,
            "predict": predict_result,
            "score": score_result,
            "output_dir": str(output_path),
        }


# ==================== Job Manager ====================

class JobManager:
    """
    Manage async jobs and retrieve results.

    Provides a simple interface for submitting tasks and checking status.
    """

    def __init__(self, celery_app: "Celery" = None):
        """Initialize job manager."""
        self.app = celery_app or globals().get("celery_app")
        if not self.app:
            raise RuntimeError("Celery app not available")

    def submit_task(
        self,
        task_name: str,
        args: tuple = (),
        kwargs: dict = None,
        queue: str = "default",
        priority: int = 5,
    ) -> str:
        """
        Submit a task for async execution.

        Args:
            task_name: Name of the task (e.g., "alphagwas.extract_variants")
            args: Positional arguments
            kwargs: Keyword arguments
            queue: Queue name (default, high, low)
            priority: Task priority (0-9, higher = more important)

        Returns:
            Job ID
        """
        task = self.app.send_task(
            task_name,
            args=args,
            kwargs=kwargs or {},
            queue=queue,
            priority=priority,
        )
        logger.info(f"Submitted task {task_name} with ID {task.id}")
        return task.id

    def get_job_status(self, job_id: str) -> JobInfo:
        """Get job status and info."""
        result = AsyncResult(job_id, app=self.app)

        status_map = {
            "PENDING": JobStatus.PENDING,
            "STARTED": JobStatus.STARTED,
            "PROGRESS": JobStatus.RUNNING,
            "SUCCESS": JobStatus.SUCCESS,
            "FAILURE": JobStatus.FAILURE,
            "REVOKED": JobStatus.REVOKED,
        }

        status = status_map.get(result.status, JobStatus.PENDING)
        info = result.info or {}

        job_info = JobInfo(
            job_id=job_id,
            task_name=result.name or "unknown",
            status=status,
            created_at=datetime.utcnow(),  # Celery doesn't track this by default
            progress=info.get("progress", 0.0) if isinstance(info, dict) else 0.0,
            message=info.get("message", "") if isinstance(info, dict) else "",
        )

        if status == JobStatus.SUCCESS:
            job_info.result = result.result
            job_info.progress = 1.0
        elif status == JobStatus.FAILURE:
            job_info.error = str(result.result) if result.result else "Unknown error"

        return job_info

    def get_job_result(self, job_id: str, timeout: float = None) -> Any:
        """
        Get job result (blocking).

        Args:
            job_id: Job ID
            timeout: Max wait time in seconds

        Returns:
            Task result

        Raises:
            TimeoutError: If timeout exceeded
            Exception: If task failed
        """
        result = AsyncResult(job_id, app=self.app)
        return result.get(timeout=timeout)

    def cancel_job(self, job_id: str) -> bool:
        """Cancel a pending or running job."""
        self.app.control.revoke(job_id, terminate=True)
        logger.info(f"Revoked job {job_id}")
        return True

    def list_active_jobs(self) -> list[str]:
        """List active job IDs."""
        inspect = self.app.control.inspect()
        active = inspect.active() or {}

        job_ids = []
        for worker_tasks in active.values():
            for task in worker_tasks:
                job_ids.append(task["id"])

        return job_ids

    def list_scheduled_jobs(self) -> list[str]:
        """List scheduled job IDs."""
        inspect = self.app.control.inspect()
        scheduled = inspect.scheduled() or {}

        job_ids = []
        for worker_tasks in scheduled.values():
            for task in worker_tasks:
                job_ids.append(task["request"]["id"])

        return job_ids

    def get_queue_stats(self) -> dict:
        """Get queue statistics."""
        inspect = self.app.control.inspect()

        return {
            "active": len(self.list_active_jobs()),
            "scheduled": len(self.list_scheduled_jobs()),
            "stats": inspect.stats() or {},
        }


# ==================== FastAPI Integration ====================

def create_queue_routes(job_manager: JobManager):
    """Create FastAPI routes for job management."""
    from fastapi import APIRouter, HTTPException

    router = APIRouter(prefix="/jobs", tags=["Jobs"])

    class SubmitJobRequest(BaseModel):
        task_name: str
        args: list = []
        kwargs: dict = {}
        queue: str = "default"
        priority: int = 5

    @router.post("/submit")
    async def submit_job(request: SubmitJobRequest):
        """Submit a new job."""
        try:
            job_id = job_manager.submit_task(
                request.task_name,
                args=tuple(request.args),
                kwargs=request.kwargs,
                queue=request.queue,
                priority=request.priority,
            )
            return {"job_id": job_id, "message": "Job submitted"}
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    @router.get("/{job_id}")
    async def get_job(job_id: str):
        """Get job status."""
        try:
            info = job_manager.get_job_status(job_id)
            return info.model_dump()
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    @router.delete("/{job_id}")
    async def cancel_job(job_id: str):
        """Cancel a job."""
        try:
            job_manager.cancel_job(job_id)
            return {"message": "Job cancelled"}
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    @router.get("/")
    async def list_jobs():
        """List all active and scheduled jobs."""
        try:
            return {
                "active": job_manager.list_active_jobs(),
                "scheduled": job_manager.list_scheduled_jobs(),
            }
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    @router.get("/stats")
    async def get_stats():
        """Get queue statistics."""
        try:
            return job_manager.get_queue_stats()
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    return router
