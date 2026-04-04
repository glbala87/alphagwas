"""
Persistent Job Storage for AlphaGWAS API.

Provides a database-backed job store that survives application restarts.
Falls back to in-memory storage if database is unavailable.
"""

import json
import logging
import os
import threading
from datetime import datetime
from typing import Any, Optional

logger = logging.getLogger(__name__)

# Try to import SQLAlchemy
try:
    from sqlalchemy import Column, DateTime, Float, Integer, String, Text, create_engine
    from sqlalchemy.ext.declarative import declarative_base
    from sqlalchemy.orm import sessionmaker

    HAS_SQLALCHEMY = True
except ImportError:
    HAS_SQLALCHEMY = False

JobStoreBase = declarative_base() if HAS_SQLALCHEMY else None


if HAS_SQLALCHEMY:

    class JobRecord(JobStoreBase):
        """Persistent job record."""

        __tablename__ = "api_jobs"

        job_id = Column(String(8), primary_key=True)
        name = Column(String(255), nullable=False)
        status = Column(String(20), nullable=False, default="pending")
        created_at = Column(String(50), nullable=False)
        updated_at = Column(String(50), nullable=False)
        progress = Column(Float, default=0.0)
        current_step = Column(String(50), nullable=True)
        message = Column(Text, nullable=True)
        results_url = Column(String(255), nullable=True)
        config_json = Column(Text, default="{}")
        input_file = Column(String(500), nullable=True)
        results_json = Column(Text, nullable=True)

        def to_dict(self) -> dict:
            return {
                "job_id": self.job_id,
                "name": self.name,
                "status": self.status,
                "created_at": self.created_at,
                "updated_at": self.updated_at,
                "progress": self.progress or 0.0,
                "current_step": self.current_step,
                "message": self.message,
                "results_url": self.results_url,
                "config": json.loads(self.config_json) if self.config_json else {},
                "input_file": self.input_file,
            }

        def get_results(self) -> Optional[dict]:
            if self.results_json:
                try:
                    return json.loads(self.results_json)
                except (json.JSONDecodeError, TypeError):
                    return None
            return None


class JobStore:
    """
    Persistent job storage with database backend.

    Falls back to in-memory dict if no database is configured.
    Thread-safe for concurrent access.
    """

    def __init__(self, db_url: Optional[str] = None, max_jobs: int = 1000):
        self._lock = threading.Lock()
        self._max_jobs = max_jobs
        self._db_url = db_url or os.getenv("ALPHAGWAS_JOB_DB_URL") or os.getenv("DATABASE_URL")
        self._engine = None
        self._Session = None

        # In-memory fallback
        self._memory_jobs: dict[str, dict] = {}
        self._memory_results: dict[str, Any] = {}
        self._use_db = False

        if self._db_url and HAS_SQLALCHEMY:
            try:
                self._engine = create_engine(self._db_url, pool_pre_ping=True)
                JobStoreBase.metadata.create_all(self._engine)
                self._Session = sessionmaker(bind=self._engine)
                self._use_db = True
                logger.info(f"Job store using database: {self._safe_url()}")
            except Exception as e:
                logger.warning(f"Failed to connect job store DB, falling back to memory: {e}")
                self._use_db = False
        else:
            logger.info("Job store using in-memory storage (set DATABASE_URL for persistence)")

    def _safe_url(self) -> str:
        if self._db_url and "@" in self._db_url:
            parts = self._db_url.split("@")
            prefix = parts[0].rsplit(":", 1)[0]
            return f"{prefix}:****@{parts[1]}"
        return self._db_url or "memory"

    # ==================== CRUD Operations ====================

    def create_job(self, job_data: dict) -> None:
        """Create a new job record."""
        job_id = job_data["job_id"]

        if self._use_db:
            session = self._Session()
            try:
                record = JobRecord(
                    job_id=job_id,
                    name=job_data.get("name", ""),
                    status=job_data.get("status", "pending"),
                    created_at=job_data.get("created_at", datetime.now().isoformat()),
                    updated_at=job_data.get("updated_at", datetime.now().isoformat()),
                    progress=job_data.get("progress", 0.0),
                    current_step=job_data.get("current_step"),
                    message=job_data.get("message"),
                    results_url=job_data.get("results_url"),
                    config_json=json.dumps(job_data.get("config", {}), default=str),
                    input_file=job_data.get("input_file"),
                )
                session.add(record)
                session.commit()
            except Exception:
                session.rollback()
                raise
            finally:
                session.close()
        else:
            with self._lock:
                self._memory_jobs[job_id] = job_data.copy()
                self._evict_if_needed()

    def get_job(self, job_id: str) -> Optional[dict]:
        """Get job by ID."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                return record.to_dict() if record else None
            finally:
                session.close()
        else:
            with self._lock:
                job = self._memory_jobs.get(job_id)
                return job.copy() if job else None

    def update_job(self, job_id: str, updates: dict) -> bool:
        """Update job fields."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                if not record:
                    return False
                for key, value in updates.items():
                    if key == "config":
                        record.config_json = json.dumps(value, default=str)
                    elif hasattr(record, key):
                        setattr(record, key, value)
                record.updated_at = datetime.now().isoformat()
                session.commit()
                return True
            except Exception:
                session.rollback()
                raise
            finally:
                session.close()
        else:
            with self._lock:
                if job_id not in self._memory_jobs:
                    return False
                self._memory_jobs[job_id].update(updates)
                self._memory_jobs[job_id]["updated_at"] = datetime.now().isoformat()
                return True

    def delete_job(self, job_id: str) -> Optional[dict]:
        """Delete job and return its data."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                if not record:
                    return None
                data = record.to_dict()
                session.delete(record)
                session.commit()
                return data
            except Exception:
                session.rollback()
                raise
            finally:
                session.close()
        else:
            with self._lock:
                return self._memory_jobs.pop(job_id, None)

    def list_jobs(
        self,
        status: Optional[str] = None,
        limit: int = 50,
    ) -> list[dict]:
        """List jobs, optionally filtered by status."""
        if self._use_db:
            session = self._Session()
            try:
                query = session.query(JobRecord)
                if status:
                    query = query.filter_by(status=status)
                query = query.order_by(JobRecord.created_at.desc()).limit(limit)
                return [r.to_dict() for r in query.all()]
            finally:
                session.close()
        else:
            with self._lock:
                jobs = list(self._memory_jobs.values())
            if status:
                jobs = [j for j in jobs if j.get("status") == status]
            jobs.sort(key=lambda x: x.get("created_at", ""), reverse=True)
            return [j.copy() for j in jobs[:limit]]

    def count_active(self) -> int:
        """Count running jobs."""
        if self._use_db:
            session = self._Session()
            try:
                return session.query(JobRecord).filter_by(status="running").count()
            finally:
                session.close()
        else:
            with self._lock:
                return sum(1 for j in self._memory_jobs.values() if j.get("status") == "running")

    def job_exists(self, job_id: str) -> bool:
        """Check if job exists."""
        if self._use_db:
            session = self._Session()
            try:
                return session.query(JobRecord).filter_by(job_id=job_id).count() > 0
            finally:
                session.close()
        else:
            with self._lock:
                return job_id in self._memory_jobs

    # ==================== Results Storage ====================

    def save_results(self, job_id: str, results: dict) -> None:
        """Save job results."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                if record:
                    record.results_json = json.dumps(results, default=str)
                    session.commit()
            except Exception:
                session.rollback()
                raise
            finally:
                session.close()
        else:
            with self._lock:
                self._memory_results[job_id] = results

    def get_results(self, job_id: str) -> Optional[dict]:
        """Get job results."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                return record.get_results() if record else None
            finally:
                session.close()
        else:
            with self._lock:
                results = self._memory_results.get(job_id)
                return results.copy() if results else None

    def delete_results(self, job_id: str) -> None:
        """Delete job results."""
        if self._use_db:
            session = self._Session()
            try:
                record = session.query(JobRecord).filter_by(job_id=job_id).first()
                if record:
                    record.results_json = None
                    session.commit()
            except Exception:
                session.rollback()
                raise
            finally:
                session.close()
        else:
            with self._lock:
                self._memory_results.pop(job_id, None)

    # ==================== Maintenance ====================

    def _evict_if_needed(self) -> None:
        """Evict oldest completed/failed jobs if over limit (memory mode only)."""
        if len(self._memory_jobs) <= self._max_jobs:
            return
        completed = [
            (jid, j) for jid, j in self._memory_jobs.items()
            if j.get("status") in ("completed", "failed")
        ]
        completed.sort(key=lambda x: x[1].get("updated_at", ""))
        to_remove = len(self._memory_jobs) - self._max_jobs
        for jid, _ in completed[:to_remove]:
            del self._memory_jobs[jid]
            self._memory_results.pop(jid, None)
