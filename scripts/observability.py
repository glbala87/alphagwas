"""
Observability Module for AlphaGWAS.

Provides structured logging, metrics collection, request tracing,
and health check capabilities for production monitoring.
"""

import contextvars
import logging
import os
import shutil
import threading
import time
import uuid
from collections import defaultdict, deque
from datetime import datetime
from typing import Any, Optional

logger = logging.getLogger(__name__)

# Context variable for request ID propagation
request_id_var: contextvars.ContextVar[str] = contextvars.ContextVar("request_id", default="")


# ---------------------------------------------------------------------------
# Structured Logging
# ---------------------------------------------------------------------------

class StructuredFormatter(logging.Formatter):
    """JSON-structured log formatter with request ID support."""

    def format(self, record: logging.LogRecord) -> str:
        import json

        log_data = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }

        # Add request ID if available
        req_id = request_id_var.get("")
        if req_id:
            log_data["request_id"] = req_id

        # Add exception info if present
        if record.exc_info and record.exc_info[0]:
            log_data["exception"] = {
                "type": record.exc_info[0].__name__,
                "message": str(record.exc_info[1]),
            }

        # Add extra fields
        for key in ("status_code", "method", "path", "duration_ms", "client_ip"):
            if hasattr(record, key):
                log_data[key] = getattr(record, key)

        return json.dumps(log_data)


def setup_structured_logging(
    level: str = "INFO",
    json_format: bool = True,
) -> logging.Logger:
    """
    Configure structured logging for the application.

    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR)
        json_format: Use JSON format (True) or plain text (False)

    Returns:
        Root logger
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Remove existing handlers
    root_logger.handlers.clear()

    handler = logging.StreamHandler()
    if json_format:
        handler.setFormatter(StructuredFormatter())
    else:
        handler.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s (%(request_id)s): %(message)s",
            defaults={"request_id": "-"},
        ))

    root_logger.addHandler(handler)
    return root_logger


# ---------------------------------------------------------------------------
# Metrics Collector
# ---------------------------------------------------------------------------

class MetricsCollector:
    """
    Thread-safe in-process metrics collector.

    Tracks counters, gauges, and histograms for application monitoring.
    Exposes metrics via get_metrics() for the /metrics endpoint.
    """

    def __init__(self, max_observations: int = 1000):
        self._lock = threading.Lock()
        self._counters: dict[str, float] = defaultdict(float)
        self._gauges: dict[str, float] = defaultdict(float)
        self._histograms: dict[str, deque] = defaultdict(lambda: deque(maxlen=max_observations))
        self._start_time = time.time()

    def inc(self, metric: str, value: float = 1.0) -> None:
        """Increment a counter."""
        with self._lock:
            self._counters[metric] += value

    def dec(self, metric: str, value: float = 1.0) -> None:
        """Decrement a gauge."""
        with self._lock:
            self._gauges[metric] -= value

    def set_gauge(self, metric: str, value: float) -> None:
        """Set a gauge to a specific value."""
        with self._lock:
            self._gauges[metric] = value

    def observe(self, metric: str, value: float) -> None:
        """Record an observation (for histogram-like metrics)."""
        with self._lock:
            self._histograms[metric].append(value)

    def get_metrics(self) -> dict[str, Any]:
        """Get all metrics as a dictionary."""
        with self._lock:
            result = {
                "uptime_seconds": round(time.time() - self._start_time, 1),
                "collected_at": datetime.utcnow().isoformat() + "Z",
                "counters": dict(self._counters),
                "gauges": dict(self._gauges),
                "histograms": {},
            }

            for name, values in self._histograms.items():
                if values:
                    sorted_vals = sorted(values)
                    n = len(sorted_vals)
                    result["histograms"][name] = {
                        "count": n,
                        "sum": round(sum(sorted_vals), 4),
                        "min": round(sorted_vals[0], 4),
                        "max": round(sorted_vals[-1], 4),
                        "avg": round(sum(sorted_vals) / n, 4),
                        "p50": round(sorted_vals[n // 2], 4),
                        "p95": round(sorted_vals[int(n * 0.95)], 4) if n >= 20 else None,
                        "p99": round(sorted_vals[int(n * 0.99)], 4) if n >= 100 else None,
                    }

            return result


# Global metrics collector instance
metrics_collector = MetricsCollector()


# ---------------------------------------------------------------------------
# Request ID Middleware
# ---------------------------------------------------------------------------

try:
    from starlette.middleware.base import BaseHTTPMiddleware
    from starlette.requests import Request
    from starlette.responses import Response

    class RequestIdMiddleware(BaseHTTPMiddleware):
        """
        Middleware that assigns a unique request ID to each request.

        - Reads X-Request-ID header or generates a new UUID
        - Stores in contextvars for log correlation
        - Adds to response headers
        - Logs request method, path, status, and duration
        """

        async def dispatch(self, request: Request, call_next) -> Response:
            req_id = request.headers.get("X-Request-ID", str(uuid.uuid4())[:8])
            request_id_var.set(req_id)

            start = time.time()
            method = request.method
            path = request.url.path

            try:
                response = await call_next(request)
                duration_ms = round((time.time() - start) * 1000, 1)

                response.headers["X-Request-ID"] = req_id
                response.headers["X-Response-Time"] = f"{duration_ms}ms"

                # Record metrics
                metrics_collector.observe("request_duration_seconds", duration_ms / 1000)
                metrics_collector.inc("requests_total")

                # Log request
                logger.info(
                    f"{method} {path} {response.status_code} {duration_ms}ms",
                    extra={
                        "method": method,
                        "path": path,
                        "status_code": response.status_code,
                        "duration_ms": duration_ms,
                        "client_ip": request.client.host if request.client else "",
                    },
                )

                return response

            except Exception as e:
                duration_ms = round((time.time() - start) * 1000, 1)
                metrics_collector.inc("requests_error_total")
                logger.error(
                    f"{method} {path} 500 {duration_ms}ms - {e}",
                    extra={"method": method, "path": path, "duration_ms": duration_ms},
                )
                raise

    HAS_MIDDLEWARE = True

except ImportError:
    HAS_MIDDLEWARE = False
    RequestIdMiddleware = None


# ---------------------------------------------------------------------------
# Health Check
# ---------------------------------------------------------------------------

class DetailedHealthCheck:
    """
    Comprehensive health check that verifies system components.

    Checks database connectivity, disk space, and memory usage.
    Returns structured health info for monitoring dashboards.
    """

    def __init__(self, db_url: Optional[str] = None):
        self._start_time = time.time()
        self._db_url = db_url or os.getenv("DATABASE_URL")

    def check(self) -> dict[str, Any]:
        """Run all health checks and return results."""
        checks = {}
        overall_status = "healthy"

        # Database check
        checks["database"] = self._check_database()
        if checks["database"]["status"] == "unhealthy":
            overall_status = "degraded"

        # Disk space check
        checks["disk"] = self._check_disk_space()
        if checks["disk"]["status"] == "unhealthy":
            overall_status = "unhealthy"

        # Memory check
        checks["memory"] = self._check_memory()
        if checks["memory"]["status"] == "unhealthy":
            overall_status = "unhealthy"

        return {
            "status": overall_status,
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "uptime_seconds": round(time.time() - self._start_time, 1),
            "checks": checks,
        }

    def _check_database(self) -> dict:
        """Check database connectivity."""
        if not self._db_url:
            return {"status": "skipped", "message": "No DATABASE_URL configured"}

        try:
            from sqlalchemy import create_engine, text

            engine = create_engine(self._db_url)
            with engine.connect() as conn:
                conn.execute(text("SELECT 1"))
            engine.dispose()
            return {"status": "healthy", "message": "Connected"}
        except ImportError:
            return {"status": "skipped", "message": "SQLAlchemy not installed"}
        except Exception as e:
            return {"status": "unhealthy", "message": str(e)}

    def _check_disk_space(self) -> dict:
        """Check available disk space in temp directory."""
        import tempfile

        try:
            usage = shutil.disk_usage(tempfile.gettempdir())
            free_gb = round(usage.free / (1024 ** 3), 2)
            pct_free = round(usage.free / usage.total * 100, 1)

            status = "healthy"
            if pct_free < 5:
                status = "unhealthy"
            elif pct_free < 15:
                status = "degraded"

            return {
                "status": status,
                "free_gb": free_gb,
                "percent_free": pct_free,
            }
        except Exception as e:
            return {"status": "unhealthy", "message": str(e)}

    def _check_memory(self) -> dict:
        """Check process memory usage."""
        try:
            import resource

            usage_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            # macOS reports in bytes, Linux in KB
            import platform
            if platform.system() == "Darwin":
                usage_mb = round(usage_kb / (1024 * 1024), 1)
            else:
                usage_mb = round(usage_kb / 1024, 1)

            status = "healthy" if usage_mb < 2048 else "degraded" if usage_mb < 4096 else "unhealthy"

            return {
                "status": status,
                "rss_mb": usage_mb,
            }
        except Exception as e:
            return {"status": "skipped", "message": str(e)}


def create_metrics_endpoint(app: Any, collector: MetricsCollector) -> None:
    """Register /metrics endpoint on a FastAPI app."""
    try:
        from fastapi.responses import JSONResponse as FastAPIJSONResponse

        @app.get("/metrics")
        async def metrics():
            return FastAPIJSONResponse(content=collector.get_metrics())

    except ImportError:
        logger.warning("FastAPI not available, skipping metrics endpoint registration")
