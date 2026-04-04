# AlphaGWAS Docker Image
# Multi-stage build for smaller final image

# ============================================
# Stage 1: Builder
# ============================================
FROM python:3.14-slim as builder

WORKDIR /build

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Create virtual environment and install dependencies
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# ============================================
# Stage 2: Runtime
# ============================================
FROM python:3.14-slim as runtime

LABEL maintainer="BalaSubramani Gattu Linga <glbala87@github>"
LABEL description="AlphaGWAS - GWAS Variant Prioritization Pipeline using AlphaGenome"
LABEL version="1.0.0"

# Install runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Set up application directory
WORKDIR /app

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /app/data/input /app/data/intermediate /app/data/output /app/config

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV ALPHAGWAS_HOME=/app

# Create non-root user for security
RUN useradd -m -u 1000 alphagwas && \
    chown -R alphagwas:alphagwas /app

USER alphagwas

# Default command
ENTRYPOINT ["python", "run_pipeline.py"]
CMD ["--help"]

# ============================================
# Health check
# ============================================
HEALTHCHECK --interval=30s --timeout=10s --start-period=10s --retries=3 \
    CMD python -c "import requests; r = requests.get('http://localhost:8000/health', timeout=5); r.raise_for_status()" 2>/dev/null || python -c "import pandas; import numpy; print('OK')" || exit 1
