# REST API

AlphaGWAS provides a FastAPI-based REST interface for programmatic access.

## Starting the API Server

```bash
# Install API dependencies
pip install alphagwas[api]

# Start server
uvicorn scripts.api:app --host 0.0.0.0 --port 8000

# Or with auto-reload for development
uvicorn scripts.api:app --reload
```

## API Documentation

Once running, access interactive documentation at:

- **Swagger UI**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

## Endpoints

### Health Check

```bash
curl http://localhost:8000/health
```

Response:
```json
{
  "status": "healthy",
  "version": "1.0.0",
  "timestamp": "2024-01-15T10:30:00"
}
```

### Submit Analysis Job

```bash
curl -X POST "http://localhost:8000/api/v1/jobs" \
  -H "Content-Type: multipart/form-data" \
  -F "file=@gwas_sumstats.tsv" \
  -F "name=my_analysis" \
  -F "pvalue_threshold=5e-8"
```

Response:
```json
{
  "job_id": "abc12345",
  "name": "my_analysis",
  "status": "running",
  "created_at": "2024-01-15T10:30:00",
  "progress": 0.0
}
```

### Check Job Status

```bash
curl http://localhost:8000/api/v1/jobs/abc12345
```

### Get Results

```bash
curl http://localhost:8000/api/v1/jobs/abc12345/results
```

### Download Results

```bash
# TSV format
curl -o results.tsv "http://localhost:8000/api/v1/jobs/abc12345/download?format=tsv"

# JSON format
curl -o results.json "http://localhost:8000/api/v1/jobs/abc12345/download?format=json"
```

### Score Individual Variants

```bash
curl -X POST "http://localhost:8000/api/v1/variants/score" \
  -H "Content-Type: application/json" \
  -d '[
    {"chromosome": "1", "position": 12345678},
    {"chromosome": "2", "position": 87654321}
  ]'
```

### List Available Tissues

```bash
curl http://localhost:8000/api/v1/tissues
```

## Python Client Example

```python
import requests

API_URL = "http://localhost:8000"

# Submit job
with open("gwas_data.tsv", "rb") as f:
    response = requests.post(
        f"{API_URL}/api/v1/jobs",
        files={"file": f},
        data={"name": "my_analysis"}
    )
job_id = response.json()["job_id"]

# Poll for completion
import time
while True:
    status = requests.get(f"{API_URL}/api/v1/jobs/{job_id}").json()
    if status["status"] == "completed":
        break
    time.sleep(5)

# Get results
results = requests.get(f"{API_URL}/api/v1/jobs/{job_id}/results").json()
print(f"Top variant: {results['top_variant']}")
```

## WebSocket Updates

For real-time job updates:

```python
import websockets
import asyncio

async def monitor_job(job_id):
    async with websockets.connect(f"ws://localhost:8000/ws/jobs/{job_id}") as ws:
        async for message in ws:
            print(message)

asyncio.run(monitor_job("abc12345"))
```

## Docker Deployment

```bash
# Build and run
docker build -t alphagwas-api .
docker run -p 8000:8000 alphagwas-api uvicorn scripts.api:app --host 0.0.0.0

# Or with docker-compose
docker-compose up api
```
