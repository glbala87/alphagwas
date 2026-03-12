# AlphaGWAS Makefile
# Convenience commands for development and deployment

.PHONY: help install test lint docker-build docker-run clean

# Default target
help:
	@echo "AlphaGWAS - Available commands:"
	@echo ""
	@echo "  make install       Install dependencies"
	@echo "  make install-dev   Install development dependencies"
	@echo "  make test          Run tests"
	@echo "  make test-cov      Run tests with coverage"
	@echo "  make lint          Run linting"
	@echo "  make format        Format code with black"
	@echo ""
	@echo "  make docker-build  Build Docker image"
	@echo "  make docker-run    Run pipeline in Docker"
	@echo "  make docker-test   Run tests in Docker"
	@echo "  make docker-shell  Interactive Docker shell"
	@echo ""
	@echo "  make clean         Clean temporary files"
	@echo "  make clean-all     Clean everything including outputs"
	@echo ""
	@echo "  make check         Verify setup and config"
	@echo "  make run           Run full pipeline"
	@echo "  make visualize     Generate visualizations"

# Installation
install:
	pip install -r requirements.txt

install-dev:
	pip install -r requirements.txt
	pip install pytest pytest-cov black flake8 mypy

# Testing
test:
	python -m pytest tests/ -v

test-cov:
	python -m pytest tests/ -v --cov=scripts --cov-report=html --cov-report=term

test-fast:
	python -m pytest tests/ -v -x --tb=short

# Linting and formatting
lint:
	flake8 scripts/ tests/ --max-line-length=100 --ignore=E501,W503

format:
	black scripts/ tests/ --line-length=100

type-check:
	mypy scripts/ --ignore-missing-imports

# Docker commands
docker-build:
	docker build -t alphagwas:latest .

docker-run:
	docker-compose run --rm alphagwas --config /app/config/config.yaml

docker-test:
	docker-compose run --rm alphagwas-test

docker-shell:
	docker-compose run --rm alphagwas-dev /bin/bash

docker-clean:
	docker-compose down --rmi local --volumes

# Pipeline commands
check:
	python run_pipeline.py --check

run:
	python run_pipeline.py --config config/config.yaml

run-step-%:
	python run_pipeline.py --config config/config.yaml --step $*

visualize:
	python run_pipeline.py --visualize

# Cleaning
clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	find . -type f -name "*.pyo" -delete 2>/dev/null || true
	find . -type f -name ".coverage" -delete 2>/dev/null || true
	rm -rf htmlcov/ .mypy_cache/ 2>/dev/null || true

clean-intermediate:
	rm -rf data/intermediate/* 2>/dev/null || true
	touch data/intermediate/.gitkeep

clean-output:
	rm -rf data/output/* 2>/dev/null || true
	touch data/output/.gitkeep

clean-all: clean clean-intermediate clean-output

# Development
dev-setup: install-dev
	@echo "Development environment ready!"
	@echo "Run 'make test' to verify setup."

# CI/CD helpers
ci-test:
	python -m pytest tests/ -v --junitxml=test-results.xml

ci-lint:
	flake8 scripts/ tests/ --max-line-length=100 --ignore=E501,W503 --format=pylint > lint-results.txt || true
