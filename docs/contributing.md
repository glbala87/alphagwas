# Contributing to AlphaGWAS

Thank you for your interest in contributing to AlphaGWAS!

## Development Setup

```bash
# Clone the repository
git clone https://github.com/glbala87/alphagwas.git
cd alphagwas

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install with development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pip install pre-commit
pre-commit install
```

## Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=scripts --cov-report=html

# Run specific test file
pytest tests/test_utils.py -v
```

## Code Style

We use the following tools for code quality:

- **black** - Code formatting
- **isort** - Import sorting
- **flake8** - Linting
- **mypy** - Type checking

Pre-commit hooks will run these automatically on commit.

### Manual formatting

```bash
# Format code
black scripts/ tests/
isort scripts/ tests/

# Check linting
flake8 scripts/ tests/

# Type check
mypy scripts/
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes
4. Run tests: `pytest tests/ -v`
5. Commit with descriptive message
6. Push to your fork
7. Open a Pull Request

## Commit Messages

Use conventional commit format:

```
feat: add new visualization module
fix: handle missing values in scoring
docs: update installation guide
test: add tests for enrichment module
refactor: simplify prediction batching
```

## Reporting Issues

Please include:

- Python version
- AlphaGWAS version
- Steps to reproduce
- Expected vs actual behavior
- Error messages/tracebacks

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
