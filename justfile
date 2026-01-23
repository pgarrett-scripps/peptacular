# List available recipes
help:
    @just --list

# Install dependencies
install:
    uv sync

# Install all dependencies
install-all:
    uv sync --all-groups

# Sync dependencies
sync:
    uv sync

# Run tests
test:
    uv run python --version
    uv run pytest tests/

# Run tests with coverage
test-cov:
    uv run pytest tests/ --cov=src/peptacular --cov-report=term-missing


# Clean build artifacts and cache
clean:
    find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
    find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
    find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
    find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
    find . -type f -name "*.pyc" -delete
    find . -type f -name "*.pyo" -delete
    find . -type f -name ".coverage" -delete

# Run linter
lint:
    uv run ruff check src/

# Format code
format:
    uv run ruff check --select I --fix src/ tests/
    uv run ruff check --select F401 --fix src/ tests/
    uv run ruff format src tests

# ty type checking
ty:
    uv run ty check src/

# Build documentation
docs:
    cd docs && uv run sphinx-build -b html . _build/html

docs-test:
    cd docs && uv run sphinx-build -b doctest . _build/doctest

# Clean documentation build
docs-clean:
    rm -rf docs/_build

# Build and open documentation
docs-open:
    just docs
    python -c "import webbrowser; webbrowser.open('file://{{justfile_directory()}}/docs/_build/html/index.html')"

# Build paper with Docker
paper:
    docker run --rm \
        --volume {{justfile_directory()}}/paper:/data \
        --user `id -u`:`id -g` \
        --env JOURNAL=joss \
        openjournals/inara

# Run lint and tests
check:
    just lint
    just test

# Clean, install, and test
all:
    just clean
    just install
    just test

# Upgrade Python syntax to 3.12+
upgrade:
    @echo "Upgrading Python syntax to 3.12+..."
    @find src tests -name "*.py" -type f -exec uv run --python-preference managed pyupgrade --py312-plus {} +
    @echo "Python syntax upgraded to 3.12+"