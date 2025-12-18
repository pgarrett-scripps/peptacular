.PHONY: help install install-dev install-docs install-all sync test test-v test-vv test-cov test-file test-watch clean lint format check all docs docs-clean docs-open

help:  ## Show this help message
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install dependencies using uv
	uv sync

install-dev:  ## Install with dev dependencies
	uv sync --group dev

install-docs:  ## Install with docs dependencies
	uv sync --group docs

install-all:  ## Install with all dependency groups
	uv sync --all-groups

sync:  ## Sync dependencies (same as install)
	uv sync

test:  ## Run tests with pytest
	uv run pytest tests/

test-v:  ## Run tests with verbose output
	uv run pytest tests/ -v

test-vv:  ## Run tests with very verbose output
	uv run pytest tests/ -vv

test-cov:  ## Run tests with coverage report
	uv run pytest tests/ --cov=src/peptacular --cov-report=term-missing

test-file:  ## Run tests for a specific file (usage: make test-file FILE=test_mod_parser.py)
	uv run pytest tests/$(FILE) -v

clean:  ## Clean up cache files and build artifacts
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name ".coverage" -delete

lint:  ## Run linter (if ruff is installed)
	uv run ruff check src/ tests/ || echo "Ruff not installed, skipping lint"

format:  ## Format code (if ruff is installed)
	uv run ruff format src/ tests/ || echo "Ruff not installed, skipping format"

docs:  ## Build documentation with Sphinx
	cd docs && uv run sphinx-build -b html . _build/html

docs-clean:  ## Clean documentation build artifacts
	rm -rf docs/_build

docs-open:  ## Build and open documentation in browser
	$(MAKE) docs
	@python -c "import webbrowser; webbrowser.open('file://$(shell pwd)/docs/_build/html/index.html')"

check:  ## Run linter and tests
	$(MAKE) lint
	$(MAKE) test

all: clean install test  ## Clean, install dependencies, and run tests
