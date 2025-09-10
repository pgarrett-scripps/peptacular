.PHONY: help install install-dev sync test test-cov test-watch lint format clean build docs run-gen

# Default target
help:
	@echo "Available commands:"
	@echo "  install      - Install production dependencies"
	@echo "  install-dev  - Install all dependencies (including dev)"
	@echo "  sync         - Sync dependencies with lock file"
	@echo "  test         - Run tests with pytest"
	@echo "  test-cov     - Run tests with coverage"
	@echo "  test-watch   - Run tests in watch mode"
	@echo "  lint         - Run linting (flake8)"
	@echo "  format       - Format code (black)"
	@echo "  clean        - Clean cache and build files"
	@echo "  build        - Build package"
	@echo "  docs         - Build documentation"
	@echo "  run-gen      - Run generate_mod_dbs.py"

# Dependency management
install:
	uv sync --no-dev

install-dev:
	uv sync --dev

sync:
	uv sync

# Testing
test:
	uv run pytest

test-cov:
	uv run pytest --cov=src --cov-report=html --cov-report=term

test-watch:
	uv run pytest --watch

# Code quality
lint:
	uv run ruff check src

format:
	uv run ruff format src

format-check:
	uv run ruff format --check src

# Maintenance
clean:
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf .pytest_cache/

# Build
build:
	uv build

# Documentation (if using sphinx)
docs:
	uv run sphinx-build -b html docs/ docs/_build/

# Project-specific commands
run-gen:
	uv run python generate_mod_dbs.py

# Development setup (run this after cloning)
dev-setup: install-dev
	@echo "Development environment ready!"