.PHONY: help install install-dev install-docs install-all sync test test-v test-vv test-cov test-file test-watch clean lint format check all docs docs-clean docs-open paper

help: 
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install: 
	uv sync

install-dev:
	uv sync --group dev

install-docs:
	uv sync --group docs

install-all:
	uv sync --all-groups

sync:
	uv sync

test:
	uv run pytest tests/

test-v:
	uv run pytest tests/ -v

test-vv:
	uv run pytest tests/ -vv

test-cov:
	uv run pytest tests/ --cov=src/peptacular --cov-report=term-missing

test-file:
	uv run pytest tests/$(FILE) -v

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name ".coverage" -delete

lint: 
	uv run ruff check src/ tests/ || echo "Ruff not installed, skipping lint"

format:
	uv run ruff format src/ tests/ || echo "Ruff not installed, skipping format"

docs:
	cd docs && uv run sphinx-build -b html . _build/html

docs-clean:
	rm -rf docs/_build

docs-open:
	$(MAKE) docs
	@python -c "import webbrowser; webbrowser.open('file://$(shell pwd)/docs/_build/html/index.html')"

paper:
	docker run --rm \
		--volume $(PWD)/paper:/data \
		--user $(shell id -u):$(shell id -g) \
		--env JOURNAL=joss \
		openjournals/inara

check:
	$(MAKE) lint
	$(MAKE) test

all: clean install test
