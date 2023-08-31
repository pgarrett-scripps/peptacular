# Contributing to Peptacular

Thank you for considering contributing to Peptacular! This document lays out the guidelines for contributing to the project, submitting bug reports, suggesting features, and more.

## Table of Contents

- [Getting Started](#getting-started)
- [Setting up Your Development Environment](#setting-up-your-development-environment)
- [Submitting a Bug Report](#submitting-a-bug-report)
- [Feature Requests](#feature-requests)
- [Pull Requests](#pull-requests)
- [Code Style](#code-style)
- [Documentation](#documentation)
- [Tests](#tests)
  
## Getting Started

1. **Fork the Repository**: Click the 'Fork' button to create a fork of the repository.
2. **Clone the Repository**: `https://github.com/pgarrett-scripps/peptacular.git`
3. **Add Upstream Remote**: `git remote add upstream https://github.com/pgarrett-scripps/peptacular.git`

## Setting up Your Development Environment

1. **Python**: Ensure you have Python >= 3.8 installed.
2. **Virtual Environment**: Optionally, set up a virtual environment.
   ```bash
   python -m venv .venv
   ```
   Activate it by running `source .venv/bin/activate` or `.venv\Scripts\activate` (Windows).
3. **Dependencies**: Install dependencies using `pip install -r requirements.txt`.

## Submitting a Bug Report

If you've found a bug, please open an issue in the repository's issue tracker and fill out the issue template for bug reports.

## Feature Requests

New features are welcomed. If you have an idea, please open an issue in the repository's issue tracker and use the feature request template.

## Pull Requests

1. **Fetch Latest Changes**: Fetch changes from the upstream master branch and merge them into your local master branch.
   ```bash
   git fetch upstream
   git merge upstream/master
   ```
2. **Create a Branch**: Create a branch for your work.
   ```bash
   git checkout -b my-feature
   ```
3. **Commit Your Changes**: Make your changes and commit them with a descriptive message.
   ```bash
   git add .
   git commit -m "Implemented XYZ feature"
   ```
4. **Push**: Push the branch to your fork.
   ```bash
   git push origin my-feature
   ```
5. **Create a Pull Request**: Open a pull request on the original repository to merge your forked branch.

## Code Style

- Follow PEP 8 guidelines for Python code.
- Use type annotations.

## Documentation

- Add comments and docstrings to your code.
- Update README.md if your changes introduce new features or alter existing ones.

## Tests

- Write tests to cover new features and changes.
- Run tests locally before submitting a pull request.
- All tests must pass before your pull request can be merged.
  
Thank you for contributing to Peptacular! We appreciate your support!
