# Contributing to OpenCEA

## Development setup

Use Python 3.10, 3.11, or 3.12.

```bash
git clone https://github.com/SHLEW06/OpenCEA.git
cd OpenCEA
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade "pip==26.1.2"
python -m pip install -e ".[dev]"
pre-commit install
```

The `dev` extra pins the direct test, lint, type, build, and publishing
tools. Pre-commit and CI use the same Ruff version and rule set.

## Running the checks

```bash
python -m pytest --cov
ruff check src tests examples scripts
ruff format --check src tests examples scripts
python -m mypy src/opencea
python -m build
python -m twine check dist/*
python scripts/check_distribution.py dist
```

CI runs the test suite on every supported Python version. It also tests
the unpacked source distribution and installs the wheel in a separate
environment before running the smoke calculation.

## Ground rules

- **The DARTH golden tests are the gate.** `tests/test_darth_reference.py` pins the engine to the published manuscript (Tables 5 and 6, to the cent and to the dollar). Never weaken or delete these; any change that moves them is a domain-behavior change and needs explicit justification.
- **Don't change numerical behavior in tooling or refactor PRs.** Discounting, within-cycle correction, dominance logic, and PSA sampling are validated; keep behavioral changes in their own clearly-labeled PRs.
- Add tests for new functionality; coverage should not drop below the configured threshold.
- Follow the existing style — Ruff handles formatting, so don't hand-format.

## Releases

OpenCEA uses semantic versioning. Maintainers must follow
[`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md). Publishing uses
PyPI Trusted Publishing and starts only when a GitHub Release is
published.
