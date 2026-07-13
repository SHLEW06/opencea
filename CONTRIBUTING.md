# Contributing to OpenCEA

## Development setup

Requires Python ≥ 3.10.

```bash
git clone https://github.com/SHLEW06/OpenCEA.git
cd OpenCEA
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install          # optional but recommended
```

## Running the checks

```bash
pytest                      # full suite (117 tests, < 5 s)
pytest --cov=opencea        # with coverage; CI enforces the fail_under
                            # threshold in pyproject.toml
ruff check src tests        # lint
ruff format src tests       # format (CI runs `ruff format --check`)
```

All of these run in CI (`.github/workflows/ci.yml`) on Python 3.10–3.12; a PR must pass lint, format check, and the full test suite with coverage.

## Ground rules

- **The DARTH golden tests are the gate.** `tests/test_darth_reference.py` pins the engine to the published manuscript (Tables 5 and 6, to the cent and to the dollar). Never weaken or delete these; any change that moves them is a domain-behavior change and needs explicit justification.
- **Don't change numerical behavior in tooling or refactor PRs.** Discounting, within-cycle correction, dominance logic, and PSA sampling are validated; keep behavioral changes in their own clearly-labeled PRs.
- Add tests for new functionality; coverage should not drop below the configured threshold.
- Follow the existing style — Ruff handles formatting, so don't hand-format.

## Releases

Versioning is semantic. Publishing runs via PyPI Trusted Publishing on GitHub release publication (see `docs/PUBLISHING.md`); maintainers only.
