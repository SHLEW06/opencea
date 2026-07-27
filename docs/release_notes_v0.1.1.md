# OpenCEA 0.1.1 release candidate

OpenCEA 0.1.1 prepares the first public package upload. It contains the
typed public API work completed after the local `v0.1.0` tag and repairs
the package validation path. The model calculations are unchanged.

## Release integrity

- Direct test, lint, type, build, and publishing tools are pinned.
  Ruff uses version `0.15.21` and the same explicit rules in local
  checks, pre-commit, and CI.
- Package metadata uses the MIT SPDX expression, includes the license
  file and PEP 561 marker, and declares support for Python 3.10 through
  3.12.
- The source distribution includes its tests, YAML fixtures, examples,
  and release documents.
- CI tests the unpacked sdist and installs the wheel in a separate
  environment. The publishing workflow repeats the release gates before
  upload.
- The installed-package example is self-contained and returns total cost
  `244.00` and total QALYs `2.44`.

## Validation target

The release candidate must pass:

```text
python -m pytest --cov
ruff check src tests examples scripts
ruff format --check src tests examples scripts
python -m mypy src/opencea
python -m build
python -m twine check dist/*
python scripts/check_distribution.py dist
```

It must also pass the tests from the unpacked sdist and the smoke
calculation after a clean wheel installation. See
`docs/RELEASE_CHECKLIST.md` for the complete procedure.

## Candidate validation

The local release gate was run on Python 3.12.13 on July 27, 2026:

- The source checkout and unpacked sdist each passed 122 tests with no
  skips. Coverage was 95.86% against the 90% threshold.
- Ruff `0.15.21` passed lint and format checks. MyPy `2.3.0` passed on
  all 10 package source files.
- Pre-commit passed every configured hook on the scoped release files.
- The pinned setuptools backend built one sdist and one
  `py3-none-any` wheel without warnings. Twine and the archive validator
  accepted both files.
- A fresh environment installed the wheel and its runtime dependencies.
  `pip check` found no broken requirements, and the smoke calculation
  returned total cost `244.00` and total QALYs `2.44`.

CI remains responsible for repeating the suite on Python 3.10 and 3.11
before a release commit can be tagged.

## Publication status

Version `0.1.1` is a release candidate. It has not been uploaded to PyPI
and no GitHub Release has been created.
