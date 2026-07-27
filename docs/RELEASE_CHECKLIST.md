# OpenCEA 0.1.1 release checklist

Use this checklist for the `0.1.1` release candidate. A maintainer should
perform the Git and PyPI steps only after reviewing the uncommitted
changes and after CI passes.

## Confirm the release state

- [ ] The working tree contains only reviewed release changes.
- [ ] `pyproject.toml`, `src/opencea/__init__.py`, and `CITATION.cff`
  all contain version `0.1.1`.
- [ ] `docs/release_notes_v0.1.1.md` and the `0.1.1` changelog entry
  describe the candidate.
- [ ] `CITATION.cff` has the actual `date-released` value added on the
  release date.
- [ ] Python support is stated as 3.10 through 3.12 in metadata,
  documentation, and CI.
- [ ] The PyPI Trusted Publisher and GitHub `pypi` environment use the
  values in `docs/PUBLISHING.md`.
- [ ] The remote does not already contain `v0.1.1`.
- [ ] The PyPI index does not already contain version `0.1.1`.

## Reproduce the checks

Create a fresh development environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade "pip==26.1.2"
python -m pip install -e ".[dev]"
```

Run the source checks:

```bash
python -m pytest --cov
ruff check src tests examples scripts
ruff format --check src tests examples scripts
python -m mypy src/opencea
pre-commit run --all-files
```

Build and inspect the release archives from a clean checkout:

```bash
python -m build
python -m twine check dist/*
python scripts/check_distribution.py dist
```

Test the unpacked sdist:

```bash
mkdir -p /tmp/opencea-sdist-0.1.1
tar -xzf dist/opencea-0.1.1.tar.gz -C /tmp/opencea-sdist-0.1.1
cd /tmp/opencea-sdist-0.1.1/opencea-0.1.1
python3.12 -m venv .venv
.venv/bin/python -m pip install --upgrade "pip==26.1.2"
.venv/bin/python -m pip install -e ".[dev]"
.venv/bin/python -m pytest --cov
```

Return to the repository, then install and test only the wheel:

```bash
python3.12 -m venv /tmp/opencea-wheel-0.1.1
/tmp/opencea-wheel-0.1.1/bin/python -m pip install --upgrade "pip==26.1.2"
/tmp/opencea-wheel-0.1.1/bin/python -m pip install dist/opencea-0.1.1-py3-none-any.whl
/tmp/opencea-wheel-0.1.1/bin/python -m pip check
/tmp/opencea-wheel-0.1.1/bin/python examples/installed_smoke.py
```

The smoke script must report version `0.1.1`, total cost `244.00`, and
total QALYs `2.44`.

## Tag and publish

- [ ] Commit the reviewed release changes.
- [ ] Push the branch and wait for every CI job, including
  `build and validate distributions`, to pass.
- [ ] Create the annotated tag `v0.1.1` on that exact commit.
- [ ] Push only that tag.
- [ ] Create a GitHub Release from `v0.1.1` and paste
  `docs/release_notes_v0.1.1.md`.
- [ ] Publish the GitHub Release.
- [ ] Confirm that the `Publish to PyPI` workflow pauses for any required
  environment approval, passes its checks, and uploads both archives.

## Verify the public release

- [ ] Open the PyPI project page and check the description, image,
  project URLs, MIT license expression, Python range, and files.
- [ ] Create another clean environment and run:

```bash
python3.12 -m venv /tmp/opencea-pypi-0.1.1
/tmp/opencea-pypi-0.1.1/bin/python -m pip install --no-cache-dir "opencea==0.1.1"
/tmp/opencea-pypi-0.1.1/bin/python -m pip check
```

- [ ] Run the README quickstart with the package installed from PyPI.
- [ ] Update the README installation text only after the public install
  succeeds.
