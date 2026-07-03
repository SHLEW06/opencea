# Publishing OpenCEA to PyPI

The `.github/workflows/publish.yml` workflow builds and uploads OpenCEA to
PyPI whenever a GitHub Release is published. It uses **PyPI Trusted
Publishing** (OpenID Connect), so no long-lived API token is stored in the
repository or in GitHub Secrets.

The one-time human steps below only need to be done once per project.

## 1. Create a PyPI account (once)

1. Go to <https://pypi.org/account/register/> and create an account.
2. Enable 2FA (required for publishing).

## 2. Register OpenCEA as a Trusted Publisher on PyPI (once)

The project name `opencea` is currently unclaimed on PyPI. You will claim
it by adding a pending trusted publisher **before** the first release, so
that the first workflow run creates the project.

1. Sign in to <https://pypi.org>.
2. Go to **Your account → Publishing → Add a new pending publisher**.
3. Fill in:
   - **PyPI project name:** `opencea`
   - **Owner:** `SHLEW06`
   - **Repository name:** `opencea`
   - **Workflow name:** `publish.yml`
   - **Environment name:** `pypi`
4. Click **Add**.

## 3. Create the `pypi` environment on GitHub (once)

1. On GitHub, go to **Settings → Environments → New environment**.
2. Name it `pypi` (must exactly match the workflow and the trusted
   publisher configuration).
3. Optionally add a required reviewer so that publish runs pause for a
   manual approval.

## 4. Publish a release

For every release after the one-time setup:

1. Make sure `pyproject.toml` and `src/opencea/__init__.py` both have the
   new version.
2. Update `CHANGELOG.md` and `CITATION.cff` (`version`, `date-released`).
3. Tag the commit: `git tag -a vX.Y.Z -m "opencea vX.Y.Z"` and push the tag.
4. On GitHub, **Releases → Draft a new release**, pick the tag, paste the
   release notes from `docs/release_notes_vX.Y.Z.md`, and click **Publish
   release**.
5. The `Publish to PyPI` workflow will run automatically, build the sdist
   and wheel, run `twine check`, and upload via Trusted Publishing.
6. After ~30 seconds the release appears on <https://pypi.org/project/opencea/>
   and `pip install opencea` works.

## Local verification before releasing

```bash
python -m build
twine check dist/*
# fresh-venv smoke:
python -m venv /tmp/opencea_smoke
/tmp/opencea_smoke/bin/pip install dist/opencea-*.whl
/tmp/opencea_smoke/bin/python -c "import opencea; print(opencea.__version__)"
```

If `twine check` fails, do not publish — fix the metadata first.
