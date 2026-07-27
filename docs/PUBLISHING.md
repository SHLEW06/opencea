# Publishing OpenCEA to PyPI

OpenCEA publishes through PyPI Trusted Publishing. GitHub exchanges an
OpenID Connect token for short-lived upload access, so the repository
does not need a PyPI API token.

Publication starts only when a GitHub Release is published. The workflow
checks that the release tag matches the version in `pyproject.toml`,
runs the release gates, builds one sdist and one wheel, inspects both
archives, and installs the wheel in a clean environment. The upload step
runs only after those checks pass.

## One-time setup

1. Create a PyPI account and enable two-factor authentication.
2. In PyPI, add a pending Trusted Publisher with these values:
   - PyPI project name: `opencea`
   - Owner: `SHLEW06`
   - Repository name: `OpenCEA`
   - Workflow name: `publish.yml`
   - Environment name: `pypi`
3. In the GitHub repository, create an environment named `pypi`.
4. Add a required reviewer to that environment if publication should
   pause for human approval.

The capitalization of the owner, repository, workflow, and environment
should match the repository and workflow configuration.

## Release procedure

Follow [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) for every release.
It covers version synchronization, local checks, archive inspection,
clean installation, tagging, publication, and post-publication
verification.

The publishing workflow must not be run from an arbitrary branch or an
untagged commit. Do not bypass a failed validation step. PyPI versions
cannot be replaced after upload, so a failed or incorrect version needs
a new version number.
