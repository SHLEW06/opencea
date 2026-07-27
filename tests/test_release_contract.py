"""Regression tests for version, package data, and release documentation."""

from __future__ import annotations

import re
import runpy
from importlib.metadata import version as installed_version
from pathlib import Path

import yaml

import opencea

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]


def test_release_versions_are_consistent() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text())
    version = project["project"]["version"]

    assert version == opencea.__version__
    assert version == installed_version("opencea")
    assert version == citation["version"]
    assert (ROOT / "docs" / f"release_notes_v{version}.md").is_file()


def test_pep561_marker_is_present() -> None:
    marker = ROOT / "src" / "opencea" / "py.typed"
    assert marker.is_file()
    assert marker.read_text().strip()


def test_installed_smoke_example(capsys) -> None:
    namespace = runpy.run_path(str(ROOT / "examples" / "installed_smoke.py"))
    namespace["main"]()

    output = capsys.readouterr().out
    assert f"OpenCEA {opencea.__version__}" in output
    assert "total_cost=244.00" in output
    assert "total_qaly=2.44" in output


def test_readme_quickstart(capsys) -> None:
    readme = (ROOT / "README.md").read_text()
    quickstart = readme.split("After installation, this calculation", 1)[1]
    snippet = quickstart.split("```python", 1)[1].split("```", 1)[0]
    exec(compile(snippet, "README.md", "exec"), {})

    assert capsys.readouterr().out.strip() == "244.0 2.4400000000000004"


def test_readme_uses_absolute_repository_links() -> None:
    readme = (ROOT / "README.md").read_text()
    targets = re.findall(r"!?\[[^]]*]\(([^)]+)\)", readme)
    assert all(target.startswith(("https://", "http://", "#")) for target in targets)
