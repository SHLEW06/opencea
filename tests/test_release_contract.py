"""Regression tests for version, package data, and release documentation."""

from __future__ import annotations

import json
import re
import runpy
from importlib.metadata import version as installed_version
from pathlib import Path

import pytest
import yaml

import opencea
from opencea import WaningSpec, breakeven_drug_price, scenario_icer

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
    quickstart = readme.split("## Five-minute example", 1)[1]
    snippet = quickstart.split("```python", 1)[1].split("```", 1)[0]
    exec(compile(snippet, "README.md", "exec"), {})

    assert capsys.readouterr().out.splitlines() == [
        "Base-case ICER: $98,900/QALY",
        "Waning-effect ICER: $206,778/QALY",
        "Sustained-effect break-even price: $6,355/year",
    ]


def test_readme_links_resolve_or_are_external() -> None:
    readme = (ROOT / "README.md").read_text()
    targets = re.findall(r"!?\[[^]]*]\(([^)]+)\)", readme)
    for target in targets:
        if target.startswith(("https://", "http://", "#")):
            continue
        local_target = target.split("#", 1)[0]
        assert (ROOT / local_target).exists(), f"broken README link: {target}"


def test_readme_presentation_assets_exist() -> None:
    assets = ROOT / "docs" / "assets"
    for name in (
        "empagliflozin-ceac.png",
        "empagliflozin-results.json",
        "empagliflozin-scenario-icers.png",
        "empagliflozin-tornado.png",
    ):
        path = assets / name
        assert path.is_file(), f"missing generated README asset: {name}"
        assert path.stat().st_size > 0


def test_generated_case_study_summary_matches_model_and_readme() -> None:
    summary = json.loads(
        (ROOT / "docs" / "assets" / "empagliflozin-results.json").read_text()
    )
    params = ROOT / "examples" / "empagliflozin_t2d.yaml"
    waning = WaningSpec(start_year=3.0, end_year=10.0)
    expected_icers = {
        "wac_sustained": scenario_icer(params),
        "net_sustained": scenario_icer(params, drug_price=4_500.0),
        "wac_waning": scenario_icer(params, waning=waning),
        "net_waning": scenario_icer(params, drug_price=4_500.0, waning=waning),
    }

    readme = (ROOT / "README.md").read_text()
    for scenario in summary["scenarios"]:
        expected = expected_icers[scenario["name"]]
        assert scenario["icer_per_qaly"] == pytest.approx(expected, abs=1e-9)
        assert f"{expected:,.0f}" in readme
        assert f"{scenario['probability_cost_effective']:.3f}" in readme

    break_even = summary["break_even"]
    expected_sustained = breakeven_drug_price(params, target_icer=100_000.0)
    expected_waning = breakeven_drug_price(
        params,
        target_icer=100_000.0,
        waning=waning,
    )
    assert break_even["sustained_annual_drug_price"] == pytest.approx(
        expected_sustained,
        abs=1e-9,
    )
    assert break_even["waning_annual_drug_price"] == pytest.approx(
        expected_waning,
        abs=1e-9,
    )
    assert f"{expected_sustained:,.0f}" in readme
    assert f"{expected_waning:,.0f}" in readme
