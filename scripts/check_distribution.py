"""Validate OpenCEA source and wheel distribution contents."""

from __future__ import annotations

import argparse
import tarfile
import zipfile
from email.parser import BytesParser
from email.policy import default
from pathlib import Path, PurePosixPath

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_NAME = "opencea"
SUPPORTED_PYTHON = (">=3.10", "<3.13")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"distribution check failed: {message}")


def project_version() -> str:
    data = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())
    return str(data["project"]["version"])


def parse_metadata(payload: bytes):
    return BytesParser(policy=default).parsebytes(payload)


def validate_metadata(payload: bytes, version: str) -> None:
    metadata = parse_metadata(payload)
    require(metadata["Name"].lower() == PACKAGE_NAME, "unexpected package name")
    require(metadata["Version"] == version, "artifact version does not match project")
    requires_python = metadata["Requires-Python"]
    require(
        all(bound in requires_python for bound in SUPPORTED_PYTHON),
        f"unexpected Requires-Python: {requires_python}",
    )
    require(metadata["License-Expression"] == "MIT", "missing MIT license expression")
    require("dev" in metadata.get_all("Provides-Extra", []), "missing dev extra")

    requirements = metadata.get_all("Requires-Dist", [])
    for dependency in ("numpy", "pandas", "pydantic", "pyyaml", "matplotlib"):
        require(
            any(item.lower().startswith(dependency) for item in requirements),
            f"missing runtime dependency: {dependency}",
        )

    project_urls = metadata.get_all("Project-URL", [])
    require(
        any("https://github.com/SHLEW06/OpenCEA" in item for item in project_urls),
        "canonical repository URL is missing",
    )


def validate_member_names(names: set[str]) -> None:
    for name in names:
        path = PurePosixPath(name)
        require(not path.is_absolute(), f"absolute archive member: {name}")
        require(".." not in path.parts, f"unsafe archive member: {name}")
        require("reports" not in path.parts, f"reports content included: {name}")
        require(".git" not in path.parts, f"Git metadata included: {name}")
        require(".env" not in path.parts, f"environment file included: {name}")
        require("__pycache__" not in path.parts, f"cache included: {name}")
        require(not name.endswith((".pyc", ".pyo")), f"bytecode included: {name}")


def validate_wheel(path: Path, version: str) -> None:
    dist_info = f"{PACKAGE_NAME}-{version}.dist-info"
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        validate_member_names(names)

        required = {
            "opencea/__init__.py",
            "opencea/py.typed",
            f"{dist_info}/METADATA",
            f"{dist_info}/WHEEL",
            f"{dist_info}/licenses/LICENSE",
        }
        missing = required - names
        require(not missing, f"wheel is missing: {sorted(missing)}")
        require(
            not any(
                PurePosixPath(name).parts[0] in {"examples", "tests"} for name in names
            ),
            "wheel contains repository-only examples or tests",
        )

        validate_metadata(archive.read(f"{dist_info}/METADATA"), version)
        wheel_metadata = archive.read(f"{dist_info}/WHEEL").decode()
        require(
            "Tag: py3-none-any" in wheel_metadata, "wheel is not platform independent"
        )
        require(archive.read("opencea/py.typed"), "PEP 561 marker is empty")


def validate_sdist(path: Path, version: str) -> None:
    root = f"{PACKAGE_NAME}-{version}"
    with tarfile.open(path, "r:gz") as archive:
        names = {member.name for member in archive.getmembers()}
        validate_member_names(names)

        required = {
            f"{root}/CHANGELOG.md",
            f"{root}/CITATION.cff",
            f"{root}/CONTRIBUTING.md",
            f"{root}/LICENSE",
            f"{root}/MANIFEST.in",
            f"{root}/PKG-INFO",
            f"{root}/README.md",
            f"{root}/docs/RELEASE_CHECKLIST.md",
            f"{root}/docs/release_notes_v{version}.md",
            f"{root}/examples/empagliflozin_t2d.yaml",
            f"{root}/examples/installed_smoke.py",
            f"{root}/examples/sick_sicker.yaml",
            f"{root}/pyproject.toml",
            f"{root}/scripts/check_distribution.py",
            f"{root}/src/opencea/py.typed",
            f"{root}/tests/test_release_contract.py",
        }
        missing = required - names
        require(not missing, f"sdist is missing: {sorted(missing)}")

        pkg_info = archive.extractfile(f"{root}/PKG-INFO")
        require(pkg_info is not None, "cannot read sdist PKG-INFO")
        validate_metadata(pkg_info.read(), version)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dist_dir",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "dist",
        help="directory containing one OpenCEA wheel and one source distribution",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    version = project_version()
    dist_dir = args.dist_dir.resolve()
    wheel = dist_dir / f"{PACKAGE_NAME}-{version}-py3-none-any.whl"
    sdist = dist_dir / f"{PACKAGE_NAME}-{version}.tar.gz"

    require(wheel.is_file(), f"missing wheel: {wheel.name}")
    require(sdist.is_file(), f"missing source distribution: {sdist.name}")
    require(
        len(list(dist_dir.glob(f"{PACKAGE_NAME}-*.whl"))) == 1,
        "distribution directory must contain exactly one OpenCEA wheel",
    )
    require(
        len(list(dist_dir.glob(f"{PACKAGE_NAME}-*.tar.gz"))) == 1,
        "distribution directory must contain exactly one OpenCEA sdist",
    )

    validate_wheel(wheel, version)
    validate_sdist(sdist, version)
    print(f"Validated {wheel.name} and {sdist.name}")


if __name__ == "__main__":
    main()
