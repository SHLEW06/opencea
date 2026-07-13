"""Shared public result contracts used by OpenCEA APIs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, TypeAlias, TypedDict

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.float64]
ParameterInput: TypeAlias = Mapping[str, Any] | str | Path


class StrategyResult(TypedDict):
    """Deterministic results returned for one evaluated strategy."""

    name: str
    trace: FloatArray
    cycle_costs: FloatArray
    cycle_qalys: FloatArray
    total_cost: float
    total_qaly: float


class CEAStrategyResult(TypedDict):
    """The subset of a strategy result consumed by ``cea_table``."""

    name: str
    total_cost: float
    total_qaly: float


class CaseStrategyResult(TypedDict):
    """Deterministic empagliflozin case-study results for one strategy."""

    total_cost: float
    total_qaly: float
    state_cost: float
    acute_event_cost: float
    trace: FloatArray


class DSAOutcome(TypedDict):
    """Cost and QALY totals expected by a DSA evaluator."""

    cost: float
    qaly: float
