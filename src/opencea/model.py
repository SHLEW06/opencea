"""Pydantic specification for a cohort state-transition cost-effectiveness model.

A cohort model is fully described by:

  - a list of mutually-exclusive health states,
  - one or more strategies (each carrying its own transition matrix and
    per-state cost / utility vectors so that treatment effects can act on
    transitions, rewards, or both),
  - an initial state distribution,
  - a cycle length and a finite time horizon (in cycles),
  - discount rates for costs and QALYs,
  - and a flag for half-cycle correction.

Validation enforces the standard integrity rules so that an invalid spec
fails fast — well before the engine runs and produces silently wrong
ICERs:

  - transition matrices are square (n_states x n_states),
  - transition rows sum to 1 (within a small tolerance),
  - probabilities are in [0, 1],
  - utilities are in [0, 1],
  - costs are non-negative,
  - vector lengths agree with the declared number of states.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Literal, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

# Tolerance for floating-point validation of probability sums / bounds.
# Loose enough to accept hand-written YAML where rows are written to a few
# decimals, tight enough to catch genuine specification errors.
_PROB_TOL = 1e-6


class Strategy(BaseModel):
    """A single decision strategy.

    Each strategy carries its own transition matrix and per-state reward
    vectors. This keeps treatment effects (on progression and / or on
    costs / utilities) localized inside the strategy rather than scattered
    through the model.
    """

    model_config = ConfigDict(extra="forbid")

    name: str = Field(
        ..., description="Human-readable strategy label, e.g. 'SoC' or 'A'."
    )
    transition_matrix: List[List[float]] = Field(
        ..., description="Square (n_states x n_states) row-stochastic matrix."
    )
    state_costs: List[float] = Field(
        ..., description="Per-cycle cost incurred by occupying each state."
    )
    state_utilities: List[float] = Field(
        ..., description="Per-cycle utility weight (QALY) for each state, in [0, 1]."
    )

    @field_validator("transition_matrix")
    @classmethod
    def _validate_transition_matrix(cls, v: List[List[float]]) -> List[List[float]]:
        arr = np.asarray(v, dtype=float)
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(f"transition_matrix must be square; got shape {arr.shape}")
        if (arr < -_PROB_TOL).any() or (arr > 1 + _PROB_TOL).any():
            raise ValueError("transition probabilities must lie in [0, 1]")
        row_sums = arr.sum(axis=1)
        if not np.allclose(row_sums, 1.0, atol=_PROB_TOL):
            bad = [
                (i, float(s))
                for i, s in enumerate(row_sums)
                if not np.isclose(s, 1.0, atol=_PROB_TOL)
            ]
            raise ValueError(
                f"transition matrix rows must sum to 1; offending rows: {bad}"
            )
        return v

    @field_validator("state_costs")
    @classmethod
    def _validate_costs(cls, v: List[float]) -> List[float]:
        arr = np.asarray(v, dtype=float)
        if (arr < -_PROB_TOL).any():
            raise ValueError("state_costs must be non-negative")
        return v

    @field_validator("state_utilities")
    @classmethod
    def _validate_utilities(cls, v: List[float]) -> List[float]:
        arr = np.asarray(v, dtype=float)
        if (arr < -_PROB_TOL).any() or (arr > 1 + _PROB_TOL).any():
            raise ValueError("state_utilities must lie in [0, 1]")
        return v

    @model_validator(mode="after")
    def _validate_dimensions(self) -> "Strategy":
        n = len(self.transition_matrix)
        if len(self.state_costs) != n:
            raise ValueError(
                f"state_costs has length {len(self.state_costs)} but transition_matrix is {n}x{n}"
            )
        if len(self.state_utilities) != n:
            raise ValueError(
                f"state_utilities has length {len(self.state_utilities)} but transition_matrix is {n}x{n}"
            )
        return self

    def n_states(self) -> int:
        return len(self.transition_matrix)


class CohortModel(BaseModel):
    """A cohort state-transition cost-effectiveness model.

    The model is time-homogeneous in stage 1 (a single transition matrix
    per strategy applies at every cycle). Time-varying transitions are a
    natural extension for a later stage.
    """

    model_config = ConfigDict(extra="forbid")

    states: List[str] = Field(..., description="Ordered names of the health states.")
    strategies: List[Strategy] = Field(
        ..., min_length=1, description="One or more comparator strategies."
    )
    initial_distribution: List[float] = Field(
        ...,
        description="Starting proportion of the cohort in each state; must sum to 1.",
    )
    cycle_length: float = Field(1.0, gt=0, description="Cycle length in years.")
    time_horizon: int = Field(..., gt=0, description="Number of cycles to simulate.")
    discount_rate_costs: float = Field(
        0.03, ge=0, lt=1, description="Annual discount rate applied to costs."
    )
    discount_rate_qalys: float = Field(
        0.03, ge=0, lt=1, description="Annual discount rate applied to QALYs."
    )
    wcc_method: Literal["simpson_1_3", "half_cycle", "none"] = Field(
        "simpson_1_3",
        description="Within-cycle-correction method. 'simpson_1_3' (DARTH default), "
        "'half_cycle' (trapezoidal), or 'none'.",
    )

    @field_validator("initial_distribution")
    @classmethod
    def _validate_initial_distribution(cls, v: List[float]) -> List[float]:
        arr = np.asarray(v, dtype=float)
        if (arr < -_PROB_TOL).any() or (arr > 1 + _PROB_TOL).any():
            raise ValueError("initial_distribution entries must lie in [0, 1]")
        if not np.isclose(arr.sum(), 1.0, atol=_PROB_TOL):
            raise ValueError(
                f"initial_distribution must sum to 1; got {float(arr.sum())}"
            )
        return v

    @model_validator(mode="after")
    def _validate_state_consistency(self) -> "CohortModel":
        n = len(self.states)
        if n < 2:
            raise ValueError("a cohort model needs at least 2 states")
        if len(self.initial_distribution) != n:
            raise ValueError(
                f"initial_distribution length {len(self.initial_distribution)} != n_states {n}"
            )
        for strat in self.strategies:
            if strat.n_states() != n:
                raise ValueError(
                    f"strategy '{strat.name}' has {strat.n_states()} states but model declares {n}"
                )
        names = [s.name for s in self.strategies]
        if len(set(names)) != len(names):
            raise ValueError(f"strategy names must be unique; got {names}")
        return self

    # -- I/O helpers ---------------------------------------------------------

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "CohortModel":
        """Load a model spec from a YAML file."""
        import yaml  # imported lazily so the core has no hard YAML runtime dep

        with open(path, "r") as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)
