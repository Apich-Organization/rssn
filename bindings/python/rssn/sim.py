"""Simulation scenarios by name, with JSON parameters.

`rssn_sim_run` takes a scenario name and a JSON parameter object; `{}` selects
the defaults of the scenarios that have them. The catalog in
:mod:`rssn.scenarios` lists every name and, for most of them, the parameters
its request type reads.
"""

from __future__ import annotations

import ctypes
import json
from typing import Any

from . import _ffi
from .errors import RssnError, Status
from .scenarios import SCENARIOS, SCENARIO_PARAMS

__all__ = ["SCENARIO_PARAMS", "SCENARIOS", "defaults", "names", "params", "run"]


def names() -> list[str]:
    """Every scenario name `run` accepts, in catalog order."""
    return [scenario["name"] for scenario in SCENARIOS]


def params(name: str) -> dict[str, dict[str, Any]]:
    """The parameters the scenario's request type reads (empty when unknown)."""
    return SCENARIO_PARAMS.get(name, {})


def defaults(name: str) -> dict[str, Any]:
    """The parameters with their defaults, where the source declares them.

    Only literals are recovered from the source; a field without a literal
    default is left out, so this is a starting point for `run`, not a complete
    schema.
    """
    return {
        field: spec["default"]
        for field, spec in params(name).items()
        if spec.get("default") is not None
    }


def run(name: str, parameters: dict[str, Any] | None = None) -> dict:
    """Run the scenario `name` and return its JSON result as a dict."""
    payload = "{}" if parameters is None else json.dumps(parameters)
    pointer = ctypes.c_void_p()
    code = _ffi.library().rssn_sim_run(
        name.encode("utf-8"), payload.encode("utf-8"), ctypes.byref(pointer)
    )
    if int(code) != int(Status.OK):
        raise RssnError(Status(int(code)), _ffi.last_error())
    return json.loads(_ffi.take_string(pointer.value))
