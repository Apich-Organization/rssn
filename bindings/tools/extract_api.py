#!/usr/bin/env python3
"""Extract the RSSN embedding surface into one machine-readable file.

The language bindings in this directory are generated from `api/rssn_api.json`,
which this script derives from the engine's own kernel manifest plus a source
scan for the simulation scenarios (the manifest does not list them and the C
ABI has no entry point that enumerates them).

Usage, from the repository root:

    python bindings/tools/extract_api.py
    python bindings/tools/extract_api.py --check   # exit 1 if the file is stale

The manifest defaults to the explicit snapshot in `bindings/api/manifest.json`
or to the output of `bazel build //:manifest`.

`--manifest` overrides the input path. Use `--library` to read the JSON
manifest from a built RSSN shared library directly. `--write-manifest`
saves the same input JSON for tools such as the Aevia API audit.
"""

from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
BINDINGS = ROOT / "bindings"
API_DIR = BINDINGS / "api"
API_JSON = API_DIR / "rssn_api.json"


# --------------------------------------------------------------------------
# dimension rendering, mirroring `kernel::signature::dim_name` but in the text
# grammar `Graph::dimension_from_text` accepts (`A*B^2/C`, `base:Name`).
# --------------------------------------------------------------------------

def _exponent_text(exponent: Any) -> str:
    if isinstance(exponent, dict):
        if "rat" in exponent:
            rat = exponent["rat"]
            num, den = rat["numer"], rat["denom"]
            return str(num) if den == 1 else f"{num}/{den}"
        if "int" in exponent:
            return str(exponent["int"])
    return str(exponent)


def _needs_parens(text: str) -> bool:
    return "/" in text or text.startswith("-")


def render_dimension(dim: Any) -> str:
    """A `DimExpr` as the engine's dimension text."""
    if dim is None:
        return "1"
    if isinstance(dim, str):
        return dim
    if "any" in dim:
        return "any"
    if "non_metric" in dim:
        return "non-metric"
    if "container" in dim:
        return "container"
    if "dimless" in dim:
        return "1"
    if "base" in dim:
        return str(dim["base"])
    if "var" in dim:
        return "?" + str(dim["var"])
    if "inverse" in dim:
        return f"1/({render_dimension(dim['inverse'])})"
    if "pow" in dim:
        base, exponent = dim["pow"]
        text = _exponent_text(exponent)
        base_text = render_dimension(base)
        if "*" in base_text:
            base_text = f"({base_text})"
        return f"{base_text}^{text}"
    if "mul" in dim:
        factors = [render_dimension(factor) for factor in dim["mul"]]
        return "*".join(factors)
    raise ValueError(f"unknown dimension expression: {dim!r}")


def _rational_text(rat: Any) -> tuple[str, str]:
    """A bignum rational encoded as `[sign, limbs]` as (numerator, denominator)."""
    def decode(part: Any) -> str:
        sign, limbs = part[0], part[1]
        value = 0
        for index, limb in enumerate(limbs):
            value += int(limb) << (32 * index)
        return str(-value if sign < 0 else value)

    if rat is None:
        return "1", "1"
    return decode(rat["numer"]), decode(rat["denom"])


def scale_record(scale: Any) -> dict[str, Any]:
    if scale is None:
        return {"numer": "1", "denom": "1", "pi": 0}
    numer, denom = _rational_text(scale)
    return {"numer": numer, "denom": denom, "pi": int(scale.get("pi", 0))}


# --------------------------------------------------------------------------
# scenario catalog: a source scan, because the manifest has no scenario list
# --------------------------------------------------------------------------

def _balanced(text: str, start: int) -> str:
    """The `{...}` block that starts at or after `start`."""
    open_at = text.index("{", start)
    depth = 0
    for index in range(open_at, len(text)):
        char = text[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return text[open_at + 1 : index]
    return ""


_RUST_TYPES = {
    "usize": "int",
    "u32": "int",
    "u64": "int",
    "i32": "int",
    "i64": "int",
    "f64": "float",
    "f32": "float",
    "bool": "bool",
    "String": "str",
}


def _json_type(rust_type: str) -> str:
    text = rust_type.strip()
    if text in _RUST_TYPES:
        return _RUST_TYPES[text]
    if text.startswith("Vec<"):
        inner = text[4:-1]
        return f"list[{_json_type(inner)}]"
    if text.startswith("Option<"):
        return _json_type(text[7:-1]) + "|null"
    if text.startswith("Array2<") or text.startswith("Array1<"):
        return "nested list"
    return text


def _literal(text: str) -> Any:
    text = text.strip()
    if text in ("true", "false"):
        return text == "true"
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    if re.fullmatch(r"-?(\d+\.\d*|\.\d+|\d+)([eE][-+]?\d+)?", text):
        return float(text)
    quoted = re.fullmatch(r'"([^"]*)"(?:\.into\(\)|\.to_owned\(\)|\.to_string\(\))?', text)
    if quoted:
        return quoted.group(1)
    if text.startswith("vec!["):
        return None
    if text in ("None",):
        return None
    if text == "true":
        return True
    return None


def _struct_fields(source: str, name: str) -> dict[str, dict[str, Any]]:
    match = re.search(rf"struct\s+{name}\s*\{{", source)
    if not match:
        return {}
    body = _balanced(source, match.start())
    fields: dict[str, dict[str, Any]] = {}
    for line in body.splitlines():
        line = line.split("//")[0].strip().rstrip(",")
        if not line or line.startswith("#[") or line.startswith("///"):
            continue
        field = re.fullmatch(r"(?:pub\s+)?([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(.+)", line)
        if not field:
            continue
        fields[field.group(1)] = {"type": _json_type(field.group(2)), "default": None}
    return fields


def _struct_defaults(source: str, name: str) -> dict[str, Any]:
    match = re.search(rf"impl\s+Default\s+for\s+{name}\s*\{{", source)
    if not match:
        return {}
    body = _balanced(source, match.start())
    defaults: dict[str, Any] = {}
    for line in body.splitlines():
        line = line.split("//")[0].strip().rstrip(",")
        field = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(.+)", line)
        if not field:
            continue
        value = _literal(field.group(2))
        if value is not None:
            defaults[field.group(1)] = value
    return defaults


def _scenario_entry(
    name: str, source: str, request: str | None, group: str, note: str = ""
) -> dict[str, Any]:
    params: dict[str, dict[str, Any]] = {}
    if request:
        params = _struct_fields(source, request)
        for field, value in _struct_defaults(source, request).items():
            if field in params:
                params[field]["default"] = value
    return {"name": name, "group": group, "params": params, "note": note}


def extract_scenarios() -> list[dict[str, Any]]:
    scenarios: list[dict[str, Any]] = []
    fixed_path = ROOT / "src" / "sim" / "scenario.rs"
    engine_path = ROOT / "src" / "sim" / "scenario_engines.rs"
    fixed = fixed_path.read_text(encoding="utf-8") if fixed_path.exists() else ""
    engines = engine_path.read_text(encoding="utf-8") if engine_path.exists() else ""

    # Fixed scenarios: the name list is `const FIXED_SCENARIOS: &[&str] = &[...]`,
    # and each request struct is `struct XxxRequest { .. }` next to its handler.
    fixed_names = re.search(r"const\s+FIXED_SCENARIOS\s*:[^=]*=\s*&\[(.*?)\];", fixed, re.S)
    handlers: dict[str, str] = {}
    for match in re.finditer(r"struct\s+(\w+Request)\s*\{", fixed):
        body = _balanced(fixed, match.start())
        flattened = re.search(r"#\[serde\(flatten\)\]\s*pub\s+(\w+)\s*:\s*([\w:]+)", body)
        if flattened:
            handlers[match.group(1)] = ("flatten", flattened.group(2))
        else:
            handlers[match.group(1)] = ("fields", "")

    # The dispatch arm says which request type a fixed scenario parses:
    # `| "ising" => { let request: IsingRequest = parse(params)?; ... }`.
    arms: dict[str, str] = {}
    for match in re.finditer(r'\|\s*"([^"]+)"\s*=>\s*\{', fixed):
        body = _balanced(fixed, match.start())
        request = re.search(r":\s*(\w+)\s*=\s*parse\(", body)
        if request:
            arms[match.group(1)] = request.group(1)

    if fixed_names:
        for raw in re.findall(r'"([^"]+)"', fixed_names.group(1)):
            scenarios.append(_scenario_entry(raw, fixed, arms.get(raw), "fixed"))

    # Engine scenarios: `engine_scenarios! { name => Req, name, ... }`.
    block = re.search(r"engine_scenarios!\s*\{(.*?)\n\}", engines, re.S)
    if block:
        for line in block.group(1).splitlines():
            line = line.split("//")[0].strip().rstrip(",")
            if not line:
                continue
            parts = [part.strip() for part in line.split("=>")]
            name = parts[0]
            request = parts[1] if len(parts) > 1 else None
            scenarios.append(_scenario_entry(name, engines, request, "engine"))
    return scenarios


# --------------------------------------------------------------------------

def read_library_manifest(path: Path) -> str:
    library = ctypes.CDLL(str(path.resolve()))
    manifest_json = library.rssn_manifest_json
    manifest_json.argtypes = []
    manifest_json.restype = ctypes.c_void_p
    release = library.rssn_string_free
    release.argtypes = [ctypes.c_void_p]
    release.restype = None
    pointer = manifest_json()
    if pointer is None:
        raise ValueError("the RSSN library returned no manifest")
    try:
        return ctypes.string_at(pointer).decode("utf-8")
    finally:
        release(pointer)


def build(manifest: dict[str, Any]) -> dict[str, Any]:
    ops = []
    for op in manifest["ops"]:
        arity = op["arity"]
        ops.append(
            {
                "name": op["name"],
                "rule_set": op["rule_set"],
                "category": op["category"] or "",
                "arity_min": int(arity["min"]),
                "arity_max": None if arity["max"] is None else int(arity["max"]),
                "doc": (op.get("doc") or "").strip(),
                "ffi": bool(op["backend_support"]["ffi"]),
                "eval": op["eval"],
                "keys": op["keys"],
                "overloads": op["overloads"],
                "overload_sources": op["overload_sources"],
                "signature_source": op["signature_source"],
                "dim_rule": op["dim_rule"],
                "dynamic_reason": op["dynamic_reason"],
                "binder": op["binder"],
                "aliases": list(op.get("aliases") or []),
                "since": op.get("since") or "",
                "stability": op.get("stability") or "",
            }
        )
    ops.sort(key=lambda entry: entry["name"])

    units = []
    for unit in manifest["units"]:
        units.append(
            {
                "name": unit["name"],
                "display": unit["display"],
                "dimension": render_dimension(unit["dimension"]),
                "scale": scale_record(unit.get("scale")),
                "offset": scale_record(unit["offset"]) if unit.get("offset") else None,
                "affine": bool(unit["affine"]),
                "aliases": list(unit.get("aliases") or []),
                "prefixable": bool(unit.get("prefixable")),
                "category": unit.get("category") or "",
                "since": unit.get("since") or "",
            }
        )

    constants = []
    for constant in manifest["constants"]:
        numer, denom = _rational_text(constant["value"])
        uncertainty = constant.get("uncertainty")
        if uncertainty is not None:
            uncertainty = "/".join(_rational_text(uncertainty))
        constants.append(
            {
                "name": constant["name"],
                "symbol": constant["symbol"],
                "numer": numer,
                "denom": denom,
                "pi": int(constant["value"].get("pi", 0)),
                "dimension": render_dimension(constant["dimension"]),
                "unit": constant.get("unit") or "",
                "exact": bool(constant.get("exact")),
                "uncertainty": uncertainty,
                "source": constant.get("source") or "",
            }
        )

    return {
        "generator": "bindings/tools/extract_api.py",
        "rssn_version": manifest["rssn_version"],
        "api_version": manifest["api_version"],
        "manifest_hash": bytes(manifest["manifest_hash"]).hex(),
        "features": manifest["features"],
        "rule_sets": [entry["name"] for entry in manifest["rule_sets"]],
        "base_dimensions": [entry["name"] for entry in manifest["base_dimensions"]],
        "value_kinds": manifest.get("value_kinds", []),
        "unit_systems": [entry["name"] for entry in manifest.get("unit_systems", [])],
        "ops": ops,
        "units": units,
        "constants": constants,
        "scenarios": extract_scenarios(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--manifest", type=Path, help="read the kernel manifest JSON")
    source.add_argument("--library", type=Path, help="read the manifest from a built RSSN shared library")
    parser.add_argument("--write-manifest", type=Path, help="save the source manifest JSON for other tools")
    parser.add_argument("--output", type=Path, default=API_JSON, help="write the binding API description here")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)

    manifest_path = args.manifest
    if manifest_path is None and args.library is None:
        candidates = (API_DIR / "manifest.json", ROOT / "bazel-bin" / "manifest" / "manifest.json")
        manifest_path = next((candidate for candidate in candidates if candidate.is_file()), None)
        if manifest_path is None:
            parser.error("no kernel manifest found; use --manifest or --library")
    try:
        if args.library is None:
            manifest_text = manifest_path.read_text(encoding="utf-8")
        else:
            manifest_text = read_library_manifest(args.library)
        payload = build(json.loads(manifest_text))
    except (OSError, AttributeError, ValueError, TypeError, KeyError) as error:
        parser.error(f"cannot read the RSSN manifest: {error}")
    if args.write_manifest is not None:
        args.write_manifest.parent.mkdir(parents=True, exist_ok=True)
        args.write_manifest.write_text(manifest_text.rstrip() + "\n", encoding="utf-8")
    output = args.output
    text = json.dumps(payload, indent=1, sort_keys=False, ensure_ascii=False) + "\n"

    if args.check:
        current = output.read_text(encoding="utf-8") if output.exists() else ""
        if current != text:
            print(f"error: {output} is stale; rerun extract_api.py", file=sys.stderr)
            return 1
        print(f"{output} is current ({len(payload['ops'])} ops)")
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding="utf-8")
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]
    print(
        f"wrote {output} ({len(payload['ops'])} ops, {len(payload['units'])} units, "
        f"{len(payload['constants'])} constants, {len(payload['scenarios'])} scenarios, "
        f"sha256:{digest})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
