#!/usr/bin/env python3
"""Run the cross-language conformance table against the Python binding.

    python bindings/tests/run_conformance.py            # verify cases.tsv
    python bindings/tests/run_conformance.py --record   # fill the expectations

`cases.tsv` holds one row per case and is shared with the Julia and C++ runners,
so the three wrappers are held to the same answers. `--record` runs the engine
and rewrites the expectation columns; verify afterwards.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import traceback

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(pathlib.Path(HERE).parents[0] / "python"))

import rssn  # noqa: E402


def read_cases(path: pathlib.Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip("\r")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        if len(fields) < len(header):
            fields += [""] * (len(header) - len(fields))
        rows.append(dict(zip(header, fields)))
    return rows


def write_cases(path: pathlib.Path, rows: list[dict[str, str]], header: list[str]) -> None:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.startswith("#")]
    lines.append("\t".join(header))
    for row in rows:
        lines.append("\t".join(row.get(column, "") for column in header))
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def config_for(row: dict[str, str]) -> rssn.Config | None:
    target = row.get("target", "").strip()
    bindings = row.get("bindings", "").strip()
    if target != "numeric" and not bindings:
        return None
    config = rssn.Config()
    if target == "numeric":
        config.numeric(float(row.get("tolerance") or 1e-12))
    for pair in filter(None, (entry.strip() for entry in bindings.split(";"))):
        name, _, value = pair.partition("=")
        config.bind(name.strip(), float(value))
    return config


def run_case(row: dict[str, str]) -> dict[str, str]:
    """Execute one case; returns the expectation columns it observed."""
    kind = row["kind"].strip()
    payload = row["payload"]
    result = {"expect_text": "", "expect_value": "", "expect_status": "ok"}

    if kind == "manifest":
        digest = rssn.manifest_hash()
        assert digest != 0, "manifest hash is zero"
        result["expect_value"] = str(digest)
        return result

    if kind == "units":
        result["expect_value"] = str(len(rssn.units.UNITS))
        return result

    if kind == "sim":
        name, _, parameters = payload.partition(" ")
        output = rssn.sim.run(name, json.loads(parameters or "{}"))
        assert isinstance(output, dict) and output, f"sim {name} returned no object"
        result["expect_text"] = ",".join(sorted(output))
        return result

    if kind == "callback":
        # A higher-order operator: the payload is `<op> <args>` and every runner
        # supplies the same residual `f(x) = x*x - 2`, so the row checks that the
        # three wrappers hand a host callback over the same way. The first field
        # of the answer record (the root) is the expectation.
        name, _, arguments = payload.partition(" ")
        values = [_literal(entry) for entry in arguments.split(",") if entry.strip()]
        answer = rssn.Session().call_value(name, *values, lambda x: x * x - 2.0)
        result["expect_text"] = f"{float(answer[0]):.10g}"
        return result

    session = rssn.Session()
    try:
        if kind == "apply":
            arguments = [float(entry) for entry in filter(None, row["bindings"].split(","))]
            result["expect_text"] = str(session.apply(payload, *arguments))
            return result
        if kind == "simplify":
            answer = session.simplify(payload)
        else:
            config = config_for(row)
            answer = session.compute(payload, config)
        result["expect_text"] = answer.text
        if answer.value is not None and row.get("target", "").strip() == "numeric":
            result["expect_value"] = repr(answer.value)
        return result
    finally:
        session.close()


def _literal(text: str) -> int | float:
    """A numeric literal from the table, keeping integers integral.

    An operator with an `Int` parameter refuses a `Float` argument, so `64` must
    not become `64.0` on the way in.
    """
    text = text.strip()
    if "." in text or "e" in text or "E" in text or text.lower() in ("inf", "nan"):
        return float(text)
    return int(text)


def expected_status_of(error: Exception) -> str:
    if isinstance(error, rssn.RssnError):
        return error.status.name.lower()
    return "exception"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true", help="rewrite the expectations")
    parser.add_argument("--cases", type=pathlib.Path, default=HERE / "cases.tsv")
    args = parser.parse_args(argv)

    rows = read_cases(args.cases)
    header = [
        "name", "kind", "payload", "target", "tolerance", "bindings",
        "expect_text", "expect_value", "expect_status",
    ]

    failures = 0
    for row in rows:
        name = row["name"]
        wanted_status = (row.get("expect_status") or "ok").strip()
        try:
            observed = run_case(row)
        except Exception as error:  # noqa: BLE001 - the status is the observation
            observed = {"expect_text": "", "expect_value": "", "expect_status": expected_status_of(error)}
            if args.record and not isinstance(error, rssn.RssnError):
                traceback.print_exc(limit=1)

        if args.record:
            row.update(observed)
            print(f"recorded {name}: {observed['expect_status']} {observed['expect_text'][:60]}")
            continue

        if observed["expect_status"] != wanted_status:
            print(f"FAIL {name}: status {observed['expect_status']} != {wanted_status}")
            failures += 1
            continue
        if wanted_status != "ok":
            print(f"ok   {name} ({wanted_status})")
            continue
        if observed["expect_text"] != row["expect_text"]:
            print(f"FAIL {name}: text {observed['expect_text']!r} != {row['expect_text']!r}")
            failures += 1
            continue
        if row["expect_value"] and observed["expect_value"]:
            want = float(row["expect_value"])
            got = float(observed["expect_value"])
            if abs(want - got) > 1e-9 * max(1.0, abs(want)):
                print(f"FAIL {name}: value {got!r} != {want!r}")
                failures += 1
                continue
        print(f"ok   {name}")

    if args.record:
        write_cases(args.cases, rows, header)
        print(f"wrote {args.cases}")
        return 0

    print(f"{len(rows) - failures}/{len(rows)} cases pass")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
