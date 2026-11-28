#!/usr/bin/env python3
"""Verify that ATTRIBUTIONS.md is up to date with cargo-about."""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n")


def main() -> int:
    repo_root = pathlib.Path(__file__).resolve().parent.parent
    config = repo_root / "about.toml"
    template = repo_root / "about.hbs"
    target = repo_root / "ATTRIBUTIONS.md"

    if not target.is_file():
        print(
            f"ERROR: {target} is missing.\n"
            "Regenerate it with:\n"
            "  cargo about generate --config about.toml about.hbs -o ATTRIBUTIONS.md",
            file=sys.stderr,
        )
        return 1

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = pathlib.Path(tmpdir) / "ATTRIBUTIONS.md"
        cmd = [
            "cargo",
            "about",
            "generate",
            "--frozen",
            "--config",
            str(config),
            str(template),
            "-o",
            str(out_path),
        ]
        result = subprocess.run(cmd, cwd=repo_root, check=False)
        if result.returncode != 0:
            print(
                f"ERROR: cargo-about failed with exit code {result.returncode}",
                file=sys.stderr,
            )
            return result.returncode

        generated = normalize(out_path.read_text(encoding="utf-8"))

    if not generated.strip():
        print("ERROR: cargo-about produced empty output.", file=sys.stderr)
        return 1

    existing = normalize(target.read_text(encoding="utf-8"))
    if generated != existing:
        print(
            f"ERROR: {target} is stale and does not match cargo-about output.\n"
            "Regenerate it with:\n"
            "  cargo about generate --config about.toml about.hbs -o ATTRIBUTIONS.md",
            file=sys.stderr,
        )
        return 1

    print(f"OK: {target} is up to date ({len(generated.splitlines())} lines).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
