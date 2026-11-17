#!/usr/bin/env python3
"""Extract release notes for a given version from CHANGELOG.md."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

HEADING_RE = re.compile(r"^##\s+\[?v?([^\s\]]+)\]?(?:\s|$)")
LINK_DEF_RE = re.compile(r"^\[[^\]]+\]:\s+\S+")


def extract_section(changelog_text: str, version: str) -> str | None:
    target = version.lstrip("v")
    lines = changelog_text.splitlines()
    collecting = False
    collected: list[str] = []

    for line in lines:
        match = HEADING_RE.match(line)
        if match:
            heading_version = match.group(1)
            if collecting:
                break
            if heading_version.lower() == target.lower():
                collecting = True
                continue
        elif collecting:
            if LINK_DEF_RE.match(line):
                break
            collected.append(line)

    if not collecting:
        return None

    body = "\n".join(collected).strip()
    return body or None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract release notes for a version from CHANGELOG.md."
    )
    parser.add_argument(
        "--version",
        required=True,
        help="Target version (with or without leading 'v').",
    )
    parser.add_argument(
        "--changelog",
        default="CHANGELOG.md",
        help="Path to CHANGELOG.md (default: CHANGELOG.md).",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Optional output file path; prints to stdout if omitted.",
    )
    args = parser.parse_args()

    version = args.version.lstrip("v")
    changelog_path = Path(args.changelog)
    notes: str | None = None

    if changelog_path.is_file():
        text = changelog_path.read_text(encoding="utf-8")
        notes = extract_section(text, version)

    if not notes:
        notes = f"Release v{version}."

    output_text = f"{notes.rstrip()}\n"
    if args.output:
        Path(args.output).write_text(output_text, encoding="utf-8")
    else:
        sys.stdout.write(output_text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
