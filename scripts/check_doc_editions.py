#!/usr/bin/env python3
"""Check that the published RSSN guides have linked English and Chinese editions."""

from pathlib import Path
import re


def main() -> None:
    docs = Path(__file__).resolve().parent.parent / "docs"
    english = {path.name.removesuffix(".md"): path for path in docs.glob("*.md")
               if not path.name.endswith(".zh.md")}
    chinese = {path.name.removesuffix(".zh.md"): path for path in docs.glob("*.zh.md")}
    if not english or english.keys() != chinese.keys():
        raise SystemExit(f"unpaired RSSN guides: {sorted(english.keys() ^ chinese.keys())}")
    for name, original in sorted(english.items()):
        translation = chinese[name]
        for page, sibling in [(original, translation), (translation, original)]:
            content = page.read_text(encoding="utf-8")
            title = content.partition("\n")[0]
            if not title.startswith("# ") or f"]({sibling.name})" not in content:
                raise SystemExit(f"{page.name}: missing title or link to {sibling.name}")
            if page == translation and not re.search(r"[\u3400-\u9fff]", title):
                raise SystemExit(f"{translation.name}: missing Chinese title")
    print(f"checked {len(english)} paired RSSN documentation editions")


if __name__ == "__main__":
    main()
