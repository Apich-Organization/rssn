"""Find the workflow spec and sources with a runfiles tree or manifest."""

import importlib.util
import json
import os
from pathlib import Path
import sys


def runfile_lookup(launcher):
    tree = Path(os.getenv("RUNFILES_DIR") or os.getenv("TEST_SRCDIR") or launcher + ".runfiles")
    manifest = Path(os.getenv("RUNFILES_MANIFEST_FILE") or launcher + ".runfiles_manifest")
    entries = {}
    if manifest.is_file():
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if line.startswith(" "):
                key, value = line[1:].split(" ", 1)
                unescape = lambda field: field.replace(r"\s", " ").replace(r"\n", "\n").replace(r"\b", "\\")
                entries[unescape(key)] = unescape(value)
            else:
                key, _, value = line.partition(" ")
                entries[key] = value or key

    def lookup(key):
        path = tree / key
        if path.exists():
            return path
        if key in entries:
            return Path(entries[key])
        # Some manifest-only runfiles trees list the files of a directory
        # artifact separately. Recover the artifact root from one member.
        prefix = key + "/"
        for member, target in entries.items():
            if member.startswith(prefix):
                depth = len(member[len(prefix):].split("/"))
                return Path(target).parents[depth - 1]
        raise FileNotFoundError(f"missing Bazel runfile {key}")

    return lookup


def main():
    _, _, spec_key, launcher, *args = sys.argv
    lookup = runfile_lookup(launcher)
    settings = json.loads(lookup(spec_key).read_text(encoding="utf-8"))
    runner = lookup(settings["runner"])
    sys.path.insert(0, str(runner.parent))
    module = importlib.util.spec_from_file_location("scientific_runner", runner)
    loaded = importlib.util.module_from_spec(module)
    module.loader.exec_module(loaded)
    loaded.main(settings, lookup, args)


if __name__ == "__main__":
    main()
