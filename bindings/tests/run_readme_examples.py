#!/usr/bin/env python3
"""Run every example the READMEs show, so the documentation cannot drift.

    python bindings/tests/run_readme_examples.py            # every language
    python bindings/tests/run_readme_examples.py python     # one of them

The READMEs are the front door of each binding, and an example that does not run
is worse than no example. This script extracts the fenced code blocks of a
README, evaluates them in order in one shared namespace (so a block may use what
an earlier one defined), and reports the blocks that fail.

A block that is *meant* to raise -- the strict-dimension conflict, a callback
that throws -- wraps the call in the language's own error handling, which is
what the README shows anyway.

C++ is the exception: a snippet is a fragment of a function, not a translation
unit. Its examples are compiled and run by `bindings/cpp/tests/workflow.cpp`,
which this script builds and runs instead.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BINDINGS = HERE.parent
ROOT = BINDINGS.parent

FENCE = re.compile(r"^```([A-Za-z+]*)\s*$")

#: The runner template for each language: a list of blocks, each evaluated at
#: top level so that state carries over and an error stays inside its block.
RUNNERS = {
    "python": """\
import sys, traceback
sys.path.insert(0, {python_path!r})   # the checkout's package, as the README assumes
blocks = {blocks!r}
failures = 0
for index, block in enumerate(blocks, 1):
    print(f"BLOCK {{index}}", flush=True)
    try:
        exec(compile(block, f"<block {{index}}>", "exec"), globals())
    except Exception as error:
        failures += 1
        print(f"FAIL {{index}}: {{type(error).__name__}}: {{error}}", flush=True)
print(f"{{len(blocks) - failures}}/{{len(blocks)}} blocks ran")
sys.exit(1 if failures else 0)
""",
    "julia": """\
push!(LOAD_PATH, {julia_path})        # the checkout's package, as the README assumes
blocks = {blocks}
failures = 0
for (index, block) in enumerate(blocks)
    println("BLOCK ", index)
    flush(stdout)
    try
        Base.include_string(Main, block)
    catch error
        global failures += 1
        println("FAIL ", index, ": ", sprint(showerror, error))
    end
end
println(length(blocks) - failures, "/", length(blocks), " blocks ran")
exit(failures == 0 ? 0 : 1)
""",
    "r": """\
grDevices::pdf(NULL)   # a plot in an example goes nowhere when this runs unattended
blocks <- {blocks}
failures <- 0L
for (index in seq_along(blocks)) {{
    cat("BLOCK", index, "\\n")
    flush(stdout())
    tryCatch({{
        eval(parse(text = blocks[[index]]), globalenv())
    }}, error = function(condition) {{
        failures <<- failures + 1L
        cat("FAIL", index, ":", conditionMessage(condition), "\\n")
    }})
}}
cat(length(blocks) - failures, "/", length(blocks), " blocks ran\\n")
quit(status = if (failures == 0L) 0L else 1L)
""",
}

LANGUAGES = {
    "python": {"readme": BINDINGS / "python" / "README.md", "fence": "python",
               "command": [sys.executable, "-u"], "suffix": ".py"},
    "julia": {"readme": BINDINGS / "julia" / "README.md", "fence": "julia",
              "command": ["julia"], "suffix": ".jl"},
    "r": {"readme": BINDINGS / "r" / "README.md", "fence": "r",
          "command": ["Rscript"], "suffix": ".R"},
}


def read_blocks(path: Path, fence: str) -> list[str]:
    """The fenced blocks of `fence` in `path`, in order."""
    blocks: list[str] = []
    current: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        match = FENCE.match(line)
        if match:
            if current is None:
                current = [] if match.group(1) == fence else None
                if match.group(1) == fence:
                    continue
            else:
                blocks.append("\n".join(current).strip("\n"))
                current = None
            continue
        if current is not None:
            current.append(line)
    if current:
        raise SystemExit(f"{path}: unterminated code fence")
    # A block that documents an installation step, or that shows a shape rather
    # than a program, says so with `# not run` on its first line.
    return [block for block in blocks
            if block.strip() and "# not run" not in block.splitlines()[0]]


def runner_for(language: str, blocks: list[str]) -> str:
    template = RUNNERS[language]
    if language == "python":
        return template.format(blocks=blocks, python_path=str(BINDINGS / "python"))
    if language == "julia":
        return template.format(
            blocks="[\n" + ",\n".join(_julia_string(block) for block in blocks) + "\n]",
            julia_path=json.dumps(str(BINDINGS / "julia" / "src")))
    return template.format(blocks="c(\n" + ",\n".join(_r_string(b) for b in blocks) + "\n)")


def _julia_string(text: str) -> str:
    """A Julia string literal for a whole block."""
    return '"""\n' + text.replace('"""', '\\"\\"\\"') + '\n"""'


def _r_string(text: str) -> str:
    """An R string literal for a whole block: `"` is escaped, `\\` is not."""
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'


def run_language(language: str, verbose: bool) -> int:
    spec = LANGUAGES[language]
    readme = spec["readme"]
    if not readme.exists():
        print(f"skip {language}: no {readme}")
        return 0
    if shutil.which(spec["command"][0]) is None:
        print(f"skip {language}: {spec['command'][0]} not on PATH")
        return 0

    blocks = read_blocks(readme, spec["fence"])
    if not blocks:
        print(f"skip {language}: no fenced examples")
        return 0

    with tempfile.TemporaryDirectory() as directory:
        script = Path(directory) / f"readme{spec['suffix']}"
        script.write_text(runner_for(language, blocks), encoding="utf-8")
        completed = subprocess.run(spec["command"] + [str(script)], capture_output=True,
                                   text=True, cwd=ROOT)
    output = completed.stdout + completed.stderr
    failures = [line for line in output.splitlines() if line.startswith("FAIL ")]
    summary = [line for line in output.splitlines() if "blocks ran" in line]
    print(f"== {language}: {readme.relative_to(ROOT)}")
    print(f"   {summary[0] if summary else 'no summary'}")
    for line in failures:
        print(f"   {line}")
    if verbose or failures:
        for line in output.splitlines():
            if not line.startswith("BLOCK "):
                print(f"     | {line}")
    return 1 if failures or completed.returncode != 0 else 0


def _import_library() -> Path | None:
    """The import library next to the shared library, as `run_all.sh` finds it.

    `RSSN_LIB` names the DLL (`librssn.so`, `librssn.dylib`); the import library
    that a link step needs sits beside it. Without `RSSN_LIB`, the build
    directories of a source checkout are tried in turn.
    """
    import os

    names = ("rssn.dll.lib", "librssn.dll.a", "rssn.lib", "librssn.a")
    environment = os.environ.get("RSSN_LIB", "")
    if environment:
        directory = Path(environment).parent
        for name in names:
            candidate = directory / name
            if candidate.exists():
                return candidate
        return None
    for directory in ("target/debug", "target-bindings/debug", "target/ffi",
                      "target/release"):
        for name in names:
            candidate = ROOT / directory / name
            if candidate.exists():
                return candidate
    return None


def run_cpp(verbose: bool) -> int:
    """The C++ examples are compiled, not extracted: build and run the tour."""
    build = ROOT / "build-bindings" / "cpp-workflow"
    library = _import_library()
    if library is None:
        print("skip cpp: no import library next to a built rssn (see bindings/README.md)")
        return 0
    configure = subprocess.run(
        ["cmake", "-S", str(BINDINGS / "cpp"), "-B", str(build), "-DRSSN_ROOT=" + str(ROOT),
         "-DRSSN_LIBRARY=" + str(library)],
        capture_output=True, text=True)
    if configure.returncode != 0:
        print("FAIL cpp: cmake configure failed")
        print(configure.stdout[-2000:], configure.stderr[-2000:])
        return 1
    compile_step = subprocess.run(["cmake", "--build", str(build), "--target", "rssn_workflow"],
                                  capture_output=True, text=True)
    if compile_step.returncode != 0:
        print("FAIL cpp: build failed")
        print(compile_step.stdout[-3000:], compile_step.stderr[-3000:])
        return 1
    executable = build / ("rssn_workflow.exe" if sys.platform == "win32" else "rssn_workflow")
    completed = subprocess.run([str(executable)], capture_output=True, text=True, cwd=ROOT)
    output = completed.stdout + completed.stderr
    print("== cpp: bindings/cpp/README.md (tests/workflow.cpp)")
    failures = [line for line in output.splitlines() if line.startswith("FAIL ")]
    print(f"   {'; '.join(failures) if failures else 'every example ran'}")
    if verbose or failures:
        for line in output.splitlines():
            print(f"     | {line}")
    return 1 if completed.returncode != 0 else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("languages", nargs="*", choices=sorted([*LANGUAGES, "cpp", "all"]),
                        default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)
    wanted = args.languages or ["all"]
    if "all" in wanted:
        wanted = sorted(LANGUAGES) + ["cpp"]

    failures = 0
    for language in wanted:
        if language == "cpp":
            failures += run_cpp(args.verbose)
        else:
            failures += run_language(language, args.verbose)
    print("every README example ran" if failures == 0 else f"{failures} language(s) failed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
