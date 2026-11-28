"""Measure prebuilt LLVM research benchmarks."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
from run import measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpu", type=int, help="pin to one CPU in the current Linux affinity mask")
    modes = parser.add_subparsers(dest="mode", required=True)

    pair = modes.add_parser("pair", help="measure and compare two compatible executables")
    pair.add_argument("baseline", type=Path)
    pair.add_argument("experiment", type=Path)
    pair.add_argument("--pairs", type=int, default=11)
    pair.add_argument("--runner", help="run both executables through a program such as Wine")

    zstd = modes.add_parser("zstd-io", help="compare four complete I/O variants")
    zstd.add_argument("study", type=Path)
    zstd.add_argument("--groups", type=int, default=21)

    execute = modes.add_parser("run", help="run a harness that already measures its variants")
    execute.add_argument("binary", type=Path)
    execute.add_argument("--runner", help="run the executable through a program such as Wine")
    execute.add_argument("--env", action="append", default=[], metavar="NAME=VALUE")

    args = sys.argv[1:]
    separator = args.index("--") if "--" in args else len(args)
    options = parser.parse_args(args[:separator])
    workload = args[separator + 1:]
    if options.cpu is not None:
        if not hasattr(os, "sched_getaffinity") or options.cpu not in os.sched_getaffinity(0):
            parser.error("--cpu must be in the current Linux process affinity mask")
        os.sched_setaffinity(0, {options.cpu})
    if options.mode == "pair":
        if options.pairs < 3:
            parser.error("--pairs must be at least 3")
        result = measure(options.baseline.resolve(strict=True), options.experiment.resolve(strict=True),
                         workload, options.pairs, [options.runner] if options.runner else ())
        print(json.dumps(result))
    elif options.mode == "zstd-io":
        if workload:
            parser.error("zstd-io takes a study directory and --groups")
        runner = Path(__file__).resolve().parent / "zstd-io/compare.py"
        subprocess.run([sys.executable, runner, options.study, "--groups", str(options.groups)], check=True)
    else:
        path = options.binary.resolve(strict=True)
        runner = [options.runner] if options.runner else []
        env = os.environ.copy()
        for entry in options.env:
            name, equals, value = entry.partition("=")
            if not name or not equals:
                parser.error("--env expects NAME=VALUE")
            env[name] = value
        subprocess.run([*runner, path, *workload], env=env, check=True)


if __name__ == "__main__":
    main()
