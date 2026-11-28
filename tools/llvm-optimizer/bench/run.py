"""Run paired benchmarks and count LLVM's annotated stack spills and reloads."""

import argparse
import collections
import json
import os
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.dont_write_bytecode = True
from callback_ir import kernel

BENCH = Path(__file__).resolve().parent
TARGET = "rssn-spill-rustfft-bench"
ANNOTATION = re.compile(r"# \d+-byte (?:Folded )?(Spill|Reload)\b")
FUNCTION = re.compile(r"# -- Begin function (\S+)")


def build(command, *, env=None, log=None):
    print("+", " ".join(map(str, command)), file=sys.stderr, flush=True)
    if log is None:
        subprocess.run(command, env=env, check=True)
    else:
        with log.open("w") as output:
            result = subprocess.run(command, env=env, stdout=output, stderr=output)
        if result.returncode:
            raise RuntimeError(f"Command failed. See {log}")


def stack_counts(asm, pattern=None):
    totals = collections.Counter({"Spill": 0, "Reload": 0})
    selected = pattern is None
    with asm.open() as stream:
        for line in stream:
            if pattern is not None:
                if function := FUNCTION.search(line):
                    selected = pattern in function[1]
                if "# -- End function" in line:
                    selected = False
            if selected and (event := ANNOTATION.search(line)):
                totals[event[1]] += 1
    return dict(totals)


def ir_build(args, root, suite, ir, runner, symbol):
    source = root / f"{suite}.ll"
    base = root / f"{suite}-off.ll"
    optimized = root / f"{suite}-transform.ll"
    tool_env = os.environ | args.llvm_env
    environment = dict(
        tool_env,
        RSSN_SPILL_MODE="off",
        RSSN_SPILL_FILTER="",
        RSSN_SPILL_ABI="off",
        RSSN_SPILL_PLAN_WIDTH=str(args.plan_width),
    )
    source.write_text(ir)
    build([args.opt, "-mcpu=native", "-passes=default<O3>", "-S", source, "-o", base], env=tool_env)
    build(
        [
            args.opt,
            "-mcpu=native",
            f"-load-pass-plugin={args.plugin}",
            "-passes=function(rssn-spill-pressure)",
            "-verify-each",
            "-pass-remarks-analysis=rssn-spill-pressure",
            "-S",
            base,
            "-o",
            optimized,
        ],
        env=environment,
        log=root / f"{suite}-remarks.log",
    )

    def emit(ir, mode):
        binary = root / f"{suite}-{mode}"
        asm = binary.with_suffix(".s")
        obj = binary.with_suffix(".o")
        opts = [args.llc, "-mcpu=native", "-O3", ir]
        build([*opts, "-o", asm], env=tool_env)
        build([*opts, "-filetype=obj", "-o", obj], env=tool_env)
        build([args.cc, "-O3", BENCH / runner, obj, "-o", binary])
        return binary, stack_counts(asm, symbol)

    return emit(base, "off"), emit(optimized, "transform")


def abi_build(args, root, mode):
    executable = root / f"abi-{mode}"
    asm, ir = executable.with_suffix(".s"), executable.with_suffix(".ll")
    env = dict(
        os.environ,
        RUSTC_BOOTSTRAP="1",
        RSSN_SPILL_MODE="off",
        RSSN_SPILL_FILTER="rssn_abi_mix",
        RSSN_SPILL_ABI=mode,
    )
    build(
        [
            args.rustc,
            "--edition=2024",
            "-Copt-level=3",
            "-Cpanic=abort",
            "-Clto=fat",
            "-Ccodegen-units=1",
            "-Ctarget-cpu=native",
            "-Cremark=rssn-call-abi",
            "-Cllvm-args=-pass-remarks-analysis=rssn-call-abi",
            "-Zverify-llvm-ir",
            "-Zverbose-asm",
            f"-Zllvm-plugins={args.plugin}",
            f"--emit=link={executable},asm={asm},llvm-ir={ir}",
            BENCH / "abi.rs",
        ],
        env=env,
        log=root / f"abi-{mode}.log",
    )
    with ir.open() as stream:
        kernel = next(line for line in stream if line.startswith("define ")
                      and "rssn_abi_mix" in line)
    if ("x86_regcallcc" in kernel) != (mode == "transform"):
        raise RuntimeError(f"ABI experiment did not produce the expected convention: {kernel}")
    return executable, stack_counts(asm, "rssn_abi_mix")


def fft_build(args, root, mode):
    target_dir = root / f"rustfft-{mode}"
    asm = root / f"rustfft-{mode}.s"
    flags = [
        "-Ctarget-cpu=native",
        f"-Zllvm-plugins={args.plugin}",
        "-Zverbose-asm",
        "-Zverify-llvm-ir",
        "-Cremark=rssn-spill-pressure",
        "-Cllvm-args=-pass-remarks-analysis=rssn-spill-pressure",
    ]
    env = dict(
        os.environ,
        RUSTC_BOOTSTRAP="1",
        RUSTC=args.rustc,
        RSSN_SPILL_MODE=mode,
        RSSN_SPILL_FILTER="",
        RSSN_SPILL_ABI="off",
        RSSN_SPILL_PLAN_WIDTH=str(args.plan_width),
        CARGO_TARGET_DIR=str(target_dir),
        CARGO_BUILD_BUILD_DIR=str(target_dir),
        CARGO_BUILD_JOBS=str(args.jobs),
        CARGO_ENCODED_RUSTFLAGS="\x1f".join(flags),
    )
    command = [
        args.cargo,
        "rustc",
        "--release",
        "--locked",
        "--manifest-path",
        BENCH / "Cargo.toml",
    ]
    if args.offline:
        command.append("--offline")
    # Cargo may store intermediate outputs in a separate build directory.
    build([*command, "--", f"--emit=asm={asm},link"], env=env, log=root / f"rustfft-{mode}.log")
    return target_dir / "release" / TARGET, stack_counts(asm, "7rustfft")


def probe(binary, arguments):
    nanos, digest = subprocess.check_output(
        [binary, *map(str, arguments)], text=True
    ).split()
    return int(nanos), digest


def measure(baseline, experiment, arguments, pairs):
    warm_base, reference = probe(baseline, arguments)
    warm_experiment, checksum = probe(experiment, arguments)
    if checksum != reference or min(warm_base, warm_experiment) <= 0:
        raise RuntimeError("warm-up outputs differ or the sample is too short")
    timings = []
    for i in range(pairs):
        binaries = [baseline, experiment] if i % 2 == 0 else [experiment, baseline]
        observations = [probe(binary, arguments) for binary in binaries]
        if {checksum for _, checksum in observations} != {reference}:
            raise RuntimeError(f"bitwise output mismatch in pair {i + 1}")
        a, b = [elapsed for elapsed, _ in observations]
        timings.append((a, b) if i % 2 == 0 else (b, a))

    ratios = [b / a for a, b in timings]
    rng = random.Random(185)
    boot = sorted(
        statistics.median(rng.choices(ratios, k=pairs)) for _ in range(4000)
    )
    return {
        "median_base_ns": statistics.median(a for a, _ in timings),
        "median_pass_ns": statistics.median(b for _, b in timings),
        "paired_runtime_change_pct": 100 * (statistics.median(ratios) - 1),
        "bootstrap_95_pct": [100 * (boot[i] - 1) for i in (100, 3899)],
        "checksum": reference,
        "samples_ns": timings,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plugin", type=Path)
    parser.add_argument("--suite", choices=("both", "all", "callback", "lanes", "abi", "rustfft"),
                        default="all")
    parser.add_argument("--work-dir", type=Path, default=BENCH.parent / "build-bench")
    parser.add_argument("--opt", default="opt-23")
    parser.add_argument("--llc", default="llc-23")
    parser.add_argument("--llvm-env", action="append", default=[], help=argparse.SUPPRESS)
    parser.add_argument("--rustc", default="rustc")
    parser.add_argument("--cargo", default="cargo")
    parser.add_argument("--cc", default="cc")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    cpus = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else ()
    parser.add_argument("--cpu", type=int, default=min(cpus) if cpus else None)
    parser.add_argument("--pairs", type=int, default=11)
    parser.add_argument("--plan-width", type=int, choices=range(1, 5), default=4)
    parser.add_argument("--fft-rounds", type=int, default=24000)
    parser.add_argument("--callback-rounds", type=int, default=4000000)
    parser.add_argument("--abi-rounds", type=int, default=60000000)
    args = parser.parse_args()
    args.llvm_env = dict(setting.split("=", 1) for setting in args.llvm_env)
    if min(args.jobs, args.fft_rounds, args.callback_rounds, args.abi_rounds) < 1 or args.pairs < 3:
        parser.error("positive jobs and rounds and at least 3 pairs are required")
    if args.cpu is not None and args.cpu not in cpus:
        parser.error("--cpu must be in the current process affinity mask")
    args.plugin = args.plugin.resolve(strict=True)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="ab-", dir=args.work_dir.resolve()))
    print(f"Artifacts in {root}", file=sys.stderr, flush=True)

    setups = []
    stacks = {}
    for suite in ("callback", "lanes", "abi", "rustfft"):
        if args.suite not in ("all", suite) and not (
            args.suite == "both" and suite in ("callback", "rustfft")
        ):
            continue
        start = time.monotonic()
        if suite == "callback":
            baseline, experiment = ir_build(
                args, root, suite, kernel(), "callback_runner.c", "rssn_callback_kernel"
            )
            cases = [(suite, [args.callback_rounds])]
        elif suite == "lanes":
            baseline, experiment = ir_build(
                args, root, suite, (BENCH.parent / "tests/fusion.ll").read_text(),
                "lane_runner.c", "pack_live_scalars",
            )
            cases = [(suite, [args.callback_rounds])]
        elif suite == "abi":
            baseline = abi_build(args, root, "off")
            experiment = abi_build(args, root, "transform")
            cases = [(suite, [args.abi_rounds])]
        else:
            baseline = fft_build(args, root, "off")
            experiment = fft_build(args, root, "transform")
            sizes = (12, 18, 27, 64, 120, 256, 512, 704, 4096)
            cases = [(f"rustfft-{n}", [n, args.fft_rounds]) for n in sizes]
        stacks[suite] = {
            "off": baseline[1],
            "transform": experiment[1],
            "total_build_seconds": time.monotonic() - start,
        }
        setups.extend((name, baseline[0], experiment[0], opts) for name, opts in cases)

    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    results = {
        "cpu": args.cpu,
        "plan_width": args.plan_width,
        "rustc": subprocess.check_output([args.rustc, "-vV"], text=True),
        "stack_annotations": stacks,
        "timings": {},
    }
    for name, baseline, experiment, opts in setups:
        results["timings"][name] = measure(baseline, experiment, opts, args.pairs)
    (root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    for name, result in results["timings"].items():
        low, high = result["bootstrap_95_pct"]
        change = result["paired_runtime_change_pct"]
        print(f"{name:<14} {change:+.2f}% runtime [{low:+.2f}%, {high:+.2f}%]")
    print(root / "results.json")


if __name__ == "__main__":
    main()
