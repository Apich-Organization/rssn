"""Compare the Rust source change and C pass in a four-way experiment."""

import argparse
import json
import random
import statistics
import subprocess
from pathlib import Path

VARIANTS = ("old-off", "old-transform", "new-off", "new-transform")
COMPARISONS = (
    ("old-off", "old-transform"),
    ("old-off", "new-off"),
    ("new-off", "new-transform"),
    ("old-off", "new-transform"),
)


def measure(study, variant, mode, path, rounds):
    output = subprocess.check_output(
        [study / variant, mode, path, str(rounds)], text=True
    )
    elapsed, digest = output.split()
    return int(elapsed), digest


def run(study, mode, path, rounds, groups, rng):
    samples = []
    expected = None
    for group in range(groups):
        start = group % len(VARIANTS)
        order = VARIANTS[start:] + VARIANTS[:start]
        if (group // len(VARIANTS)) % 2:
            order = order[::-1]
        measurements = {}
        for variant in order:
            elapsed, digest = measure(study, variant, mode, path, rounds)
            if expected is None:
                expected = digest
            if digest != expected:
                raise RuntimeError(f"{variant} reads {path} as {digest}, expected {expected}")
            measurements[variant] = elapsed
        samples.append(measurements)

    effects = {}
    for base, other in COMPARISONS:
        differences = [100 * (row[other] / row[base] - 1) for row in samples]
        bootstrap = sorted(
            statistics.median(rng.choices(differences, k=groups)) for _ in range(5000)
        )
        effects[f"{base}/{other}"] = [
            statistics.median(differences),
            bootstrap[125],
            bootstrap[4874],
        ]
    return dict(rounds=rounds, checksum=expected, effects=effects, samples_ns=samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path, help="directory with the four binaries and generated inputs")
    parser.add_argument("--groups", type=int, default=21, help="number of matched sample groups")
    args = parser.parse_args()
    if args.groups < 1:
        parser.error("--groups must be positive")
    study = args.study.resolve()
    fixtures = Path(__file__).resolve().parents[4] / "tests/fixtures"
    inputs = (
        ("fixture_parquet", "parquet", fixtures / "parquet/codec_zstd.parquet", 650),
        ("fixture_hdf5", "hdf5", fixtures / "hdf5/filters_plugins.h5", 1800),
        ("many_pages", "parquet", study / "generated.parquet", 25),
        ("many_chunks", "hdf5", study / "generated.h5", 50),
    )
    rng = random.Random(314159)
    for name, mode, path, rounds in inputs:
        result = run(study, mode, path, rounds, args.groups, rng)
        print(json.dumps(dict(name=name, **result)), flush=True)


if __name__ == "__main__":
    main()
