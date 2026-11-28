"""Extract the Tape column evaluator and the device workloads it can also run."""

import argparse
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
from extract import BACKEND, RSSN, enum, item


def extract(tape_source):
    tape = tape_source.read_text()
    backend = (BACKEND / "mod.rs").read_text()
    tests = (BACKEND / "device/tests.rs").read_text()
    ops = (RSSN / "src/graph/op.rs").read_text()
    parts = [
        "type Reg = u32;",
        re.search(r"^pub type EvalFn.*;", ops, re.M)[0],
        enum(tape, "Cmp"),
        enum(backend, "Intrinsic"),
        "#[derive(Copy, Clone)]\n" + item(tape, "enum Inst {"),
        item(tape, "struct Tape {"),
        item(tape, "trait Scalar:"),
        item(tape, "trait FromExtern:"),
        "impl Intrinsic {\n" + item(tape, "fn apply(") + "\n}",
        item(tape, "impl Scalar for f64 {"),
        item(tape, "impl FromExtern for f64 {"),
        item(tape, "fn int_pow("),
        re.search(r"pub const COLUMN_BLOCK.*;", tape)[0],
        item(tape, "const fn bit("),
        *(item(tape, "fn " + name) for name in
          ("compare_f64(", "gather_slot(", "with_operands<", "apply<T")),
        *("#[inline]\n" + item(tape, "fn " + name + "(") for name in
          ("eval_unary_slice", "eval_binary_slice", "eval_select_slice")),
    ]
    if "pub(super) fn eval_powi_slice" in tape:
        batch = item(tape, "fn eval_powi_slice(").replace("pub(super) ", "", 1)
        parts.append("#[inline]\n" + batch)
    column = item(tape, "fn eval_columns_range(")
    parts.append("impl Tape {\n    #[inline(never)]\n" + column + "\n}")
    for case in ("everything", "smooth"):
        original = item(tests, "fn " + case + "(")
        instructions = original[original.index("vec!["):original.index("];") + 1]
        values = re.search(
            r"tape\(insts, (\d+), (\d+), (vec!\[[\d, ]+\])\)", original,
        )
        parts.append(f"""fn {case}() -> Tape {{
    Tape {{ insts: {instructions}, args: vec![], inputs: {values[1]},
           params: {values[2]}, outputs: {values[3]} }}
}}""")
    return "\n\n".join(parts) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--tape-source", type=Path, default=BACKEND / "tape.rs")
    args = parser.parse_args()
    args.output.write_text(extract(args.tape_source))
