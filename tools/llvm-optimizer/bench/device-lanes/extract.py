"""Extract the RSSN CPU lane interpreter and its two device test workloads."""

import argparse
import re
from pathlib import Path

RSSN = Path(__file__).resolve().parents[4]
BACKEND = RSSN / "src/backend"


def item(source, declaration):
    match = re.search(
        r"^(?P<indent> *)(?:pub(?:\([^\n]*?\))? )?" + re.escape(declaration),
        source, re.M,
    )
    if match is None:
        raise ValueError(f"missing Rust item {declaration}")
    indent = match["indent"]
    end = source.index("\n" + indent + "}", match.end()) + len(indent) + 2
    return source[match.start():end]


def enum(source, name):
    return "#[derive(Clone, Copy, Debug, PartialEq)]\n" + item(source, "enum " + name + " {")


def test_program(source, name):
    original = item(source, f"fn {name}(")
    start = original.index("vec![")
    end = original.index("];", start) + 1
    ops = original[start:end].replace("Inst::", "KernelOp::")
    outputs = re.search(r"vec!\[([\d, ]+)\]\)\)\.unwrap\(\)", original)[1]
    return f"""fn {name}() -> KernelIr {{
    KernelIr {{ ops: {ops}, inputs: 2, outputs: vec![{outputs}] }}
}}"""


def extract(cpu_source):
    cpu = cpu_source.read_text()
    tape = (BACKEND / "tape.rs").read_text()
    kernel = (BACKEND / "device/kernel.rs").read_text()
    module = (BACKEND / "mod.rs").read_text()
    tests = (BACKEND / "device/tests.rs").read_text()
    lanes = re.search(r"pub const LANES: usize = \d+;", cpu)[0]
    helpers = cpu[cpu.index("#[inline(always)]\nfn unary("):
                  cpu.index("/// Executes the program for rows")]
    executor = item(cpu, "fn exec_block(")
    if "eval_powi_slice(" in executor:
        batch = item(tape, "fn eval_powi_slice(").replace("pub(super) ", "", 1)
        helpers += "#[inline]\n" + batch + "\n"
    return "\n\n".join([
        "type Reg = u32;",
        enum(tape, "Cmp"),
        enum(module, "Intrinsic"),
        enum(kernel, "KernelOp"),
        "pub struct KernelIr { ops: Vec<KernelOp>, inputs: usize, outputs: Vec<Reg> }",
        "impl Intrinsic {\n" + item(tape, "fn apply(") + "\n}",
        item(tape, "fn int_pow("),
        *(item(kernel, "fn " + name + "(") for name in
          ("select_is_true", "cmp_value", "gather_slot")),
        lanes,
        helpers,
        "#[inline(never)]\npub " + executor,
        test_program(tests, "everything"),
        test_program(tests, "smooth"),
    ]) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cpu-source", type=Path, default=BACKEND / "device/cpu.rs")
    args = parser.parse_args()
    args.output.write_text(extract(args.cpu_source))
