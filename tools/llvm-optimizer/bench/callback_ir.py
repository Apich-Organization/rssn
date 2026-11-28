"""Generate a late-IR SIMD kernel with independent data live across an ABI call."""

import argparse

LANES = 4


def kernel(vectors=12):
    ir = [
        'target triple = "x86_64-unknown-linux-gnu"',
        "declare double @rssn_gain(i64) nounwind nofree nosync willreturn memory(none)",
        "",
        "define void @rssn_callback_kernel(ptr noundef noalias readonly %input,",
        "                                  ptr noundef noalias writeonly %output, i64 %key) {",
        "entry:",
    ]
    for i in range(vectors):
        ir += [
            f"  %p{i} = getelementptr double, ptr %input, i64 {i * LANES}",
            f"  %v{i} = load <4 x double>, ptr %p{i}, align 8",
        ]
    ir += [
        "  %gain = call double @rssn_gain(i64 %key)",
        "  %lane = insertelement <4 x double> poison, double %gain, i64 0",
        "  %scale = shufflevector <4 x double> %lane, <4 x double> poison, <4 x i32> zeroinitializer",
    ]
    for i in reversed(range(vectors)):
        ir += [
            f"  %product{i} = fmul <4 x double> %v{i}, %scale",
            f"  %q{i} = getelementptr double, ptr %output, i64 {i * LANES}",
            f"  store <4 x double> %product{i}, ptr %q{i}, align 8",
        ]
    return "\n".join(ir + ["  ret void", "}", ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vectors", type=int, default=12, choices=range(4, 17))
    print(kernel(parser.parse_args().vectors), end="")
