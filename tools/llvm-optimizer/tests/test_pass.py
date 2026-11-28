import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(TESTS.parent / "bench"))
from callback_ir import kernel


def run(args, *, input=None, env=None):
    result = subprocess.run(args, input=input, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"{args}\n{result.stderr}")
    return result.stdout, result.stderr


def body(ir, name):
    pattern = r"^define [^\n]*?@" + re.escape(name) + r"\([^{}]*\{(.*?)^}"
    return re.search(pattern, ir, re.M | re.S)[1]


def shifted_tail(ir):
    function = body(ir, "tail_zero_numerator")
    for suffix in ("", ".r"):
        function = function.replace(
            f"  %score{suffix} = fsub double %product{suffix}, %x\n",
            f"  %decay{suffix} = fmul double %scale, %x\n"
            f"  %score{suffix} = fsub double %product{suffix}, %decay{suffix}\n"
            f"  %adjusted{suffix} = fsub double %score{suffix}, %shift\n",
        ).replace(
            f"@llvm.exp.f64(double %score{suffix})", f"@llvm.exp.f64(double %adjusted{suffix})",
        )
    signature = ("define double @tail_shifted(double %seed, double %c, double %den, "
                 "double %scale, double %shift, i32 %rounds)")
    return ir + "\n" + signature + " {" + function + "}\n"


def hoist_unit_base(ir):
    ir = ir.replace("left:\n  %b0 = fsub double 1.0, %t\n", "left:\n")
    return ir.replace("  %odd = trunc", "  %b0 = fsub double 1.0, %t\n  %odd = trunc")


def ir_tests(opt, llc, plugin):
    env = dict(
        os.environ, RSSN_SPILL_MODE="off", RSSN_SPILL_FILTER="",
        RSSN_SPILL_PLAN_WIDTH="4", RSSN_SPILL_ABI="off", RSSN_SPILL_LIBM="off",
        RSSN_SPILL_LIBM_REUSE="off",
        RSSN_LIBM_UNIT_LOG="off",
        RSSN_LIBM_TAIL_CUTOFF="off",
        RSSN_POW_EXPONENT_PACK="off",
        RSSN_LOOP_EXIT_COMPRESSION="off",
        RSSN_LOOP_PREHEADER_PRESSURE="off",
        RSSN_FP_ENDPOINT_GUARDS="off",
        RSSN_FP_CLASS_INTERVAL="off",
        RSSN_SCALAR_LOAD_PACK="off",
        RSSN_SWITCH_LOAD_REUSE="off",
        RSSN_FP_REGISTER_BANK="off",
        RSSN_POW_UNIT_BASE="off",
    )

    def transform(ir, passes="rssn-spill-pressure", *, settings=env, cpu="haswell",
                  extra=()):
        return run(
            [
                opt,
                f"-mcpu={cpu}",
                *extra,
                f"-load-pass-plugin={plugin}",
                f"-passes={passes}"
                if "<" in passes or passes.startswith("rssn-call-abi")
                or passes in ("rssn-release", "verify")
                else f"-passes=function({passes})",
                "-verify-each",
                "-pass-remarks-analysis=rssn-.*",
                "-S",
                "-o",
                "-",
                "-",
            ],
            input=ir,
            env=settings,
        )

    memory, remarks = transform((TESTS / "memory.ll").read_text())
    assert body(memory, "merge_disjoint").count("load double") == 1, memory
    for name in ("alias", "atomic", "volatile", "lifetime", "call_pressure", "unknown"):
        assert body(memory, f"keep_{name}").count("load double") == 2, name
    assert body(memory, "keep_metadata").count("load ptr") == 2, memory
    assert "merged loads 1" in remarks, remarks

    cfg, remarks = transform((TESTS / "cfg.ll").read_text())
    assert body(cfg, "split_defined").count("fneg") == 2, cfg
    assert body(cfg, "split_defined").index("call i64 @gate") < body(
        cfg, "split_defined"
    ).index("fneg"), cfg
    assert "extra copies 1" in remarks, remarks
    assert body(cfg, "keep_undef_operand").count("fneg") == 1, cfg
    assert body(cfg, "keep_phi_edge").index("fneg") < body(
        cfg, "keep_phi_edge"
    ).index("call i64 @gate"), cfg
    assert "SIMD live through 1 call sites 2 -> 2" in next(
        line for line in remarks.splitlines() if "keep_phi_edge " in line
    ), remarks
    moved = body(cfg, "sink_single_successor")
    assert moved.count("fneg") == 1 and moved.index("left:") < moved.index("fneg"), moved

    fusion, remarks = transform((TESTS / "fusion.ll").read_text())
    packed = body(fusion, "pack_contiguous")
    assert packed.count("load <2 x double>") == 1 and "!noundef" in packed, packed
    halves = body(fusion, "concat_halves")
    assert halves.count("load <4 x double>") == 1 and "shufflevector" not in halves, halves
    for name in ("gap", "swapped", "clobbered", "volatile_lane", "lanes_clobbered",
                 "mixed_scalar_ops"):
        assert "load <2 x double>" not in body(fusion, f"keep_{name}"), name
    assert body(fusion, "merge_equivalent_addresses").count("load double") == 1
    # An IR carrier with only extracts can disappear in SelectionDAG.
    # Verify the ABI spill pair on target assembly as well as in the IR.
    original = (TESTS / "fusion.ll").read_text()
    for ir, expected in ((original, 2), (fusion, 1)):
        asm, _ = run([llc, "-O3", "-mcpu=haswell", "-o", "-", "-"], input=ir)
        fn = re.search(r"^pack_live_scalars:.*?^\.Lfunc_end", asm, re.M | re.S)[0]
        assert fn.count("byte Spill") == expected, fn
        assert fn.count("Reload") == expected, fn
    carrier = body(fusion, "pack_live_scalars")
    assert carrier.count("load <2 x double>") == 1, carrier
    assert carrier.index("load <2 x double>") < carrier.index("call double"), carrier
    assert carrier[carrier.index("call double"):].count("extractelement") == 2, carrier
    assert "fmul <2 x double>" in carrier, carrier
    assert "lane groups 1, source loads fused 2" in remarks, remarks
    narrow, _ = transform((TESTS / "fusion.ll").read_text(), cpu="core2")
    assert body(narrow, "concat_halves").count("load <2 x double>") == 2, narrow

    joint_ir = (TESTS / "joint.ll").read_text()
    grouped, remarks = transform(joint_ir)
    shared = body(grouped, "shared_joint")
    assert shared[:shared.index("call double")].count("load double") == 0, shared
    assert "joint plans 1, joint members 3" in remarks, remarks
    single, _ = transform(joint_ir, settings=dict(env, RSSN_SPILL_PLAN_WIDTH="1"))
    shared = body(single, "shared_joint")
    assert shared[:shared.index("call double")].count("load double") == 2, shared
    for ir, expected in ((single, 2), (grouped, 0)):
        asm, _ = run([llc, "-O3", "-mcpu=haswell", "-o", "-", "-"], input=ir)
        fn = re.search(r"^shared_joint:.*?^\.Lfunc_end", asm, re.M | re.S)[0]
        assert fn.count("byte Spill") == expected, fn
        assert fn.count("Reload") == expected, fn
    _, avx = transform(joint_ir, "rssn-spill-pressure-report")
    _, sse = transform(joint_ir, "rssn-spill-pressure-report", cpu="core2")
    assert "1 call sites 1 -> 1" in next(
        line for line in avx.splitlines() if "width_probe " in line
    ), avx
    assert "1 call sites 2 -> 2" in next(
        line for line in sse.splitlines() if "width_probe " in line
    ), sse

    frames_ir = (TESTS / "private-frame.ll").read_text()
    frames, remarks = transform(frames_ir)
    for name in ("private_frame", "private_frame_gep"):
        fn = body(frames, name)
        assert fn.index("call double @observer") < fn.index("load double"), fn
    for name in ("unknown_capture", "escaped_frame", "frame_visible_to_call",
                 "expired_frame", "frame_exposed_implicitly"):
        fn = body(frames, name)
        assert fn.index("load double") < fn.index("call double"), fn
    assert "private call boundaries 1" in remarks, remarks
    # localescape is an implicit frame escape; remove that negative fixture
    # from codegen because it needs the platform's actual unwinding model.
    for ir, expected in ((frames_ir, 1), (frames, 0)):
        ir = ir[:ir.index("declare void @llvm.localescape")]
        asm, _ = run([llc, "-O3", "-mcpu=haswell", "-o", "-", "-"], input=ir)
        fn = re.search(r"^private_frame:.*?^\.Lfunc_end", asm, re.M | re.S)[0]
        assert fn.count("byte Spill") == expected, fn
        assert fn.count("Reload") == expected, fn

    abi_ir = (TESTS / "call-abi.ll").read_text()
    abi, remarks = transform(abi_ir, "rssn-call-abi")
    assert "define internal x86_regcallcc void @leaf_result" in abi, abi
    assert "call x86_regcallcc void @leaf_result" in abi, abi
    assert "sret(%Result)" in abi and "range(i64 0, 1152921504606846976)" in abi, abi
    assert "define internal i64 @address_taken" in abi, abi
    assert "define internal i64 @nonleaf" in abi, abi
    assert "define internal i64 @leaf_forward" in abi, abi
    assert "nonleaf callee" in remarks and "save traffic exceed benefit" in remarks, remarks
    report, remarks = transform(abi_ir, "rssn-call-abi-report")
    assert "define internal fastcc void @leaf_result" in report, report
    assert "result eligible" in remarks, remarks
    for ir, expected in ((abi_ir, True), (abi, False)):
        asm, _ = run([llc, "-O3", "-mcpu=haswell", "-o", "-", "-"], input=ir)
        fn = re.search(r"^leaf_result:.*?^\.Lfunc_end", asm, re.M | re.S)[0]
        assert ("(%rsp)" in fn) == expected, fn

    indirect_fields = """
target triple = "x86_64-pc-windows-msvc"
define internal fastcc double @read_fields(ptr readonly captures(none) %env) {
  %b = getelementptr inbounds ptr, ptr %env, i64 1
  %pa = load ptr, ptr %env
  %pb = load ptr, ptr %b
  %x = load double, ptr %pa
  %y = load double, ptr %pb
  %sum = fadd double %x, %y
  ret double %sum
}
define double @field_caller(ptr %a, ptr %b) {
  %env = alloca [2 x ptr]
  %next = getelementptr inbounds ptr, ptr %env, i64 1
  store ptr %a, ptr %env
  store ptr %b, ptr %next
  %r = call fastcc double @read_fields(ptr %env)
  ret double %r
}
"""
    for target in ("x86_64-pc-windows-msvc", "x86_64-unknown-linux-gnu"):
        ir = indirect_fields.replace("x86_64-pc-windows-msvc", target)
        report, remarks = transform(ir, "rssn-call-abi-report")
        assert "indirect scalar fields 2, scalar loads 2, direct sites 1" in remarks, remarks
        assert "caller stack environments 1" in remarks, remarks
        unchanged, _ = transform(ir, "rssn-call-abi")
        assert report == unchanged, report
        for rejected in (ir.replace(" readonly", ""),
                         ir.replace("load double, ptr %pb", "load volatile double, ptr %pb")):
            _, remarks = transform(rejected, "rssn-call-abi-report")
            assert "indirect scalar fields" not in remarks, remarks

    vector_ir = (TESTS / "libm-vector.ll").read_text()
    report, remarks = transform(vector_ir, "rssn-libm-pressure-report")
    assert "uniform candidate" in remarks and "libm.evaluate" not in report, remarks
    with_veclib, remarks = transform(vector_ir, "rssn-libm-pressure",
                                    extra=("-vector-library=LIBMVEC",))
    assert "libm.evaluate" not in with_veclib and not remarks, remarks
    for target in ("x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc"):
        source = vector_ir.replace("x86_64-unknown-linux-gnu", target)
        guarded, remarks = transform(source, "rssn-libm-pressure")
        assert "uniform guarded llvm.exp.v2f64" in remarks, remarks
        assert "false_edge guarded llvm.log.v4f64" in remarks, remarks
        assert '!{!"branch_weights", i32 80, i32 20}' in guarded, guarded
        for name in ("uniform", "false_edge"):
            fn = body(guarded, name)
            assert fn.index("freeze i1") < fn.index("br i1") < fn.index("call") < fn.index("phi"), fn
        for name in ("nonuniform", "shared", "unprofitable", "strict"):
            assert "libm.evaluate" not in body(guarded, f"keep_{name}"), name
        for name, branches in (("shared_branch", 1), ("dependent_branches", 2)):
            fn = body(guarded, name)
            assert fn.count("br i1") == branches and fn.count("phi <2 x double>") == 2, fn
        asm, _ = run([llc, "-O3", "-mcpu=haswell", "-o", "-", "-"], input=guarded)
        fn = asm.split("uniform:", 1)[1].split("false_edge:", 1)[0]
        assert re.search(r"\bj[a-z]+\s+.*\n(?s:.*?)\bcallq?\s+.*\bexp", fn), fn

    reuse = (TESTS / "libm-reuse.ll").read_text()
    _, remarks = transform(reuse, "rssn-libm-reuse-report")
    assert "mirrored_pair candidate" in remarks, remarks
    for target in ("x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc"):
        source = reuse.replace("x86_64-unknown-linux-gnu", target)
        revised, remarks = transform(source, "rssn-libm-reuse")
        fn = body(revised, "mirrored_pair")
        assert "guarded reuse" in remarks and "2 matched consumers" in remarks, remarks
        assert fn.count("freeze double") == 2 and "icmp ne i64" in fn, fn
        assert fn.count("call double @sinh") == 2, fn
        cold = fn.split("\nexp.recompute:", 1)[1].split("br label", 1)[0]
        assert all(x in cold for x in ("@llvm.exp", "fadd double", "fdiv double")), cold
        assert "phi double [ %d, %right ], [ %r, %exp.recompute ]" in fn, fn
    for modified in (
        reuse.replace("%bx = fmul double %b, 2.000000e+00",
                      "%bx = fmul double %b, 3.000000e+00"),
        reuse.replace("%r = fdiv double 2.000000e+00, %q",
                      "%r = fdiv double %q, 2.000000e+00"),
        reuse.replace("%e = call double @llvm.exp", "%e = call nofpclass(nan) double @llvm.exp"),
        reuse.replace("%r = fdiv double 2.000000e+00, %q",
                      "%r = fdiv double 2.000000e+00, %q, !fpmath !0")
             + "!0 = !{float 2.0}\n",
        reuse.replace("@mirrored_pair(double %x) {", "@mirrored_pair(double %x) strictfp {"),
    ):
        revised, _ = transform(modified, "rssn-libm-reuse")
        assert "exp.recompute" not in revised, revised

    pow_ir = (TESTS / "pow-exponent-pack.ll").read_text()
    packed, remarks = transform(pow_ir, "rssn-pow-exponent-pack")
    assert "packed two pow exponents at 2 windows" in remarks, remarks
    fn = body(packed, "pack_pow_exponents")
    assert fn.count("freeze <2 x double>") == 1, fn
    assert fn.count('asm "vunpckhpd $1, $1, $0"') == 2, fn
    assert fn.count('asm "", "=x,0"') == 2, fn
    lowered, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"],
                     input=packed)
    assert lowered.count("vunpckhpd") == 2, lowered
    # A full-width YMM value spills independently of the preserved XMM budget.
    wide_pow = re.sub(r"double (%a[0-9])", r"<4 x double> \1", pow_ir)
    wide_pow = re.sub(
        r"declare double @measure\([^)]*\)",
        "declare double @measure(double, double, " + ", ".join(["<4 x double>"] * 10) + ")",
        wide_pow,
    )
    for source, cpu in (
        (pow_ir, "core2"),
        (pow_ir.replace("x86_64-pc-windows-msvc", "x86_64-unknown-linux-gnu"), "haswell"),
        (pow_ir.replace("ret double %sum", "%late = fadd double %sum, %e0\n  ret double %late"),
         "haswell"),
        (wide_pow, "haswell"),
    ):
        unchanged, _ = transform(source, "rssn-pow-exponent-pack", cpu=cpu)
        assert "asm" not in unchanged, unchanged
    shared_budget = (TESTS / "pow-exponent-budget.ll").read_text()
    packed, remarks = transform(shared_budget, "rssn-pow-exponent-pack")
    assert packed.count("freeze <2 x double>") == 1, remarks
    # An extra live scalar leaves a second excess XMM slot to remove.
    more_pressure = shared_budget.replace("double %a4)", "double %a4, double %a5)")
    more_pressure = more_pressure.replace(
        "declare double @consume(double, double, double, double, double, double)",
        "declare double @consume(double, double, double, double, double, double, double)",
    )
    packed, remarks = transform(more_pressure, "rssn-pow-exponent-pack")
    assert packed.count("freeze <2 x double>") == 2, remarks

    guarded_pair = pow_ir.replace("declare double @llvm.pow", "declare void @opaque()\n\ndeclare double @llvm.pow")
    for block, low, high in (("positive", "a", "b"), ("negative", "c", "d")):
        guarded_pair = guarded_pair.replace(
            f"  %{high} = call double @llvm.pow",
            f"  br i1 %skip, label %{block}.skip, label %{block}.eval, !prof !1\n\n"
            f"{block}.skip:\n  br label %{block}.done\n\n{block}.eval:\n"
            f"  %{high} = call double @llvm.pow",
        )
        end = guarded_pair.index("  br label %join", guarded_pair.index(f"{block}.eval:"))
        guarded_pair = (guarded_pair[:end] + f"  br label %{block}.done\n\n{block}.done:\n"
                        f"  %{high}.result = phi double [ 1.0, %{block}.skip ], [ %{high}, %{block}.eval ]\n"
                        + guarded_pair[end:])
        guarded_pair = (guarded_pair.replace(f"[ %{low}, %{block} ]", f"[ %{low}, %{block}.done ]")
                        .replace(f"[ %{high}, %{block} ]", f"[ %{high}.result, %{block}.done ]"))
    guarded_pair = guarded_pair.replace("i32 %n,", "i32 %n, i1 %skip,", 1)
    guarded_pair += '\n!1 = !{!"branch_weights", i32 1, i32 100}\n'
    packed, remarks = transform(guarded_pair, "rssn-pow-exponent-pack")
    assert "2 windows, cross-edge windows 2, overlapping calls 4" in remarks, remarks
    assert packed.count('asm "vunpckhpd') == 2, packed
    for rejected in (
        guarded_pair.replace("i32 1, i32 100}", "i32 100, i32 1}"),
        guarded_pair.replace("positive.eval:\n", "positive.eval:\n  call void @opaque()\n"),
    ):
        unchanged, _ = transform(rejected, "rssn-pow-exponent-pack")
        assert "asm" not in unchanged, unchanged

    bank_ir = (TESTS / "fp-register-bank.ll").read_text()
    report, remarks = transform(bank_ir, "rssn-fp-register-bank-report")
    assert "can carry one stack argument in the GPR bank" in remarks, remarks
    assert "can carry one FP temporary in the GPR bank" in remarks, remarks
    assert "fp.bank" not in report, report
    banked, remarks = transform(bank_ir, "rssn-fp-register-bank")
    fn = body(banked, "bank_exponents")
    assert fn.count('asm "vmovq $1, $0", "=x,r"') == 4, fn
    assert "freeze double %exponent" in fn and "freeze double %e1" not in fn, fn
    assert fn.count("@llvm.pow.f64") == 4, fn
    lowered, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"],
                     input=banked)
    assert lowered.count("vmovq") >= 4, lowered
    int_args = ", ".join(f"i64 %g{i}" for i in range(6))
    crowded = bank_ir.replace("double %a9) {", f"double %a9, {int_args}) {{")
    crowded = crowded.replace("double, double, double, double)\n",
                               "double, double, double, double, "
                               + ", ".join(["i64"] * 6) + ")\n", 1)
    crowded = crowded.replace("double %a9)\n", f"double %a9, {int_args})\n")
    wide = re.sub(r"double (%a[0-9])", r"<4 x double> \1", bank_ir)
    wide = re.sub(
        r"declare double @observe\([^)]*\)",
        "declare double @observe(double, double, " + ", ".join(["<4 x double>"] * 10) + ")",
        wide,
    )
    for rejected, cpu in (
        (bank_ir, "core2"),
        (bank_ir.replace("x86_64-pc-windows-msvc", "x86_64-unknown-linux-gnu"), "haswell"),
        (bank_ir.replace("double %a9) {", "double %a9) strictfp {"), "haswell"),
        (bank_ir.replace("  %exponent = fneg", "  call void asm sideeffect \"\", \"\"()\n"
                         "  %exponent = fneg"), "haswell"),
        (crowded, "haswell"),
        (wide, "haswell"),
    ):
        unchanged, _ = transform(rejected, "rssn-fp-register-bank", cpu=cpu)
        assert "fp.bank" not in unchanged, unchanged

    unit_ir = (TESTS / "pow-unit-base.ll").read_text()
    report, remarks = transform(unit_ir, "rssn-pow-unit-base-report")
    assert "can guard 3 powers per node in 2 blocks" in remarks, remarks
    assert "pow.unit." not in report, report
    unit, remarks = transform(unit_ir, "rssn-pow-unit-base")
    assert "guarded 3 powers per node in 2 blocks" in remarks, remarks
    standalone, standalone_remarks = transform(
        unit_ir, "rssn-release", settings=dict(env, RSSN_POW_UNIT_BASE="transform"),
    )
    assert standalone == unit and standalone_remarks == remarks
    untouched, _ = transform(unit_ir, "verify")
    disabled, disabled_remarks = transform(unit_ir, "rssn-release")
    assert disabled == untouched and not disabled_remarks
    fn = body(unit, "unit_bases")
    assert fn.count("pow.unit.cutoff =") == 1, fn
    assert fn.count("freeze double %t") == 1, fn
    assert fn.count("= phi double [") >= 3 and fn.count("@llvm.pow.f64") == 6, fn
    guarded, remarks = transform(hoist_unit_base(unit_ir), "rssn-pow-unit-base")
    assert "guarded 3 powers per node in 2 blocks" in remarks, remarks
    assert "%b0 = fsub double 1.000000e+00, %t" in guarded, guarded
    assert re.search(r"pow\.unit\.evaluate:.*?fsub double 1\.000000e\+00, %unit\.defined", guarded, re.S), guarded
    for rejected in (
        unit_ir.replace("i32 %rounds) {", "i32 %rounds) strictfp {"),
        unit_ir.replace("i32 %rounds) {", "i32 %rounds) denormal_fpenv(preservesign) {"),
        unit_ir.replace("%b0 = fsub", "%b0 = fsub fast"),
        unit_ir.replace("  %sy = fmul", "  call void @observe()\n  %sy = fmul"),
        unit_ir.replace("%s0 = fmul double %initial, %q0",
                        "%s0 = fmul double %initial, %q0\n"
                        "  %q0.extra = fadd double %q0, %s0"),
        unit_ir.replace("%ty = fmul double %m, %y", "%ty = fmul double %m, %x"),
    ):
        unchanged, _ = transform(rejected, "rssn-pow-unit-base")
        assert "pow.unit." not in unchanged, unchanged

    unit_log_ir = (TESTS / "libm-unit-log.ll").read_text()
    report, remarks = transform(unit_log_ir, "rssn-libm-unit-log-report")
    assert "can cache a unit-node result for 2 log/exp chains" in remarks, remarks
    assert "unit.log." not in report, report
    guarded, remarks = transform(unit_log_ir, "rssn-libm-unit-log")
    fn = body(guarded, "unit_log_chain")
    assert "cached a unit-node result for 2 log/exp chains" in remarks, remarks
    assert fn.count("unit.log.cached =") == 1, fn
    assert fn.count("icmp eq i64") == 2 and "unit.log.node = freeze" not in fn, fn
    assert fn.count("@llvm.log.f64") == 2 and fn.count("@llvm.exp.f64") == 5, fn
    already_guarded = unit_log_ir.replace(
        "  %lx = call", "  %one.x = fcmp oeq double %x, 1.0\n"
        "  br i1 %one.x, label %join, label %x.nonunit\n\nx.nonunit:\n  %lx = call",
    ).replace(
        "  %ly = call", "  %one.y = fcmp oeq double 1.0, %y\n"
        "  br i1 %one.y, label %join, label %y.nonunit\n\ny.nonunit:\n  %ly = call",
    ).replace(
        "[ %ex, %left ], [ %ey, %right ]",
        "[ 0.0, %left ], [ %ex, %x.nonunit ], [ 0.0, %right ], [ %ey, %y.nonunit ]",
    )
    unchanged, _ = transform(already_guarded, "rssn-libm-unit-log")
    assert "unit.log.cached" not in unchanged, unchanged
    bits_guarded = already_guarded.replace(
        "%one.x = fcmp oeq double %x, 1.0",
        "%bits.x = bitcast double %x to i64\n  %one.x = icmp eq i64 %bits.x, 4607182418800017408",
    ).replace(
        "%one.y = fcmp oeq double 1.0, %y",
        "%bits.y = bitcast double %y to i64\n  %one.y = icmp eq i64 4607182418800017408, %bits.y",
    )
    unchanged, _ = transform(bits_guarded, "rssn-libm-unit-log")
    assert "unit.log.cached" not in unchanged, unchanged
    unit_path = already_guarded.replace("label %join, label %x.nonunit", "label %x.nonunit, label %join")
    unit_path = unit_path.replace("label %join, label %y.nonunit", "label %y.nonunit, label %join")
    guarded, _ = transform(unit_path, "rssn-libm-unit-log")
    assert "unit.log.cached" in guarded, guarded
    relaxed, _ = transform(
        unit_log_ir.replace("call noundef double @llvm.exp", "call double @llvm.exp"),
        "rssn-libm-unit-log",
    )
    assert relaxed.count("freeze double %x") == 1 and "fcmp oeq" in relaxed, relaxed
    for rejected in (
        unit_log_ir.replace("i32 %rounds) {", "i32 %rounds) strictfp {"),
        unit_log_ir.replace("i32 %rounds) {", "i32 %rounds) denormal_fpenv(preservesign) {"),
        unit_log_ir.replace("%cx = fmul", "%cx = fmul fast"),
        unit_log_ir.replace("%lx = call double @llvm.log.f64", "%lx = call double @log")
        + "\ndeclare double @log(double)\n",
        unit_log_ir.replace("  %cx = fmul", "  call void @observe(double %lx)\n  %cx = fmul"),
        unit_log_ir.replace("%x = fadd double %a, 1.0", "%x = fadd double %a, 2.0"),
    ):
        unchanged, _ = transform(rejected, "rssn-libm-unit-log")
        assert "unit.log.cached" not in unchanged, unchanged

    tail_ir = (TESTS / "libm-tail-cutoff.ll").read_text()
    report, remarks = transform(tail_ir, "rssn-libm-tail-cutoff-report")
    assert remarks.count("can place a shared cutoff before 2 log chains") == 2, remarks
    assert "exp.tail.cutoff" not in report, report
    tail, remarks = transform(tail_ir, "rssn-libm-tail-cutoff")
    assert remarks.count("placed a shared cutoff before 2 log chains") == 2, remarks
    existing = body(tail, "tail_existing")
    zero = body(tail, "tail_zero_numerator")
    assert existing.count("fcmp olt double %x, %exp.tail.cutoff") == 2, existing
    assert "fcmp olt double %score, -7.450000e+02" in existing, existing
    assert zero.count("fdiv double %exp.tail.result") == 2, zero
    assert zero.count("phi double [ 0.000000e+00") >= 2, zero
    shifted, remarks = transform(shifted_tail(tail_ir), "rssn-libm-tail-cutoff")
    assert remarks.count("placed a shared cutoff before 2 log chains") == 3, remarks
    shifted = body(shifted, "tail_shifted")
    assert shifted.count("fdiv double %exp.tail.result") == 2, shifted
    assert "%adjusted = fsub double %score" in shifted, shifted
    for order in ("rssn-fp-class-interval,rssn-libm-tail-cutoff",
                  "rssn-libm-tail-cutoff,rssn-fp-class-interval"):
        composed, remarks = transform(shifted_tail(tail_ir), order)
        assert remarks.count("placed a shared cutoff before 2 log chains") == 3, remarks
        assert "@llvm.is.fpclass" not in body(composed, "tail_shifted").split("header:")[1].split("valid:")[0], composed
    for rejected in (
        tail_ir.replace("i32 %rounds) {", "i32 %rounds) strictfp {"),
        tail_ir.replace("i32 %rounds) {", "i32 %rounds) denormal_fpenv(preservesign) {"),
        tail_ir.replace("double %weight, i32 615", "double %weight, i32 3"),
    ):
        unchanged, _ = transform(rejected, "rssn-libm-tail-cutoff")
        assert "exp.tail.cutoff" not in unchanged, unchanged
    for rejected in (
        tail_ir.replace("[ 0.0, %left.score ]", "[ 1.0, %left.score ]"),
        tail_ir.replace("left.score:\n", "left.score:\n  call void @observe()\n"),
        tail_ir.replace("%decay = fmul", "%decay = fmul fast"),
        tail_ir.replace("%x = call noundef double @llvm.exp.f64(double %t)",
                        "%x = load double, ptr @external_node", 1)
        + "\n@external_node = external global double\n",
    ):
        unchanged, _ = transform(rejected, "rssn-libm-tail-cutoff")
        assert "exp.tail.cutoff" not in body(unchanged, "tail_existing"), unchanged

    exit_ir = (TESTS / "loop-exit-compression.ll").read_text()
    _, remarks = transform(exit_ir, "rssn-loop-exit-compression-report")
    assert "can move 2 exit operations" in remarks, remarks
    early, remarks = transform(exit_ir, "rssn-loop-exit-compression")
    assert "2 fewer live SIMD values" in remarks, remarks
    assert "%coefficient = fdiv" in body(early, "compress_exit").split("loop:", 1)[0], early
    lowered, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"],
                     input=early)
    assert "vdivsd" in lowered
    for kept in (
        exit_ir.replace("i32 %n) {", "i32 %n) strictfp {"),
        exit_ir.replace("ret double %result",
                        "call void @touch(double %a, double %b, double %c)\n  ret double %result"),
        exit_ir.replace("%product = fmul double %b, %c", "%product = fmul double %b, %y"),
    ):
        original, _ = transform(kept, "rssn-loop-exit-compression")
        assert "%coefficient = fdiv" in body(original, "compress_exit").split("exit:", 1)[1]

    switch_reads = (TESTS / "switch-load-reuse.ll").read_text()
    original, _ = transform(switch_reads, "verify")
    report, remarks = transform(switch_reads, "rssn-switch-load-reuse-report")
    assert report == original and "3 switch arms covers 1 narrower reads" in remarks
    rewritten, remarks = transform(switch_reads, "rssn-switch-load-reuse")
    switch_fn = body(rewritten, "hash_mode")
    assert "reused" in remarks and switch_fn.count("load ") == 1, rewritten
    assert "load b64, ptr %source, align 1" in switch_fn, switch_fn
    assert "extractelement <2 x b32>" in switch_fn, switch_fn
    assert switch_fn.index("load ") < switch_fn.index("switch i32"), switch_fn
    run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"],
        input=rewritten)
    for kept in (
        switch_reads.replace("%w2 = load i64, ptr %source", "%w2 = load i64, ptr %out"),
        switch_reads.replace("%w1 = load i64", "store i64 7, ptr %source\n  %w1 = load i64"),
        switch_reads.replace("%part = load i32", "%part = load volatile i32"),
        switch_reads.replace("%w2 = load i64", "%w2 = load atomic i64")
                    .replace("%w2 = load atomic i64, ptr %source, align 1",
                             "%w2 = load atomic i64, ptr %source unordered, align 8"),
        switch_reads.replace("%w0 = load i64", "fence acquire\n  %w0 = load i64"),
    ):
        retained, _ = transform(kept, "rssn-switch-load-reuse")
        assert "switch.bytes" not in retained, retained

    preheader = (TESTS / "loop-preheader-pressure.ll").read_text()
    original, _ = transform(preheader, "verify")
    report, remarks = transform(preheader, "rssn-loop-preheader-pressure-report")
    assert report == original and "placements 2" in remarks and "eligible" in remarks, remarks
    rewritten, remarks = transform(preheader, "rssn-loop-preheader-pressure")
    fn = body(rewritten, "dispatch")
    assert "applied" in remarks, remarks
    assert fn.index("arm:") < fn.index("%lane =") < fn.index("%wide =") < fn.index("@consume"), fn
    run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"], input=rewritten)
    costs = []
    for vector_ir in (preheader, preheader.replace("<4 x ", "<8 x ")):
        _, remarks = transform(vector_ir, "rssn-loop-preheader-pressure-report",
                               cpu="skylake-avx512")
        match = re.search(r"estimated call pressure ([0-9.e+]+) ->", remarks)
        assert match, remarks
        costs.append(float(match[1]))
    assert 0 < costs[0] < costs[1], costs
    edge_use = preheader.replace("  %out = phi double", "  %edge = phi <4 x i64> "
                                 "[ %wide, %arm ], [ zeroinitializer, %call ]\n"
                                 "  %out = phi double").replace(
        "  %next = add", "  call void @consume(<4 x i64> %edge, i64 %i)\n  %next = add")
    rewritten, _ = transform(edge_use, "rssn-loop-preheader-pressure")
    fn = body(rewritten, "dispatch")
    assert fn.index("arm:") < fn.index("%wide =") < fn.index("latch:"), fn
    costly = preheader.replace("i64 noundef %v,", "<4 x i64> noundef %v, <4 x i64> noundef %w,")
    costly = "\n".join(
        "  %wide = add <4 x i64> %v, %w" if line.startswith("  %wide =") else line
        for line in costly.splitlines() if not line.startswith("  %lane ="))
    unchanged, remarks = transform(costly, "rssn-loop-preheader-pressure")
    assert "placements 1" in remarks and "rejected" in remarks, remarks
    assert "%wide =" in body(unchanged, "dispatch").split("header:", 1)[0], unchanged
    for retained in (
        preheader.replace("i64 noundef %v", "i64 %v"),
        preheader.replace("i32 1, i32 1000", "i32 1000, i32 1"),
        preheader.replace("x86_64-pc-windows-msvc", "x86_64-unknown-linux-gnu"),
        preheader.replace("  ret double %out", "  call void @consume(<4 x i64> %wide, i64 0)\n"
                          "  ret double %out"),
    ):
        unchanged, _ = transform(retained, "rssn-loop-preheader-pressure")
        assert "%wide =" in body(unchanged, "dispatch").split("header:", 1)[0], unchanged

    endpoint = (TESTS / "fp-endpoint-guards.ll").read_text()
    report, remarks = transform(endpoint, "rssn-fp-endpoint-guards-report")
    assert "can merge endpoint guards, 3 exp chain operations guarded" in remarks, remarks
    assert "weight.nonzero.finite" not in report, report
    guarded, remarks = transform(endpoint, "rssn-fp-endpoint-guards")
    fn = body(guarded, "endpoint_node")
    assert "merged endpoint guards" in remarks, remarks
    assert fn.index("weight.nonzero.finite") < fn.index("lower:") < fn.index("@llvm.exp"), fn
    assert "freeze i1 %left" in fn and "%zero = fcmp" not in fn, fn
    report, remarks = transform(endpoint, "rssn-fp-class-interval-report")
    assert "mapped_weight can compact a finite nonzero class check" in remarks, remarks
    assert "@llvm.is.fpclass" in body(report, "mapped_weight"), report
    classified, _ = transform(endpoint, "rssn-fp-class-interval")
    fn = body(classified, "mapped_weight")
    assert "%weight.zero.or.nonfinite = icmp uge i64" in fn, fn
    assert "@llvm.is.fpclass" not in fn, fn
    carried = body(classified, "carried_magnitude")
    assert "%weight.zero.or.nonfinite = icmp uge i64" in carried, carried
    unproved = endpoint.replace("%initial = call double @llvm.fabs.f64(double %first)",
                                "%initial = fneg double %first")
    unchanged, _ = transform(unproved, "rssn-fp-class-interval")
    assert "@llvm.is.fpclass" in body(unchanged, "carried_magnitude"), unchanged
    inverse = endpoint.replace("%node = call double @llvm.exp.f64(double %arg)",
                               "%node = fadd double %arg, 2.0")
    inverse = inverse.replace("i32 615)", "i32 408)")
    changed, _ = transform(inverse, "rssn-fp-class-interval")
    assert "weight.nonzero.finite = icmp ult i64" in body(changed, "mapped_weight"), changed
    for rejected in (
        endpoint.replace("double %seed) {", "double %seed) strictfp {"),
        endpoint.replace("i32 615)", "i32 519)"),
        endpoint.replace("double %weight, i32 615)", "double %factor, i32 615)"),
    ):
        unchanged, _ = transform(rejected, "rssn-fp-class-interval")
        assert "@llvm.is.fpclass" in body(unchanged, "mapped_weight"), unchanged
    for rejected in (
        endpoint.replace("i32 %n) {", "i32 %n) strictfp {"),
        endpoint.replace("i32 %n) {", "i32 %n) denormal_fpenv(preservesign) {"),
        endpoint.replace("%positive = fcmp ugt", "%positive = fcmp ogt"),
        endpoint.replace("%scaled = fmul", "call void @observe()\n  %scaled = fmul"),
        endpoint.replace("[ 0.000000e+00, %upper ]", "[ 1.000000e+00, %upper ]"),
        endpoint.replace("%distance, 5.000000e-01", "%distance, -5.000000e-01"),
    ):
        unchanged, _ = transform(rejected, "rssn-fp-endpoint-guards")
        assert "weight.nonzero.finite" not in unchanged, unchanged
    for tail in ("%ret = fadd double %result, %distance",
                 "%ret = select i1 %valid, double %result, double 0.000000e+00"):
        still_used = endpoint.replace("ret double %result", tail + "\n  ret double %ret")
        without_sink, remarks = transform(still_used, "rssn-fp-endpoint-guards")
        assert "merged endpoint guards, 0 exp chain operations guarded" in remarks, remarks
        fn = body(without_sink, "endpoint_node")
        assert fn.index("@llvm.exp") < fn.index("weight.nonzero.finite"), fn

    scalar_packs = (TESTS / "scalar-load-pack.ll").read_text()
    packed, remarks = transform(scalar_packs, "rssn-scalar-load-pack")
    assert "9 scalar loads into 3 packs of 4 bytes" in remarks, remarks
    assert packed.count("load <4 x i8>") == 3, packed
    assert packed.count("bitcast <2 x i8>") == 3, packed
    assert packed.count("load i32") == 3, packed
    report, remarks = transform(scalar_packs, "rssn-scalar-load-pack-report")
    assert "9 scalar loads into 3 packs of 4 bytes" in remarks, remarks
    assert "load <4 x i8>" not in report, report
    win_packs = scalar_packs.replace("x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc")
    win_packed, remarks = transform(win_packs, "rssn-scalar-load-pack")
    assert "12 scalar loads into 3 packs of 8 bytes" in remarks, remarks
    assert win_packed.count("load <8 x i8>") == 3, win_packed
    win_records = win_packs
    for field in "abc":
        win_records = win_records.replace(f"load i32, ptr %r{field}.high, align 4",
                                          "trunc i64 %size to i32")
    win_records = win_records.replace("[8 x i8]", "[4 x i8]")
    record_packed, remarks = transform(win_records, "rssn-scalar-load-pack")
    assert "9 scalar loads into 3 packs of 4 bytes" in remarks, remarks
    assert "independent reads" in remarks
    assert record_packed.count("load <4 x i8>") == 3, record_packed
    dependent_records = "\n".join(
        line for line in win_records.splitlines()
        if not line.startswith(("  %rb =", "  %rc =", "  %rb.high =", "  %rc.high ="))
    )
    dependent_records = dependent_records.replace(
        "  %b0 = load",
        "  %index = zext i16 %a0 to i64\n"
        "  %rb = getelementptr inbounds [4 x i8], ptr %b, i64 %index\n"
        "  %rc = getelementptr inbounds [4 x i8], ptr %c, i64 %index\n"
        "  %b0 = load",
    )
    analyzed, remarks = transform(dependent_records, "rssn-scalar-load-pack")
    assert "dependent reads, below rewrite threshold" in remarks, remarks
    assert "load <4 x i8>" not in analyzed, analyzed
    params = ", ".join(f"i64 %live{i}" for i in range(24))
    pressured_records = dependent_records.replace("i64 %size)", f"i64 %size, {params})", 1)
    live_types = ", ".join(["i64"] * 24)
    pressured_records = pressured_records.replace(
        "declare void @consume_fields",
        f"declare void @consume_live({live_types})\ndeclare void @consume_fields",
    )
    pressured_records = pressured_records.replace(
        "%next = add", f"call void @consume_live({params})\n  %next = add",
    )
    record_packed, remarks = transform(pressured_records, "rssn-scalar-load-pack")
    assert "dependent reads, eligible" in remarks, remarks
    assert record_packed.count("load <4 x i8>") == 3, record_packed
    other_records = win_records.replace("[4 x i8]", "[8 x i8]")
    other_packed, _ = transform(other_records, "rssn-scalar-load-pack")
    assert "load <4 x i8>" not in other_packed, other_packed
    for blocked in (
        scalar_packs.replace("%a2 = load", "store i8 7, ptr %ra.byte, align 1\n  %a2 = load"),
        scalar_packs.replace("%a2 = load", "fence acquire\n  %a2 = load"),
        scalar_packs.replace("%a0 = load i16", "%a0 = load volatile i16")
                    .replace("%c0 = load i16", "%c0 = load volatile i16"),
        scalar_packs.replace("ptr %ra, i64 3", "ptr %ra, i64 4")
                    .replace("ptr %rc, i64 3", "ptr %rc, i64 4"),
    ):
        unchanged, _ = transform(blocked, "rssn-scalar-load-pack")
        assert "load <4 x i8>" not in unchanged, unchanged
    lowered, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs", "-o", "-", "-"],
                     input=packed)
    assert "field_loop" in lowered, lowered

    win_pressure = (TESTS / "win64-call-pressure.ll").read_text()
    win_report, remarks = transform(win_pressure, "rssn-spill-pressure-report")
    assert "Win64 ABI call sites 1, over preservation budget 1" in remarks, remarks
    assert "at selected site GP 1 / 8, XMM-sized 2 / 10, native wide SIMD 1" in remarks, remarks
    assert "Win64 possible libm sites 1, over preservation budget 1" in remarks, remarks
    assert "at selected site GP 2 / 8, XMM-sized 1 / 10, native wide SIMD 1" in remarks, remarks
    win_transform, diagnostics = transform(win_pressure, "rssn-spill-pressure")
    assert win_report == win_transform, win_report
    assert "Win64" not in diagnostics, diagnostics
    win_byte_pressure = win_pressure.replace(
        "%saved = load i64, ptr %p, align 8",
        "%saved.raw = load b64, ptr %p, align 8")
    win_byte_pressure = win_byte_pressure.replace(
        "%bits = add i64 %saved, %r",
        "%saved = bitcast b64 %saved.raw to i64\n  %bits = add i64 %saved, %r")
    _, remarks = transform(win_byte_pressure, "rssn-spill-pressure-report")
    assert "at selected site GP 1 / 8, XMM-sized 2 / 10, native wide SIMD 1" in remarks, remarks
    assert "at selected site GP 2 / 8, XMM-sized 1 / 10, native wide SIMD 1" in remarks, remarks

    # An outlined Rust closure uses fastcc and the Win64 preservation mask.
    fast_win = win_pressure.replace("win64cc i64 @opaque", "fastcc i64 @opaque")
    _, remarks = transform(fast_win, "rssn-spill-pressure-report")
    assert "Win64 ABI call sites 1, over preservation budget 1" in remarks, remarks
    mir, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs",
                  "-stop-after=virtregrewriter", "-o", "-", "-"], input=fast_win)
    assert re.search(r"CALL64pcrel32 @opaque, csr_win64\b", mir), mir
    _, remarks = transform(win_pressure, "rssn-spill-pressure-report", cpu="x86-64")
    assert "native wide SIMD 0, unclassified SIMD peak 2" in remarks, remarks

    stack_args = win_pressure.split("declare", 1)[0] + """
declare fastcc i64 @stack_target(double, i64, ptr, double, i64, double)
define fastcc i64 @win_stack_probe(double %a, i64 %b, ptr %c, double %d,
                                    i64 %e, double %f, i64 %unused) {
  %r = notail call fastcc i64 @stack_target(double %a, i64 %b, ptr %c, double %d,
                                           i64 %e, double %f)
  %v = add i64 %r, %b
  ret i64 %v
}
"""
    stack_report, remarks = transform(stack_args, "rssn-spill-pressure-report")
    assert "Win64 incoming scalar stack slots 3, used parameters 2" in remarks, remarks
    stack_call = "Win64 calls with scalar stack slots 1, at selected site slots 2"
    assert f"{stack_call}, callee stack_target" in remarks, remarks
    stack_off, _ = transform(stack_args)
    assert stack_report == stack_off, stack_report
    mir, _ = run([llc, "-O3", "-mcpu=haswell", "-verify-machineinstrs",
                  "-stop-after=virtregrewriter", "-o", "-", "-"], input=stack_args)
    assert "ADJCALLSTACKDOWN64 48" in mir, mir
    assert re.search(r"CALL64pcrel32 @stack_target, csr_win64\b", mir), mir
    for rejected in (
        stack_args.replace("double %f", "<2 x double> %f")
                  .replace("i64, double)", "i64, <2 x double>)"),
        stack_args.replace("ptr %c", "ptr byval(i64) %c")
                  .replace("i64, ptr,", "i64, ptr byval(i64),"),
    ):
        _, remarks = transform(rejected, "rssn-spill-pressure-report")
        assert "scalar stack slots" not in remarks, remarks

    raw = kernel()
    optimized, remarks = transform(raw)
    call = optimized.index("call double @rssn_gain")
    assert optimized[:call].count("load <4 x double>") == 0, remarks
    assert optimized.count("load <4 x double>") == 12, optimized
    assert optimized.index("%p0 =") > call, optimized
    for attr in ("nosync", "nofree", "nounwind"):
        blocked, _ = transform(raw.replace(f" {attr}", ""))
        assert blocked[: blocked.index("call double @rssn_gain")].count(
            "load <4 x double>"
        ) == 12, attr
    report, remarks = transform(raw, "rssn-spill-pressure-report")
    assert report[: report.index("call double @rssn_gain")].count(
        "load <4 x double>"
    ) == 12, report
    assert "SIMD live through 1 call sites 12 -> 12" in remarks, remarks
    _, remarks = transform(raw, "default<O3>", settings=dict(env, RSSN_SPILL_MODE="report"))
    assert remarks.count("rssn_callback_kernel GP peak") == 1, remarks
    _, remarks = transform(raw, "default<O2>", settings=dict(env, RSSN_SPILL_MODE="report"))
    assert remarks == "", remarks
    _, remarks = transform(
        raw,
        settings=dict(env, RSSN_SPILL_FILTER="no-matching-function"),
    )
    assert remarks == "", remarks


def edge_restore_tests(llc, plugin):
    mir = (TESTS / "edge-restore.mir").read_text()
    env = dict(os.environ, RSSN_SPILL_FILTER="")

    def run_mir(source, pipeline, extra=()):
        return run([
            llc, f"-load-pass-plugin={plugin}", "-x=mir", f"-run-pass={pipeline}",
            *extra,
            "-pass-remarks-analysis=rssn-(edge-restore|machine-stack-traffic)",
            "-verify-machineinstrs",
            "-o", "-", "-",
        ], input=source, env=env)

    base, _ = run_mir(mir, "none")
    report, remarks = run_mir(mir, "rssn-machine-stack-traffic")
    assert report == base
    assert "calls 2, spill reads 2 (2 after calls, 0 wider than 8 bytes, 16 known bytes)" in remarks
    assert "spill writes 1, immutable fixed-stack reads 0 (0 after calls)" in remarks
    assert "constant-pool reads 2 (2 after calls)" in remarks
    folded = mir.replace("\nstack:\n", "\nfixedStack:\n  - { id: 0, offset: 40, size: 8, alignment: 8, isImmutable: true }\nstack:\n")
    folded = folded.replace(
        "    $rax = MOV64rr $r10",
        "    CMP64mi8 %stack.0, 1, $noreg, 0, $noreg, 0, implicit-def dead $eflags :: (load (s64) from %stack.0)\n"
        "    $rax = MOV64rm %fixed-stack.0, 1, $noreg, 0, $noreg :: (load (s64) from %fixed-stack.0)",
    )
    _, remarks = run_mir(folded, "rssn-machine-stack-traffic")
    assert "spill reads 3 (2 after calls, 0 wider than 8 bytes, 24 known bytes)" in remarks
    assert "immutable fixed-stack reads 1 (0 after calls)" in remarks
    reported, remarks = run_mir(mir, "rssn-edge-restore-report")
    assert base == reported and "can share 2 restores on 2 incoming edges" in remarks
    changed, remarks = run_mir(mir, "rssn-edge-restore")
    assert "shared 2 restores on 2 incoming edges" in remarks
    assert "in destination" in remarks and changed.count("  bb.") == base.count("  bb.")
    assert changed.count("from constant-pool") == 1
    assert changed.count("from %stack.0") == 1
    run([llc, "-x=mir", "-O3", "-start-after=virtregrewriter",
         "-verify-machineinstrs", "-filetype=obj", "-o", os.devnull, "-"], input=changed)

    # Only one branch needs the restored registers. The same rule covers an
    # explicit branch target and a single-predecessor fallthrough.
    start, remaining = mir.split("  bb.2:\n", 1)
    _, join = remaining.split("  bb.3:\n", 1)
    private = (start + "  bb.2:\n    successors: %bb.4\n\n    JMP_1 %bb.4\n\n"
               + "  bb.3:\n" + join)
    join_body, tail = join.split("  bb.4:\n", 1)
    fallthrough = (
        start.replace("    JCC_1 %bb.3, 4, implicit $eflags\n    JMP_1 %bb.4",
                      "    JCC_1 %bb.4, 5, implicit $eflags")
        + "  bb.3:\n" + join_body
        + "  bb.2:\n    successors: %bb.4\n\n    JMP_1 %bb.4\n\n"
        + "  bb.4:\n" + tail
    )
    for case in (private, fallthrough):
        before, _ = run_mir(case, "none")
        observed, remarks = run_mir(case, "rssn-edge-restore-report")
        assert observed == before and "can sink 2 restores on 1 incoming edge" in remarks
        after, remarks = run_mir(case, "rssn-edge-restore")
        assert "sank 2 restores on 1 incoming edge, in destination" in remarks
        assert after.count("  bb.") == before.count("  bb.")
        assert after.index("from %stack.0") > after.index("JCC_1", after.index("CALL64pcrel32"))

    # Another incoming path must retain its own registers. The unconditional
    # jump in the layout predecessor makes a shared fallthrough block legal.
    other_pred = (
        "  bb.4:\n    successors: %bb.3\n\n"
        "    $r10 = MOV64ri 1\n    $r11 = MOV64ri 2\n    JMP_1 %bb.3"
    )
    mixed = mir.replace("  bb.4:\n    $rax = MOV32ri64 0\n    RET64 implicit $rax",
                        other_pred)
    changed, remarks = run_mir(mixed, "rssn-edge-restore")
    assert "on fallthrough" in remarks
    assert changed.count("from constant-pool") == 1
    assert "MOV64ri 1" in changed and "MOV64ri 2" in changed
    elf = mixed.replace("x86_64-pc-windows-msvc", "x86_64-unknown-linux-gnu")
    _, remarks = run_mir(elf, "rssn-edge-restore", ("-basic-block-sections=all",))
    assert "with extra jump" in remarks
    noreturn = elf.replace(
        "  declare i32 @opaque()", "  declare i32 @opaque()\n  declare void @die() noreturn nounwind",
    ).replace(
        "  bb.3:\n", "  bb.5:\n    successors:\n\n"
        "    CALL64m $rip, 1, $noreg, target-flags(x86-gotpcrel) @die, $noreg, csr_64,"
        " implicit $rsp, implicit-def $rsp :: (load (s64) from got)\n\n"
        "  bb.3:\n",
    )
    _, remarks = run_mir(noreturn, "rssn-edge-restore")
    assert "on fallthrough" in remarks
    indirect = noreturn.replace("target-flags(x86-gotpcrel) @die", "@die")
    _, remarks = run_mir(indirect, "rssn-edge-restore")
    assert "with extra jump" in remarks

    # The extra path now falls through to the join from its layout
    # predecessor. A shared block can only be placed out of line.
    cold_penalty = mixed.replace(other_pred + "\n", "").replace(
        "  bb.3:\n", other_pred.removesuffix("    JMP_1 %bb.3") + "\n  bb.3:\n",
    )
    changed, remarks = run_mir(cold_penalty, "rssn-edge-restore")
    assert "with extra jump" in remarks and changed.count("from %stack.0") == 1
    hot_exit = cold_penalty.replace(
        "%bb.3(0x20000000), %bb.4(0x60000000)",
        "%bb.3(0x79000000), %bb.4(0x07000000)",
    )
    alternate_alias = cold_penalty.replace(
        "  bb.4:\n    successors: %bb.3\n\n",
        "  bb.4:\n    successors: %bb.3\n    liveins: $r10d\n\n"
        "    TEST32rr $r10d, $r10d, implicit-def dead $eflags\n",
    )
    overwrite = cold_penalty.replace(
        "    TEST32rr $eax, $eax", "    MOV64mi32 %stack.0, 1, $noreg, 0, $noreg, 1"
        " :: (store (s64) into %stack.0)\n    TEST32rr $eax, $eax",
    )
    for source in (hot_exit, alternate_alias, overwrite):
        expected, _ = run_mir(source, "none")
        actual, remarks = run_mir(source, "rssn-edge-restore")
        assert actual == expected and not remarks, remarks


def portable_target_tests(opt, llc, plugin, backends):
    targets = {
        "aarch64": ("aarch64-unknown-linux-gnu", "aarch64-pc-windows-msvc",
                    "aarch64-apple-darwin", "aarch64-unknown-freebsd"),
        "arm": ("armv7-unknown-linux-gnueabihf", "thumbv7a-pc-windows-msvc"),
        "x86": ("i686-unknown-linux-gnu", "i686-pc-windows-msvc"),
        "x86-64": ("x86_64-unknown-linux-gnu", "x86_64-pc-windows-msvc",
                   "x86_64-apple-darwin", "x86_64-unknown-freebsd"),
    }
    passes = (
        ("pow-unit-base", "guarded 3 powers per node in 2 blocks"),
        ("fp-endpoint-guards", "merged endpoint guards"),
        ("libm-unit-log", "cached a unit-node result for 2 log/exp chains"),
    )
    for name, expected in passes:
        original = (TESTS / f"{name}.ll").read_text()
        for backend, triples in targets.items():
            if backend not in backends:
                print(f"Skipping {name} for {backend}: backend absent from this LLVM SDK")
                continue
            for target in triples:
                source = original.replace("x86_64-unknown-linux-gnu", target)
                changed, remarks = run([
                    opt, f"-load-pass-plugin={plugin}", f"-passes=function(rssn-{name})",
                    f"-pass-remarks-analysis=rssn-{name}", "-verify-each", "-S", "-o", "-", "-",
                ], input=source, env=dict(os.environ, RSSN_SPILL_FILTER=""))
                assert expected in remarks, (target, remarks)
                run([llc, "-O3", "-verify-machineinstrs", "-filetype=obj", "-o", os.devnull, "-"],
                    input=changed)


def rustc_tests(rustc, plugin, linux_x86):
    with tempfile.TemporaryDirectory(prefix="rssn-spill-verify-") as tmp:
        probe = Path(tmp) / "probe.rs"
        probe.write_text(
            """
#[inline(never)]
#[unsafe(no_mangle)]
pub extern "C" fn rssn_pressure_probe(a: f64, b: f64) -> f64 {
    std::hint::black_box(a) + 2.0 * std::hint::black_box(b)
}
fn main() { println!("{}", rssn_pressure_probe(2., 4.)); }
"""
        )
        env = dict(
            os.environ,
            RUSTC_BOOTSTRAP="1",
            RSSN_SPILL_MODE="report",
            RSSN_SPILL_FILTER="rssn_pressure_probe",
            RSSN_SPILL_ABI="off",
        )
        for lto in ("off", "thin", "fat"):
            exe = Path(tmp) / (lto + (".exe" if os.name == "nt" else ""))
            _, remarks = run(
                [
                    rustc,
                    "--edition=2024",
                    "-Copt-level=3",
                    "-Ccodegen-units=1",
                    f"-Clto={lto}",
                    "-Cremark=rssn-spill-pressure",
                    "-Cllvm-args=-pass-remarks-analysis=rssn-spill-pressure",
                    "-Zverify-llvm-ir",
                    f"-Zllvm-plugins={plugin}",
                    str(probe),
                    "-o",
                    str(exe),
                ],
                env=env,
            )
            if linux_x86:
                assert remarks.count("rssn_pressure_probe GP peak") == 1, (lto, remarks)
            assert run([str(exe)])[0].strip() == "10", lto


def reuse_runtime_test(opt, llc, rustc, plugin):
    with tempfile.TemporaryDirectory(prefix="rssn-libm-reuse-") as tmp:
        tmp = Path(tmp)
        runner = tmp / "runner.rs"
        runner.write_text("""
use std::sync::atomic::{AtomicUsize, Ordering};
static CALLS: AtomicUsize = AtomicUsize::new(0);
unsafe extern "C" { fn mirrored_pair(x: f64) -> f64; }

#[unsafe(no_mangle)]
pub extern "C" fn sinh(x: f64) -> f64 {
    CALLS.fetch_add(1, Ordering::Relaxed);
    if x < -1.5 { 2.0 * x } else { x }
}

fn main() {
    for x in [1.0, 2.0, -2.0] {
        let value = unsafe { mirrored_pair(std::hint::black_box(x)) };
        println!("{:016x}", value.to_bits());
    }
    assert_eq!(CALLS.load(Ordering::Relaxed), 6);
}
""")
        base = (TESTS / "libm-reuse.ll").read_text()
        other, _ = run([
            opt, f"-load-pass-plugin={plugin}", "-passes=function(rssn-libm-reuse)",
            "-verify-each", "-S", "-o", "-", "-",
        ], input=base, env=dict(os.environ, RSSN_SPILL_FILTER=""))
        results = []
        for mode, ir in (("base", base), ("reuse", other)):
            obj, exe = tmp / f"{mode}.o", tmp / mode
            run([llc, "-O3", "-filetype=obj", "-o", str(obj), "-"], input=ir)
            run([rustc, "--edition=2024", "-O", f"-Clink-arg={obj}",
                 str(runner), "-o", str(exe)])
            results.append(run([str(exe)])[0])
        assert results[0] == results[1], results


def endpoint_runtime_test(opt, llc, rustc, plugin):
    with tempfile.TemporaryDirectory(prefix="rssn-endpoint-") as tmp:
        tmp = Path(tmp)
        runner = tmp / "class-boundaries.rs"
        runner.write_text("""
unsafe extern "C" {
    fn endpoint_node(weight: f64, arg: f64, left: bool, n: i32) -> f64;
    fn mapped_weight(factor: f64, seed: f64) -> u32;
    fn carried_magnitude(first: f64, candidate: f64, n: u32) -> u32;
}

fn main() {
    for weight in [
        0.0, -0.0, f64::from_bits(1), f64::from_bits(0x8000_0000_0000_0001),
        f64::MIN_POSITIVE, -f64::MIN_POSITIVE, 1.0, -1.0, f64::MAX, -f64::MAX,
        f64::INFINITY, f64::NEG_INFINITY, f64::NAN, -f64::NAN,
        f64::from_bits(0x7ff0_0000_0000_0001), f64::from_bits(0xfff0_0000_0000_0001),
    ] {
        for exponent in [f64::NEG_INFINITY, -30.0, 0.0, 36.0, 38.0, 709.0,
                         710.0, f64::INFINITY, f64::NAN] {
            let count = unsafe { mapped_weight(weight, exponent) };
            println!("{count}");
            println!("{}", unsafe { carried_magnitude(weight, exponent, 3) });
            for left in [true, false] {
                let output = unsafe { endpoint_node(weight, exponent, left, 3) };
                println!("{:016x}", output.to_bits());
            }
        }
    }
}
""")
        base = (TESTS / "fp-endpoint-guards.ll").read_text()
        changed, _ = run([
            opt, f"-load-pass-plugin={plugin}",
            "-passes=function(rssn-fp-endpoint-guards,rssn-fp-class-interval)",
            "-verify-each", "-S", "-o", "-", "-",
        ], input=base, env=dict(os.environ, RSSN_SPILL_FILTER=""))
        results = []
        for mode, ir in (("base", base), ("guarded", changed)):
            obj, exe = tmp / f"{mode}.o", tmp / mode
            run([llc, "-O3", "-verify-machineinstrs", "-filetype=obj",
                 "-o", str(obj), "-"], input=ir)
            run([rustc, "--edition=2024", "-O", f"-Clink-arg={obj}",
                 str(runner), "-o", str(exe)])
            results.append(run([str(exe)])[0])
        assert results[0] == results[1], results


def pow_unit_runtime_test(opt, llc, rustc, plugin):
    with tempfile.TemporaryDirectory(prefix="rssn-unit-pow-") as tmp:
        tmp = Path(tmp)
        runner = tmp / "pow-boundaries.rs"
        runner.write_text("""
unsafe extern "C" {
    fn unit_bases(t: f64, x: f64, y: f64, e0: f64, e1: f64, e2: f64,
                  initial: f64, rounds: i32) -> f64;
}

fn main() {
    let cutoff = 2.0_f64.powi(-55).to_bits();
    let snan = f64::from_bits(0x7ff0_0000_0000_0001);
    for t in [
        0.0, -0.0, f64::from_bits(1), f64::MIN_POSITIVE,
        f64::from_bits(cutoff - 1), f64::from_bits(cutoff),
        f64::from_bits(cutoff + 1), 2.0_f64.powi(-56), 2.0_f64.powi(-54),
        -2.0_f64.powi(-55), 0.5, f64::INFINITY, f64::NAN, snan,
    ] {
        for (x, y) in [
            (1.0, -1.0), (0.4, 0.99), (f64::from_bits(2.0_f64.to_bits() - 1), -1.0),
            (3.0, -4.0), (0.0, f64::from_bits(1)), (f64::MAX, 1.0),
            (f64::INFINITY, 1.0), (f64::NAN, snan),
        ] {
            for exponent in [0.0, -2.25, 1.0, f64::INFINITY,
                             f64::NEG_INFINITY, f64::NAN, snan] {
                for initial in [-0.0, 0.625, snan] {
                    for rounds in [1, 2] {
                        let output = unsafe {
                            unit_bases(t, x, y, exponent, -0.25, exponent, initial, rounds)
                        };
                        println!("{:016x}", output.to_bits());
                    }
                }
            }
        }
    }
}
""")
        base = (TESTS / "pow-unit-base.ll").read_text()
        variants = {"base": base}
        for mode, ir in (("guarded", base), ("shared", hoist_unit_base(base))):
            changed, _ = run([
                opt, f"-load-pass-plugin={plugin}", "-passes=function(rssn-pow-unit-base)",
                "-verify-each", "-S", "-o", "-", "-",
            ], input=ir, env=dict(os.environ, RSSN_SPILL_FILTER=""))
            assert "pow.unit.cutoff" in changed, changed
            variants[mode] = changed
        results = []
        for mode, ir in variants.items():
            obj, exe = tmp / f"{mode}.o", tmp / mode
            run([llc, "-O3", "-verify-machineinstrs", "-filetype=obj",
                 "-o", str(obj), "-"], input=ir)
            run([rustc, "--edition=2024", "-O", f"-Clink-arg={obj}",
                 str(runner), "-o", str(exe)])
            results.append(run([str(exe)])[0].splitlines())
        assert all(len(outputs) == 4704 for outputs in results)
        mismatch = next(((i, values) for i, values in enumerate(zip(*results))
                         if len(set(values)) != 1), None)
        assert mismatch is None, mismatch


def unit_log_runtime_test(opt, llc, rustc, plugin):
    with tempfile.TemporaryDirectory(prefix="rssn-unit-log-") as tmp:
        tmp = Path(tmp)
        runner = tmp / "log-boundaries.rs"
        runner.write_text("""
unsafe extern "C" {
    fn unit_log_chain(seed: f64, c: f64, rounds: i32) -> f64;
}

fn main() {
    let snan = f64::from_bits(0x7ff0_0000_0000_0001);
    let tie = -53.0 * std::f64::consts::LN_2;
    let tiebits = tie.to_bits();
    for seed in [
        f64::NEG_INFINITY, -746.0, -40.0, f64::from_bits(tiebits + 1), tie,
        f64::from_bits(tiebits - 1), -1.0, -0.0, 0.0, 710.0,
        f64::INFINITY, f64::NAN, snan,
    ] {
        for coefficient in [
            f64::NEG_INFINITY, -3.25, -0.0, 0.0, f64::from_bits(1),
            0.3, f64::MAX, f64::INFINITY, f64::NAN, snan,
        ] {
            for rounds in [1, 2] {
                let y = unsafe { unit_log_chain(seed, coefficient, rounds) };
                println!("{:016x}", y.to_bits());
            }
        }
    }
}
""")
        base = (TESTS / "libm-unit-log.ll").read_text()
        changed, _ = run([
            opt, f"-load-pass-plugin={plugin}", "-passes=function(rssn-libm-unit-log)",
            "-verify-each", "-S", "-o", "-", "-",
        ], input=base, env=dict(os.environ, RSSN_SPILL_FILTER=""))
        assert "unit.log.cached" in changed, changed
        results = []
        for mode, ir in (("base", base), ("guarded", changed)):
            obj, exe = tmp / f"{mode}.o", tmp / mode
            run([llc, "-O3", "-verify-machineinstrs", "-filetype=obj",
                 "-o", str(obj), "-"], input=ir)
            run([rustc, "--edition=2024", "-O", f"-Clink-arg={obj}",
                 str(runner), "-o", str(exe)])
            results.append(run([str(exe)])[0].splitlines())
        assert len(results[0]) == len(results[1]) == 260
        mismatch = next(((i, a, b) for i, (a, b) in enumerate(zip(*results)) if a != b), None)
        assert mismatch is None, mismatch


def libm_tail_runtime_test(opt, llc, rustc, plugin):
    with tempfile.TemporaryDirectory(prefix="rssn-libm-tail-") as tmp:
        tmp = Path(tmp)
        runner = tmp / "tail-boundaries.rs"
        runner.write_text("""
unsafe extern "C" {
    fn tail_existing(seed: f64, c: f64, scale: f64, rounds: i32) -> f64;
    fn tail_zero_numerator(seed: f64, c: f64, den: f64, rounds: i32) -> f64;
    fn tail_shifted(seed: f64, c: f64, den: f64, scale: f64, shift: f64, rounds: i32) -> f64;
}

fn main() {
    // LLVM can choose a different NaN payload if register allocation changes
    // the operand placement for commutative operations.
    let bits = |x: f64| if x.is_nan() { f64::NAN.to_bits() } else { x.to_bits() };
    let edge = 746.0_f64.ln();
    let snan = f64::from_bits(0x7ff0_0000_0000_0001);
    for u in [
        f64::NEG_INFINITY, -746.0, -1.0, 0.0, 745.0_f64.ln(),
        f64::from_bits(edge.to_bits() - 1), edge, f64::from_bits(edge.to_bits() + 1),
        14.0, 700.0, 710.0, f64::INFINITY, snan,
    ] {
        for c in [
            f64::NEG_INFINITY, -2.25, -0.0, 0.0, 0.25, 2.0,
            f64::MAX, f64::INFINITY, f64::NAN, snan,
        ] {
            for v in [
                -1.0, -0.0, 0.0, f64::from_bits(1), f64::MIN_POSITIVE,
                0.001, 1.0, f64::MAX, f64::INFINITY, f64::NAN,
            ] {
                for n in [1, 2] {
                    let a = unsafe { tail_existing(u, c, v, n) };
                    let b = unsafe { tail_zero_numerator(u, c, v, n) };
                    println!("{:016x} {:016x}", bits(a), bits(b));
                }
            }
        }
    }
    for (u, c, a) in [
        (edge, 0.0, 1.0), (edge + 0.25, -1.0, 1.0),
        (20.0, 1.0, 0.001), (30.0, 1.0e12, 1.0),
        (700.0, 1.0e300, 1.0), (700.0, f64::MAX, f64::MAX),
        (710.0, -0.0, 1.0), (f64::NAN, -1.0, 1.0),
        (12.0, f64::NEG_INFINITY, 1.0), (40.0, -1.0, f64::INFINITY),
    ] {
        for g in [-1.0e308, -1000.0, -746.0, -0.25, -0.0, 0.0, 0.25,
                  746.0, 1.0e308, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, snan] {
            for den in [-1.0, 0.0, 1.0] {
                for n in [1, 2] {
                    println!("{:016x}", bits(unsafe { tail_shifted(u, c, den, a, g, n) }));
                }
            }
        }
    }
}
""")
        base = shifted_tail((TESTS / "libm-tail-cutoff.ll").read_text())
        changed, _ = run([
            opt, f"-load-pass-plugin={plugin}", "-passes=function(rssn-libm-tail-cutoff)",
            "-verify-each", "-S", "-o", "-", "-",
        ], input=base, env=dict(os.environ, RSSN_SPILL_FILTER=""))
        assert changed.count("exp.tail.cutoff =") == 3, changed
        results = []
        for mode, ir in (("base", base), ("cutoff", changed)):
            obj, exe = tmp / f"{mode}.o", tmp / mode
            run([llc, "-O3", "-verify-machineinstrs", "-filetype=obj",
                 "-o", str(obj), "-"], input=ir)
            run([rustc, "--edition=2024", "-O", f"-Clink-arg={obj}",
                 str(runner), "-o", str(exe)])
            results.append(run([str(exe)])[0].splitlines())
        assert len(results[0]) == len(results[1]) == 3380
        mismatch = next(((i, a, b) for i, (a, b) in enumerate(zip(*results)) if a != b), None)
        assert mismatch is None, mismatch


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("plugin", type=Path)
    parser.add_argument("--opt", default="opt-23")
    parser.add_argument("--llc", default="llc-23")
    parser.add_argument("--rustc", default="rustc")
    parser.add_argument("--tests", choices=("all", "ir", "rustc"), default="all")
    args = parser.parse_args()
    plugin_path = str(args.plugin.resolve())
    if args.tests in ("all", "ir"):
        version, _ = run([args.llc, "--version"])
        backends = set(re.findall(r"^\s+([\w-]+)\s+-", version, re.M))
        portable_target_tests(args.opt, args.llc, plugin_path, backends)
        if "x86-64" in backends:
            ir_tests(args.opt, args.llc, plugin_path)
            edge_restore_tests(args.llc, plugin_path)
    if args.tests in ("all", "rustc"):
        version, _ = run([args.rustc, "-vV"])
        linux_x86 = "host: x86_64-unknown-linux-gnu" in version.splitlines()
        rustc_tests(args.rustc, plugin_path, linux_x86)
        if linux_x86:
            reuse_runtime_test(args.opt, args.llc, args.rustc, plugin_path)
            endpoint_runtime_test(args.opt, args.llc, args.rustc, plugin_path)
            pow_unit_runtime_test(args.opt, args.llc, args.rustc, plugin_path)
            unit_log_runtime_test(args.opt, args.llc, args.rustc, plugin_path)
            libm_tail_runtime_test(args.opt, args.llc, args.rustc, plugin_path)
