"""RSSN artifact, test and API maintenance workflows."""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import tomllib


GENERATED = [
    "rssn/src/rules/numerical_abi/contract_generated.rs",
    "rssn/src/rules/numerical_abi/structured_generated.rs",
    "rssn/src/rules/numerical_abi/generated.rs",
    "rssn/src/rules/numerical_abi/codec_generated.rs",
]


def format_sources(build, check):
    command = ["rustup", "run", "nightly", "cargo", "fmt", "--all"]
    if check:
        build.run([*command, "--check"])
    else:
        with build.refresh(("rssn/**/*.rs",)):
            build.run(command)


TIMING_EXAMPLES = {
    "backends": "backend_timing",
    "engine-timing": "engine_timing",
}

# A trailing :: selects a suite. Other filters select one exact libtest case.
BENCHMARK_TESTS = {
    "extract": (None, "graph::extract::tests::extraction_hot_path_timing_report"),
    "eval-cost": (None, "graph::extract::tests::eval_cost_benchmark_report"),
    "cp-als": (None, "kernels::tensor::tests::test_cp_cost_20x20x20_rank3"),
    "groebner-elimination": (None, "rules::poly::groebner::tests::hard_examples_report_basis_sizes_and_times"),
    "fft-plans": ("kernels", "transforms::composite_and_power_of_two_plans_return_finite_values"),
    "parallel-kernels": ("kernels", "par_sites::"),
    "parallel-sim": ("sim", "par_sites::"),
}


def run_cargo_target(build, workflow, args):
    parser = argparse.ArgumentParser(description=f"Run an RSSN {workflow}")
    if workflow == "benchmark":
        parser.epilog = "Timing workloads: " + ", ".join((*TIMING_EXAMPLES, *BENCHMARK_TESTS)) + "."
    parser.add_argument("name")
    parser.add_argument("--profile")
    # The inner separator passes application or benchmark flags through as-is.
    separator = args.index("--") if "--" in args else len(args)
    options, target_args = parser.parse_known_args(args[:separator])
    target_args.extend(args[separator + 1:])
    build.assert_runnable()
    name = options.name
    profile = options.profile
    env = None
    if workflow == "example":
        command = ["run", "--example", name]
    elif name in TIMING_EXAMPLES:
        command = ["run", "--example", TIMING_EXAMPLES[name]]
        profile = profile or "tools"
    elif name in BENCHMARK_TESTS:
        target, filter_name = BENCHMARK_TESTS[name]
        exact = [] if filter_name.endswith("::") else ["--exact"]
        command = ["test", *(["--test", target] if target else ["--lib"])]
        target_args = [filter_name, *exact, "--nocapture", *target_args]
        env = {"RUST_TEST_THREADS": build.env.get("RUST_TEST_THREADS", "1")}
    else:
        command = ["bench", "--bench", name]
    profile_args = ["--profile", profile] if profile else []
    build.cargo([*command, *profile_args, *(["--", *target_args] if target_args else [])], env=env)


def benchmark_ci(build, args):
    parser = argparse.ArgumentParser(description="Run the RSSN CI benchmark suite")
    parser.add_argument("--full", action="store_true", help="collect Criterion and timing reports")
    options = parser.parse_args(args)
    engine_args = ["--noplot"] if options.full else ["--test"]
    run_cargo_target(build, "benchmark", ["--profile", "tools", "engine", "--", *engine_args])
    if options.full:
        for name in (*TIMING_EXAMPLES, *BENCHMARK_TESTS):
            run_cargo_target(build, "benchmark", [name])


def rustdoc(build):
    build.python("scripts/check_doc_editions.py")
    build.cargo(["doc", "--workspace", "--no-deps", "--all-features"])


def rust_examples(build):
    build.assert_native()
    build.cargo(["build", "--examples"])
    directory = build.cargo_output("debug") / "examples"
    suffix = ".exe" if "windows" in build.target else ""
    for source in sorted((build.root / "examples").glob("*.rs")):
        if source.stem not in TIMING_EXAMPLES.values():
            build.run([directory / (source.stem + suffix)])


def generate_contracts(build, mode):
    for family in ("contracts", "structured", "bindings", "codecs"):
        build.tool("rssn-api-gen", f"--{family}-{mode}", project="rssn")


def library_names(target):
    if "windows-msvc" in target:
        return "rssn.dll", "rssn.lib", ["rssn.dll.lib"]
    if "windows-gnu" in target:
        return "rssn.dll", "librssn.a", ["librssn.dll.a"]
    if "apple" in target:
        return "librssn.dylib", "librssn.a", []
    return "librssn.so", "librssn.a", []


def ffi(build, kind, profile=None):
    shared, static, imports = library_names(build.target)
    types, files = {
        "shared": ("cdylib", [shared, *imports]),
        "static": ("staticlib", [static]),
        "both": ("staticlib,cdylib", [shared, *imports, static]),
    }[kind]
    # One rustc invocation when both C artifacts are requested. The shared-only
    # consumers never pay to assemble or copy the static archive.
    overrides = None
    if kind == "both":
        profile = profile or ("release" if build.release else "ffi")
        flags = ["--profile", profile, "--all-features"]
        if profile == "ffi":
            overrides = {
                "CARGO_PROFILE_FFI_LTO": build.lto,
                "CARGO_PROFILE_FFI_CODEGEN_UNITS": build.codegen_units,
                "CARGO_PROFILE_FFI_INCREMENTAL": "false",
            }
    else:
        flags = build.profile(["--all-features"])
    late_ffi = (build.root / "src/lib.rs", *files) if build.standalone_llvm else None
    build.cargo(["rustc", "--lib", *flags, "--crate-type", types], env=overrides, late_ffi=late_ffi)
    build.copy_artifacts("ffi", build.cargo_output("debug" if profile == "dev" else profile), files)


# Portable SDK sources and their examples. Native binaries are already
# shipped at the archive root, so local build and install output paths stay out.
SDK_PACKAGE_FILES = (
    "bindings/README.md",
    "bindings/api/*.md",
    "bindings/api/*.txt",
    "bindings/cpp/CMakeLists.txt",
    "bindings/cpp/README.md",
    "bindings/cpp/include/**/*.hpp",
    "bindings/cpp/tests/*.cpp",
    "bindings/examples/README.md",
    "bindings/examples/cpp/*.cpp",
    "bindings/examples/julia/*.jl",
    "bindings/examples/python/*.py",
    "bindings/examples/r/*.R",
    "bindings/julia/Project.toml",
    "bindings/julia/README.md",
    "bindings/julia/src/*.jl",
    "bindings/julia/test/*.jl",
    "bindings/python/pyproject.toml",
    "bindings/python/README.md",
    "bindings/python/setup.py",
    "bindings/python/rssn/*.py",
    "bindings/python/rssn/py.typed",
    "bindings/python/tests/*.py",
    "bindings/r/DESCRIPTION",
    "bindings/r/NAMESPACE",
    "bindings/r/README.md",
    "bindings/r/R/*.R",
    "bindings/r/inst/examples/*.R",
    "bindings/r/man/*.Rd",
    "bindings/r/src/*.c",
    "bindings/r/tests/*.R",
    "bindings/tests/cases.tsv",
    "bindings/tests/*.cpp",
    "bindings/tests/*.py",
    "bindings/tests/*.R",
    "bindings/tools/*.py",
)


def package_ffi(build, args, *, from_binaries=False):
    parser = argparse.ArgumentParser(description="Package the Bazel-built RSSN C distribution")
    parser.add_argument("--dist-name")
    parser.add_argument("--dist-dir", type=Path, default=build.workspace / "dist")
    parser.add_argument("--version")
    if from_binaries:
        parser.add_argument("--native-dir", type=Path, required=True)
        parser.add_argument("--headers-dir", type=Path, required=True)
    options = parser.parse_args(args)
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    out = caller / options.dist_dir
    if from_binaries:
        artifacts = (caller / options.native_dir).resolve()
        headers = [(caller / options.headers_dir / name).resolve() for name in ("rssn.h", "rssn.hpp")]
    else:
        artifacts, *headers = build.extra
    sources = build.root
    version = options.version or tomllib.loads((sources / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    suffix = "-dev" if build.settings["build"] == ["dev"] else ""
    name = options.dist_name or f"rssn-{version}-{build.target}{suffix}"
    if name in (".", "..") or not re.fullmatch(r"[A-Za-z0-9_.+-]+", name):
        raise ValueError(f"invalid distribution name: {name!r}")
    shared, static, imports = library_names(build.target)
    entries = {library: artifacts / library for library in [shared, static, *imports]}
    entries.update({header.name: header for header in headers})
    names = (
        "LICENSE", "ARCHITECTURE.md", "ATTRIBUTIONS.md",
        "CHANGELOG.md", "CODE_OF_CONDUCT.md", "CONTRIBUTING.md", "SECURITY.md",
        "examples/README.md",
    )
    entries.update({entry: sources / entry for entry in names})
    entries.update({p.name: p for p in sources.glob("scripts/ffi-package/README*.md")})
    for pattern in ("docs/*.md", "examples/*.rs", "examples/c/*.c", *SDK_PACKAGE_FILES):
        entries.update({p.relative_to(sources).as_posix(): p for p in sources.glob(pattern) if p.is_file()})
    checkout = build.checkout("rssn")
    source_ref = "v0.3"
    if (checkout / ".git").exists() and shutil.which("git", path=build.env.get("PATH")):
        revision = build.run(["git", "rev-parse", "--verify", "HEAD"], cwd=checkout,
                             capture_output=True, text=True, check=False)
        if revision.returncode == 0:
            source_ref = revision.stdout.strip()
    repo = "https://github.com/Apich-Organization/rssn"
    out.mkdir(parents=True, exist_ok=True)
    archive = out / f"{name}.tar.gz"
    with tempfile.TemporaryDirectory(dir=out) as temporary:
        stage = Path(temporary) / name
        for entry, source in entries.items():
            destination = stage / entry
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.suffix == ".md":
                text = source.read_text(encoding="utf-8")
                text = text.replace("@RSSN_SOURCE@", f"{repo}/tree/{source_ref}")
                text = re.sub(r"(\]\()(?:\.\./)?((?:src|tools)/[^\s)]+)",
                              lambda m: f"{m[1]}{repo}/blob/{source_ref}/{m[2]}", text)
                text = text.replace("](README.md#testing)", f"]({repo}/blob/{source_ref}/README.md#testing)")
                destination.write_text(text, encoding="utf-8", newline="\n")
            else:
                shutil.copyfile(source, destination)

        def normalize(info):
            info.uid = info.gid = info.mtime = 0
            info.uname = info.gname = ""
            info.mode = 0o755 if info.isdir() else 0o644
            return info

        compressed = Path(temporary) / archive.name
        with compressed.open("wb") as output, gzip.GzipFile(fileobj=output, filename="", mode="wb", mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w", format=tarfile.PAX_FORMAT) as tar:
                tar.add(stage, arcname=name, filter=normalize)
        with compressed.open("rb") as output:
            digest = hashlib.file_digest(output, "sha256").hexdigest()
        checksum = compressed.with_suffix(".gz.sha256")
        checksum.write_text(f"{digest}  {archive.name}\n", encoding="ascii", newline="\n")
        if from_binaries or build.settings["build"] == ["validate"]:
            verify_ffi_package(build, compressed)
        compressed.replace(archive)
        checksum.replace(archive.with_suffix(".gz.sha256"))
    print(archive)


def verify_ffi_clients(build, *cmake_options):
    build.assert_native()
    generator = {
        "i686-pc-windows-msvc": ["-A", "Win32"],
        "x86_64-pc-windows-msvc": ["-A", "x64"],
        "aarch64-pc-windows-msvc": ["-A", "ARM64"],
    }.get(build.target, [])
    with tempfile.TemporaryDirectory(prefix="rssn-ffi-check-", dir=build.scratch) as temp:
        build.run(["cmake", "-S", "tests/ffi/package", "-B", temp,
                   "-DCMAKE_BUILD_TYPE=Release", *generator, *cmake_options])
        # Avoid linking several copies of the static library at once.
        build.run(["cmake", "--build", temp, "--config", "Release", "--parallel", "1"])
        build.run(["ctest", "--test-dir", temp, "--build-config", "Release", "--output-on-failure"])


def verify_ffi_package(build, archive):
    verify_ffi_clients(build, f"-DRSSN_ARCHIVE={archive}")


def test_ffi_clients(build, directory, linkage):
    verify_ffi_clients(build, f"-DRSSN_LIBRARY_DIR={directory}", f"-DRSSN_LINKAGE={linkage}")


def plugin_build(build):
    build.need_rust()
    sdk = build.config["sdk"]
    sdk_version = build.config["sdk_info"].get("version")
    if not sdk or not sdk_version:
        raise ValueError("set --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config for the matching LLVM SDK")
    rust_llvm = build.config["llvm_version"]
    if sdk_version != rust_llvm:
        raise ValueError(f"rustc uses LLVM {rust_llvm}; {sdk} reports {sdk_version}")
    source = build.root / "tools/llvm-optimizer"
    cache = build.state / "cmake"
    host_flags = []
    if import_library := build.config["env"].get("RSSN_LLVM_IMPORT_LIBRARY"):
        host_flags.append(f"-DRSSN_LLVM_IMPORT_LIBRARY={import_library}")
    build.run([
        "cmake", "-S", source, "-B", cache, "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_TESTING=OFF", f"-DRSSN_RUSTC={build.rustc_exe}", f"-DRSSN_LLVM_CONFIG={sdk}",
        *host_flags,
    ])
    build.run(["cmake", "--build", cache, "--config", "Release", "--parallel", build.jobs])
    module, = build.outputs
    for directory in (cache / "Release", cache):
        for suffix in (Path(module).suffix, ".dylib", ".so", ".dll"):
            for prefix in ("", "lib"):
                file = directory / (prefix + "RssnLLVMOptimizer" + suffix)
                if file.is_file():
                    build.write_output(module, file)
                    return
    raise ValueError(f"CMake built no RssnLLVMOptimizer loadable module in {cache}")


def llvm_benchmark(build, args):
    build.assert_native()
    build.need_rust()
    bench = build.root / "tools/llvm-optimizer/bench"
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    llvm_env = [f"--llvm-env={name}={value}" for name, value in build.llvm_loader_env().items()]
    build.python(bench / "run.py", build.plugin,
                 "--work-dir", build.state / "runs", "--opt", build.sdk_tool("opt"),
                 "--llc", build.sdk_tool("llc"), "--rustc", build.rustc_exe,
                 "--cargo", build.cargo_exe, *llvm_env, *args, cwd=caller)


C_BENCHMARKS = ("zstd", "zstd-fse", "win-appell")


def llvm_benchmark_environment(parser, entries):
    overrides = {}
    for entry in entries:
        key, equals, value = entry.partition("=")
        if not equals or not re.fullmatch(r"[A-Z][A-Z0-9_]*", key):
            parser.error("--env expects NAME=VALUE")
        overrides[key] = value
    return overrides


def llvm_zstd_sources(build, bench):
    with (build.root / "Cargo.lock").open("rb") as lock:
        locked = next(pkg["version"] for pkg in tomllib.load(lock)["package"] if pkg["name"] == "zstd-sys")
    # Use Cargo's configured registry or source replacement for the locked
    # research package, without resolving the full RSSN workspace.
    build.need_rust()
    manifest = bench / "zstd-io/Cargo.toml"
    result = build.run([build.cargo_exe, "metadata", "--locked", "--format-version", "1",
                        "--manifest-path", manifest, "--filter-platform", build.target],
                       env={"RUSTC": str(build.rustc_exe)}, capture_output=True, text=True)
    package = next(pkg for pkg in json.loads(result.stdout)["packages"] if pkg["name"] == "zstd-sys")
    if package["version"] != locked:
        raise ValueError(f"LLVM zstd study pins {package['version']}, but the RSSN workspace pins {locked}")
    return Path(package["manifest_path"]).parent / "zstd/lib"


def llvm_zstd_benchmark(build, options, bench, library, out, tool_env, clang):
    if options.case == "zstd-fse":
        unit = library / "common/fse_decompress.c"
        optimization = "rssn-scalar-load-pack"
        function = "FSE_decompress_wksp_bmi2"
    elif options.study == "double-fast":
        unit = library / "compress/zstd_double_fast.c"
        optimization = "rssn-switch-load-reuse"
        function = "ZSTD_fillDoubleHashTable"
    else:
        unit = library / "decompress/zstd_decompress_block.c"
        optimization = "rssn-scalar-load-pack"
        function = "ZSTD_decompressSequences"
    if not unit.is_file():
        raise ValueError(f"{library} is missing the zstd C source {unit.name}")

    # The selected unit is compiled to O3 IR once. Both versions share the
    # same archive and the exact same harness object.
    include = ["-I", library, "-I", library / "common"]
    exe_suffix = ".exe" if build.target.endswith("-windows-gnu") else ""
    target_flags = [f"--target={build.target}"]
    compile_flags = [*target_flags, "-O3", "-DNDEBUG", f"-march={options.cpu}",
                     "-DZSTD_DISABLE_ASM=1", "-DZSTD_LIB_DEPRECATED=0", *include, *options.cflag]
    base = out / "base.ll"
    transformed = out / "transform.ll"
    build.run([clang, *compile_flags, "-S", "-emit-llvm", unit, "-o", base], env=tool_env)
    build.run([build.sdk_tool("opt"), f"-load-pass-plugin={build.plugin}", f"-passes=function({optimization}),verify",
               f"-pass-remarks-analysis={optimization}", "-S", base, "-o", transformed],
              env={**tool_env, "RSSN_SPILL_MODE": "off", "RSSN_SPILL_FILTER": function})

    # Recreate the archive from exact members so an earlier run cannot leave
    # a second copy of the studied translation unit in the link.
    archive = out / "libcommon.a"
    with tempfile.TemporaryDirectory(dir=out) as staging:
        objects = []
        for source in sorted(path for family in ("common", "compress", "decompress")
                             for path in (library / family).glob("*.c") if path != unit):
            obj = Path(staging) / f"{source.stem}.o"
            namespace = ["-DXXH_NAMESPACE=ZSTD_"] if source.stem == "xxhash" else []
            build.run([clang, *compile_flags, *namespace, "-ffunction-sections", "-fdata-sections",
                       "-c", source, "-o", obj], env=tool_env)
            objects.append(obj)
        new_archive = Path(staging) / archive.name
        build.run([build.sdk_tool("llvm-ar"), "rcs", new_archive, *objects], env=tool_env)
        new_archive.replace(archive)

    harness = out / "harness.o"
    build.run([clang, *target_flags, "-O2", *include, *options.cflag, "-c",
               bench / f"{options.case}.c", "-o", harness], env=tool_env)
    link_flags = ["-no-pie"] if not exe_suffix else []
    for name, ir in (("base", base), ("transform", transformed)):
        obj = out / f"{name}.o"
        build.run([build.sdk_tool("llc"), "-O3", "-verify-machineinstrs", "-filetype=obj", ir, "-o", obj],
                  env=tool_env)
        build.run([clang, *target_flags, *options.cflag, *link_flags, harness, obj, archive,
                   "-o", out / (name + exe_suffix)], env=tool_env)


def llvm_windows_benchmark(build, parser, options, bench, caller, out, tool_env, clang):
    symbols = {
        "appell": "rssn_de_appell",
        "lerch-small": "rssn_de_small",
        "lerch-positive": "rssn_de_positive",
        "lerch-negative": "rssn_de_negative",
    }
    symbol = symbols[options.kernel]
    recipes = options.variant or [
        "base=verify",
        "pack=function(rssn-pow-exponent-pack)",
        "reuse=function(rssn-libm-reuse)",
        "both=function(rssn-libm-reuse,rssn-pow-exponent-pack)",
    ]
    if options.kernel != "appell" and not options.variant:
        parser.error("Lerch runner variants require --variant NAME=PASSES")
    if build.target != "x86_64-pc-windows-msvc":
        parser.error("win-appell requires --define=rust_target=x86_64-pc-windows-msvc")
    variants = {}
    for recipe in recipes:
        name, equals, passes = recipe.partition("=")
        if not equals or not passes or not re.fullmatch(r"[a-z][a-z0-9-]*", name) or name in variants:
            parser.error("--variant expects a unique NAME=PASSES")
        variants[name] = passes
    if variants.get("base") != "verify" or len(variants) not in (2, 4):
        parser.error("provide two or four variants including base=verify")

    lld_link = options.lld_link or build.sdk_tool("lld-link")
    target = f"--target={build.target}"
    define = [] if options.kernel == "appell" else [f"-DRSSN_WIN_{options.kernel.upper().replace('-', '_')}=1"]
    work = out / "work"
    work.mkdir(exist_ok=True)
    if options.ir:
        source = (caller / options.ir).resolve(strict=True)
    else:
        source = work / "full.ll"
        llvm_benchmark_build(build, [
            "de", "--out-dir", str(work / "cargo"), f"--rustflag=-Ctarget-cpu={options.cpu}",
            *[f"--env={entry}" for entry in options.env],
            "--", f"--emit=llvm-ir={source},link", "-Zno-link",
        ])
    extract = work / "kernel.ll"
    dependencies = []
    if options.kernel == "lerch-negative":
        dependencies = ["--rfunc=gamma_numerical$|statrs8function5gamma(8ln_gamma|5gamma)$"]
    build.run([build.sdk_tool("llvm-extract"), "-S", "--recursive", f"--func={symbol}",
               *dependencies, source, "-o", extract], env=tool_env)
    runner = work / "runner.obj"
    build.run([clang, target, "-O2", "-ffreestanding", *define, "-c", bench / "de/win_appell.c",
               "-o", runner], env=tool_env)
    fltused = work / "fltused.obj"
    build.run([clang, target, "-x", "c", "-c", "-", "-o", fltused],
              env=tool_env, input=b"int _fltused;\n")

    for library in ("math", "kernel32"):
        build.run([build.sdk_tool("llvm-dlltool"), "-m", "i386:x86-64",
                   "-d", bench / "de" / f"{library}.def", "-l", work / f"{library}.lib"], env=tool_env)
    for name, passes in variants.items():
        ir = out / f"{name}.ll"
        obj = out / f"{name}.obj"
        build.run([build.sdk_tool("opt"), f"-mcpu={options.cpu}", f"-load-pass-plugin={build.plugin}",
                   "-S", "-verify-each", f"-passes={passes}", extract, "-o", ir], env=tool_env)
        build.run([build.sdk_tool("llc"), "-O3", f"-mcpu={options.cpu}", "-verify-machineinstrs", "-filetype=obj",
                   ir, "-o", obj], env=tool_env)
        build.run([lld_link, "/dll", "/noentry", "/nodefaultlib", f"/export:{symbol}",
                   f"/out:{out / (name + '.dll')}", obj, fltused, work / "math.lib"], env=tool_env)
    build.run([lld_link, "/entry:mainCRTStartup", "/subsystem:console", "/nodefaultlib",
               f"/out:{out / 'appell-suite.exe'}", runner, work / "kernel32.lib"], env=tool_env)


def llvm_c_benchmark_build(build, args):
    parser = argparse.ArgumentParser(description="Build paired LLVM C research benchmarks")
    cases = parser.add_subparsers(dest="case", required=True)
    for case in C_BENCHMARKS:
        harness = cases.add_parser(case)
        harness.add_argument("--out-dir", required=True, type=Path)
        harness.add_argument("--env", action="append", default=[], metavar="NAME=VALUE")
        harness.add_argument("--clang", help="Clang command when installed separately from the LLVM SDK")
        harness.add_argument("--cpu", default="x86-64-v3", help="LLVM x86-64 CPU level")
        if case == "win-appell":
            harness.add_argument("--ir", type=Path, help="use existing Win64 fat LTO IR for the DE study")
            harness.add_argument("--kernel", choices=("appell", "lerch-small", "lerch-positive", "lerch-negative"),
                                 default="appell")
            harness.add_argument("--variant", action="append", metavar="NAME=PASSES")
            harness.add_argument("--lld-link", help="COFF linker when installed separately from the LLVM SDK")
        else:
            harness.add_argument("--zstd-lib", type=Path, help="zstd/lib from the locked zstd-sys dependency")
            harness.add_argument("--cflag", action="append", default=[], help="Extra compiler and linker argument")
            if case == "zstd":
                harness.add_argument("--study", choices=("decode", "double-fast"), default="decode")
    options = parser.parse_args(args)
    if not build.config["sdk_info"].get("bindir"):
        parser.error("set --repo_env=RSSN_LLVM_CONFIG=/path/to/llvm-config for the matching LLVM SDK")
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    bench = build.root / "tools/llvm-optimizer/bench"
    out = (caller / options.out_dir).resolve()
    overrides = llvm_benchmark_environment(parser, options.env)
    tool_env = {**build.llvm_loader_env(), **overrides}
    clang = options.clang or build.sdk_tool("clang")
    out.mkdir(parents=True, exist_ok=True)
    if options.case == "win-appell":
        llvm_windows_benchmark(build, parser, options, bench, caller, out, tool_env, clang)
    else:
        native_linux = build.target == build.rust_host and build.target.startswith("x86_64-") and "-linux-" in build.target
        windows_gnu = options.case == "zstd" and build.target == "x86_64-pc-windows-gnu"
        if not (native_linux or windows_gnu):
            parser.error("zstd C studies require x86-64 Linux or a Windows GNU target for the zstd harness")
        library = (caller / options.zstd_lib).resolve() if options.zstd_lib else llvm_zstd_sources(build, bench)
        llvm_zstd_benchmark(build, options, bench, library, out, tool_env, clang)


def llvm_count_exp_build(build, args):
    parser = argparse.ArgumentParser(description="Build a Linux libm call-counting preload probe")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--env", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--mode", choices=("exp", "pow", "log"), default="exp")
    options = parser.parse_args(args)
    build.assert_native()
    if "-linux-" not in build.target:
        parser.error("the libm preload probe requires Linux")
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    out = (caller / options.out_dir).resolve()
    overrides = llvm_benchmark_environment(parser, options.env)
    out.mkdir(parents=True, exist_ok=True)
    mode_define = {"exp": [], "pow": ["-DRSSN_COUNT_POW=1"], "log": ["-DRSSN_COUNT_LOG=1"]}
    probe = {"exp": "count_exp", "pow": "count_pow", "log": "count_exp_log"}[options.mode]
    source = build.root / "tools/llvm-optimizer/bench/de/count_exp.c"
    build.run([overrides.get("CC") or build.env.get("CC") or "cc", "-shared", "-fPIC", "-O2", *mode_define[options.mode],
               source, "-ldl", "-o", out / f"{probe}.so"], env=overrides)


def llvm_benchmark_ci(build, args):
    parser = argparse.ArgumentParser(description="Run the Linux LLVM benchmark suite")
    parser.add_argument("--reports-dir", required=True, type=Path)
    options = parser.parse_args(args)
    build.assert_native()
    if not (build.target.startswith("x86_64-") and "-linux-" in build.target):
        parser.error("the LLVM CI suite requires x86-64 Linux")
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    reports = (caller / options.reports_dir).resolve()
    llvm_benchmark(build, ["--suite", "all", "--pairs", "3", "--work-dir", str(reports)])

    rssn = build.root
    fixtures = rssn / "tests/fixtures"
    hdf5 = fixtures / "hdf5/filters_plugins.h5"
    dictionary = reports / "sample.dict"
    with (fixtures / "hdf5/dense_nested.h5").open("rb") as source:
        source.seek(-65536, os.SEEK_END)
        dictionary.write_bytes(source.read())
    bench = rssn / "tools/llvm-optimizer/bench"
    for name, build_args, workload in (
        ("decode", ["zstd"], ["100", fixtures / "parquet/codec_zstd.parquet", hdf5]),
        ("double-fast", ["zstd", "--study", "double-fast"], ["--compress-min4", "100", hdf5, dictionary]),
        ("fse", ["zstd-fse"], [bench / "fixtures/zstd-fse.hex", "25000"]),
    ):
        out = reports / name
        llvm_c_benchmark_build(build, [*build_args, "--out-dir", str(out)])
        build.python(bench / "compare.py", "pair", out / "base", out / "transform", "--pairs", "3",
                     "--output", out / "results.json", "--", *workload)
    for mode in ("exp", "pow", "log"):
        llvm_count_exp_build(build, ["--mode", mode, "--out-dir", str(reports / "libm")])


def llvm_benchmark_build(build, args):
    if args and args[0] == "count-exp":
        llvm_count_exp_build(build, args[1:])
        return
    parser = argparse.ArgumentParser(description="Build an LLVM research benchmark")
    parser.add_argument("case", choices=("example", "de", "zstd-io", "lanes", "tape", "count-exp"))
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--profile", choices=("dev", "release"), default="release")
    parser.add_argument("--example", default="numeric_kernels")
    parser.add_argument("--source", type=Path, help="older CPU or Tape source")
    parser.add_argument("--lib", action="store_true", help="build the DE library")
    parser.add_argument("--rustflag", action="append", default=[])
    parser.add_argument("--env", action="append", default=[], metavar="NAME=VALUE")
    separator = args.index("--") if "--" in args else len(args)
    options = parser.parse_args(args[:separator])
    rustc_args = args[separator + 1:]
    if options.source and options.case not in ("lanes", "tape"):
        parser.error("--source applies to lanes and tape")
    if options.lib and options.case != "de":
        parser.error("--lib applies to de")
    bench = build.root / "tools/llvm-optimizer/bench"
    caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
    out = (caller / options.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    overrides = llvm_benchmark_environment(parser, options.env)
    env = {"RUSTC_BOOTSTRAP": "1", "RUSTC": str(build.rustc_exe),
           "CARGO_BUILD_BUILD_DIR": str(out),
           "CARGO_ENCODED_RUSTFLAGS": "\x1f".join(options.rustflag), **overrides}
    build.need_rust()
    if options.case in ("lanes", "tape"):
        lane = options.case == "lanes"
        source = out / "kernel.rs"
        script = "extract.py" if lane else "extract_tape.py"
        source_flag = "--cpu-source" if lane else "--tape-source"
        extraction = [source_flag, (caller / options.source).resolve()] if options.source else []
        build.python(bench / "device-lanes" / script, source, *extraction)
        env["LANE_SOURCE" if lane else "TAPE_SOURCE"] = str(source)
        entry = bench / "device-lanes" / ("main.rs" if lane else "tape.rs")
        suffix = ".exe" if "windows" in build.target else ""
        build.run([build.rustc_exe, "--edition=2024", "--target", build.target,
                   "-Copt-level=3", "-Cpanic=abort", *options.rustflag, *rustc_args,
                   entry, "-o", out / (options.case + suffix)], env=env)
        return
    manifest = build.root / "Cargo.toml" if options.case == "example" else bench / options.case / "Cargo.toml"
    if options.case == "de":
        source = out / "kernels.rs"
        build.python(bench / "de/extract.py", source)
        env["RSSN_DE_SOURCE"] = str(source)
        target = ["--lib"] if options.lib else ["--bin", "rssn-de-study"]
    elif options.case == "example":
        target = ["--example", options.example]
    else:
        target = ["--bin", "rssn-io-codec-study"]
    target_flag = ["--target", build.target] if build.settings["target"] else []
    extra = ["--", *rustc_args] if rustc_args else []
    build.run([build.cargo_exe, "rustc", "--locked", "--jobs", build.jobs,
               "--manifest-path", manifest, "--target-dir", out, "--profile",
               options.profile, *target_flag, *target, *extra], env=env)


def llvm_test(build, tests):
    build.assert_native()
    plugin, = build.extra
    # Extracted SDKs can have a non-system location for their shared LLVM library.
    # Put this on the test environment, never the release compiler's environment.
    build.python("tools/llvm-optimizer/tests/test_pass.py", plugin, "--tests", tests,
                 "--opt", build.sdk_tool("opt"), "--llc", build.sdk_tool("llc"),
                 "--rustc", build.rustc_exe, env=build.llvm_loader_env())


def standalone_release_probe(build, root, freestanding):
    source, value = root / "probe.rs", root / "value.rs"
    probe = (
        '#[unsafe(no_mangle)]\n'
        'pub extern "C" fn rssn_pressure_probe(a: f64, b: f64) -> f64 {\n'
        '    core::hint::black_box(a) + include!("value.rs") * core::hint::black_box(b)\n}\n'
    )
    extra = []
    if freestanding:
        # Exercise rustc's actual COFF link and custom linker arguments even
        # when the build host has no redistributable MSVC CRT libraries.
        probe = '#![no_std]\n' + probe + (
            '#[panic_handler]\nfn panic(_: &core::panic::PanicInfo<\'_>) -> ! { loop {} }\n'
            '#[unsafe(no_mangle)]\npub static _fltused: i32 = 0;\n'
        )
        extra = ["--", "-Clink-arg=/NOENTRY", "-Clink-arg=/NODEFAULTLIB"]
    source.write_text(probe, encoding="utf-8")
    value.write_text("2.0", encoding="utf-8")
    args = ["rustc", "--release", "--all-features", "--lib", "--crate-type", "staticlib,cdylib",
            "--message-format=json-render-diagnostics", *extra]
    shared, static, imports = library_names(build.target)
    libraries = [shared, static, *imports]
    expected = set(libraries)
    for wanted_fresh in (False, True, False):
        result = build.cargo(args, cwd=root, late_ffi=(source, *libraries), stdout=subprocess.PIPE, text=True)
        artifacts = [message for line in result.stdout.splitlines()
                     if (message := json.loads(line))["reason"] == "compiler-artifact"]
        if len(artifacts) != 1 or artifacts[0]["fresh"] != wanted_fresh:
            raise ValueError(f"expected one Cargo artifact with fresh={wanted_fresh}: {artifacts}")
        files = {Path(file).name for file in artifacts[0]["filenames"] if Path(file).is_file()}
        if not expected <= files:
            raise ValueError(f"standalone test is missing Cargo artifacts: {expected - files}")
        if wanted_fresh:
            value.write_text("3.0", encoding="utf-8")

    if build.rust_host == build.target:
        library = build.cargo_output() / shared
        build.run([build.config["python"], "-c", (
            "import ctypes, sys; f = ctypes.CDLL(sys.argv[1]).rssn_pressure_probe; "
            "f.argtypes = (ctypes.c_double, ctypes.c_double); f.restype = ctypes.c_double; "
            "assert f(2.0, 4.0) == 14.0"
        ), library], cwd=root)


def release_pass_test(build):
    # The spill-pressure diagnostic is only part of the Linux x86-64 contract.
    check_remark = build.target == "x86_64-unknown-linux-gnu"
    with tempfile.TemporaryDirectory(dir=build.scratch) as temporary:
        root = Path(temporary)
        freestanding = (build.standalone_llvm and build.target.endswith("-windows-msvc") and
                        build.rust_host != build.target)
        panic = '[profile.release]\npanic = "abort"\n' if freestanding else ""
        (root / "Cargo.toml").write_text(
            '[package]\nname = "rssn"\nversion = "0.1.0"\nedition = "2024"\n'
            '[lib]\npath = "probe.rs"\n[workspace]\n' + panic, encoding="utf-8",
        )
        (root / "Cargo.lock").write_text(
            'version = 4\n[[package]]\nname = "rssn"\nversion = "0.1.0"\n', encoding="utf-8",
        )
        if build.standalone_llvm:
            standalone_release_probe(build, root, freestanding)
            return
        (root / "probe.rs").write_text(
            '#[inline(never)] #[unsafe(no_mangle)]\n'
            'pub extern "C" fn rssn_pressure_probe(a: f64, b: f64) -> f64 {\n'
            '    std::hint::black_box(a) + 2.0 * std::hint::black_box(b)\n}\n'
            '#[test] fn pressure_test() { assert_eq!(rssn_pressure_probe(2.0, 4.0), 10.0); }\n', encoding="utf-8",
        )
        for crate_type, lto in (("rlib", "off"), ("cdylib", build.lto), ("staticlib", build.lto)):
            result = build.cargo([
                "rustc", "--release", "--all-features", "--lib", "--crate-type", crate_type, "--",
                "-Cremark=rssn-spill-pressure", "-Cllvm-args=-pass-remarks-analysis=rssn-spill-pressure",
            ], cwd=root, lto=lto, stderr=subprocess.PIPE, text=True, check=False, env={
                "RSSN_SPILL_MODE": "report",
                "RSSN_SPILL_FILTER": "rssn_pressure_probe",
            })
            if result.returncode or (check_remark and "rssn_pressure_probe GP peak" not in result.stderr):
                raise ValueError(f"could not validate {crate_type} with LTO={lto}:\n{result.stderr}")
        if build.target == build.rust_host or build.target_runner:
            build.cargo(["test", "--release"], cwd=root)


def binding_suites(build, root, env):
    print(f"library: {env['RSSN_LIB']}", flush=True)
    python_root = root / "bindings/python"
    python = build.config["python"]
    julia = ["julia", "--project=bindings/julia"]
    build.run([python, "-m", "unittest", "discover", "-s", "tests"], cwd=python_root, env=env)
    build.python("bindings/tests/run_conformance.py", cwd=root, env=env)
    build.run([*julia, "bindings/julia/test/runtests.jl"], cwd=root, env=env)

    python_path = python_root.as_posix()
    if build.env.get("PYTHONPATH"):
        python_path += os.pathsep + build.env["PYTHONPATH"]
    build.python("bindings/examples/python/exact_tensor.py", cwd=root, env={**env, "PYTHONPATH": python_path})
    build.run([*julia, "bindings/examples/julia/exact_tensor.jl"], cwd=root, env=env)

    build.run(["Rscript", "bindings/tests/run_conformance.R"], cwd=root, env=env)
    with tempfile.TemporaryDirectory(dir=build.scratch) as temporary:
        names_file = Path(temporary) / "operators.txt"
        with names_file.open("wb") as output:
            build.run([
                python, "-c", r'import rssn; print("\n".join(op["name"] for op in rssn.manifest()["ops"]))',
            ], cwd=python_root, env=env, stdout=output)
        build.run(["Rscript", "bindings/r/tests/test-rssn.R"], cwd=root,
                  env={**env, "RSSN_CONFORMANCE_OP_NAMES": str(names_file)})
    for script in ("bindings/r/tests/test-json.R", "bindings/r/tests/test-surface.R",
                   "bindings/examples/r/exact_tensor.R"):
        build.run(["Rscript", script], cwd=root, env=env)

    cmake_root = build.state / "cpp-release"
    build.run(["cmake", "-S", "bindings/cpp", "-B", cmake_root, "-U", "RSSN_LIBRARY"], cwd=root, env=env)
    build.run(["cmake", "--build", cmake_root, "--config", "Release"], cwd=root, env=env)
    build.run(["ctest", "--test-dir", cmake_root, "-C", "Release", "--output-on-failure"], cwd=root, env=env)
    build.python("bindings/tests/run_readme_examples.py", "python", "julia", "r", cwd=root, env=env)
    print("all four language suites passed")


def bindings(build, check_only, library=None):
    build.assert_native()
    if library is None:
        directory, = build.extra
        shared, _, _ = library_names(build.target)
        library = directory / shared
    if not library.is_file():
        raise ValueError(f"RSSN shared library does not exist at {library}")
    if check_only:
        build.python("bindings/tools/extract_api.py", "--library", library)
        build.python("bindings/tools/generate.py", "--check")
        return
    r_library = build.state / "r-library"
    r_library.mkdir(parents=True, exist_ok=True)
    env = {"RSSN_LIB": str(library), "R_LIBS_USER": str(r_library)}
    # R writes compiled objects beside its sources during installation.
    with tempfile.TemporaryDirectory(dir=build.scratch) as temporary:
        r_source = Path(temporary) / "rssn"
        shutil.copytree(build.root / "bindings/r", r_source)
        build.run(["R", "CMD", "INSTALL", f"--library={r_library}", r_source], env=env)
    binding_suites(build, build.root, env)


def manifest(build):
    result = build.cargo(["run", "--bin", "signature-audit", "--all-features", "--", "--json"], stdout=-1)
    path = build.outputs["manifest.json"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(result.stdout)


def refresh_bindings(build, args, *, from_manifest=False):
    build.assert_native()
    parser = argparse.ArgumentParser(description="Regenerate the RSSN language bindings")
    parser.add_argument("--check", action="store_true")
    if from_manifest:
        parser.add_argument("--manifest", required=True, type=Path)
    options = parser.parse_args(args)
    if from_manifest:
        caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
        manifest_file = (caller / options.manifest).resolve()
    else:
        (manifest_file,) = build.extra
    # Binding generator chooses these output paths itself. Keep the transient
    # API description in the build mirror, and publish the generated veneers.
    patterns = [
        "rssn/bindings/python/**/*.py", "rssn/bindings/julia/**/*.jl",
        "rssn/bindings/r/**/*.R", "rssn/bindings/r/NAMESPACE", "rssn/bindings/r/**/*.Rd",
        "rssn/bindings/cpp/**/*.hpp",
    ]
    with build.refresh(patterns):
        build.python("bindings/tools/extract_api.py", "--manifest", manifest_file)
        build.python("bindings/tools/generate.py", *(["--check"] if options.check else []))


def dispatch(build, workflow, args):
    if workflow in ("fmt", "fmt-check"):
        format_sources(build, workflow == "fmt-check")
    elif workflow in ("example", "benchmark"):
        run_cargo_target(build, workflow, args)
    elif workflow == "benchmark_ci":
        benchmark_ci(build, args)
    elif workflow in ("signature_audit", "doc_audit"):
        build.assert_runnable()
        tool = "signature-audit" if workflow == "signature_audit" else "doc_gen_check"
        build.tool(tool, *args, binary=True)
    elif workflow in ("docs", "rustdoc"):
        rustdoc(build)
        if workflow == "docs":
            rust_examples(build)
    elif workflow == "rust-examples":
        rust_examples(build)
    elif workflow == "publish-crate":
        parser = argparse.ArgumentParser(description="Publish the checked-out RSSN crate")
        parser.add_argument("--dry-run", action="store_true")
        options = parser.parse_args(args)
        # Cargo uses the checkout's Git metadata and package inclusion rules.
        flags = ["--dry-run"] if options.dry_run else []
        build.cargo(["publish", *flags], cwd=build.checkout("rssn"), target=False)
    elif workflow == "llvm-benchmark":
        llvm_benchmark(build, args)
    elif workflow == "llvm-benchmark-ci":
        llvm_benchmark_ci(build, args)
    elif workflow == "llvm-benchmark-build":
        llvm_benchmark_build(build, args)
    elif workflow == "llvm-c-benchmark-build":
        llvm_c_benchmark_build(build, args)
    elif workflow == "llvm-benchmark-compare":
        caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
        build.python(build.root / "tools/llvm-optimizer/bench/compare.py", *args, cwd=caller)
    elif workflow == "llvm-build":
        plugin_build(build)
    elif workflow in ("llvm-test", "llvm-test-tools", "llvm-test-rustc"):
        llvm_test(build, {"llvm-test": "all", "llvm-test-tools": "ir", "llvm-test-rustc": "rustc"}[workflow])
    elif workflow == "release-pass-test":
        release_pass_test(build)
        if build.outputs:
            result = build.outputs["release-pass.ok"]
            result.parent.mkdir(parents=True, exist_ok=True)
            result.write_text("LLVM release probe passed\n", encoding="utf-8")
    elif workflow == "ffi":
        ffi(build, *build.settings["build"])
    elif workflow in ("package-ffi", "package-ffi-from-binaries"):
        package_ffi(build, args, from_binaries=workflow == "package-ffi-from-binaries")
    elif workflow == "ffi-package-smoke":
        parser = argparse.ArgumentParser(description="Test the C and Python clients of an existing FFI archive")
        parser.add_argument("--archive", type=Path, required=True)
        archive = parser.parse_args(args).archive
        caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
        verify_ffi_package(build, caller / archive)
    elif workflow == "ffi-native-clients":
        directory, = build.extra
        linkage, = build.settings["build"]
        test_ffi_clients(build, directory, linkage)
    elif workflow == "ffi-native-clients-from-library":
        parser = argparse.ArgumentParser(description="Run native FFI client tests against existing libraries")
        parser.add_argument("--library-dir", type=Path, required=True)
        parser.add_argument("--linkage", choices=("shared", "static"), required=True)
        options = parser.parse_args(args)
        caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
        test_ffi_clients(build, (caller / options.library_dir).resolve(), options.linkage)
    elif workflow == "rustlib":
        # Whole-program LTO is deferred for a standalone rlib. Use the ordinary
        # O3 pipeline so this artifact already contains the plugin's machine code.
        build.cargo(["build", "--lib", *build.profile()], lto="off")
        build.write_output("librssn.rlib", build.cargo_output() / "librssn.rlib")
    elif workflow == "check":
        build.cargo(["check", *(args or ["--workspace", "--all-targets"]), *build.profile()])
    elif workflow == "manifest":
        build.assert_native()
        manifest(build)
    elif workflow == "clippy":
        build.cargo(["clippy", "--workspace", "--all-targets", "--all-features", "--", "-D", "warnings"])
    elif workflow == "rust-tests":
        build.assert_runnable()
        options = build.settings["build"]
        if build.release:
            options = [*build.profile(), *[option for option in options if option != "--all-features"]]
        build.rust_workspace_tests(options, args)
    elif workflow == "gpu-test":
        build.assert_native()
        result = build.cargo([
            "test", "--workspace", "--all-features", "--lib", "gpu_", "--", "--nocapture",
        ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        print(result.stdout, end="")
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, result.args)
        if "no wgpu adapter; skipping" in result.stdout or "gpu test device:" not in result.stdout:
            raise ValueError("the GPU test did not execute on a wgpu adapter")
    elif workflow == "generated-sources-check":
        generate_contracts(build, "check")
    elif workflow == "contracts":
        build.assert_native()
        generate_contracts(build, "check")
        build.tool("signature-audit", binary=True)
        build.tool("doc_gen_check", binary=True)
    elif workflow in ("bindings", "bindings-check"):
        bindings(build, workflow == "bindings-check")
    elif workflow == "bindings-from-library":
        parser = argparse.ArgumentParser(description="Check and test bindings against an existing RSSN library")
        parser.add_argument("--library", type=Path, required=True)
        inputs = parser.parse_args(args)
        caller = Path(build.env.get("BUILD_WORKING_DIRECTORY", str(build.workspace)))
        library = (caller / inputs.library).resolve()
        bindings(build, True, library)
        bindings(build, False, library)
    elif workflow == "refresh-contracts":
        with build.refresh(GENERATED):
            generate_contracts(build, "write")
    elif workflow in ("headers-check", "refresh-headers", "headers"):
        build.assert_native()
        command = ["run", "-p", "rssn-api-gen", "--features", "headers", "--bin", "rssn-headers", "--"]
        if workflow == "headers-check":
            build.cargo([*command, "--check"])
        elif workflow == "headers":
            build.cargo([*command, "--write"])
            for name in build.outputs:
                build.write_output(name, build.root / name)
        else:
            with build.refresh(["rssn/rssn.h", "rssn/rssn.hpp"]):
                build.cargo([*command, "--write"])
    elif workflow in ("refresh-bindings", "refresh-bindings-from-manifest"):
        refresh_bindings(build, args, from_manifest=workflow.endswith("-from-manifest"))
    elif workflow == "witness":
        build.assert_native()
        build.cargo(["run", "--profile", "tools", "--bin", "signature-mine", "--", "--check"])
    elif workflow == "refresh-witness":
        build.assert_native()
        with tempfile.TemporaryDirectory(dir=build.scratch) as temporary:
            trace = str(Path(temporary) / "sorts.tsv")
            build.cargo(["test", "--lib"], env={"RSSN_SORT_TRACE": trace})
            with build.refresh(["rssn/signatures.witnessed.json"]):
                build.cargo(["run", "--profile", "tools", "--bin", "signature-mine", "--", trace])
    else:
        raise ValueError(f"unknown RSSN workflow: {workflow}")
