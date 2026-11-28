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


def ffi(build, kind):
    shared, static, imports = library_names(build.target)
    types, files = {
        "shared": ("cdylib", [shared, *imports]),
        "static": ("staticlib", [static]),
        "both": ("staticlib,cdylib", [shared, *imports, static]),
    }[kind]
    # One rustc invocation when both C artifacts are requested. The shared-only
    # consumers never pay to assemble or copy the static archive.
    profile = None
    overrides = None
    if kind == "both":
        profile = "release" if build.release else "ffi"
        flags = ["--profile", profile, "--all-features"]
        if not build.release:
            overrides = {
                "CARGO_PROFILE_FFI_LTO": build.lto,
                "CARGO_PROFILE_FFI_CODEGEN_UNITS": build.codegen_units,
                "CARGO_PROFILE_FFI_INCREMENTAL": "false",
            }
    else:
        flags = build.profile(["--all-features"])
    late_ffi = (build.root / "src/lib.rs", *files) if build.standalone_llvm else None
    build.cargo(["rustc", "--lib", *flags, "--crate-type", types], env=overrides, late_ffi=late_ffi)
    build.copy_artifacts("ffi", build.cargo_output(profile), files)


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


def package_ffi(build, args, *, artifacts=None, headers=None, from_binaries=False):
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
    sources = build.root
    version = options.version or tomllib.loads((sources / "Cargo.toml").read_text(encoding="utf-8"))["package"]["version"]
    suffix = "-llvm" if build.settings["build"] == ["llvm"] else ""
    name = options.dist_name or f"rssn-{version}-{build.target}{suffix}"
    if name in (".", "..") or not re.fullmatch(r"[A-Za-z0-9_.+-]+", name):
        raise ValueError(f"invalid distribution name: {name!r}")
    if artifacts is None:
        artifacts, *headers = build.extra
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


def legacy_ffi_dev(build, args):
    parser = argparse.ArgumentParser(description="Build the development FFI package for legacy clients")
    parser.add_argument("--profile", choices=("dev", "test"), required=True)
    parser.add_argument("--package", action="store_true")
    parser.add_argument("--dist-dir")
    parser.add_argument("--dist-name")
    parser.add_argument("--version")
    options, cargo_args = parser.parse_known_args(args)
    build.cargo(["rustc", "--lib", "--profile", options.profile,
                 "--crate-type", "staticlib,cdylib", *cargo_args])
    source = build.cargo_output("debug")
    root = build.workspace / Path(build.env.get("CARGO_TARGET_DIR", "target"))
    if build.settings["target"]:
        root /= build.target
    root /= "debug"
    root.mkdir(parents=True, exist_ok=True)
    shared, static, imports = library_names(build.target)
    names = (shared, static, *imports)
    with tempfile.TemporaryDirectory(dir=root) as temp:
        for name in names:
            shutil.copy2(source / name, Path(temp) / name)
        for name in names:
            (Path(temp) / name).replace(root / name)
    if options.package:
        package_args = []
        for name in ("dist_dir", "dist_name", "version"):
            if value := getattr(options, name):
                package_args.extend(("--" + name.replace("_", "-"), value))
        package_ffi(build, package_args, artifacts=root,
                    headers=[build.root / "rssn.h", build.root / "rssn.hpp"])


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


def binding_suites(build, root, env, *, require_all, has_r):
    print(f"library: {env.get('RSSN_LIB') or 'automatic discovery'}", flush=True)
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

    readmes = ["python", "julia"]
    if has_r:
        build.run(["Rscript", "bindings/tests/run_conformance.R"], cwd=root, env=env)
        with tempfile.TemporaryDirectory(dir=build.scratch) as temporary:
            names_file = Path(temporary) / "operators.txt"
            with names_file.open("wb") as output:
                build.run([
                    python, "-c", 'import rssn; print("\\n".join(op["name"] for op in rssn.manifest()["ops"]))',
                ], cwd=python_root, env=env, stdout=output)
            build.run(["Rscript", "bindings/r/tests/test-rssn.R"], cwd=root,
                      env={**env, "RSSN_CONFORMANCE_OP_NAMES": str(names_file)})
        for script in ("bindings/r/tests/test-json.R", "bindings/r/tests/test-surface.R",
                       "bindings/examples/r/exact_tensor.R"):
            build.run(["Rscript", script], cwd=root, env=env)
        readmes.append("r")
    else:
        print("R suite skipped. Install R and the rssn package to include it.", flush=True)

    cmake_root = build.state / ("cpp-release" if require_all else "cpp-local")
    cmake_options = ["-U", "RSSN_LIBRARY"] if require_all else []
    build.run(["cmake", "-S", "bindings/cpp", "-B", cmake_root, *cmake_options], cwd=root, env=env)
    build.run(["cmake", "--build", cmake_root, "--config", "Release"], cwd=root, env=env)
    build.run(["ctest", "--test-dir", cmake_root, "-C", "Release", "--output-on-failure"], cwd=root, env=env)
    build.python("bindings/tests/run_readme_examples.py", *readmes, cwd=root, env=env)
    print("all four language suites passed" if has_r else "Python, Julia and C++ suites passed (R skipped)")


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
    binding_suites(build, build.root, env, require_all=True, has_r=True)


def legacy_bindings(build, args):
    build.assert_native()
    parser = argparse.ArgumentParser(description="Run binding suites from a source checkout")
    parser.add_argument("--require-all", action="store_true")
    options = parser.parse_args(args)
    root = build.checkout("rssn")
    library = build.env.get("RSSN_LIB")
    if options.require_all and (not library or not (root / library).is_file()):
        parser.error("set RSSN_LIB to the shared library file for all four bindings")
    env = {"RSSN_LIB": str(root / library)} if library and (root / library).exists() else {}
    rscript = shutil.which("Rscript", path=build.env.get("PATH"))
    has_r = bool(rscript and build.run([
        rscript, "-e", 'quit(status = if (requireNamespace("rssn", quietly = TRUE)) 0L else 1L)',
    ], cwd=root, env=env, check=False).returncode == 0)
    if options.require_all and not has_r:
        parser.error("install R and the rssn package from bindings/r for full conformance")
    binding_suites(build, root, env, require_all=options.require_all, has_r=has_r)


def manifest(build):
    result = build.cargo(["run", "--bin", "signature-audit", "--all-features", "--", "--json"], stdout=-1)
    path = build.outputs["manifest.json"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(result.stdout)


def refresh_bindings(build, args):
    build.assert_native()
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
        build.python("bindings/tools/generate.py", *args)


def dispatch(build, workflow, args):
    if workflow == "llvm-build":
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
        ffi(build, build.settings["build"][0])
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
    elif workflow == "ffi-legacy-dev":
        legacy_ffi_dev(build, args)
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
        build.cargo(["test", *options, *args])
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
    elif workflow == "bindings-legacy":
        legacy_bindings(build, args)
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
    elif workflow == "refresh-bindings":
        refresh_bindings(build, args)
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
