"""Local Bazel action adapter. Cargo remains the Rust compiler scheduler."""

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True

_TOOL_ENVIRONMENT = re.compile(
    r"(?:CARGO_TARGET_[A-Z0-9_]+_(?:LINKER|RUNNER)|"
    r"(?:(?:HOST|TARGET)_)?(?:CC|CXX|AR|RANLIB|CFLAGS|CXXFLAGS|ARFLAGS)(?:_[A-Za-z0-9_-]+)?|"
    r"PKG_CONFIG(?:_[A-Za-z0-9_-]+)?|CMAKE_[A-Z0-9_]+|"
    r"SDKROOT|INCLUDE|LIB|LIBPATH|MACOSX_DEPLOYMENT_TARGET|IPHONEOS_DEPLOYMENT_TARGET)"
)


def same(left, right):
    return right.is_file() and left.read_bytes() == right.read_bytes()


def write_atomic(path, data, mode=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        handle.write(data)
    temporary = Path(handle.name)
    if mode is None:
        mode = stat.S_IMODE(path.stat().st_mode) if path.is_file() else 0o644
    try:
        temporary.chmod(mode)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def lock(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / ".lock").open("a+b") as handle:
        if os.name == "nt":
            import msvcrt

            while True:
                handle.seek(0)
                try:
                    msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as error:
                    if error.errno not in (13, 36):
                        raise
                    time.sleep(0.2)
        else:
            import fcntl

            fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def command(args, cwd, env, *, check=True, **kwargs):
    args = list(args)
    if os.name == "nt" and args:
        head = Path(args[0])
        if not head.suffix and not head.is_absolute():
            args[0] = shutil.which(str(args[0]), path=env.get("PATH") if env else None) or args[0]
    print("+ " + subprocess.list2cmdline([str(arg) for arg in args]), flush=True)
    return subprocess.run(args, cwd=cwd, env=env, check=check, **kwargs)


def rust_cpu_flags(target, baseline):
    if baseline == "distribution":
        if target.startswith("x86_64-"):
            return ["-Ctarget-cpu=x86-64-v3"]
        if target.startswith("aarch64-"):
            return ["-Ctarget-feature=+neon"]
    # Keep the Rust target CPU for the 32-bit x86 ABI and for non-native tools.
    if target.startswith(("i686-", "wasm")):
        return []
    return ["-Ctarget-cpu=generic"]


def release_pass_modes(target):
    modes = {
        "RSSN_SPILL_MODE": "transform",
        "RSSN_SPILL_FILTER": "",
        "RSSN_SPILL_PLAN_WIDTH": "4",
        "RSSN_SPILL_ABI": "off",
        "RSSN_SPILL_LIBM": "off",
        "RSSN_SPILL_LIBM_REUSE": "off",
        "RSSN_LIBM_TAIL_CUTOFF": "transform",
        "RSSN_LIBM_UNIT_LOG": "transform",
        "RSSN_POW_EXPONENT_PACK": "off",
        "RSSN_LOOP_EXIT_COMPRESSION": "transform",
        "RSSN_LOOP_PREHEADER_PRESSURE": "transform",
        "RSSN_FP_ENDPOINT_GUARDS": "transform",
        "RSSN_FP_CLASS_INTERVAL": "transform",
        "RSSN_SCALAR_LOAD_PACK": "transform",
        "RSSN_SWITCH_LOAD_REUSE": "transform",
        "RSSN_FP_REGISTER_BANK": "off",
        "RSSN_POW_UNIT_BASE": "transform",
    }
    if target.startswith("x86_64-") and "-linux-" in target and not target.endswith("gnux32"):
        modes.update(RSSN_SPILL_ABI="transform", RSSN_SPILL_LIBM="transform",
                     RSSN_SPILL_LIBM_REUSE="transform")
    elif target == "x86_64-pc-windows-msvc":
        modes.update(RSSN_SPILL_LIBM="transform", RSSN_SPILL_LIBM_REUSE="transform",
                     RSSN_POW_EXPONENT_PACK="transform", RSSN_FP_REGISTER_BANK="transform")
    return modes


class Build:
    def __init__(self, settings, lookup):
        self.settings = settings
        self.lookup = lookup
        self.config = json.loads(lookup(settings["config"]).read_text(encoding="utf-8"))
        self.workspace = Path(self.config["workspace"])
        self.inputs = {key: lookup(value) for key, value in settings["sources"].items()}
        self.extra = [lookup(value).absolute() for value in settings["extra"]]
        self.outputs = {key: Path(value).absolute() for key, value in settings.get("outputs", {}).items()}
        self.project = settings["project"]
        self.plugin = lookup(settings["plugin"]).absolute() if settings["plugin"] else None
        self.env = os.environ.copy()
        # Avoid hidden Cargo profile, Rust flags, and wrappers from the caller.
        for key in list(self.env):
            if key.startswith("CARGO_PROFILE_") or key in (
                "CARGO_ENCODED_RUSTFLAGS", "RUSTFLAGS", "RUSTDOCFLAGS",
                "RUSTC_WRAPPER", "RUSTC_WORKSPACE_WRAPPER", "CARGO_INCREMENTAL",
                "CARGO_BUILD_TARGET", "CARGO_BUILD_BUILD_DIR", "DEV", "RUSTC_BOOTSTRAP",
                "SOURCE_DATE_EPOCH", "AEVIA_BLESS", "RSSN_SORT_TRACE",
            ) or key.startswith(("VERGEN_", "RSSN_SPILL_", "RSSN_BAZEL_LATE_", "RUSTUP_TOOLCHAIN")) or key.endswith("_RUSTFLAGS"):
                self.env.pop(key)
        self.env.update(self.config["env"])
        self.env["VERGEN_IDEMPOTENT"] = "1"
        self.env["RUST_MIN_STACK"] = "67108864" if self.project == "aevia" else "8388608"
        self.env["CARGO_TERM_COLOR"] = "never"
        self.scratch = Path(os.getenv("TEST_TMPDIR", tempfile.gettempdir()))
        installed = self.config["rust"]
        tool_bin = Path(installed.get("bin", "."))
        suffix = ".exe" if self.config["windows"] else ""
        self.cargo_exe = Path(installed.get("cargo") or (tool_bin / ("cargo" + suffix)))
        self.rustc_exe = Path(installed.get("rustc") or (tool_bin / ("rustc" + suffix)))
        self.rust_host = next(
            (line[6:] for line in installed.get("version", "").splitlines() if line.startswith("host: ")), ""
        )
        self.target = settings["target"] or self.rust_host
        if self.target and not re.fullmatch(r"[0-9A-Za-z_-]+", self.target):
            raise ValueError(f"expected a Rust target triple, got {self.target!r}")
        target_configs = self.config.get("target_env", {})
        if not isinstance(target_configs, dict):
            raise ValueError("RSSN_TARGET_ENV must map Rust target triples to environment objects")
        target_environment = target_configs.get(self.target, {})
        if not isinstance(target_environment, dict):
            raise ValueError(f"target environment for {self.target} must be an object")
        for variable, value in target_environment.items():
            if not _TOOL_ENVIRONMENT.fullmatch(variable) or not isinstance(value, str):
                raise ValueError(f"{variable}: expected a string-valued target toolchain variable")
        if settings["workflow"] != "llvm-build":
            self.env.update(target_environment)
        key_prefix = "CARGO_TARGET_" + self.target.replace("-", "_").upper()
        self.target_runner = self.env.get(key_prefix + "_RUNNER")
        self.release = self.plugin is not None
        self.cpu_baseline = settings.get("cpu_baseline", "portable")
        if self.cpu_baseline not in ("portable", "distribution"):
            raise ValueError("rssn_cpu_baseline must be portable or distribution")
        self.lto = settings["lto"]
        self.codegen_units = settings.get("codegen_units", "16")
        if not self.codegen_units.isdecimal() or int(self.codegen_units) < 1:
            raise ValueError("rssn_codegen_units must be a positive integer")
        if self.lto not in ("thin", "off", "fat"):
            raise ValueError("rssn_lto must be thin, off or fat")
        self.standalone_llvm = settings.get("standalone_llvm", False)
        if self.standalone_llvm and (
            not (self.target.endswith(("-windows-msvc", "-freebsd", "-apple-darwin")) or
                 "-linux-gnu" in self.target or "-linux-musl" in self.target) or
            self.lto != "fat" or
            self.codegen_units != "1" or not self.release
        ):
            raise ValueError("rssn_llvm_standalone requires the LLVM plugin, a Windows MSVC, Linux, macOS or FreeBSD target, fat LTO and one codegen unit")
        threads = settings.get("frontend_threads", "1")
        if not threads.isdecimal() or int(threads) < 1:
            raise ValueError("rssn_frontend_threads must be a positive integer")
        self.frontend_threads = int(threads)
        jobs = settings.get("jobs", "")
        if jobs and (not jobs.isdecimal() or int(jobs) < 1):
            raise ValueError("rssn_jobs must be a positive integer")
        self.jobs = jobs or str(len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1)
        root = Path(self.config["env"].get("BAZEL_SCI_CACHE", self.workspace / ".bazel-cache"))
        tools = dict(self.config["env"])
        for name in ("BAZEL_SCI_CACHE", "PYTHON", "RSSN_SCCACHE", "RSSN_TARGET_ENV"):
            tools.pop(name, None)
        if settings["workflow"] != "llvm-build":
            # Changing SDK selection alone must not duplicate development and
            # header caches. Release Cargo units key on the plugin digest.
            for name in ("CMAKE_GENERATOR", "RSSN_FETCH_LLVM_SDK", "RSSN_LLVM_CONFIG", "RSSN_LLVM_IMPORT_LIBRARY"):
                tools.pop(name, None)
            if target_environment:
                tools["target_environment"] = target_environment
        else:
            tools["llvm"] = [self.config["sdk"], self.config["sdk_info"]]
        key = json.dumps(
            [str(self.workspace), self.project, settings["scope"], installed, tools],
            sort_keys=True,
        )
        self.state = (root if root.is_absolute() else self.workspace / root) / hashlib.sha256(key.encode()).hexdigest()[:20]
        self.mirror = self.state / "src"
        self.root = self.mirror / self.project
        self.target_dir = self.state / "cargo"
        if self.frontend_threads != 1 or self.standalone_llvm:
            self.wrapper_source = self.lookup(settings["wrapper"]).absolute()
            self.wrapper_digest = hashlib.sha256(self.wrapper_source.read_bytes()).hexdigest()[:20]
        if self.frontend_threads != 1:
            self.target_dir /= f"frontend-{self.frontend_threads}-{self.wrapper_digest}"

    def need_rust(self):
        if not self.rust_host or not self.rustc_exe.is_file() or not self.cargo_exe.is_file():
            raise ValueError("Install a Rust toolchain satisfying Cargo.toml and run bazel sync --configure")

    def stage(self):
        # Keep paths and unchanged mtimes stable for Cargo and rustc. Remove deleted
        # inputs to ensure removed modules or fixtures cannot mask a build failure.
        index = self.state / "inputs.json"
        old = set(json.loads(index.read_text(encoding="utf-8"))) if index.is_file() else set()
        for key, source in self.inputs.items():
            destination = self.mirror / key
            mode = 0o755 if source.stat().st_mode & stat.S_IXUSR else 0o644
            if not same(source, destination):
                write_atomic(destination, source.read_bytes(), mode)
            elif os.name != "nt" and stat.S_IMODE(destination.stat().st_mode) != mode:
                destination.chmod(mode)
        for key in old - self.inputs.keys():
            (self.mirror / key).unlink(missing_ok=True)
        if old != self.inputs.keys():
            write_atomic(index, json.dumps(sorted(self.inputs)).encode())

    def compiler_wrapper(self):
        directory = self.state / "wrappers" / self.wrapper_digest
        name = "rustc-wrapper.exe" if self.config["windows"] else "rustc-wrapper"
        wrapper = directory / name
        if not wrapper.is_file():
            directory.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=directory) as temporary:
                output = Path(temporary) / name
                self.run([self.rustc_exe, "--edition=2024", self.wrapper_source, "-o", output])
                output.replace(wrapper)
        return wrapper

    def sdk_tool(self, name):
        bindir = Path(self.config["sdk_info"]["bindir"])
        suffix = ".exe" if self.config["windows"] else ""
        major = self.config["llvm_major"]
        for candidate in (bindir / f"{name}{suffix}", bindir / f"{name}-{major}{suffix}"):
            if candidate.is_file():
                return candidate
        raise ValueError(f"LLVM SDK is missing {name} in {bindir}")

    def llvm_loader_env(self):
        sdk = self.config["sdk_info"]
        if self.config["windows"]:
            return {"PATH": sdk["bindir"] + os.pathsep + self.env["PATH"]}
        loader = "DYLD_LIBRARY_PATH" if "apple" in self.rust_host else "LD_LIBRARY_PATH"
        return {loader: sdk["libdir"]} if sdk.get("libdir") else {}

    def cargo(self, args, *, project=None, env=None, target=True, cwd=None, lto=None, late_ffi=None, **kwargs):
        self.need_rust()
        subcommand, *rest = args
        if late_ffi is not None and (
            not self.standalone_llvm or subcommand != "rustc" or lto not in (None, "fat")
        ):
            raise ValueError("the standalone LLVM pipeline requires cargo rustc with fat LTO")
        local = dict(self.env, CARGO_TARGET_DIR=str(self.target_dir), RUSTC=str(self.rustc_exe))
        # A user-wide build-dir override must not escape the locked build state.
        local["CARGO_BUILD_BUILD_DIR"] = str(self.target_dir)
        cache = self.config["env"].get("RSSN_SCCACHE")
        if self.frontend_threads != 1 or late_ffi is not None:
            local.update(
                RUSTC_WRAPPER=str(self.compiler_wrapper()),
                RSSN_BAZEL_FRONTEND_THREADS=str(self.frontend_threads),
            )
            local.pop("RSSN_BAZEL_SCCACHE", None)
            if cache:
                local["RSSN_BAZEL_SCCACHE"] = cache
        elif cache:
            local["RUSTC_WRAPPER"] = cache
        flags = rust_cpu_flags(self.target, self.cpu_baseline) if target else ["-Ctarget-cpu=generic"]
        if target and "-Ctarget-cpu=x86-64-v3" in flags:
            # cc-rs consumes these for target C/C++ dependencies, independently
            # of host build scripts and proc macros.
            native_flag = "/arch:AVX2" if self.target.endswith("-msvc") else "-march=x86-64-v3"
            for language in ("CFLAGS", "CXXFLAGS"):
                key = f"{language}_{self.target}"
                alternate = f"{language}_{self.target.replace('-', '_')}"
                existing = local.get(key, local.get(alternate, ""))
                local[key] = " ".join(filter(None, (existing, native_flag)))
        if self.release and target:
            modes = release_pass_modes(self.target)
            modes.update({k: v for k, v in (env or {}).items() if k in modes})
            content = self.plugin.read_bytes()
            mode_bytes = json.dumps(modes, sort_keys=True).encode()
            digest = hashlib.sha256(mode_bytes + content).hexdigest()
            plugin = self.state / "plugins" / digest[:20] / self.plugin.name
            if not same(self.plugin, plugin):
                write_atomic(plugin, content)
            if late_ffi is None:
                # Explicit --target keeps target flags off host build scripts and proc macros.
                flags.append(f"-Zllvm-plugins={plugin}")
            else:
                adapter = Path(__file__).resolve()
                llvm = {name: str(self.sdk_tool(name)) for name in ("opt", "llc")}
                source, *libraries = late_ffi
                local.update(
                    RSSN_BAZEL_LATE_PYTHON=self.config["python"],
                    RSSN_BAZEL_LATE_ADAPTER=str(adapter),
                    RSSN_BAZEL_LATE_SOURCE=str(source),
                    RSSN_BAZEL_LATE_LIBRARIES=json.dumps(libraries),
                    RSSN_BAZEL_LATE_PLUGIN=str(plugin),
                    RSSN_BAZEL_LATE_TOOLS=json.dumps(llvm),
                    RSSN_BAZEL_LATE_LOADER=json.dumps(self.llvm_loader_env()),
                )
                # Cargo tracks explicit rustc arguments. Include the external pipeline
                # in the final crate fingerprint without invalidating its dependencies.
                identity = [digest, self.wrapper_digest, self.config["sdk_info"], llvm, libraries]
                fingerprint = hashlib.sha256(adapter.read_bytes() + json.dumps(identity, sort_keys=True).encode()).hexdigest()
                if "--" not in rest:
                    rest.append("--")
                rest.append(f"--rssn-bazel-late={fingerprint}")
            local.update(
                RUSTC_BOOTSTRAP="1",
                CARGO_PROFILE_RELEASE_OPT_LEVEL="3",
                CARGO_PROFILE_RELEASE_CODEGEN_UNITS=self.codegen_units,
                CARGO_PROFILE_RELEASE_LTO=lto or self.lto,
                CARGO_PROFILE_RELEASE_INCREMENTAL="false",
                **modes,
            )
        else:
            # Keep source locations without paying for full debug information
            # in both the Rust test binaries and the interpreter.
            local.update(CARGO_PROFILE_DEV_DEBUG="line-tables-only", CARGO_PROFILE_TEST_DEBUG="line-tables-only")
        local["CARGO_ENCODED_RUSTFLAGS"] = "\x1f".join(flags)
        if env:
            local.update(env)
        if subcommand == "test" and "windows-msvc" in self.target:
            local["CARGO_PROFILE_TEST_DEBUG"] = "0"
        options = ["--locked", "--jobs", self.jobs]
        if target:
            options.extend(["--target", self.target])
        if self.settings["timings"]:
            options.append("--timings")
        return command([self.cargo_exe, subcommand, *options, *rest], cwd or (self.mirror / project if project else self.root), local, **kwargs)

    def rust_workspace_tests(self, options, arguments, *, env=None):
        # Test settings are encoded in Bazel's workflow, while --test_arg
        # arrives separately. The bounded mode accepts either entrypoint.
        flags = [*options, *arguments]
        batches = [flag for flag in flags if flag.startswith("--batches=")]
        if not batches:
            self.cargo(["test", *flags], env=env)
            return
        if len(batches) != 1:
            raise ValueError("Pass --batches once")
        batch = batches[0]
        flags.remove(batch)
        self.bounded_rust_tests(flags, int(batch.removeprefix("--batches=")))

    def bounded_rust_tests(self, options, batch_size):
        """Run the full workspace with bounded test-executable storage."""
        if options != ["--workspace", "--all-features"] or self.release:
            raise ValueError("--batches is available on //:rust_all_features_test in the test profile")
        if not 1 <= batch_size <= 16:
            raise ValueError("--batches must be between 1 and 16")
        self.need_rust()
        graph = self.run(
            [self.cargo_exe, "metadata", "--locked", "--no-deps", "--format-version", "1"],
            capture_output=True, text=True, env={"RUSTC": str(self.rustc_exe)},
        )
        metadata = json.loads(graph.stdout)
        members = set(metadata["workspace_members"])
        targets = [target for package in metadata["packages"] if package["id"] in members
                   for target in package["targets"]]
        integration = sorted({target["name"] for target in targets
                              if target["kind"] == ["test"] and target["test"]})
        examples = sorted({target["name"] for target in targets if target["kind"] == ["example"]})
        unit_binaries = sorted({target["name"] for target in targets
                                if target["kind"] == ["bin"] and target["test"]})
        successful = True
        target_dir = self.target_dir.resolve()

        def run(subcommand, flags, *, discard_example=False):
            nonlocal successful
            result = self.cargo(
                [subcommand, *options, "--message-format=json-render-diagnostics", *flags],
                stdout=subprocess.PIPE, text=True, check=False,
                env={"CARGO_INCREMENTAL": "0"},
            )
            for line in result.stdout.splitlines():
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    print(line)
                    continue
                if not isinstance(event, dict):
                    continue
                if event.get("reason") == "compiler-message":
                    if rendered := event["message"].get("rendered"):
                        print(rendered, end="")
                    continue
                if event.get("reason") != "compiler-artifact":
                    continue
                executable = event.get("executable")
                if not executable:
                    continue
                if not (event["profile"]["test"] or discard_example and event["target"]["kind"] == ["example"]):
                    # One of the integration tests may run a normal package binary.
                    continue
                # The compiler can leave both its hashed output and a hard link
                # named after the target. Remove only the exact executable inode.
                installed = Path(executable).resolve()
                if not installed.is_relative_to(target_dir):
                    raise ValueError(f"Unexpected test executable outside the Bazel cache: {installed}")
                if installed.is_file():
                    for name in set(event["filenames"]) | {executable}:
                        path = Path(name).resolve()
                        if (path != installed and path.is_relative_to(target_dir)
                                and path.exists() and path.samefile(installed)):
                            path.unlink()
                    installed.unlink()
            successful = result.returncode == 0 and successful

        for start in range(0, len(integration), batch_size):
            tests = integration[start:start + batch_size]
            run("test", ["--no-fail-fast", *[arg for name in tests for arg in ("--test", name)]])
        run("test", ["--no-fail-fast", "--lib", *[arg for name in unit_binaries for arg in ("--bin", name)]])
        run("test", ["--no-fail-fast", "--doc"])
        for start in range(0, len(examples), batch_size):
            group = examples[start:start + batch_size]
            run("build", ["--profile", "test", *[arg for name in group for arg in ("--example", name)]],
                discard_example=True)
        if not successful:
            raise subprocess.CalledProcessError(1, f"{self.project} Rust test batches")

    def python(self, script, *args, env=None, cwd=None):
        return command([self.config["python"], script, *args], cwd or self.root, dict(self.env, **(env or {})))

    def run(self, args, *, env=None, cwd=None, **kwargs):
        return command(args, cwd or self.root, dict(self.env, **(env or {})), **kwargs)

    def profile(self, normal=()):
        return ["--release", "--all-features"] if self.release else list(normal)

    def cargo_output(self, profile=None):
        return self.target_dir / self.target / (profile or ("release" if self.release else "debug"))

    def write_output(self, name, source):
        destination = self.outputs[name]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    def copy_artifacts(self, name, src, files):
        directory = self.outputs[name]
        directory.mkdir(parents=True, exist_ok=True)
        for filename in files:
            shutil.copy2(src / filename, directory / filename)

    def tool(self, package, *args, project=None, binary=False, env=None):
        switches = ["--bin" if binary else "-p", package]
        return self.cargo(["run", *switches, "--", *args], project=project, env=env)

    def assert_runnable(self):
        if self.target != self.rust_host and not self.target_runner:
            raise ValueError(f"cross-running {self.target} requires a Cargo target runner in RSSN_TARGET_ENV")

    def assert_native(self):
        if self.target != self.rust_host:
            raise ValueError(f"this workflow executes native artifacts; remove rust_target={self.target}")

    def checkout(self, project):
        source = self.inputs[f"{project}/Cargo.toml"].resolve()
        if not source.is_file():
            raise ValueError(f"missing checkout for {project}")
        return source.parent

    @contextmanager
    def refresh(self, patterns):
        # Publish only outputs selected by the caller and detect changes made by
        # another editor while the generators were running.
        before = {}
        for pattern in patterns:
            before.update({p: p.read_bytes() for p in self.mirror.glob(pattern) if p.is_file()})
        try:
            yield
            changed = []
            for pattern in patterns:
                for artifact in self.mirror.glob(pattern):
                    if not artifact.is_file() or artifact in changed:
                        continue
                    data = artifact.read_bytes()
                    if before.get(artifact) == data:
                        continue
                    project, *rest = artifact.relative_to(self.mirror).parts
                    dest = self.checkout(project).joinpath(*rest)
                    current = dest.read_bytes() if dest.is_file() else None
                    if current not in (before.get(artifact), data):
                        raise ValueError(f"{dest} changed while generation was running; rerun after reviewing the change")
                    changed.append(artifact)
            for artifact in changed:
                project, *rest = artifact.relative_to(self.mirror).parts
                dest = self.checkout(project).joinpath(*rest)
                write_atomic(dest, artifact.read_bytes())
                print(f"Updated {dest}")
        finally:
            for pattern in patterns:
                for artifact in self.mirror.glob(pattern):
                    if artifact.is_file() and artifact not in before:
                        artifact.unlink()


def main(settings, lookup, args):
    try:
        build = Build(settings, lookup)
        if build.project == "aevia":
            path = build.extra.pop(0)
            import importlib.util

            module = importlib.util.spec_from_file_location("aevia_tasks", path)
            tasks = importlib.util.module_from_spec(module)
            module.loader.exec_module(tasks)
        else:
            import rssn_tasks as tasks

        if settings["workflow"] == "llvm-build":
            # Keep CMake's working tree separate from all Cargo snapshots.
            build.state = build.state / "llvm"
            build.mirror = build.state / "src"
            build.root = build.mirror / "rssn"
        with lock(build.state):
            build.stage()
            tasks.dispatch(build, settings["workflow"], args)
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode or 1)
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"Bazel scientific build: {error}", file=sys.stderr)
        sys.exit(1)


def llvm_link(rustc, arguments):
    source = os.path.normcase(os.path.abspath(os.environ["RSSN_BAZEL_LATE_SOURCE"]))
    inputs = [i for i, arg in enumerate(arguments) if os.path.normcase(os.path.abspath(arg)) == source]
    if len(inputs) != 1 or "--out-dir" not in arguments:
        raise ValueError("expected a single crate source and Cargo's --out-dir for standalone LLVM")
    output_index = arguments.index("--out-dir") + 1
    output = Path(arguments[output_index])
    llvm = json.loads(os.environ["RSSN_BAZEL_LATE_TOOLS"])
    llvm_env = dict(os.environ, **json.loads(os.environ["RSSN_BAZEL_LATE_LOADER"]))
    plugin = os.environ["RSSN_BAZEL_LATE_PLUGIN"]
    libraries = json.loads(os.environ["RSSN_BAZEL_LATE_LIBRARIES"])
    # Match the Unix targets' PIC defaults for both shared libraries and archives.
    relocation = ["-relocation-model=pic"] if any(name.endswith((".so", ".dylib", ".a")) for name in libraries) else []

    with tempfile.TemporaryDirectory(dir=output, prefix="llvm-") as directory:
        stage = Path(directory)
        arguments[output_index] = directory
        # Preserve Cargo's jobserver file descriptors on Unix. The temporary
        # output directory also scopes the object paths in rustc's link record.
        subprocess.run([rustc, *arguments, "-Zno-link", "-Csave-temps=yes", "--emit=llvm-bc"],
                       check=True, close_fds=os.name == "nt")
        records, modules = list(stage.glob("*.rlink")), list(stage.glob("*.rcgu.bc"))
        if len(records) != 1 or len(modules) != 1 or not modules[0].with_suffix(".o").is_file():
            raise ValueError("fat LTO must leave one rustc link record and one codegen unit with its object")
        obj = modules[0].with_suffix(".o")
        optimized = stage / "optimized.bc"
        subprocess.run([llvm["opt"], f"-load-pass-plugin={plugin}", "-passes=rssn-release,verify",
                        modules[0], "-o", optimized], check=True, env=llvm_env)
        subprocess.run([llvm["llc"], "-O3", "-filetype=obj", "-function-sections", "-data-sections", *relocation,
                        optimized, "-o", obj], check=True, env=llvm_env)

        # Output paths come from the .rlink. The linking options come from the
        # new rustc session, so keep the original --extern, -L and -C arguments.
        arguments[inputs[0]] = str(records[0])
        subprocess.run([rustc, *arguments, "-Zlink-only"], check=True, close_fds=os.name == "nt")
        # The link record keeps Cargo's crate name and extra-filename suffix.
        stem = records[0].stem
        files = [stage / name.replace("rssn", stem, 1) for name in libraries]
        files.append(stage / (stem + ".d"))
        missing = [file.name for file in files if not file.is_file()]
        if missing:
            raise ValueError("standalone LLVM linking did not produce " + ", ".join(missing))
        pdb = stage / (stem + ".pdb")
        if pdb.is_file():
            files.append(pdb)
        for file in files:
            file.replace(output / file.name)


if __name__ == "__main__":
    if len(sys.argv) > 3 and sys.argv[1] == "llvm-link":
        try:
            llvm_link(sys.argv[2], sys.argv[3:])
        except subprocess.CalledProcessError as error:
            sys.exit(error.returncode or 1)
        except (OSError, ValueError, KeyError) as error:
            sys.exit(f"RSSN standalone LLVM: {error}")
        sys.exit(0)
    if len(sys.argv) != 3 or sys.argv[1] != "action":
        raise SystemExit("usage: runner.py action SPEC")
    main(json.loads(Path(sys.argv[2]).read_text(encoding="utf-8")), Path, [])
