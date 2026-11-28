"""Local Cargo/Rust/LLVM tool discovery, shared by both Bazel modules."""

def _rust_tool(ctx, name, override, windows):
    if override:
        return ctx.path(override)
    if windows:
        # An SSH session on Windows may be unable to execute the Rustup symlink.
        rustup = ctx.which("rustup")
        if rustup:
            result = ctx.execute([str(rustup), "which", name], working_directory = str(ctx.workspace_root))
            if result.return_code == 0:
                executable = ctx.path(result.stdout.strip())
                if executable.exists:
                    return executable
    return ctx.which(name)

# Official LLVM release SDKs for hosted Unix release jobs. The selected
# compiler's LLVM version, operating system and architecture determine the asset.
_SDK_ASSETS = {
    ("23.1.5", "linux", "x86_64"): ("Linux-X64", "c3107384835f788768f614db610e87fba7077f125ca950eb02a78bab0a308065"),
    ("23.1.5", "linux", "aarch64"): ("Linux-ARM64", "830ceb70edf63849ce66e97bf85310a3bf2bdd04c860a9f2fa4ef1560e290cef"),
    ("23.1.5", "mac os x", "aarch64"): ("macOS-ARM64", "503de30ef90311e0b143623801b19af2250adac6d69faed65ad0051170202e3b"),
}

def _release_sdk(ctx, llvm_version):
    arch = ctx.os.arch
    if arch == "amd64":
        arch = "x86_64"
    if arch == "arm64":
        arch = "aarch64"
    asset = _SDK_ASSETS.get((llvm_version, ctx.os.name, arch))
    if not asset:
        return None
    platform, sha256 = asset
    filename = "LLVM-%s-%s.tar.zst" % (llvm_version, platform)
    tar = ctx.which("tar")
    if not tar or not ctx.which("zstd"):
        fail("Install tar and zstd to fetch the LLVM SDK, or set RSSN_LLVM_CONFIG")
    archive = ctx.path(filename)
    ctx.download(
        url = "https://github.com/llvm/llvm-project/releases/download/llvmorg-%s/%s" % (llvm_version, filename),
        sha256 = sha256,
        output = archive,
    )
    root = ctx.path("sdk")
    ctx.file("sdk/.bazel-sdk")
    extract = ctx.execute([str(tar), "--use-compress-program=zstd", "-xf", str(archive), "-C", str(root)], timeout = 900)
    if extract.return_code:
        fail("Could not extract %s: %s" % (filename, extract.stderr))
    ctx.delete(archive)
    for directory in [root] + root.readdir():
        config = directory.get_child("bin").get_child("llvm-config")
        if config.exists:
            return config
    fail("The pinned LLVM SDK %s does not contain bin/llvm-config" % filename)

def _host_impl(ctx):
    env = {name: ctx.getenv(name) for name in [
        "ANDROID_NDK_HOME",
        "BAZEL_SCI_CACHE",
        "CARGO",
        "CARGO_HOME",
        "CMAKE_GENERATOR",
        "CARGO_TARGET_AARCH64_APPLE_DARWIN_LINKER",
        "CARGO_TARGET_AARCH64_PC_WINDOWS_MSVC_LINKER",
        "CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER",
        "CARGO_TARGET_X86_64_APPLE_DARWIN_LINKER",
        "CARGO_TARGET_X86_64_PC_WINDOWS_MSVC_LINKER",
        "CC",
        "CC_aarch64_pc_windows_msvc",
        "CC_aarch64_unknown_linux_gnu",
        "CC_x86_64_pc_windows_msvc",
        "CFLAGS",
        "CFLAGS_aarch64_pc_windows_msvc",
        "CFLAGS_aarch64_unknown_linux_gnu",
        "CFLAGS_x86_64_pc_windows_msvc",
        "CXX",
        "HOME",
        "INCLUDE",
        "LIB",
        "LIBPATH",
        "LOCALAPPDATA",
        "MACOSX_DEPLOYMENT_TARGET",
        "PATH",
        "PYTHON",
        "RSSN_FETCH_LLVM_SDK",
        "RSSN_LLVM_CONFIG",
        "RSSN_LLVM_IMPORT_LIBRARY",
        "RSSN_TARGET_ENV",
        "RSSN_SCCACHE",
        "RUSTC",
        "SDKROOT",
        "SYSTEMROOT",
        "USERPROFILE",
    ]}
    windows = "windows" in ctx.os.name
    python = None
    candidates = [ctx.path(env["PYTHON"])] if env["PYTHON"] else [ctx.which("python3"), ctx.which("python")]
    for candidate in candidates:
        if candidate:
            version = ctx.execute([str(candidate), "-c", "import sys; sys.exit(sys.version_info < (3, 11))"])
            if version.return_code == 0:
                python = candidate
                break
    if not python:
        fail("Install Python 3.11+ or pass --repo_env=PYTHON=/path/to/python")
    rustc = _rust_tool(ctx, "rustc", env["RUSTC"], windows)
    cargo = _rust_tool(ctx, "cargo", env["CARGO"], windows)
    rust = {}
    if rustc and cargo:
        version = ctx.execute([str(rustc), "-vV"], working_directory = str(ctx.workspace_root))
        sysroot = ctx.execute([str(rustc), "--print", "sysroot"], working_directory = str(ctx.workspace_root))
        cargo_version = ctx.execute([str(cargo), "-V"], working_directory = str(ctx.workspace_root))
        if version.return_code == 0 and sysroot.return_code == 0 and cargo_version.return_code == 0:
            rust = {
                "version": version.stdout,
                "cargo_version": cargo_version.stdout,
                # Freeze rustup's resolution so an action in the mirror uses the same toolchain.
                "bin": str(ctx.path(sysroot.stdout.strip()).get_child("bin")),
                "rustc": str(rustc) if env["RUSTC"] else "",
                "cargo": str(cargo) if env["CARGO"] else "",
            }
    llvm_version = ""
    for line in rust.get("version", "").splitlines():
        if line.startswith("LLVM version: "):
            llvm_version = line[len("LLVM version: "):]
    llvm_major = llvm_version.split(".")[0]
    sdk = env["RSSN_LLVM_CONFIG"]
    if not sdk and env["RSSN_FETCH_LLVM_SDK"] == "1":
        sdk = _release_sdk(ctx, llvm_version)
    sdk = sdk or (ctx.which("llvm-config-" + llvm_major) if llvm_major else None) or ctx.which("llvm-config")
    sdk_info = {}
    if sdk:
        for option in ["version", "includedir", "bindir", "libdir", "cxxflags"]:
            result = ctx.execute([str(sdk), "--" + option])
            if result.return_code == 0:
                sdk_info[option] = result.stdout.strip()

    # CMake uses the default module suffix (.so on ELF and Darwin, .dll on Windows).
    plugin_filename = "RssnLLVMOptimizer.dll" if windows else "RssnLLVMOptimizer.so"
    target_env = {}
    if env["RSSN_TARGET_ENV"]:
        target_env = json.decode(ctx.read(env["RSSN_TARGET_ENV"]))
    config = {
        "workspace": str(ctx.workspace_root),
        "rust": rust,
        "cargo": str(cargo) if cargo else "",
        "python": str(python),
        "windows": windows,
        "sdk": str(sdk) if sdk else "",
        "sdk_info": sdk_info,
        "llvm_major": llvm_major,
        "llvm_version": llvm_version,
        "target_env": target_env,
        "env": {k: v for k, v in env.items() if v != None},
    }
    ctx.file("config.json", json.encode_indent(config) + "\n")
    ctx.file("settings.bzl", "PYTHON = %r\nWINDOWS = %r\nPLUGIN_FILENAME = %r\n" % (str(python), windows, plugin_filename))
    ctx.file("BUILD.bazel", 'exports_files(["config.json", "settings.bzl"], visibility = ["//visibility:public"])\n')

_scientific_host = repository_rule(
    implementation = _host_impl,
    local = True,
    configure = True,
)

def _host_extension(_):
    _scientific_host(name = "scientific_host")

host = module_extension(implementation = _host_extension)
