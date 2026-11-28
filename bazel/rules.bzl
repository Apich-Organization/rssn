"""Bazel actions and local workflow launchers for the scientific workspaces."""

load("@scientific_host//:settings.bzl", "PLUGIN_FILENAME", "PYTHON", "WINDOWS")

LLVM_PLUGIN_OUTPUT = PLUGIN_FILENAME

def _source_map(ctx, runtime):
    sources = {}
    for target, project in [(ctx.attr.rssn, "rssn"), (ctx.attr.aevia, "aevia")]:
        if not target:
            continue
        for f in target.files.to_list():
            rel = f.short_path
            if rel.startswith("../"):
                rel = rel.split("/", 2)[2]
            sources[project + "/" + rel] = _runfile(ctx, f) if runtime else f.path
    return sources

def _runfile(ctx, f):
    if f.short_path.startswith("../"):
        return f.short_path[3:]
    return ctx.workspace_name + "/" + f.short_path

def _path(ctx, f, runtime):
    return _runfile(ctx, f) if runtime else f.path

def _settings(ctx, runtime):
    release_distribution = ctx.attr.workflow == "release-pass-test" or (
        ctx.attr.workflow == "ffi" and ctx.attr.arguments == ["both"]
    )
    return {
        "build": ctx.attr.arguments,
        "config": _path(ctx, ctx.file._config, runtime),
        "extra": [_path(ctx, f, runtime) for f in ctx.files.data],
        "plugin": _path(ctx, ctx.file.plugin, runtime) if ctx.file.plugin else "",
        "project": ctx.attr.project,
        "scope": ctx.attr.scope,
        "sources": _source_map(ctx, runtime),
        "target": ctx.var.get("rust_target", ""),
        "cpu_baseline": ctx.var.get("rssn_cpu_baseline", "portable"),
        "lto": ctx.var.get("rssn_lto", "fat" if release_distribution else "thin"),
        "codegen_units": ctx.var.get("rssn_codegen_units", "1" if release_distribution else "16"),
        "standalone_llvm": ctx.var.get("rssn_llvm_standalone") == "1" if release_distribution else False,
        "jobs": ctx.var.get("rssn_jobs", ""),
        "frontend_threads": ctx.var.get("rssn_frontend_threads", "1"),
        "wrapper": _path(ctx, ctx.file._wrapper, runtime),
        "workflow": ctx.attr.workflow,
        "timings": ctx.var.get("rssn_timings") == "1",
    }

def _cargo_artifact_impl(ctx):
    spec = ctx.actions.declare_file(ctx.label.name + ".action.json")
    settings = _settings(ctx, False)
    outputs = {}
    for name in ctx.attr.outputs:
        path = ctx.label.name + "/" + name
        if name in ctx.attr.directories:
            outputs[name] = ctx.actions.declare_directory(path)
        else:
            outputs[name] = ctx.actions.declare_file(path)
    if ctx.attr.workflow == "llvm-build":
        # Native CMake input is independent of the requested Rust target or LTO.
        settings.update(target = "", cpu_baseline = "portable", lto = "thin", timings = False, frontend_threads = "1", codegen_units = "16")
    if ctx.attr.workflow == "headers":
        settings["target"] = ""
    settings["outputs"] = {k: f.path for k, f in outputs.items()}
    ctx.actions.write(spec, json.encode(settings))
    srcs = ctx.files.rssn + ctx.files.aevia + ctx.files.data + [
        ctx.file._runner,
        ctx.file._tasks,
        ctx.file._config,
        ctx.file._wrapper,
        spec,
    ]
    if ctx.file.plugin:
        srcs.append(ctx.file.plugin)
    ctx.actions.run(
        executable = PYTHON,
        arguments = [ctx.file._runner.path, "action", spec.path],
        inputs = depset(srcs),
        outputs = outputs.values(),
        mnemonic = "Scientific" + ctx.attr.mnemonic,
        progress_message = ctx.attr.progress + " %{label}",
        # Cargo, CMake and the SDK use persistent local directories and host tools.
        execution_requirements = {"local": "1", "no-remote": "1"},
        use_default_shell_env = True,
    )
    return [DefaultInfo(files = depset(outputs.values()))]

def _shell_quote(value):
    return "'" + value.replace("'", "'\"'\"'") + "'"

def _workflow_impl(ctx):
    spec = ctx.actions.declare_file(ctx.label.name + ".workflow.json")
    settings = _settings(ctx, True)
    settings["tasks"] = _runfile(ctx, ctx.file._tasks)
    settings["runner"] = _runfile(ctx, ctx.file._runner)
    ctx.actions.write(spec, json.encode(settings))
    executable = ctx.actions.declare_file(ctx.label.name + (".bat" if WINDOWS else ".sh"))
    spec_key = _runfile(ctx, spec)

    # A minimal bootstrap locates the manifest or symlink tree for both run and test.
    bootstrap = _runfile(ctx, ctx.file._bootstrap)
    if WINDOWS:
        content = ('@echo off\r\n"%s" -c "%s" "%s" "%s" "%%~f0" %%*\r\n' %
                   (PYTHON, _BOOTSTRAP_COMMAND, bootstrap, spec_key))
    else:
        content = "#!/bin/sh\nexec %s -c %s %s %s \"$0\" \"$@\"\n" % (
            _shell_quote(PYTHON),
            _shell_quote(_BOOTSTRAP_COMMAND),
            _shell_quote(bootstrap),
            _shell_quote(spec_key),
        )
    ctx.actions.write(executable, content, is_executable = True)
    files = ctx.files.rssn + ctx.files.aevia + ctx.files.data + [
        spec,
        ctx.file._runner,
        ctx.file._tasks,
        ctx.file._config,
        ctx.file._bootstrap,
        ctx.file._wrapper,
    ]
    if ctx.file.plugin:
        files.append(ctx.file.plugin)
    return [DefaultInfo(executable = executable, runfiles = ctx.runfiles(files = files))]

_BOOTSTRAP_COMMAND = "; ".join([
    "import os,runpy,sys",
    "from pathlib import Path",
    "key, launcher = sys.argv[1], sys.argv[3]",
    "directory = os.getenv('RUNFILES_DIR') or os.getenv('TEST_SRCDIR') or launcher+'.runfiles'",
    "manifest = os.getenv('RUNFILES_MANIFEST_FILE') or launcher+'.runfiles_manifest'",
    "path = Path(directory, key)",
    "entry = None if path.is_file() else next(line for line in Path(manifest).read_text(encoding='utf-8').splitlines() if line.lstrip().startswith(key+' '))",
    "value = entry.lstrip().partition(' ')[2] if entry else str(path)",
    "value = value.replace(chr(92)+'s',' ').replace(chr(92)+'n',chr(10)).replace(chr(92)+'b',chr(92)) if entry and entry.startswith(' ') else value",
    "runpy.run_path(value, run_name='__main__')",
])

_attrs = {
    "aevia": attr.label(),
    "arguments": attr.string_list(),
    "data": attr.label_list(allow_files = True),
    "plugin": attr.label(allow_single_file = True, cfg = "exec"),
    "project": attr.string(default = "rssn"),
    "rssn": attr.label(),
    "scope": attr.string(),
    "workflow": attr.string(mandatory = True),
    "_bootstrap": attr.label(default = Label("//bazel:bootstrap.py"), allow_single_file = True),
    "_config": attr.label(default = Label("@scientific_host//:config.json"), allow_single_file = True),
    "_runner": attr.label(default = Label("//bazel:runner.py"), allow_single_file = True),
    "_tasks": attr.label(default = Label("//bazel:rssn_tasks.py"), allow_single_file = True),
    "_wrapper": attr.label(default = Label("//bazel:rustc_wrapper.rs"), allow_single_file = True),
}

cargo_artifact = rule(
    implementation = _cargo_artifact_impl,
    attrs = dict(_attrs, outputs = attr.string_list(mandatory = True), directories = attr.string_list(), mnemonic = attr.string(default = "Cargo"), progress = attr.string(default = "Building")),
)
scientific_run = rule(implementation = _workflow_impl, executable = True, attrs = _attrs)
scientific_test = rule(implementation = _workflow_impl, test = True, attrs = _attrs)
