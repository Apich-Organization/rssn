//! The default build script is standalone and works in archives and worktrees.
use std::env;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;
use std::process::Output;

fn run(
    script: &Path,
    root: &Path,
    out: &Path,
    vars: &[(&str, &str)],
) -> Output {
    let mut command = Command::new(script);
    command
        .env("CARGO_MANIFEST_DIR", root)
        .env("OUT_DIR", out)
        .env("TARGET", "test-target")
        .env("OPT_LEVEL", "2")
        .env("DEBUG", "false")
        .env_remove("DEV")
        .env_remove("SOURCE_DATE_EPOCH")
        .env_remove("VERGEN_IDEMPOTENT");
    for name in [
        "VERGEN_BUILD_DATE",
        "VERGEN_BUILD_TIMESTAMP",
        "VERGEN_GIT_SHA",
    ] {
        command.env_remove(name);
    }
    for (name, value) in vars {
        command.env(name, value);
    }
    command.output().expect("run build script")
}
fn text(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).unwrap()
}
fn git(
    root: &Path,
    args: &[&str],
) -> bool {
    Command::new("git")
        .current_dir(root)
        .args([
            "-c",
            "core.hooksPath=/dev/null",
            "-c",
            "commit.gpgsign=false",
        ])
        .args(args)
        .env("GIT_AUTHOR_NAME", "RSSN test")
        .env("GIT_AUTHOR_EMAIL", "test@example.invalid")
        .env("GIT_COMMITTER_NAME", "RSSN test")
        .env("GIT_COMMITTER_EMAIL", "test@example.invalid")
        .output()
        .is_ok_and(|result| result.status.success())
}
#[test]
fn standalone_build_script_supports_overlay_archives_and_worktrees() {
    let temp: PathBuf = env::temp_dir().join(format!("rssn-build-script-{}", std::process::id()));
    let _ = fs::remove_dir_all(&temp);
    let archive = temp.join("archive");
    let out = temp.join("out");
    fs::create_dir_all(&archive).unwrap();
    fs::create_dir_all(&out).unwrap();
    let executable = temp.join(format!("build-script{}", env::consts::EXE_SUFFIX));
    let rustc = env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("build.rs");
    let compiled = Command::new(rustc)
        .arg("--edition=2024")
        .arg(source)
        .arg("-o")
        .arg(&executable)
        .output()
        .unwrap();
    assert!(
        compiled.status.success(),
        "{}",
        String::from_utf8_lossy(&compiled.stderr)
    );

    let vars = [
        ("SOURCE_DATE_EPOCH", "1671809360"),
        ("VERGEN_IDEMPOTENT", "1"),
    ];
    let missing = run(&executable, &archive, &out, &vars);
    assert!(!missing.status.success());
    assert!(
        String::from_utf8_lossy(&missing.stderr)
            .contains("required source file signatures.overlay.json")
    );
    assert!(!out.join("signatures.overlay.json").exists());
    let baseline =
        fs::read(Path::new(env!("CARGO_MANIFEST_DIR")).join("signatures.overlay.json")).unwrap();
    fs::write(archive.join("signatures.overlay.json"), &baseline).unwrap();
    let first = run(&executable, &archive, &out, &vars);
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert_eq!(
        fs::read(out.join("signatures.overlay.json")).unwrap(),
        baseline
    );
    let printed = text(&first);
    assert!(printed.contains("cargo:rustc-env=VERGEN_BUILD_DATE=2022-12-23\n"));
    assert!(
        printed.contains("cargo:rustc-env=VERGEN_BUILD_TIMESTAMP=2022-12-23T15:29:20.000000000Z\n")
    );
    assert!(printed.contains("cargo:rustc-env=VERGEN_GIT_SHA=VERGEN_IDEMPOTENT_OUTPUT\n"));
    assert!(printed.contains("cargo:rerun-if-changed=signatures.overlay.json\n"));
    assert!(printed.contains("cargo:rerun-if-env-changed=DEV\n"));
    assert_eq!(printed, text(&run(&executable, &archive, &out, &vars)));

    let overlay = "{ \"ops\": { \"test\": [] } }\n";
    fs::write(archive.join("signatures.overlay.json"), overlay).unwrap();
    let second = run(
        &executable,
        &archive,
        &out,
        &[
            ("SOURCE_DATE_EPOCH", "1671809360"),
            ("VERGEN_GIT_SHA", "archive-sha"),
        ],
    );
    assert!(second.status.success());
    assert_eq!(
        fs::read_to_string(out.join("signatures.overlay.json")).unwrap(),
        overlay
    );
    assert!(text(&second).contains("cargo:rustc-env=VERGEN_GIT_SHA=archive-sha\n"));
    // A stale OUT_DIR copy must not hide a subsequently missing source file.
    fs::remove_file(archive.join("signatures.overlay.json")).unwrap();
    let missing = run(&executable, &archive, &out, &vars);
    assert!(!missing.status.success());
    fs::write(archive.join("signatures.overlay.json"), &baseline).unwrap();
    let headers = run(&executable, &archive, &out, &[("DEV", "1")]);
    assert!(!headers.status.success());
    assert!(String::from_utf8_lossy(&headers.stderr).contains("--features dev-headers"));

    // Git is optional in source archives, but when present worktree-specific
    // HEAD and the common branch reference must both invalidate the build.
    let repo = temp.join("repository");
    fs::create_dir_all(&repo).unwrap();
    fs::write(repo.join("signatures.overlay.json"), &baseline).unwrap();
    if git(&repo, &["init", "-q"])
        && git(&repo, &["add", "signatures.overlay.json"])
        && git(&repo, &["commit", "-q", "--allow-empty", "-m", "initial"])
    {
        let worktree = temp.join("worktree");
        let path = worktree.to_str().unwrap();
        assert!(git(
            &repo,
            &["worktree", "add", "-q", "-b", "rssn-build-test", path]
        ));
        let before = run(&executable, &worktree, &out, &[]);
        assert!(before.status.success());
        let before = text(&before);
        assert!(before.contains("cargo:rerun-if-changed=") && before.contains("/HEAD"));
        assert!(before.contains("refs/heads/rssn-build-test"));
        assert!(git(
            &worktree,
            &["commit", "-q", "--allow-empty", "-m", "second"]
        ));
        let after = text(&run(&executable, &worktree, &out, &[]));
        let sha = |text: &str| {
            text.lines()
                .find(|line| line.starts_with("cargo:rustc-env=VERGEN_GIT_SHA="))
                .unwrap()
                .to_owned()
        };
        assert_ne!(sha(&before), sha(&after));
    }
    fs::remove_dir_all(temp).unwrap();
}
