//! Lean build metadata and optional C/C++ header generation.
use std::env;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::time::SystemTime;
use std::time::UNIX_EPOCH;
const KEYS: &[&str] = &[
    "VERGEN_BUILD_DATE",
    "VERGEN_BUILD_TIMESTAMP",
    "VERGEN_GIT_BRANCH",
    "VERGEN_GIT_SHA",
    "VERGEN_GIT_COMMIT_AUTHOR_NAME",
    "VERGEN_GIT_COMMIT_AUTHOR_EMAIL",
    "VERGEN_GIT_COMMIT_MESSAGE",
    "VERGEN_GIT_DESCRIBE",
    "VERGEN_RUSTC_SEMVER",
    "VERGEN_RUSTC_LLVM_VERSION",
    "VERGEN_RUSTC_HOST_TRIPLE",
    "VERGEN_CARGO_OPT_LEVEL",
    "VERGEN_CARGO_DEBUG",
    "VERGEN_CARGO_TARGET_TRIPLE",
    "VERGEN_SYSINFO_NAME",
    "VERGEN_SYSINFO_OS_VERSION",
    "VERGEN_SYSINFO_KERNEL_VERSION",
    "VERGEN_SYSINFO_CPU_BRAND",
    "VERGEN_SYSINFO_CPU_CORE_COUNT",
    "VERGEN_SYSINFO_TOTAL_MEMORY",
];
fn command(
    program: &str,
    args: &[&str],
) -> Option<String> {
    let output = Command::new(program).args(args).output().ok()?;
    output
        .status
        .success()
        .then(|| {
            String::from_utf8_lossy(&output.stdout)
                .trim()
                .replace(['\r', '\n'], " ")
        })
        .filter(|v| !v.is_empty())
}
fn git(args: &[&str]) -> Option<String> {
    command("git", args)
}
fn date(
    seconds: u64,
    nanos: u32,
) -> Option<(String, String)> {
    // Gregorian civil date, valid for the reproducible Unix epoch range.
    let days = i64::try_from(seconds / 86400).ok()?;
    let z = days.checked_add(719468)?;
    let era = z.div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    let y = y + i64::from(m <= 2);
    if !(0..=9999).contains(&y) {
        return None;
    }
    let day = format!("{y:04}-{m:02}-{d:02}");
    let tod = seconds % 86400;
    let timestamp = format!(
        "{day}T{:02}:{:02}:{:02}.{nanos:09}Z",
        tod / 3600,
        tod / 60 % 60,
        tod % 60
    );
    Some((day, timestamp))
}
fn os_release(key: &str) -> Option<String> {
    fs::read_to_string("/etc/os-release")
        .ok()?
        .lines()
        .find_map(|line| {
            line.strip_prefix(key)
                .and_then(|line| line.strip_prefix('='))
                .map(|v| v.trim_matches('"').to_owned())
        })
}
fn proc_field(
    path: &str,
    key: &str,
) -> Option<String> {
    fs::read_to_string(path).ok()?.lines().find_map(|line| {
        line.split_once(':')
            .filter(|(k, _)| k.trim() == key)
            .map(|(_, v)| v.trim().to_owned())
    })
}
fn memory_size(mut bytes: u64) -> String {
    let units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"];
    let mut index = 0;
    while bytes >= 1024 && index < units.len() - 1 {
        bytes /= 1024;
        index += 1;
    }
    format!("{bytes} {}", units[index])
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = env::var("CARGO_MANIFEST_DIR")?;
    env::set_current_dir(&root)?;
    let out = env::var("OUT_DIR")?;
    let overlay = Path::new("signatures.overlay.json");
    let target = Path::new(&out).join("signatures.overlay.json");
    fs::copy(overlay, target).map_err(|error| {
        std::io::Error::new(
            error.kind(),
            format!("required source file signatures.overlay.json: {error}; restore it from Git"),
        )
    })?;
    println!("cargo:rerun-if-changed=signatures.overlay.json");
    // The witnessed signatures are mined, not written (`signature-mine`), and
    // versioned so that the manifest of a clean clone does not depend on
    // having run the test suite. Unlike the overlay they are evidence, not
    // source: a tree without the file still builds, with the derived
    // signatures only.
    let witnessed = Path::new("signatures.witnessed.json");
    let target = Path::new(&out).join("signatures.witnessed.json");
    if witnessed.exists() {
        fs::copy(witnessed, target)?;
    } else {
        fs::write(target, "{ \"version\": 1, \"ops\": {} }\n")?;
    }
    println!("cargo:rerun-if-changed=signatures.witnessed.json");
    for name in ["SOURCE_DATE_EPOCH", "VERGEN_IDEMPOTENT", "DEV", "RUSTC"] {
        println!("cargo:rerun-if-env-changed={name}");
    }
    for key in KEYS {
        println!("cargo:rerun-if-env-changed={key}");
    }
    for arg in ["HEAD", "packed-refs"] {
        if let Some(path) = git(&["rev-parse", "--git-path", arg]) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    if let Some(head) = git(&["symbolic-ref", "-q", "HEAD"])
        && let Some(path) = git(&["rev-parse", "--git-path", &head])
    {
        println!("cargo:rerun-if-changed={path}");
    }
    let idempotent = env::var_os("VERGEN_IDEMPOTENT").is_some();
    let time = match env::var("SOURCE_DATE_EPOCH") {
        | Ok(value) => {
            Some(
                date(value.parse::<u64>()?, 0)
                    .ok_or("SOURCE_DATE_EPOCH is outside the supported calendar range")?,
            )
        },
        | Err(env::VarError::NotPresent) if !idempotent => {
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .ok()
                .and_then(|value| date(value.as_secs(), value.subsec_nanos()))
        },
        | Err(env::VarError::NotPresent) => None,
        | Err(error) => return Err(error.into()),
    };
    let rustc = env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    let rustc_info = Command::new(rustc)
        .arg("-vV")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).into_owned())
        .unwrap_or_default();
    let compiler = |key: &str| {
        rustc_info
            .lines()
            .find_map(|line| line.strip_prefix(key).map(str::to_owned))
    };
    let placeholder = if idempotent {
        "VERGEN_IDEMPOTENT_OUTPUT"
    } else {
        "unknown"
    };
    for &key in KEYS {
        let value = env::var(key)
            .ok()
            .or_else(|| {
                match key {
                    | "VERGEN_BUILD_DATE" => time.as_ref().map(|t| t.0.clone()),
                    | "VERGEN_BUILD_TIMESTAMP" => time.as_ref().map(|t| t.1.clone()),
                    | "VERGEN_CARGO_OPT_LEVEL" => env::var("OPT_LEVEL").ok(),
                    | "VERGEN_CARGO_DEBUG" => env::var("DEBUG").ok(),
                    | "VERGEN_CARGO_TARGET_TRIPLE" => env::var("TARGET").ok(),
                    | "VERGEN_RUSTC_SEMVER" => compiler("release: "),
                    | "VERGEN_RUSTC_LLVM_VERSION" => compiler("LLVM version: "),
                    | "VERGEN_RUSTC_HOST_TRIPLE" => compiler("host: "),
                    | _ if idempotent => None,
                    | "VERGEN_GIT_BRANCH" => git(&["rev-parse", "--abbrev-ref", "HEAD"]),
                    | "VERGEN_GIT_SHA" => git(&["rev-parse", "HEAD"]),
                    | "VERGEN_GIT_COMMIT_AUTHOR_NAME" => git(&["log", "-1", "--format=%an"]),
                    | "VERGEN_GIT_COMMIT_AUTHOR_EMAIL" => git(&["log", "-1", "--format=%ae"]),
                    | "VERGEN_GIT_COMMIT_MESSAGE" => git(&["log", "-1", "--format=%s"]),
                    | "VERGEN_GIT_DESCRIBE" => git(&["describe", "--always"]),
                    | "VERGEN_SYSINFO_NAME" => os_release("NAME"),
                    | "VERGEN_SYSINFO_OS_VERSION" => {
                        command("uname", &["-s"])
                            .zip(os_release("NAME"))
                            .map(|(kernel, release)| format!("{kernel} ({release})"))
                    },
                    | "VERGEN_SYSINFO_KERNEL_VERSION" => command("uname", &["-r"]),
                    | "VERGEN_SYSINFO_CPU_BRAND" => proc_field("/proc/cpuinfo", "model name"),
                    | "VERGEN_SYSINFO_CPU_CORE_COUNT" => {
                        std::thread::available_parallelism()
                            .ok()
                            .map(|n| n.to_string())
                    },
                    | "VERGEN_SYSINFO_TOTAL_MEMORY" => {
                        proc_field("/proc/meminfo", "MemTotal")
                            .and_then(|value| value.split_whitespace().next()?.parse::<u64>().ok())
                            .and_then(|kb| kb.checked_mul(1024))
                            .map(memory_size)
                    },
                    | _ => None,
                }
            })
            .unwrap_or_else(|| placeholder.into());
        println!("cargo:rustc-env={key}={}", value.replace(['\r', '\n'], " "));
    }
    if env::var("DEV").as_deref() == Ok("1") {
        headers()?;
    }
    Ok(())
}
#[cfg(not(feature = "dev-headers"))]
fn headers() -> Result<(), Box<dyn std::error::Error>> {
    Err("DEV=1 requires --features dev-headers".into())
}
#[cfg(feature = "dev-headers")]
fn headers() -> Result<(), Box<dyn std::error::Error>> {
    let root = env::var("CARGO_MANIFEST_DIR")?;
    cbindgen::generate(&root)?.write_to_file(Path::new(&root).join("rssn.h"));
    let config = cbindgen::Config {
        language: cbindgen::Language::Cxx,
        namespace: Some("rssn".into()),
        ..cbindgen::Config::from_file("cbindgen.toml")?
    };
    cbindgen::Builder::new()
        .with_crate(&root)
        .with_config(config)
        .generate()?
        .write_to_file(Path::new(&root).join("rssn.hpp"));
    for path in ["src", "cbindgen.toml", "Cargo.toml"] {
        println!("cargo:rerun-if-changed={path}");
    }
    Ok(())
}
