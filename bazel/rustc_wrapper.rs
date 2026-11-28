//! Apply RSSN-specific compiler options to the Cargo units that request them.

use std::env;
use std::ffi::OsStr;
use std::process::Command;
use std::process::{
    self,
};

fn run() -> Result<i32, Box<dyn std::error::Error>> {
    let mut invocation = env::args_os().skip(1);
    let rustc = invocation.next().ok_or("Cargo did not supply rustc")?;
    let mut arguments: Vec<_> = invocation.collect();
    let late = arguments.iter().position(|arg| {
        arg.to_str()
            .is_some_and(|arg| arg.starts_with("--rssn-bazel-late="))
    });
    let compiling_rssn = arguments
        .windows(2)
        .any(|pair| pair[0] == "--crate-name" && pair[1] == "rssn")
        || arguments.iter().any(|arg| arg == "--crate-name=rssn");

    let mut command = if let Some(marker) = late {
        arguments.remove(marker);
        let mut command = Command::new(
            env::var_os("RSSN_BAZEL_LATE_PYTHON").ok_or("missing LLVM adapter Python")?,
        );
        command
            .arg(env::var_os("RSSN_BAZEL_LATE_ADAPTER").ok_or("missing LLVM adapter")?)
            .args([OsStr::new("llvm-link"), &rustc]);
        command
    } else if let Some(cache) = env::var_os("RSSN_BAZEL_SCCACHE") {
        let mut command = Command::new(cache);
        command.arg(&rustc);
        command
    } else {
        Command::new(&rustc)
    };
    command.args(arguments);

    let threads = env::var("RSSN_BAZEL_FRONTEND_THREADS")?;
    if compiling_rssn && threads != "1" {
        if env::var_os("RUSTC_BOOTSTRAP").as_deref() != Some(OsStr::new("1")) {
            command.env("RUSTC_BOOTSTRAP", "rssn");
        }
        command.arg(format!("-Zthreads={threads}"));
    }
    Ok(command.status()?.code().unwrap_or(1))
}

fn main() {
    match run() {
        | Ok(code) => process::exit(code),
        | Err(error) => {
            eprintln!("RSSN rustc wrapper: {error}");
            process::exit(1);
        },
    }
}
