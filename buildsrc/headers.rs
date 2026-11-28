//! Shared C/C++ header generation for Cargo and the standalone generator.

use std::path::Path;

pub fn generate(
    root: &Path
) -> Result<[(&'static str, cbindgen::Bindings); 2], Box<dyn std::error::Error>> {
    let config = cbindgen::Config::from_file(root.join("cbindgen.toml"))?;
    let cpp_config = cbindgen::Config {
        language: cbindgen::Language::Cxx,
        ..config.clone()
    };
    let generate = |config| {
        cbindgen::Builder::new()
            .with_src(root.join("src/ffi/mod.rs"))
            // Preserve the existing opaque Ball and Wavelet declarations.
            .with_src(root.join("src/graph/store.rs"))
            .with_src(root.join("src/kernels/wavelets.rs"))
            .with_config(config)
            .generate()
    };
    Ok([
        ("rssn.h", generate(config)?),
        ("rssn.hpp", generate(cpp_config)?),
    ])
}

/// Write the canonical header, formatting the C++ declarations before publication.
pub fn header_bytes(
    name: &str,
    bindings: &cbindgen::Bindings,
) -> std::io::Result<Vec<u8>> {
    use std::io::Write;
    use std::process::Command;
    use std::process::Stdio;

    let mut generated = Vec::new();
    bindings.write(&mut generated);
    if !name.ends_with(".hpp") {
        return Ok(generated);
    }
    let formatter = std::env::var_os("CLANG_FORMAT").unwrap_or_else(|| "clang-format".into());
    let mut child = Command::new(formatter)
        .args(["--style=LLVM", &format!("--assume-filename={name}")])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    let mut input = child.stdin.take().expect("clang-format stdin is piped");
    // Read the formatted header while supplying input to avoid filling either pipe.
    let writer = std::thread::spawn(move || input.write_all(&generated));
    let output = child.wait_with_output();
    writer
        .join()
        .map_err(|_| std::io::Error::other("clang-format input thread failed"))??;
    let output = output?;
    if !output.status.success() {
        return Err(std::io::Error::other(format!(
            "clang-format failed: {}",
            output.status
        )));
    }
    Ok(output.stdout)
}
