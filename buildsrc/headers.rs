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
