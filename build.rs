use std::env;

use vergen_gitcl::Build;
use vergen_gitcl::Cargo;
use vergen_gitcl::Emitter;
use vergen_gitcl::Gitcl;
use vergen_gitcl::Rustc;
use vergen_gitcl::Sysinfo;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    const ENV_VAR_NAME: &str = "DEV";

    println!("cargo:rustc-env=VERGEN_SYSINFO_KERNEL_VERSION=unknown");

    println!("cargo:rustc-env=VERGEN_SYSINFO_KERNEL_VERSION=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_BRANCH=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_COMMIT_AUTHOR_EMAIL=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_COMMIT_AUTHOR_NAME=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_COMMIT_MESSAGE=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_DESCRIBE=unknown");
    println!("cargo:rustc-env=VERGEN_GIT_SHA=unknown");

    let mut emitter = Emitter::default();

    let build = Build::all_build();
    let gitcl = Gitcl::all_git();
    let si = Sysinfo::all_sysinfo();
    let rustc = Rustc::all_rustc();
    let cargo = Cargo::all_cargo();
    emitter
        .add_instructions(&gitcl)?
        .add_instructions(&build)?
        .add_instructions(&si)?
        .add_instructions(&rustc)?
        .add_instructions(&cargo)?
        .emit()?;

    // The signature overlay is a hand-maintained JSON file and `*.json` is
    // ignored by git, so a clean checkout does not contain it; `include_str!`
    // has no "optional" form and the crate would not compile. Materialise it
    // in `$OUT_DIR` first — unconditionally, because this has nothing to do
    // with the `DEV` gate below.
    write_overlay()?;

    println!("cargo:rerun-if-changed=signatures.overlay.json");

    // 1. Attempt to get the environment variable's value.
    match env::var(ENV_VAR_NAME) {
        // 2. If the variable is successfully retrieved (Ok(value)), check if the value is "1".
        | Ok(value) if value == "1" => {
            println!(
                "cargo:warning=Condition met: {} is set to '1'.",
                ENV_VAR_NAME
            );

            // 3. If the value is "1", call the function and use '?' to propagate errors.
            generate_headers()?;
        },
        // 4. If the variable is set to any other value, or not set (Err), skip.
        | _ => {
            println!(
                "cargo:warning=Skipping header generation. Set {}='1' to enable.",
                ENV_VAR_NAME
            );
        },
    }

    Ok(())
}

fn generate_headers() -> Result<(), Box<dyn std::error::Error>> {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR")?;

    // Generate C header using cbindgen.toml
    match cbindgen::generate(&crate_dir) {
        | Ok(bindings) => {
            bindings.write_to_file("rssn.h");

            println!("cargo:warning=Generated rssn.h");
        },
        | Err(e) => {
            println!(
                "cargo:warning=Failed \
                 to generate C \
                 bindings: {:?}",
                e
            );

            println!("cargo:warning=Continuing build without C header generation");
        },
    }

    // Generate C++ header with custom config
    let cpp_config = cbindgen::Config {
        language: cbindgen::Language::Cxx,
        namespace: Some("rssn".to_string()),
        ..cbindgen::Config::from_file("cbindgen.toml").unwrap_or_default()
    };

    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(cpp_config)
        .generate()
    {
        | Ok(bindings) => {
            bindings.write_to_file("rssn.hpp");

            println!("cargo:warning=Generated rssn.hpp");
        },
        | Err(e) => {
            println!(
                "cargo:warning=Failed \
                 to generate C++ \
                 bindings: {:?}",
                e
            );

            println!("cargo:warning=Continuing build without C++ header generation");
        },
    }

    println!("cargo:rerun-if-changed=src/");

    println!("cargo:rerun-if-changed=cbindgen.toml");

    println!("cargo:rerun-if-changed=Cargo.toml");

    Ok(())
}

/// Copies `signatures.overlay.json` into `$OUT_DIR`, or writes an empty one.
///
/// The overlay is the one place a human corrects an inferred signature, so it
/// has to exist for the crate to build at all: an empty overlay is a valid
/// overlay, a missing one is not.
fn write_overlay() -> Result<(), Box<dyn std::error::Error>> {
    const EMPTY: &str = "{\"ops\": {}}\n";
    let source = std::path::Path::new("signatures.overlay.json");
    let out_dir = env::var("OUT_DIR")?;
    let target = std::path::Path::new(&out_dir).join("signatures.overlay.json");
    if source.is_file() {
        std::fs::copy(source, &target)?;
        println!("cargo:warning=Using signatures.overlay.json from the crate root");
    } else {
        std::fs::write(&target, EMPTY)?;
        println!(
            "cargo:warning=signatures.overlay.json is missing; using an empty signature overlay"
        );
    }
    Ok(())
}
