use std::env;
use std::io;

use vergen_gitcl::BuildBuilder;
use vergen_gitcl::CargoBuilder;
use vergen_gitcl::Emitter;
use vergen_gitcl::GitclBuilder;
use vergen_gitcl::RustcBuilder;
use vergen_gitcl::SysinfoBuilder;

fn main() -> Result<
    (),
    Box<dyn std::error::Error>,
> {

    const ENV_VAR_NAME: &str = "DEV";

    let mut emitter =
        Emitter::default();

    emitter.add_instructions(
        &BuildBuilder::all_build()?,
    )?;

    emitter.add_instructions(
        &CargoBuilder::all_cargo()?,
    )?;

    emitter.add_instructions(
        &GitclBuilder::all_git()?,
    )?;

    emitter.add_instructions(
        &RustcBuilder::all_rustc()?,
    )?;

    emitter.add_instructions(
        &SysinfoBuilder::all_sysinfo()?,
    )?;

    emitter.emit()?;

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

fn generate_headers() -> Result<
    (),
    Box<dyn std::error::Error>,
> {

    let crate_dir = std::env::var(
        "CARGO_MANIFEST_DIR",
    )?;

    // Generate C header using cbindgen.toml
    match cbindgen::generate(&crate_dir)
    {
        | Ok(bindings) => {

            bindings.write_to_file(
                "rssn.h",
            );

            println!(
                "cargo:warning=Generated rssn.h"
            );
        },
        | Err(e) => {

            println!(
                "cargo:warning=Failed \
                 to generate C \
                 bindings: {:?}",
                e
            );

            println!(
                "cargo:warning=Continuing build without C header generation"
            );
        },
    }

    // Generate C++ header with custom config
    let cpp_config = cbindgen::Config {
        language:
            cbindgen::Language::Cxx,
        namespace: Some(
            "rssn".to_string(),
        ),
        ..cbindgen::Config::from_file(
            "cbindgen.toml",
        )
        .unwrap_or_default()
    };

    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(cpp_config)
        .generate()
    {
        | Ok(bindings) => {

            bindings.write_to_file(
                "rssn.hpp",
            );

            println!(
                "cargo:warning=Generated rssn.hpp"
            );
        },
        | Err(e) => {

            println!(
                "cargo:warning=Failed \
                 to generate C++ \
                 bindings: {:?}",
                e
            );

            println!(
                "cargo:warning=Continuing build without C++ header generation"
            );
        },
    }

    println!(
        "cargo:rerun-if-changed=src/"
    );

    println!(
        "cargo:rerun-if-changed=cbindgen.toml"
    );

    Ok(())
}
