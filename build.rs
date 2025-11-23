use vergen_gitcl::{
    BuildBuilder, CargoBuilder, Emitter, GitclBuilder, RustcBuilder, SysinfoBuilder,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut emitter = Emitter::default();

    emitter.add_instructions(&BuildBuilder::all_build()?)?;
    emitter.add_instructions(&CargoBuilder::all_cargo()?)?;
    emitter.add_instructions(&GitclBuilder::all_git()?)?;
    emitter.add_instructions(&RustcBuilder::all_rustc()?)?;
    emitter.add_instructions(&SysinfoBuilder::all_sysinfo()?)?;

    emitter.emit()?;

    // Generate C and C++ headers using cbindgen
    generate_headers()?;

    Ok(())
}

fn generate_headers() -> Result<(), Box<dyn std::error::Error>> {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR")?;
    
    // Generate C header using cbindgen.toml
    match cbindgen::generate(&crate_dir) {
        Ok(bindings) => {
            bindings.write_to_file("rssn.h");
            println!("cargo:warning=Generated rssn.h");
        }
        Err(e) => {
            println!("cargo:warning=Failed to generate C bindings: {:?}", e);
            println!("cargo:warning=Continuing build without C header generation");
        }
    }
    
    // Generate C++ header with custom config
    let cpp_config = cbindgen::Config {
        language: cbindgen::Language::Cxx,
        namespace: Some("rssn".to_string()),
        ..cbindgen::Config::from_file("cbindgen.toml")
            .unwrap_or_default()
    };
    
    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(cpp_config)
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file("rssn.hpp");
            println!("cargo:warning=Generated rssn.hpp");
        }
        Err(e) => {
            println!("cargo:warning=Failed to generate C++ bindings: {:?}", e);
            println!("cargo:warning=Continuing build without C++ header generation");
        }
    }
    
    println!("cargo:rerun-if-changed=src/");
    println!("cargo:rerun-if-changed=cbindgen.toml");
    
    Ok(())
}
