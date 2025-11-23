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
    
    // Generate C header
    let c_config = cbindgen::Config {
        language: cbindgen::Language::C,
        cpp_compat: true,
        include_guard: Some("RSSN_H".to_string()),
        namespace: None,
        documentation: true,
        documentation_style: cbindgen::DocumentationStyle::C,
        export: cbindgen::ExportConfig {
            include: vec!["rssn".to_string()],
            exclude: vec![],
            ..Default::default()
        },
        parse: cbindgen::ParseConfig {
            parse_deps: false,
            include: None,
            exclude: vec![],
            clean: false,
            extra_bindings: vec![],
            expand: cbindgen::ParseExpandConfig {
                crates: vec![],
                all_features: false,
                default_features: true,
                features: vec![],
            },
        },
        ..Default::default()
    };
    
    match cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_config(c_config)
        .generate()
    {
        Ok(bindings) => {
            bindings.write_to_file("rssn.h");
            println!("cargo:warning=Generated rssn.h");
        }
        Err(e) => {
            println!("cargo:warning=Failed to generate C bindings: {:?}", e);
            println!("cargo:warning=Continuing build without C header generation");
        }
    }
    
    // Generate C++ header
    let cpp_config = cbindgen::Config {
        language: cbindgen::Language::Cxx,
        cpp_compat: true,
        include_guard: Some("RSSN_HPP".to_string()),
        namespace: Some("rssn".to_string()),
        documentation: true,
        documentation_style: cbindgen::DocumentationStyle::Doxy,
        export: cbindgen::ExportConfig {
            include: vec!["rssn".to_string()],
            exclude: vec![],
            ..Default::default()
        },
        parse: cbindgen::ParseConfig {
            parse_deps: false,
            include: None,
            exclude: vec![],
            clean: false,
            extra_bindings: vec![],
            expand: cbindgen::ParseExpandConfig {
                crates: vec![],
                all_features: false,
                default_features: true,
                features: vec![],
            },
        },
        ..Default::default()
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
    
    Ok(())
}
