//! Bootstrap RSSN bindings from source syntax and reviewed semantic policies.
use rssn_api_model::{contract_codegen, numerical_codegen, structured_codegen};
use std::path::{Path, PathBuf};

fn generate(
    root: &Path,
    policy: &Path,
    structured: bool,
) -> Result<String, Box<dyn std::error::Error>> {
    let source = std::fs::read_to_string(policy)?;
    Ok(if structured {
        structured_codegen::generate(root, &source)?
    } else {
        numerical_codegen::generate(root, &source)?
    })
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let command = args
        .next()
        .ok_or("expected --bindings-{write,check,report} or --structured-{write,check,report}")?;
    let root = args.next().map_or_else(
        || PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../.."),
        PathBuf::from,
    );
    if args.next().is_some() {
        return Err("usage: rssn-api-gen COMMAND [RSSN_ROOT]".into());
    }
    if let Some(action) = command.strip_prefix("--contracts-") {
        let policy = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("contracts.tsv");
        let target = root.join("src/rules/numerical_abi/contract_generated.rs");
        let result = contract_codegen::generate(&std::fs::read_to_string(policy)?)?;
        match action {
            | "report" => print!("{result}"),
            | "write" => {
                std::fs::write(&target, result)?;
                println!("wrote {}", target.display());
            },
            | "check" => {
                if std::fs::read_to_string(&target)? != result {
                    return Err("generated contracts are stale".into());
                }
                println!("generated contracts verified");
            },
            | _ => return Err("unknown generator action".into()),
        }
        return Ok(());
    }
    let (structured, action) = if let Some(action) = command.strip_prefix("--structured-") {
        (true, action)
    } else if let Some(action) = command.strip_prefix("--bindings-") {
        (false, action)
    } else {
        return Err("unknown generator command".into());
    };
    let policy_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).to_path_buf();
    let policy = policy_dir.join(if structured {
        "structured-policy.tsv"
    } else {
        "numerical-policy.tsv"
    });
    let target = root.join(if structured {
        "src/rules/numerical_abi/structured_generated.rs"
    } else {
        "src/rules/numerical_abi/generated.rs"
    });
    let result = generate(&root, &policy, structured)?;
    match action {
        | "report" => print!("{result}"),
        | "write" => {
            std::fs::write(&target, result)?;
            println!("wrote {}", target.display());
        },
        | "check" => {
            if std::fs::read_to_string(&target)? != result {
                return Err("generated RSSN source is stale; review report and run write".into());
            }
            println!("generated RSSN source verified");
        },
        | _ => return Err("unknown generator action".into()),
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn generation_bootstraps_without_generated_rssn_files() {
        let root =
            std::env::temp_dir().join(format!("rssn-source-bootstrap-{}", std::process::id()));
        let kernels = root.join("src/kernels");
        std::fs::create_dir_all(&kernels).unwrap();
        std::fs::write(kernels.join("sample.rs"), "/// A documented scalar.\npub fn scale(x: f64) -> f64 { x }\npub struct Report { pub value: f64, pub ok: bool }\n").unwrap();
        let contract_dir = root.join("tools/rssn-api-gen");
        std::fs::create_dir_all(&contract_dir).unwrap();
        std::fs::write(
            contract_dir.join("contracts.tsv"),
            "Dimensionless\t-\t-\t-\t-\tDimensionless.\n",
        )
        .unwrap();
        let numerical = root.join("numerical.tsv");
        std::fs::write(
            &numerical,
            "crate::kernels::sample::scale\tnum_scale\tnumerical_special\tDimensionless\tfinite\n",
        )
        .unwrap();
        let structured = root.join("structured.tsv");
        std::fs::write(
            &structured,
            "crate::kernels::sample::Report\tREPORT\tbound\tD,E\tvalue=D;ok=-\n",
        )
        .unwrap();
        assert!(!root.join("src/rules/numerical_abi/generated.rs").exists());
        assert!(
            generate(&root, &numerical, false)
                .unwrap()
                .contains("num_scale")
        );
        assert!(
            generate(&root, &structured, true)
                .unwrap()
                .contains("record!(crate::kernels::sample::Report")
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
