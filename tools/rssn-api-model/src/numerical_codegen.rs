//! Generate declared numerical adapters from Rust syntax plus reviewed policy.
//! The policy contains semantic choices; it never describes Rust types.
use quote::ToTokens;
use std::collections::BTreeSet;
use std::path::Path;
use syn::FnArg;
use syn::Item;
use syn::ReturnType;
use syn::Visibility;

fn compact(tokens: impl ToTokens) -> String {
    tokens
        .to_token_stream()
        .to_string()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}
fn function_shape(
    f: &syn::ItemFn,
    output: &str,
) -> Result<&'static str, String> {
    let sig = &f.sig;
    if f.attrs
        .iter()
        .any(|attr| attr.path().is_ident("cfg") || attr.path().is_ident("cfg_attr"))
        || !matches!(f.vis, Visibility::Public(_))
        || !sig.generics.params.is_empty()
        || sig.generics.where_clause.is_some()
        || sig.asyncness.is_some()
        || sig.unsafety.is_some()
        || sig.abi.is_some()
        || sig.variadic.is_some()
    {
        return Err(format!(
            "{} is not an eligible public monomorphic Rust function",
            sig.ident
        ));
    }
    let args: Vec<_> = sig
        .inputs
        .iter()
        .map(|arg| {
            match arg {
                | FnArg::Typed(arg) => Ok(compact(&arg.ty)),
                | _ => Err(format!("{} has a receiver", sig.ident)),
            }
        })
        .collect::<Result<_, _>>()?;
    let result = match &sig.output {
        | ReturnType::Type(_, ty) => compact(ty),
        | ReturnType::Default => "()".into(),
    };
    match (
        args.iter()
            .map(String::as_str)
            .collect::<Vec<_>>()
            .as_slice(),
        result.as_str(),
        output,
    ) {
        | (["f64"], "f64", "finite") => Ok("Scalar"),
        | (["f64"], "f64", "allow_infinity") => Ok("ScalarAllowInfinity"),
        | (["&[f64]"], "f64", "finite") => Ok("VectorScalar"),
        | (["&[f64]"], "Result<f64,String>", "finite") => Ok("VectorScalarFallible"),
        | (["&[f64]"], "Vec<f64>", "finite") => Ok("VectorVector"),
        | (["&mut[f64]"], "f64", "finite") => Ok("MutableVectorScalar"),
        | (["&[f64]", "f64"], "f64", "finite") => Ok("VectorScalarScalar"),
        | (["&[f64]", "f64"], "Vec<f64>", "finite") => Ok("VectorScalarVector"),
        | (["&[f64]"], "Result<Vec<f64>,String>", "finite") => Ok("VectorVectorFallible"),
        | (["&[f64]", "&[f64]"], "Result<f64,String>", "finite") => Ok("VectorsScalar"),
        | (["&[f64]", "&[f64]"], "Result<Vec<f64>,String>", "finite") => Ok("VectorsVector"),
        | _ => {
            Err(format!(
                "{} has unsupported signature or output policy: ({}) -> {result}; {output}",
                sig.ident,
                args.join(", ")
            ))
        },
    }
}
fn docs(f: &syn::ItemFn) -> Result<String, String> {
    let mut lines = Vec::new();
    for attr in &f.attrs {
        if !attr.path().is_ident("doc") {
            continue;
        }
        if let syn::Meta::NameValue(value) = &attr.meta
            && let syn::Expr::Lit(syn::ExprLit {
                lit: syn::Lit::Str(text),
                ..
            }) = &value.value
        {
            let text = text.value();
            let line = text.trim();
            if line.is_empty() {
                if !lines.is_empty() {
                    break;
                }
            } else {
                lines.push(line.to_owned());
            }
        }
    }
    if lines.is_empty() {
        Err(format!("{} has no Rust documentation", f.sig.ident))
    } else {
        Ok(lines.join(" "))
    }
}
fn ident(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .enumerate()
            .all(|(i, c)| c == '_' || c.is_ascii_alphabetic() || (i > 0 && c.is_ascii_digit()))
}
/// Generate a stable checked-in Rust fragment. Each policy path must resolve to
/// exactly one public free function with an explicitly supported concrete ABI.
pub fn generate(
    root: &Path,
    policy: &str,
) -> Result<String, String> {
    let contract_path = root.join("tools/rssn-api-gen/contracts.tsv");
    let contract_source = std::fs::read_to_string(&contract_path)
        .map_err(|error| format!("{}: {error}", contract_path.display()))?;
    let contracts = crate::contract_codegen::simple_summaries(&contract_source)?;
    let mut rows = Vec::new();
    let mut names = BTreeSet::new();
    let mut paths = BTreeSet::new();
    for (line, text) in policy.lines().enumerate() {
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = text.split('\t').collect();
        let (path, name, namespace, dimension, output, abi_doc) = match fields.as_slice() {
            | [path, name, namespace, dimension, output] => {
                (*path, *name, *namespace, *dimension, *output, None)
            },
            | [path, name, namespace, dimension, output, doc] if !doc.is_empty() => {
                (*path, *name, *namespace, *dimension, *output, Some(*doc))
            },
            | _ => {
                return Err(format!(
                    "policy line {} must have five fields plus an optional legacy ABI description",
                    line + 1
                ));
            },
        };
        let parts: Vec<_> = path.split("::").collect();
        if parts.len() != 4
            || parts[0] != "crate"
            || parts[1] != "kernels"
            || !parts[2..].iter().all(|p| ident(p))
            || !ident(name)
            || !ident(namespace)
        {
            return Err(format!("invalid policy identifiers on line {}", line + 1));
        }
        let contract = contracts
            .get(dimension)
            .ok_or_else(|| format!("unknown dimension policy {dimension} on line {}", line + 1))?;
        if contract.structured_result || contract.nested_paths {
            return Err(format!(
                "dimension policy {dimension} needs a structured adapter"
            ));
        }
        if !names.insert(name) || !paths.insert(path) {
            return Err(format!("duplicate export or source on line {}", line + 1));
        }
        let source = root.join("src/kernels").join(format!("{}.rs", parts[2]));
        let text =
            std::fs::read_to_string(&source).map_err(|e| format!("{}: {e}", source.display()))?;
        let parsed = syn::parse_file(&text).map_err(|e| format!("{}: {e}", source.display()))?;
        let mut matches = parsed.items.iter().filter_map(|item| {
            match item {
                | Item::Fn(f) if f.sig.ident == parts[3] => Some(f),
                | _ => None,
            }
        });
        let function = matches
            .next()
            .ok_or_else(|| format!("missing public Rust function {path}"))?;
        if matches.next().is_some() {
            return Err(format!("ambiguous Rust function {path}"));
        }
        let adapter = function_shape(function, output)?;
        if contract
            .arguments
            .iter()
            .any(|index| *index >= function.sig.inputs.len())
        {
            return Err(format!(
                "dimension policy {dimension} refers to a missing argument of {path}"
            ));
        }
        let doc = match abi_doc {
            | Some(doc) => doc.to_string(),
            | None => docs(function)?,
        };
        rows.push((path, name, namespace, dimension, output, adapter, doc));
    }
    rows.sort_by_key(|row| row.0);
    let mut out = String::from(
        "// Generated by rssn-api-gen; edit numerical-policy.tsv or Rust source docs.\n\
        // Signatures come from Rust syntax and are checked by typed adapter pointers.\n\
        macro_rules! define_numerical_exports {\n\
        ($name:ident; $($manual:expr,)*) => {\n\
        /// Numerical exports with reviewed dimensional policies.\n\
        pub static $name: &[Export] = &[$($manual,)*\n",
    );
    for (path, name, namespace, dimension, output, adapter, doc) in rows {
        out.push_str(&format!("    Export {{ rust_path: {path:?}, name: {name:?}, namespace: {namespace:?}, category: {namespace:?}, doc: {doc:?}, stability: Stability::Experimental, allows_infinite_output: {}, dimension_policy: DimensionPolicy::{dimension}, adapter: Adapter::{adapter}({path}) }},\n", output == "allow_infinity"));
    }
    out.push_str("];\n};\n}\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn f(source: &str) -> syn::ItemFn {
        syn::parse_str(source).unwrap()
    }
    #[test]
    fn signature_detection_rejects_generic_unsupported_or_unaudited_outputs() {
        assert_eq!(
            function_shape(&f("pub fn f(x: &[f64]) -> f64 { 0.0 }"), "finite").unwrap(),
            "VectorScalar"
        );
        assert!(function_shape(&f("pub fn f<T>(x: T) -> f64 { 0.0 }"), "finite").is_err());
        assert_eq!(
            function_shape(&f("pub fn f(x: &mut [f64]) -> f64 { 0.0 }"), "finite").unwrap(),
            "MutableVectorScalar"
        );
        assert!(function_shape(&f("pub fn f(x: &mut [u8]) -> f64 { 0.0 }"), "finite").is_err());
        assert!(
            function_shape(&f("pub fn f(x: &[f64]) -> f64 { 0.0 }"), "allow_infinity").is_err()
        );
    }
    #[test]
    fn generated_output_is_deterministic_and_stale_policies_fail() {
        let root = std::env::temp_dir().join(format!("rssn-bindings-audit-{}", std::process::id()));
        let dir = root.join("src/kernels");
        std::fs::create_dir_all(&dir).unwrap();
        let policies = root.join("tools/rssn-api-gen");
        std::fs::create_dir_all(&policies).unwrap();
        std::fs::write(policies.join("contracts.tsv"), "Dimensionless\t-\t-\t-\t-\tDimensionless.\nPreserve\tD=0\t-\tD\t-\tPreserve.\nNovel\tD=0\t-\tD^2\t-\tNovel.\nRecord\tD=0\t-\t-\tx=D\tRecord.\nMissing\tD=1\t-\tD\t-\tMissing.\n").unwrap();
        std::fs::write(dir.join("fixture.rs"), "/// Square.\npub fn square(x:f64)->f64{x*x}\n/// Magnitude.\npub fn magnitude(x:&[f64])->f64{0.0}\n").unwrap();
        let a = "crate::kernels::fixture::square\tnum_square\tnumerical_fixture\tDimensionless\tfinite\n";
        let b = "crate::kernels::fixture::magnitude\tnum_magnitude\tnumerical_fixture\tPreserve\tfinite\n";
        assert_eq!(
            generate(&root, &format!("{a}{b}")).unwrap(),
            generate(&root, &format!("{b}{a}")).unwrap()
        );
        let legacy = "crate::kernels::fixture::square\tnum_square\tnumerical_fixture\tDimensionless\tfinite\tStable ABI documentation.\n";
        assert!(
            generate(&root, legacy)
                .unwrap()
                .contains("Stable ABI documentation.")
        );
        assert!(
            generate(
                &root,
                "crate::kernels::fixture::square\tnum_square\tnumerical_fixture\tNovel\tfinite\n"
            )
            .unwrap()
            .contains("DimensionPolicy::Novel")
        );
        assert!(
            generate(
                &root,
                "crate::kernels::fixture::square\tnum_square\tnumerical_fixture\tRecord\tfinite\n"
            )
            .is_err()
        );
        assert!(
            generate(
                &root,
                "crate::kernels::fixture::square\tnum_square\tnumerical_fixture\tMissing\tfinite\n"
            )
            .is_err()
        );
        assert!(generate(&root, "crate::kernels::fixture::missing\tnum_missing\tnumerical_fixture\tPreserve\tfinite\n").is_err());
        std::fs::write(
            dir.join("fixture.rs"),
            "/// Changed.\npub fn square(x:&mut [u8])->f64{0.0}\n",
        )
        .unwrap();
        assert!(generate(&root, a).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }
    #[test]
    fn rustdoc_not_manifest_supplies_description() {
        assert_eq!(
            docs(&f(
                "/// First line.\n/// Second line.\n///\n/// More details.\npub fn f(x:f64)->f64{x}"
            ))
            .unwrap(),
            "First line. Second line."
        );
        assert!(docs(&f("pub fn f(x:f64)->f64{x}")).is_err());
    }
}
