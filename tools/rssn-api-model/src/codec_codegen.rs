//! Source-derived composable codecs for reviewed, concrete kernel functions.
use quote::ToTokens;
use std::collections::BTreeSet;
use std::path::Path;
use syn::FnArg;
use syn::Item;
use syn::ReturnType;
use syn::Type;
use syn::Visibility;

fn compact(value: impl ToTokens) -> String {
    value
        .to_token_stream()
        .to_string()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}
fn ident(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .enumerate()
            .all(|(i, c)| c == '_' || c.is_ascii_alphabetic() || (i > 0 && c.is_ascii_digit()))
}
fn qualify(
    ty: &Type,
    module: &str,
    complex_alias: &BTreeSet<String>,
) -> String {
    match ty {
        | Type::Reference(r) => {
            format!(
                "&{}{}",
                if r.mutability.is_some() {
                    "mut "
                } else {
                    ""
                },
                qualify(&r.elem, module, complex_alias)
            )
        },
        | Type::Slice(s) => format!("[{}]", qualify(&s.elem, module, complex_alias)),
        | Type::Array(a) => {
            format!(
                "[{};{}]",
                qualify(&a.elem, module, complex_alias),
                compact(&a.len)
            )
        },
        | Type::Tuple(t) => {
            format!(
                "({},)",
                t.elems
                    .iter()
                    .map(|t| qualify(t, module, complex_alias))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        },
        | Type::Path(path) => {
            let text = compact(ty);
            match text.as_str() {
                | "Mat" | "dense::Mat" => "crate::kernels::dense::Mat".into(),
                | "String" => "String".into(),
                | "BigInt" => "num_bigint::BigInt".into(),
                | "Complex" | "Complex64" | "Complex<f64>"
                    if complex_alias.contains(text.trim_end_matches("<f64>")) =>
                {
                    "num_complex::Complex64".into()
                },
                | _ => {
                    if path.path.segments.len() == 1 {
                        let last = &path.path.segments[0];
                        if let syn::PathArguments::AngleBracketed(args) = &last.arguments {
                            let args = args
                                .args
                                .iter()
                                .map(|a| {
                                    match a {
                                        | syn::GenericArgument::Type(t) => {
                                            qualify(t, module, complex_alias)
                                        },
                                        | _ => compact(a),
                                    }
                                })
                                .collect::<Vec<_>>()
                                .join(",");
                            return format!("{}<{args}>", last.ident);
                        }
                        if text.chars().next().is_some_and(char::is_uppercase) {
                            return format!("crate::{module}::{text}");
                        }
                    }
                    text
                },
            }
        },
        | _ => compact(ty),
    }
}
fn collect_imports(
    tree: &syn::UseTree,
    prefix: &str,
    imports: &mut BTreeSet<(String, String)>,
) {
    match tree {
        | syn::UseTree::Path(path) => {
            collect_imports(&path.tree, &format!("{prefix}{}::", path.ident), imports)
        },
        | syn::UseTree::Name(name) => {
            imports.insert((name.ident.to_string(), format!("{prefix}{}", name.ident)));
        },
        | syn::UseTree::Rename(name) => {
            imports.insert((name.rename.to_string(), format!("{prefix}{}", name.ident)));
        },
        | syn::UseTree::Group(group) => {
            for item in &group.items {
                collect_imports(item, prefix, imports);
            }
        },
        | syn::UseTree::Glob(_) => {},
    }
}

fn owned(
    ty: &Type,
    module: &str,
    complex_alias: &BTreeSet<String>,
) -> Result<(String, bool), String> {
    if let Type::Reference(reference) = ty {
        if reference.mutability.is_some() {
            return Err("mutable arguments require a dedicated adapter".into());
        }
        if compact(&reference.elem) == "str" {
            return Ok(("String".into(), true));
        }
        if let Type::Slice(slice) = reference.elem.as_ref() {
            return Ok((
                format!("Vec<{}>", qualify(&slice.elem, module, complex_alias)),
                true,
            ));
        }
        return Ok((qualify(&reference.elem, module, complex_alias), true));
    }
    Ok((qualify(ty, module, complex_alias), false))
}
fn output(
    ty: &Type,
    module: &str,
    complex_alias: &BTreeSet<String>,
) -> Result<(String, bool), String> {
    if let Type::Path(path) = ty {
        if path.path.is_ident("Self") {
            return Err("Self output requires an associated adapter".into());
        }
        if let Some(segment) = path.path.segments.last()
            && segment.ident == "Result"
            && let syn::PathArguments::AngleBracketed(args) = &segment.arguments
            && let Some(syn::GenericArgument::Type(inner)) = args.args.first()
        {
            return Ok((
                format!("Result<{}, OpError>", qualify(inner, module, complex_alias)),
                true,
            ));
        }
    }
    Ok((qualify(ty, module, complex_alias), false))
}
/// Generate codec bindings. Rust checks each generated call and codec trait;
/// policies supply names and dimensions, never duplicate Rust signatures.
pub fn generate(
    root: &Path,
    policy: &str,
) -> Result<String, String> {
    let contract_source = std::fs::read_to_string(root.join("tools/rssn-api-gen/contracts.tsv"))
        .map_err(|e| e.to_string())?;
    let contracts = crate::contract_codegen::simple_summaries(&contract_source)?;
    let mut names = BTreeSet::new();
    let mut paths = BTreeSet::new();
    let mut rows = Vec::new();
    for (line, text) in policy.lines().enumerate() {
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = text.split('\t').collect();
        let [path, name, namespace, dimension] = fields.as_slice() else {
            return Err(format!(
                "codec policy line {} needs path, name, namespace, dimension",
                line + 1
            ));
        };
        let parts: Vec<_> = path.split("::").collect();
        if parts.len() != 4
            || parts[0] != "crate"
            || !matches!(parts[1], "kernels" | "backend")
            || !parts[2..].iter().all(|p| ident(p))
            || !ident(name)
            || !ident(namespace)
            || !ident(dimension)
            || !names.insert(*name)
            || !paths.insert(*path)
        {
            return Err(format!(
                "invalid or duplicate codec policy line {}",
                line + 1
            ));
        }
        if !matches!(
            *namespace,
            "numerical_vector" | "numerical_stats" | "numerical_special" | "numerical_matrix"
        ) {
            return Err(format!(
                "{path} uses a numerical namespace that is not installed"
            ));
        }
        let contract = contracts
            .get(*dimension)
            .ok_or_else(|| format!("{path} has unknown dimensional policy {dimension}"))?;
        let module = parts[2];
        let public_modules = std::fs::read_to_string(root.join(format!("src/{}/mod.rs", parts[1])))
            .map_err(|e| e.to_string())?;
        let module_items = syn::parse_file(&public_modules).map_err(|e| e.to_string())?;
        if !module_items.items.iter().any(|item| matches!(item,Item::Mod(m) if m.ident == module && matches!(m.vis,Visibility::Public(_)) && !m.attrs.iter().any(|a| a.path().is_ident("cfg") || a.path().is_ident("cfg_attr")))) { return Err(format!("{path} parent is not an unconditional public kernel module")); }
        let source = std::fs::read_to_string(root.join(format!("src/{}/{module}.rs", parts[1])))
            .map_err(|e| e.to_string())?;
        let file = syn::parse_file(&source).map_err(|e| e.to_string())?;
        let functions: Vec<_> = file
            .items
            .iter()
            .filter_map(|item| {
                match item {
                    | Item::Fn(f) if f.sig.ident == parts[3] => Some(f),
                    | _ => None,
                }
            })
            .collect();
        let [f] = functions.as_slice() else {
            return Err(format!("missing or ambiguous public function {path}"));
        };
        let mut imports = BTreeSet::new();
        for item in &file.items {
            if let Item::Use(item) = item {
                collect_imports(&item.tree, "", &mut imports);
            }
        }
        let complex_alias: BTreeSet<String> = imports
            .iter()
            .filter(|(_, target)| {
                matches!(
                    target.as_str(),
                    "num_complex::Complex64" | "num_complex::Complex"
                )
            })
            .map(|(name, _)| name.clone())
            .collect();
        let signature_text = compact(&f.sig);
        if file.items.iter().any(|item| matches!(item, Item::Type(alias) if alias.ident == "Complex" && signature_text.contains("Complex"))) {
            return Err(format!("{path} uses a local Complex alias; a dedicated adapter must define its representation"));
        }
        let sig = &f.sig;
        if !matches!(f.vis, Visibility::Public(_))
            || !sig.generics.params.is_empty()
            || sig.generics.where_clause.is_some()
            || sig.asyncness.is_some()
            || sig.unsafety.is_some()
            || sig.abi.is_some()
            || sig.variadic.is_some()
            || f.attrs
                .iter()
                .any(|a| a.path().is_ident("cfg") || a.path().is_ident("cfg_attr"))
        {
            return Err(format!("{path} requires a dedicated adapter"));
        }
        if contract
            .arguments
            .iter()
            .any(|index| *index >= sig.inputs.len())
        {
            return Err(format!(
                "{path} dimensional policy refers to a missing argument"
            ));
        }
        let mut args = Vec::new();
        let mut call = Vec::new();
        for (index, arg) in sig.inputs.iter().enumerate() {
            let FnArg::Typed(arg) = arg else {
                return Err(format!("{path} has a receiver"));
            };
            let (mut ty, borrow) =
                owned(&arg.ty, &format!("{}::{module}", parts[1]), &complex_alias)?;
            let nested_slices = ty == "Vec<&[f64]>";
            if nested_slices {
                ty = "Vec<Vec<f64>>".into();
            }
            args.push(format!("arg{index}: {ty}"));
            call.push(if nested_slices {
                format!("&arg{index}.iter().map(Vec::as_slice).collect::<Vec<_>>()")
            } else {
                format!("{}arg{index}", if borrow { "&" } else { "" })
            });
        }
        let ReturnType::Type(_, ret) = &sig.output else {
            return Err(format!("{path} has no output"));
        };
        let (ret, fallible) = output(ret, &format!("{}::{module}", parts[1]), &complex_alias)?;
        let doc = f
            .attrs
            .iter()
            .filter_map(|attr| {
                if attr.path().is_ident("doc")
                    && let syn::Meta::NameValue(v) = &attr.meta
                    && let syn::Expr::Lit(syn::ExprLit {
                        lit: syn::Lit::Str(s),
                        ..
                    }) = &v.value
                {
                    Some(s.value().trim().to_string())
                } else {
                    None
                }
            })
            .take_while(|s| !s.is_empty())
            .collect::<Vec<_>>()
            .join(" ");
        if doc.is_empty() {
            return Err(format!("{path} has no documentation"));
        }
        let error_mapping = if fallible {
            let ReturnType::Type(_, original) = &sig.output else {
                unreachable!()
            };
            let Type::Path(original) = original.as_ref() else {
                unreachable!()
            };
            let syn::PathArguments::AngleBracketed(arguments) = &original
                .path
                .segments
                .last()
                .ok_or("missing Result segment")?
                .arguments
            else {
                unreachable!()
            };
            match arguments.args.last() {
                | Some(syn::GenericArgument::Type(Type::Path(error)))
                    if error.path.is_ident("String") =>
                {
                    ".map_err(OpError::Domain)"
                },
                | Some(syn::GenericArgument::Type(Type::Path(error)))
                    if error
                        .path
                        .segments
                        .last()
                        .is_some_and(|segment| segment.ident == "DenseError") =>
                {
                    ".map_err(dense_error)"
                },
                | _ => ".map_err(|error| OpError::Domain(error.to_string()))",
            }
        } else {
            ""
        };
        let expression = format!("{path}({}){error_mapping}", call.join(", "));
        rows.push((
            *path,
            *name,
            *namespace,
            *dimension,
            args.join(", "),
            ret,
            expression,
            doc,
        ));
    }
    rows.sort_by_key(|r| r.0);
    let mut out = String::from(
        "// Generated by rssn-api-gen --codecs-write; edit codec-policy.tsv and kernel sources.\n",
    );
    for (path, name, namespace, dimension, args, ret, expression, doc) in &rows {
        out.push_str(&format!("binding!({name}, {name:?}, {namespace:?}, {doc:?}, {:?}, DimensionPolicy::{dimension}, ({args}) -> {ret}, |{args}| {{ {expression} }});\n",path));
    }
    out.push_str("pub(super) static GENERATED_CODEC_EXPORTS: &[&CodecExport] = &[\n");
    for (_, name, ..) in &rows {
        out.push_str(&format!("    &{name}::EXPORT,\n"));
    }
    out.push_str("];\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn codecs_derive_arguments_results_and_reject_private_paths() {
        let root = std::env::temp_dir().join(format!("rssn-codecs-{}", std::process::id()));
        std::fs::create_dir_all(root.join("src/kernels")).unwrap();
        std::fs::create_dir_all(root.join("tools/rssn-api-gen")).unwrap();
        std::fs::write(
            root.join("tools/rssn-api-gen/contracts.tsv"),
            "Dimensionless\t-\t-\t-\t-\tRaw.\n",
        )
        .unwrap();
        std::fs::write(root.join("src/kernels/mod.rs"), "pub mod example;").unwrap();
        std::fs::write(root.join("src/kernels/example.rs"),"/// Test.\npub fn fit(x: &[f64], n: usize) -> Result<Vec<f64>, String> { Ok(x.to_vec()) }\nfn private(x:f64)->f64{x}").unwrap();
        let generated = generate(
            &root,
            "crate::kernels::example::fit\tnum_fit\tnumerical_stats\tDimensionless",
        )
        .unwrap();
        assert!(generated.contains("arg0: Vec<f64>, arg1: usize"));
        assert!(generated.contains("fit(&arg0, arg1).map_err(OpError::Domain)"));
        assert!(
            generate(
                &root,
                "crate::kernels::example::private\tnum_private\tnumerical_stats\tDimensionless"
            )
            .is_err()
        );
        std::fs::write(root.join("src/kernels/mod.rs"), "mod example;").unwrap();
        assert!(
            generate(
                &root,
                "crate::kernels::example::fit\tnum_fit\tnumerical_stats\tDimensionless"
            )
            .is_err()
        );
        std::fs::write(root.join("src/kernels/mod.rs"), "pub mod example;").unwrap();
        std::fs::write(root.join("src/kernels/example.rs"), "pub type Complex = (f64, f64);\n/// Tuple transform.\npub fn fft(x: &[Complex]) -> Vec<Complex> { x.to_vec() }").unwrap();
        let ambiguous = generate(
            &root,
            "crate::kernels::example::fft\tnum_fft\tnumerical_vector\tDimensionless",
        )
        .unwrap_err();
        assert!(ambiguous.contains("local Complex alias"), "{ambiguous}");
        std::fs::write(root.join("src/kernels/example.rs"), "use num_complex::Complex64 as Complex;\n/// Complex transform.\npub fn fft(x: &[Complex]) -> Vec<Complex> { x.to_vec() }").unwrap();
        let generated = generate(
            &root,
            "crate::kernels::example::fft\tnum_fft\tnumerical_vector\tDimensionless",
        )
        .unwrap();
        assert!(generated.contains("arg0: Vec<num_complex::Complex64>"));
        assert_eq!(
            qualify(
                &syn::parse_str::<Type>("Complex").unwrap(),
                "kernels::example",
                &BTreeSet::new()
            ),
            "crate::kernels::example::Complex"
        );
        std::fs::remove_dir_all(root).unwrap();
    }
}
