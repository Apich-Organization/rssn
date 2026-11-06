//! Draft structured codec generator; no result type becomes callable merely
//! because its public Rust fields can be encoded.
use quote::ToTokens;
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};
use syn::{Fields, Item, Visibility};

fn compact(tokens: impl ToTokens) -> String {
    tokens
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
// Parse reviewed formulas as an expression tree, not an enumerated list of
// special cases. Every name must occur in the contract variable declaration.
pub(crate) fn formula(
    source: &str,
    variables: &BTreeSet<&str>,
) -> Result<String, String> {
    fn integer(expr: &syn::Expr) -> Option<i64> {
        match expr {
            | syn::Expr::Lit(syn::ExprLit {
                lit: syn::Lit::Int(n),
                ..
            }) => n.base10_parse().ok(),
            | syn::Expr::Unary(syn::ExprUnary {
                op: syn::UnOp::Neg(_),
                expr,
                ..
            }) => integer(expr)?.checked_neg(),
            | _ => None,
        }
    }
    fn emit(
        expr: &syn::Expr,
        variables: &BTreeSet<&str>,
    ) -> Result<String, String> {
        let prefix = "super::dimensions::DimFormula";
        match expr {
            | syn::Expr::Paren(value) => emit(&value.expr, variables),
            | syn::Expr::Lit(syn::ExprLit {
                lit: syn::Lit::Int(n),
                ..
            }) if n.base10_parse::<u64>().ok() == Some(1) => Ok(format!("{prefix}::One")),
            | syn::Expr::Path(path)
                if path.path.segments.len() == 1
                    && variables.contains(path.path.segments[0].ident.to_string().as_str()) =>
            {
                Ok(format!(
                    "{prefix}::Var({:?})",
                    path.path.segments[0].ident.to_string()
                ))
            },
            | syn::Expr::Binary(binary)
                if matches!(binary.op, syn::BinOp::Mul(_) | syn::BinOp::Div(_)) =>
            {
                let left = emit(&binary.left, variables)?;
                let mut right = emit(&binary.right, variables)?;
                if matches!(binary.op, syn::BinOp::Div(_)) {
                    right = format!("{prefix}::Pow(&{right}, -1, 1)");
                }
                Ok(format!("{prefix}::Mul(&[{left}, {right}])"))
            },
            | syn::Expr::Binary(binary) if matches!(binary.op, syn::BinOp::BitXor(_)) => {
                let power = integer(&binary.right).ok_or("power must be an integer")?;
                Ok(format!(
                    "{prefix}::Pow(&{}, {power}, 1)",
                    emit(&binary.left, variables)?
                ))
            },
            | syn::Expr::Call(call) => {
                let syn::Expr::Path(name) = call.func.as_ref() else {
                    return Err("invalid formula function".into());
                };
                if name.path.is_ident("sqrt") && call.args.len() == 1 {
                    Ok(format!(
                        "{prefix}::Pow(&{}, 1, 2)",
                        emit(&call.args[0], variables)?
                    ))
                } else if name.path.is_ident("pow") && call.args.len() == 3 {
                    let numerator =
                        integer(&call.args[1]).ok_or("power numerator must be an integer")?;
                    let denominator = integer(&call.args[2])
                        .filter(|d| *d > 0)
                        .ok_or("power denominator must be positive")?;
                    Ok(format!(
                        "{prefix}::Pow(&{}, {numerator}, {denominator})",
                        emit(&call.args[0], variables)?
                    ))
                } else {
                    Err("unsupported formula function or arity".into())
                }
            },
            | _ => Err("unsupported formula or unbound dimension variable".into()),
        }
    }
    let expr = syn::parse_str::<syn::Expr>(source).map_err(|error| error.to_string())?;
    emit(&expr, variables)
}
/// Generate exact public record codecs plus auditable per-field policies.
/// Bound rows must cover every numeric leaf; unbound rows are codecs only and
/// can never claim a callable dimensional contract.
pub fn generate(
    root: &Path,
    policy: &str,
) -> Result<String, String> {
    let mut contracts = BTreeSet::new();
    let mut rows = Vec::new();
    for (line, text) in policy.lines().enumerate() {
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = text.split('\t').collect();
        let [path, contract, status, vars, field_policy] = fields.as_slice() else {
            return Err(format!(
                "line {} needs type, contract, status, variables and field policy",
                line + 1
            ));
        };
        let mut variables = BTreeSet::new();
        if *vars != "-" {
            for variable in vars.split(',') {
                if !ident(variable) || !variables.insert(variable) {
                    return Err(format!("invalid or duplicate variable {variable}"));
                }
            }
        }
        if !ident(contract) || !contracts.insert(contract.to_ascii_uppercase()) {
            return Err(format!("duplicate or invalid contract {contract}"));
        }
        let parts: Vec<_> = path.split("::").collect();
        if parts.len() != 4
            || parts[0] != "crate"
            || parts[1] != "kernels"
            || !parts[2..].iter().all(|p| ident(p))
            || !matches!(*status, "bound" | "unbound")
        {
            return Err(format!("invalid type or status on line {}", line + 1));
        }
        let source = root.join("src/kernels").join(format!("{}.rs", parts[2]));
        let text =
            std::fs::read_to_string(&source).map_err(|e| format!("{}: {e}", source.display()))?;
        let file = syn::parse_file(&text).map_err(|e| format!("{}: {e}", source.display()))?;
        let mut matches = file.items.iter().filter_map(|item| match item {
            | Item::Struct(item) if item.ident == parts[3] => Some(item),
            | _ => None,
        });
        let item = matches
            .next()
            .ok_or_else(|| format!("missing public struct {path}"))?;
        if matches.next().is_some()
            || !matches!(item.vis, Visibility::Public(_))
            || !item.generics.params.is_empty()
            || item.generics.where_clause.is_some()
            || item
                .attrs
                .iter()
                .any(|a| a.path().is_ident("cfg") || a.path().is_ident("cfg_attr"))
        {
            return Err(format!("unsupported struct {path}"));
        }
        let Fields::Named(fields) = &item.fields else {
            return Err(format!("{path} has non-named fields"));
        };
        let mut policy_fields = BTreeMap::new();
        for token in field_policy.split(';') {
            let (name, dim) = token
                .split_once('=')
                .ok_or_else(|| format!("invalid field policy {token}"))?;
            if !ident(name) || policy_fields.insert(name, dim).is_some() {
                return Err(format!("duplicate or invalid field {name}"));
            }
        }
        let mut encoded = Vec::new();
        for field in &fields.named {
            if !matches!(field.vis, Visibility::Public(_))
                || field
                    .attrs
                    .iter()
                    .any(|a| a.path().is_ident("cfg") || a.path().is_ident("cfg_attr"))
            {
                return Err(format!("{path} contains private or conditional field"));
            }
            let name = field
                .ident
                .as_ref()
                .ok_or_else(|| format!("unnamed field in {path}"))?
                .to_string();
            let ty = compact(&field.ty);
            let (codec, numeric) = match ty.as_str() {
                | "f64" | "Vec<f64>" => (ty, true),
                | "usize" | "bool" | "Vec<usize>" => (ty, false),
                | "Mat" if parts[2] == "dense" => ("crate::kernels::dense::Mat".to_string(), true),
                | _ => return Err(format!("{path}.{name} has unsupported field type {ty}")),
            };
            let dim = policy_fields
                .remove(name.as_str())
                .ok_or_else(|| format!("{path}.{name} has no dimension policy"))?;
            let valid = if !numeric {
                dim == "-"
            } else if *status == "unbound" {
                dim == "?"
            } else {
                formula(dim, &variables).is_ok()
            };
            if !valid {
                return Err(format!(
                    "{path}.{name} has invalid {status} dimensional policy {dim}"
                ));
            }
            encoded.push((name, codec, dim.to_string()));
        }
        if !policy_fields.is_empty() {
            return Err(format!("{path} contains stale field policies"));
        }
        rows.push((*path, *contract, *status, variables, encoded));
    }
    rows.sort_by_key(|row| (row.0, row.1));
    let mut out = String::from(
        "// Generated structured record codecs; source policy does not register operators.\n",
    );
    let mut policy_rows = Vec::new();
    let mut emitted_types = BTreeSet::new();
    for (path, contract, status, variables, fields) in rows {
        out.push_str(&format!(
            "// {path}: {status}; {}\n",
            fields
                .iter()
                .map(|(n, _, d)| format!("{n}={d}"))
                .collect::<Vec<_>>()
                .join(";")
        ));
        if emitted_types.insert(path) {
            out.push_str(&format!(
                "record!({path} {{ {} }});\n",
                fields
                    .iter()
                    .map(|(n, t, _)| format!("{n}: {t}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        let name = contract.to_ascii_uppercase();
        let policies = fields
            .iter()
            .map(|(name, _, dimension)| {
                let value = match dimension.as_str() {
                    | "-" => "super::dimensions::FieldDimension::NonNumeric".to_string(),
                    | "?" => "super::dimensions::FieldDimension::Unbound".to_string(),
                    | _ => format!(
                        "super::dimensions::FieldDimension::Formula(&{})",
                        formula(dimension, &variables).expect("validated formula")
                    ),
                };
                format!("({name:?}, {value})")
            })
            .collect::<Vec<_>>()
            .join(", ");
        out.push_str(&format!("pub(in crate::rules) const STRUCTURED_{name}: &[(&str, super::dimensions::FieldDimension)] = &[{policies}];\n"));
        policy_rows.push(format!("super::dimensions::StructuredPolicy {{ rust_path: {path:?}, contract: {contract:?}, bound: {}, fields: STRUCTURED_{name} }}", status == "bound"));
    }
    out.push_str("/// Source-derived record schemas and reviewed dimensional status.\n");
    out.push_str(&format!(
        "pub static STRUCTURED_POLICIES: &[super::dimensions::StructuredPolicy] = &[{}];\n",
        policy_rows.join(", ")
    ));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dimensional_formulas_are_composable_and_reject_unbound_or_invalid_powers() {
        let variables = BTreeSet::from(["D", "E", "velocity"]);
        let expr = formula("sqrt(D*E)/pow(velocity,1,3)", &variables).unwrap();
        assert!(expr.contains("::Mul("));
        assert!(expr.contains("::Pow("));
        assert!(expr.contains("Var(\"D\")"));
        assert!(formula("D*unknown", &variables).is_err());
        assert!(formula("pow(D,1,0)", &variables).is_err());
        assert!(formula("pow(D,1,-2)", &variables).is_err());
    }
    #[test]
    fn dimensions_cover_every_public_field_and_schema_drift_fails() {
        let root =
            std::env::temp_dir().join(format!("rssn-structured-audit-{}", std::process::id()));
        let dir = root.join("src/kernels");
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("sample.rs"),
            "pub struct Report { pub value: f64, pub ok: bool }\n",
        )
        .unwrap();
        let policy = "crate::kernels::sample::Report\tREPORT\tbound\tD,E\tvalue=D;ok=-\n";
        let generated = generate(&root, policy).unwrap();
        assert!(
            generated.contains("record!(crate::kernels::sample::Report { value: f64, ok: bool })")
        );
        assert!(generated.contains("DimFormula::Var(\"D\")"));
        let two =
            format!("{policy}crate::kernels::sample::Report\tOTHER\tbound\tD,E\tvalue=1;ok=-\n");
        let multiple = generate(&root, &two).unwrap();
        assert_eq!(
            multiple
                .matches("record!(crate::kernels::sample::Report")
                .count(),
            1
        );
        assert!(multiple.contains("STRUCTURED_REPORT"));
        assert!(multiple.contains("STRUCTURED_OTHER"));
        let collision = format!("{policy}crate::kernels::sample::Report\treport\tbound\tD,E\tvalue=1;ok=-\n");
        assert!(generate(&root, &collision).is_err());
        assert!(
            generate(
                &root,
                "crate::kernels::sample::Report\tREPORT\tbound\tD,E\tvalue=?;ok=-\n"
            )
            .is_err()
        );
        assert!(
            generate(
                &root,
                "crate::kernels::sample::Report\tREPORT\tbound\tD,E\tvalue=D\n"
            )
            .is_err()
        );
        std::fs::write(
            dir.join("sample.rs"),
            "pub struct Report { pub value: f64, pub ok: bool, pub count: usize }\n",
        )
        .unwrap();
        assert!(generate(&root, policy).is_err());
        std::fs::write(
            dir.join("sample.rs"),
            "pub struct Report { pub value: f64, ok: bool }\n",
        )
        .unwrap();
        assert!(generate(&root, policy).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }
}
