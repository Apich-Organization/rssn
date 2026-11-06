//! Dimensional policies for numerical exports.
//!
//! Each policy generates parameter and result signatures, including record fields.
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::signature::TypeParam;
use crate::kernel::signature::TypeParamKind;

/// A composable, exact dimensional expression over named input witnesses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DimFormula {
    /// The dimensionless dimension.
    One,
    /// A dimension bound by an input path.
    Var(&'static str),
    /// The product of any number of dimensions.
    Mul(&'static [Self]),
    /// An exact rational power; generated denominators are positive.
    Pow(&'static Self, i64, i64),
}
impl DimFormula {
    pub(in crate::rules) fn instantiate(
        &self,
        vars: &[(&str, &DimExpr)],
    ) -> DimExpr {
        match self {
            | Self::One => DimExpr::Dimless,
            | Self::Var(name) => {
                vars.iter()
                    .find(|(key, _)| key == name)
                    .expect("validated formula variable")
                    .1
                    .clone()
            },
            | Self::Mul(parts) => {
                DimExpr::Mul(parts.iter().map(|part| part.instantiate(vars)).collect())
            },
            | Self::Pow(base, n, d) => base.instantiate(vars).pow(Exp::Rat(Rational::new(*n, *d))),
        }
    }
}
/// An explicitly reviewed field formula or a nonnumeric/unbound field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FieldDimension {
    /// The field has the stated exact dimensional formula.
    Formula(&'static DimFormula),
    /// A boolean or integer status/count, not a physical quantity.
    NonNumeric,
    /// No operator-specific dimension policy has been reviewed.
    Unbound,
}
/// Reviewed dimensional policy for a generated Rust result codec.
#[derive(Clone, Copy, Debug)]
pub struct StructuredPolicy {
    /// Public Rust result type.
    pub rust_path: &'static str,
    /// Named operation-specific interpretation of this Rust result type.
    pub contract: &'static str,
    /// Whether all numeric fields have explicit formulas.
    pub bound: bool,
    /// Exhaustive fields and their reviewed dimensions.
    pub fields: &'static [(&'static str, FieldDimension)],
}

/// One reviewed dimension witness at a numerical argument and record path.
#[derive(Clone, Copy, Debug)]
pub struct InputBinding {
    /// Contract variable name.
    pub name: &'static str,
    /// Positional argument index.
    pub argument: usize,
    /// Named record fields, from outermost to innermost.
    pub path: &'static [&'static str],
}
/// A declared dimension on an additional numerical input.
#[derive(Clone, Copy, Debug)]
pub struct InputConstraint {
    /// Positional argument index.
    pub argument: usize,
    /// Nested record field path.
    pub path: &'static [&'static str],
    /// Required dimension.
    pub formula: &'static DimFormula,
}
/// A dimension assigned to a nested result record path.
#[derive(Clone, Copy, Debug)]
pub struct OutputBinding {
    /// Named result fields, from outermost to innermost.
    pub path: &'static [&'static str],
    /// The result dimension.
    pub formula: &'static DimFormula,
}
/// One operation's reviewed input and result dimensional relationship.
#[derive(Clone, Copy, Debug)]
pub struct DimContract {
    /// Variables bound by input witnesses.
    pub inputs: &'static [InputBinding],
    /// Other dimensioned inputs and their required formulas.
    pub constraints: &'static [InputConstraint],
    /// Exhaustive named result fields.
    pub outputs: &'static [(&'static str, FieldDimension)],
    /// Additional nested result paths for structured codecs.
    pub output_paths: &'static [OutputBinding],
    /// Formula for a scalar or tensor result, if physically dimensioned.
    pub result: Option<&'static DimFormula>,
    /// Require all unbound arguments to carry no physical dimension.
    pub other_args_dimensionless: bool,
}
fn path_type_mut<'a>(
    sig: &'a mut Signature,
    arg: usize,
    path: &[&str],
) -> Option<&'a mut TypeExpr> {
    let mut current = sig.params.get_mut(arg)?;
    for name in path {
        let TypeExpr::Record(fields) = current else {
            return None;
        };
        let index = fields.iter().position(|(field, _)| field == name)?;
        if fields.iter().filter(|(field, _)| field == name).count() != 1 {
            return None;
        }
        current = &mut fields[index].1;
    }
    Some(current)
}
fn numeric_leaf(ty: &TypeExpr) -> bool {
    match ty {
        | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => numeric_leaf(elem),
        | TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex => true,
        | _ => false,
    }
}
fn nonnumeric_leaf(ty: &TypeExpr) -> bool {
    match ty {
        | TypeExpr::Tensor { elem, .. } => nonnumeric_leaf(elem),
        | TypeExpr::Int | TypeExpr::Bool => true,
        | _ => false,
    }
}
/// Derive a manifest signature from the same contract used at runtime.
pub(in crate::rules) fn apply_contract(
    mut sig: Signature,
    contract: &DimContract,
) -> Signature {
    let mut names = std::collections::BTreeSet::new();
    let mut paths = std::collections::BTreeSet::new();
    for input in contract.inputs {
        assert!(
            names.insert(input.name),
            "duplicate contract input variable"
        );
        assert!(
            paths.insert((input.argument, input.path)),
            "duplicate contract input path"
        );
    }
    for check in contract.constraints {
        assert!(
            paths.insert((check.argument, check.path)),
            "duplicate contract input or constraint path"
        );
    }
    let vars: Vec<_> = contract
        .inputs
        .iter()
        .map(|input| (input.name, DimExpr::Var(input.name.into())))
        .collect();
    let refs: Vec<_> = vars.iter().map(|(name, dim)| (*name, dim)).collect();
    for (input, (_, dimension)) in contract.inputs.iter().zip(&vars) {
        sig.type_params.push(TypeParam {
            name: input.name.into(),
            kind: TypeParamKind::Dim,
            bounds: vec![],
        });
        let target =
            path_type_mut(&mut sig, input.argument, input.path).expect("contract input path");
        assert!(
            numeric_leaf(target),
            "contract input must be an unannotated numeric leaf"
        );
        attach(target, dimension.clone());
    }
    for check in contract.constraints {
        let target =
            path_type_mut(&mut sig, check.argument, check.path).expect("contract constraint path");
        assert!(
            numeric_leaf(target),
            "contract constraint must target an unannotated numeric leaf"
        );
        let dim = check.formula.instantiate(&refs);
        if dim != DimExpr::Dimless {
            attach(target, dim);
        }
    }
    let mut output_names = std::collections::BTreeSet::new();
    if !contract.outputs.is_empty() {
        let TypeExpr::Record(fields) = &mut sig.ret else {
            panic!("contract outputs require a record");
        };
        assert_eq!(
            fields.len(),
            contract.outputs.len(),
            "contract must cover each result field"
        );
        for (name, field_policy) in contract.outputs {
            assert!(
                output_names.insert(*name),
                "duplicate contract result field"
            );
            let mut matches = fields.iter_mut().filter(|(field, _)| field == name);
            let target = &mut matches.next().expect("contract result field").1;
            assert!(matches.next().is_none(), "duplicate schema result field");
            match field_policy {
                | FieldDimension::Formula(formula) => {
                    assert!(
                        numeric_leaf(target),
                        "contract formula requires numeric result"
                    );
                    let dim = formula.instantiate(&refs);
                    if dim != DimExpr::Dimless {
                        attach(target, dim);
                    }
                },
                | FieldDimension::NonNumeric => {
                    assert!(
                        nonnumeric_leaf(target),
                        "contract status field must be integer or bool"
                    );
                },
                | FieldDimension::Unbound => panic!("unbound contract result"),
            }
        }
    }
    let mut nested_paths = std::collections::BTreeSet::new();
    for output in contract.output_paths {
        assert!(
            nested_paths.insert(output.path),
            "duplicate contract output path"
        );
        assert!(!output.path.is_empty(), "empty nested contract output path");
        assert!(
            !output_names.contains(output.path[0]),
            "overlapping contract output paths"
        );
        let mut current = &mut sig.ret;
        for name in output.path {
            let TypeExpr::Record(fields) = current else {
                panic!("contract output path requires a record");
            };
            let index = fields
                .iter()
                .position(|(field, _)| field == name)
                .expect("contract output field");
            assert_eq!(
                fields.iter().filter(|(field, _)| field == name).count(),
                1,
                "duplicate contract output field"
            );
            current = &mut fields[index].1;
        }
        assert!(
            numeric_leaf(current),
            "contract output must target an unannotated numeric leaf"
        );
        let dimension = output.formula.instantiate(&refs);
        if dimension != DimExpr::Dimless {
            attach(current, dimension);
        }
    }
    if let Some(result) = contract.result {
        assert!(
            numeric_leaf(&sig.ret),
            "contract scalar result must be numeric"
        );
        let dim = result.instantiate(&refs);
        if dim != DimExpr::Dimless {
            attach(&mut sig.ret, dim);
        }
    }
    sig
}

fn attach(
    ty: &mut TypeExpr,
    dim: DimExpr,
) {
    match ty {
        | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => attach(elem, dim),
        | TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex => {
            let scalar = std::mem::replace(ty, TypeExpr::Never);
            *ty = TypeExpr::Dim(Box::new(scalar), dim);
        },
        | _ => {},
    }
}
fn field(
    ty: &mut TypeExpr,
    name: &str,
    dim: DimExpr,
) {
    if let TypeExpr::Record(fields) = ty
        && let Some((_, t)) = fields.iter_mut().find(|(field, _)| field == name)
    {
        attach(t, dim);
    }
}
pub(in crate::rules) fn apply_record_fields(
    ty: &mut TypeExpr,
    fields: &[(&str, FieldDimension)],
    d: &DimExpr,
    e: &DimExpr,
) {
    for (name, formula) in fields {
        let dim = match formula {
            | FieldDimension::Formula(formula) => formula.instantiate(&[("D", d), ("E", e)]),
            | FieldDimension::NonNumeric => continue,
            | FieldDimension::Unbound => {
                panic!("unbound structured dimensional policy cannot be registered")
            },
        };
        if dim != DimExpr::Dimless {
            field(ty, name, dim);
        }
    }
}

include!("contract_generated.rs");

/// Apply a descriptor's policy to its shape-derived signature.
#[must_use]
pub fn apply(
    sig: Signature,
    policy: DimensionPolicy,
) -> Signature {
    apply_contract(sig, contract_for(policy))
}

#[cfg(test)]
mod contract_validation_tests {
    use super::*;
    fn signature() -> Signature {
        Signature::new(
            vec![TypeExpr::Float],
            TypeExpr::Record(vec![("value".into(), TypeExpr::Float)]),
        )
    }
    #[test]
    fn stale_dimensionless_paths_and_duplicate_bindings_are_rejected() {
        let bad_input = DimContract {
            inputs: &[],
            constraints: &[InputConstraint {
                argument: 1,
                path: &[],
                formula: &DimFormula::One,
            }],
            outputs: &[],
            output_paths: &[],
            result: None,
            other_args_dimensionless: false,
        };
        assert!(std::panic::catch_unwind(|| apply_contract(signature(), &bad_input)).is_err());
        let bad_output = DimContract {
            inputs: &[],
            constraints: &[],
            outputs: &[("missing", FieldDimension::Formula(&DimFormula::One))],
            output_paths: &[],
            result: None,
            other_args_dimensionless: false,
        };
        assert!(std::panic::catch_unwind(|| apply_contract(signature(), &bad_output)).is_err());
        let duplicate = DimContract {
            inputs: &[
                InputBinding {
                    name: "D",
                    argument: 0,
                    path: &[],
                },
                InputBinding {
                    name: "E",
                    argument: 0,
                    path: &[],
                },
            ],
            constraints: &[],
            outputs: &[],
            output_paths: &[],
            result: None,
            other_args_dimensionless: false,
        };
        assert!(std::panic::catch_unwind(|| apply_contract(signature(), &duplicate)).is_err());
    }
}
