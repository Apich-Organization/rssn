//! Runtime checking and value witnesses for numerical dimensional policies.
use super::dimensions::DimensionPolicy;
use crate::kernel::call::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::DynValue;
use crate::kernel::value::Value;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use std::collections::BTreeMap;

type Dimension = BTreeMap<String, BigRational>;
pub(in crate::rules) fn canonical(expr: &DimExpr) -> Result<Dimension, OpError> {
    crate::kernel::overload::canonical_concrete_dimension(expr).ok_or_else(|| {
        OpError::Unsupported {
            backend: "numerical".into(),
            feature: "unknown or non-scalar dimension".into(),
        }
    })
}
pub(in crate::rules) fn expression(dim: &Dimension) -> Result<DimExpr, OpError> {
    let mut parts = Vec::new();
    for (name, power) in dim {
        let base = DimExpr::Base(name.clone());
        if power.is_one() {
            parts.push(base);
            continue;
        }
        let numer = power
            .numer()
            .to_i64()
            .ok_or_else(|| OpError::Domain("dimension exponent overflows ABI".into()))?;
        let denom = power
            .denom()
            .to_i64()
            .ok_or_else(|| OpError::Domain("dimension exponent overflows ABI".into()))?;
        parts.push(base.pow(Exp::Rat(Rational::new(numer, denom))));
    }
    Ok(match parts.len() {
        | 0 => DimExpr::Dimless,
        | 1 => parts.remove(0),
        | _ => DimExpr::Mul(parts),
    })
}
fn dimension(value: &Value) -> Result<Dimension, OpError> {
    fn schema_dimension(ty: &TypeExpr) -> Result<Dimension, OpError> {
        match ty {
            | TypeExpr::Dim(_, d) | TypeExpr::Point(d) => canonical(d),
            | TypeExpr::Tensor { elem, .. } => schema_dimension(elem),
            | _ => Ok(Dimension::new()),
        }
    }
    match value {
        | Value::Quantity(q) => canonical(&q.unit.dimension),
        | Value::Dyn(d) => canonical(&d.dim),
        | Value::Variant(v) if v.tag == "None" && v.payload.is_empty() => Ok(Dimension::new()),
        | Value::Variant(v) if v.tag == "Some" && v.payload.len() == 1 => dimension(&v.payload[0]),
        | Value::Tensor(t) => {
            let declared = schema_dimension(&t.elem)?;
            for item in &t.data {
                if matches!(
                    item,
                    Value::Quantity(_) | Value::Dyn(_) | Value::Tensor(_) | Value::Variant(_)
                ) && dimension(item)? != declared
                {
                    return Err(OpError::Domain(
                        "tensor element dimension disagrees with its declared element type".into(),
                    ));
                }
            }
            Ok(declared)
        },
        | _ => Ok(Dimension::new()),
    }
}

fn substitute(
    expr: &DimExpr,
    vars: &BTreeMap<&str, Dimension>,
) -> Result<DimExpr, OpError> {
    match expr {
        | DimExpr::Var(name) => {
            expression(vars.get(name.as_str()).ok_or_else(|| {
                OpError::Unsupported {
                    backend: "numerical".into(),
                    feature: format!("unbound dimension {name}"),
                }
            })?)
        },
        | DimExpr::Mul(parts) => {
            let items = parts
                .iter()
                .map(|p| substitute(p, vars))
                .collect::<Result<Vec<_>, _>>()?;
            expression(&canonical(&DimExpr::Mul(items))?)
        },
        | DimExpr::Pow(base, exp) => {
            let base = substitute(base, vars)?;
            expression(&canonical(&DimExpr::Pow(Box::new(base), exp.clone()))?)
        },
        | DimExpr::Inverse(base) => {
            let base = substitute(base, vars)?;
            expression(&canonical(&DimExpr::Inverse(Box::new(base)))?)
        },
        | other => Ok(other.clone()),
    }
}
fn instantiated_type(
    ty: &TypeExpr,
    vars: &BTreeMap<&str, Dimension>,
) -> Result<TypeExpr, OpError> {
    Ok(match ty {
        | TypeExpr::Dim(inner, dim) => {
            let inner = instantiated_type(inner, vars)?;
            let dim = substitute(dim, vars)?;
            if dim == DimExpr::Dimless {
                inner
            } else {
                TypeExpr::Dim(Box::new(inner), dim)
            }
        },
        | TypeExpr::Tensor { rank, shape, elem } => {
            TypeExpr::Tensor {
                rank: rank.clone(),
                shape: shape.clone(),
                elem: Box::new(instantiated_type(elem, vars)?),
            }
        },
        | TypeExpr::Record(fields) => {
            TypeExpr::Record(
                fields
                    .iter()
                    .map(|(name, ty)| Ok((name.clone(), instantiated_type(ty, vars)?)))
                    .collect::<Result<_, OpError>>()?,
            )
        },
        | TypeExpr::Variant(variants) => {
            TypeExpr::Variant(
                variants
                    .iter()
                    .map(|(tag, payload)| {
                        Ok((
                            tag.clone(),
                            payload
                                .iter()
                                .map(|ty| instantiated_type(ty, vars))
                                .collect::<Result<_, OpError>>()?,
                        ))
                    })
                    .collect::<Result<_, OpError>>()?,
            )
        },
        | other => other.clone(),
    })
}
fn scalar_representation_matches(
    value: &Value,
    ty: &TypeExpr,
) -> bool {
    match value {
        | Value::Dyn(dynamic) => scalar_representation_matches(&dynamic.value, ty),
        | Value::Quantity(quantity) => scalar_representation_matches(&quantity.value, ty),
        | _ => {
            matches!(
                (ty, value),
                (TypeExpr::Int, Value::Int(_))
                    | (
                        TypeExpr::Real,
                        Value::Int(_) | Value::Real(_) | Value::Float(_)
                    )
                    | (TypeExpr::Float, Value::Float(_))
                    | (
                        TypeExpr::Complex,
                        Value::Int(_) | Value::Real(_) | Value::Float(_) | Value::Complex(_)
                    )
                    | (TypeExpr::Bool, Value::Bool(_))
                    | (TypeExpr::String, Value::String(_))
            )
        },
    }
}
fn schema_representation_matches(
    expected: &TypeExpr,
    actual: &TypeExpr,
) -> bool {
    use crate::kernel::signature::NatExpr;
    match (expected, actual) {
        | (TypeExpr::Dim(e, ed), TypeExpr::Dim(a, ad)) => {
            canonical(ed)
                .ok()
                .zip(canonical(ad).ok())
                .is_some_and(|(left, right)| left == right)
                && schema_representation_matches(e, a)
        },
        | (TypeExpr::Dim(inner, _), other) => schema_representation_matches(inner, other),
        | (other, TypeExpr::Dim(inner, dim)) => {
            canonical(dim).is_ok_and(|dimension| dimension.is_empty())
                && schema_representation_matches(other, inner)
        },
        | (
            TypeExpr::Tensor {
                rank: NatExpr::Lit(er),
                shape: ea,
                elem: e,
            },
            TypeExpr::Tensor {
                rank: NatExpr::Lit(ar),
                shape: aa,
                elem: a,
            },
        ) => er == ar && (ea.is_none() || ea == aa) && schema_representation_matches(e, a),
        | (TypeExpr::Record(a), TypeExpr::Record(b)) => {
            let left: std::collections::BTreeSet<_> = a.iter().map(|(name, _)| name).collect();
            let right: std::collections::BTreeSet<_> = b.iter().map(|(name, _)| name).collect();
            left.len() == a.len()
                && right.len() == b.len()
                && left == right
                && a.iter().all(|(name, ty)| {
                    b.iter().filter(|(n, _)| n == name).count() == 1
                        && b.iter()
                            .any(|(n, other)| n == name && schema_representation_matches(ty, other))
                })
        },
        | (TypeExpr::Variant(expected), TypeExpr::Variant(actual)) => {
            let expected_tags: std::collections::BTreeSet<_> =
                expected.iter().map(|(tag, _)| tag).collect();
            let actual_tags: std::collections::BTreeSet<_> =
                actual.iter().map(|(tag, _)| tag).collect();
            !expected.is_empty()
                && !actual.is_empty()
                && expected_tags.len() == expected.len()
                && actual_tags.len() == actual.len()
                && actual.iter().all(|(tag, payload)| {
                    expected
                        .iter()
                        .find(|(name, _)| name == tag)
                        .is_some_and(|(_, schema)| {
                            schema.len() == payload.len()
                                && schema.iter().zip(payload).all(|(expected, actual)| {
                                    schema_representation_matches(expected, actual)
                                })
                        })
                })
        },
        | (TypeExpr::Real, TypeExpr::Int | TypeExpr::Real | TypeExpr::Float)
        | (
            TypeExpr::Complex,
            TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex,
        )
        | (TypeExpr::Int, TypeExpr::Int)
        | (TypeExpr::Float, TypeExpr::Float)
        | (TypeExpr::Bool, TypeExpr::Bool)
        | (TypeExpr::String, TypeExpr::String) => true,
        | _ => false,
    }
}
fn decorate(
    value: Value,
    ty: &TypeExpr,
) -> Result<Value, OpError> {
    match (value, ty) {
        | (value, TypeExpr::Dim(inner, dim)) => {
            if !scalar_representation_matches(&value, inner) {
                return Err(OpError::Domain(
                    "numerical result has an incompatible scalar representation".into(),
                ));
            }
            if matches!(value, Value::Dyn(_) | Value::Quantity(_)) {
                if dimension(&value)? != canonical(dim)? {
                    return Err(OpError::Domain(
                        "numerical result has a conflicting dimension witness".into(),
                    ));
                }
                return Ok(value);
            }
            if canonical(dim)?.is_empty() {
                return Ok(value);
            }
            Ok(Value::Dyn(DynValue {
                value: Box::new(value),
                dim: dim.clone(),
            }))
        },
        | (
            value @ (Value::Dyn(_) | Value::Quantity(_)),
            ty @ (TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex),
        ) => {
            if !scalar_representation_matches(&value, ty) {
                return Err(OpError::Domain(
                    "numerical result has an incompatible scalar representation".into(),
                ));
            }
            if !dimension(&value)?.is_empty() {
                return Err(OpError::Domain(
                    "numerical result has a conflicting dimension witness".into(),
                ));
            }
            Ok(value)
        },
        | (Value::Tensor(mut tensor), TypeExpr::Tensor { rank, shape, elem }) => {
            use crate::kernel::signature::NatExpr;
            if !matches!(rank, NatExpr::Lit(n) if usize::try_from(*n).ok() == Some(tensor.shape.len()))
                || tensor
                    .shape
                    .iter()
                    .try_fold(1usize, |a, n| a.checked_mul(*n))
                    != Some(tensor.data.len())
            {
                return Err(OpError::Domain(
                    "numerical tensor result has invalid rank or data length".into(),
                ));
            }
            if let Some(expected) = shape
                && (expected.len() != tensor.shape.len() || expected.iter().zip(&tensor.shape).any(|(e,a)| !matches!(e,NatExpr::Lit(n) if usize::try_from(*n).ok() == Some(*a)))) {
                return Err(OpError::Domain("numerical tensor result has incompatible extents".into()));
            }
            if !schema_representation_matches(elem, &tensor.elem) {
                return Err(OpError::Domain(
                    "numerical tensor result has incompatible element schema".into(),
                ));
            }
            if let TypeExpr::Dim(_, declared) = &tensor.elem {
                let expected = match elem.as_ref() {
                    | TypeExpr::Dim(_, dim) => canonical(dim)?,
                    | _ => Dimension::new(),
                };
                if canonical(declared)? != expected {
                    return Err(OpError::Domain(
                        "numerical tensor result has a conflicting element dimension".into(),
                    ));
                }
            }
            tensor.data = tensor
                .data
                .into_iter()
                .map(|value| decorate(value, elem))
                .collect::<Result<_, _>>()?;
            tensor.elem = (**elem).clone();
            Ok(Value::Tensor(tensor))
        },
        | (Value::Record(mut record), TypeExpr::Record(fields)) => {
            let declared: std::collections::BTreeSet<_> =
                fields.iter().map(|(name, _)| name).collect();
            let actual: std::collections::BTreeSet<_> =
                record.fields.iter().map(|(name, _)| name).collect();
            if declared.len() != fields.len()
                || actual.len() != record.fields.len()
                || actual != declared
            {
                return Err(OpError::Domain(
                    "numerical result record fields do not match the declared schema".into(),
                ));
            }
            for (name, value) in &mut record.fields {
                let ty = fields
                    .iter()
                    .find(|(field, _)| field == name)
                    .map(|(_, ty)| ty)
                    .ok_or_else(|| {
                        OpError::Internal(format!("undeclared numerical result field {name}"))
                    })?;
                *value = decorate(value.clone(), ty)?;
            }
            Ok(Value::Record(record))
        },
        | (Value::Variant(mut variant), TypeExpr::Variant(variants)) => {
            let tags: std::collections::BTreeSet<_> = variants.iter().map(|(tag, _)| tag).collect();
            if tags.len() != variants.len() {
                return Err(OpError::Domain(
                    "numerical variant result schema contains duplicate tags".into(),
                ));
            }
            let schema = variants
                .iter()
                .find(|(tag, _)| tag == &variant.tag)
                .map(|(_, payload)| payload)
                .ok_or_else(|| {
                    OpError::Domain("numerical variant result has an undeclared tag".into())
                })?;
            if schema.len() != variant.payload.len() {
                return Err(OpError::Domain(
                    "numerical variant result payload has incorrect arity".into(),
                ));
            }
            variant.payload = variant
                .payload
                .into_iter()
                .zip(schema)
                .map(|(value, ty)| decorate(value, ty))
                .collect::<Result<_, _>>()?;
            Ok(Value::Variant(variant))
        },
        | (Value::Opaque(value), TypeExpr::Opaque(kind)) if &value.kind == kind => {
            Ok(Value::Opaque(value))
        },
        | (value, ty) if scalar_representation_matches(&value, ty) => Ok(value),
        | (_, _) => {
            Err(OpError::Domain(
                "numerical result does not match its declared schema".into(),
            ))
        },
    }
}
fn path_value<'a>(
    args: &'a [Value],
    arg: usize,
    path: &[&str],
) -> Result<&'a Value, OpError> {
    let mut value = args
        .get(arg)
        .ok_or_else(|| OpError::Domain(format!("missing contract argument {arg}")))?;
    for name in path {
        let Value::Record(record) = value else {
            return Err(OpError::Domain(format!(
                "contract path {name} requires a record"
            )));
        };
        let mut matching = record.fields.iter().filter(|(field, _)| field == name);
        value = &matching
            .next()
            .ok_or_else(|| OpError::Domain(format!("missing contract field {name}")))?
            .1;
        if matching.next().is_some() {
            return Err(OpError::Domain(format!("duplicate contract field {name}")));
        }
    }
    Ok(value)
}
fn numeric_element(ty: &TypeExpr) -> bool {
    match ty {
        | TypeExpr::Dim(inner, _) | TypeExpr::Tensor { elem: inner, .. } => numeric_element(inner),
        | TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex => true,
        | _ => false,
    }
}
/// Bind and check the exact input paths described by a dimensional contract.
pub(super) fn prepare_contract(
    contract: &super::dimensions::DimContract,
    args: &[Value],
) -> Result<BTreeMap<&'static str, Dimension>, OpError> {
    let mut vars = BTreeMap::new();
    for input in contract.inputs {
        let value = path_value(args, input.argument, input.path)?;
        if !matches!(
            value,
            Value::Tensor(_)
                | Value::Int(_)
                | Value::Float(_)
                | Value::Real(_)
                | Value::Complex(_)
                | Value::Dyn(_)
                | Value::Quantity(_)
        ) {
            return Err(OpError::Domain(
                "contract dimension requires a numeric value".into(),
            ));
        }
        if let Value::Tensor(tensor) = value {
            let elem = match &tensor.elem {
                | TypeExpr::Dim(inner, _) => inner.as_ref(),
                | other => other,
            };
            if !numeric_element(elem) {
                return Err(OpError::Domain(
                    "contract dimension requires a numeric tensor element".into(),
                ));
            }
        }
        let dim = dimension(value)?;
        if vars.insert(input.name, dim).is_some() {
            return Err(OpError::Internal(format!(
                "duplicate contract variable {}",
                input.name
            )));
        }
    }
    if contract.other_args_dimensionless {
        for (index, value) in args.iter().enumerate() {
            if !contract.inputs.iter().any(|input| input.argument == index)
                && !contract
                    .constraints
                    .iter()
                    .any(|check| check.argument == index)
                && !dimension(value)?.is_empty()
            {
                return Err(OpError::Domain(format!(
                    "argument {index} must be dimensionless"
                )));
            }
        }
    }
    let expressions = vars
        .iter()
        .map(|(name, dim)| Ok((*name, expression(dim)?)))
        .collect::<Result<Vec<_>, OpError>>()?;
    let refs = expressions
        .iter()
        .map(|(name, expr)| (*name, expr))
        .collect::<Vec<_>>();
    for check in contract.constraints {
        let value = path_value(args, check.argument, check.path)?;
        if matches!(value, Value::Variant(v) if v.tag == "None" && v.payload.is_empty()) {
            continue;
        }
        let actual = dimension(value)?;
        let expected = canonical(&check.formula.instantiate(&refs))?;
        if actual != expected {
            return Err(OpError::Domain(format!(
                "argument {} has incompatible contract dimension",
                check.argument
            )));
        }
    }
    Ok(vars)
}
/// Validate the input dimensions before running a numerical kernel.
pub(super) fn prepare(
    policy: DimensionPolicy,
    args: &[Value],
) -> Result<BTreeMap<&'static str, Dimension>, OpError> {
    prepare_contract(super::dimensions::contract_for(policy), args)
}
/// Resolve explicitly declared shape-dependent output powers using validated input tensors.
pub(super) fn finish_shaped(
    result_type: &TypeExpr,
    vars: &BTreeMap<&str, Dimension>,
    args: &[Value],
    result: Value,
) -> Result<Value, OpError> {
    fn resolve_dim(
        dim: &mut DimExpr,
        args: &[Value],
    ) -> Result<(), OpError> {
        match dim {
            | DimExpr::Pow(base, exponent) => {
                resolve_dim(base, args)?;
                if let Exp::Symbol(symbol) = exponent
                    && let Some(source) = symbol.strip_prefix("shape:")
                {
                    let parts = source
                        .split(':')
                        .map(str::parse::<i64>)
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|_| {
                            OpError::Internal("invalid shape dimension expression".into())
                        })?;
                    if parts.len() != 3 {
                        return Err(OpError::Internal("invalid shape power arity".into()));
                    }
                    let arg = usize::try_from(parts[0])
                        .map_err(|_| OpError::Internal("invalid shape argument".into()))?;
                    let axis = usize::try_from(parts[1])
                        .map_err(|_| OpError::Internal("invalid shape axis".into()))?;
                    let Some(Value::Tensor(tensor)) = args.get(arg) else {
                        return Err(OpError::Domain(
                            "shape-dependent units require a tensor argument".into(),
                        ));
                    };
                    let extent = tensor.shape.get(axis).copied().ok_or_else(|| {
                        OpError::Domain("shape dimension axis is out of range".into())
                    })?;
                    let power = i64::try_from(extent)
                        .ok()
                        .and_then(|n| n.checked_mul(parts[2]))
                        .ok_or_else(|| {
                            OpError::Domain("shape dimension exponent overflows".into())
                        })?;
                    *exponent = Exp::Rat(Rational::new(power, 1));
                }
            },
            | DimExpr::Mul(parts) => {
                for part in parts {
                    resolve_dim(part, args)?;
                }
            },
            | DimExpr::Inverse(base) => resolve_dim(base, args)?,
            | _ => {},
        }
        Ok(())
    }
    fn resolve_type(
        ty: &mut TypeExpr,
        args: &[Value],
    ) -> Result<(), OpError> {
        match ty {
            | TypeExpr::Dim(inner, dim) => {
                resolve_dim(dim, args)?;
                resolve_type(inner, args)?;
            },
            | TypeExpr::Tensor { elem, .. } | TypeExpr::Sparse { elem, .. } => {
                resolve_type(elem, args)?;
            },
            | TypeExpr::Record(fields) => {
                for (_, ty) in fields {
                    resolve_type(ty, args)?;
                }
            },
            | TypeExpr::Variant(variants) => {
                for (_, payload) in variants {
                    for ty in payload {
                        resolve_type(ty, args)?;
                    }
                }
            },
            | _ => {},
        }
        Ok(())
    }
    let mut result_type = result_type.clone();
    resolve_type(&mut result_type, args)?;
    decorate(result, &instantiated_type(&result_type, vars)?)
}
/// Attach the instantiated dimensional witnesses to a kernel result.
pub(super) fn finish(
    sig: &crate::kernel::signature::Signature,
    vars: &BTreeMap<&str, Dimension>,
    result: Value,
) -> Result<Value, OpError> {
    decorate(result, &instantiated_type(&sig.ret, vars)?)
}

#[cfg(test)]
mod contract_tests {
    use super::*;
    use crate::kernel::value::RecordValue;
    #[test]
    fn malformed_nested_result_record_is_not_decorated() {
        let schema = TypeExpr::Record(vec![(
            "outer".into(),
            TypeExpr::Record(vec![
                ("a".into(), TypeExpr::Float),
                ("b".into(), TypeExpr::Float),
            ]),
        )]);
        let malformed = Value::Record(RecordValue {
            fields: vec![(
                "outer".into(),
                Value::Record(RecordValue {
                    fields: vec![
                        ("a".into(), Value::Float(1.0)),
                        ("a".into(), Value::Float(2.0)),
                    ],
                }),
            )],
        });
        assert!(decorate(malformed, &schema).is_err());
    }
    #[test]
    fn malformed_scalars_and_tensor_shapes_fail_closed() {
        use crate::kernel::signature::NatExpr;
        use crate::kernel::value::TensorValue;
        let scalar = TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base("Length".into()));
        assert!(decorate(Value::Bool(true), &scalar).is_err());
        assert!(decorate(Value::Bool(true), &TypeExpr::Float).is_err());
        let schema = TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Lit(2)]),
            elem: Box::new(TypeExpr::Float),
        };
        let value = |shape, data, elem| Value::Tensor(TensorValue { shape, data, elem });
        assert!(
            decorate(
                value(
                    vec![1, 2],
                    vec![Value::Float(1.0), Value::Float(2.0)],
                    TypeExpr::Float
                ),
                &schema
            )
            .is_err()
        );
        assert!(
            decorate(
                value(vec![2], vec![Value::Float(1.0)], TypeExpr::Float),
                &schema
            )
            .is_err()
        );
        assert!(
            decorate(
                value(
                    vec![2],
                    vec![Value::Float(1.0), Value::Float(2.0)],
                    TypeExpr::Bool
                ),
                &schema
            )
            .is_err()
        );
        let empty = TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Lit(0)]),
            elem: Box::new(TypeExpr::Float),
        };
        assert!(decorate(value(vec![0], vec![], TypeExpr::Bool), &empty).is_err());
    }
}

#[cfg(test)]
mod variant_contract_tests {
    use super::*;
    use crate::kernel::value::VariantValue;

    fn tagged(
        tag: &str,
        payload: Vec<Value>,
    ) -> Value {
        Value::Variant(VariantValue {
            tag: tag.into(),
            payload,
        })
    }

    #[test]
    fn option_payloads_receive_instantiated_dimensions_without_accepting_malformed_tags() {
        let schema = TypeExpr::Variant(vec![
            ("None".into(), vec![]),
            (
                "Some".into(),
                vec![TypeExpr::Float.with_dim(DimExpr::Var("D".into()))],
            ),
        ]);
        let length = canonical(&DimExpr::Base("Length".into())).unwrap();
        let schema = instantiated_type(&schema, &BTreeMap::from([("D", length)])).unwrap();
        assert!(decorate(tagged("None", vec![]), &schema).is_ok());
        let value = decorate(tagged("Some", vec![Value::Float(2.0)]), &schema).unwrap();
        let Value::Variant(value) = value else {
            panic!("expected variant")
        };
        let Value::Dyn(value) = &value.payload[0] else {
            panic!("expected dimensional payload")
        };
        assert_eq!(value.dim, DimExpr::Base("Length".into()));
        for invalid in [
            tagged("Some", vec![]),
            tagged("None", vec![Value::Float(2.0)]),
            tagged("Unknown", vec![]),
            tagged("Some", vec![Value::Bool(true)]),
            tagged(
                "Some",
                vec![Value::Dyn(DynValue {
                    value: Box::new(Value::Float(2.0)),
                    dim: DimExpr::Base("Time".into()),
                })],
            ),
        ] {
            assert!(decorate(invalid, &schema).is_err());
        }
        let malformed = TypeExpr::Variant(vec![("None".into(), vec![]), ("None".into(), vec![])]);
        assert!(decorate(tagged("None", vec![]), &malformed).is_err());
    }
}
