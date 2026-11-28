//! Explicit dimensional, budgeted adapters for adaptive stochastic paths.
//! Lévy noise is dimensionless at an explicit reference-time scale, avoiding
//! dimension exponents inferred from an arbitrary floating-point alpha.
use super::CallCtx;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpError;
use super::RecordValue;
use super::RefCell;
use super::RuleError;
use super::Signature;
use super::ToPrimitive;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::canonical;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::effectful;
use super::common::positive_budget;
use super::common::register;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::state_type;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::expression;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::Exp;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Rational;
use crate::kernel::value::TensorValue;
use crate::sim::stochastic::JumpDistribution;
use crate::sim::stochastic::JumpSpec;
use crate::sim::stochastic::Method;
use crate::sim::stochastic::Options;
use crate::sim::stochastic::SdeError;
use crate::sim::stochastic::Structure;
use crate::sim::stochastic::{
    self,
};
use std::cell::Cell;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Family {
    Brownian,
    Levy,
    Jump,
}
#[derive(Clone, Copy)]
struct Adapter {
    family: Family,
    general: bool,
}
impl Adapter {
    const fn base(self) -> usize {
        if self.general { 13 } else { 11 }
    }

    const fn callback_start(self) -> usize {
        self.base()
            + if matches!(self.family, Family::Brownian) {
                0
            } else {
                2
            }
    }

    const fn arity(self) -> usize {
        self.callback_start()
            + if matches!(self.family, Family::Jump) {
                3
            } else {
                2
            }
    }
}
// The mark extent is independently anchored by the constant mark argument;
// it must never be conflated with the state extent N or noise channel count.
fn mark_type() -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: Some(vec![NatExpr::Var("M".into())]),
        elem: Box::new(TypeExpr::Float),
    }
}
fn dimension_product(
    a: DimExpr,
    b: DimExpr,
) -> Result<DimExpr, OpError> {
    expression(&canonical(&DimExpr::Mul(vec![a, b]))?)
}
fn check_coefficients(
    ty: &TypeExpr,
    state: usize,
    noise: usize,
    general: bool,
) -> Result<DimExpr, OpError> {
    if !general {
        return vector_dimension(ty, Some(state));
    }
    let TypeExpr::Tensor {
        rank: NatExpr::Lit(2),
        shape,
        elem,
    } = ty
    else {
        return Err(OpError::Domain(
            "general diffusion must return a rank-two real tensor".into(),
        ));
    };
    if let Some(shape) = shape
        && !matches!(shape.as_slice(),[NatExpr::Lit(n),NatExpr::Lit(m)] if *n==state as u64&&*m==noise as u64)
    {
        return Err(OpError::Domain(
            "general diffusion schema must have concrete state x noise extents".into(),
        ));
    }
    scalar_dimension(elem)
}
fn coefficient_values(
    value: &Value,
    dimension: &DimExpr,
    state: usize,
    noise: usize,
    general: bool,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<Vec<f64>, OpError> {
    if !general {
        return vector(value, dimension, Some(state), index, ctx);
    }
    let Value::Tensor(t) = value else {
        return Err(OpError::Domain(
            "general diffusion returned a non-tensor".into(),
        ));
    };
    if t.shape != [state, noise] || state.checked_mul(noise) != Some(t.data.len()) {
        return Err(OpError::Domain(
            "general diffusion returned inconsistent state x noise storage".into(),
        ));
    }
    if canonical(&check_coefficients(
        &ctx.value_type(value)?,
        state,
        noise,
        true,
    )?)? != canonical(dimension)?
    {
        return Err(OpError::Domain(
            "diffusion coefficient dimensions disagree with their contract".into(),
        ));
    }
    t.data
        .iter()
        .map(|value| {
            if matches!(value, Value::Dyn(_) | Value::Quantity(_)) {
                check_dimension(ctx, value, dimension, index)?;
            }
            numeric(value, index)
        })
        .collect()
}
impl EvalDyn for Adapter {
    fn signature(&self) -> Signature {
        let d = DimExpr::Var("D".into());
        let time = DimExpr::Var("T".into());
        let rate = DimExpr::Mul(vec![d.clone(), DimExpr::Inverse(Box::new(time.clone()))]);
        let diffusion = if self.family == Family::Levy {
            d.clone()
        } else {
            DimExpr::Mul(vec![
                d.clone(),
                time.clone().pow(Exp::Rat(Rational::new(-1, 2))),
            ])
        };
        let coefficient = if self.general {
            TypeExpr::Tensor {
                rank: NatExpr::Lit(2),
                shape: None,
                elem: Box::new(dimensional_type(diffusion)),
            }
        } else {
            state_type(diffusion)
        };
        let mut params = vec![
            dimensional_type(time.clone()),
            dimensional_type(time.clone()),
            dimensional_type(time.clone()),
            state_type(d.clone()),
            dimensional_type(d.clone()),
            TypeExpr::Float,
            TypeExpr::Int,
            TypeExpr::Int,
            TypeExpr::Int,
            TypeExpr::Int,
            TypeExpr::Int,
        ];
        if self.general {
            params.extend([TypeExpr::Int, TypeExpr::Bool]);
        }
        match self.family {
            | Family::Brownian => {},
            | Family::Levy => params.extend([TypeExpr::Float, TypeExpr::Float]),
            | Family::Jump => {
                params.extend([
                    dimensional_type(DimExpr::Inverse(Box::new(time.clone()))),
                    mark_type(),
                ]);
            },
        }
        params.extend([
            TypeExpr::Function(
                vec![dimensional_type(time.clone()), state_type(d.clone())],
                Box::new(state_type(rate)),
            ),
            TypeExpr::Function(
                vec![dimensional_type(time.clone()), state_type(d.clone())],
                Box::new(coefficient),
            ),
        ]);
        if self.family == Family::Jump {
            params.push(TypeExpr::Function(
                vec![
                    dimensional_type(time.clone()),
                    state_type(d.clone()),
                    mark_type(),
                ],
                Box::new(state_type(d.clone())),
            ));
        }
        let mut fields = vec![
            ("times".into(), tensor_type(1, time)),
            ("states".into(), tensor_type(2, d)),
            ("evaluations".into(), TypeExpr::Int),
            ("steps".into(), TypeExpr::Int),
            ("rejected".into(), TypeExpr::Int),
        ];
        if self.family == Family::Jump {
            fields.extend([
                ("total_jumps".into(), TypeExpr::Int),
                (
                    "jump_counts".into(),
                    TypeExpr::Tensor {
                        rank: NatExpr::Lit(1),
                        shape: None,
                        elem: Box::new(TypeExpr::Int),
                    },
                ),
            ]);
        }
        let mut signature = Signature::new(params, TypeExpr::Record(fields));
        signature.type_params = vec![
            TypeParam {
                name: "D".into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            },
            TypeParam {
                name: "T".into(),
                kind: TypeParamKind::Dim,
                bounds: vec![],
            },
            TypeParam {
                name: "N".into(),
                kind: TypeParamKind::Nat,
                bounds: vec![],
            },
        ];
        if self.family == Family::Jump {
            signature.type_params.push(TypeParam {
                name: "M".into(),
                kind: TypeParamKind::Nat,
                bounds: vec![],
            });
        }
        effectful(signature)
    }

    #[allow(clippy::too_many_lines)] // preflight, shared callback lifetime and final error precedence are one audited transaction
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        if args.len() != self.arity() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(self.arity() as u8),
                got: args.len(),
            });
        }
        let at = self.callback_start();
        let schema = self.signature();
        let (drift, drift_sig) =
            CallbackRun::new(ctx, &args[at], at, 2, schema.params[at].clone())?;
        let (diffusion, diff_sig) =
            CallbackRun::new(ctx, &args[at + 1], at + 1, 2, schema.params[at + 1].clone())?;
        let time_dimension = scalar_dimension(&drift_sig.params[0])?;
        let state_dimension = vector_dimension(&drift_sig.params[1], None)?;
        let initial = vector(&args[3], &state_dimension, None, 3, ctx)?;
        if initial.is_empty() {
            return Err(OpError::Domain("SDE initial state must be nonempty".into()));
        }
        let n = initial.len();
        let noise = if self.general {
            positive_budget(&args[11], 11, "noise channels")?
        } else {
            n
        };
        let commutative = if self.general {
            let Value::Bool(value) = args[12] else {
                return Err(OpError::Domain("commutative must be Bool".into()));
            };
            value
        } else {
            true
        };
        if n.checked_mul(noise).is_none_or(|v| v > 8_388_608) {
            return Err(OpError::Domain(
                "SDE coefficient storage exceeds the finite work budget".into(),
            ));
        }
        let rate_dimension = dimension_product(
            state_dimension.clone(),
            DimExpr::Inverse(Box::new(time_dimension.clone())),
        )?;
        let diffusion_dimension = if self.family == Family::Levy {
            state_dimension.clone()
        } else {
            dimension_product(
                state_dimension.clone(),
                time_dimension.clone().pow(Exp::Rat(Rational::new(-1, 2))),
            )?
        };
        vector_dimension(&drift_sig.params[1], Some(n))?;
        if canonical(&vector_dimension(&drift_sig.ret, Some(n))?)? != canonical(&rate_dimension)?
            || canonical(&scalar_dimension(&diff_sig.params[0])?)? != canonical(&time_dimension)?
            || canonical(&vector_dimension(&diff_sig.params[1], Some(n))?)?
                != canonical(&state_dimension)?
            || canonical(&check_coefficients(&diff_sig.ret, n, noise, self.general)?)?
                != canonical(&diffusion_dimension)?
        {
            return Err(OpError::Domain(
                "SDE drift/diffusion callback shape or dimensions violate the explicit contract"
                    .into(),
            ));
        }
        for i in [0usize, 1, 2] {
            check_dimension(ctx, &args[i], &time_dimension, i)?;
        }
        check_dimension(ctx, &args[4], &state_dimension, 4)?;
        check_dimension(ctx, &args[5], &DimExpr::Dimless, 5)?;
        let t0 = numeric(&args[0], 0)?;
        let t1 = numeric(&args[1], 1)?;
        let reference = numeric(&args[2], 2)?;
        if reference <= 0.0 || t1 <= t0 {
            return Err(OpError::Domain(
                "SDE requires positive reference time and increasing endpoints".into(),
            ));
        }
        let atol = numeric(&args[4], 4)?;
        let rtol = numeric(&args[5], 5)?;
        let cells = positive_budget(&args[6], 6, "finest noise cells")?;
        let steps = positive_budget(&args[7], 7, "proposal/accepted-step budget")?;
        let budget = positive_budget(&args[8], 8, "callback budget")?;
        let Value::Int(seed) = &args[9] else {
            return Err(OpError::Domain(
                "SDE seed must be an exact nonnegative Int".into(),
            ));
        };
        let seed = seed
            .to_u64()
            .ok_or_else(|| OpError::Domain("SDE seed must fit u64".into()))?;
        let Value::Int(method) = &args[10] else {
            return Err(OpError::Domain("SDE method must be an exact Int".into()));
        };
        let method = match method.to_u32() {
            | Some(0) => Method::EulerMaruyama,
            | Some(1) => Method::Milstein,
            | Some(2) => Method::MilsteinStratonovich,
            | Some(3) => Method::HeunStratonovich,
            | Some(4) => Method::SrkPlaten,
            | Some(5) => Method::MilsteinNoncommutative,
            | Some(6) => Method::MilsteinNoncommutativeStratonovich,
            | Some(7) => Method::DriftImplicitEuler,
            | Some(8) => Method::ThetaMaruyama,
            | Some(9) => Method::DriftImplicitMilstein,
            | _ => return Err(OpError::Domain("SDE method must be in 0..=9".into())),
        };
        let structure = if self.general {
            Structure::General { noise, commutative }
        } else {
            Structure::Diagonal
        };
        let options = Options::adaptive(
            method,
            structure,
            seed,
            0.0,
            (t1 - t0) / reference,
            rtol,
            atol,
            steps,
            cells,
        );
        let jump = if self.family == Family::Jump {
            Some(CallbackRun::new(
                ctx,
                &args[at + 2],
                at + 2,
                3,
                schema.params[at + 2].clone(),
            )?)
        } else {
            None
        };
        let mark = if self.family == Family::Jump {
            vector(
                &args[self.base() + 1],
                &DimExpr::Dimless,
                None,
                self.base() + 1,
                ctx,
            )?
        } else {
            vec![]
        };
        if let Some((_, sig)) = &jump
            && (canonical(&scalar_dimension(&sig.params[0])?)? != canonical(&time_dimension)?
                || canonical(&vector_dimension(&sig.params[1], Some(n))?)?
                    != canonical(&state_dimension)?
                || canonical(&vector_dimension(&sig.params[2], Some(mark.len()))?)?
                    != canonical(&DimExpr::Dimless)?
                || canonical(&vector_dimension(&sig.ret, Some(n))?)?
                    != canonical(&state_dimension)?)
        {
            return Err(OpError::Domain("jump callback must map (physical time,state,dimensionless mark) to a state increment".into()));
        }
        let used = Cell::new(0usize);
        let budget_error = RefCell::new(None);
        let allowed = || {
            if drift.failed()
                || diffusion.failed()
                || jump.as_ref().is_some_and(|(run, _)| run.failed())
                || budget_error.borrow().is_some()
            {
                return false;
            }
            if used.get() >= budget {
                *budget_error.borrow_mut() = Some(OpError::BudgetExceeded);
                return false;
            }
            used.set(used.get() + 1);
            true
        };
        let state_value = |y: &[f64]| {
            Value::Tensor(TensorValue {
                shape: vec![n],
                data: y.iter().copied().map(Value::Float).collect(),
                elem: dimensional_type(state_dimension.clone()),
            })
        };
        let a = |u: f64, y: &[f64]| {
            if !allowed() {
                return Err(SdeError::Invalid("callback budget/failure".into()));
            }
            drift
                .evaluate(
                    &[
                        dimensional_float(t0 + reference * u, time_dimension.clone()),
                        state_value(y),
                    ],
                    |value| {
                        let mut values = vector(&value, &rate_dimension, Some(n), at, ctx)?;
                        for value in &mut values {
                            *value *= reference;
                        }
                        if values.iter().any(|v| !v.is_finite()) {
                            return Err(OpError::Domain(
                                "reference time overflows the normalized SDE drift".into(),
                            ));
                        }
                        Ok(values)
                    },
                )
                .ok_or_else(|| SdeError::Invalid("drift callback failed".into()))
        };
        let g = |u: f64, y: &[f64]| {
            if !allowed() {
                return Err(SdeError::Invalid("callback budget/failure".into()));
            }
            diffusion
                .evaluate(
                    &[
                        dimensional_float(t0 + reference * u, time_dimension.clone()),
                        state_value(y),
                    ],
                    |value| {
                        let mut values = coefficient_values(
                            &value,
                            &diffusion_dimension,
                            n,
                            noise,
                            self.general,
                            at + 1,
                            ctx,
                        )?;
                        if self.family != Family::Levy {
                            for value in &mut values {
                                *value *= reference.sqrt();
                            }
                        }
                        if values.iter().any(|v| !v.is_finite()) {
                            return Err(OpError::Domain(
                                "reference time overflows the normalized SDE diffusion".into(),
                            ));
                        }
                        Ok(values)
                    },
                )
                .ok_or_else(|| SdeError::Invalid("diffusion callback failed".into()))
        };
        let j = |u: f64, y: &[f64], mark: &[f64]| {
            if !allowed() {
                return Err(SdeError::Invalid("callback budget/failure".into()));
            }
            let Some((jump, _)) = &jump else {
                return Err(SdeError::Invalid("missing jump callback".into()));
            };
            jump.evaluate(
                &[
                    dimensional_float(t0 + reference * u, time_dimension.clone()),
                    state_value(y),
                    Value::Tensor(TensorValue {
                        shape: vec![mark.len()],
                        data: mark.iter().copied().map(Value::Float).collect(),
                        elem: TypeExpr::Float,
                    }),
                ],
                |value| vector(&value, &state_dimension, Some(n), at + 2, ctx),
            )
            .ok_or_else(|| SdeError::Invalid("jump callback failed".into()))
        };
        let mut jump_info = None;
        let outcome = match self.family {
            | Family::Brownian => stochastic::integrate(a, g, &initial, &options),
            | Family::Levy => {
                for i in [self.base(), self.base() + 1] {
                    check_dimension(ctx, &args[i], &DimExpr::Dimless, i)?;
                }
                stochastic::integrate_levy(
                    a,
                    g,
                    &initial,
                    numeric(&args[self.base()], self.base())?,
                    numeric(&args[self.base() + 1], self.base() + 1)?,
                    &options,
                )
            },
            | Family::Jump => {
                let inverse = DimExpr::Inverse(Box::new(time_dimension.clone()));
                check_dimension(ctx, &args[self.base()], &inverse, self.base())?;
                let spec = JumpSpec {
                    intensity: numeric(&args[self.base()], self.base())? * reference,
                    distribution: JumpDistribution::Constant(mark),
                    compensated: false,
                };
                stochastic::integrate_jump_diffusion(a, g, j, &initial, &spec, &options).map(
                    |result| {
                        jump_info = Some((result.total_jumps, result.jump_counts));
                        result.path
                    },
                )
            },
        };
        drift.finish()?;
        diffusion.finish()?;
        if let Some((jump, _)) = jump {
            jump.finish()?;
        }
        if let Some(error) = budget_error.into_inner() {
            return Err(error);
        }
        let mut path = outcome.map_err(|error| {
            match error {
                | SdeError::StepLimit(message) => OpError::NotReduced(message),
                | other => OpError::Domain(other.to_string()),
            }
        })?;
        for time in &mut path.t {
            *time = t0 + reference * *time;
        }
        if let Some(last) = path.t.last_mut() {
            *last = t1;
        }
        let rows = path.t.len();
        let mut fields = vec![
            (
                "times".into(),
                output_tensor(vec![rows], path.t, time_dimension)?,
            ),
            (
                "states".into(),
                output_tensor(
                    vec![rows, n],
                    path.y.into_iter().flatten().collect(),
                    state_dimension,
                )?,
            ),
            ("evaluations".into(), Value::int(used.get())),
            ("steps".into(), Value::int(path.steps)),
            ("rejected".into(), Value::int(path.rejected)),
        ];
        if let Some((total, counts)) = jump_info {
            fields.extend([
                ("total_jumps".into(), Value::int(total)),
                (
                    "jump_counts".into(),
                    Value::Tensor(TensorValue {
                        shape: vec![counts.len()],
                        data: counts.into_iter().map(Value::int).collect(),
                        elem: TypeExpr::Int,
                    }),
                ),
            ]);
        }
        Ok(Value::Record(RecordValue { fields }))
    }
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (family, stem, path) in [
        (
            Family::Brownian,
            "num_sde_adaptive",
            "crate::sim::stochastic::integrate",
        ),
        (
            Family::Levy,
            "num_sde_levy_adaptive",
            "crate::sim::stochastic::integrate_levy",
        ),
        (
            Family::Jump,
            "num_sde_jump_adaptive",
            "crate::sim::stochastic::integrate_jump_diffusion",
        ),
    ] {
        for general in [false, true] {
            let evaluator = Adapter { family, general };
            let name = format!("{stem}_{}", if general { "general" } else { "diag" });
            register(
                installer,
                &name,
                path,
                "Adaptive stochastic path with exact declared callback dimensions and finite budgets. Arguments: t0,t1,positive reference_time,y0,state atol,dimensionless rtol,positive finest cells,proposal/path step budget,callback budget,nonnegative u64 seed,method code 0..9; general variants then take noise count and commutative Bool. Levy variants next take alpha,beta and use dimensionless stable noise at the reference-time scale (diffusion has state units). Jump variants next take intensity (inverse time),dimensionless constant mark vector with independent extent M, then drift,diffusion,jump callbacks sharing that M; other families finish with drift,diffusion callbacks. Brownian diffusion has state/sqrt(time) units. Method codes follow SdeMethod declaration order. Rejections reuse a bounded dyadic noise realization; resolution or work exhaustion is an error. Stable alpha<2 supports only Ito Euler/implicit Euler/theta, while alpha=2 uses the correctly scaled Brownian schemes. Jump times/marks are fixed before rejection and continuous segments stop exactly at each arrival. Local error estimates are heuristic, not global certificates.",
                evaluator.arity() as u8,
                evaluator,
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod mark_extent_tests {
    use super::Adapter;
    use super::EvalDyn;
    use super::Family;
    use super::NatExpr;
    use super::TypeExpr;
    use super::TypeParamKind;

    #[test]
    fn jump_mark_extent_is_anchored_independently_of_state() {
        for general in [false, true] {
            let adapter = Adapter {
                family: Family::Jump,
                general,
            };
            let signature = adapter.signature();
            let mark = &signature.params[adapter.base() + 1];
            let TypeExpr::Tensor { shape, .. } = mark else {
                panic!("mark must be a tensor");
            };
            assert_eq!(shape, &Some(vec![NatExpr::Var("M".into())]));
            let TypeExpr::Function(params, _) = &signature.params[adapter.callback_start() + 2]
            else {
                panic!("jump must be a callback");
            };
            assert_eq!(&params[2], mark);
            let TypeExpr::Tensor { shape, .. } = &params[1] else {
                panic!("state must be a tensor");
            };
            assert_eq!(shape, &Some(vec![NatExpr::Var("N".into())]));
            assert!(
                signature
                    .type_params
                    .iter()
                    .any(|p| p.name == "M" && p.kind == TypeParamKind::Nat)
            );
        }
    }
}
