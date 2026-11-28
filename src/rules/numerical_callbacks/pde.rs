//! Dimensioned adaptive PDE exports, with explicit coordinate/time/field contracts.
use super::Arc;
use super::Arity;
use super::BackendSupport;
use super::CallCtx;
use super::Determinism;
use super::DimExpr;
use super::EvalDyn;
use super::Installer;
use super::OpDescriptor;
use super::OpError;
use super::Purity;
use super::RecordValue;
use super::RefCell;
use super::RuleError;
use super::Signature;
use super::Stability;
use super::ThreadSafety;
use super::ToPrimitive;
use super::TypeExpr;
use super::TypeParam;
use super::TypeParamKind;
use super::Value;
use super::canonical;
use super::check_dimension;
use super::common::CallbackRun;
use super::common::positive_budget;
use super::dimensional_float;
use super::dimensional_type;
use super::evolution::output_tensor;
use super::evolution::tensor_type;
use super::evolution::vector;
use super::evolution::vector_dimension;
use super::numeric;
use super::scalar_dimension;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::NatExpr;
use crate::kernel::value::VariantValue;
use crate::sim::fdm_nd::Boundary3D;
use crate::sim::fdm_nd::BoundaryCondition;
use crate::sim::fdm_nd::FdmError;
use crate::sim::fdm_nd::TimeStepping;
use crate::sim::fdm_nd::adaptive::PdeAdaptiveOptions;
use crate::sim::fdm_nd::adaptive::PdeAdaptiveSolution;
use crate::sim::fdm_nd::adaptive::PdeTerminalEvent;
use crate::sim::fdm_nd::adaptive::{
    self,
};
use crate::sim::pde_api::FvmScheme;
use crate::sim::pde_api::{
    self,
};
use std::cell::Cell;

#[derive(Clone, Copy)]
enum Scheme {
    Heat(TimeStepping),
    Wave,
    Advection,
    Fem(f64),
    FemSource(f64),
    Fvm(FvmScheme),
    Dg,
    DgNd,
    Spectral,
}
#[derive(Clone, Copy)]
struct Pde {
    scheme: Scheme,
    event: bool,
}
fn ratio(
    a: DimExpr,
    b: DimExpr,
) -> DimExpr {
    DimExpr::Mul(vec![a, DimExpr::Inverse(Box::new(b))])
}
fn record(value: &Value) -> Result<&RecordValue, OpError> {
    if let Value::Record(r) = value {
        Ok(r)
    } else {
        Err(OpError::Domain(
            "PDE controls and boundary data must be records".into(),
        ))
    }
}
fn field<'a>(
    record: &'a RecordValue,
    name: &str,
) -> Result<&'a Value, OpError> {
    record
        .get(name)
        .ok_or_else(|| OpError::Domain(format!("PDE record is missing {name}")))
}
fn same(
    a: &DimExpr,
    b: &DimExpr,
) -> Result<(), OpError> {
    if canonical(a)? == canonical(b)? {
        Ok(())
    } else {
        Err(OpError::Domain(
            "PDE coordinate, time, field, coefficient, or callback dimensions disagree".into(),
        ))
    }
}
fn scalar(
    value: &Value,
    dim: &DimExpr,
    index: usize,
    ctx: &CallCtx<'_>,
) -> Result<f64, OpError> {
    check_dimension(ctx, value, dim, index)?;
    numeric(value, index)
}
fn ints(
    value: &Value,
    index: usize,
) -> Result<Vec<usize>, OpError> {
    let Value::Tensor(t) = value else {
        return Err(OpError::Domain(
            "PDE grid extents and boundary kinds must be integer vectors".into(),
        ));
    };
    if t.elem != TypeExpr::Int {
        return Err(OpError::Domain(
            "PDE grid extents and boundary kinds must be dimensionless integers".into(),
        ));
    }
    if t.shape.len() != 1 || t.shape[0] != t.data.len() {
        return Err(OpError::Domain("malformed PDE integer vector".into()));
    }
    t.data
        .iter()
        .map(|v| {
            let Value::Int(n) = v else {
                return Err(OpError::Domain(
                    "PDE integer vector contains a noninteger".into(),
                ));
            };
            n.to_usize().ok_or_else(|| {
                OpError::Domain(format!("argument {index} contains an invalid integer"))
            })
        })
        .collect()
}
fn integer_vector() -> TypeExpr {
    TypeExpr::Tensor {
        rank: NatExpr::Lit(1),
        shape: None,
        elem: Box::new(TypeExpr::Int),
    }
}
fn controls_type(
    t: &DimExpr,
    u: &DimExpr,
) -> TypeExpr {
    TypeExpr::Record(vec![
        ("t0".into(), dimensional_type(t.clone())),
        ("t1".into(), dimensional_type(t.clone())),
        ("rtol".into(), TypeExpr::Float),
        ("atol".into(), dimensional_type(u.clone())),
        ("h0".into(), dimensional_type(t.clone())),
        ("h_min".into(), dimensional_type(t.clone())),
        ("h_max".into(), dimensional_type(t.clone())),
        ("max_steps".into(), TypeExpr::Int),
        ("max_rejections".into(), TypeExpr::Int),
        ("max_step_calls".into(), TypeExpr::Int),
        ("max_evaluations".into(), TypeExpr::Int),
    ])
}
fn boundary_type(
    l: &DimExpr,
    u: &DimExpr,
) -> TypeExpr {
    TypeExpr::Record(vec![
        ("kind".into(), integer_vector()),
        ("value".into(), tensor_type(1, u.clone())),
        ("flux".into(), tensor_type(1, ratio(u.clone(), l.clone()))),
        ("alpha".into(), tensor_type(1, DimExpr::Dimless)),
        ("beta".into(), tensor_type(1, l.clone())),
        ("gamma".into(), tensor_type(1, u.clone())),
    ])
}
fn controls(
    value: &Value,
    t: &DimExpr,
    u: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<(PdeAdaptiveOptions, usize), OpError> {
    let r = record(value)?;
    let n = |name| positive_budget(field(r, name)?, 0, name);
    Ok((
        PdeAdaptiveOptions {
            t0: scalar(field(r, "t0")?, t, 0, ctx)?,
            t1: scalar(field(r, "t1")?, t, 0, ctx)?,
            rtol: scalar(field(r, "rtol")?, &DimExpr::Dimless, 0, ctx)?,
            atol: scalar(field(r, "atol")?, u, 0, ctx)?,
            h0: scalar(field(r, "h0")?, t, 0, ctx)?,
            h_min: scalar(field(r, "h_min")?, t, 0, ctx)?,
            h_max: scalar(field(r, "h_max")?, t, 0, ctx)?,
            max_steps: n("max_steps")?,
            max_rejections: n("max_rejections")?,
            max_step_calls: n("max_step_calls")?,
        },
        n("max_evaluations")?,
    ))
}
fn boundary(
    value: &Value,
    dim: usize,
    l: &DimExpr,
    u: &DimExpr,
    ctx: &CallCtx<'_>,
) -> Result<Boundary3D, OpError> {
    let r = record(value)?;
    let kind = ints(field(r, "kind")?, 0)?;
    let values = vector(field(r, "value")?, u, None, 0, ctx)?;
    let flux = vector(
        field(r, "flux")?,
        &ratio(u.clone(), l.clone()),
        None,
        0,
        ctx,
    )?;
    let alpha = vector(field(r, "alpha")?, &DimExpr::Dimless, None, 0, ctx)?;
    let beta = vector(field(r, "beta")?, l, None, 0, ctx)?;
    let gamma = vector(field(r, "gamma")?, u, None, 0, ctx)?;
    if kind.len() != 1 && kind.len() != 2 * dim {
        return Err(OpError::Domain(
            "boundary vectors need one broadcast entry or two faces per active axis".into(),
        ));
    }
    if [&values, &flux, &alpha, &beta, &gamma]
        .iter()
        .any(|v| v.len() != kind.len())
    {
        return Err(OpError::Domain(
            "all boundary vectors must have equal lengths".into(),
        ));
    }
    let mut faces = Vec::new();
    for face in 0..6 {
        if face >= 2 * dim {
            faces.push(BoundaryCondition::Periodic);
            continue;
        }
        let i = if kind.len() == 1 { 0 } else { face };
        faces.push(match kind[i] {
            | 0 => BoundaryCondition::Dirichlet(values[i]),
            | 1 => BoundaryCondition::Neumann(flux[i]),
            | 2 => {
                BoundaryCondition::Robin {
                    alpha: alpha[i],
                    beta: beta[i],
                    gamma: gamma[i],
                }
            },
            | 3 => BoundaryCondition::Periodic,
            | _ => {
                return Err(OpError::Domain(
                    "boundary kind is 0 Dirichlet, 1 Neumann, 2 Robin, or 3 periodic".into(),
                ));
            },
        });
    }
    Ok(Boundary3D {
        x_min: faces[0].clone(),
        x_max: faces[1].clone(),
        y_min: faces[2].clone(),
        y_max: faces[3].clone(),
        z_min: faces[4].clone(),
        z_max: faces[5].clone(),
    })
}
fn kernel_error(error: FdmError) -> OpError {
    match error{FdmError::NotConverged=>OpError::NotReduced("adaptive PDE tolerance, linear solve, minimum step, or attempt/rejection budget exhausted".into()),other=>OpError::Domain(other.to_string())}
}
struct Budget {
    limit: usize,
    count: Cell<usize>,
    halted: Cell<bool>,
    error: RefCell<Option<OpError>>,
}
impl Budget {
    const fn new(limit: usize) -> Self {
        Self {
            limit,
            count: Cell::new(0),
            halted: Cell::new(false),
            error: RefCell::new(None),
        }
    }

    fn invoke<T>(
        &self,
        run: &CallbackRun<'_, '_>,
        args: &[Value],
        decode: impl FnOnce(Value) -> Result<T, OpError>,
    ) -> Option<T> {
        if self.halted.get() {
            return None;
        }
        if self.count.get() >= self.limit {
            self.halted.set(true);
            *self.error.borrow_mut() = Some(OpError::NotReduced(
                "PDE callback evaluation budget exhausted".into(),
            ));
            return None;
        }
        self.count.set(self.count.get() + 1);
        let value = run.evaluate(args, decode);
        if value.is_none() {
            self.halted.set(true);
        }
        value
    }

    fn finish(&self) -> Result<(), OpError> {
        self.error.borrow_mut().take().map_or(Ok(()), Err)
    }
}
fn callback_scalar(
    signature: &Signature,
    argument: &DimExpr,
    result: &DimExpr,
) -> Result<(), OpError> {
    same(&scalar_dimension(&signature.params[0])?, argument)?;
    same(&scalar_dimension(&signature.ret)?, result)
}
impl Pde {
    const fn initial_index(self) -> usize {
        match self.scheme {
            | Scheme::Fem(_) | Scheme::FemSource(_) => 6,
            | Scheme::Fvm(_) => 1,
            | Scheme::Dg | Scheme::DgNd | Scheme::Spectral => 4,
            | _ => 3,
        }
    }

    const fn control_index(self) -> usize {
        match self.scheme {
            | Scheme::Wave => 7,
            | Scheme::Fem(_) | Scheme::FemSource(_) => 8,
            | Scheme::Fvm(_) => 3,
            | _ => 5,
        }
    }

    const fn base_count(self) -> usize {
        match self.scheme {
            | Scheme::Wave => 8,
            | Scheme::Fem(_) => 9,
            | Scheme::FemSource(_) => 10,
            | Scheme::Spectral => 7,
            | _ => 6,
        }
    }
}
fn matrix_dimension(
    ty: &TypeExpr,
    rows: usize,
    cols: usize,
) -> Result<DimExpr, OpError> {
    let TypeExpr::Tensor {
        rank: NatExpr::Lit(2),
        shape,
        elem,
    } = ty
    else {
        return Err(OpError::Domain(
            "FEM source coordinate input must be a rank-two tensor".into(),
        ));
    };
    if let Some(shape) = shape
        && shape != &vec![NatExpr::Lit(rows as u64), NatExpr::Lit(cols as u64)]
    {
        return Err(OpError::Domain(
            "FEM source coordinate shape must match nodal coordinates".into(),
        ));
    }
    scalar_dimension(elem)
}
impl EvalDyn for Pde {
    fn signature(&self) -> Signature {
        let t = DimExpr::Var("T".into());
        let l = DimExpr::Var("L".into());
        let u = DimExpr::Var("U".into());
        let rate = ratio(u.clone(), t.clone());
        let speed = ratio(l.clone(), t.clone());
        let diffusivity = ratio(l.clone().pow(2_i64), t.clone());
        let control = controls_type(&t, &u);
        let bc = boundary_type(&l, &u);
        let state = tensor_type(1, u.clone());
        let mut args = match self.scheme {
            | Scheme::Heat(_) => {
                vec![
                    integer_vector(),
                    dimensional_type(l.clone()),
                    dimensional_type(diffusivity),
                    state.clone(),
                    bc,
                    control,
                ]
            },
            | Scheme::Wave => {
                vec![
                    integer_vector(),
                    dimensional_type(l.clone()),
                    dimensional_type(speed),
                    state.clone(),
                    tensor_type(1, rate.clone()),
                    dimensional_type(rate.clone()),
                    bc,
                    control,
                ]
            },
            | Scheme::Advection => {
                vec![
                    integer_vector(),
                    dimensional_type(l.clone()),
                    tensor_type(1, speed),
                    state.clone(),
                    bc,
                    control,
                ]
            },
            | Scheme::Fem(_) | Scheme::FemSource(_) => {
                let mut params = vec![
                    integer_vector(),
                    tensor_type(1, l.clone()),
                    tensor_type(1, l.clone()),
                    TypeExpr::Int,
                    dimensional_type(diffusivity),
                    dimensional_type(ratio(DimExpr::Dimless, t.clone())),
                    state.clone(),
                    bc,
                    control,
                ];
                if matches!(self.scheme, Scheme::FemSource(_)) {
                    params.push(TypeExpr::Function(
                        vec![dimensional_type(t.clone()), tensor_type(2, l.clone())],
                        Box::new(tensor_type(1, rate.clone())),
                    ));
                }
                params
            },
            | Scheme::Fvm(_) => {
                vec![
                    dimensional_type(l.clone()),
                    state.clone(),
                    TypeExpr::Bool,
                    control,
                    TypeExpr::Function(
                        vec![dimensional_type(u.clone())],
                        Box::new(dimensional_type(DimExpr::Mul(vec![
                            u.clone(),
                            speed.clone(),
                        ]))),
                    ),
                    TypeExpr::Function(
                        vec![dimensional_type(u.clone())],
                        Box::new(dimensional_type(speed)),
                    ),
                ]
            },
            | Scheme::Dg => {
                vec![
                    dimensional_type(l.clone()),
                    TypeExpr::Int,
                    TypeExpr::Int,
                    dimensional_type(speed),
                    state.clone(),
                    control,
                ]
            },
            | Scheme::DgNd => {
                vec![
                    integer_vector(),
                    tensor_type(1, l.clone()),
                    TypeExpr::Int,
                    tensor_type(1, speed),
                    state.clone(),
                    control,
                ]
            },
            | Scheme::Spectral => {
                vec![
                    integer_vector(),
                    tensor_type(1, l.clone()),
                    dimensional_type(diffusivity),
                    dimensional_type(ratio(DimExpr::Dimless, t.clone())),
                    state.clone(),
                    control,
                    TypeExpr::Function(
                        vec![dimensional_type(t.clone()), state.clone()],
                        Box::new(tensor_type(1, rate.clone())),
                    ),
                ]
            },
        };
        if self.event {
            args.push(TypeExpr::Int);
            let mut input = vec![dimensional_type(t.clone()), state];
            if matches!(self.scheme, Scheme::Wave) {
                input.push(tensor_type(1, rate.clone()));
            }
            args.push(TypeExpr::Function(input, Box::new(TypeExpr::Float)));
        }
        let mut signature = Signature::new(
            args,
            TypeExpr::Record(vec![
                ("times".into(), tensor_type(1, t.clone())),
                ("states".into(), tensor_type(2, u)),
                ("velocities".into(), tensor_type(2, rate)),
                ("coordinates".into(), tensor_type(2, l)),
                ("accepted".into(), TypeExpr::Int),
                ("rejected".into(), TypeExpr::Int),
                ("step_calls".into(), TypeExpr::Int),
                ("evaluations".into(), TypeExpr::Int),
                ("error_ratio".into(), TypeExpr::Float),
                (
                    "event_time".into(),
                    TypeExpr::Variant(vec![
                        ("None".into(), vec![]),
                        ("Some".into(), vec![dimensional_type(t)]),
                    ]),
                ),
            ]),
        );
        signature.type_params = ["T", "L", "U"]
            .into_iter()
            .map(|name| {
                TypeParam {
                    name: name.into(),
                    kind: TypeParamKind::Dim,
                    bounds: vec![],
                }
            })
            .collect();
        if self.event
            || matches!(
                self.scheme,
                Scheme::Fvm(_) | Scheme::Spectral | Scheme::FemSource(_)
            )
        {
            super::common::effectful(signature)
        } else {
            signature
        }
    }

    // Keep each reviewed dimensional preflight beside its selected solver dispatch.
    #[allow(clippy::too_many_lines)]
    fn eval(
        &self,
        args: &[Value],
        ctx: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let count = self.base_count() + if self.event { 2 } else { 0 };
        if args.len() != count {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(count as u8),
                got: args.len(),
            });
        }
        let ci = self.control_index();
        let control = record(&args[ci])?;
        let t = scalar_dimension(&ctx.value_type(field(control, "t1")?)?)?;
        let u = vector_dimension(&ctx.value_type(&args[self.initial_index()])?, None)?;
        let l = match self.scheme {
            | Scheme::Fem(_) | Scheme::FemSource(_) | Scheme::Spectral | Scheme::DgNd => {
                vector_dimension(&ctx.value_type(&args[1])?, None)?
            },
            | _ => {
                scalar_dimension(&ctx.value_type(
                    &args[usize::from(!matches!(self.scheme, Scheme::Dg | Scheme::Fvm(_)))],
                )?)?
            },
        };
        let (options, limit) = controls(&args[ci], &t, &u, ctx)?;
        let budget = Budget::new(limit);
        let initial = vector(
            &args[self.initial_index()],
            &u,
            None,
            self.initial_index(),
            ctx,
        )?;
        let n = initial.len();
        let speed = ratio(l.clone(), t.clone());
        let rate = ratio(u.clone(), t.clone());
        let diffusivity = ratio(l.clone().pow(2_i64), t.clone());
        let mut event_run = None;
        let mut event_direction = 0_i8;
        if self.event {
            let i = self.base_count();
            let Value::Int(direction) = &args[i] else {
                return Err(OpError::Domain(
                    "event direction must be -1, 0, or 1".into(),
                ));
            };
            event_direction = direction
                .to_i8()
                .filter(|d| (-1..=1).contains(d))
                .ok_or_else(|| OpError::Domain("event direction must be -1, 0, or 1".into()))?;
            let arity = 2 + u8::from(matches!(self.scheme, Scheme::Wave));
            let (run, sig) = CallbackRun::new(
                ctx,
                &args[i + 1],
                i + 1,
                arity,
                self.signature().params[i + 1].clone(),
            )?;
            same(&scalar_dimension(&sig.params[0])?, &t)?;
            same(&vector_dimension(&sig.params[1], Some(n))?, &u)?;
            if arity == 3 {
                same(&vector_dimension(&sig.params[2], Some(n))?, &rate)?;
            }
            same(&scalar_dimension(&sig.ret)?, &DimExpr::Dimless)?;
            event_run = Some(run);
        }
        let mut event_fn = |time: f64, state: &[f64]| -> Result<f64, FdmError> {
            let run = event_run
                .as_ref()
                .ok_or_else(|| FdmError::Domain("missing PDE event callback".into()))?;
            let mut values = vec![
                dimensional_float(time, t.clone()),
                output_tensor(vec![n], state[..n].to_vec(), u.clone())
                    .map_err(|e| FdmError::Domain(e.to_string()))?,
            ];
            if matches!(self.scheme, Scheme::Wave) {
                values.push(
                    output_tensor(vec![n], state[n..].to_vec(), rate.clone())
                        .map_err(|e| FdmError::Domain(e.to_string()))?,
                );
            }
            budget
                .invoke(run, &values, |value| {
                    scalar(&value, &DimExpr::Dimless, 0, ctx)
                })
                .ok_or_else(|| FdmError::Domain("PDE event callback failed".into()))
        };
        let mut terminal = PdeTerminalEvent {
            direction: event_direction,
            evaluate: &mut event_fn,
        };
        let event = if self.event {
            Some(&mut terminal)
        } else {
            None
        };
        let mut coordinates = Vec::new();
        let mut callbacks = Vec::new();
        let solution = match self.scheme {
            | Scheme::Heat(_) | Scheme::Wave | Scheme::Advection => {
                let extents = ints(&args[0], 0)?;
                if !(1..=3).contains(&extents.len()) {
                    return Err(OpError::Domain(
                        "PDE spatial rank must be 1, 2, or 3".into(),
                    ));
                }
                let h = scalar(&args[1], &l, 1, ctx)?;
                let bc_index = if matches!(self.scheme, Scheme::Wave) {
                    6
                } else {
                    4
                };
                let bc = boundary(&args[bc_index], extents.len(), &l, &u, ctx)?;
                let total = extents
                    .iter()
                    .try_fold(1usize, |n, d| n.checked_mul(*d))
                    .ok_or_else(|| OpError::Domain("PDE grid size overflows".into()))?;
                if total != n {
                    return Err(OpError::Domain(
                        "PDE initial length differs from the extent product".into(),
                    ));
                }
                for i in 0..n {
                    let mut offset = i;
                    for d in 0..3 {
                        if d < extents.len() {
                            coordinates.push((offset % extents[d]) as f64 * h);
                            offset /= extents[d];
                        } else {
                            coordinates.push(0.0);
                        }
                    }
                }
                match self.scheme {
                    | Scheme::Heat(method) => {
                        adaptive::heat_adaptive(
                            &extents,
                            h,
                            scalar(&args[2], &diffusivity, 2, ctx)?,
                            &initial,
                            &bc,
                            method,
                            &options,
                            event,
                        )
                    },
                    | Scheme::Wave => {
                        adaptive::wave_adaptive(
                            &extents,
                            h,
                            scalar(&args[2], &speed, 2, ctx)?,
                            &initial,
                            &vector(&args[4], &rate, Some(n), 4, ctx)?,
                            scalar(&args[5], &rate, 5, ctx)?,
                            &bc,
                            &options,
                            event,
                        )
                    },
                    | _ => {
                        adaptive::advection_adaptive(
                            &extents,
                            h,
                            &vector(&args[2], &speed, Some(extents.len()), 2, ctx)?,
                            &initial,
                            &bc,
                            &options,
                            event,
                        )
                    },
                }
            },
            | Scheme::Fem(theta) | Scheme::FemSource(theta) => {
                use crate::sim::boundary::MarkerBcs;
                use crate::sim::fem_simplex::FemProblem;
                use crate::sim::fem_simplex::FemSpace;
                use crate::sim::fem_simplex::SimplexMesh;
                use crate::sim::fem_simplex::marker;
                let cells = ints(&args[0], 0)?;
                if !(1..=3).contains(&cells.len()) || cells.contains(&0) {
                    return Err(OpError::Domain(
                        "FEM needs one to three positive cell counts".into(),
                    ));
                }
                let lo = vector(&args[1], &l, Some(cells.len()), 1, ctx)?;
                let hi = vector(&args[2], &l, Some(cells.len()), 2, ctx)?;
                if lo.iter().zip(&hi).any(|(a, b)| a >= b) {
                    return Err(OpError::Domain(
                        "FEM upper bounds must exceed lower bounds".into(),
                    ));
                }
                let degree = positive_budget(&args[3], 3, "FEM degree")?;
                let vertices = cells.iter().try_fold(1usize, |total, c| {
                    c.checked_add(1).and_then(|n| total.checked_mul(n))
                });
                if !(1..=2).contains(&degree) || vertices.is_none_or(|v| v > n) {
                    return Err(OpError::Domain(
                        "FEM degree must be 1 or 2 and initial data must cover the mesh vertices"
                            .into(),
                    ));
                }
                let bc = boundary(&args[7], cells.len(), &l, &u, ctx)?;
                let mesh = match cells.len() {
                    | 1 => SimplexMesh::interval(cells[0], lo[0], hi[0]),
                    | 2 => {
                        SimplexMesh::rectangle(cells[0], cells[1], [lo[0], lo[1]], [hi[0], hi[1]])
                    },
                    | _ => {
                        SimplexMesh::boxed(
                            [cells[0], cells[1], cells[2]],
                            [lo[0], lo[1], lo[2]],
                            [hi[0], hi[1], hi[2]],
                        )
                    },
                };
                let space = FemSpace::new(mesh, degree).map_err(OpError::Domain)?;
                if space.dofs.len() != n {
                    return Err(OpError::Domain("FEM initial field must contain one entry per vertex/edge degree of freedom".into()));
                }
                coordinates = space.dofs.iter().flatten().copied().collect();
                let faces = [bc.x_min, bc.x_max, bc.y_min, bc.y_max, bc.z_min, bc.z_max];
                let markers = [
                    marker::XMIN,
                    marker::XMAX,
                    marker::YMIN,
                    marker::YMAX,
                    marker::ZMIN,
                    marker::ZMAX,
                ];
                let mut bcs = MarkerBcs::new();
                for i in 0..2 * cells.len() {
                    bcs = bcs.set(
                        markers[i],
                        pde_api::to_rssn_bc(&faces[i]).map_err(kernel_error)?,
                    );
                }
                let source = |_: f64, _: [f64; 3]| 0.0;
                let problem = FemProblem {
                    k: scalar(&args[4], &diffusivity, 4, ctx)?,
                    c: scalar(&args[5], &ratio(DimExpr::Dimless, t.clone()), 5, ctx)?,
                    f: &source,
                    bcs: &bcs,
                };
                if matches!(self.scheme, Scheme::FemSource(_)) {
                    let (run, sig) =
                        CallbackRun::new(ctx, &args[9], 9, 2, self.signature().params[9].clone())?;
                    same(&scalar_dimension(&sig.params[0])?, &t)?;
                    same(&matrix_dimension(&sig.params[1], n, 3)?, &l)?;
                    same(&vector_dimension(&sig.ret, Some(n))?, &rate)?;
                    let nodal_source = |time: f64, coords: &[[f64; 3]]| {
                        budget
                            .invoke(
                                &run,
                                &[
                                    dimensional_float(time, t.clone()),
                                    output_tensor(
                                        vec![n, 3],
                                        coords.iter().flatten().copied().collect(),
                                        l.clone(),
                                    )
                                    .unwrap_or(Value::Float(f64::NAN)),
                                ],
                                |v| vector(&v, &rate, Some(n), 0, ctx),
                            )
                            .ok_or_else(|| FdmError::Domain("FEM source callback failed".into()))
                    };
                    let result = pde_api::adaptive::fem_heat_source_adaptive(
                        &space,
                        &problem,
                        &initial,
                        theta,
                        nodal_source,
                        &options,
                        event,
                    );
                    callbacks.push(run);
                    result
                } else {
                    pde_api::adaptive::fem_heat_adaptive(
                        &space, &problem, &initial, theta, &options, event,
                    )
                }
            },
            | Scheme::Fvm(scheme) => {
                let length = scalar(&args[0], &l, 0, ctx)?;
                let Value::Bool(periodic) = &args[2] else {
                    return Err(OpError::Domain("FVM periodic flag must be Boolean".into()));
                };
                let (flux, fsig) =
                    CallbackRun::new(ctx, &args[4], 4, 1, self.signature().params[4].clone())?;
                let (waves, wsig) =
                    CallbackRun::new(ctx, &args[5], 5, 1, self.signature().params[5].clone())?;
                let flux_dimension = DimExpr::Mul(vec![u.clone(), speed.clone()]);
                callback_scalar(&fsig, &u, &flux_dimension)?;
                callback_scalar(&wsig, &u, &speed)?;
                let f = |value| {
                    budget
                        .invoke(&flux, &[dimensional_float(value, u.clone())], |result| {
                            scalar(&result, &flux_dimension, 0, ctx)
                        })
                        .unwrap_or(f64::NAN)
                };
                let g = |value| {
                    budget
                        .invoke(&waves, &[dimensional_float(value, u.clone())], |result| {
                            scalar(&result, &speed, 0, ctx)
                        })
                        .unwrap_or(f64::NAN)
                };
                for i in 0..n {
                    coordinates.extend([(i as f64 + 0.5) * length / n as f64, 0.0, 0.0]);
                }
                let result = pde_api::adaptive::fvm_adaptive(
                    length, &initial, *periodic, scheme, f, g, &options, event,
                );
                callbacks.push(flux);
                callbacks.push(waves);
                result
            },
            | Scheme::Dg => {
                let length = scalar(&args[0], &l, 0, ctx)?;
                let elements = positive_budget(&args[1], 1, "DG elements")?;
                let degree = positive_budget(&args[2], 2, "DG degree")?;
                if !length.is_finite()
                    || length <= 0.0
                    || degree.checked_add(1).and_then(|p| elements.checked_mul(p)) != Some(n)
                {
                    return Err(OpError::Domain("DG initial data must match elements*(degree+1) on a positive finite domain".into()));
                }
                let mesh = crate::sim::dg::DgMesh::line(elements, 0.0, length, degree);
                let (nodes, _) = crate::sim::dg::lgl(degree);
                for e in 0..elements {
                    for xi in &nodes {
                        coordinates.extend([
                            (e as f64 + f64::midpoint(1.0, *xi)) * length / elements as f64,
                            0.0,
                            0.0,
                        ]);
                    }
                }
                pde_api::adaptive::dg_advection_adaptive(
                    &mesh,
                    scalar(&args[3], &speed, 3, ctx)?,
                    &initial,
                    &options,
                    event,
                )
            },
            | Scheme::DgNd => {
                let cells = ints(&args[0], 0)?;
                if !(1..=2).contains(&cells.len()) || cells.contains(&0) {
                    return Err(OpError::Domain(
                        "DG needs one or two positive element counts".into(),
                    ));
                }
                let lengths = vector(&args[1], &l, Some(cells.len()), 1, ctx)?;
                let velocities = vector(&args[3], &speed, Some(cells.len()), 3, ctx)?;
                let degree = positive_budget(&args[2], 2, "DG degree")?;
                let elements = cells.iter().try_fold(1usize, |n, d| n.checked_mul(*d));
                let nodes = degree
                    .checked_add(1)
                    .and_then(|p| p.checked_pow(cells.len() as u32));
                if lengths.iter().any(|v| *v <= 0.0)
                    || elements.and_then(|e| nodes.and_then(|p| e.checked_mul(p))) != Some(n)
                {
                    return Err(OpError::Domain(
                        "DG nodal shape or physical box is invalid".into(),
                    ));
                }
                let mesh = if cells.len() == 1 {
                    crate::sim::dg::DgMesh::line(cells[0], 0.0, lengths[0], degree)
                } else {
                    crate::sim::dg::DgMesh::rect(
                        [cells[0], cells[1]],
                        [0.0, 0.0],
                        [lengths[0], lengths[1]],
                        degree,
                    )
                };
                let (xi, _) = crate::sim::dg::lgl(degree);
                let p = degree + 1;
                for e in 0..elements.unwrap_or(0) {
                    for node in 0..nodes.unwrap_or(0) {
                        let x = (e % cells[0]) as f64 + f64::midpoint(1.0, xi[node % p]);
                        let y = if cells.len() == 2 {
                            ((e / cells[0]) as f64 + f64::midpoint(1.0, xi[node / p])) * lengths[1]
                                / cells[1] as f64
                        } else {
                            0.0
                        };
                        coordinates.extend([x * lengths[0] / cells[0] as f64, y, 0.0]);
                    }
                }
                pde_api::adaptive::dg_advection_nd_adaptive(
                    &mesh,
                    [velocities[0], velocities.get(1).copied().unwrap_or(0.0)],
                    &initial,
                    &options,
                    event,
                )
            },
            | Scheme::Spectral => {
                let extents = ints(&args[0], 0)?;
                let lengths = vector(&args[1], &l, Some(extents.len()), 1, ctx)?;
                if extents.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(n) {
                    return Err(OpError::Domain(
                        "spectral initial data does not match the extent product".into(),
                    ));
                }
                if lengths.iter().any(|v| *v <= 0.0) {
                    return Err(OpError::Domain("spectral lengths must be positive".into()));
                }
                let (run, sig) =
                    CallbackRun::new(ctx, &args[6], 6, 2, self.signature().params[6].clone())?;
                same(&scalar_dimension(&sig.params[0])?, &t)?;
                same(&vector_dimension(&sig.params[1], Some(n))?, &u)?;
                same(&vector_dimension(&sig.ret, Some(n))?, &rate)?;
                let nonlinearity = |time, state: &[f64]| {
                    budget
                        .invoke(
                            &run,
                            &[
                                dimensional_float(time, t.clone()),
                                output_tensor(vec![n], state.to_vec(), u.clone())
                                    .unwrap_or(Value::Float(f64::NAN)),
                            ],
                            |value| vector(&value, &rate, Some(n), 0, ctx),
                        )
                        .unwrap_or_else(|| vec![f64::NAN; n])
                };
                let grid = crate::sim::spectral::SpectralGrid::new(&extents, &lengths)
                    .map_err(OpError::Domain)?;
                for i in 0..n {
                    coordinates.extend(grid.coord(i));
                }
                let result = pde_api::adaptive::spectral_etdrk4_adaptive(
                    &extents,
                    &lengths,
                    scalar(&args[2], &diffusivity, 2, ctx)?,
                    scalar(&args[3], &ratio(DimExpr::Dimless, t.clone()), 3, ctx)?,
                    &initial,
                    nonlinearity,
                    &options,
                    event,
                );
                callbacks.push(run);
                result
            },
        };
        for callback in callbacks {
            callback.finish()?;
        }
        if let Some(callback) = event_run {
            callback.finish()?;
        }
        budget.finish()?;
        let result = solution.map_err(kernel_error)?;
        encode(
            result,
            n,
            matches!(self.scheme, Scheme::Wave),
            coordinates,
            &t,
            &l,
            &u,
            budget.count.get(),
        )
    }
}
#[allow(clippy::too_many_arguments)]
fn encode(
    result: PdeAdaptiveSolution,
    n: usize,
    wave: bool,
    coordinates: Vec<f64>,
    t: &DimExpr,
    l: &DimExpr,
    u: &DimExpr,
    evaluations: usize,
) -> Result<Value, OpError> {
    let samples = result.times.len();
    let mut values = Vec::new();
    let mut velocities = Vec::new();
    for state in result.states {
        if state.len() != if wave { 2 * n } else { n } {
            return Err(OpError::Domain("PDE result state length changed".into()));
        }
        values.extend_from_slice(&state[..n]);
        if wave {
            velocities.extend_from_slice(&state[n..]);
        }
    }
    let event = Value::Variant(VariantValue {
        tag: if result.event_time.is_some() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: result
            .event_time
            .into_iter()
            .map(|v| dimensional_float(v, t.clone()))
            .collect(),
    });
    Ok(Value::Record(RecordValue {
        fields: vec![
            (
                "times".into(),
                output_tensor(vec![samples], result.times, t.clone())?,
            ),
            (
                "states".into(),
                output_tensor(vec![samples, n], values, u.clone())?,
            ),
            (
                "velocities".into(),
                output_tensor(
                    if wave {
                        vec![samples, n]
                    } else {
                        vec![0, 0]
                    },
                    velocities,
                    ratio(u.clone(), t.clone()),
                )?,
            ),
            (
                "coordinates".into(),
                output_tensor(vec![n, 3], coordinates, l.clone())?,
            ),
            ("accepted".into(), Value::int(result.accepted)),
            ("rejected".into(), Value::int(result.rejected)),
            ("step_calls".into(), Value::int(result.step_calls)),
            ("evaluations".into(), Value::int(evaluations)),
            ("error_ratio".into(), Value::Float(result.error_ratio)),
            ("event_time".into(), event),
        ],
    }))
}

pub(super) fn install(installer: &mut Installer<'_>) -> Result<(), RuleError> {
    for (name, path, scheme) in [
        (
            "num_pde_heat_explicit_adaptive",
            "crate::sim::fdm_nd::adaptive::heat_adaptive",
            Scheme::Heat(TimeStepping::Explicit),
        ),
        (
            "num_pde_heat_be_adaptive",
            "crate::sim::fdm_nd::adaptive::heat_adaptive",
            Scheme::Heat(TimeStepping::BackwardEuler),
        ),
        (
            "num_pde_heat_cn_adaptive",
            "crate::sim::fdm_nd::adaptive::heat_adaptive",
            Scheme::Heat(TimeStepping::CrankNicolson),
        ),
        (
            "num_pde_heat_adi_adaptive",
            "crate::sim::fdm_nd::adaptive::heat_adaptive",
            Scheme::Heat(TimeStepping::Adi),
        ),
        (
            "num_pde_wave_adaptive",
            "crate::sim::fdm_nd::adaptive::wave_adaptive",
            Scheme::Wave,
        ),
        (
            "num_pde_advection_adaptive",
            "crate::sim::fdm_nd::adaptive::advection_adaptive",
            Scheme::Advection,
        ),
        (
            "num_pde_fem_be_adaptive",
            "crate::sim::pde_api::adaptive::fem_heat_adaptive",
            Scheme::Fem(1.0),
        ),
        (
            "num_pde_fem_cn_adaptive",
            "crate::sim::pde_api::adaptive::fem_heat_adaptive",
            Scheme::Fem(0.5),
        ),
        (
            "num_pde_fem_be_source_adaptive",
            "crate::sim::pde_api::adaptive::fem_heat_source_adaptive",
            Scheme::FemSource(1.0),
        ),
        (
            "num_pde_fem_cn_source_adaptive",
            "crate::sim::pde_api::adaptive::fem_heat_source_adaptive",
            Scheme::FemSource(0.5),
        ),
        (
            "num_pde_dg_nd_adaptive",
            "crate::sim::pde_api::adaptive::dg_advection_nd_adaptive",
            Scheme::DgNd,
        ),
        (
            "num_pde_fvm_upwind_adaptive",
            "crate::sim::pde_api::adaptive::fvm_adaptive",
            Scheme::Fvm(FvmScheme::Upwind),
        ),
        (
            "num_pde_fvm_lax_friedrichs_adaptive",
            "crate::sim::pde_api::adaptive::fvm_adaptive",
            Scheme::Fvm(FvmScheme::LaxFriedrichs),
        ),
        (
            "num_pde_fvm_muscl_minmod_adaptive",
            "crate::sim::pde_api::adaptive::fvm_adaptive",
            Scheme::Fvm(FvmScheme::MusclMinmod),
        ),
        (
            "num_pde_fvm_muscl_van_leer_adaptive",
            "crate::sim::pde_api::adaptive::fvm_adaptive",
            Scheme::Fvm(FvmScheme::MusclVanLeer),
        ),
        (
            "num_pde_dg_adaptive",
            "crate::sim::pde_api::adaptive::dg_advection_adaptive",
            Scheme::Dg,
        ),
        (
            "num_pde_spectral_etdrk4_adaptive",
            "crate::sim::pde_api::adaptive::spectral_etdrk4_adaptive",
            Scheme::Spectral,
        ),
    ] {
        for event in [false, true] {
            let pde = Pde { scheme, event };
            let name = if event {
                format!("{name}_event")
            } else {
                name.into()
            };
            let callbacks = event
                || matches!(
                    scheme,
                    Scheme::Fvm(_) | Scheme::Spectral | Scheme::FemSource(_)
                );
            let signature = pde.signature();
            let count = signature.params.len() as u8;
            let mut descriptor=OpDescriptor::new(&name,Arity::Fixed(count)).rust_path(path).category("numerical_callback")
                .doc("Tolerance-controlled PDE evolution with explicit coordinate L, time T, and field U dimensions. Diffusivity is L^2/T; velocity L/T; reaction 1/T; Neumann U/L; Robin alpha dimensionless, beta L, gamma U. Grid/state lengths and boundary pairing/corner data are checked. Controls: t0,t1,rtol,atol,h0,h_min,h_max,max_steps,max_rejections,max_step_calls,max_evaluations. Boundary record: kind,value,flux,alpha,beta,gamma vectors (one broadcast face or two per active axis). Only accepted states are returned; step doubling controls temporal error, independently of spatial error. Optional trailing direction and dimensionless terminal-event callback; wave events receive displacement and velocity separately. FEM uses a fixed simplex mesh; source variants take a vectorized (time,coordinates)->nodal-rate callback with a consistent mass load. FVM preserves periodic/transmissive flux methods; classic LF's timestep-dependent viscosity may prevent strict temporal convergence. DG uses periodic SSPRK3 in one or two dimensions; Fourier uses ETDRK4 and a dimensioned nonlinear rate callback. Budgets and callback failures are explicit errors.")
                .stability(Stability::Experimental).purity(if callbacks{Purity::Effectful}else{Purity::Pure}).determinism(Determinism::EnvDependent).thread_safety(ThreadSafety::SessionBound)
                .backend(BackendSupport{interpreter:true,jit:false,ffi:true,gpu:false}).signature(signature).eval_dyn(Arc::new(pde));
            if callbacks {
                descriptor = descriptor.higher_order();
            }
            let op = installer.op(descriptor)?;
            installer.dim(op, crate::graph::dim::DimRule::Container);
        }
    }
    Ok(())
}
