//! Native geometric and algebraic multigrid hierarchies with operator reference domains.
use super::super::native_object_api::Maybe;
use super::Array;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::LIMIT;
use super::NativeType;
use super::NewObject;
use super::Normalized;
use super::OpError;
use super::Output;
use super::RuleError;
use super::api;
use super::array;
use super::dim;
use super::err;
use super::object_binding;
use super::product;
use super::scale;
use super::sparse::Csr;
use super::sparse::inverse;
use super::sparse::result;
use super::sparse::vector;
use crate::kernel::signature::DimExpr;
use crate::kernels::amg_general as ag;
use crate::kernels::multigrid as mg;
use crate::kernels::reuse::f64_scratch;
use ndarray::ArrayD;
use ndarray::IxDyn;
struct Amg<T> {
    data: T,
    n: usize,
    scale: f64,
}
macro_rules! hierarchy {
    ($t:path,$kind:literal) => {
        impl NativeType for Amg<$t> {
            const KIND: &'static str = $kind;
            const MAX_BYTES: usize = 16 * 1024 * 1024;

            fn validate(&self) -> Result<(), OpError> {
                if self.n == 0 || self.n > 512 || !self.scale.is_finite() || self.scale <= 0.0 {
                    return Err(err("AMG hierarchy size/reference budget exceeded"));
                }
                Ok(())
            }

            fn retained_bytes(&self) -> usize {
                size_of::<Self>() + self.data.heap_bytes()
            }
        }
    };
}
hierarchy!(mg::Amg, "rssn/native/amg/v1");
hierarchy!(ag::AmgGeneral, "rssn/native/amg-general/v1");
const GEOMETRIC_MAX_SIDE: usize = 127;
const GEOMETRIC_MAX_POINTS: usize = 32_768;
struct Geometric {
    data: mg::GeometricMg,
    scale: f64,
    workspace: mg::GeometricTransferWorkspace,
}

// Above side 31, the next 3D grid exceeds the hierarchy budget. Each
// direction needs nf²*nc and nf*nc² slots, with the two buffers exchanged
// by interpolation. After both directions each buffer needs at most
// 31²*15 slots. The single 2D buffer needs at most 127*63, and 1D uses none.
const _: () = {
    const MAX_3D_SIDE: usize = 31;
    const MAX_TRANSFER_SLOTS: usize = 2 * MAX_3D_SIDE.pow(2) * (MAX_3D_SIDE / 2);
    assert!((2 * MAX_3D_SIDE + 1).pow(3) > GEOMETRIC_MAX_POINTS);
    assert!(GEOMETRIC_MAX_SIDE * (GEOMETRIC_MAX_SIDE / 2) <= MAX_TRANSFER_SLOTS);
    assert!(size_of::<Geometric>() + MAX_TRANSFER_SLOTS * size_of::<f64>() <= Geometric::MAX_BYTES);
};

impl Geometric {
    const fn workspace_fits(&self) -> bool {
        self.workspace.heap_bytes() <= Self::MAX_BYTES - size_of::<Self>()
    }
}

impl NativeType for Geometric {
    const KIND: &'static str = "rssn/native/geometric-mg/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.size() > GEOMETRIC_MAX_POINTS
            || !self.scale.is_finite()
            || self.scale <= 0.0
            || !self.workspace_fits()
        {
            return Err(err("geometric hierarchy size/reference budget exceeded"));
        }
        Ok(())
    }
}

fn geometric_level_size(
    hierarchy: &Geometric,
    level: usize,
) -> Result<usize, OpError> {
    hierarchy
        .data
        .level_size(level)
        .ok_or_else(|| err("multigrid level out of range"))
}

#[derive(Clone, Copy)]
enum GridTransfer {
    Restrict,
    Prolong,
}

fn geometric_transfer(
    hierarchy: &mut Geometric,
    level: usize,
    values: &Array,
    dimension: DimExpr,
    direction: GridTransfer,
) -> Result<NewObject<Array>, OpError> {
    let fine = geometric_level_size(hierarchy, level)?;
    let coarse = geometric_level_size(hierarchy, level + 1)?;
    let (source_len, destination_len) = match direction {
        | GridTransfer::Restrict => (fine, coarse),
        | GridTransfer::Prolong => (coarse, fine),
    };
    let input = vector(values, source_len)?;
    let mut output = vec![0.0; destination_len];
    let transferred = match direction {
        | GridTransfer::Restrict => {
            hierarchy
                .data
                .restrict_into(level, input, &mut output, &mut hierarchy.workspace)
        },
        | GridTransfer::Prolong => {
            hierarchy
                .data
                .prolong_into(level, input, &mut output, &mut hierarchy.workspace)
        },
    };
    // The arena reserves the existing maximum at creation. Check actual
    // capacity before releasing exclusive access, keeping only workspaces
    // that fit. An oversized transient buffer still serves the current call.
    if !hierarchy.workspace_fits() {
        hierarchy.workspace = mg::GeometricTransferWorkspace::default();
    }
    transferred.map_err(err)?;
    result(output, values.scale, dimension)
}

object_binding!(
    mg_geometric_level_size,
    "num_linear_geometric_level_size",
    "crate::kernels::multigrid::GeometricMg::level_size",
    "Returns the number of interior samples at a zero-based geometric level. Level zero is the finest grid.",
    true,
    (hierarchy: Handle<Geometric>, level: usize) -> usize,
    |ctx, hierarchy, level| hierarchy.with(ctx, |a| geometric_level_size(a, level))
);

object_binding!(
    mg_geometric_restrict,
    "num_linear_geometric_restrict",
    "crate::kernels::multigrid::GeometricMg::restrict_into",
    "Restricts a rank-one array from the specified fine level to its next coarser level by full weighting. The result keeps the input array's physical dimension and reference scale.",
    true,
    (hierarchy: Handle<Geometric>, level: usize, values: Handle<Array>) -> NewObject<Array>,
    |ctx, hierarchy, level, values| {
        let dimension = dim(ctx, &values.value)?;
        hierarchy.with_mut(ctx, |a| values.with(ctx, |u| {
            geometric_transfer(a, level, u, dimension, GridTransfer::Restrict)
        }))
    }
);

object_binding!(
    mg_geometric_prolong,
    "num_linear_geometric_prolong",
    "crate::kernels::multigrid::GeometricMg::prolong_into",
    "Interpolates a rank-one array from the next coarser level to the specified fine level, with homogeneous Dirichlet boundaries. The result keeps the input array's physical dimension and reference scale.",
    true,
    (hierarchy: Handle<Geometric>, level: usize, values: Handle<Array>) -> NewObject<Array>,
    |ctx, hierarchy, level, values| {
        let dimension = dim(ctx, &values.value)?;
        hierarchy.with_mut(ctx, |a| values.with(ctx, |u| {
            geometric_transfer(a, level, u, dimension, GridTransfer::Prolong)
        }))
    }
);

struct MgResult {
    data: mg::MgResult,
    scale: f64,
}
impl NativeType for MgResult {
    const KIND: &'static str = "rssn/native/multigrid-result/v1";
    const MAX_BYTES: usize = LIMIT * 8 + 65536;

    fn validate(&self) -> Result<(), OpError> {
        if self.data.x.len() > LIMIT
            || self.data.history.len() > 2049
            || self
                .data
                .x
                .iter()
                .chain(&self.data.history)
                .any(|v| !v.is_finite())
            || !self.scale.is_finite()
            || self.scale <= 0.0
        {
            return Err(err("invalid finite multigrid result or work budget"));
        }
        Ok(())
    }
}
fn cycle(name: &str) -> Result<mg::Cycle, OpError> {
    match name {
        | "V" => Ok(mg::Cycle::V),
        | "W" => Ok(mg::Cycle::W),
        | _ => Err(err("cycle must be V or W")),
    }
}

object_binding!(
    mg_amg_general_workspace_len,
    "num_linear_amg_general_workspace_len",
    "crate::kernels::amg_general::AmgGeneral::workspace_len",
    "Reports the number of reusable Float slots needed by one cycle of the general AMG hierarchy. A terminal-only hierarchy requires zero slots.",
    true,
    (hierarchy: Handle<Amg<ag::AmgGeneral>>) -> usize,
    |ctx, hierarchy| hierarchy.with(ctx, |a| Ok(a.data.workspace_len()))
);

object_binding!(
    mg_amg_general_precondition_many,
    "num_linear_amg_general_precondition_many",
    "crate::kernels::amg_general::AmgGeneral::precondition_with_workspace",
    "Applies the general AMG preconditioner from zero to each row of a rank-two RHS array. Rows share one hierarchy and one reusable work buffer. The result keeps the input shape and order, with reference scale and dimension equal to the RHS divided by the operator.",
    true,
    (hierarchy: Handle<Amg<ag::AmgGeneral>>, rhs: Handle<Array>) -> NewObject<Array>,
    |ctx, hierarchy, rhs| {
        let result_dimension = product(dim(ctx, &rhs.value)?, inverse(dim(ctx, &hierarchy.value)?)?)?;
        hierarchy.with(ctx, |a| rhs.with(ctx, |b| {
            let [rows, cols] = b.data.shape() else {
                return Err(err("AMG batch right-hand sides must have rank two"));
            };
            if *cols != a.n {
                return Err(err("AMG batch row length must match the operator"));
            }
            let b_values = b.data.as_standard_layout();
            let values = b_values.as_slice().ok_or_else(|| err("AMG batch array layout is invalid"))?;
            let mut out = vec![0.0; values.len()];
            let workspace_len = a.data.workspace_len();
            let mut work = f64_scratch(workspace_len);
            work.resize(workspace_len, 0.0);
            for (source, destination) in values.chunks_exact(a.n).zip(out.chunks_exact_mut(a.n)) {
                a.data.precondition_with_workspace(source, destination, &mut work).map_err(err)?;
            }
            let data = ArrayD::from_shape_vec(IxDyn(&[*rows, *cols]), out).map_err(err)?;
            let reference = api::scale_without_underflow(b.scale, &[], &[a.scale])?;
            Ok(array(data, reference, result_dimension))
        }))
    }
);

fn controls(
    name: &str,
    smoother: &str,
    pre: usize,
    post: usize,
    tolerance: f64,
    cycles: usize,
) -> Result<mg::MgOptions, OpError> {
    if pre > 8 || post > 8 || cycles == 0 || cycles > 2048 || tolerance <= 0.0 {
        return Err(err("multigrid control work budget exceeded"));
    }
    let smoother = match smoother {
        | "Jacobi" => mg::Smoother::Jacobi,
        | "RedBlackGaussSeidel" => mg::Smoother::RedBlackGaussSeidel,
        | _ => return Err(err("invalid native geometric smoother")),
    };
    Ok(mg::MgOptions {
        cycle: cycle(name)?,
        pre,
        post,
        smoother,
        tol: tolerance,
        max_cycles: cycles,
    })
}
fn setup(
    a: &Csr,
    theta: f64,
    coarse: usize,
    levels: usize,
    sweeps: usize,
    near: &[Vec<f64>],
) -> Result<(), OpError> {
    if a.data.rows == 0
        || a.data.rows != a.data.cols
        || a.data.rows > 512
        || !(0.0..=1.0).contains(&theta)
        || coarse == 0
        || coarse > 512
        || levels == 0
        || levels > 10
        || sweeps > 8
        || near.len() > 16
        || near.iter().any(|v| v.len() != a.data.rows)
    {
        return Err(err("AMG setup shape/options/work budget exceeded"));
    }
    Ok(())
}
object_binding!(mg_amg_build,"num_linear_amg_build","crate::kernels::multigrid::Amg::build","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Csr>,theta:f64,coarse_size:usize,max_levels:usize,sweeps:usize,cycle_name:String,near_nullspace:Vec<Vec<f64>>,smooth:bool)->NewObject<Amg<mg::Amg>>,|ctx,a,theta,coarse_size,max_levels,sweeps,cycle_name,near_nullspace,smooth|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{setup(a,theta,coarse_size,max_levels,sweeps,&near_nullspace)?;let options=mg::AmgOptions{theta,coarse_size,max_levels,sweeps,cycle:cycle(&cycle_name)?,near_nullspace,smooth_prolongator:smooth};let data=mg::Amg::build(&a.data,&options).map_err(|e|err(format!("{e:?}")))?;Ok(NewObject::new(Amg{data,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(mg_amg_general_build,"num_linear_amg_general_build","crate::kernels::amg_general::AmgGeneral::build","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Csr>,kind:String,theta:f64,coarse_size:usize,max_levels:usize,sweeps:usize,cycle_name:String,near_nullspace:Vec<Vec<f64>>)->NewObject<Amg<ag::AmgGeneral>>,|ctx,a,kind,theta,coarse_size,max_levels,sweeps,cycle_name,near_nullspace|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|{setup(a,theta,coarse_size,max_levels,sweeps,&near_nullspace)?;let kind=match kind.as_str(){"RugeStuben"=>ag::AmgKind::RugeStuben,"SmoothedAggregation"=>ag::AmgKind::SmoothedAggregation,_=>return Err(err("invalid AMG coarsening kind"))};let options=ag::AmgGeneralOptions{kind,theta,coarse_size,max_levels,sweeps,cycle:cycle(&cycle_name)?,near_nullspace};let data=ag::AmgGeneral::build(&a.data,&options).map_err(|e|err(format!("{e:?}")))?;Ok(NewObject::new(Amg{data,n:a.data.rows,scale:a.scale},vec![d]))})});
object_binding!(mg_amg_num_levels,"num_linear_amg_num_levels","crate::kernels::multigrid::Amg::num_levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.num_levels()))});
object_binding!(mg_amg_operator_complexity,"num_linear_amg_operator_complexity","crate::kernels::multigrid::Amg::operator_complexity","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.operator_complexity()))});
object_binding!(mg_amg_precondition,"num_linear_amg_precondition","crate::kernels::multigrid::Amg::precondition","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.precondition(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_vcycle,"num_linear_amg_vcycle","crate::kernels::multigrid::Amg::vcycle","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>,initial:Vec<f64>)->NewObject<Array>,|ctx,a,rhs,initial|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.len()!=a.n{return Err(err("initial iterate length mismatch"));}let mut out=initial;a.data.vcycle(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_solve,"num_linear_amg_solve","crate::kernels::multigrid::Amg::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<mg::Amg>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;if tolerance<=0.0||max_cycles==0||max_cycles>2048{return Err(err("AMG solve control/work budget invalid"));}a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.as_ref().is_some_and(|v|v.len()!=a.n){return Err(err("initial iterate length mismatch"));}Ok(NewObject::new(MgResult{data:a.data.solve(input,initial.as_deref(),tolerance,max_cycles),scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_amg_general_num_levels,"num_linear_amg_general_num_levels","crate::kernels::amg_general::AmgGeneral::num_levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.num_levels()))});
object_binding!(mg_amg_general_operator_complexity,"num_linear_amg_general_operator_complexity","crate::kernels::amg_general::AmgGeneral::operator_complexity","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.operator_complexity()))});
object_binding!(mg_amg_general_precondition,"num_linear_amg_general_precondition","crate::kernels::amg_general::AmgGeneral::precondition","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>)->NewObject<Array>,|ctx,a,rhs|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;let mut out=vec![0.0;a.n];a.data.precondition(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_general_vcycle,"num_linear_amg_general_vcycle","crate::kernels::amg_general::AmgGeneral::vcycle","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>,initial:Vec<f64>)->NewObject<Array>,|ctx,a,rhs,initial|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.len()!=a.n{return Err(err("initial iterate length mismatch"));}let mut out=initial;a.data.vcycle(input,&mut out);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
object_binding!(mg_amg_general_solve,"num_linear_amg_general_solve","crate::kernels::amg_general::AmgGeneral::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Amg<ag::AmgGeneral>>,rhs:Handle<Array>,initial:Option<Vec<f64>>,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;if tolerance<=0.0||max_cycles==0||max_cycles>2048{return Err(err("AMG solve control/work budget invalid"));}a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.n)?;if initial.as_ref().is_some_and(|v|v.len()!=a.n){return Err(err("initial iterate length mismatch"));}Ok(NewObject::new(MgResult{data:a.data.solve(input,initial.as_deref(),tolerance,max_cycles),scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_multigrid_mean_factor,"num_linear_multigrid_mean_factor","crate::kernels::multigrid::MgResult::mean_factor","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<MgResult>)->f64,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.mean_factor()))});
object_binding!(mg_multigrid_result,"num_linear_multigrid_result","","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<MgResult>)->Dynamic,|ctx,a|{let _=ctx;let d=dim(ctx,&a.value)?;a.with(ctx,|a|Ok(Dynamic(api::record(vec![("solution",a.data.x.clone().encode()?),("reference",api::quantity(a.scale,d)?),("history",a.data.history.clone().encode()?),("converged",a.data.converged.encode()?)]))))});
object_binding!(
    mg_geometric_new,
    "num_linear_geometric_new",
    "crate::kernels::multigrid::GeometricMg::new",
    "Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",
    false,
    (dimension: usize, n: usize, shift: Normalized<f64>, reference: Input) -> Maybe<NewObject<Geometric>>,
    |_, dimension, n, shift, reference| {
        if dimension == 0
            || dimension > 3
            || n > GEOMETRIC_MAX_SIDE
            || n.checked_pow(dimension as u32).is_none_or(|size| size > GEOMETRIC_MAX_POINTS)
        {
            return Err(err("geometric hierarchy extent/work budget exceeded"));
        }
        let (scale, dimension_ref) = api::reference(&reference.0)?;
        Ok(Maybe(mg::GeometricMg::new(dimension, n, shift.0).map(|data| {
            let hierarchy = Geometric {
                data,
                scale,
                workspace: mg::GeometricTransferWorkspace::default(),
            };
            NewObject::new(hierarchy, vec![dimension_ref])
        })))
    }
);
object_binding!(mg_geometric_levels,"num_linear_geometric_levels","crate::kernels::multigrid::GeometricMg::levels","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.levels()))});
object_binding!(mg_geometric_size,"num_linear_geometric_size","crate::kernels::multigrid::GeometricMg::size","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>)->usize,|ctx,a|{let _=ctx;a.with(ctx,|a|Ok(a.data.size()))});
object_binding!(
    mg_geometric_apply,
    "num_linear_geometric_apply",
    "crate::kernels::multigrid::GeometricMg::apply",
    "Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",
    true,
    (a: Handle<Geometric>, level: usize, values: Handle<Array>) -> NewObject<Array>,
    |ctx, a, level, values| {
        let d = product(dim(ctx, &a.value)?, dim(ctx, &values.value)?)?;
        a.with(ctx, |a| values.with(ctx, |u| {
            let size = geometric_level_size(a, level)?;
            let input = vector(u, size)?;
            let mut out = vec![0.0; size];
            a.data.apply(level, input, &mut out);
            result(out, scale(a.scale, u.scale)?, d)
        }))
    }
);
object_binding!(mg_geometric_residual_norm,"num_linear_geometric_residual_norm","crate::kernels::multigrid::GeometricMg::residual_norm","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,values:Vec<f64>,rhs:Handle<Array>)->Dynamic,|ctx,a,values,rhs|{let _=ctx;let d=dim(ctx,&rhs.value)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if values.len()!=input.len(){return Err(err("initial iterate length mismatch"));}Ok(Dynamic(api::quantity(a.data.residual_norm(&values,input)*r.scale,d)?))}))});
object_binding!(mg_geometric_solve,"num_linear_geometric_solve","crate::kernels::multigrid::GeometricMg::solve","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let data=a.data.solve(input,initial.as_deref(),&options);Ok(NewObject::new(MgResult{data,scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_geometric_fmg,"num_linear_geometric_fmg","crate::kernels::multigrid::GeometricMg::fmg","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<MgResult>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let data=a.data.fmg(input,max_cycles,&options);Ok(NewObject::new(MgResult{data,scale:scale(r.scale,1.0/a.scale)?},vec![d]))}))});
object_binding!(mg_geometric_cycle_once,"num_linear_geometric_cycle_once","crate::kernels::multigrid::GeometricMg::cycle_once","Native multigrid hierarchy in an explicit positive operator reference. RHS array references determine solution units RHS/operator. Initial iterates are normalized in that solution reference. Native cycle choices and convergence histories are preserved; setup/cycle/storage bounds are enforced.",true,(a:Handle<Geometric>,rhs:Handle<Array>,initial:Option<Vec<f64>>,cycle_name:String,smoother:String,pre:usize,post:usize,tolerance:f64,max_cycles:usize)->NewObject<Array>,|ctx,a,rhs,initial,cycle_name,smoother,pre,post,tolerance,max_cycles|{let _=ctx;let d=product(dim(ctx,&rhs.value)?,inverse(dim(ctx,&a.value)?)?)?;let options=controls(&cycle_name,&smoother,pre,post,tolerance,max_cycles)?;a.with(ctx,|a|rhs.with(ctx,|r|{let input=vector(r,a.data.size())?;if initial.as_ref().is_some_and(|v|v.len()!=input.len()){return Err(err("initial iterate length mismatch"));}let mut out=initial.unwrap_or_else(||vec![0.0;input.len()]);a.data.cycle_once(&mut out,input,&options);result(out,scale(r.scale,1.0/a.scale)?,d)}))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(i, BINDINGS, "numerical_matrix")
}
static BINDINGS: &[&api::Binding] = &[
    &mg_amg_build::BINDING,
    &mg_amg_general_build::BINDING,
    &mg_amg_num_levels::BINDING,
    &mg_amg_operator_complexity::BINDING,
    &mg_amg_precondition::BINDING,
    &mg_amg_vcycle::BINDING,
    &mg_amg_solve::BINDING,
    &mg_amg_general_num_levels::BINDING,
    &mg_amg_general_operator_complexity::BINDING,
    &mg_amg_general_precondition::BINDING,
    &mg_amg_general_precondition_many::BINDING,
    &mg_amg_general_workspace_len::BINDING,
    &mg_amg_general_vcycle::BINDING,
    &mg_amg_general_solve::BINDING,
    &mg_multigrid_mean_factor::BINDING,
    &mg_multigrid_result::BINDING,
    &mg_geometric_new::BINDING,
    &mg_geometric_levels::BINDING,
    &mg_geometric_size::BINDING,
    &mg_geometric_level_size::BINDING,
    &mg_geometric_restrict::BINDING,
    &mg_geometric_prolong::BINDING,
    &mg_geometric_apply::BINDING,
    &mg_geometric_residual_norm::BINDING,
    &mg_geometric_solve::BINDING,
    &mg_geometric_fmg::BINDING,
    &mg_geometric_cycle_once::BINDING,
];

#[cfg(test)]
mod transfer_tests {
    use super::*;

    #[test]
    fn geometric_workspace_retains_both_largest_transfer_buffers() {
        let mut hierarchy = Geometric {
            data: mg::GeometricMg::new(3, 31, 0.0).unwrap_or_else(|| panic!("expected hierarchy")),
            scale: 1.0,
            workspace: mg::GeometricTransferWorkspace::default(),
        };
        for level in (0..hierarchy.data.levels() - 1).rev() {
            let n =
                geometric_level_size(&hierarchy, level).unwrap_or_else(|error| panic!("{error}"));
            let values = Array {
                data: ArrayD::from_elem(IxDyn(&[n]), 1.0),
                scale: 1.0,
            };
            let coarse = geometric_transfer(
                &mut hierarchy,
                level,
                &values,
                DimExpr::Dimless,
                GridTransfer::Restrict,
            )
            .unwrap_or_else(|error| panic!("{error}"));
            geometric_transfer(
                &mut hierarchy,
                level,
                &coarse.value,
                DimExpr::Dimless,
                GridTransfer::Prolong,
            )
            .unwrap_or_else(|error| panic!("{error}"));
        }
        assert!(hierarchy.workspace.heap_bytes() >= 2 * 31 * 31 * 15 * size_of::<f64>());
        assert!(hierarchy.workspace_fits());
    }
}
