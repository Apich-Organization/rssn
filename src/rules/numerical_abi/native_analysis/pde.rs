//! Native PDE objects with explicit length/time/value reference frames.
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernels::pde as p;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::rules::numerical_abi::native_backend::CompiledSnapshot;
use crate::rules::numerical_abi::native_backend::compiled_snapshot;
use crate::rules::numerical_abi::native_object_api::Dynamic;
use crate::rules::numerical_abi::native_object_api::Handle;
use crate::rules::numerical_abi::native_object_api::Input;
use crate::rules::numerical_abi::native_object_api::NativeArg;
use crate::rules::numerical_abi::native_object_api::NativeType;
use crate::rules::numerical_abi::native_object_api::NewObject;
use crate::rules::numerical_abi::native_object_api::object_binding;
use crate::rules::numerical_abi::native_object_api::{
    self as api,
};
#[derive(Clone)]
struct Frame {
    length: f64,
    time: f64,
    value: f64,
    dims: Vec<DimExpr>,
}
impl Frame {
    fn new(
        l: &Input,
        t: &Input,
        v: &Input,
    ) -> Result<Self, OpError> {
        let (length, dl) = api::reference(&l.0)?;
        let (time, dt) = api::reference(&t.0)?;
        let (value, dv) = api::reference(&v.0)?;
        Ok(Self {
            length,
            time,
            value,
            dims: vec![dl, dt, dv],
        })
    }

    fn dims(&self) -> Vec<DimExpr> {
        self.dims.clone()
    }

    fn check(&self) -> Result<(), OpError> {
        if [self.length, self.time, self.value]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || self.dims.len() != 3
        {
            return Err(api::error("invalid PDE reference frame"));
        }
        Ok(())
    }

    fn equal(
        &self,
        other: &Self,
    ) -> Result<(), OpError> {
        for (a, b) in self.dims.iter().zip(&other.dims) {
            api::same_dimension(a, b)?;
        }
        if ![self.length, self.time, self.value]
            .into_iter()
            .zip([other.length, other.time, other.value])
            .all(|(a, b)| crate::kernels::exact_eq(a, b))
        {
            return Err(api::error(
                "PDE grid and boundary must use the same reference magnitudes",
            ));
        }
        Ok(())
    }

    fn unit(
        &self,
        length: i64,
        time: i64,
        value: i64,
    ) -> DimExpr {
        DimExpr::Mul(
            self.dims
                .iter()
                .zip([length, time, value])
                .map(|(d, p)| DimExpr::Pow(Box::new(d.clone()), Exp::Rat(Rational::new(p, 1))))
                .collect(),
        )
    }

    fn measured(
        &self,
        v: &Input,
        l: i64,
        t: i64,
        u: i64,
    ) -> Result<f64, OpError> {
        api::same_dimension(&self.unit(l, t, u), &api::dimension(&v.0)?)?;
        let x = super::numeric(&v.0, 0)?;
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        for (scale, exponent) in [(self.length, l), (self.time, t), (self.value, u)] {
            for _ in 0..exponent.unsigned_abs() {
                if exponent < 0 {
                    numerator.push(scale);
                } else {
                    denominator.push(scale);
                }
            }
        }
        let x = api::scale_product(x, &numerator, &denominator)?;
        if !x.is_finite() {
            return Err(api::error("PDE normalization is not finite"));
        }
        Ok(x)
    }
}
#[derive(Clone, Copy)]
enum Grid {
    D1(p::Grid1D),
    D2(p::Grid2D),
    D3(p::Grid3D),
}
#[derive(Clone)]
struct GridObject {
    native: Grid,
    frame: Frame,
}
impl GridObject {
    const fn rank(&self) -> usize {
        match self.native {
            | Grid::D1(_) => 1,
            | Grid::D2(_) => 2,
            | Grid::D3(_) => 3,
        }
    }

    const fn count(&self) -> usize {
        match self.native {
            | Grid::D1(g) => g.n,
            | Grid::D2(g) => g.nx * g.ny,
            | Grid::D3(g) => g.nx * g.ny * g.nz,
        }
    }
}
impl NativeType for GridObject {
    const KIND: &'static str = "rssn/native/analysis/pde/grid/v1";
    const MAX_BYTES: usize = 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.check()?;
        if self.count() > 4096 || self.count() < 3 {
            return Err(api::error("PDE grid must have3..4096 total nodes"));
        }
        Ok(())
    }
}
#[derive(Clone)]
struct BoundaryObject {
    native: p::Boundary,
    frame: Frame,
    snapshot: Option<CompiledSnapshot>,
}
impl NativeType for BoundaryObject {
    const KIND: &'static str = "rssn/native/analysis/pde/boundary/v1";
    const MAX_BYTES: usize = 2048;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.check()
    }
}
struct FieldObject {
    native: p::PdeSolution,
    frame: Frame,
    work: Option<(usize, usize)>,
}
impl NativeType for FieldObject {
    const KIND: &'static str = "rssn/native/analysis/pde/field/v1";
    const MAX_BYTES: usize = 8 * 1024 * 1024;

    fn validate(&self) -> Result<(), OpError> {
        self.frame.check()?;
        let values = match &self.native {
            | p::PdeSolution::D1(v) => &v.u,
            | p::PdeSolution::D2(v) => &v.u,
            | p::PdeSolution::D3(v) => &v.u,
        };
        if values.len() > 4096 || values.iter().any(|v| !v.is_finite()) {
            return Err(api::error(
                "PDE field storage is nonfinite or exceeds grid budget",
            ));
        }
        Ok(())
    }
}
fn vector(
    values: Vec<f64>,
    d: DimExpr,
) -> Result<Value, OpError> {
    if values.iter().any(|x| !x.is_finite()) {
        return Err(api::error("physical PDE output is not finite"));
    }
    Ok(Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), d),
        data: values.into_iter().map(Value::Float).collect(),
    }))
}
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn program(
    ctx: &CallCtx<'_>,
    handle: &Input,
    frame: &Frame,
    rank: usize,
    out: DimExpr,
    budget: usize,
) -> Result<CompiledSnapshot, OpError> {
    let f = compiled_snapshot(ctx, &handle.0, budget)?;
    if f.input_dimensions.len() != rank {
        return Err(api::error(
            "compiled PDE function arity differs from grid rank",
        ));
    }
    for d in &f.input_dimensions {
        api::same_dimension(&frame.dims[0], d)?;
    }
    api::same_dimension(&out, &f.output_dimension)?;
    Ok(f)
}
object_binding!(grid_new,"num_analysis_pde_grid","","Construct an actual native1D/2D/3D rectangular grid from physical bounds and explicit positive length/time/value references. Shape product is bounded before native allocations.",false,(lower:Input,upper:Input,shape:Vec<usize>,length:Input,time:Input,value:Input)->NewObject<GridObject>,|_ctx:&CallCtx<'_>,lower:Input,upper:Input,shape:Vec<usize>,length:Input,time:Input,value:Input|{
 let frame=Frame::new(&length,&time,&value)?;api::same_dimension(&frame.dims[0],&api::dimension(&lower.0)?)?;api::same_dimension(&frame.dims[0],&api::dimension(&upper.0)?)?;let lo=<Vec<f64> as Arg>::decode(&lower.0,0)?;let hi=<Vec<f64> as Arg>::decode(&upper.0,1)?;
 if !(1..=3).contains(&shape.len())||lo.len()!=shape.len()||hi.len()!=shape.len()||shape.iter().any(|n|*n<3)||shape.iter().try_fold(1usize,|a,b|a.checked_mul(*b)).is_none_or(|n|n>4096)||lo.iter().zip(&hi).any(|(a,b)|a>=b){return Err(api::error("invalid rectangular PDE bounds/shape or grid allocation budget"));}
 let lo=lo.into_iter().map(|x|x/frame.length).collect::<Vec<_>>();let hi=hi.into_iter().map(|x|x/frame.length).collect::<Vec<_>>();let native=match shape.len(){1=>Grid::D1(p::Grid1D{x0:lo[0],x1:hi[0],n:shape[0]}),2=>Grid::D2(p::Grid2D{x0:lo[0],x1:hi[0],y0:lo[1],y1:hi[1],nx:shape[0],ny:shape[1]}),_=>Grid::D3(p::Grid3D{x0:lo[0],x1:hi[0],y0:lo[1],y1:hi[1],z0:lo[2],z1:hi[2],nx:shape[0],ny:shape[1],nz:shape[2]})};Ok(NewObject::new(GridObject{native,frame:frame.clone()},frame.dims()))
});
object_binding!(grid_square,"num_analysis_pde_grid_unit_square","crate::kernels::pde::Grid2D::unit_square","Native unit-square grid interpreted through explicit physical references.",false,(n:usize,length:Input,time:Input,value:Input)->NewObject<GridObject>,|_ctx:&CallCtx<'_>,n:usize,length:Input,time:Input,value:Input|{if !(3..=64).contains(&n){return Err(api::error("unit-square side needs3..64 nodes"));}let frame=Frame::new(&length,&time,&value)?;Ok(NewObject::new(GridObject{native:Grid::D2(p::Grid2D::unit_square(n)),frame:frame.clone()},frame.dims()))});
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn boundary_constant(
    grid: GridObject,
    value: Input,
    neumann: bool,
) -> Result<NewObject<BoundaryObject>, OpError> {
    let number = grid
        .frame
        .measured(&value, if neumann { -1 } else { 0 }, 0, 1)?;
    let native = if neumann {
        p::Boundary::neumann_const(number)
    } else {
        p::Boundary::dirichlet_const(number)
    };
    Ok(NewObject::new(
        BoundaryObject {
            native,
            frame: grid.frame.clone(),
            snapshot: None,
        },
        grid.frame.dims(),
    ))
}
object_binding!(dirichlet_const,"num_analysis_pde_dirichlet_const","crate::kernels::pde::Boundary::dirichlet_const","Constant measured Dirichlet data in the retained grid reference frame.",false,(grid:Handle<GridObject>,value:Input)->NewObject<BoundaryObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,value:Input|boundary_constant(grid.with(ctx,|g|Ok(g.clone()))?,value,false));
object_binding!(neumann_const,"num_analysis_pde_neumann_const","crate::kernels::pde::Boundary::neumann_const","Constant measured outward normal gradient has field/length dimension.",false,(grid:Handle<GridObject>,value:Input)->NewObject<BoundaryObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,value:Input|boundary_constant(grid.with(ctx,|g|Ok(g.clone()))?,value,true));
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn boundary_program(
    ctx: &CallCtx<'_>,
    grid: Handle<GridObject>,
    function: Input,
    budget: usize,
    neumann: bool,
) -> Result<NewObject<BoundaryObject>, OpError> {
    let grid = grid.with(ctx, |g| Ok(g.clone()))?;
    let frame = grid.frame;
    let f = program(
        ctx,
        &function,
        &frame,
        2,
        frame.unit(if neumann { -1 } else { 0 }, 0, 1),
        budget,
    )?;
    let retained = f.clone();
    let l = frame.length;
    let scale = if neumann {
        frame.value / frame.length
    } else {
        frame.value
    };
    let sample = move |x: f64, y: f64| retained.sample(&[x * l, y * l]) / scale;
    let native = if neumann {
        p::Boundary::neumann(sample)
    } else {
        p::Boundary::dirichlet(sample)
    };
    Ok(NewObject::dependent(
        BoundaryObject {
            native,
            frame: frame.clone(),
            snapshot: Some(f),
        },
        frame.dims(),
        vec![function.0],
    ))
}
object_binding!(dirichlet,"num_analysis_pde_dirichlet","crate::kernels::pde::Boundary::dirichlet","Retain an immutable compiled two-coordinate physical boundary program, with lifetime evaluation budget and dependency-tracked program storage.",false,(grid:Handle<GridObject>,function:Input,budget:usize)->NewObject<BoundaryObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,function:Input,budget:usize|boundary_program(ctx,grid,function,budget,false));
object_binding!(neumann,"num_analysis_pde_neumann","crate::kernels::pde::Boundary::neumann","Retain an immutable compiled outward-gradient boundary program; field/length output and physical tangential coordinates are checked.",false,(grid:Handle<GridObject>,function:Input,budget:usize)->NewObject<BoundaryObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,function:Input,budget:usize|boundary_program(ctx,grid,function,budget,true));
fn boundary_list(
    ctx: &CallCtx<'_>,
    input: &Input,
    grid: &GridObject,
) -> Result<Vec<BoundaryObject>, OpError> {
    let Value::Tensor(t) = &input.0 else {
        return Err(api::error(
            "PDE boundaries require a vector of boundary handles",
        ));
    };
    if t.shape != [2 * grid.rank()] || t.data.len() != 2 * grid.rank() {
        return Err(api::error("PDE face count differs from grid dimension"));
    }
    t.data
        .iter()
        .map(|v| {
            let handle = <Handle<BoundaryObject> as NativeArg>::decode(v, 0, ctx)?;
            let value = handle.with(ctx, |b| Ok(b.clone()))?;
            grid.frame.equal(&value.frame)?;
            Ok(value)
        })
        .collect()
}
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // One normalized frame/validation path serves genuine native dimensional variants.
#[allow(
    clippy::needless_pass_by_value,
    reason = "The typed native ABI transfers decoded arguments by value; this shared dispatcher mirrors its explicit source signature."
)]
fn run(
    ctx: &CallCtx<'_>,
    kind: usize,
    grid: Handle<GridObject>,
    bc: Input,
    coefficient: Input,
    end: Input,
    steps: usize,
    method: String,
    relative: f64,
    budget: usize,
    primary: Input,
    secondary: Input,
    dispatch: bool,
) -> Result<NewObject<FieldObject>, OpError> {
    let grid = grid.with(ctx, |g| Ok(g.clone()))?;
    let rank = kind % 3 + 1;
    if grid.rank() != rank {
        return Err(api::error("PDE entry dimension differs from grid"));
    }
    let group = kind / 3;
    if group > 4 {
        return Err(api::error("unknown native PDE equation family"));
    }
    if steps == 0
        || steps > 1_000_000
        || steps
            .checked_mul(grid.count())
            .is_none_or(|n| n > 4_000_000)
    {
        return Err(api::error("PDE step times grid work budget exceeded"));
    }
    let frame = grid.frame.clone();
    let end = frame.measured(&end, 0, 1, 0)?;
    if !end.is_finite() || end < 0.0 {
        return Err(api::error("PDE end time must be nonnegative"));
    }
    let boundaries = boundary_list(ctx, &bc, &grid)?;
    let values = boundaries
        .iter()
        .map(|b| b.native.clone())
        .collect::<Vec<_>>();
    let b1: [p::Boundary; 2] = if rank == 1 {
        values
            .clone()
            .try_into()
            .map_err(|_| api::error("one-dimensional face count"))?
    } else {
        std::array::from_fn(|_| p::Boundary::dirichlet_const(0.0))
    };
    let b2: [p::Boundary; 4] = if rank == 2 {
        values
            .clone()
            .try_into()
            .map_err(|_| api::error("two-dimensional face count"))?
    } else {
        std::array::from_fn(|_| p::Boundary::dirichlet_const(0.0))
    };
    let b3: [p::Boundary; 6] = if rank == 3 {
        values
            .try_into()
            .map_err(|_| api::error("three-dimensional face count"))?
    } else {
        std::array::from_fn(|_| p::Boundary::dirichlet_const(0.0))
    };
    if !relative.is_finite() || relative <= 0.0 {
        return Err(api::error("relative elliptic tolerance must be positive"));
    }
    let solver = match method.as_str() {
        | "SparseLu" => p::PoissonSolver::SparseLu,
        | "Cg" => p::PoissonSolver::Cg { tol: relative },
        | "MeanFree" => p::PoissonSolver::MeanFree { tol: relative },
        | _ if group == 2 || group == 3 => {
            return Err(api::error("Poisson solver must be SparseLu,Cg,orMeanFree"));
        },
        | _ => p::PoissonSolver::SparseLu,
    };
    let scheme = match method.as_str() {
        | "LaxWendroff" => p::AdvectionScheme::LaxWendroff,
        | "Upwind" => p::AdvectionScheme::Upwind,
        | _ if group == 4 => {
            return Err(api::error("advection scheme must be Upwind or LaxWendroff"));
        },
        | _ => p::AdvectionScheme::Upwind,
    };
    let initial = if group == 3 {
        None
    } else {
        Some(program(
            ctx,
            &primary,
            &frame,
            rank,
            frame.unit(if group == 2 { -2 } else { 0 }, 0, 1),
            budget,
        )?)
    };
    let velocity = if group == 1 {
        Some(program(
            ctx,
            &secondary,
            &frame,
            rank,
            frame.unit(0, -1, 1),
            budget,
        )?)
    } else {
        None
    };
    let factor = if group == 2 {
        frame.value / (frame.length * frame.length)
    } else {
        frame.value
    };
    let sample = |x: &[f64]| {
        initial.as_ref().map_or(0.0, |f| {
            f.sample(&x.iter().map(|v| v * frame.length).collect::<Vec<_>>()) / factor
        })
    };
    let speed = |x: &[f64]| {
        velocity.as_ref().map_or(0.0, |f| {
            f.sample(&x.iter().map(|v| v * frame.length).collect::<Vec<_>>())
                / (frame.value / frame.time)
        })
    };
    let f1 = |x| sample(&[x]);
    let f2 = |x, y| sample(&[x, y]);
    let f3 = |x, y, z| sample(&[x, y, z]);
    let v1 = |x| speed(&[x]);
    let v2 = |x, y| speed(&[x, y]);
    let v3 = |x, y, z| speed(&[x, y, z]);
    let c = if group <= 1 {
        frame.measured(&coefficient, if group == 0 { 2 } else { 1 }, -1, 0)?
    } else {
        0.0
    };
    let a = if group == 4 {
        api::same_dimension(&frame.unit(1, -1, 0), &api::dimension(&coefficient.0)?)?;
        let a = <Vec<f64> as Arg>::decode(&coefficient.0, 2)?;
        if a.len() != rank {
            return Err(api::error(
                "advection velocity count differs from grid rank",
            ));
        }
        a.into_iter()
            .map(|v| v / (frame.length / frame.time))
            .collect::<Vec<_>>()
    } else {
        vec![0.0; rank]
    };
    let native = match (group, grid.native, dispatch) {
        | (0, Grid::D1(g), false) => {
            p::heat_1d(c, &g, &b1, &f1, end, steps).map(p::PdeSolution::D1)
        },
        | (0, Grid::D2(g), false) => {
            p::heat_2d(c, &g, &b2, &f2, end, steps).map(p::PdeSolution::D2)
        },
        | (0, Grid::D3(g), false) => {
            p::heat_3d(c, &g, &b3, &f3, end, steps).map(p::PdeSolution::D3)
        },
        | (1, Grid::D1(g), false) => {
            p::wave_1d(c, &g, &b1, &f1, &v1, end, steps).map(p::PdeSolution::D1)
        },
        | (1, Grid::D2(g), false) => {
            p::wave_2d(c, &g, &b2, &f2, &v2, end, steps).map(p::PdeSolution::D2)
        },
        | (1, Grid::D3(g), false) => {
            p::wave_3d(c, &g, &b3, &f3, &v3, end, steps).map(p::PdeSolution::D3)
        },
        | (2, Grid::D1(g), false) => p::poisson_1d(&g, &b1, &f1, solver).map(p::PdeSolution::D1),
        | (2, Grid::D2(g), false) => p::poisson_2d(&g, &b2, &f2, solver).map(p::PdeSolution::D2),
        | (2, Grid::D3(g), false) => p::poisson_3d(&g, &b3, &f3, solver).map(p::PdeSolution::D3),
        | (3, Grid::D1(g), _) => p::laplace_1d(&g, &b1, solver).map(p::PdeSolution::D1),
        | (3, Grid::D2(g), _) => p::laplace_2d(&g, &b2, solver).map(p::PdeSolution::D2),
        | (3, Grid::D3(g), _) => p::laplace_3d(&g, &b3, solver).map(p::PdeSolution::D3),
        | (4, Grid::D1(g), false) => {
            p::advection_1d(a[0], &g, scheme, &f1, end, steps).map(p::PdeSolution::D1)
        },
        | (4, Grid::D2(g), false) => {
            p::advection_2d([a[0], a[1]], &g, scheme, &f2, end, steps).map(p::PdeSolution::D2)
        },
        | (4, Grid::D3(g), false) => {
            p::advection_3d([a[0], a[1], a[2]], &g, scheme, &f3, end, steps).map(p::PdeSolution::D3)
        },
        | (0, Grid::D1(g), true) => {
            p::pde_solver(&p::PdeProblem::Heat1D {
                alpha: c,
                grid: g,
                bc: b1,
                u0: &f1,
                t_end: end,
                steps,
            })
        },
        | (0, Grid::D2(g), true) => {
            p::pde_solver(&p::PdeProblem::Heat2D {
                alpha: c,
                grid: g,
                bc: b2,
                u0: &f2,
                t_end: end,
                steps,
            })
        },
        | (0, Grid::D3(g), true) => {
            p::pde_solver(&p::PdeProblem::Heat3D {
                alpha: c,
                grid: g,
                bc: b3,
                u0: &f3,
                t_end: end,
                steps,
            })
        },
        | (1, Grid::D1(g), true) => {
            p::pde_solver(&p::PdeProblem::Wave1D {
                c,
                grid: g,
                bc: b1,
                u0: &f1,
                v0: &v1,
                t_end: end,
                steps,
            })
        },
        | (1, Grid::D2(g), true) => {
            p::pde_solver(&p::PdeProblem::Wave2D {
                c,
                grid: g,
                bc: b2,
                u0: &f2,
                v0: &v2,
                t_end: end,
                steps,
            })
        },
        | (1, Grid::D3(g), true) => {
            p::pde_solver(&p::PdeProblem::Wave3D {
                c,
                grid: g,
                bc: b3,
                u0: &f3,
                v0: &v3,
                t_end: end,
                steps,
            })
        },
        | (2, Grid::D1(g), true) => {
            p::pde_solver(&p::PdeProblem::Poisson1D {
                grid: g,
                bc: b1,
                f: &f1,
                solver,
            })
        },
        | (2, Grid::D2(g), true) => {
            p::pde_solver(&p::PdeProblem::Poisson2D {
                grid: g,
                bc: b2,
                f: &f2,
                solver,
            })
        },
        | (2, Grid::D3(g), true) => {
            p::pde_solver(&p::PdeProblem::Poisson3D {
                grid: g,
                bc: b3,
                f: &f3,
                solver,
            })
        },
        | (4, Grid::D1(g), true) => {
            p::pde_solver(&p::PdeProblem::Advection1D {
                a: a[0],
                grid: g,
                scheme,
                u0: &f1,
                t_end: end,
                steps,
            })
        },
        | (4, Grid::D2(g), true) => {
            p::pde_solver(&p::PdeProblem::Advection2D {
                a: [a[0], a[1]],
                grid: g,
                scheme,
                u0: &f2,
                t_end: end,
                steps,
            })
        },
        | (4, Grid::D3(g), true) => {
            p::pde_solver(&p::PdeProblem::Advection3D {
                a: [a[0], a[1], a[2]],
                grid: g,
                scheme,
                u0: &f3,
                t_end: end,
                steps,
            })
        },
        | _ => return Err(api::error("invalid PDE dispatch")),
    };
    if let Some(f) = initial {
        f.finish()?;
    }
    if let Some(f) = velocity {
        f.finish()?;
    }
    for b in boundaries {
        if let Some(f) = b.snapshot {
            f.finish()?;
        }
    }
    let native = native.map_err(api::error)?;
    Ok(NewObject::new(
        FieldObject {
            native,
            frame: frame.clone(),
            work: None,
        },
        frame.dims(),
    ))
}

object_binding!(heat_1d,"num_analysis_pde_heat_1d","crate::kernels::pde::heat_1d","Native heat_1d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,0,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(heat_2d,"num_analysis_pde_heat_2d","crate::kernels::pde::heat_2d","Native heat_2d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,1,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(heat_3d,"num_analysis_pde_heat_3d","crate::kernels::pde::heat_3d","Native heat_3d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,2,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(wave_1d,"num_analysis_pde_wave_1d","crate::kernels::pde::wave_1d","Native wave_1d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,3,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(wave_2d,"num_analysis_pde_wave_2d","crate::kernels::pde::wave_2d","Native wave_2d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,4,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(wave_3d,"num_analysis_pde_wave_3d","crate::kernels::pde::wave_3d","Native wave_3d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,5,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(poisson_1d,"num_analysis_pde_poisson_1d","crate::kernels::pde::poisson_1d","Native poisson_1d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,6,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(poisson_2d,"num_analysis_pde_poisson_2d","crate::kernels::pde::poisson_2d","Native poisson_2d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,7,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(poisson_3d,"num_analysis_pde_poisson_3d","crate::kernels::pde::poisson_3d","Native poisson_3d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,8,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(laplace_1d,"num_analysis_pde_laplace_1d","crate::kernels::pde::laplace_1d","Native laplace_1d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,9,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(laplace_2d,"num_analysis_pde_laplace_2d","crate::kernels::pde::laplace_2d","Native laplace_2d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,10,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(laplace_3d,"num_analysis_pde_laplace_3d","crate::kernels::pde::laplace_3d","Native laplace_3d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,11,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(advection_1d,"num_analysis_pde_advection_1d","crate::kernels::pde::advection_1d","Native advection_1d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,12,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(advection_2d,"num_analysis_pde_advection_2d","crate::kernels::pde::advection_2d","Native advection_2d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,13,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(advection_3d,"num_analysis_pde_advection_3d","crate::kernels::pde::advection_3d","Native advection_3d in the grid physical reference frame. Boundary handles share that frame. Heat coefficient L^2/T,wave/advection L/T; primary compiled output V for initial data or V/L^2 for Poisson; secondary V/T only for wave. Method SparseLu/Cg/MeanFree for elliptic,Upwind/LaxWendroff for advection; relative>0. Dedicated Laplace ignores primary; other unused coefficients/secondary are explicit placeholder Values.",false,(grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|run(ctx,14,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,false));

object_binding!(dispatch,"num_analysis_pde_solver","crate::kernels::pde::pde_solver","Native PdeProblem dispatch. Kind0..2heat1D..3D,3..5wave,6..8Poisson,12..14advection; use dedicated Laplace endpoints for9..11. Reference and compiled-function contracts match dedicated native solvers.",false,(kind:usize,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input)->NewObject<FieldObject>,|ctx:&CallCtx<'_>,kind:usize,grid:Handle<GridObject>,boundaries:Input,coefficient:Input,end:Input,steps:usize,method:String,relative:f64,budget:usize,primary:Input,secondary:Input|{if (9..=11).contains(&kind){return Err(api::error("PdeProblem has no Laplace variant; use dedicated Laplace APIs"));}run(ctx,kind,grid,boundaries,coefficient,end,steps,method,relative,budget,primary,secondary,true)});

object_binding!(grid1_dx,"num_analysis_pde_grid1_dx","crate::kernels::pde::Grid1D::dx","Physical native grid spacing in the retained length domain.",false,(grid:Handle<GridObject>)->Dynamic,|ctx:&CallCtx<'_>,grid:Handle<GridObject>|grid.with(ctx,|grid|{let Grid::D1(g)=grid.native else{return Err(api::error("grid rank mismatch"));};Ok(Dynamic(api::quantity(g.dx()*grid.frame.length,grid.frame.dims[0].clone())?))}));

object_binding!(grid2_dx,"num_analysis_pde_grid2_dx","crate::kernels::pde::Grid2D::dx","Physical native grid spacing in the retained length domain.",false,(grid:Handle<GridObject>)->Dynamic,|ctx:&CallCtx<'_>,grid:Handle<GridObject>|grid.with(ctx,|grid|{let Grid::D2(g)=grid.native else{return Err(api::error("grid rank mismatch"));};Ok(Dynamic(api::quantity(g.dx()*grid.frame.length,grid.frame.dims[0].clone())?))}));

object_binding!(grid2_dy,"num_analysis_pde_grid2_dy","crate::kernels::pde::Grid2D::dy","Physical native grid spacing in the retained length domain.",false,(grid:Handle<GridObject>)->Dynamic,|ctx:&CallCtx<'_>,grid:Handle<GridObject>|grid.with(ctx,|grid|{let Grid::D2(g)=grid.native else{return Err(api::error("grid rank mismatch"));};Ok(Dynamic(api::quantity(g.dy()*grid.frame.length,grid.frame.dims[0].clone())?))}));

object_binding!(grid1_nodes,"num_analysis_pde_grid1_nodes","crate::kernels::pde::Grid1D::nodes","Native1D grid nodes restored to measured coordinates.",false,(grid:Handle<GridObject>)->Dynamic,|ctx:&CallCtx<'_>,grid:Handle<GridObject>|grid.with(ctx,|grid|{let Grid::D1(g)=grid.native else{return Err(api::error("grid needs rank1"));};Ok(Dynamic(vector(g.nodes().into_iter().map(|x|x*grid.frame.length).collect(),grid.frame.dims[0].clone())?))}));
object_binding!(field2_at,"num_analysis_pde_field2_at","crate::kernels::pde::Field2D::at","Checked native2D field sample retains measured ordinate units.",false,(field:Handle<FieldObject>,i:usize,j:usize)->Dynamic,|ctx:&CallCtx<'_>,field:Handle<FieldObject>,i:usize,j:usize|field.with(ctx,|field|{let p::PdeSolution::D2(f)=&field.native else{return Err(api::error("field needs rank2"));};
if i>=f.nx||j>=f.ny{return Err(api::error("field index out of bounds"));}Ok(Dynamic(api::quantity(f.at(i,j)*field.frame.value,field.frame.dims[2].clone())?))}));
object_binding!(field2_point,"num_analysis_pde_field2_point","crate::kernels::pde::Field2D::point","Checked native2D field node coordinates retain measured length units.",false,(field:Handle<FieldObject>,i:usize,j:usize)->Dynamic,|ctx:&CallCtx<'_>,field:Handle<FieldObject>,i:usize,j:usize|field.with(ctx,|field|{let p::PdeSolution::D2(f)=&field.native else{return Err(api::error("field needs rank2"));};
if i>=f.nx||j>=f.ny{return Err(api::error("field index out of bounds"));}let(x,y)=f.point(i,j);Ok(Dynamic(api::record(vec![("item0",api::quantity(x*field.frame.length,field.frame.dims[0].clone())?),("item1",api::quantity(y*field.frame.length,field.frame.dims[0].clone())?)])))}));
object_binding!(field_data,"num_analysis_pde_field_data","","Retrieve native field samples,shape,and physical snapshot time with preserved measurement units.",false,(field:Handle<FieldObject>)->Dynamic,|ctx:&CallCtx<'_>,field:Handle<FieldObject>|field.with(ctx,|field|{let(values,shape,t)=match &field.native{p::PdeSolution::D1(f)=>(f.u.clone(),vec![f.u.len()],f.t),p::PdeSolution::D2(f)=>(f.u.clone(),vec![f.ny,f.nx],f.t),p::PdeSolution::D3(f)=>(f.u.clone(),vec![f.nz,f.ny,f.nx],f.t)};Ok(Dynamic(api::record(vec![("values",vector(values.into_iter().map(|v|v*field.frame.value).collect(),field.frame.dims[2].clone())?),("shape",Output::encode(shape)?),("accepted_steps",Output::encode(field.work.map(|w|w.0))?),("evaluations",Output::encode(field.work.map(|w|w.1))?),("time",api::quantity(t*field.frame.time,field.frame.dims[1].clone())?)])))}));

object_binding!(field1_max_error,"num_analysis_pde_field1_max_error","crate::kernels::pde::Field1D::max_error","Native maximum nodal error against an owned physical exact-solution program; returns measured ordinate error.",false,(field:Handle<FieldObject>,exact:Input,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,field:Handle<FieldObject>,exact:Input,budget:usize|field.with(ctx,|field|{let p::PdeSolution::D1(f)=&field.native else{return Err(api::error("field dimension mismatch"));};let exact=program(ctx,&exact,&field.frame,1,field.frame.dims[2].clone(),budget)?;let l=field.frame.length;let v=field.frame.value;let error=f.max_error(|x|exact.sample(&[x*l])/v);exact.finish()?;Ok(Dynamic(api::quantity(error*v,field.frame.dims[2].clone())?))}));

object_binding!(field2_max_error,"num_analysis_pde_field2_max_error","crate::kernels::pde::Field2D::max_error","Native maximum nodal error against an owned physical exact-solution program; returns measured ordinate error.",false,(field:Handle<FieldObject>,exact:Input,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,field:Handle<FieldObject>,exact:Input,budget:usize|field.with(ctx,|field|{let p::PdeSolution::D2(f)=&field.native else{return Err(api::error("field dimension mismatch"));};let exact=program(ctx,&exact,&field.frame,2,field.frame.dims[2].clone(),budget)?;let l=field.frame.length;let v=field.frame.value;let error=f.max_error(|x,y|exact.sample(&[x*l,y*l])/v);exact.finish()?;Ok(Dynamic(api::quantity(error*v,field.frame.dims[2].clone())?))}));
#[path = "pde/mol.rs"]
mod mol;

object_binding!(grid3_dx,"num_analysis_pde_grid3_dx","crate::kernels::pde::Grid3D::dx","Actual native3D grid method; physical spacings retain the grid length reference.",false,(grid:Handle<GridObject>)->Dynamic,|ctx,grid|grid.with(ctx,|grid|{let Grid::D3(g)=grid.native else{return Err(api::error("grid needs rank3"));};Ok(Dynamic(api::quantity(api::scale_product(g.dx(),&[grid.frame.length],&[])?,grid.frame.dims[0].clone())?))}));
object_binding!(grid3_dy,"num_analysis_pde_grid3_dy","crate::kernels::pde::Grid3D::dy","Actual native3D grid method; physical spacings retain the grid length reference.",false,(grid:Handle<GridObject>)->Dynamic,|ctx,grid|grid.with(ctx,|grid|{let Grid::D3(g)=grid.native else{return Err(api::error("grid needs rank3"));};Ok(Dynamic(api::quantity(api::scale_product(g.dy(),&[grid.frame.length],&[])?,grid.frame.dims[0].clone())?))}));
object_binding!(grid3_dz,"num_analysis_pde_grid3_dz","crate::kernels::pde::Grid3D::dz","Actual native3D grid method; physical spacings retain the grid length reference.",false,(grid:Handle<GridObject>)->Dynamic,|ctx,grid|grid.with(ctx,|grid|{let Grid::D3(g)=grid.native else{return Err(api::error("grid needs rank3"));};Ok(Dynamic(api::quantity(api::scale_product(g.dz(),&[grid.frame.length],&[])?,grid.frame.dims[0].clone())?))}));
object_binding!(grid3_len,"num_analysis_pde_grid3_len","crate::kernels::pde::Grid3D::len","Actual native3D grid method; physical spacings retain the grid length reference.",false,(grid:Handle<GridObject>)->usize,|ctx,grid|grid.with(ctx,|grid|{let Grid::D3(g)=grid.native else{return Err(api::error("grid needs rank3"));};Ok(g.len())}));
object_binding!(grid3_is_empty,"num_analysis_pde_grid3_is_empty","crate::kernels::pde::Grid3D::is_empty","Actual native3D grid method; physical spacings retain the grid length reference.",false,(grid:Handle<GridObject>)->bool,|ctx,grid|grid.with(ctx,|grid|{let Grid::D3(g)=grid.native else{return Err(api::error("grid needs rank3"));};Ok(g.is_empty())}));
object_binding!(grid_cube,"num_analysis_pde_grid_unit_cube","crate::kernels::pde::Grid3D::unit_cube","Native unit-cube grid interpreted through explicit physical references, bounded before native allocation.",false,(n:usize,length:Input,time:Input,value:Input)->NewObject<GridObject>,|_ctx,n,length,time,value|{if !(3..=16).contains(&n){return Err(api::error("unit-cube side needs3..16nodes"));}let frame=Frame::new(&length,&time,&value)?;Ok(NewObject::new(GridObject{native:Grid::D3(p::Grid3D::unit_cube(n)),frame:frame.clone()},frame.dims()))});
object_binding!(field3_at,"num_analysis_pde_field3_at","crate::kernels::pde::Field3D::at","Checked native3D field sample with physical field units.",false,(field:Handle<FieldObject>,i:usize,j:usize,k:usize)->Dynamic,|ctx,field,i,j,k|field.with(ctx,|field|{let p::PdeSolution::D3(f)=&field.native else{return Err(api::error("field needs rank3"));};
if i>=f.nx||j>=f.ny||k>=f.nz{return Err(api::error("field3 index out of range"));}Ok(Dynamic(api::quantity(api::scale_product(f.at(i,j,k),&[field.frame.value],&[])?,field.frame.dims[2].clone())?))}));
object_binding!(field3_point,"num_analysis_pde_field3_point","crate::kernels::pde::Field3D::point","Checked actual native3D node coordinates in physical length units.",false,(field:Handle<FieldObject>,i:usize,j:usize,k:usize)->Dynamic,|ctx,field,i,j,k|field.with(ctx,|field|{let p::PdeSolution::D3(f)=&field.native else{return Err(api::error("field needs rank3"));};
if i>=f.nx||j>=f.ny||k>=f.nz{return Err(api::error("field3 index out of range"));}let(x,y,z)=f.point(i,j,k);Ok(Dynamic(vector([x,y,z].into_iter().map(|v|api::scale_product(v,&[field.frame.length],&[])).collect::<Result<_,_>>()?,field.frame.dims[0].clone())?))}));
object_binding!(field3_max_error,"num_analysis_pde_field3_max_error","crate::kernels::pde::Field3D::max_error","Actual native3D maximum nodal error against an immutable physical exact-solution program, with original evaluation failures preserved.",false,(field:Handle<FieldObject>,exact:Input,budget:usize)->Dynamic,|ctx,field,exact,budget|field.with(ctx,|field|{let p::PdeSolution::D3(f)=&field.native else{return Err(api::error("field needs rank3"));};let exact=program(ctx,&exact,&field.frame,3,field.frame.dims[2].clone(),budget)?;let l=field.frame.length;let v=field.frame.value;let error=f.max_error(|x,y,z|{let physical=[x,y,z].map(|x|api::scale_product(x,&[l],&[]).unwrap_or(f64::NAN));api::scale_product(exact.sample(&physical),&[],&[v]).unwrap_or(f64::NAN)});exact.finish()?;Ok(Dynamic(api::quantity(api::scale_product(error,&[v],&[])?,field.frame.dims[2].clone())?))}));
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &grid3_dx::BINDING,
            &grid3_dy::BINDING,
            &grid3_dz::BINDING,
            &grid3_len::BINDING,
            &grid3_is_empty::BINDING,
            &grid_cube::BINDING,
            &field3_at::BINDING,
            &field3_point::BINDING,
            &field3_max_error::BINDING,
            &grid_new::BINDING,
            &grid_square::BINDING,
            &dirichlet_const::BINDING,
            &neumann_const::BINDING,
            &dirichlet::BINDING,
            &neumann::BINDING,
            &heat_1d::BINDING,
            &heat_2d::BINDING,
            &heat_3d::BINDING,
            &wave_1d::BINDING,
            &wave_2d::BINDING,
            &wave_3d::BINDING,
            &poisson_1d::BINDING,
            &poisson_2d::BINDING,
            &poisson_3d::BINDING,
            &laplace_1d::BINDING,
            &laplace_2d::BINDING,
            &laplace_3d::BINDING,
            &advection_1d::BINDING,
            &advection_2d::BINDING,
            &advection_3d::BINDING,
            &dispatch::BINDING,
            &grid1_dx::BINDING,
            &grid2_dx::BINDING,
            &grid2_dy::BINDING,
            &grid1_nodes::BINDING,
            &field2_at::BINDING,
            &field2_point::BINDING,
            &field_data::BINDING,
            &field1_max_error::BINDING,
            &field2_max_error::BINDING,
            mol::FIXED,
            mol::TIMED,
        ],
        "numerical_vector",
    )
}
