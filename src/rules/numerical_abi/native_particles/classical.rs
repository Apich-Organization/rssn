#![allow(clippy::wildcard_imports)] // Native classical entry points with coherent SI contracts and bounded typed callbacks.
use super::engine;
use super::*;
use crate::rules::numerical_abi::codecs::Arg;
use crate::sim::classical as c;
use std::cell::Cell;
use std::cell::RefCell;
const HBAR: f64 = 1.054_571_817e-34;
pub(super) fn si() -> Frame {
    Frame {
        values: [1.0; 4],
        dims: [
            DimExpr::base("Mass"),
            DimExpr::base("Length"),
            DimExpr::base("Time"),
            DimExpr::Mul(vec![DimExpr::base("Current"), DimExpr::base("Time")]),
        ],
        physical: true,
    }
}
impl NativeType for c::Particle3D {
    const KIND: &'static str = "rssn/native/sim/classical-particle-si/v1";
    const MAX_BYTES: usize = 128;

    fn validate(&self) -> Result<(), OpError> {
        if self.mass <= 0.0 {
            return Err(api::error("classical particle mass must bepositive"));
        }
        engine::finite([self.mass, self.x, self.y, self.z, self.vx, self.vy, self.vz])
    }
}
struct Particles(Vec<c::Particle3D>);
impl NativeArg for Particles {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Opaque(c::Particle3D::KIND.into()))
    }

    fn decode(
        v: &Value,
        _: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = v else {
            return Err(api::error("particle handle vector expected"));
        };
        if t.shape != [t.data.len()] || t.data.is_empty() || t.data.len() > 512 {
            return Err(api::error(
                "classical particles neednonempty rank1<=512handles",
            ));
        }
        Ok(Self(
            t.data
                .iter()
                .map(|v| {
                    let h = <Handle<c::Particle3D> as NativeArg>::decode(v, 0, ctx)?;
                    h.with(ctx, |p| Ok(*p))
                })
                .collect::<Result<_, _>>()?,
        ))
    }
}
fn particle_data(p: &c::Particle3D) -> Result<Value, OpError> {
    let f = si();
    Ok(api::record(vec![
        ("mass", f.output(p.mass, MASS, false)?),
        ("position", f.array(vec![p.x, p.y, p.z], LENGTH)?),
        ("velocity", f.array(vec![p.vx, p.vy, p.vz], VELOCITY)?),
    ]))
}
fn history(
    data: &[Vec<f64>],
    f: &Frame,
) -> Result<Dynamic, OpError> {
    if data.iter().any(|r| r.len() != 6) {
        return Err(api::error("native particle history haswrongstatesize"));
    }
    let pos = data.iter().map(|r| [r[0], r[1], r[2]]).collect();
    let vel = data.iter().map(|r| [r[3], r[4], r[5]]).collect();
    Ok(Dynamic(api::record(vec![
        ("positions", engine::packed(f, pos, LENGTH)?),
        ("velocities", engine::packed(f, vel, VELOCITY)?),
    ])))
}
struct Stop;
struct Gate<'a, 'b> {
    ctx: &'a CallCtx<'b>,
    left: Cell<usize>,
    failure: RefCell<Option<OpError>>,
}
impl<'a, 'b> Gate<'a, 'b> {
    fn new(
        ctx: &'a CallCtx<'b>,
        budget: usize,
    ) -> Result<Self, OpError> {
        if budget == 0 || budget > 1_000_000 {
            return Err(api::error("classical callback budget must be1..1000000"));
        }
        Ok(Self {
            ctx,
            left: Cell::new(budget),
            failure: RefCell::new(None),
        })
    }

    fn call<T>(
        &self,
        work: impl FnOnce() -> Result<T, OpError>,
    ) -> T {
        let result = self.left.get().checked_sub(1).map_or_else(
            || Err(api::error("native classical callback budget exhausted")),
            |left| {
                self.left.set(left);
                work()
            },
        );
        match result {
            | Ok(v) => v,
            | Err(e) => {
                *self.failure.borrow_mut() = Some(e);
                std::panic::resume_unwind(Box::new(Stop))
            },
        }
    }

    fn run<T>(
        &self,
        work: impl FnOnce() -> Result<T, String>,
    ) -> Result<T, OpError> {
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(work));
        if let Some(e) = self.failure.borrow_mut().take() {
            return Err(e);
        }
        match result {
            | Ok(v) => v.map_err(api::error),
            | Err(e) => std::panic::resume_unwind(e),
        }
    }
}
struct ForceFn(CallbackRef);
struct Potential<const D: usize>(CallbackRef);
struct Initial<const VECTOR: bool> {
    handle: CallbackRef,
    dimension: DimExpr,
}
fn scalar_dimension(t: &TypeExpr) -> Result<DimExpr, OpError> {
    match t {
        | TypeExpr::Float => Ok(DimExpr::Dimless),
        | TypeExpr::Dim(t, d) if **t == TypeExpr::Float => Ok(d.clone()),
        | _ => Err(api::error("initial callback result must beconcreteFloat")),
    }
}
fn callback(
    v: &Value,
    c: &CallCtx<'_>,
    arity: usize,
) -> Result<CallbackRef, OpError> {
    let Value::Function(h) = v else {
        return Err(api::error("typed callback expected"));
    };
    if c.callback_signature(*h)?.params.len() != arity {
        return Err(api::error("callback arity differsfromnativefunction"));
    }
    Ok(*h)
}
impl NativeArg for ForceFn {
    fn ty() -> TypeExpr {
        let f = si();
        TypeExpr::Function(
            vec![
                TypeExpr::Float.with_dim(f.dim(TIME).expect("SI dimension")),
                TypeExpr::Record(vec![
                    (
                        "position".into(),
                        api::tensor_type(
                            1,
                            TypeExpr::Float.with_dim(f.dim(LENGTH).expect("SI dimension")),
                        ),
                    ),
                    (
                        "velocity".into(),
                        api::tensor_type(
                            1,
                            TypeExpr::Float.with_dim(f.dim(VELOCITY).expect("SI dimension")),
                        ),
                    ),
                ]),
            ],
            Box::new(api::tensor_type(
                1,
                TypeExpr::Float.with_dim(f.dim(FORCE).expect("SI dimension")),
            )),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, 2)?))
    }
}
impl<const D: usize> NativeArg for Potential<D> {
    fn ty() -> TypeExpr {
        TypeExpr::Function(
            vec![TypeExpr::Float.with_dim(DimExpr::base("Length")); D],
            Box::new(TypeExpr::Float.with_dim(si().dim(ENERGY).expect("SI dimension"))),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, D)?))
    }
}
impl<const VECTOR: bool> NativeArg for Initial<VECTOR> {
    fn ty() -> TypeExpr {
        let x = TypeExpr::Float.with_dim(DimExpr::base("Length"));
        TypeExpr::Function(
            vec![if VECTOR {
                api::tensor_type(1, x)
            } else {
                x
            }],
            Box::new(TypeExpr::Float.with_dim(DimExpr::Var("U".into()))),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let handle = callback(v, c, 1)?;
        let dimension = scalar_dimension(&c.callback_signature(handle)?.ret)?;
        Ok(Self { handle, dimension })
    }
}
fn initial_value<const V: bool>(
    g: &Gate<'_, '_>,
    f: &Initial<V>,
    x: &[f64],
) -> f64 {
    g.call(|| {
        let x = if V {
            si().array(x.to_vec(), LENGTH)?
        } else {
            si().output(x[0], LENGTH, false)?
        };
        let out = g.ctx.invoke_callback(f.handle, &[x])?;
        api::same_dimension(&api::dimension(&out)?, &f.dimension)?;
        crate::rules::numerical_abi::numeric(&out, 0)
    })
}
fn measured_history(
    rows: Vec<Vec<f64>>,
    d: DimExpr,
) -> Result<Dynamic, OpError> {
    let n = rows.len();
    let m = rows.first().map_or(0, Vec::len);
    if rows.iter().any(|r| r.len() != m) {
        return Err(api::error("history outputragged"));
    }
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![n, m],
        data: rows
            .into_iter()
            .flatten()
            .map(|v| api::quantity(v, d.clone()))
            .collect::<Result<_, _>>()?,
        elem: TypeExpr::Float.with_dim(d),
    })))
}
fn grid_work(
    nodes: &[usize],
    steps: usize,
) -> Result<usize, OpError> {
    let n = nodes
        .iter()
        .try_fold(1usize, |a, n| a.checked_mul(*n))
        .filter(|n| *n <= 262144)
        .ok_or_else(|| api::error("classical grid allocationbound exceeded"))?;
    if nodes.iter().any(|n| *n < 3)
        || steps == 0
        || steps > 100000
        || n.checked_mul(steps + 1).is_none_or(|n| n > 1_048_576)
    {
        return Err(api::error("classical grid/history work bound exceeded"));
    }
    Ok(n)
}
object_binding!(classical_particle_new,"num_native_classical_particle_new","crate::sim::classical::Particle3D::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(mass:Input,position:Input,velocity:Input)->NewObject<c::Particle3D>,|_ctx:&CallCtx<'_>,mass:Input,position:Input,velocity:Input|{let f=si();let m=f.scalar(&mass.0,MASS,false)?;let p=engine::array::<3>(&f,&position.0,LENGTH)?;let v=engine::array::<3>(&f,&velocity.0,VELOCITY)?;Ok(NewObject::new(c::Particle3D::new(m,p[0],p[1],p[2],v[0],v[1],v[2]),f.dims.to_vec()))});
object_binding!(classical_particle_kinetic_energy,"num_native_classical_particle_kinetic_energy","crate::sim::classical::Particle3D::kinetic_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particle:Handle<c::Particle3D>)->Dynamic,|ctx:&CallCtx<'_>,particle:Handle<c::Particle3D>|{particle.with(ctx,|p|Ok(Dynamic(si().output(p.kinetic_energy(),ENERGY,false)?)))});
object_binding!(classical_particle_momentum,"num_native_classical_particle_momentum","crate::sim::classical::Particle3D::momentum","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particle:Handle<c::Particle3D>)->Dynamic,|ctx:&CallCtx<'_>,particle:Handle<c::Particle3D>|{particle.with(ctx,|p|Ok(Dynamic(si().output(p.momentum(),MOMENTUM,false)?)))});
object_binding!(classical_particle_data,"num_native_classical_particle_data","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particle:Handle<c::Particle3D>)->Dynamic,|ctx:&CallCtx<'_>,particle:Handle<c::Particle3D>|{particle.with(ctx,|p|Ok(Dynamic(particle_data(p)?)))});
object_binding!(classical_simulate_particle_motion,"num_native_classical_simulate_particle_motion","crate::sim::classical::simulate_particle_motion","Actual native RK4 particle motion. Force receives seconds plus a record with separately measured position and velocity; force results are Newtons. Positive mass/ordered finite times and step/callback budgets are checked, and first callback failure aborts immediately.",false,(force:ForceFn,mass:Input,position:Input,velocity:Input,time_range:Input,steps:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,force:ForceFn,mass:Input,position:Input,velocity:Input,time_range:Input,steps:usize,callback_budget:usize|{let f=si();let m=f.scalar(&mass.0,MASS,false)?;let p=engine::array::<3>(&f,&position.0,LENGTH)?;let v=engine::array::<3>(&f,&velocity.0,VELOCITY)?;let t=engine::array::<2>(&f,&time_range.0,TIME)?;if m<=0.0||steps==0||steps>100000||t[0]>=t[1]{return Err(api::error("particle mass/time/steps invalid"));}let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||c::simulate_particle_motion(|time,state|gate.call(||{let args=[f.output(time,TIME,false)?,api::record(vec![("position",f.array(state[..3].to_vec(),LENGTH)?),("velocity",f.array(state[3..].to_vec(),VELOCITY)?)])];let value=ctx.invoke_callback(force.0,&args)?;engine::array::<3>(&f,&value,FORCE)}),m,(p[0],p[1],p[2]),(v[0],v[1],v[2]),(t[0],t[1]),steps))?;history(&(out),&f)});
object_binding!(classical_simulate_ising_model,"num_native_classical_simulate_ising_model","crate::sim::classical::simulate_ising_model","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(size:usize,temperature:Input,exchange_reference:Input,steps:usize)->Dynamic,|_ctx:&CallCtx<'_>,size:usize,temperature:Input,exchange_reference:Input,steps:usize|{if !(2..=128).contains(&size)||steps>50_000_000{return Err(api::error("Ising lattice/proposal bounds exceeded"));}let f=si();let kt=f.scalar(&temperature.0,ENERGY,true)?;let(j,d)=api::reference(&exchange_reference.0)?;api::same_dimension(&d,&f.dim(ENERGY)?)?;if kt<0.0{return Err(api::error("temperature must benonnegative"));}let out=c::simulate_ising_model(size,kt/j,steps);let data=out.into_iter().flatten().map(|s|Value::int(i64::from(s))).collect();Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![size,size],data,elem:TypeExpr::Int})))});
object_binding!(classical_solve_1d_schrodinger,"num_native_classical_solve_1d_schrodinger","crate::sim::classical::solve_1d_schrodinger","Native normalized Schrodinger operator -onehalfLaplacian+V, mapped to physical SI using positive mass and coordinate reference with energy=hbar^2/(mass*length^2). Potential callbacks receive physical coordinates and return energy. Eigenvectors are dimensionless discrete coefficients; no continuum wavefunction normalization is fabricated. Dense eigensystem limited256nodes.",false,(potential:Potential<1>,ranges:Input,nodes:[usize;1],mass:Input,coordinate_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,potential:Potential<1>,ranges:Input,nodes:[usize;1],mass:Input,coordinate_reference:Input,callback_budget:usize|{let f=si();let m=f.scalar(&mass.0,MASS,false)?;let(l,ld)=api::reference(&coordinate_reference.0)?;api::same_dimension(&ld,&f.dim(LENGTH)?)?;let n=nodes.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).filter(|n|*n<=256).ok_or_else(||api::error("dense Schrodinger eigensystem budget256nodes"))?;if m<=0.0||nodes.iter().any(|n|*n<3){return Err(api::error("Schrodinger needspositive mass and>=3nodesperaxis"));}let Value::Tensor(t)=&ranges.0 else{return Err(api::error("coordinate range matrix required"));};
if t.shape!=[1,2]{return Err(api::error("Schrodinger ranges shape must be[D,2]"));}api::same_dimension(&api::dimension(&ranges.0)?,&ld)?;let range=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if range.iter().any(|r|r[0]>=r[1]){return Err(api::error("Schrodinger ranges mustincrease"));}let e=api::scale_product(HBAR,&[HBAR],&[m,l,l])?;if e<=0.0{return Err(api::error("Schrodinger physical energy scale not representable"));}let gate=Gate::new(ctx,callback_budget)?;let potential_call=|x:&[f64]|gate.call(||{let args=x.iter().map(|x|f.output(*x*l,LENGTH,false)).collect::<Result<Vec<_>,_>>()?;let out=ctx.invoke_callback(potential.0,&args)?;Ok(f.scalar(&out,ENERGY,false)?/e)});let(energies,vectors)=gate.run(||c::solve_1d_schrodinger(|x|potential_call(&[x]),(range[0][0]/l,range[0][1]/l),nodes[0]))?;let energies=f.array(energies.into_iter().map(|x|x*e).collect(),ENERGY)?;let eigenvectors=Value::Tensor(TensorValue{shape:vec![vectors.rows(),vectors.cols()],data:vectors.into_data().into_iter().map(Value::Float).collect(),elem:TypeExpr::Float});let _=n;Ok(Dynamic(api::record(vec![("energies",energies),("discrete_eigenvectors",eigenvectors),("coordinate_reference",coordinate_reference.0)])))});
object_binding!(classical_solve_2d_schrodinger,"num_native_classical_solve_2d_schrodinger","crate::sim::classical::solve_2d_schrodinger","Native normalized Schrodinger operator -onehalfLaplacian+V, mapped to physical SI using positive mass and coordinate reference with energy=hbar^2/(mass*length^2). Potential callbacks receive physical coordinates and return energy. Eigenvectors are dimensionless discrete coefficients; no continuum wavefunction normalization is fabricated. Dense eigensystem limited256nodes.",false,(potential:Potential<2>,ranges:Input,nodes:[usize;2],mass:Input,coordinate_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,potential:Potential<2>,ranges:Input,nodes:[usize;2],mass:Input,coordinate_reference:Input,callback_budget:usize|{let f=si();let m=f.scalar(&mass.0,MASS,false)?;let(l,ld)=api::reference(&coordinate_reference.0)?;api::same_dimension(&ld,&f.dim(LENGTH)?)?;let n=nodes.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).filter(|n|*n<=256).ok_or_else(||api::error("dense Schrodinger eigensystem budget256nodes"))?;if m<=0.0||nodes.iter().any(|n|*n<3){return Err(api::error("Schrodinger needspositive mass and>=3nodesperaxis"));}let Value::Tensor(t)=&ranges.0 else{return Err(api::error("coordinate range matrix required"));};
if t.shape!=[2,2]{return Err(api::error("Schrodinger ranges shape must be[D,2]"));}api::same_dimension(&api::dimension(&ranges.0)?,&ld)?;let range=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if range.iter().any(|r|r[0]>=r[1]){return Err(api::error("Schrodinger ranges mustincrease"));}let e=api::scale_product(HBAR,&[HBAR],&[m,l,l])?;if e<=0.0{return Err(api::error("Schrodinger physical energy scale not representable"));}let gate=Gate::new(ctx,callback_budget)?;let potential_call=|x:&[f64]|gate.call(||{let args=x.iter().map(|x|f.output(*x*l,LENGTH,false)).collect::<Result<Vec<_>,_>>()?;let out=ctx.invoke_callback(potential.0,&args)?;Ok(f.scalar(&out,ENERGY,false)?/e)});let(energies,vectors)=gate.run(||c::solve_2d_schrodinger(|x,y|potential_call(&[x,y]),(range[0][0]/l,range[0][1]/l,range[1][0]/l,range[1][1]/l),(nodes[0],nodes[1])))?;let energies=f.array(energies.into_iter().map(|x|x*e).collect(),ENERGY)?;let eigenvectors=Value::Tensor(TensorValue{shape:vec![vectors.rows(),vectors.cols()],data:vectors.into_data().into_iter().map(Value::Float).collect(),elem:TypeExpr::Float});let _=n;Ok(Dynamic(api::record(vec![("energies",energies),("discrete_eigenvectors",eigenvectors),("coordinate_reference",coordinate_reference.0)])))});
object_binding!(classical_solve_3d_schrodinger,"num_native_classical_solve_3d_schrodinger","crate::sim::classical::solve_3d_schrodinger","Native normalized Schrodinger operator -onehalfLaplacian+V, mapped to physical SI using positive mass and coordinate reference with energy=hbar^2/(mass*length^2). Potential callbacks receive physical coordinates and return energy. Eigenvectors are dimensionless discrete coefficients; no continuum wavefunction normalization is fabricated. Dense eigensystem limited256nodes.",false,(potential:Potential<3>,ranges:Input,nodes:[usize;3],mass:Input,coordinate_reference:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,potential:Potential<3>,ranges:Input,nodes:[usize;3],mass:Input,coordinate_reference:Input,callback_budget:usize|{let f=si();let m=f.scalar(&mass.0,MASS,false)?;let(l,ld)=api::reference(&coordinate_reference.0)?;api::same_dimension(&ld,&f.dim(LENGTH)?)?;let n=nodes.iter().try_fold(1usize,|a,n|a.checked_mul(*n)).filter(|n|*n<=256).ok_or_else(||api::error("dense Schrodinger eigensystem budget256nodes"))?;if m<=0.0||nodes.iter().any(|n|*n<3){return Err(api::error("Schrodinger needspositive mass and>=3nodesperaxis"));}let Value::Tensor(t)=&ranges.0 else{return Err(api::error("coordinate range matrix required"));};
if t.shape!=[3,2]{return Err(api::error("Schrodinger ranges shape must be[D,2]"));}api::same_dimension(&api::dimension(&ranges.0)?,&ld)?;let range=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if range.iter().any(|r|r[0]>=r[1]){return Err(api::error("Schrodinger ranges mustincrease"));}let e=api::scale_product(HBAR,&[HBAR],&[m,l,l])?;if e<=0.0{return Err(api::error("Schrodinger physical energy scale not representable"));}let gate=Gate::new(ctx,callback_budget)?;let potential_call=|x:&[f64]|gate.call(||{let args=x.iter().map(|x|f.output(*x*l,LENGTH,false)).collect::<Result<Vec<_>,_>>()?;let out=ctx.invoke_callback(potential.0,&args)?;Ok(f.scalar(&out,ENERGY,false)?/e)});let(energies,vectors)=gate.run(||c::solve_3d_schrodinger(|x,y,z|potential_call(&[x,y,z]),(range[0][0]/l,range[0][1]/l,range[1][0]/l,range[1][1]/l,range[2][0]/l,range[2][1]/l),(nodes[0],nodes[1],nodes[2])))?;let energies=f.array(energies.into_iter().map(|x|x*e).collect(),ENERGY)?;let eigenvectors=Value::Tensor(TensorValue{shape:vec![vectors.rows(),vectors.cols()],data:vectors.into_data().into_iter().map(Value::Float).collect(),elem:TypeExpr::Float});let _=n;Ok(Dynamic(api::record(vec![("energies",energies),("discrete_eigenvectors",eigenvectors),("coordinate_reference",coordinate_reference.0)])))});
object_binding!(classical_solve_heat_equation_1d_crank_nicolson,"num_native_classical_solve_heat_equation_1d_crank_nicolson","crate::sim::classical::solve_heat_equation_1d_crank_nicolson","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial:Initial<false>,diffusivity:Input,ranges:Input,nodes:[usize;1],time_range:Input,steps:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,initial:Initial<false>,diffusivity:Input,ranges:Input,nodes:[usize;1],time_range:Input,steps:usize,callback_budget:usize|{let f=si();let alpha=f.scalar(&diffusivity.0,[0,2,-1,0],false)?;let t=engine::array::<2>(&f,&time_range.0,TIME)?;if alpha<0.0||t[0]>=t[1]{return Err(api::error("diffusion/time range invalid"));}let _n=grid_work(&nodes,steps)?;api::same_dimension(&api::dimension(&ranges.0)?,&f.dim(LENGTH)?)?;let r=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if r.len()!=nodes.len()||r.iter().any(|r|r.len()!=2||r[0]>=r[1])||nodes.is_empty()||nodes.len()>3{return Err(api::error("heat range/node shape invalid"));}let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||c::solve_heat_equation_1d_crank_nicolson(&|x|initial_value(&gate,&initial,&[x]),alpha,(r[0][0],r[0][1]),nodes[0],(t[0],t[1]),steps))?;measured_history(out,initial.dimension)});
object_binding!(classical_solve_heat_equation_2d_crank_nicolson,"num_native_classical_solve_heat_equation_2d_crank_nicolson","crate::sim::classical::solve_heat_equation_2d_crank_nicolson","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:[usize;2],time_range:Input,steps:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:[usize;2],time_range:Input,steps:usize,callback_budget:usize|{let f=si();let alpha=f.scalar(&diffusivity.0,[0,2,-1,0],false)?;let t=engine::array::<2>(&f,&time_range.0,TIME)?;if alpha<0.0||t[0]>=t[1]{return Err(api::error("diffusion/time range invalid"));}let _n=grid_work(&nodes,steps)?;api::same_dimension(&api::dimension(&ranges.0)?,&f.dim(LENGTH)?)?;let r=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if r.len()!=nodes.len()||r.iter().any(|r|r.len()!=2||r[0]>=r[1])||nodes.is_empty()||nodes.len()>3{return Err(api::error("heat range/node shape invalid"));}let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||c::solve_heat_equation_2d_crank_nicolson(&|x,y|initial_value(&gate,&initial,&[x,y]),alpha,[(r[0][0],r[0][1]),(r[1][0],r[1][1])],nodes,(t[0],t[1]),steps))?;measured_history(out,initial.dimension)});
object_binding!(classical_solve_heat_equation_3d_crank_nicolson,"num_native_classical_solve_heat_equation_3d_crank_nicolson","crate::sim::classical::solve_heat_equation_3d_crank_nicolson","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:[usize;3],time_range:Input,steps:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:[usize;3],time_range:Input,steps:usize,callback_budget:usize|{let f=si();let alpha=f.scalar(&diffusivity.0,[0,2,-1,0],false)?;let t=engine::array::<2>(&f,&time_range.0,TIME)?;if alpha<0.0||t[0]>=t[1]{return Err(api::error("diffusion/time range invalid"));}let _n=grid_work(&nodes,steps)?;api::same_dimension(&api::dimension(&ranges.0)?,&f.dim(LENGTH)?)?;let r=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if r.len()!=nodes.len()||r.iter().any(|r|r.len()!=2||r[0]>=r[1])||nodes.is_empty()||nodes.len()>3{return Err(api::error("heat range/node shape invalid"));}let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||c::solve_heat_equation_3d_crank_nicolson(&|x,y,z|initial_value(&gate,&initial,&[x,y,z]),alpha,[(r[0][0],r[0][1]),(r[1][0],r[1][1]),(r[2][0],r[2][1])],nodes,(t[0],t[1]),steps))?;measured_history(out,initial.dimension)});
object_binding!(classical_solve_heat_equation_nd_crank_nicolson,"num_native_classical_solve_heat_equation_nd_crank_nicolson","crate::sim::classical::solve_heat_equation_nd_crank_nicolson","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:Vec<usize>,time_range:Input,steps:usize,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,initial:Initial<true>,diffusivity:Input,ranges:Input,nodes:Vec<usize>,time_range:Input,steps:usize,callback_budget:usize|{let f=si();let alpha=f.scalar(&diffusivity.0,[0,2,-1,0],false)?;let t=engine::array::<2>(&f,&time_range.0,TIME)?;if alpha<0.0||t[0]>=t[1]{return Err(api::error("diffusion/time range invalid"));}let _n=grid_work(&nodes,steps)?;api::same_dimension(&api::dimension(&ranges.0)?,&f.dim(LENGTH)?)?;let r=<Vec<Vec<f64>> as Arg>::decode(&ranges.0,0)?;if r.len()!=nodes.len()||r.iter().any(|r|r.len()!=2||r[0]>=r[1])||nodes.is_empty()||nodes.len()>3{return Err(api::error("heat range/node shape invalid"));}let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||c::solve_heat_equation_nd_crank_nicolson(&|x|initial_value(&gate,&initial,x),alpha,&r.iter().map(|r|(r[0],r[1])).collect::<Vec<_>>(),&nodes,(t[0],t[1]),steps))?;measured_history(out,initial.dimension)});
object_binding!(classical_solve_wave_equation_1d,"num_native_classical_solve_wave_equation_1d","crate::sim::classical::solve_wave_equation_1d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial_u:Input,initial_ut:Input,spacing:Input,speed:Input,dt:Input,steps:usize)->Dynamic,|_ctx:&CallCtx<'_>,initial_u:Input,initial_ut:Input,spacing:Input,speed:Input,dt:Input,steps:usize|{let f=si();let d=api::dimension(&initial_u.0)?;let rate=DimExpr::Mul(vec![d.clone(),DimExpr::Inverse(Box::new(f.dim(TIME)?))]);api::same_dimension(&rate,&api::dimension(&initial_ut.0)?)?;let u=<Vec<f64> as Arg>::decode(&initial_u.0,0)?;let ut=<Vec<f64> as Arg>::decode(&initial_ut.0,0)?;let nodes=[u.len()];let n=grid_work(&nodes,steps)?;if u.len()!=n||ut.len()!=n{return Err(api::error("wave initialvectors shape mismatch"));}let c=f.scalar(&speed.0,VELOCITY,false)?;let dt=f.scalar(&dt.0,TIME,false)?;if c<0.0||dt<=0.0{return Err(api::error("wave speed/time step invalid"));}let h=f.scalar(&spacing.0,LENGTH,false)?;if h<=0.0{return Err(api::error("wave spacing must bepositive"));}let out=c::solve_wave_equation_1d(&u,&ut,c,h,dt,steps).map_err(api::error)?;measured_history(out,d)});
object_binding!(classical_solve_wave_equation_2d,"num_native_classical_solve_wave_equation_2d","crate::sim::classical::solve_wave_equation_2d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial_u:Input,initial_ut:Input,nodes:[usize;2],spacing:Input,speed:Input,dt:Input,steps:usize)->Dynamic,|_ctx:&CallCtx<'_>,initial_u:Input,initial_ut:Input,nodes:[usize;2],spacing:Input,speed:Input,dt:Input,steps:usize|{let f=si();let d=api::dimension(&initial_u.0)?;let rate=DimExpr::Mul(vec![d.clone(),DimExpr::Inverse(Box::new(f.dim(TIME)?))]);api::same_dimension(&rate,&api::dimension(&initial_ut.0)?)?;let u=<Vec<f64> as Arg>::decode(&initial_u.0,0)?;let ut=<Vec<f64> as Arg>::decode(&initial_ut.0,0)?;let n=grid_work(&nodes,steps)?;if u.len()!=n||ut.len()!=n{return Err(api::error("wave initialvectors shape mismatch"));}let c=f.scalar(&speed.0,VELOCITY,false)?;let dt=f.scalar(&dt.0,TIME,false)?;if c<0.0||dt<=0.0{return Err(api::error("wave speed/time step invalid"));}let h=engine::array::<2>(&f,&spacing.0,LENGTH)?;if h.iter().any(|h|*h<=0.0){return Err(api::error("wave spacings mustbepositive"));}let out=c::solve_wave_equation_2d(&u,&ut,nodes,h,c,dt,steps).map_err(api::error)?;measured_history(out,d)});
object_binding!(classical_solve_wave_equation_3d,"num_native_classical_solve_wave_equation_3d","crate::sim::classical::solve_wave_equation_3d","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(initial_u:Input,initial_ut:Input,nodes:[usize;3],spacing:Input,speed:Input,dt:Input,steps:usize)->Dynamic,|_ctx:&CallCtx<'_>,initial_u:Input,initial_ut:Input,nodes:[usize;3],spacing:Input,speed:Input,dt:Input,steps:usize|{let f=si();let d=api::dimension(&initial_u.0)?;let rate=DimExpr::Mul(vec![d.clone(),DimExpr::Inverse(Box::new(f.dim(TIME)?))]);api::same_dimension(&rate,&api::dimension(&initial_ut.0)?)?;let u=<Vec<f64> as Arg>::decode(&initial_u.0,0)?;let ut=<Vec<f64> as Arg>::decode(&initial_ut.0,0)?;let n=grid_work(&nodes,steps)?;if u.len()!=n||ut.len()!=n{return Err(api::error("wave initialvectors shape mismatch"));}let c=f.scalar(&speed.0,VELOCITY,false)?;let dt=f.scalar(&dt.0,TIME,false)?;if c<0.0||dt<=0.0{return Err(api::error("wave speed/time step invalid"));}let h=engine::array::<3>(&f,&spacing.0,LENGTH)?;if h.iter().any(|h|*h<=0.0){return Err(api::error("wave spacings mustbepositive"));}let out=c::solve_wave_equation_3d(&u,&ut,nodes,h,c,dt,steps).map_err(api::error)?;measured_history(out,d)});
object_binding!(classical_projectile_motion_with_drag,"num_native_classical_projectile_motion_with_drag","crate::sim::classical::projectile_motion_with_drag","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(parameters:Input)->Dynamic,|_ctx:&CallCtx<'_>,parameters:Input|{let f=si();let read=|key,p|f.scalar(api::field(&parameters.0,key)?,p,false);let v0=read("speed",VELOCITY)?;let angle=crate::rules::numerical_abi::numeric(api::field(&parameters.0,"angle")?,0)?;api::normalized(api::field(&parameters.0,"angle")?)?;let mass=read("mass",MASS)?;let drag_coeff=crate::rules::numerical_abi::numeric(api::field(&parameters.0,"drag_coefficient")?,0)?;api::normalized(api::field(&parameters.0,"drag_coefficient")?)?;let area=read("area",[0,2,0,0])?;let air_density=read("air_density",[1,-3,0,0])?;let dt=read("dt",TIME)?;let max_time=read("max_time",TIME)?;if mass<=0.0||v0<0.0||drag_coeff<0.0||area<0.0||air_density<0.0||dt<=0.0||max_time<0.0||max_time/dt>100000.0{return Err(api::error("projectile physical controls/step budget invalid"));}let out=c::projectile_motion_with_drag(c::ProjectileParams{v0,angle,mass,drag_coeff,area,air_density,dt,max_time});let mut times=vec![];let mut pos=vec![];let mut vel=vec![];for(t,x,y,vx,vy)in out{times.push(t);pos.push([x,y]);vel.push([vx,vy]);}Ok(Dynamic(api::record(vec![("time",f.array(times,TIME)?),("positions",engine::packed(&f,pos,LENGTH)?),("velocities",engine::packed(&f,vel,VELOCITY)?)])))});
object_binding!(classical_simulate_n_body,"num_native_classical_simulate_n_body","crate::sim::classical::simulate_n_body","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particles:Particles,dt:Input,steps:usize,gravity_constant:Input)->Dynamic,|_ctx:&CallCtx<'_>,particles:Particles,dt:Input,steps:usize,gravity_constant:Input|{let f=si();let dt=f.scalar(&dt.0,TIME,false)?;let g=f.scalar(&gravity_constant.0,[-1,3,-2,0],false)?;if dt<=0.0||steps>100000||particles.0.len().checked_mul(particles.0.len()).and_then(|n|n.checked_mul(steps)).is_none_or(|n|n>50_000_000)||particles.0.len().checked_mul(steps+1).is_none_or(|n|n>131072){return Err(api::error("classical Nbody time/work/history bounds invalid"));}let history=c::simulate_n_body(particles.0,dt,steps,g);let mut data=vec![];for row in history{let values=row.iter().map(particle_data).collect::<Result<Vec<_>,_>>()?;data.push(Value::Tensor(TensorValue{shape:vec![values.len()],data:values,elem:api::dynamic_type()}));}Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![data.len()],data,elem:api::dynamic_type()})))});
object_binding!(classical_gravitational_potential_energy,"num_native_classical_gravitational_potential_energy","crate::sim::classical::gravitational_potential_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particles:Particles,gravity_constant:Input)->Dynamic,|_ctx:&CallCtx<'_>,particles:Particles,gravity_constant:Input|{let f=si();Ok(Dynamic(f.output(c::gravitational_potential_energy(&particles.0,f.scalar(&gravity_constant.0,[-1,3,-2,0],false)?),ENERGY,false)?))});
object_binding!(classical_total_kinetic_energy,"num_native_classical_total_kinetic_energy","crate::sim::classical::total_kinetic_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(particles:Particles)->Dynamic,|_ctx:&CallCtx<'_>,particles:Particles|{let f=si();Ok(Dynamic(f.output(c::total_kinetic_energy(&particles.0),ENERGY,false)?))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &classical_particle_new::BINDING,
            &classical_particle_kinetic_energy::BINDING,
            &classical_particle_momentum::BINDING,
            &classical_particle_data::BINDING,
            &classical_simulate_particle_motion::BINDING,
            &classical_simulate_ising_model::BINDING,
            &classical_solve_1d_schrodinger::BINDING,
            &classical_solve_2d_schrodinger::BINDING,
            &classical_solve_3d_schrodinger::BINDING,
            &classical_solve_heat_equation_1d_crank_nicolson::BINDING,
            &classical_solve_heat_equation_2d_crank_nicolson::BINDING,
            &classical_solve_heat_equation_3d_crank_nicolson::BINDING,
            &classical_solve_heat_equation_nd_crank_nicolson::BINDING,
            &classical_solve_wave_equation_1d::BINDING,
            &classical_solve_wave_equation_2d::BINDING,
            &classical_solve_wave_equation_3d::BINDING,
            &classical_projectile_motion_with_drag::BINDING,
            &classical_simulate_n_body::BINDING,
            &classical_gravitational_potential_energy::BINDING,
            &classical_total_kinetic_energy::BINDING,
        ],
        "numerical_vector",
    )
}
