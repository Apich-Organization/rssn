#![allow(clippy::wildcard_imports)] // Native split-step wave solvers with hbar and independent physical coordinate/time references.
mod linear;
mod models;
use super::*;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::sim::split_step::SplitOrder;
use crate::sim::split_step::SplitStep;
use num_complex::Complex64;
use std::cell::Cell;
use std::cell::RefCell;
const HBAR: f64 = 1.054_571_817e-34;
struct Wave {
    native: SplitStep,
    length: f64,
    time: f64,
    energy: f64,
    physical: bool,
}
impl NativeType for Wave {
    const KIND: &'static str = "rssn/native/sim/split-step/v1";
    const MAX_BYTES: usize = 8192;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.grid.size() > 65536
            || self.native.mass <= 0.0
            || [self.length, self.time, self.energy, self.native.mass]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(api::error("split-step grid/reference/mass invalid"));
        }
        Ok(())
    }
}
impl Wave {
    fn dim(
        &self,
        kind: &str,
    ) -> DimExpr {
        if !self.physical {
            return DimExpr::Dimless;
        }
        match kind {
            | "length" => DimExpr::base("Length"),
            | "time" => DimExpr::base("Time"),
            | "energy" => energy_dimension(),
            | "psi" => self.samples().dimension(),
            | _ => DimExpr::Dimless,
        }
    }

    fn read(
        &self,
        v: &Value,
        kind: &str,
        scale: f64,
    ) -> Result<f64, OpError> {
        api::same_dimension(&api::dimension(v)?, &self.dim(kind))?;
        api::scale_product(crate::rules::numerical_abi::numeric(v, 0)?, &[], &[scale])
    }

    fn out(
        &self,
        v: f64,
        kind: &str,
        scale: f64,
    ) -> Result<Value, OpError> {
        api::quantity(api::scale_product(v, &[scale], &[])?, self.dim(kind))
    }

    fn samples(&self) -> WaveSamples {
        WaveSamples {
            dimensions: self.native.grid.dim,
            size: self.native.grid.size(),
            reference: if self.physical {
                Some(self.length)
            } else {
                None
            },
        }
    }

    fn psi_factor(&self) -> Result<f64, OpError> {
        self.samples().factor()
    }

    fn psi(
        &self,
        value: &Input,
    ) -> Result<Vec<Complex64>, OpError> {
        self.samples().decode(value)
    }

    fn psi_output(
        &self,
        values: Vec<Complex64>,
    ) -> Result<Dynamic, OpError> {
        self.samples().encode(values)
    }

    fn coupling(
        &self,
        v: &Value,
    ) -> Result<f64, OpError> {
        let dim = DimExpr::Mul(vec![
            self.dim("energy"),
            self.dim("length")
                .pow(Exp::Rat(Rational::int(self.native.grid.dim as i64))),
        ]);
        api::same_dimension(&api::dimension(v)?, &dim)?;
        let l = vec![self.length; self.native.grid.dim];
        let mut den = vec![self.energy];
        den.extend(l);
        api::scale_product(crate::rules::numerical_abi::numeric(v, 0)?, &[], &den)
    }

    fn point(
        &self,
        x: [f64; 3],
    ) -> Result<Value, OpError> {
        let d = self.dim("length");
        Ok(Value::Tensor(TensorValue {
            shape: vec![3],
            data: x
                .into_iter()
                .map(|v| self.out(v, "length", self.length))
                .collect::<Result<_, _>>()?,
            elem: TypeExpr::Float.with_dim(d),
        }))
    }
}
fn energy_dimension() -> DimExpr {
    DimExpr::Mul(vec![
        DimExpr::base("Mass"),
        DimExpr::base("Length").pow(2),
        DimExpr::base("Time").pow(-2),
    ])
}

struct WaveSamples {
    dimensions: usize,
    size: usize,
    reference: Option<f64>,
}

impl WaveSamples {
    fn dimension(&self) -> DimExpr {
        if self.reference.is_some() {
            DimExpr::base("Length").pow(Exp::Rat(Rational::new(-(self.dimensions as i64), 2)))
        } else {
            DimExpr::Dimless
        }
    }

    fn factor(&self) -> Result<f64, OpError> {
        let f = self
            .reference
            .map_or(1.0, |length| length.powf(self.dimensions as f64 / 2.0));
        if f <= 0.0 || !f.is_finite() {
            return Err(api::error(
                "wavefunction reference amplitude not representable",
            ));
        }
        Ok(f)
    }

    fn decode(
        &self,
        v: &Input,
    ) -> Result<Vec<Complex64>, OpError> {
        self.decode_value(&v.0)
    }

    fn decode_value(
        &self,
        value: &Value,
    ) -> Result<Vec<Complex64>, OpError> {
        api::same_dimension(&api::dimension(value)?, &self.dimension())?;
        let Value::Tensor(t) = value else {
            return Err(api::error("wavefunction tensor expected"));
        };
        if t.shape != [self.size] || t.data.len() != self.size {
            return Err(api::error(
                "wavefunction must be a flat tensor matching the grid",
            ));
        }
        let factor = self.factor()?;
        t.data
            .iter()
            .map(|v| {
                let c = complex(v)? * factor;
                if !c.re.is_finite() || !c.im.is_finite() {
                    return Err(api::error("wavefunction data must be finite"));
                }
                Ok(c)
            })
            .collect()
    }

    fn encode(
        &self,
        v: Vec<Complex64>,
    ) -> Result<Dynamic, OpError> {
        let factor = self.factor()?;
        let data = v
            .into_iter()
            .map(|z| {
                let z = z / factor;
                if !z.re.is_finite() || !z.im.is_finite() {
                    return Err(api::error("native wavefunction result is nonfinite"));
                }
                Output::encode(z)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Dynamic(Value::Tensor(TensorValue {
            shape: vec![data.len()],
            data,
            elem: TypeExpr::Complex.with_dim(self.dimension()),
        })))
    }
}

fn complex(v: &Value) -> Result<Complex64, OpError> {
    match v {
        | Value::Dyn(d) => complex(&d.value),
        | Value::Quantity(q) => {
            let unit = Value::Quantity(crate::kernel::value::QuantityValue {
                value: Box::new(Value::Float(1.0)),
                unit: q.unit.clone(),
            });
            Ok(complex(&q.value)? * crate::rules::numerical_abi::numeric(&unit, 0)?)
        },
        | _ => <Complex64 as Arg>::decode(v, 0),
    }
}
fn sizes(
    n: &[usize],
    lengths: &[f64],
) -> Result<(), OpError> {
    if n.is_empty()
        || n.len() > 3
        || n.len() != lengths.len()
        || n.iter().any(|n| *n < 4 || !n.is_power_of_two())
        || n.iter()
            .try_fold(1usize, |a, n| a.checked_mul(*n))
            .is_none_or(|n| n > 65536)
        || lengths.iter().any(|v| *v <= 0.0 || !v.is_finite())
    {
        return Err(api::error(
            "spectral grid needs1..3poweroftwoaxes>=4,<=65536points,positivefinite lengths",
        ));
    }
    Ok(())
}
fn order(v: &str) -> Result<SplitOrder, OpError> {
    match v {
        | "Strang" => Ok(SplitOrder::Strang),
        | "Yoshida4" => Ok(SplitOrder::Yoshida4),
        | _ => Err(api::error("split order mustbeStrang/Yoshida4")),
    }
}
struct SampleFn(CallbackRef);
struct PotentialFn(CallbackRef);
struct AbsorberFn(CallbackRef);
fn variable(
    n: &str,
    t: TypeExpr,
) -> TypeExpr {
    t.with_dim(DimExpr::Var(n.into()))
}
impl NativeArg for SampleFn {
    fn ty() -> TypeExpr {
        TypeExpr::Function(
            vec![api::tensor_type(1, variable("X", TypeExpr::Float))],
            Box::new(variable("P", TypeExpr::Complex)),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, 1)?))
    }
}
impl NativeArg for PotentialFn {
    fn ty() -> TypeExpr {
        TypeExpr::Function(
            vec![
                variable("T", TypeExpr::Float),
                api::tensor_type(1, variable("X", TypeExpr::Float)),
            ],
            Box::new(variable("E", TypeExpr::Float)),
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
impl NativeArg for AbsorberFn {
    fn ty() -> TypeExpr {
        TypeExpr::Function(
            vec![api::tensor_type(1, variable("X", TypeExpr::Float))],
            Box::new(variable("E", TypeExpr::Float)),
        )
    }

    fn decode(
        v: &Value,
        _: usize,
        c: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        Ok(Self(callback(v, c, 1)?))
    }
}
fn callback(
    v: &Value,
    c: &CallCtx<'_>,
    arity: usize,
) -> Result<CallbackRef, OpError> {
    let Value::Function(h) = v else {
        return Err(api::error("wave callback expected"));
    };
    if c.callback_signature(*h)?.params.len() != arity {
        return Err(api::error("wave callback arity mismatch"));
    }
    Ok(*h)
}
struct Stop;
struct Gate<'a, 'b> {
    ctx: &'a CallCtx<'b>,
    remaining: Cell<usize>,
    failure: RefCell<Option<OpError>>,
}
impl<'a, 'b> Gate<'a, 'b> {
    fn new(
        ctx: &'a CallCtx<'b>,
        budget: usize,
    ) -> Result<Self, OpError> {
        if budget == 0 || budget > 10_000_000 {
            return Err(api::error("wave callback budget mustbe1..10million"));
        }
        Ok(Self {
            ctx,
            remaining: Cell::new(budget),
            failure: RefCell::new(None),
        })
    }

    fn invoke<T>(
        &self,
        work: impl FnOnce() -> Result<T, OpError>,
    ) -> T {
        let result = self.remaining.get().checked_sub(1).map_or_else(
            || Err(api::error("wave callback budget exhausted")),
            |left| {
                self.remaining.set(left);
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
        work: impl FnOnce() -> T,
    ) -> Result<T, OpError> {
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(work));
        if let Some(e) = self.failure.borrow_mut().take() {
            return Err(e);
        }
        match result {
            | Ok(v) => Ok(v),
            | Err(e) => std::panic::resume_unwind(e),
        }
    }

    fn potential(
        &self,
        w: &Wave,
        h: CallbackRef,
        t: f64,
        x: [f64; 3],
    ) -> f64 {
        self.invoke(|| {
            let value = self
                .ctx
                .invoke_callback(h, &[w.out(t, "time", w.time)?, w.point(x)?])?;
            w.read(&value, "energy", w.energy)
        })
    }
}
fn evolution(
    size: usize,
    steps: usize,
) -> Result<(), OpError> {
    if steps > 100000 || steps.checked_mul(size).is_none_or(|n| n > 50_000_000) {
        return Err(api::error("splitstep grid-step work limit exceeded"));
    }
    Ok(())
}
struct PhysicalWaveReferences {
    length: f64,
    time: f64,
    mass: f64,
    energy: f64,
}

impl PhysicalWaveReferences {
    fn new(
        length: &Value,
        time: &Value,
        mass: &Value,
    ) -> Result<Self, OpError> {
        let (length, length_dimension) = api::reference(length)?;
        let (time, time_dimension) = api::reference(time)?;
        api::same_dimension(&length_dimension, &DimExpr::base("Length"))?;
        api::same_dimension(&time_dimension, &DimExpr::base("Time"))?;
        api::same_dimension(&api::dimension(mass)?, &DimExpr::base("Mass"))?;
        let mass = crate::rules::numerical_abi::numeric(mass, 0)?;
        let mass = api::scale_product(mass, &[length, length], &[HBAR, time])?;
        let energy = api::scale_product(HBAR, &[], &[time])?;
        if mass <= 0.0 || energy <= 0.0 {
            return Err(api::error(
                "quantum mass and reference scales must be positive and representable",
            ));
        }
        Ok(Self {
            length,
            time,
            mass,
            energy,
        })
    }

    fn lengths(
        &self,
        lengths: &Value,
    ) -> Result<Vec<f64>, OpError> {
        api::same_dimension(&api::dimension(lengths)?, &DimExpr::base("Length"))?;
        let mut lengths = <Vec<f64> as Arg>::decode(lengths, 0)?;
        for length in &mut lengths {
            *length /= self.length;
        }
        Ok(lengths)
    }
}

object_binding!(split_centered,"num_native_split_centered","crate::sim::split_step::SplitStep::centered","Actual native split-step solver in1D/2D/3D with independent physical coordinate/time references, SI mass and hbar. Potential/energy scale=hbar/time_reference; normalizedmass=mL²/(hbarT). Wavefunctions carry Length^(-dimension/2), preserving dimensionless norm; contact coupling has Energy*Length^dimension units.",false,(nodes:Vec<usize>,lengths:Input,mass:Input,coordinate_reference:Input,time_reference:Input)->NewObject<Wave>,|_ctx:&CallCtx<'_>,nodes:Vec<usize>,lengths:Input,mass:Input,coordinate_reference:Input,time_reference:Input|{
    let refs = PhysicalWaveReferences::new(&coordinate_reference.0, &time_reference.0, &mass.0)?;
    let lengths = refs.lengths(&lengths.0)?;
    sizes(&nodes, &lengths)?;
    let native = SplitStep::centered(&nodes, &lengths, refs.mass).map_err(api::error)?;
    Ok(NewObject::new(Wave {native, length: refs.length, time: refs.time, energy: refs.energy, physical: true}, vec![DimExpr::base("Length"), DimExpr::base("Time")]))
});
object_binding!(split_centered_normalized,"num_native_split_centered_normalized","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(nodes:Vec<usize>,lengths:Vec<f64>,mass:f64)->NewObject<Wave>,|_ctx:&CallCtx<'_>,nodes:Vec<usize>,lengths:Vec<f64>,mass:f64|{sizes(&nodes,&lengths)?;if mass<=0.0{return Err(api::error("normalized mass mustbepositive"));}Ok(NewObject::normalized(Wave{native:SplitStep::centered(&nodes,&lengths,mass).map_err(api::error)?,length:1.0,time:1.0,energy:1.0,physical:false}))});
object_binding!(split_x,"num_native_split_x","crate::sim::split_step::SplitStep::x","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,index:usize|{solver.with(ctx,|w|{if index>=w.native.grid.size(){return Err(api::error("spectral point index outsidegrid"));}Ok(Dynamic(w.point(w.native.x(index))?))})});
object_binding!(split_norm2,"num_native_split_norm2","crate::sim::split_step::SplitStep::norm2","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input)->f64,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input|{solver.with(ctx,|w|Ok(w.native.norm2(&w.psi(&psi)?)))});
object_binding!(split_kinetic_energy,"num_native_split_kinetic_energy","crate::sim::split_step::SplitStep::kinetic_energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input|{solver.with(ctx,|w|Ok(Dynamic(w.out(w.native.kinetic_energy(&w.psi(&psi)?),"energy",w.energy)?)))});
object_binding!(split_sample,"num_native_split_sample","crate::sim::split_step::SplitStep::sample","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,function:SampleFn,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,function:SampleFn,callback_budget:usize|{solver.with(ctx,|w|{let gate=Gate::new(ctx,callback_budget)?;let factor=w.psi_factor()?;let values=gate.run(||w.native.sample(|x|gate.invoke(||{let out=ctx.invoke_callback(function.0,&[w.point(x)?])?;api::same_dimension(&api::dimension(&out)?,&w.dim("psi"))?;Ok(complex(&out)?*factor)})))?;w.psi_output(values)})});
object_binding!(split_energy,"num_native_split_energy","crate::sim::split_step::SplitStep::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,time:Input,coupling:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,time:Input,coupling:Input,callback_budget:usize|{solver.with(ctx,|w|{let psi=w.psi(&psi)?;let t=w.read(&time.0,"time",w.time)?;let g=w.coupling(&coupling.0)?;let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||w.native.energy(&psi,&|t,x|gate.potential(w,potential.0,t,x),t,g))?;Ok(Dynamic(w.out(out,"energy",w.energy)?))})});
object_binding!(split_chemical_potential,"num_native_split_chemical_potential","crate::sim::split_step::SplitStep::chemical_potential","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,time:Input,coupling:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,time:Input,coupling:Input,callback_budget:usize|{solver.with(ctx,|w|{let psi=w.psi(&psi)?;let t=w.read(&time.0,"time",w.time)?;let g=w.coupling(&coupling.0)?;let gate=Gate::new(ctx,callback_budget)?;let out=gate.run(||w.native.chemical_potential(&psi,&|t,x|gate.potential(w,potential.0,t,x),t,g))?;Ok(Dynamic(w.out(out,"energy",w.energy)?))})});
object_binding!(split_rotate,"num_native_split_rotate","crate::sim::split_step::SplitStep::rotate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,angle:f64)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,angle:f64|{solver.with(ctx,|w|{if w.native.grid.dim<2||angle.abs()>=std::f64::consts::FRAC_PI_2{return Err(api::error("native shear rotation needs>=2D and |angle|<pi/2"));}let mut psi=w.psi(&psi)?;w.native.rotate(&mut psi,angle);w.psi_output(psi)})});
object_binding!(split_evolve,"num_native_split_evolve","crate::sim::split_step::SplitStep::evolve","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize|{solver.with(ctx,|w|{evolution(w.native.grid.size(),steps)?;let psi=w.psi(&psi)?;let g=w.coupling(&coupling.0)?;let dt=w.read(&dt.0,"time",w.time)?;if dt<=0.0{return Err(api::error("wave timestep must bepositive"));}let gate=Gate::new(ctx,callback_budget)?;let t=w.read(&initial_time.0,"time",w.time)?;let order=order(&scheme)?;let out=gate.run(||w.native.evolve(&psi,&|t,x|gate.potential(w,potential.0,t,x),g,t,dt,steps,order))?;w.psi_output(out)})});
object_binding!(split_evolve_rotating,"num_native_split_evolve_rotating","crate::sim::split_step::SplitStep::evolve_rotating","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,angular_velocity:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,angular_velocity:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize|{solver.with(ctx,|w|{evolution(w.native.grid.size(),steps)?;let psi=w.psi(&psi)?;let g=w.coupling(&coupling.0)?;let dt=w.read(&dt.0,"time",w.time)?;if dt<=0.0{return Err(api::error("wave timestep must bepositive"));}let gate=Gate::new(ctx,callback_budget)?;let t=w.read(&initial_time.0,"time",w.time)?;let order=order(&scheme)?;api::same_dimension(&api::dimension(&angular_velocity.0)?,&w.dim("time").pow(-1))?;let omega=crate::rules::numerical_abi::numeric(&angular_velocity.0,0)?*w.time;if w.native.grid.dim<2||omega*dt!=0.0&&omega.abs()*dt>=1.0{return Err(api::error("rotating splitstep needs>=2D and bounded rotationpersubstep"));}let out=gate.run(||w.native.evolve_rotating(&psi,&|t,x|gate.potential(w,potential.0,t,x),g,omega,t,dt,steps,order))?;w.psi_output(out)})});
object_binding!(split_evolve_absorbing,"num_native_split_evolve_absorbing","crate::sim::split_step::SplitStep::evolve_absorbing","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,absorber:AbsorberFn,coupling:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,absorber:AbsorberFn,coupling:Input,initial_time:Input,dt:Input,steps:usize,scheme:String,callback_budget:usize|{solver.with(ctx,|w|{evolution(w.native.grid.size(),steps)?;let psi=w.psi(&psi)?;let g=w.coupling(&coupling.0)?;let dt=w.read(&dt.0,"time",w.time)?;if dt<=0.0{return Err(api::error("wave timestep must bepositive"));}let gate=Gate::new(ctx,callback_budget)?;let t=w.read(&initial_time.0,"time",w.time)?;let order=order(&scheme)?;let absorption=|x|gate.invoke(||{let out=ctx.invoke_callback(absorber.0,&[w.point(x)?])?;let value=w.read(&out,"energy",w.energy)?;if value<0.0{return Err(api::error("absorbing potential must benonnegative"));}Ok(value)});let out=gate.run(||w.native.evolve_absorbing(&psi,&|t,x|gate.potential(w,potential.0,t,x),&absorption,g,t,dt,steps,order))?;w.psi_output(out)})});
object_binding!(split_ground_state,"num_native_split_ground_state","crate::sim::split_step::SplitStep::ground_state","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,dt:Input,steps:usize,tolerance:Input,callback_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,solver:Handle<Wave>,psi:Input,potential:PotentialFn,coupling:Input,dt:Input,steps:usize,tolerance:Input,callback_budget:usize|{solver.with(ctx,|w|{evolution(w.native.grid.size(),steps)?;let psi=w.psi(&psi)?;let g=w.coupling(&coupling.0)?;let dt=w.read(&dt.0,"time",w.time)?;if dt<=0.0{return Err(api::error("wave timestep must bepositive"));}let gate=Gate::new(ctx,callback_budget)?;let tolerance=w.read(&tolerance.0,"energy",w.energy)?;if tolerance<=0.0||steps==0||w.native.norm2(&psi)<=0.0{return Err(api::error("groundstate needspositive tolerance,stepbudget,andnonzero guess"));}let out=gate.run(||w.native.ground_state(&psi,&|t,x|gate.potential(w,potential.0,t,x),g,dt,steps,tolerance))?;let psi=w.psi_output(out.psi)?.0;Ok(Dynamic(api::record(vec![("psi",psi),("energy",w.out(out.energy,"energy",w.energy)?),("chemical_potential",w.out(out.mu,"energy",w.energy)?),("steps",Value::int(out.steps)),("converged",Value::Bool(out.converged))])))})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    linear::install(i)?;
    models::install(i)?;
    api::install(
        i,
        &[
            &split_centered::BINDING,
            &split_centered_normalized::BINDING,
            &split_x::BINDING,
            &split_norm2::BINDING,
            &split_kinetic_energy::BINDING,
            &split_sample::BINDING,
            &split_energy::BINDING,
            &split_chemical_potential::BINDING,
            &split_rotate::BINDING,
            &split_evolve::BINDING,
            &split_evolve_rotating::BINDING,
            &split_evolve_absorbing::BINDING,
            &split_ground_state::BINDING,
        ],
        "numerical_vector",
    )
}
