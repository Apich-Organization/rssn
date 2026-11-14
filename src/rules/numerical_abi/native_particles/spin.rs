#![allow(clippy::wildcard_imports)] // Native spin and quantum Monte Carlo use explicit exchange-energy reference and Kelvin/kBT conversion.
use super::*;
use crate::kernels::random::Rng;
use crate::rules::numerical_abi::codecs::Arg;
use crate::rules::numerical_abi::codecs::Output;
use crate::sim::models_ising_statistical as ising;
use crate::sim::spin_models as spin;
use crate::sim::sse;
struct Lattice<const D: usize>(spin::Lattice<D>);
struct On<const D: usize, const N: usize> {
    native: spin::OnModel<D, N>,
    frame: Frame,
}
struct Potts<const D: usize> {
    native: spin::Potts<D>,
    frame: Frame,
}
struct Temper<const D: usize, const N: usize> {
    native: spin::ParallelTempering<D, N>,
    frame: Frame,
}
struct PottsTemper<const D: usize> {
    native: spin::PottsTempering<D>,
    frame: Frame,
}
fn sites<const D: usize>(
    dims: [usize; D],
    periodic: bool,
) -> Result<usize, OpError> {
    if dims.iter().any(|n| *n < if periodic { 3 } else { 1 }) {
        return Err(api::error("spin lattice side below boundary requirement"));
    }
    dims.iter()
        .try_fold(1usize, |a, n| a.checked_mul(*n))
        .filter(|n| *n <= 4096)
        .ok_or_else(|| api::error("spin lattice site allocation bound4096exceeded"))
}
fn finite(v: impl IntoIterator<Item = f64>) -> Result<(), OpError> {
    if v.into_iter().any(|x| !x.is_finite()) {
        return Err(api::error("spin state/result is nonfinite"));
    }
    Ok(())
}
fn beta(
    f: &Frame,
    v: &Value,
) -> Result<f64, OpError> {
    let b = f.scalar(v, [-1, -2, 2, 0], false)?;
    if b < 0.0 || b > 100.0 {
        return Err(api::error("normalized inverse-temperature must be0..100"));
    }
    Ok(b)
}
fn betas(
    f: &Frame,
    v: &Value,
) -> Result<Vec<f64>, OpError> {
    let b = f.vector(v, [-1, -2, 2, 0])?;
    if b.len() < 2 || b.len() > 32 || b.iter().any(|b| *b < 0.0 || *b > 100.0) {
        return Err(api::error(
            "tempering ladder needs2..32inverse temperatures in0..100normalized range",
        ));
    }
    Ok(b)
}
fn update(
    kind: &str,
    delta: f64,
) -> Result<spin::Update, OpError> {
    match kind {
        | "Metropolis" if delta > 0.0 && delta <= 4.0 => Ok(spin::Update::Metropolis(delta)),
        | "HeatBath" => Ok(spin::Update::HeatBath),
        | "Wolff" => Ok(spin::Update::Wolff),
        | "SwendsenWang" => Ok(spin::Update::SwendsenWang),
        | _ => {
            Err(api::error(
                "spin update must beMetropolis(delta in(0,4])/HeatBath/Wolff/SwendsenWang",
            ))
        },
    }
}
fn work(
    sites: usize,
    replicas: usize,
    rounds: usize,
) -> Result<(), OpError> {
    if rounds > 1_000_000
        || sites
            .checked_mul(replicas)
            .and_then(|n| n.checked_mul(rounds))
            .is_none_or(|n| n > 50_000_000)
    {
        return Err(api::error("Monte Carlo site-replica-round budget exceeded"));
    }
    Ok(())
}
fn wl_controls(
    final_f: f64,
    flatness: f64,
    budget: usize,
) -> Result<(), OpError> {
    if final_f <= 0.0
        || final_f >= 1.0
        || flatness <= 0.0
        || flatness >= 1.0
        || budget == 0
        || budget > 50_000_000
    {
        return Err(api::error(
            "Wang-Landau needs0<ln_f_final,flatness<1 and1..50millionproposal budget",
        ));
    }
    Ok(())
}
fn dos(
    f: &Frame,
    rows: Vec<(f64, f64)>,
) -> Result<Dynamic, OpError> {
    let (e, g): (Vec<_>, Vec<_>) = rows.into_iter().unzip();
    Ok(Dynamic(api::record(vec![
        ("energy", f.array(e, ENERGY)?),
        ("log_density", Output::encode(g)?),
    ])))
}
impl<const D: usize> NativeType for Lattice<D> {
    const KIND: &'static str = match D {
        | 1 => "rssn/native/sim/spin/lattice/d1/v1",
        | 2 => "rssn/native/sim/spin/lattice/d2/v1",
        | 3 => "rssn/native/sim/spin/lattice/d3/v1",
        | _ => "rssn/native/sim/spin/unsupported/v1",
    };
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.0.sites() > 4096 {
            return Err(api::error("lattice storage bound exceeded"));
        }
        Ok(())
    }
}
impl<const D: usize, const N: usize> NativeType for On<D, N> {
    const KIND: &'static str = match (D, N) {
        | (1, 1) => "rssn/native/sim/spin/on-model/d1/n1/v1",
        | (1, 2) => "rssn/native/sim/spin/on-model/d1/n2/v1",
        | (1, 3) => "rssn/native/sim/spin/on-model/d1/n3/v1",
        | (1, 4) => "rssn/native/sim/spin/on-model/d1/n4/v1",
        | (2, 1) => "rssn/native/sim/spin/on-model/d2/n1/v1",
        | (2, 2) => "rssn/native/sim/spin/on-model/d2/n2/v1",
        | (2, 3) => "rssn/native/sim/spin/on-model/d2/n3/v1",
        | (2, 4) => "rssn/native/sim/spin/on-model/d2/n4/v1",
        | (3, 1) => "rssn/native/sim/spin/on-model/d3/n1/v1",
        | (3, 2) => "rssn/native/sim/spin/on-model/d3/n2/v1",
        | (3, 3) => "rssn/native/sim/spin/on-model/d3/n3/v1",
        | (3, 4) => "rssn/native/sim/spin/on-model/d3/n4/v1",
        | _ => "rssn/native/sim/spin/unsupported/v1",
    };
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.lattice.sites() > 4096
            || self.native.spins.len() != self.native.lattice.sites()
        {
            return Err(api::error("O(N) spin shape/count invalid"));
        }
        finite(self.native.spins.iter().flatten().copied())
    }
}
impl<const D: usize> NativeType for Potts<D> {
    const KIND: &'static str = match D {
        | 1 => "rssn/native/sim/spin/potts/d1/v1",
        | 2 => "rssn/native/sim/spin/potts/d2/v1",
        | 3 => "rssn/native/sim/spin/potts/d3/v1",
        | _ => "rssn/native/sim/spin/unsupported/v1",
    };
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.lattice.sites() > 4096 || self.native.q > 256 {
            return Err(api::error("Potts lattice/q bound exceeded"));
        }
        Ok(())
    }
}
impl<const D: usize, const N: usize> NativeType for Temper<D, N> {
    const KIND: &'static str = match (D, N) {
        | (1, 1) => "rssn/native/sim/spin/tempering/d1/n1/v1",
        | (1, 2) => "rssn/native/sim/spin/tempering/d1/n2/v1",
        | (1, 3) => "rssn/native/sim/spin/tempering/d1/n3/v1",
        | (1, 4) => "rssn/native/sim/spin/tempering/d1/n4/v1",
        | (2, 1) => "rssn/native/sim/spin/tempering/d2/n1/v1",
        | (2, 2) => "rssn/native/sim/spin/tempering/d2/n2/v1",
        | (2, 3) => "rssn/native/sim/spin/tempering/d2/n3/v1",
        | (2, 4) => "rssn/native/sim/spin/tempering/d2/n4/v1",
        | (3, 1) => "rssn/native/sim/spin/tempering/d3/n1/v1",
        | (3, 2) => "rssn/native/sim/spin/tempering/d3/n2/v1",
        | (3, 3) => "rssn/native/sim/spin/tempering/d3/n3/v1",
        | (3, 4) => "rssn/native/sim/spin/tempering/d3/n4/v1",
        | _ => "rssn/native/sim/spin/unsupported/v1",
    };
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.replicas.len() > 32 {
            return Err(api::error("replica count bound exceeded"));
        }
        Ok(())
    }
}
impl<const D: usize> NativeType for PottsTemper<D> {
    const KIND: &'static str = match D {
        | 1 => "rssn/native/sim/spin/potts-tempering/d1/v1",
        | 2 => "rssn/native/sim/spin/potts-tempering/d2/v1",
        | 3 => "rssn/native/sim/spin/potts-tempering/d3/v1",
        | _ => "rssn/native/sim/spin/unsupported/v1",
    };
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.replicas.len() > 32 {
            return Err(api::error("replica count bound exceeded"));
        }
        Ok(())
    }
}
struct OnMeasure(CallbackRef);
struct PottsMeasure(CallbackRef);
fn measure_type(potts: bool) -> TypeExpr {
    let e = DimExpr::Var("E".into());
    let fields = vec![
        ("energy".into(), TypeExpr::Float.with_dim(e.clone())),
        ("magnetization".into(), TypeExpr::Float),
        (
            "spins".into(),
            api::tensor_type(
                if potts { 1 } else { 2 },
                if potts {
                    TypeExpr::Int
                } else {
                    TypeExpr::Float
                },
            ),
        ),
        (
            "inverse_temperature".into(),
            TypeExpr::Float.with_dim(DimExpr::Inverse(Box::new(e.clone()))),
        ),
        ("coupling".into(), TypeExpr::Float.with_dim(e)),
    ];
    TypeExpr::Function(
        vec![TypeExpr::Int, api::tensor_type(1, TypeExpr::Record(fields))],
        Box::new(TypeExpr::Bool),
    )
}
macro_rules! measure_arg {
    ($t:ident,$potts:expr) => {
        impl NativeArg for $t {
            fn ty() -> TypeExpr {
                measure_type($potts)
            }

            fn decode(
                v: &Value,
                _: usize,
                c: &CallCtx<'_>,
            ) -> Result<Self, OpError> {
                let Value::Function(h) = v else {
                    return Err(api::error("expected replica measurement callback"));
                };
                if c.callback_signature(*h)?.params.len() != 2 {
                    return Err(api::error(
                        "measurement callback needsround andreplica snapshots",
                    ));
                }
                Ok(Self(*h))
            }
        }
    };
}
measure_arg!(OnMeasure, false);
measure_arg!(PottsMeasure, true);
struct StopMeasure;
fn measurement(
    ctx: &CallCtx<'_>,
    h: CallbackRef,
    round: usize,
    snapshots: Vec<Value>,
    ty: TypeExpr,
) -> Result<(), OpError> {
    let input = Value::Tensor(TensorValue {
        shape: vec![snapshots.len()],
        data: snapshots,
        elem: ty,
    });
    match ctx.invoke_callback(h, &[Value::int(round), input])? {
        | Value::Bool(true) => Ok(()),
        | Value::Bool(false) => Err(api::error("replica measurement requested cancellation")),
        | _ => Err(api::error("measurement must return Bool acknowledgement")),
    }
}
fn on_snapshot<const D: usize, const N: usize>(
    m: &spin::OnModel<D, N>,
    f: &Frame,
) -> Result<Value, OpError> {
    let spins = Value::Tensor(TensorValue {
        shape: vec![m.spins.len(), N],
        data: m
            .spins
            .iter()
            .flatten()
            .copied()
            .map(Value::Float)
            .collect(),
        elem: TypeExpr::Float,
    });
    Ok(api::record(vec![
        ("energy", f.output(m.energy(), ENERGY, false)?),
        (
            "magnetization",
            Value::Float(m.abs_magnetisation_per_site()),
        ),
        ("spins", spins),
        (
            "inverse_temperature",
            f.output(m.beta, [-1, -2, 2, 0], false)?,
        ),
        ("coupling", f.output(m.j, ENERGY, false)?),
    ]))
}
fn potts_snapshot<const D: usize>(
    m: &spin::Potts<D>,
    f: &Frame,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("energy", f.output(m.energy(), ENERGY, false)?),
        ("magnetization", Value::Float(m.order_parameter())),
        ("spins", Output::encode(m.spins.clone())?),
        (
            "inverse_temperature",
            f.output(m.beta, [-1, -2, 2, 0], false)?,
        ),
        ("coupling", f.output(m.j, ENERGY, false)?),
    ]))
}
fn measure_schema(
    ctx: &CallCtx<'_>,
    h: CallbackRef,
) -> Result<TypeExpr, OpError> {
    match &ctx.callback_signature(h)?.params[1] {
        | TypeExpr::Tensor { elem, .. } => Ok((**elem).clone()),
        | _ => {
            Err(api::error(
                "measurement replica parameter must be rank1record tensor",
            ))
        },
    }
}
object_binding!(spin_lattice_new_d1,"num_native_spin_lattice_new_d1","crate::sim::spin_models::Lattice<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(dims:[usize;1],periodic:bool)->NewObject<Lattice<1>>,|_ctx:&CallCtx<'_>,dims:[usize;1],periodic:bool|{sites(dims,periodic)?;Ok(NewObject::normalized(Lattice(spin::Lattice::new(dims,periodic).map_err(api::error)?)))});
object_binding!(spin_lattice_sites_d1,"num_native_spin_lattice_sites_d1","crate::sim::spin_models::Lattice<D>::sites","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<1>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<1>>|{lattice.with(ctx,|s|Ok(s.0.sites()))});
object_binding!(spin_lattice_bonds_d1,"num_native_spin_lattice_bonds_d1","crate::sim::spin_models::Lattice<D>::bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<1>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<1>>|{lattice.with(ctx,|s|Ok(s.0.bonds()))});
object_binding!(spin_lattice_diagonal_bonds_d1,"num_native_spin_lattice_diagonal_bonds_d1","crate::sim::spin_models::Lattice<D>::diagonal_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<1>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<1>>|{lattice.with(ctx,|s|Ok(s.0.diagonal_bonds()))});
object_binding!(spin_lattice_neighbours_d1,"num_native_spin_lattice_neighbours_d1","crate::sim::spin_models::Lattice<D>::neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<1>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<1>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_lattice_second_neighbours_d1,"num_native_spin_lattice_second_neighbours_d1","crate::sim::spin_models::Lattice<D>::second_neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<1>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<1>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.second_neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_potts_new_d1,"num_native_spin_potts_new_d1","crate::sim::spin_models::Potts<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,q:usize,inverse_temperature:Input,seed:u64)->NewObject<Potts<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,q:usize,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::Potts::new(dims,periodic,q,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Potts{native,frame:f},dim))});
object_binding!(spin_potts_energy_d1,"num_native_spin_potts_energy_d1","crate::sim::spin_models::Potts<D>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_potts_order_parameter_d1,"num_native_spin_potts_order_parameter_d1","crate::sim::spin_models::Potts<D>::order_parameter","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{Ok(s.native.order_parameter())})});
object_binding!(spin_potts_metropolis_sweep_d1,"num_native_spin_potts_metropolis_sweep_d1","crate::sim::spin_models::Potts<D>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{Ok(s.native.metropolis_sweep())})});
object_binding!(spin_potts_heat_bath_sweep_d1,"num_native_spin_potts_heat_bath_sweep_d1","crate::sim::spin_models::Potts<D>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep();Ok(true)})});
object_binding!(spin_potts_wolff_step_d1,"num_native_spin_potts_wolff_step_d1","crate::sim::spin_models::Potts<D>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts Wolff needsferromagnetic coupling"));}Ok(s.native.wolff_step())})});
object_binding!(spin_potts_swendsen_wang_step_d1,"num_native_spin_potts_swendsen_wang_step_d1","crate::sim::spin_models::Potts<D>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts cluster updatesneedferromagnetic coupling"));}Ok(s.native.swendsen_wang_step())})});
object_binding!(spin_potts_update_d1,"num_native_spin_potts_update_d1","crate::sim::spin_models::Potts<D>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<1>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<1>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{if matches!(u,spin::Update::Wolff|spin::Update::SwendsenWang)&&s.native.j<0.0{return Err(api::error("cluster updatesneedferromagnetic coupling"));}s.native.update(u);s.validate()?;Ok(true)})});
object_binding!(spin_potts_tempering_new_d1,"num_native_spin_potts_tempering_new_d1","crate::sim::spin_models::PottsTempering<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64)->NewObject<PottsTemper<1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::PottsTempering::new(dims,periodic,q,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(PottsTemper{native,frame:f},dim))});
object_binding!(spin_on_new_d1_n1,"num_native_spin_on_new_d1_n1","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<1,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d1_n1,"num_native_spin_on_cold_start_d1_n1","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d1_n1,"num_native_spin_on_energy_d1_n1","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d1_n1,"num_native_spin_on_magnetisation_d1_n1","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->[f64;1],|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d1_n1,"num_native_spin_on_abs_magnetisation_per_site_d1_n1","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d1_n1,"num_native_spin_on_heat_bath_sweep_d1_n1","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d1_n1,"num_native_spin_on_overrelax_sweep_d1_n1","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d1_n1,"num_native_spin_on_wolff_step_d1_n1","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d1_n1,"num_native_spin_on_swendsen_wang_step_d1_n1","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d1_n1,"num_native_spin_on_local_field_d1_n1","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d1_n1,"num_native_spin_on_metropolis_sweep_d1_n1","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d1_n1,"num_native_spin_on_update_d1_n1","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d1_n1,"num_native_spin_on_couplings_d1_n1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,1>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,1>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d1_n1,"num_native_spin_tempering_new_d1_n1","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<1,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d1_n1,"num_native_spin_tempering_run_d1_n1","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<1,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<1,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d1_n1,"num_native_spin_wang_landau_on_d1_n1","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<1,1>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d1_n1,"num_native_spin_replica_exchange_wang_landau_d1_n1","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<1,1>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d1_n2,"num_native_spin_on_new_d1_n2","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<1,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d1_n2,"num_native_spin_on_cold_start_d1_n2","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d1_n2,"num_native_spin_on_energy_d1_n2","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d1_n2,"num_native_spin_on_magnetisation_d1_n2","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->[f64;2],|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d1_n2,"num_native_spin_on_abs_magnetisation_per_site_d1_n2","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d1_n2,"num_native_spin_on_heat_bath_sweep_d1_n2","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d1_n2,"num_native_spin_on_overrelax_sweep_d1_n2","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d1_n2,"num_native_spin_on_wolff_step_d1_n2","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d1_n2,"num_native_spin_on_swendsen_wang_step_d1_n2","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d1_n2,"num_native_spin_on_local_field_d1_n2","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d1_n2,"num_native_spin_on_metropolis_sweep_d1_n2","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d1_n2,"num_native_spin_on_update_d1_n2","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d1_n2,"num_native_spin_on_couplings_d1_n2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,2>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,2>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d1_n2,"num_native_spin_tempering_new_d1_n2","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<1,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d1_n2,"num_native_spin_tempering_run_d1_n2","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<1,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<1,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d1_n2,"num_native_spin_wang_landau_on_d1_n2","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<1,2>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d1_n2,"num_native_spin_replica_exchange_wang_landau_d1_n2","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<1,2>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d1_n3,"num_native_spin_on_new_d1_n3","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<1,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d1_n3,"num_native_spin_on_cold_start_d1_n3","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d1_n3,"num_native_spin_on_energy_d1_n3","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d1_n3,"num_native_spin_on_magnetisation_d1_n3","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->[f64;3],|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d1_n3,"num_native_spin_on_abs_magnetisation_per_site_d1_n3","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d1_n3,"num_native_spin_on_heat_bath_sweep_d1_n3","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d1_n3,"num_native_spin_on_overrelax_sweep_d1_n3","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d1_n3,"num_native_spin_on_wolff_step_d1_n3","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d1_n3,"num_native_spin_on_swendsen_wang_step_d1_n3","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d1_n3,"num_native_spin_on_local_field_d1_n3","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d1_n3,"num_native_spin_on_metropolis_sweep_d1_n3","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d1_n3,"num_native_spin_on_update_d1_n3","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d1_n3,"num_native_spin_on_couplings_d1_n3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,3>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,3>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d1_n3,"num_native_spin_tempering_new_d1_n3","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<1,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d1_n3,"num_native_spin_tempering_run_d1_n3","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<1,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<1,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d1_n3,"num_native_spin_wang_landau_on_d1_n3","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<1,3>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d1_n3,"num_native_spin_replica_exchange_wang_landau_d1_n3","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<1,3>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d1_n4,"num_native_spin_on_new_d1_n4","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<1,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d1_n4,"num_native_spin_on_cold_start_d1_n4","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d1_n4,"num_native_spin_on_energy_d1_n4","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d1_n4,"num_native_spin_on_magnetisation_d1_n4","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->[f64;4],|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d1_n4,"num_native_spin_on_abs_magnetisation_per_site_d1_n4","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d1_n4,"num_native_spin_on_heat_bath_sweep_d1_n4","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d1_n4,"num_native_spin_on_overrelax_sweep_d1_n4","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d1_n4,"num_native_spin_on_wolff_step_d1_n4","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d1_n4,"num_native_spin_on_swendsen_wang_step_d1_n4","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d1_n4,"num_native_spin_on_local_field_d1_n4","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d1_n4,"num_native_spin_on_metropolis_sweep_d1_n4","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d1_n4,"num_native_spin_on_update_d1_n4","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d1_n4,"num_native_spin_on_couplings_d1_n4","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<1,4>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<1,4>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d1_n4,"num_native_spin_tempering_new_d1_n4","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<1,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d1_n4,"num_native_spin_tempering_run_d1_n4","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<1,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<1,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d1_n4,"num_native_spin_wang_landau_on_d1_n4","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<1,4>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d1_n4,"num_native_spin_replica_exchange_wang_landau_d1_n4","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<1,4>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_potts_tempering_run_d1,"num_native_spin_potts_tempering_run_d1","crate::sim::spin_models::PottsTempering<D>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<PottsTemper<1>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<PottsTemper<1>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|potts_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}
if let Err(e)=result{std::panic::resume_unwind(e);}s.validate()?;Ok(true)})});
object_binding!(spin_enumerate_ising_dos_d1,"num_native_spin_enumerate_ising_dos_d1","crate::sim::spin_models::enumerate_ising_dos","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool|{if sites(dims,periodic)?>20{return Err(api::error("exact Ising enumeration limitedto20sites"));}let f=frame(ctx,&f)?;let out=spin::enumerate_ising_dos(dims,periodic).map_err(api::error)?;let(e,c):(Vec<_>,Vec<_>)=out.into_iter().unzip();Ok(Dynamic(api::record(vec![("energy",f.array(e.into_iter().map(|e|e as f64).collect(),ENERGY)?),("count",Output::encode(c)?)])))});
object_binding!(spin_wang_landau_ising_d1,"num_native_spin_wang_landau_ising_d1","crate::sim::spin_models::wang_landau_ising","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("Wang Landau lattice bound256"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_ising::<1>(dims,periodic,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out.into_iter().map(|(e,g)|(e as f64,g)).collect())});
object_binding!(spin_lattice_new_d2,"num_native_spin_lattice_new_d2","crate::sim::spin_models::Lattice<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(dims:[usize;2],periodic:bool)->NewObject<Lattice<2>>,|_ctx:&CallCtx<'_>,dims:[usize;2],periodic:bool|{sites(dims,periodic)?;Ok(NewObject::normalized(Lattice(spin::Lattice::new(dims,periodic).map_err(api::error)?)))});
object_binding!(spin_lattice_sites_d2,"num_native_spin_lattice_sites_d2","crate::sim::spin_models::Lattice<D>::sites","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<2>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<2>>|{lattice.with(ctx,|s|Ok(s.0.sites()))});
object_binding!(spin_lattice_bonds_d2,"num_native_spin_lattice_bonds_d2","crate::sim::spin_models::Lattice<D>::bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<2>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<2>>|{lattice.with(ctx,|s|Ok(s.0.bonds()))});
object_binding!(spin_lattice_diagonal_bonds_d2,"num_native_spin_lattice_diagonal_bonds_d2","crate::sim::spin_models::Lattice<D>::diagonal_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<2>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<2>>|{lattice.with(ctx,|s|Ok(s.0.diagonal_bonds()))});
object_binding!(spin_lattice_neighbours_d2,"num_native_spin_lattice_neighbours_d2","crate::sim::spin_models::Lattice<D>::neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<2>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<2>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_lattice_second_neighbours_d2,"num_native_spin_lattice_second_neighbours_d2","crate::sim::spin_models::Lattice<D>::second_neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<2>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<2>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.second_neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_potts_new_d2,"num_native_spin_potts_new_d2","crate::sim::spin_models::Potts<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,q:usize,inverse_temperature:Input,seed:u64)->NewObject<Potts<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,q:usize,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::Potts::new(dims,periodic,q,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Potts{native,frame:f},dim))});
object_binding!(spin_potts_energy_d2,"num_native_spin_potts_energy_d2","crate::sim::spin_models::Potts<D>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_potts_order_parameter_d2,"num_native_spin_potts_order_parameter_d2","crate::sim::spin_models::Potts<D>::order_parameter","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{Ok(s.native.order_parameter())})});
object_binding!(spin_potts_metropolis_sweep_d2,"num_native_spin_potts_metropolis_sweep_d2","crate::sim::spin_models::Potts<D>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{Ok(s.native.metropolis_sweep())})});
object_binding!(spin_potts_heat_bath_sweep_d2,"num_native_spin_potts_heat_bath_sweep_d2","crate::sim::spin_models::Potts<D>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep();Ok(true)})});
object_binding!(spin_potts_wolff_step_d2,"num_native_spin_potts_wolff_step_d2","crate::sim::spin_models::Potts<D>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts Wolff needsferromagnetic coupling"));}Ok(s.native.wolff_step())})});
object_binding!(spin_potts_swendsen_wang_step_d2,"num_native_spin_potts_swendsen_wang_step_d2","crate::sim::spin_models::Potts<D>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts cluster updatesneedferromagnetic coupling"));}Ok(s.native.swendsen_wang_step())})});
object_binding!(spin_potts_update_d2,"num_native_spin_potts_update_d2","crate::sim::spin_models::Potts<D>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<2>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<2>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{if matches!(u,spin::Update::Wolff|spin::Update::SwendsenWang)&&s.native.j<0.0{return Err(api::error("cluster updatesneedferromagnetic coupling"));}s.native.update(u);s.validate()?;Ok(true)})});
object_binding!(spin_potts_tempering_new_d2,"num_native_spin_potts_tempering_new_d2","crate::sim::spin_models::PottsTempering<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64)->NewObject<PottsTemper<2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::PottsTempering::new(dims,periodic,q,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(PottsTemper{native,frame:f},dim))});
object_binding!(spin_on_new_d2_n1,"num_native_spin_on_new_d2_n1","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<2,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d2_n1,"num_native_spin_on_cold_start_d2_n1","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d2_n1,"num_native_spin_on_energy_d2_n1","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d2_n1,"num_native_spin_on_magnetisation_d2_n1","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->[f64;1],|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d2_n1,"num_native_spin_on_abs_magnetisation_per_site_d2_n1","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d2_n1,"num_native_spin_on_heat_bath_sweep_d2_n1","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d2_n1,"num_native_spin_on_overrelax_sweep_d2_n1","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d2_n1,"num_native_spin_on_wolff_step_d2_n1","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d2_n1,"num_native_spin_on_swendsen_wang_step_d2_n1","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d2_n1,"num_native_spin_on_local_field_d2_n1","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d2_n1,"num_native_spin_on_metropolis_sweep_d2_n1","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d2_n1,"num_native_spin_on_update_d2_n1","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d2_n1,"num_native_spin_on_couplings_d2_n1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,1>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,1>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d2_n1,"num_native_spin_tempering_new_d2_n1","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<2,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d2_n1,"num_native_spin_tempering_run_d2_n1","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<2,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<2,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d2_n1,"num_native_spin_wang_landau_on_d2_n1","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<2,1>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d2_n1,"num_native_spin_replica_exchange_wang_landau_d2_n1","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<2,1>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d2_n2,"num_native_spin_on_new_d2_n2","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<2,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d2_n2,"num_native_spin_on_cold_start_d2_n2","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d2_n2,"num_native_spin_on_energy_d2_n2","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d2_n2,"num_native_spin_on_magnetisation_d2_n2","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->[f64;2],|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d2_n2,"num_native_spin_on_abs_magnetisation_per_site_d2_n2","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d2_n2,"num_native_spin_on_heat_bath_sweep_d2_n2","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d2_n2,"num_native_spin_on_overrelax_sweep_d2_n2","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d2_n2,"num_native_spin_on_wolff_step_d2_n2","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d2_n2,"num_native_spin_on_swendsen_wang_step_d2_n2","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d2_n2,"num_native_spin_on_local_field_d2_n2","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d2_n2,"num_native_spin_on_metropolis_sweep_d2_n2","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d2_n2,"num_native_spin_on_update_d2_n2","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d2_n2,"num_native_spin_on_couplings_d2_n2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,2>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,2>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d2_n2,"num_native_spin_tempering_new_d2_n2","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<2,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d2_n2,"num_native_spin_tempering_run_d2_n2","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<2,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<2,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d2_n2,"num_native_spin_wang_landau_on_d2_n2","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<2,2>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d2_n2,"num_native_spin_replica_exchange_wang_landau_d2_n2","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<2,2>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d2_n3,"num_native_spin_on_new_d2_n3","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<2,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d2_n3,"num_native_spin_on_cold_start_d2_n3","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d2_n3,"num_native_spin_on_energy_d2_n3","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d2_n3,"num_native_spin_on_magnetisation_d2_n3","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->[f64;3],|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d2_n3,"num_native_spin_on_abs_magnetisation_per_site_d2_n3","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d2_n3,"num_native_spin_on_heat_bath_sweep_d2_n3","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d2_n3,"num_native_spin_on_overrelax_sweep_d2_n3","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d2_n3,"num_native_spin_on_wolff_step_d2_n3","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d2_n3,"num_native_spin_on_swendsen_wang_step_d2_n3","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d2_n3,"num_native_spin_on_local_field_d2_n3","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d2_n3,"num_native_spin_on_metropolis_sweep_d2_n3","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d2_n3,"num_native_spin_on_update_d2_n3","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d2_n3,"num_native_spin_on_couplings_d2_n3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,3>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,3>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d2_n3,"num_native_spin_tempering_new_d2_n3","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<2,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d2_n3,"num_native_spin_tempering_run_d2_n3","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<2,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<2,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d2_n3,"num_native_spin_wang_landau_on_d2_n3","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<2,3>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d2_n3,"num_native_spin_replica_exchange_wang_landau_d2_n3","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<2,3>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d2_n4,"num_native_spin_on_new_d2_n4","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<2,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d2_n4,"num_native_spin_on_cold_start_d2_n4","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d2_n4,"num_native_spin_on_energy_d2_n4","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d2_n4,"num_native_spin_on_magnetisation_d2_n4","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->[f64;4],|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d2_n4,"num_native_spin_on_abs_magnetisation_per_site_d2_n4","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d2_n4,"num_native_spin_on_heat_bath_sweep_d2_n4","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d2_n4,"num_native_spin_on_overrelax_sweep_d2_n4","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d2_n4,"num_native_spin_on_wolff_step_d2_n4","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d2_n4,"num_native_spin_on_swendsen_wang_step_d2_n4","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d2_n4,"num_native_spin_on_local_field_d2_n4","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d2_n4,"num_native_spin_on_metropolis_sweep_d2_n4","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d2_n4,"num_native_spin_on_update_d2_n4","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d2_n4,"num_native_spin_on_couplings_d2_n4","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<2,4>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<2,4>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d2_n4,"num_native_spin_tempering_new_d2_n4","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<2,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d2_n4,"num_native_spin_tempering_run_d2_n4","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<2,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<2,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d2_n4,"num_native_spin_wang_landau_on_d2_n4","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<2,4>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d2_n4,"num_native_spin_replica_exchange_wang_landau_d2_n4","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<2,4>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_potts_tempering_run_d2,"num_native_spin_potts_tempering_run_d2","crate::sim::spin_models::PottsTempering<D>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<PottsTemper<2>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<PottsTemper<2>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|potts_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}
if let Err(e)=result{std::panic::resume_unwind(e);}s.validate()?;Ok(true)})});
object_binding!(spin_enumerate_ising_dos_d2,"num_native_spin_enumerate_ising_dos_d2","crate::sim::spin_models::enumerate_ising_dos","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool|{if sites(dims,periodic)?>20{return Err(api::error("exact Ising enumeration limitedto20sites"));}let f=frame(ctx,&f)?;let out=spin::enumerate_ising_dos(dims,periodic).map_err(api::error)?;let(e,c):(Vec<_>,Vec<_>)=out.into_iter().unzip();Ok(Dynamic(api::record(vec![("energy",f.array(e.into_iter().map(|e|e as f64).collect(),ENERGY)?),("count",Output::encode(c)?)])))});
object_binding!(spin_wang_landau_ising_d2,"num_native_spin_wang_landau_ising_d2","crate::sim::spin_models::wang_landau_ising","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("Wang Landau lattice bound256"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_ising::<2>(dims,periodic,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out.into_iter().map(|(e,g)|(e as f64,g)).collect())});
object_binding!(spin_lattice_new_d3,"num_native_spin_lattice_new_d3","crate::sim::spin_models::Lattice<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(dims:[usize;3],periodic:bool)->NewObject<Lattice<3>>,|_ctx:&CallCtx<'_>,dims:[usize;3],periodic:bool|{sites(dims,periodic)?;Ok(NewObject::normalized(Lattice(spin::Lattice::new(dims,periodic).map_err(api::error)?)))});
object_binding!(spin_lattice_sites_d3,"num_native_spin_lattice_sites_d3","crate::sim::spin_models::Lattice<D>::sites","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<3>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<3>>|{lattice.with(ctx,|s|Ok(s.0.sites()))});
object_binding!(spin_lattice_bonds_d3,"num_native_spin_lattice_bonds_d3","crate::sim::spin_models::Lattice<D>::bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<3>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<3>>|{lattice.with(ctx,|s|Ok(s.0.bonds()))});
object_binding!(spin_lattice_diagonal_bonds_d3,"num_native_spin_lattice_diagonal_bonds_d3","crate::sim::spin_models::Lattice<D>::diagonal_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<3>>)->usize,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<3>>|{lattice.with(ctx,|s|Ok(s.0.diagonal_bonds()))});
object_binding!(spin_lattice_neighbours_d3,"num_native_spin_lattice_neighbours_d3","crate::sim::spin_models::Lattice<D>::neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<3>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<3>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_lattice_second_neighbours_d3,"num_native_spin_lattice_second_neighbours_d3","crate::sim::spin_models::Lattice<D>::second_neighbours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(lattice:Handle<Lattice<3>>,index:usize)->Neighbours,|ctx:&CallCtx<'_>,lattice:Handle<Lattice<3>>,index:usize|{lattice.with(ctx,|s|{if index>=s.0.sites(){return Err(api::error("spin lattice index outside lattice"));}Ok(Neighbours(s.0.second_neighbours(index).iter().map(|i|(*i!=usize::MAX).then_some(*i)).collect()))})});
object_binding!(spin_potts_new_d3,"num_native_spin_potts_new_d3","crate::sim::spin_models::Potts<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,q:usize,inverse_temperature:Input,seed:u64)->NewObject<Potts<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,q:usize,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::Potts::new(dims,periodic,q,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Potts{native,frame:f},dim))});
object_binding!(spin_potts_energy_d3,"num_native_spin_potts_energy_d3","crate::sim::spin_models::Potts<D>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_potts_order_parameter_d3,"num_native_spin_potts_order_parameter_d3","crate::sim::spin_models::Potts<D>::order_parameter","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{Ok(s.native.order_parameter())})});
object_binding!(spin_potts_metropolis_sweep_d3,"num_native_spin_potts_metropolis_sweep_d3","crate::sim::spin_models::Potts<D>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->f64,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{Ok(s.native.metropolis_sweep())})});
object_binding!(spin_potts_heat_bath_sweep_d3,"num_native_spin_potts_heat_bath_sweep_d3","crate::sim::spin_models::Potts<D>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep();Ok(true)})});
object_binding!(spin_potts_wolff_step_d3,"num_native_spin_potts_wolff_step_d3","crate::sim::spin_models::Potts<D>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts Wolff needsferromagnetic coupling"));}Ok(s.native.wolff_step())})});
object_binding!(spin_potts_swendsen_wang_step_d3,"num_native_spin_potts_swendsen_wang_step_d3","crate::sim::spin_models::Potts<D>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>|{model.with_mut(ctx,|s|{if s.native.j<0.0{return Err(api::error("Potts cluster updatesneedferromagnetic coupling"));}Ok(s.native.swendsen_wang_step())})});
object_binding!(spin_potts_update_d3,"num_native_spin_potts_update_d3","crate::sim::spin_models::Potts<D>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<Potts<3>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<Potts<3>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{if matches!(u,spin::Update::Wolff|spin::Update::SwendsenWang)&&s.native.j<0.0{return Err(api::error("cluster updatesneedferromagnetic coupling"));}s.native.update(u);s.validate()?;Ok(true)})});
object_binding!(spin_potts_tempering_new_d3,"num_native_spin_potts_tempering_new_d3","crate::sim::spin_models::PottsTempering<D>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64)->NewObject<PottsTemper<3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,q:usize,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;if !(2..=256).contains(&q){return Err(api::error("Potts q must be2..256"));}let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::PottsTempering::new(dims,periodic,q,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(PottsTemper{native,frame:f},dim))});
object_binding!(spin_on_new_d3_n1,"num_native_spin_on_new_d3_n1","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<3,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d3_n1,"num_native_spin_on_cold_start_d3_n1","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d3_n1,"num_native_spin_on_energy_d3_n1","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d3_n1,"num_native_spin_on_magnetisation_d3_n1","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->[f64;1],|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d3_n1,"num_native_spin_on_abs_magnetisation_per_site_d3_n1","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d3_n1,"num_native_spin_on_heat_bath_sweep_d3_n1","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d3_n1,"num_native_spin_on_overrelax_sweep_d3_n1","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d3_n1,"num_native_spin_on_wolff_step_d3_n1","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d3_n1,"num_native_spin_on_swendsen_wang_step_d3_n1","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d3_n1,"num_native_spin_on_local_field_d3_n1","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d3_n1,"num_native_spin_on_metropolis_sweep_d3_n1","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d3_n1,"num_native_spin_on_update_d3_n1","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d3_n1,"num_native_spin_on_couplings_d3_n1","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,1>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,1>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d3_n1,"num_native_spin_tempering_new_d3_n1","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<3,1>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d3_n1,"num_native_spin_tempering_run_d3_n1","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<3,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<3,1>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d3_n1,"num_native_spin_wang_landau_on_d3_n1","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<3,1>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d3_n1,"num_native_spin_replica_exchange_wang_landau_d3_n1","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<3,1>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d3_n2,"num_native_spin_on_new_d3_n2","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<3,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d3_n2,"num_native_spin_on_cold_start_d3_n2","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d3_n2,"num_native_spin_on_energy_d3_n2","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d3_n2,"num_native_spin_on_magnetisation_d3_n2","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->[f64;2],|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d3_n2,"num_native_spin_on_abs_magnetisation_per_site_d3_n2","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d3_n2,"num_native_spin_on_heat_bath_sweep_d3_n2","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d3_n2,"num_native_spin_on_overrelax_sweep_d3_n2","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d3_n2,"num_native_spin_on_wolff_step_d3_n2","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d3_n2,"num_native_spin_on_swendsen_wang_step_d3_n2","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d3_n2,"num_native_spin_on_local_field_d3_n2","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d3_n2,"num_native_spin_on_metropolis_sweep_d3_n2","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d3_n2,"num_native_spin_on_update_d3_n2","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d3_n2,"num_native_spin_on_couplings_d3_n2","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,2>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,2>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d3_n2,"num_native_spin_tempering_new_d3_n2","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<3,2>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d3_n2,"num_native_spin_tempering_run_d3_n2","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<3,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<3,2>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d3_n2,"num_native_spin_wang_landau_on_d3_n2","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<3,2>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d3_n2,"num_native_spin_replica_exchange_wang_landau_d3_n2","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<3,2>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d3_n3,"num_native_spin_on_new_d3_n3","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<3,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d3_n3,"num_native_spin_on_cold_start_d3_n3","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d3_n3,"num_native_spin_on_energy_d3_n3","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d3_n3,"num_native_spin_on_magnetisation_d3_n3","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->[f64;3],|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d3_n3,"num_native_spin_on_abs_magnetisation_per_site_d3_n3","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d3_n3,"num_native_spin_on_heat_bath_sweep_d3_n3","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d3_n3,"num_native_spin_on_overrelax_sweep_d3_n3","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d3_n3,"num_native_spin_on_wolff_step_d3_n3","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d3_n3,"num_native_spin_on_swendsen_wang_step_d3_n3","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d3_n3,"num_native_spin_on_local_field_d3_n3","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d3_n3,"num_native_spin_on_metropolis_sweep_d3_n3","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d3_n3,"num_native_spin_on_update_d3_n3","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d3_n3,"num_native_spin_on_couplings_d3_n3","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,3>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,3>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d3_n3,"num_native_spin_tempering_new_d3_n3","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<3,3>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d3_n3,"num_native_spin_tempering_run_d3_n3","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<3,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<3,3>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d3_n3,"num_native_spin_wang_landau_on_d3_n3","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<3,3>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d3_n3,"num_native_spin_replica_exchange_wang_landau_d3_n3","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<3,3>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_on_new_d3_n4,"num_native_spin_on_new_d3_n4","crate::sim::spin_models::OnModel<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64)->NewObject<On<3,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperature:Input,seed:u64|{sites(dims,periodic)?;let f=frame(ctx,&f)?;let beta=beta(&f,&inverse_temperature.0)?;let native=spin::OnModel::new(dims,periodic,beta,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(On{native,frame:f},dim))});
object_binding!(spin_on_cold_start_d3_n4,"num_native_spin_on_cold_start_d3_n4","crate::sim::spin_models::OnModel<D,N>::cold_start","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{s.native.cold_start();Ok(true)})});
object_binding!(spin_on_energy_d3_n4,"num_native_spin_on_energy_d3_n4","crate::sim::spin_models::OnModel<D,N>::energy","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{Ok(Dynamic(s.frame.output(s.native.energy(),ENERGY,false)?))})});
object_binding!(spin_on_magnetisation_d3_n4,"num_native_spin_on_magnetisation_d3_n4","crate::sim::spin_models::OnModel<D,N>::magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->[f64;4],|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{Ok(s.native.magnetisation())})});
object_binding!(spin_on_abs_magnetisation_per_site_d3_n4,"num_native_spin_on_abs_magnetisation_per_site_d3_n4","crate::sim::spin_models::OnModel<D,N>::abs_magnetisation_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{Ok(s.native.abs_magnetisation_per_site())})});
object_binding!(spin_on_heat_bath_sweep_d3_n4,"num_native_spin_on_heat_bath_sweep_d3_n4","crate::sim::spin_models::OnModel<D,N>::heat_bath_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{s.native.heat_bath_sweep().map_err(api::error)?;Ok(true)})});
object_binding!(spin_on_overrelax_sweep_d3_n4,"num_native_spin_on_overrelax_sweep_d3_n4","crate::sim::spin_models::OnModel<D,N>::overrelax_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{s.native.overrelax_sweep();Ok(true)})});
object_binding!(spin_on_wolff_step_d3_n4,"num_native_spin_on_wolff_step_d3_n4","crate::sim::spin_models::OnModel<D,N>::wolff_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{s.native.wolff_step().map_err(api::error)})});
object_binding!(spin_on_swendsen_wang_step_d3_n4,"num_native_spin_on_swendsen_wang_step_d3_n4","crate::sim::spin_models::OnModel<D,N>::swendsen_wang_step","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>)->usize,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>|{model.with_mut(ctx,|s|{s.native.swendsen_wang_step().map_err(api::error)})});
object_binding!(spin_on_local_field_d3_n4,"num_native_spin_on_local_field_d3_n4","crate::sim::spin_models::OnModel<D,N>::local_field_public","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>,index:usize)->Dynamic,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>,index:usize|{model.with(ctx,|s|{if index>=s.native.lattice.sites(){return Err(api::error("local field indexoutside lattice"));}Ok(Dynamic(s.frame.array(s.native.local_field_public(index).to_vec(),ENERGY)?))})});
object_binding!(spin_on_metropolis_sweep_d3_n4,"num_native_spin_on_metropolis_sweep_d3_n4","crate::sim::spin_models::OnModel<D,N>::metropolis_sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>,delta:f64)->f64,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>,delta:f64|{if delta<=0.0||delta>4.0{return Err(api::error("Metropolis proposal width must bein(0,4]"));}model.with_mut(ctx,|s|Ok(s.native.metropolis_sweep(delta)))});
object_binding!(spin_on_update_d3_n4,"num_native_spin_on_update_d3_n4","crate::sim::spin_models::OnModel<D,N>::update","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>,kind:String,proposal:f64)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>,kind:String,proposal:f64|{let u=update(&kind,proposal)?;model.with_mut(ctx,|s|{s.native.update(u).map_err(api::error)?;s.validate()?;Ok(true)})});
object_binding!(spin_on_couplings_d3_n4,"num_native_spin_on_couplings_d3_n4","","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(model:Handle<On<3,4>>,coupling:Input,diagonal_coupling:Input,field:Input)->bool,|ctx:&CallCtx<'_>,model:Handle<On<3,4>>,coupling:Input,diagonal_coupling:Input,field:Input|{model.with_mut(ctx,|s|{let j=s.frame.scalar(&coupling.0,ENERGY,false)?;let j2=s.frame.scalar(&diagonal_coupling.0,ENERGY,false)?;let h=s.frame.scalar(&field.0,ENERGY,false)?;if j.abs()>100.0||j2.abs()>100.0||h.abs()>100.0{return Err(api::error("normalized coupling field bounds exceed100"));}s.native.j=j;s.native.j2=j2;s.native.field=h;Ok(true)})});
object_binding!(spin_tempering_new_d3_n4,"num_native_spin_tempering_new_d3_n4","crate::sim::spin_models::ParallelTempering<D,N>::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64)->NewObject<Temper<3,4>>,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,inverse_temperatures:Input,seed:u64|{let count=sites(dims,periodic)?;let f=frame(ctx,&f)?;let betas=betas(&f,&inverse_temperatures.0)?;work(count,betas.len(),1)?;let native=spin::ParallelTempering::new(dims,periodic,&betas,seed).map_err(api::error)?;let dim=f.dims.to_vec();Ok(NewObject::new(Temper{native,frame:f},dim))});
object_binding!(spin_tempering_run_d3_n4,"num_native_spin_tempering_run_d3_n4","crate::sim::spin_models::ParallelTempering<D,N>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<Temper<3,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<Temper<3,4>>,rounds:usize,kind:String,proposal:f64,measure:OnMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|on_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}match result{Ok(out)=>out.map_err(api::error)?,Err(e)=>std::panic::resume_unwind(e)}s.validate()?;Ok(true)})});
object_binding!(spin_wang_landau_on_d3_n4,"num_native_spin_wang_landau_on_d3_n4","crate::sim::spin_models::wang_landau_on","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,window:Input,bins:usize,delta:f64,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{let n=sites(dims,periodic)?;if n>256||!(2..=4096).contains(&bins)||delta<=0.0||delta>4.0{return Err(api::error("WL lattice/bin/proposal bounds invalid"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let e=engine::array::<2>(&f,&window.0,ENERGY)?;if e[0]>=e[1]{return Err(api::error("WL energy window must increase"));}let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_on::<3,4>(dims,periodic,(e[0],e[1],bins),delta,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_replica_exchange_wang_landau_d3_n4,"num_native_spin_replica_exchange_wang_landau_d3_n4","crate::sim::spin_models::replica_exchange_wang_landau","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,options:Input,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("REWL lattice sitebound256"));}let f=frame(ctx,&f)?;let options=rewl_options(&f,&options)?;wl_controls(options.ln_f_final,options.flatness,proposal_budget)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::replica_exchange_wang_landau::<3,4>(dims,periodic,&options)).map_err(api::error)?;dos(&f,out)});
object_binding!(spin_potts_tempering_run_d3,"num_native_spin_potts_tempering_run_d3","crate::sim::spin_models::PottsTempering<D>::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(tempering:Handle<PottsTemper<3>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure)->bool,|ctx:&CallCtx<'_>,tempering:Handle<PottsTemper<3>>,rounds:usize,kind:String,proposal:f64,measure:PottsMeasure|{let u=update(&kind,proposal)?;let schema=measure_schema(ctx,measure.0)?;tempering.with_mut(ctx,|s|{let n=s.native.replicas.first().map_or(0,|r|r.lattice.sites());work(n,s.native.replicas.len(),rounds)?;let mut failure=None;let frame=s.frame.clone();let result=std::panic::catch_unwind(std::panic::AssertUnwindSafe(||s.native.run(rounds,u,|round,replicas|{let snapshot=replicas.iter().map(|m|potts_snapshot(m,&frame)).collect::<Result<Vec<_>,_>>();let result=snapshot.and_then(|snapshot|measurement(ctx,measure.0,round,snapshot,schema.clone()));if let Err(e)=result{failure=Some(e);std::panic::resume_unwind(Box::new(StopMeasure));}})));if let Some(e)=failure{return Err(e);}
if let Err(e)=result{std::panic::resume_unwind(e);}s.validate()?;Ok(true)})});
object_binding!(spin_enumerate_ising_dos_d3,"num_native_spin_enumerate_ising_dos_d3","crate::sim::spin_models::enumerate_ising_dos","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool|{if sites(dims,periodic)?>20{return Err(api::error("exact Ising enumeration limitedto20sites"));}let f=frame(ctx,&f)?;let out=spin::enumerate_ising_dos(dims,periodic).map_err(api::error)?;let(e,c):(Vec<_>,Vec<_>)=out.into_iter().unzip();Ok(Dynamic(api::record(vec![("energy",f.array(e.into_iter().map(|e|e as f64).collect(),ENERGY)?),("count",Output::encode(c)?)])))});
object_binding!(spin_wang_landau_ising_d3,"num_native_spin_wang_landau_ising_d3","crate::sim::spin_models::wang_landau_ising","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,ln_f_final:f64,flatness:f64,seed:u64,proposal_budget:usize|{if sites(dims,periodic)?>256{return Err(api::error("Wang Landau lattice bound256"));}wl_controls(ln_f_final,flatness,proposal_budget)?;let f=frame(ctx,&f)?;let out=spin::with_wang_landau_budget(proposal_budget,||spin::wang_landau_ising::<3>(dims,periodic,ln_f_final,flatness,seed)).map_err(api::error)?;dos(&f,out.into_iter().map(|(e,g)|(e as f64,g)).collect())});
fn rewl_options(
    f: &Frame,
    v: &Input,
) -> Result<spin::RewlOptions, OpError> {
    let range = engine::array::<2>(f, api::field(&v.0, "range")?, ENERGY)?;
    let integer = |n| <usize as Arg>::decode(api::field(&v.0, n)?, 0);
    let number = |n| {
        api::normalized(api::field(&v.0, n)?)?;
        crate::rules::numerical_abi::numeric(api::field(&v.0, n)?, 0)
    };
    let o = spin::RewlOptions {
        range: (range[0], range[1]),
        nbins: integer("bins")?,
        windows: integer("windows")?,
        overlap: number("overlap")?,
        walkers: integer("walkers")?,
        delta: number("proposal")?,
        ln_f_final: number("ln_f_final")?,
        flatness: number("flatness")?,
        sweeps_per_exchange: integer("sweeps_per_exchange")?,
        seed: <u64 as Arg>::decode(api::field(&v.0, "seed")?, 0)?,
    };
    if o.range.0 >= o.range.1
        || o.nbins < 2
        || o.nbins > 4096
        || o.windows == 0
        || o.windows > 16
        || o.walkers == 0
        || o.walkers > 16
        || o.windows * o.walkers > 64
        || !(0.0..1.0).contains(&o.overlap)
        || o.delta <= 0.0
        || o.delta > 4.0
        || o.sweeps_per_exchange == 0
        || o.sweeps_per_exchange > 10000
    {
        return Err(api::error("REWL options exceed domain/allocation bounds"));
    }
    Ok(o)
}
object_binding!(spin_heat_bath_direction_n1,"num_native_spin_heat_bath_direction_n1","crate::sim::spin_models::heat_bath_direction","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input)->[f64;1],|ctx:&CallCtx<'_>,f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input|{let f=frame(ctx,&f)?;let field=engine::array::<1>(&f,&field.0,ENERGY)?;let b=beta(&f,&inverse_temperature.0)?;if field.iter().any(|x|x.abs()>100.0){return Err(api::error("normalized heatbath field bound100"));}rng.with_mut(ctx,|r|spin::heat_bath_direction(r,&field,b).map_err(api::error))});
object_binding!(spin_heat_bath_direction_n2,"num_native_spin_heat_bath_direction_n2","crate::sim::spin_models::heat_bath_direction","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input)->[f64;2],|ctx:&CallCtx<'_>,f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input|{let f=frame(ctx,&f)?;let field=engine::array::<2>(&f,&field.0,ENERGY)?;let b=beta(&f,&inverse_temperature.0)?;if field.iter().any(|x|x.abs()>100.0){return Err(api::error("normalized heatbath field bound100"));}rng.with_mut(ctx,|r|spin::heat_bath_direction(r,&field,b).map_err(api::error))});
object_binding!(spin_heat_bath_direction_n3,"num_native_spin_heat_bath_direction_n3","crate::sim::spin_models::heat_bath_direction","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input)->[f64;3],|ctx:&CallCtx<'_>,f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input|{let f=frame(ctx,&f)?;let field=engine::array::<3>(&f,&field.0,ENERGY)?;let b=beta(&f,&inverse_temperature.0)?;if field.iter().any(|x|x.abs()>100.0){return Err(api::error("normalized heatbath field bound100"));}rng.with_mut(ctx,|r|spin::heat_bath_direction(r,&field,b).map_err(api::error))});
object_binding!(spin_heat_bath_direction_n4,"num_native_spin_heat_bath_direction_n4","crate::sim::spin_models::heat_bath_direction","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input)->[f64;4],|ctx:&CallCtx<'_>,f:Handle<Frame>,rng:Handle<Rng>,field:Input,inverse_temperature:Input|{let f=frame(ctx,&f)?;let field=engine::array::<4>(&f,&field.0,ENERGY)?;let b=beta(&f,&inverse_temperature.0)?;if field.iter().any(|x|x.abs()>100.0){return Err(api::error("normalized heatbath field bound100"));}rng.with_mut(ctx,|r|spin::heat_bath_direction(r,&field,b).map_err(api::error))});
object_binding!(spin_canonical_from_dos,"num_native_spin_canonical_from_dos","crate::sim::spin_models::canonical_from_dos","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,energies:Input,log_density:Vec<f64>,inverse_temperature:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,energies:Input,log_density:Vec<f64>,inverse_temperature:Input|{let f=frame(ctx,&f)?;let e=f.vector(&energies.0,ENERGY)?;if e.is_empty()||e.len()!=log_density.len(){return Err(api::error("densityofstates arrays must match nonempty"));}let b=beta(&f,&inverse_temperature.0)?;let rows=e.into_iter().zip(log_density).collect::<Vec<_>>();let(z,e,c)=spin::canonical_from_dos(&rows,b);Ok(Dynamic(api::record(vec![("log_partition",Output::encode(z)?),("mean_energy",f.output(e,ENERGY,false)?),("heat_capacity_over_kB",Output::encode(c)?)])))});
object_binding!(spin_reweight,"num_native_spin_reweight","crate::sim::spin_models::reweight","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,energies:Input,observable:Input,beta0:Input,beta1:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,energies:Input,observable:Input,beta0:Input,beta1:Input|{let f=frame(ctx,&f)?;let e=f.vector(&energies.0,ENERGY)?;let d=api::dimension(&observable.0)?;let o=<Vec<f64> as Arg>::decode(&observable.0,0)?;if e.is_empty()||e.len()!=o.len(){return Err(api::error("reweight arrays mustmatch nonempty"));}let b0=beta(&f,&beta0.0)?;let b1=beta(&f,&beta1.0)?;Ok(Dynamic(api::quantity(spin::reweight(&e,&o,b0,b1),d)?))});
object_binding!(spin_mean_and_error,"num_native_spin_mean_and_error","crate::sim::spin_models::mean_and_error","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(values:Input)->Dynamic,|_ctx:&CallCtx<'_>,values:Input|{let d=api::dimension(&values.0)?;let x=<Vec<f64> as Arg>::decode(&values.0,0)?;if x.is_empty(){return Err(api::error("mean error needsnonempty sample"));}let(a,b)=spin::mean_and_error(&x);Ok(Dynamic(api::record(vec![("mean",api::quantity(a,d.clone())?),("error",api::quantity(b,d)?)])))});
object_binding!(spin_binder_cumulant,"num_native_spin_binder_cumulant","crate::sim::spin_models::binder_cumulant","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(values:Vec<f64>)->f64,|_ctx:&CallCtx<'_>,values:Vec<f64>|{if values.is_empty()||values.iter().all(|x|*x==0.0){return Err(api::error("Binder cumulant needsnonzero magnetization variance sample"));}Ok(spin::binder_cumulant(&values))});
object_binding!(spin_susceptibility,"num_native_spin_susceptibility","crate::sim::spin_models::susceptibility","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,magnetization:Vec<f64>,inverse_temperature:Input,sites:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,magnetization:Vec<f64>,inverse_temperature:Input,sites:usize|{if magnetization.is_empty()||sites==0{return Err(api::error("susceptibility needsnonempty sample andpositivesites"));}let f=frame(ctx,&f)?;let b=beta(&f,&inverse_temperature.0)?;Ok(Dynamic(f.output(spin::susceptibility(&magnetization,b,sites),[-1,-2,2,0],false)?))});
object_binding!(spin_specific_heat,"num_native_spin_specific_heat","crate::sim::spin_models::specific_heat","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,energies:Input,inverse_temperature:Input,sites:usize)->f64,|ctx:&CallCtx<'_>,f:Handle<Frame>,energies:Input,inverse_temperature:Input,sites:usize|{let f=frame(ctx,&f)?;let e=f.vector(&energies.0,ENERGY)?;if e.is_empty()||sites==0{return Err(api::error("specificheat needsnonempty sample andpositivesites"));}let b=beta(&f,&inverse_temperature.0)?;Ok(spin::specific_heat(&e,b,sites))});
object_binding!(spin_binder_crossing,"num_native_spin_binder_crossing","crate::sim::spin_models::binder_crossing","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(temperatures:Input,small:Vec<f64>,large:Vec<f64>)->Dynamic,|_ctx:&CallCtx<'_>,temperatures:Input,small:Vec<f64>,large:Vec<f64>|{let d=api::dimension(&temperatures.0)?;let t=<Vec<f64> as Arg>::decode(&temperatures.0,0)?;if t.len()<2||t.len()!=small.len()||t.len()!=large.len(){return Err(api::error("Binder crossing arrays mustmatch andcontain>=2points"));}let out=spin::binder_crossing(&t,&small,&large).map(|v|api::quantity(v,d)).transpose()?;Ok(Dynamic(Value::Variant(crate::kernel::value::VariantValue{tag:if out.is_some(){"Some"}else{"None"}.into(),payload:out.into_iter().collect()})))});
object_binding!(spin_fit_exponent,"num_native_spin_fit_exponent","crate::sim::spin_models::fit_exponent","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(sizes:Input,values:Input,size_reference:Input,value_reference:Input)->f64,|_ctx:&CallCtx<'_>,sizes:Input,values:Input,size_reference:Input,value_reference:Input|{let(sr,sd)=api::reference(&size_reference.0)?;let(vr,vd)=api::reference(&value_reference.0)?;api::same_dimension(&sd,&api::dimension(&sizes.0)?)?;api::same_dimension(&vd,&api::dimension(&values.0)?)?;let x=<Vec<f64> as Arg>::decode(&sizes.0,0)?.into_iter().map(|x|x/sr).collect::<Vec<_>>();let y=<Vec<f64> as Arg>::decode(&values.0,0)?.into_iter().map(|x|x/vr).collect::<Vec<_>>();if x.len()<2||x.len()!=y.len()||x.iter().chain(&y).any(|x|*x<=0.0)||x.iter().all(|a|(*a-x[0]).abs()<=0.0){return Err(api::error("exponent fit needspositive matched values anddistinctsizes"));}Ok(spin::fit_exponent(&x,&y))});
object_binding!(spin_collapse_variable,"num_native_spin_collapse_variable","crate::sim::spin_models::collapse_variable","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(temperature:Input,critical_temperature:Input,nu:f64,size:Input,temperature_reference:Input,size_reference:Input)->f64,|_ctx:&CallCtx<'_>,temperature:Input,critical_temperature:Input,nu:f64,size:Input,temperature_reference:Input,size_reference:Input|{let(tr,td)=api::reference(&temperature_reference.0)?;let(lr,ld)=api::reference(&size_reference.0)?;api::same_dimension(&td,&api::dimension(&temperature.0)?)?;api::same_dimension(&td,&api::dimension(&critical_temperature.0)?)?;api::same_dimension(&ld,&api::dimension(&size.0)?)?;let t=crate::rules::numerical_abi::numeric(&temperature.0,0)?/tr;let tc=crate::rules::numerical_abi::numeric(&critical_temperature.0,0)?/tr;let l=crate::rules::numerical_abi::numeric(&size.0,0)?/lr;if nu<=0.0||l<=0.0{return Err(api::error("scaling exponent/reference sizes must bepositive"));}Ok(spin::collapse_variable(t,tc,nu,l))});
object_binding!(spin_onsager_tc,"num_native_spin_onsager_tc","crate::sim::spin_models::onsager_tc","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,coupling:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,coupling:Input|{let f=frame(ctx,&f)?;let j=f.scalar(&coupling.0,ENERGY,false)?;if j<=0.0{return Err(api::error("exact ferromagnetic benchmark requirespositivecoupling"));}Ok(Dynamic(f.output(spin::onsager_tc(j),ENERGY,true)?))});
object_binding!(spin_onsager_energy_per_site,"num_native_spin_onsager_energy_per_site","crate::sim::spin_models::onsager_energy_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,coupling:Input,temperature:Input)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,coupling:Input,temperature:Input|{let f=frame(ctx,&f)?;let j=f.scalar(&coupling.0,ENERGY,false)?;if j<=0.0{return Err(api::error("exact ferromagnetic benchmark requirespositivecoupling"));}let t=f.scalar(&temperature.0,ENERGY,true)?;if t<=0.0{return Err(api::error("benchmark temperature mustbepositive"));}Ok(Dynamic(f.output(spin::onsager_energy_per_site(t,j),ENERGY,false)?))});
object_binding!(spin_onsager_magnetisation,"num_native_spin_onsager_magnetisation","crate::sim::spin_models::onsager_magnetisation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,coupling:Input,temperature:Input)->f64,|ctx:&CallCtx<'_>,f:Handle<Frame>,coupling:Input,temperature:Input|{let f=frame(ctx,&f)?;let j=f.scalar(&coupling.0,ENERGY,false)?;if j<=0.0{return Err(api::error("exact ferromagnetic benchmark requirespositivecoupling"));}let t=f.scalar(&temperature.0,ENERGY,true)?;if t<=0.0{return Err(api::error("benchmark temperature mustbepositive"));}Ok(spin::onsager_magnetisation(t,j))});
object_binding!(spin_ising1d_energy_per_site,"num_native_spin_ising1d_energy_per_site","crate::sim::spin_models::ising1d_energy_per_site","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,coupling:Input,inverse_temperature:Input,n:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,coupling:Input,inverse_temperature:Input,n:usize|{let f=frame(ctx,&f)?;let j=f.scalar(&coupling.0,ENERGY,false)?;if j<=0.0{return Err(api::error("exact ferromagnetic benchmark requirespositivecoupling"));}let b=beta(&f,&inverse_temperature.0)?;if n==0||n>4096{return Err(api::error("Ising chain size mustbe1..4096"));}Ok(Dynamic(f.output(spin::ising1d_energy_per_site(b,j,n),ENERGY,false)?))});
struct Neighbours(Vec<Option<usize>>);
impl api::NativeOutput for Neighbours {
    fn ty() -> TypeExpr {
        api::tensor_type(1, <Option<usize> as Output>::ty())
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.0.len()],
            data: self
                .0
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: <Option<usize> as Output>::ty(),
        }))
    }
}
struct Sse {
    native: Option<sse::SseXxz>,
    frame: Frame,
}
impl NativeType for Sse {
    const KIND: &'static str = "rssn/native/sim/sse-xxz/v1";
    const MAX_BYTES: usize = 8_388_608;

    fn validate(&self) -> Result<(), OpError> {
        if self.native.is_none() {
            return Err(api::error("SSE state invalidated after budget exhaustion"));
        }
        self.frame.validate()
    }
}
fn sse_sample(v: sse::SseSample) -> Value {
    api::record(vec![
        ("operators", Value::int(v.n_ops)),
        ("magnetization", Value::Float(v.magnetisation)),
        ("staggered", Value::Float(v.staggered)),
    ])
}
fn sse_result(
    v: sse::SseResult,
    f: &Frame,
) -> Result<Value, OpError> {
    Ok(api::record(vec![
        ("energy", f.output(v.energy, ENERGY, false)?),
        ("energy_error", f.output(v.energy_error, ENERGY, false)?),
        ("heat_capacity_over_kB", Output::encode(v.specific_heat)?),
        ("magnetization_squared", Output::encode(v.magnetisation_sq)?),
        ("staggered_squared", Output::encode(v.staggered_sq)?),
        ("mean_operators", Output::encode(v.mean_ops)?),
    ]))
}
fn ising_parameters(
    f: &Frame,
    width: usize,
    height: usize,
    temperature: &Input,
    steps: usize,
) -> Result<ising::IsingParameters, OpError> {
    if width == 0 || height == 0 || width.checked_mul(height).is_none_or(|n| n > 4096) {
        return Err(api::error("Ising grid allocation bound4096"));
    }
    work(width * height, 1, steps)?;
    let temperature = f.scalar(&temperature.0, ENERGY, true)?;
    if temperature < 0.0 {
        return Err(api::error("Ising temperature must benonnegative"));
    }
    Ok(ising::IsingParameters {
        width,
        height,
        temperature,
        mc_steps: steps,
    })
}
fn ising_observables(
    v: &ising::IsingObservables,
    f: &Frame,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("energy", f.output(v.energy, ENERGY, false)?),
        ("magnetization", Output::encode(v.magnetization)?),
        ("heat_capacity_over_kB", Output::encode(v.specific_heat)?),
        (
            "susceptibility",
            f.output(v.susceptibility, [-1, -2, 2, 0], false)?,
        ),
        ("binder", Output::encode(v.binder)?),
        ("mean_cluster_fraction", Output::encode(v.mean_cluster)?),
    ])))
}
fn lattice_observables(
    v: &ising::LatticeObservables,
    f: &Frame,
) -> Result<Dynamic, OpError> {
    Ok(Dynamic(api::record(vec![
        ("energy", f.output(v.energy, ENERGY, false)?),
        ("magnetization", Output::encode(v.magnetization)?),
        ("heat_capacity_over_kB", Output::encode(v.specific_heat)?),
        (
            "susceptibility",
            f.output(v.susceptibility, [-1, -2, 2, 0], false)?,
        ),
        ("binder", Output::encode(v.binder)?),
    ])))
}
object_binding!(sse_chain_bonds,"num_native_sse_chain_bonds","crate::sim::sse::chain_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize,periodic:bool)->Vec<[usize;2]>,|_ctx:&CallCtx<'_>,n:usize,periodic:bool|{if !(2..=4096).contains(&n){return Err(api::error("SSE chain size mustbe2..4096"));}Ok(sse::chain_bonds(n,periodic).into_iter().map(|(a,b)|[a,b]).collect())});
object_binding!(sse_square_bonds,"num_native_sse_square_bonds","crate::sim::sse::square_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(width:usize,height:usize,periodic:bool)->Vec<[usize;2]>,|_ctx:&CallCtx<'_>,width:usize,height:usize,periodic:bool|{if width<2||height<2||width.checked_mul(height).is_none_or(|n|n>4096){return Err(api::error("SSE square allocation bound4096"));}Ok(sse::square_bonds(width,height,periodic).into_iter().map(|(a,b)|[a,b]).collect())});
object_binding!(sse_new,"num_native_sse_new","crate::sim::sse::SseXxz::new","Actual native XXZ with transverse Jxy=one exchange-energy reference, unit spin hbar=1, physical inverse energy temperature and bounded operator storage. Exhausting an explicit work budget invalidates the mutable state rather than exposing a partial directed-loop update.",false,(f:Handle<Frame>,n_sites:usize,bonds:Vec<[usize;2]>,anisotropy:Input,inverse_temperature:Input,seed:u64)->NewObject<Sse>,|ctx:&CallCtx<'_>,f:Handle<Frame>,n_sites:usize,bonds:Vec<[usize;2]>,anisotropy:Input,inverse_temperature:Input,seed:u64|{if !(2..=512).contains(&n_sites)||bonds.is_empty()||bonds.len()>4096{return Err(api::error("SSE sites/bonds exceed boundednonempty contract"));}let f=frame(ctx,&f)?;let jz=f.scalar(&anisotropy.0,ENERGY,false)?;let beta=beta(&f,&inverse_temperature.0)?;if jz.abs()>4.0||beta<=0.0{return Err(api::error("SSE needspositivebeta and |Jz/Jxy|<=4"));}let native=sse::SseXxz::new(n_sites,bonds.into_iter().map(|[a,b]|(a,b)).collect(),jz,beta,seed).map_err(api::error)?;let dims=f.dims.to_vec();Ok(NewObject::new(Sse{native:Some(native),frame:f},dims))});
object_binding!(sse_n_bonds,"num_native_sse_n_bonds","crate::sim::sse::SseXxz::n_bonds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Sse>)->usize,|ctx:&CallCtx<'_>,simulation:Handle<Sse>|{simulation.with(ctx,|s|Ok(s.native.as_ref().ok_or_else(||api::error("SSE state invalidated"))?.n_bonds()))});
object_binding!(sse_constant,"num_native_sse_constant","crate::sim::sse::SseXxz::constant","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Sse>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Sse>|{simulation.with(ctx,|s|Ok(Dynamic(s.frame.output(s.native.as_ref().ok_or_else(||api::error("SSE state invalidated"))?.constant(),ENERGY,false)?)))});
object_binding!(sse_sample,"num_native_sse_sample","crate::sim::sse::SseXxz::sample","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Sse>)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Sse>|{simulation.with(ctx,|s|Ok(Dynamic(sse_sample(s.native.as_ref().ok_or_else(||api::error("SSE state invalidated"))?.sample()))))});
object_binding!(sse_sweep,"num_native_sse_sweep","crate::sim::sse::SseXxz::sweep","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Sse>,operation_budget:usize)->bool,|ctx:&CallCtx<'_>,simulation:Handle<Sse>,operation_budget:usize|{if operation_budget==0||operation_budget>50_000_000{return Err(api::error("SSE operation budget mustbe1..50million"));}simulation.with_mut(ctx,|s|{let result=sse::with_native_budget(operation_budget,65536,||s.native.as_mut().expect("validated SSE handle").sweep());if let Err(e)=result{s.native=None;return Err(api::error(e));}Ok(true)})});
object_binding!(sse_run,"num_native_sse_run","crate::sim::sse::SseXxz::run","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(simulation:Handle<Sse>,thermalization:usize,sweeps:usize,operation_budget:usize)->Dynamic,|ctx:&CallCtx<'_>,simulation:Handle<Sse>,thermalization:usize,sweeps:usize,operation_budget:usize|{if thermalization>100000||sweeps==0||sweeps>100000||operation_budget==0||operation_budget>50_000_000{return Err(api::error("SSE sweep/work budget invalid"));}simulation.with_mut(ctx,|s|{let result=sse::with_native_budget(operation_budget,65536,||s.native.as_mut().expect("validated SSE handle").run(thermalization,sweeps));match result{Ok(v)=>Ok(Dynamic(sse_result(v,&s.frame)?)),Err(e)=>{s.native=None;Err(api::error(e))}}})});
object_binding!(ising_run_ising_simulation,"num_native_ising_run_ising_simulation","crate::sim::models_ising_statistical::run_ising_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,width:usize,height:usize,temperature:Input,steps:usize)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,width:usize,height:usize,temperature:Input,steps:usize|{let f=frame(ctx,&f)?;let params=ising_parameters(&f,width,height,&temperature,steps)?;let(spins,magnetization)=ising::run_ising_simulation(&params);Ok(Dynamic(api::record(vec![("spins",Output::encode(spins.into_iter().map(i64::from).collect::<Vec<_>>())?),("magnetization",Output::encode(magnetization)?)])))});
object_binding!(ising_run_ising_simulation_seeded,"num_native_ising_run_ising_simulation_seeded","crate::sim::models_ising_statistical::run_ising_simulation_seeded","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,width:usize,height:usize,temperature:Input,steps:usize,seed:u64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,width:usize,height:usize,temperature:Input,steps:usize,seed:u64|{let f=frame(ctx,&f)?;let params=ising_parameters(&f,width,height,&temperature,steps)?;let(spins,magnetization)=ising::run_ising_simulation_seeded(&params,seed);Ok(Dynamic(api::record(vec![("spins",Output::encode(spins.into_iter().map(i64::from).collect::<Vec<_>>())?),("magnetization",Output::encode(magnetization)?)])))});
object_binding!(ising_run_wolff_simulation,"num_native_ising_run_wolff_simulation","crate::sim::models_ising_statistical::run_wolff_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,width:usize,height:usize,temperature:Input,thermalization:usize,sweeps:usize,seed:u64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,width:usize,height:usize,temperature:Input,thermalization:usize,sweeps:usize,seed:u64|{let f=frame(ctx,&f)?;let total=thermalization.checked_add(sweeps).ok_or_else(||api::error("sweep count overflow"))?;let p=ising_parameters(&f,width,height,&temperature,total)?;if width<3||height<3||p.temperature<=0.0||sweeps<2{return Err(api::error("Wolff needsperiodicsides>=3,positiveT,>=2samples"));}ising_observables(&(ising::run_wolff_simulation(width,height,p.temperature,thermalization,sweeps,seed)),&f)});
object_binding!(ising_run_lattice_simulation_d1,"num_native_ising_run_lattice_simulation_d1","crate::sim::models_ising_statistical::run_lattice_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;1],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;1],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64|{let n=sites(dims,periodic)?;let f=frame(ctx,&f)?;let t=f.scalar(&temperature.0,ENERGY,true)?;if t<=0.0||sweeps<2{return Err(api::error("lattice simulation needspositiveT and>=2samples"));}work(n,1,thermalization.checked_add(sweeps).ok_or_else(||api::error("sweep count overflow"))?)?;let model=match model.as_str(){"Ising"=>ising::SpinModel::Ising,"Xy"=>ising::SpinModel::Xy,"Heisenberg"=>ising::SpinModel::Heisenberg,"Potts" if (2..=256).contains(&q)=>ising::SpinModel::Potts(q),_=>return Err(api::error("model mustbeIsing/Xy/Heisenberg/Potts(q2..256)"))};let algorithm=match update(&algorithm,proposal)?{spin::Update::Metropolis(x)=>ising::SpinAlgorithm::Metropolis(x),spin::Update::HeatBath=>ising::SpinAlgorithm::HeatBath,spin::Update::Wolff=>ising::SpinAlgorithm::Wolff,spin::Update::SwendsenWang=>ising::SpinAlgorithm::SwendsenWang};lattice_observables(&(ising::run_lattice_simulation(dims,periodic,model,algorithm,t,(thermalization,sweeps),seed).map_err(api::error)?),&f)});
object_binding!(ising_run_lattice_simulation_d2,"num_native_ising_run_lattice_simulation_d2","crate::sim::models_ising_statistical::run_lattice_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;2],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;2],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64|{let n=sites(dims,periodic)?;let f=frame(ctx,&f)?;let t=f.scalar(&temperature.0,ENERGY,true)?;if t<=0.0||sweeps<2{return Err(api::error("lattice simulation needspositiveT and>=2samples"));}work(n,1,thermalization.checked_add(sweeps).ok_or_else(||api::error("sweep count overflow"))?)?;let model=match model.as_str(){"Ising"=>ising::SpinModel::Ising,"Xy"=>ising::SpinModel::Xy,"Heisenberg"=>ising::SpinModel::Heisenberg,"Potts" if (2..=256).contains(&q)=>ising::SpinModel::Potts(q),_=>return Err(api::error("model mustbeIsing/Xy/Heisenberg/Potts(q2..256)"))};let algorithm=match update(&algorithm,proposal)?{spin::Update::Metropolis(x)=>ising::SpinAlgorithm::Metropolis(x),spin::Update::HeatBath=>ising::SpinAlgorithm::HeatBath,spin::Update::Wolff=>ising::SpinAlgorithm::Wolff,spin::Update::SwendsenWang=>ising::SpinAlgorithm::SwendsenWang};lattice_observables(&(ising::run_lattice_simulation(dims,periodic,model,algorithm,t,(thermalization,sweeps),seed).map_err(api::error)?),&f)});
object_binding!(ising_run_lattice_simulation_d3,"num_native_ising_run_lattice_simulation_d3","crate::sim::models_ising_statistical::run_lattice_simulation","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(f:Handle<Frame>,dims:[usize;3],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64)->Dynamic,|ctx:&CallCtx<'_>,f:Handle<Frame>,dims:[usize;3],periodic:bool,model:String,q:usize,algorithm:String,proposal:f64,temperature:Input,thermalization:usize,sweeps:usize,seed:u64|{let n=sites(dims,periodic)?;let f=frame(ctx,&f)?;let t=f.scalar(&temperature.0,ENERGY,true)?;if t<=0.0||sweeps<2{return Err(api::error("lattice simulation needspositiveT and>=2samples"));}work(n,1,thermalization.checked_add(sweeps).ok_or_else(||api::error("sweep count overflow"))?)?;let model=match model.as_str(){"Ising"=>ising::SpinModel::Ising,"Xy"=>ising::SpinModel::Xy,"Heisenberg"=>ising::SpinModel::Heisenberg,"Potts" if (2..=256).contains(&q)=>ising::SpinModel::Potts(q),_=>return Err(api::error("model mustbeIsing/Xy/Heisenberg/Potts(q2..256)"))};let algorithm=match update(&algorithm,proposal)?{spin::Update::Metropolis(x)=>ising::SpinAlgorithm::Metropolis(x),spin::Update::HeatBath=>ising::SpinAlgorithm::HeatBath,spin::Update::Wolff=>ising::SpinAlgorithm::Wolff,spin::Update::SwendsenWang=>ising::SpinAlgorithm::SwendsenWang};lattice_observables(&(ising::run_lattice_simulation(dims,periodic,model,algorithm,t,(thermalization,sweeps),seed).map_err(api::error)?),&f)});
object_binding!(ising_phase_transition_scenario,"num_native_ising_phase_transition_scenario","crate::sim::models_ising_statistical::simulate_ising_phase_transition_scenario","IO: execute the original native hard-coded phase-transition scenario, create its caller-selected output directory and write CSV/NPY files. Explicit work acknowledgement requires205million site updates. Coordinates/temperature are the original reduced-unit demonstration; no physical values are silently erased.",true,(output_directory:String,site_update_budget:usize)->bool,|_ctx:&CallCtx<'_>,output_directory:String,site_update_budget:usize|{const REQUIRED:usize=41*50*50*2000;if output_directory.is_empty()||output_directory.len()>4096||site_update_budget<REQUIRED{return Err(api::error("native phase scenario writesCSV/NPY and requires explicit>=205million siteupdate budget"));}ising::simulate_ising_phase_transition_scenario(std::path::Path::new(&output_directory)).map_err(api::error)?;Ok(true)});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &sse_chain_bonds::BINDING,
            &sse_square_bonds::BINDING,
            &sse_new::BINDING,
            &sse_n_bonds::BINDING,
            &sse_constant::BINDING,
            &sse_sample::BINDING,
            &sse_sweep::BINDING,
            &sse_run::BINDING,
            &ising_run_ising_simulation::BINDING,
            &ising_run_ising_simulation_seeded::BINDING,
            &ising_run_wolff_simulation::BINDING,
            &ising_run_lattice_simulation_d1::BINDING,
            &ising_run_lattice_simulation_d2::BINDING,
            &ising_run_lattice_simulation_d3::BINDING,
            &ising_phase_transition_scenario::BINDING,
            &spin_lattice_new_d1::BINDING,
            &spin_lattice_sites_d1::BINDING,
            &spin_lattice_bonds_d1::BINDING,
            &spin_lattice_diagonal_bonds_d1::BINDING,
            &spin_lattice_neighbours_d1::BINDING,
            &spin_lattice_second_neighbours_d1::BINDING,
            &spin_potts_new_d1::BINDING,
            &spin_potts_energy_d1::BINDING,
            &spin_potts_order_parameter_d1::BINDING,
            &spin_potts_metropolis_sweep_d1::BINDING,
            &spin_potts_heat_bath_sweep_d1::BINDING,
            &spin_potts_wolff_step_d1::BINDING,
            &spin_potts_swendsen_wang_step_d1::BINDING,
            &spin_potts_update_d1::BINDING,
            &spin_potts_tempering_new_d1::BINDING,
            &spin_on_new_d1_n1::BINDING,
            &spin_on_cold_start_d1_n1::BINDING,
            &spin_on_energy_d1_n1::BINDING,
            &spin_on_magnetisation_d1_n1::BINDING,
            &spin_on_abs_magnetisation_per_site_d1_n1::BINDING,
            &spin_on_heat_bath_sweep_d1_n1::BINDING,
            &spin_on_overrelax_sweep_d1_n1::BINDING,
            &spin_on_wolff_step_d1_n1::BINDING,
            &spin_on_swendsen_wang_step_d1_n1::BINDING,
            &spin_on_local_field_d1_n1::BINDING,
            &spin_on_metropolis_sweep_d1_n1::BINDING,
            &spin_on_update_d1_n1::BINDING,
            &spin_on_couplings_d1_n1::BINDING,
            &spin_tempering_new_d1_n1::BINDING,
            &spin_tempering_run_d1_n1::BINDING,
            &spin_wang_landau_on_d1_n1::BINDING,
            &spin_replica_exchange_wang_landau_d1_n1::BINDING,
            &spin_on_new_d1_n2::BINDING,
            &spin_on_cold_start_d1_n2::BINDING,
            &spin_on_energy_d1_n2::BINDING,
            &spin_on_magnetisation_d1_n2::BINDING,
            &spin_on_abs_magnetisation_per_site_d1_n2::BINDING,
            &spin_on_heat_bath_sweep_d1_n2::BINDING,
            &spin_on_overrelax_sweep_d1_n2::BINDING,
            &spin_on_wolff_step_d1_n2::BINDING,
            &spin_on_swendsen_wang_step_d1_n2::BINDING,
            &spin_on_local_field_d1_n2::BINDING,
            &spin_on_metropolis_sweep_d1_n2::BINDING,
            &spin_on_update_d1_n2::BINDING,
            &spin_on_couplings_d1_n2::BINDING,
            &spin_tempering_new_d1_n2::BINDING,
            &spin_tempering_run_d1_n2::BINDING,
            &spin_wang_landau_on_d1_n2::BINDING,
            &spin_replica_exchange_wang_landau_d1_n2::BINDING,
            &spin_on_new_d1_n3::BINDING,
            &spin_on_cold_start_d1_n3::BINDING,
            &spin_on_energy_d1_n3::BINDING,
            &spin_on_magnetisation_d1_n3::BINDING,
            &spin_on_abs_magnetisation_per_site_d1_n3::BINDING,
            &spin_on_heat_bath_sweep_d1_n3::BINDING,
            &spin_on_overrelax_sweep_d1_n3::BINDING,
            &spin_on_wolff_step_d1_n3::BINDING,
            &spin_on_swendsen_wang_step_d1_n3::BINDING,
            &spin_on_local_field_d1_n3::BINDING,
            &spin_on_metropolis_sweep_d1_n3::BINDING,
            &spin_on_update_d1_n3::BINDING,
            &spin_on_couplings_d1_n3::BINDING,
            &spin_tempering_new_d1_n3::BINDING,
            &spin_tempering_run_d1_n3::BINDING,
            &spin_wang_landau_on_d1_n3::BINDING,
            &spin_replica_exchange_wang_landau_d1_n3::BINDING,
            &spin_on_new_d1_n4::BINDING,
            &spin_on_cold_start_d1_n4::BINDING,
            &spin_on_energy_d1_n4::BINDING,
            &spin_on_magnetisation_d1_n4::BINDING,
            &spin_on_abs_magnetisation_per_site_d1_n4::BINDING,
            &spin_on_heat_bath_sweep_d1_n4::BINDING,
            &spin_on_overrelax_sweep_d1_n4::BINDING,
            &spin_on_wolff_step_d1_n4::BINDING,
            &spin_on_swendsen_wang_step_d1_n4::BINDING,
            &spin_on_local_field_d1_n4::BINDING,
            &spin_on_metropolis_sweep_d1_n4::BINDING,
            &spin_on_update_d1_n4::BINDING,
            &spin_on_couplings_d1_n4::BINDING,
            &spin_tempering_new_d1_n4::BINDING,
            &spin_tempering_run_d1_n4::BINDING,
            &spin_wang_landau_on_d1_n4::BINDING,
            &spin_replica_exchange_wang_landau_d1_n4::BINDING,
            &spin_potts_tempering_run_d1::BINDING,
            &spin_enumerate_ising_dos_d1::BINDING,
            &spin_wang_landau_ising_d1::BINDING,
            &spin_lattice_new_d2::BINDING,
            &spin_lattice_sites_d2::BINDING,
            &spin_lattice_bonds_d2::BINDING,
            &spin_lattice_diagonal_bonds_d2::BINDING,
            &spin_lattice_neighbours_d2::BINDING,
            &spin_lattice_second_neighbours_d2::BINDING,
            &spin_potts_new_d2::BINDING,
            &spin_potts_energy_d2::BINDING,
            &spin_potts_order_parameter_d2::BINDING,
            &spin_potts_metropolis_sweep_d2::BINDING,
            &spin_potts_heat_bath_sweep_d2::BINDING,
            &spin_potts_wolff_step_d2::BINDING,
            &spin_potts_swendsen_wang_step_d2::BINDING,
            &spin_potts_update_d2::BINDING,
            &spin_potts_tempering_new_d2::BINDING,
            &spin_on_new_d2_n1::BINDING,
            &spin_on_cold_start_d2_n1::BINDING,
            &spin_on_energy_d2_n1::BINDING,
            &spin_on_magnetisation_d2_n1::BINDING,
            &spin_on_abs_magnetisation_per_site_d2_n1::BINDING,
            &spin_on_heat_bath_sweep_d2_n1::BINDING,
            &spin_on_overrelax_sweep_d2_n1::BINDING,
            &spin_on_wolff_step_d2_n1::BINDING,
            &spin_on_swendsen_wang_step_d2_n1::BINDING,
            &spin_on_local_field_d2_n1::BINDING,
            &spin_on_metropolis_sweep_d2_n1::BINDING,
            &spin_on_update_d2_n1::BINDING,
            &spin_on_couplings_d2_n1::BINDING,
            &spin_tempering_new_d2_n1::BINDING,
            &spin_tempering_run_d2_n1::BINDING,
            &spin_wang_landau_on_d2_n1::BINDING,
            &spin_replica_exchange_wang_landau_d2_n1::BINDING,
            &spin_on_new_d2_n2::BINDING,
            &spin_on_cold_start_d2_n2::BINDING,
            &spin_on_energy_d2_n2::BINDING,
            &spin_on_magnetisation_d2_n2::BINDING,
            &spin_on_abs_magnetisation_per_site_d2_n2::BINDING,
            &spin_on_heat_bath_sweep_d2_n2::BINDING,
            &spin_on_overrelax_sweep_d2_n2::BINDING,
            &spin_on_wolff_step_d2_n2::BINDING,
            &spin_on_swendsen_wang_step_d2_n2::BINDING,
            &spin_on_local_field_d2_n2::BINDING,
            &spin_on_metropolis_sweep_d2_n2::BINDING,
            &spin_on_update_d2_n2::BINDING,
            &spin_on_couplings_d2_n2::BINDING,
            &spin_tempering_new_d2_n2::BINDING,
            &spin_tempering_run_d2_n2::BINDING,
            &spin_wang_landau_on_d2_n2::BINDING,
            &spin_replica_exchange_wang_landau_d2_n2::BINDING,
            &spin_on_new_d2_n3::BINDING,
            &spin_on_cold_start_d2_n3::BINDING,
            &spin_on_energy_d2_n3::BINDING,
            &spin_on_magnetisation_d2_n3::BINDING,
            &spin_on_abs_magnetisation_per_site_d2_n3::BINDING,
            &spin_on_heat_bath_sweep_d2_n3::BINDING,
            &spin_on_overrelax_sweep_d2_n3::BINDING,
            &spin_on_wolff_step_d2_n3::BINDING,
            &spin_on_swendsen_wang_step_d2_n3::BINDING,
            &spin_on_local_field_d2_n3::BINDING,
            &spin_on_metropolis_sweep_d2_n3::BINDING,
            &spin_on_update_d2_n3::BINDING,
            &spin_on_couplings_d2_n3::BINDING,
            &spin_tempering_new_d2_n3::BINDING,
            &spin_tempering_run_d2_n3::BINDING,
            &spin_wang_landau_on_d2_n3::BINDING,
            &spin_replica_exchange_wang_landau_d2_n3::BINDING,
            &spin_on_new_d2_n4::BINDING,
            &spin_on_cold_start_d2_n4::BINDING,
            &spin_on_energy_d2_n4::BINDING,
            &spin_on_magnetisation_d2_n4::BINDING,
            &spin_on_abs_magnetisation_per_site_d2_n4::BINDING,
            &spin_on_heat_bath_sweep_d2_n4::BINDING,
            &spin_on_overrelax_sweep_d2_n4::BINDING,
            &spin_on_wolff_step_d2_n4::BINDING,
            &spin_on_swendsen_wang_step_d2_n4::BINDING,
            &spin_on_local_field_d2_n4::BINDING,
            &spin_on_metropolis_sweep_d2_n4::BINDING,
            &spin_on_update_d2_n4::BINDING,
            &spin_on_couplings_d2_n4::BINDING,
            &spin_tempering_new_d2_n4::BINDING,
            &spin_tempering_run_d2_n4::BINDING,
            &spin_wang_landau_on_d2_n4::BINDING,
            &spin_replica_exchange_wang_landau_d2_n4::BINDING,
            &spin_potts_tempering_run_d2::BINDING,
            &spin_enumerate_ising_dos_d2::BINDING,
            &spin_wang_landau_ising_d2::BINDING,
            &spin_lattice_new_d3::BINDING,
            &spin_lattice_sites_d3::BINDING,
            &spin_lattice_bonds_d3::BINDING,
            &spin_lattice_diagonal_bonds_d3::BINDING,
            &spin_lattice_neighbours_d3::BINDING,
            &spin_lattice_second_neighbours_d3::BINDING,
            &spin_potts_new_d3::BINDING,
            &spin_potts_energy_d3::BINDING,
            &spin_potts_order_parameter_d3::BINDING,
            &spin_potts_metropolis_sweep_d3::BINDING,
            &spin_potts_heat_bath_sweep_d3::BINDING,
            &spin_potts_wolff_step_d3::BINDING,
            &spin_potts_swendsen_wang_step_d3::BINDING,
            &spin_potts_update_d3::BINDING,
            &spin_potts_tempering_new_d3::BINDING,
            &spin_on_new_d3_n1::BINDING,
            &spin_on_cold_start_d3_n1::BINDING,
            &spin_on_energy_d3_n1::BINDING,
            &spin_on_magnetisation_d3_n1::BINDING,
            &spin_on_abs_magnetisation_per_site_d3_n1::BINDING,
            &spin_on_heat_bath_sweep_d3_n1::BINDING,
            &spin_on_overrelax_sweep_d3_n1::BINDING,
            &spin_on_wolff_step_d3_n1::BINDING,
            &spin_on_swendsen_wang_step_d3_n1::BINDING,
            &spin_on_local_field_d3_n1::BINDING,
            &spin_on_metropolis_sweep_d3_n1::BINDING,
            &spin_on_update_d3_n1::BINDING,
            &spin_on_couplings_d3_n1::BINDING,
            &spin_tempering_new_d3_n1::BINDING,
            &spin_tempering_run_d3_n1::BINDING,
            &spin_wang_landau_on_d3_n1::BINDING,
            &spin_replica_exchange_wang_landau_d3_n1::BINDING,
            &spin_on_new_d3_n2::BINDING,
            &spin_on_cold_start_d3_n2::BINDING,
            &spin_on_energy_d3_n2::BINDING,
            &spin_on_magnetisation_d3_n2::BINDING,
            &spin_on_abs_magnetisation_per_site_d3_n2::BINDING,
            &spin_on_heat_bath_sweep_d3_n2::BINDING,
            &spin_on_overrelax_sweep_d3_n2::BINDING,
            &spin_on_wolff_step_d3_n2::BINDING,
            &spin_on_swendsen_wang_step_d3_n2::BINDING,
            &spin_on_local_field_d3_n2::BINDING,
            &spin_on_metropolis_sweep_d3_n2::BINDING,
            &spin_on_update_d3_n2::BINDING,
            &spin_on_couplings_d3_n2::BINDING,
            &spin_tempering_new_d3_n2::BINDING,
            &spin_tempering_run_d3_n2::BINDING,
            &spin_wang_landau_on_d3_n2::BINDING,
            &spin_replica_exchange_wang_landau_d3_n2::BINDING,
            &spin_on_new_d3_n3::BINDING,
            &spin_on_cold_start_d3_n3::BINDING,
            &spin_on_energy_d3_n3::BINDING,
            &spin_on_magnetisation_d3_n3::BINDING,
            &spin_on_abs_magnetisation_per_site_d3_n3::BINDING,
            &spin_on_heat_bath_sweep_d3_n3::BINDING,
            &spin_on_overrelax_sweep_d3_n3::BINDING,
            &spin_on_wolff_step_d3_n3::BINDING,
            &spin_on_swendsen_wang_step_d3_n3::BINDING,
            &spin_on_local_field_d3_n3::BINDING,
            &spin_on_metropolis_sweep_d3_n3::BINDING,
            &spin_on_update_d3_n3::BINDING,
            &spin_on_couplings_d3_n3::BINDING,
            &spin_tempering_new_d3_n3::BINDING,
            &spin_tempering_run_d3_n3::BINDING,
            &spin_wang_landau_on_d3_n3::BINDING,
            &spin_replica_exchange_wang_landau_d3_n3::BINDING,
            &spin_on_new_d3_n4::BINDING,
            &spin_on_cold_start_d3_n4::BINDING,
            &spin_on_energy_d3_n4::BINDING,
            &spin_on_magnetisation_d3_n4::BINDING,
            &spin_on_abs_magnetisation_per_site_d3_n4::BINDING,
            &spin_on_heat_bath_sweep_d3_n4::BINDING,
            &spin_on_overrelax_sweep_d3_n4::BINDING,
            &spin_on_wolff_step_d3_n4::BINDING,
            &spin_on_swendsen_wang_step_d3_n4::BINDING,
            &spin_on_local_field_d3_n4::BINDING,
            &spin_on_metropolis_sweep_d3_n4::BINDING,
            &spin_on_update_d3_n4::BINDING,
            &spin_on_couplings_d3_n4::BINDING,
            &spin_tempering_new_d3_n4::BINDING,
            &spin_tempering_run_d3_n4::BINDING,
            &spin_wang_landau_on_d3_n4::BINDING,
            &spin_replica_exchange_wang_landau_d3_n4::BINDING,
            &spin_potts_tempering_run_d3::BINDING,
            &spin_enumerate_ising_dos_d3::BINDING,
            &spin_wang_landau_ising_d3::BINDING,
            &spin_heat_bath_direction_n1::BINDING,
            &spin_heat_bath_direction_n2::BINDING,
            &spin_heat_bath_direction_n3::BINDING,
            &spin_heat_bath_direction_n4::BINDING,
            &spin_canonical_from_dos::BINDING,
            &spin_reweight::BINDING,
            &spin_mean_and_error::BINDING,
            &spin_binder_cumulant::BINDING,
            &spin_susceptibility::BINDING,
            &spin_specific_heat::BINDING,
            &spin_binder_crossing::BINDING,
            &spin_fit_exponent::BINDING,
            &spin_collapse_variable::BINDING,
            &spin_onsager_tc::BINDING,
            &spin_onsager_energy_per_site::BINDING,
            &spin_onsager_magnetisation::BINDING,
            &spin_ising1d_energy_per_site::BINDING,
        ],
        "numerical_vector",
    )
}
