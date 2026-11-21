//! Mixed-dimension one-dimensional compressible gas calls.
use super::CallCtx;
use super::DimExpr;
use super::Dynamic;
use super::Handle;
use super::Input;
use super::Installer;
use super::NativeArg;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::Output;
use super::RuleError;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::api;
use super::err;
use super::mechanics::REFERENCE_DOMAINS;
use super::mechanics::References;
use super::object_binding;
use super::structured_ops::checked_snapshot;
use crate::rules::numerical_abi::numeric;
use crate::sim::gas_dynamics as gas;

const GAS_DOMAINS: [[i64; 4]; 3] = [[-3, 0, 1, 0], [1, -1, 0, 0], [-1, -2, 1, 0]];
const GAS_COMPONENTS: [&str; 3] = ["density", "velocity", "pressure"];

fn gas_domain(component: usize) -> DimExpr {
    DimExpr::Mul(
        REFERENCE_DOMAINS
            .into_iter()
            .zip(GAS_DOMAINS[component])
            .filter(|&(_, power)| power != 0)
            .map(|(name, power)| DimExpr::base(name).pow(power))
            .collect(),
    )
}

struct GasVelocity(f64);

impl NativeArg for GasVelocity {
    fn ty() -> TypeExpr {
        TypeExpr::Float.with_dim(gas_domain(1))
    }

    fn decode(
        value: &Value,
        index: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::same_dimension(&api::dimension(value)?, &gas_domain(1))?;
        numeric(value, index).map(Self)
    }
}

struct GasVelocities(Vec<f64>);

impl NativeArg for GasVelocities {
    fn ty() -> TypeExpr {
        api::tensor_type(1, GasVelocity::ty())
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let batch = api::MeasuredVector::<{ super::LIMIT }>::decode(value, index, ctx)?;
        api::same_dimension(&batch.dimension, &gas_domain(1))?;
        Ok(Self(batch.values))
    }
}

struct RiemannState<const BATCH: bool>(Value);

impl<const BATCH: bool> api::NativeOutput for RiemannState<BATCH> {
    fn ty() -> TypeExpr {
        TypeExpr::Record(
            GAS_COMPONENTS
                .into_iter()
                .enumerate()
                .map(|(component, name)| {
                    let quantity = TypeExpr::Float.with_dim(gas_domain(component));
                    let ty = if BATCH {
                        api::tensor_type(1, quantity)
                    } else {
                        quantity
                    };
                    (name.into(), ty)
                })
                .collect(),
        )
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        api::bounded(&self.0)?;
        Ok(self.0)
    }
}

struct PreparedRiemann {
    solution: gas::ExactRiemann,
    references: [(f64, DimExpr); 3],
}

impl NativeType for PreparedRiemann {
    const KIND: &'static str = "rssn/native/gas-exact-riemann/v1";
    const MAX_BYTES: usize = References::MAX_BYTES;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

const fn components(state: gas::Primitive) -> [f64; 3] {
    [state.rho, state.u, state.p]
}

fn state(v: &[f64]) -> Result<gas::Primitive, OpError> {
    if v.len() != 3 || v[0] <= 0.0 || v[2] <= 0.0 {
        return Err(err(
            "normalized primitive must be [positive density, velocity, positive pressure]",
        ));
    }
    Ok(gas::Primitive {
        rho: v[0],
        u: v[1],
        p: v[2],
    })
}
fn gamma(g: f64) -> Result<(), OpError> {
    if g <= 1.0 || g > 10.0 {
        return Err(err("gas gamma must be in (1,10]"));
    }
    Ok(())
}
fn limiter(s: &str) -> Result<gas::Limiter, OpError> {
    match s {
        | "None" => Ok(gas::Limiter::None),
        | "Minmod" => Ok(gas::Limiter::Minmod),
        | "VanLeer" => Ok(gas::Limiter::VanLeer),
        | _ => Err(err("limiter must be None, Minmod or VanLeer")),
    }
}
fn output(
    f: &References,
    x: Vec<f64>,
    w: Vec<gas::Primitive>,
) -> Result<Dynamic, OpError> {
    let refs = [
        f.derived(GAS_DOMAINS[0])?,
        f.derived(GAS_DOMAINS[1])?,
        f.derived(GAS_DOMAINS[2])?,
    ];
    Ok(Dynamic(api::record(vec![
        ("coordinates", x.encode()?),
        (
            "states",
            w.into_iter()
                .map(|w| vec![w.rho, w.u, w.p])
                .collect::<Vec<_>>()
                .encode()?,
        ),
        (
            "coordinate_reference",
            api::quantity(f.scales[0], f.dims[0].clone())?,
        ),
        (
            "density_reference",
            api::quantity(refs[0].0, refs[0].1.clone())?,
        ),
        (
            "velocity_reference",
            api::quantity(refs[1].0, refs[1].1.clone())?,
        ),
        (
            "pressure_reference",
            api::quantity(refs[2].0, refs[2].1.clone())?,
        ),
    ])))
}
object_binding!(gas_hllc,"num_grid_gas_hllc","crate::sim::gas_dynamics::hllc_flux","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64|{gamma(adiabatic)?;frame.with(ctx,|f|{let flux=gas::hllc_flux(state(&left)?,state(&right)?,adiabatic);let powers=[[-2,-1,1,0],[-1,-2,1,0],[0,-3,1,0]];let mut fields=Vec::new();for (i,name) in ["mass_flux","momentum_flux","energy_flux"].into_iter().enumerate(){let(s,d)=f.derived(powers[i])?;fields.push((name,api::quantity(flux[i]*s,d)?));}Ok(Dynamic(api::record(fields)))})});
object_binding!(gas_exact,"num_grid_gas_exact","crate::sim::gas_dynamics::exact_riemann","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64,similarity:f64)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,left:Vec<f64>,right:Vec<f64>,adiabatic:f64,similarity:f64|{gamma(adiabatic)?;frame.with(ctx,|f|{let w=gas::exact_riemann(state(&left)?,state(&right)?,adiabatic,similarity);output(f,vec![],vec![w])})});
object_binding!(
    gas_prepare, "num_grid_gas_prepare_riemann", "crate::sim::gas_dynamics::ExactRiemann::new",
    "Prepares the exact ideal-gas Riemann solution for repeated sampling. Each normalized input vector is [density, velocity, pressure], and gas density and pressure must be positive. The mechanical frame fixes the physical references of all three quantities.",
    true,
    (frame: Handle<References>, left: Vec<f64>, right: Vec<f64>, adiabatic: f64) -> NewObject<PreparedRiemann>,
    |ctx, frame, left, right, adiabatic| {
        gamma(adiabatic)?;
        let left = state(&left)?;
        let right = state(&right)?;
        frame.with(ctx, |frame| {
            let references = [frame.derived(GAS_DOMAINS[0])?, frame.derived(GAS_DOMAINS[1])?, frame.derived(GAS_DOMAINS[2])?];
            let dims = references.iter().map(|(_, dim)| dim.clone()).collect();
            let solution = gas::ExactRiemann::new(left, right, adiabatic);
            Ok(NewObject::new(PreparedRiemann { solution, references }, dims))
        })
    }
);
object_binding!(
    gas_sample, "num_grid_gas_sample_riemann", "crate::sim::gas_dynamics::ExactRiemann::sample",
    "Samples a prepared exact Riemann solution at a measured self-similar velocity x/t. Returns the measured density, velocity and pressure. Within a vacuum, density and pressure are zero and velocity is x/t.",
    false,
    (plan: Handle<PreparedRiemann>, velocity: GasVelocity) -> RiemannState<false>,
    |ctx, plan, velocity| plan.with(ctx, |plan| {
        let state = plan.solution.sample(api::scale_product(velocity.0, &[], &[plan.references[1].0])?);
        let fields = GAS_COMPONENTS.into_iter().zip(components(state)).zip(&plan.references)
            .map(|((name, value), (scale, dim))| Ok((name, api::quantity(api::scale_product(value, &[*scale], &[])?, dim.clone())?)))
            .collect::<Result<_, OpError>>()?;
        Ok(RiemannState(api::record(fields)))
    })
);
object_binding!(
    gas_sample_batch, "num_grid_gas_sample_riemann_batch", "crate::sim::gas_dynamics::ExactRiemann::sample_batch",
    "Samples up to 65536 measured self-similar velocities x/t. The record contains measured density, velocity and pressure vectors in the input order. Independent samples can use the active CPU pool.",
    false,
    (plan: Handle<PreparedRiemann>, velocities: GasVelocities) -> RiemannState<true>,
    |ctx, plan, velocities| plan.with(ctx, |plan| {
        let mut queries = velocities.0;
        let count = queries.len();
        for value in &mut queries {
            *value = api::scale_product(*value, &[], &[plan.references[1].0])?;
        }
        let mut columns = std::array::from_fn::<_, 3, _>(|_| Vec::with_capacity(count));
        for state in plan.solution.sample_batch(&queries) {
            for ((column, value), (scale, _)) in columns.iter_mut().zip(components(state)).zip(&plan.references) {
                column.push(Value::Float(api::scale_product(value, &[*scale], &[])?));
            }
        }
        let fields = GAS_COMPONENTS.into_iter().zip(columns).zip(&plan.references).map(|((name, data), (_, dim))| {
            (name, Value::Tensor(TensorValue { shape: vec![count], elem: TypeExpr::Float.with_dim(dim.clone()), data }))
        }).collect();
        Ok(RiemannState(api::record(fields)))
    })
);
object_binding!(gas_sod,"num_grid_gas_sod","crate::sim::gas_dynamics::sod_shock_tube","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,cells:usize,time:f64,reconstruction:String,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,cells:usize,time:f64,reconstruction:String,budget:usize|{if !(6..=1024).contains(&cells)||time<0.0||budget==0||budget>100000{return Err(err("gas grid/time/step budget invalid"));}let lim=limiter(&reconstruction)?;frame.with(ctx,|f|{let(x,rho,p)=crate::sim::models_geodesic_relativity::with_native_budget(budget,||gas::sod_shock_tube(cells,time,lim)).map_err(err)?;Ok(Dynamic(api::record(vec![("coordinates",x.encode()?),("density",rho.encode()?),("pressure",p.encode()?),("coordinate_reference",api::quantity(f.scales[0],f.dims[0].clone())?),("density_reference",{let(s,d)=f.derived([-3,0,1,0])?;api::quantity(s,d)?}),("pressure_reference",{let(s,d)=f.derived([-1,-2,1,0])?;api::quantity(s,d)?})])))})});
object_binding!(gas_solve,"num_grid_gas_solve","crate::sim::gas_dynamics::solve_euler_1d","Native ideal-gas method in an explicit mechanical reference frame. Numeric primitives are normalized density Mass/Length^3, velocity Length/Time and pressure Mass/(Length Time^2); numeric times and coordinates use the frame references. Distinct units remain separate in returned records. Native timesteps and callback evaluations have explicit budgets.",true,(frame:Handle<References>,density:Input,velocity:Input,pressure:Input,lo:f64,hi:f64,cells:usize,time:f64,adiabatic:f64,cfl:f64,reconstruction:String,boundary:String,budget:usize)->Dynamic,|ctx:&CallCtx<'_>,frame:Handle<References>,density:Input,velocity:Input,pressure:Input,lo:f64,hi:f64,cells:usize,time:f64,adiabatic:f64,cfl:f64,reconstruction:String,boundary:String,budget:usize|{gamma(adiabatic)?;if !(6..=1024).contains(&cells)||hi<=lo||time<0.0||cfl<=0.0||cfl>=1.0||budget==0||budget>100000{return Err(err("gas domain/time/grid/CFL/budget invalid"));}let lim=limiter(&reconstruction)?;let bc=match boundary.as_str(){"Reflective"=>gas::Boundary::Reflective,"Transmissive"=>gas::Boundary::Transmissive,_=>return Err(err("gas boundary must be Reflective or Transmissive"))};frame.with(ctx,|f|{let(rs,rd)=f.derived([-3,0,1,0])?;let(vs,vd)=f.derived([1,-1,0,0])?;let(ps,pd)=f.derived([-1,-2,1,0])?;let inputs=[f.dims[0].clone()];let rho=checked_snapshot(ctx,&density,cells*8,&inputs,&rd)?;let vel=checked_snapshot(ctx,&velocity,cells*8,&inputs,&vd)?;let pre=checked_snapshot(ctx,&pressure,cells*8,&inputs,&pd)?;let invalid=std::sync::atomic::AtomicBool::new(false);let result=crate::sim::models_geodesic_relativity::with_native_budget(budget,||gas::solve_euler_1d(|x|{let w=gas::Primitive{rho:rho.sample(&[x*f.scales[0]])/rs,u:vel.sample(&[x*f.scales[0]])/vs,p:pre.sample(&[x*f.scales[0]])/ps};
if !w.rho.is_finite()||w.rho<=0.0||!w.u.is_finite()||!w.p.is_finite()||w.p<=0.0{invalid.store(true,std::sync::atomic::Ordering::Relaxed);}w},lo,hi,cells,time,adiabatic,cfl,lim,bc));rho.finish()?;vel.finish()?;pre.finish()?;if invalid.load(std::sync::atomic::Ordering::Relaxed){return Err(err("gas initial callback returned invalid primitive"));}let(x,w)=result.map_err(err)?;output(f,x,w)})});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &gas_hllc::BINDING,
            &gas_exact::BINDING,
            &gas_prepare::BINDING,
            &gas_sample::BINDING,
            &gas_sample_batch::BINDING,
            &gas_sod::BINDING,
            &gas_solve::BINDING,
        ],
        "numerical_vector",
    )
}
