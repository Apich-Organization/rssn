//! Measured native Schrödinger calls and reusable fixed-potential phases.
use super::Arg;
use super::CallCtx;
use super::Complex64;
use super::DimExpr;
use super::Dynamic;
use super::ENERGY;
use super::Handle;
use super::Input;
use super::Installer;
use super::MASS;
use super::MAX_GRID_POINTS;
use super::NativeArg;
use super::NativeType;
use super::NewObject;
use super::OpError;
use super::RuleError;
use super::TIME;
use super::TensorValue;
use super::TypeExpr;
use super::Value;
use super::WaveSamples;
use super::api;
use super::bounded;
use super::directory;
use super::energy_dimension;
use super::object_binding;
use super::output_density;
use super::read;
use super::si_dimension;
use crate::rules::numerical_abi::native_object_api::MeasuredVector;
use crate::rules::numerical_abi::native_object_api::NativeOutput;
use crate::sim::models_schrodinger_quantum::SchrodingerParameters;
use crate::sim::models_schrodinger_quantum::SchrodingerPlan2D;
use crate::sim::models_schrodinger_quantum::{
    self as quantum,
};

fn wave_dimension() -> DimExpr {
    DimExpr::base("Length").pow(-1)
}

struct Potential(Vec<f64>);

impl NativeArg for Potential {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Float.with_dim(energy_dimension()))
    }

    fn decode(
        value: &Value,
        index: usize,
        ctx: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        api::bounded(value)?;
        let vector = <MeasuredVector<MAX_GRID_POINTS> as NativeArg>::decode(value, index, ctx)?;
        api::same_dimension(&vector.dimension, &energy_dimension())?;
        Ok(Self(vector.values))
    }
}

const fn wave_samples(size: usize) -> WaveSamples {
    WaveSamples {
        dimensions: 2,
        size,
        reference: Some(1.0),
    }
}

struct InitialWave(Vec<Complex64>);

impl NativeArg for InitialWave {
    fn ty() -> TypeExpr {
        api::tensor_type(1, TypeExpr::Complex.with_dim(wave_dimension()))
    }

    fn decode(
        value: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(api::error(
                "Schrodinger wavefunction must be a flat complex tensor",
            ));
        };
        let size = tensor.data.len();
        if size > MAX_GRID_POINTS {
            return Err(api::error(
                "Schrodinger wavefunction exceeds the grid budget",
            ));
        }
        api::bounded(value)?;
        Ok(Self(wave_samples(size).decode_value(value)?))
    }
}

struct TerminalWave {
    time: f64,
    psi: Vec<Complex64>,
}

impl NativeOutput for TerminalWave {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            (
                "time".into(),
                TypeExpr::Float.with_dim(DimExpr::base("Time")),
            ),
            ("final_psi".into(), InitialWave::ty()),
        ])
    }

    fn encode(
        self,
        _: &CallCtx<'_>,
    ) -> Result<Value, OpError> {
        let psi = wave_samples(self.psi.len()).encode(self.psi)?.0;
        Ok(api::record(vec![
            ("time", api::quantity(self.time, DimExpr::base("Time"))?),
            ("final_psi", psi),
        ]))
    }
}

struct PreparedSchrodinger {
    plan: SchrodingerPlan2D,
    nodes: [usize; 2],
    dt: f64,
}

impl NativeType for PreparedSchrodinger {
    const KIND: &'static str = "rssn/native/sim/schrodinger-plan-2d/v1";
    const MAX_BYTES: usize = 3 * MAX_GRID_POINTS * size_of::<Complex64>() + 4096;

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.plan.heap_bytes()
    }
}

fn parameters(
    [nx, ny]: [usize; 2],
    [lengths, dt, hbar, mass]: [&Input; 4],
    steps: usize,
    potential: Vec<f64>,
) -> Result<SchrodingerParameters, OpError> {
    api::same_dimension(&api::dimension(&lengths.0)?, &DimExpr::base("Length"))?;
    let [lx, ly] = <[f64; 2] as Arg>::decode(&lengths.0, 0)?;
    Ok(SchrodingerParameters {
        nx,
        ny,
        lx,
        ly,
        time_steps: steps,
        dt: read(dt, TIME)?,
        hbar: read(hbar, [1, 2, -1, 0])?,
        mass: read(mass, MASS)?,
        potential,
    })
}

object_binding!(
    plan_new, "num_native_quantum_schrodinger_plan_new",
    "crate::sim::models_schrodinger_quantum::SchrodingerPlan2D::new",
    "Prepares periodic two-dimensional Strang split-step evolution for repeated wavefunctions. Supply positive lengths, mass and reduced Planck constant in physical units, a signed measured time increment, and one energy sample per grid point in x-fast order. Every grid extent is positive. A zero increment gives the identity.",
    false,
    (nodes: [usize; 2], lengths: Input, dt: Input, hbar: Input, mass: Input, potential: Potential) -> NewObject<PreparedSchrodinger>,
    |_, nodes, lengths, dt, hbar, mass, potential| {
        bounded(nodes, 1, 1)?;
        let params = parameters(nodes, [&lengths, &dt, &hbar, &mass], 0, potential.0)?;
        let plan = SchrodingerPlan2D::new(&params).map_err(api::error)?;
        Ok(NewObject::new(PreparedSchrodinger { plan, nodes, dt: params.dt }, vec![wave_dimension()]))
    }
);

object_binding!(
    plan_evolve, "num_native_quantum_schrodinger_plan_evolve",
    "crate::sim::models_schrodinger_quantum::SchrodingerPlan2D::evolve",
    "Evolves a complex wavefunction with inverse-length amplitudes on a prepared periodic grid. The input and final_psi are flattened with x varying fastest. Returns the final wavefunction and time equal to steps multiplied by the signed time increment. Zero steps preserve the initial wave. Uniform potentials can evolve directly to the final state.",
    false,
    (plan: Handle<PreparedSchrodinger>, initial_psi: InitialWave, steps: usize) -> TerminalWave,
    |ctx, plan, initial_psi, steps| plan.with(ctx, |prepared| {
        bounded(prepared.nodes, steps.max(1), 1)?;
        let time = prepared.dt * steps as f64;
        if !time.is_finite() {
            return Err(api::error("Schrodinger final time is not representable"));
        }
        let psi = prepared.plan.evolve(&initial_psi.0, steps, None, |_, _| Ok(()))
            .map_err(api::error)?;
        Ok(TerminalWave { time, psi })
    })
);

object_binding!(
    run, "num_native_quantum_run_schrodinger",
    "crate::sim::models_schrodinger_quantum::run_schrodinger_simulation",
    "Advances the two-dimensional periodic Schrödinger model with a positive time increment. Input lengths, mass, reduced Planck constant and potential use SI dimensions. Returns inverse-length final amplitudes, inverse-square-length density after steps 1, 11, 21 and so on, and their measured times. Samples have the x coordinate varying fastest.",
    false,
    (nodes: [usize; 2], lengths: Input, dt: Input, steps: usize, hbar: Input, mass: Input, potential: Input, initial_psi: Input) -> Dynamic,
    |_, nodes, lengths, dt, steps, hbar, mass, potential, initial_psi| {
        let count = bounded(nodes, steps, 2)?;
        if count.checked_mul(steps.div_ceil(10)).is_none_or(|n| n > api::VALUE_LIMIT) {
            return Err(api::error("Schrodinger density history exceeds the output budget"));
        }
        api::same_dimension(&api::dimension(&potential.0)?, &si_dimension(ENERGY)?)?;
        let potential = <Vec<f64> as Arg>::decode(&potential.0, 0)?;
        let params = parameters(nodes, [&lengths, &dt, &hbar, &mass], steps, potential)?;
        if params.dt <= 0.0 {
            return Err(api::error("Schrodinger sampled timestep must be positive"));
        }
        let samples = wave_samples(count);
        let mut psi = samples.decode(&initial_psi)?;
        let frames = quantum::run_schrodinger_simulation(&params, &mut psi).map_err(api::error)?;
        let n = frames.len();
        let density = output_density(frames.into_iter().flatten(), vec![n, nodes[1], nodes[0]])?;
        let time_dim = DimExpr::base("Time");
        let times = (0..steps).step_by(10)
            .map(|i| api::quantity((i + 1) as f64 * params.dt, time_dim.clone()))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Dynamic(api::record(vec![
            ("density_frames", density),
            ("final_psi", samples.encode(psi)?.0),
            ("snapshot_times", Value::Tensor(TensorValue {
                shape: vec![times.len()], data: times, elem: TypeExpr::Float.with_dim(time_dim),
            })),
        ])))
    }
);

object_binding!(
    double_slit, "num_native_quantum_double_slit_scenario",
    "crate::sim::models_schrodinger_quantum::simulate_double_slit_scenario",
    "IO: execute the original reduced-unit 256 by 256, 300-step double-slit demonstration and write the NPY density. Provide an explicit grid-update budget.",
    true,
    (output_directory: String, grid_update_budget: usize) -> bool,
    |_, output_directory, grid_update_budget| {
        directory(&output_directory, grid_update_budget, 256 * 256 * 300)?;
        quantum::simulate_double_slit_scenario(std::path::Path::new(&output_directory)).map_err(api::error)?;
        Ok(true)
    }
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &run::BINDING,
            &double_slit::BINDING,
            &plan_new::BINDING,
            &plan_evolve::BINDING,
        ],
        "numerical_vector",
    )
}
