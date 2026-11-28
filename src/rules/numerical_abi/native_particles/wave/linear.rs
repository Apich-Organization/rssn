//! Native adapters for a fixed linear Schrödinger propagator.
use super::*;
use crate::sim::split_step::LinearSplitPlan;

const MAX_POINTS: usize = 65_536;

struct PreparedLinear {
    plan: LinearSplitPlan,
    samples: WaveSamples,
}

impl NativeType for PreparedLinear {
    const KIND: &'static str = "rssn/native/sim/linear-split-plan/v1";
    const MAX_BYTES: usize = 5 * MAX_POINTS * size_of::<Complex64>() + 8192;

    fn validate(&self) -> Result<(), OpError> {
        self.samples.factor()?;
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        size_of::<Self>() + self.plan.heap_bytes()
    }
}

struct Potential {
    values: Vec<f64>,
    dimension: DimExpr,
}

impl NativeArg for Potential {
    fn ty() -> TypeExpr {
        api::tensor_type(1, variable("E", TypeExpr::Float))
    }

    fn decode(
        value: &Value,
        _: usize,
        _: &CallCtx<'_>,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(api::error("split-step potential must be a flat tensor"));
        };
        if tensor.shape != [tensor.data.len()] || tensor.data.len() > MAX_POINTS {
            return Err(api::error(
                "split-step potential exceeds the shape or work budget",
            ));
        }
        Ok(Self {
            dimension: api::dimension(value)?,
            values: <Vec<f64> as Arg>::decode(value, 0)?,
        })
    }
}

impl Potential {
    fn normalized(
        mut self,
        size: usize,
        physical: bool,
        scale: f64,
    ) -> Result<Vec<f64>, OpError> {
        if self.values.len() != size {
            return Err(api::error(
                "split-step potential must contain one sample per point",
            ));
        }
        let dimension = if physical {
            energy_dimension()
        } else {
            DimExpr::Dimless
        };
        api::same_dimension(&self.dimension, &dimension)?;
        if physical {
            for value in &mut self.values {
                *value = api::scale_without_underflow(*value, &[], &[scale])?;
            }
        }
        Ok(self.values)
    }
}

fn timestep(
    value: &Value,
    physical: bool,
    scale: f64,
) -> Result<f64, OpError> {
    let dimension = if physical {
        DimExpr::base("Time")
    } else {
        DimExpr::Dimless
    };
    api::same_dimension(&api::dimension(value)?, &dimension)?;
    let value = crate::rules::numerical_abi::numeric(value, 0)?;
    api::scale_without_underflow(value, &[], &[scale])
}

object_binding!(
    linear_plan_new,
    "num_native_split_linear_plan_new",
    "crate::sim::split_step::LinearSplitPlan::new",
    "Prepares a fixed linear Schrödinger propagator on a periodic grid with one to three axes, odd or even. Pass lengths and coordinate reference in Length units, mass in Mass units, flattened real potential in Energy units, and both timestep and time reference in Time units. Samples have the x axis varying fastest. Use Strang or Yoshida4, a negative timestep for inverse evolution, or zero for the identity.",
    false,
    (nodes: Vec<usize>, lengths: Input, mass: Input, potential: Potential, dt: Input,
        coordinate_reference: Input, time_reference: Input, scheme: String) -> NewObject<PreparedLinear>,
    |_, nodes, lengths, mass, potential, dt, coordinate_reference, time_reference, scheme| {
        let size = nodes.iter().try_fold(1_usize, |size, &n| size.checked_mul(n));
        let Some(size) = size.filter(|size| (1..=MAX_POINTS).contains(size)) else {
            return Err(api::error("linear split-step grid must have 1 to 65536 points"));
        };
        if !(1..=3).contains(&nodes.len()) {
            return Err(api::error("linear split-step grid must have one to three axes"));
        }
        let refs = PhysicalWaveReferences::new(&coordinate_reference.0, &time_reference.0, &mass.0)?;
        let samples = WaveSamples { dimensions: nodes.len(), size, reference: Some(refs.length) };
        let plan = LinearSplitPlan::new(
            &nodes,
            &refs.lengths(&lengths.0)?,
            refs.mass,
            &potential.normalized(size, true, refs.energy)?,
            timestep(&dt.0, true, refs.time)?,
            order(&scheme)?,
        ).map_err(api::error)?;
        let dimension = samples.dimension();
        Ok(NewObject::new(PreparedLinear { plan, samples }, vec![dimension]))
    }
);

object_binding!(
    prepare_linear,
    "num_native_split_prepare_linear",
    "crate::sim::split_step::SplitStep::prepare_linear",
    "Prepares a fixed linear propagator using an existing centered split-step solver's grid, mass and references. Use one real potential sample per grid point, a signed timestep, and Strang or Yoshida4. A normalized solver uses dimensionless values.",
    false,
    (solver: Handle<Wave>, potential: Potential, dt: Input, scheme: String) -> NewObject<PreparedLinear>,
    |ctx, solver, potential, dt, scheme| solver.with(ctx, |w| {
        let samples = w.samples();
        let plan = w.native.prepare_linear(
            &potential.normalized(samples.size, w.physical, w.energy)?,
            timestep(&dt.0, w.physical, w.time)?,
            order(&scheme)?,
        ).map_err(api::error)?;
        let dimension = samples.dimension();
        Ok(NewObject::new(PreparedLinear { plan, samples }, vec![dimension]))
    })
);

object_binding!(
    linear_plan_size,
    "num_native_split_linear_plan_size",
    "crate::sim::split_step::LinearSplitPlan::size",
    "Returns the number of points in the grid of a prepared linear split-step propagator.",
    true,
    (plan: Handle<PreparedLinear>) -> usize,
    |ctx, plan| plan.with(ctx, |prepared| Ok(prepared.plan.size()))
);

object_binding!(
    linear_plan_evolve,
    "num_native_split_linear_plan_evolve",
    "crate::sim::split_step::LinearSplitPlan::evolve",
    "Reuses the phases and Fourier plans to evolve a flat complex wavefunction. Physical samples carry Length^(-dimensions/2) and the result keeps that dimension. A zero step count or zero timestep gives the identity.",
    false,
    (plan: Handle<PreparedLinear>, psi: Input, steps: usize) -> Dynamic,
    |ctx, plan, psi, steps| plan.with(ctx, |prepared| {
        evolution(prepared.plan.size(), steps)?;
        let psi = prepared.samples.decode(&psi)?;
        let result = prepared.plan.evolve(&psi, steps).map_err(api::error)?;
        prepared.samples.encode(result)
    })
);

pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &linear_plan_new::BINDING,
            &prepare_linear::BINDING,
            &linear_plan_size::BINDING,
            &linear_plan_evolve::BINDING,
        ],
        "numerical_vector",
    )
}
