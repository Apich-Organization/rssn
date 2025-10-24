use crate::compute::state::State;
use crate::compute::computation::ComputationProgress;

pub trait Computable {
    fn compute(&self, state: &mut State, progress: &mut ComputationProgress) -> Result<(), String>;
}
