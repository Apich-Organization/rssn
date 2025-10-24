use crate::compute::computation::ComputationProgress;
use crate::compute::state::State;

pub trait Computable {
    fn compute(&self, state: &mut State, progress: &mut ComputationProgress) -> Result<(), String>;
}
