use crate::compute::computation::ComputationProgress;
use crate::compute::state::State;

/// A trait for objects that can perform a computation.
///
/// Implementors of this trait represent a unit of work that can be executed
/// by the compute engine.

pub trait Computable {
    /// Performs the computation.
    ///
    /// # Arguments
    /// * `state` - The current state of the computation.
    /// * `progress` - A mutable reference to update the progress.
    ///
    /// # Returns
    /// * `Result<(), String>` - `Ok(())` if successful, or `Err(String)` with an error message.

    fn compute(
        &self,
        state: &mut State,
        progress: &mut ComputationProgress,
    ) -> Result<(), String>;
}
