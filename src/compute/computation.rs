use std::sync::Arc;
use std::sync::Condvar;
use std::sync::Mutex;
use std::sync::atomic::AtomicBool;

use serde::Deserialize;
use serde::Serialize;

use crate::compute::state::State;
use crate::symbolic::core::Expr;

// Using String for now, but could be a more complex type
/// The type of the result of a computation. Currently a String.

pub type Value = String;

/// The status of a computation.
#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    PartialEq,
    Eq,
)]

pub enum ComputationStatus {
    /// The computation is pending execution.
    Pending,
    /// The computation is currently running.
    Running,
    /// The computation is paused.
    Paused,
    /// The computation has completed successfully.
    Completed,
    /// The computation failed with an error message.
    Failed(String),
}

/// Represents the progress of a computation.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub struct ComputationProgress {
    /// The percentage of completion (0.0 to 100.0).
    pub percentage: f32,
    /// A description of the current step.
    pub description: String,
}

/// Represents a computation task.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub struct Computation {
    /// A unique identifier for the computation.
    pub id: String,
    /// The expression being computed.
    pub expr: Arc<Expr>,
    /// The current status of the computation.
    pub status: ComputationStatus,
    /// The current progress of the computation.
    pub progress: ComputationProgress,
    /// The result of the computation, if available.
    pub result: Option<Value>,
    /// The state associated with the computation.
    pub state: State,
    /// Synchronization primitives for pausing/resuming.
    #[serde(
        skip,
        default = "default_pause"
    )]
    pub pause:
        Arc<(Mutex<bool>, Condvar)>,
    /// A signal to cancel the computation.
    #[serde(
        skip,
        default = "default_cancel_signal"
    )]
    pub cancel_signal: Arc<AtomicBool>,
}

pub(crate) fn default_pause()
-> Arc<(Mutex<bool>, Condvar)> {

    Arc::new((
        Mutex::new(false),
        Condvar::new(),
    ))
}

fn default_cancel_signal()
-> Arc<AtomicBool> {

    Arc::new(AtomicBool::new(
        false,
    ))
}
