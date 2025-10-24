use crate::symbolic::core::Expr;
use serde::{Serialize, Deserialize};
use std::sync::{Arc, Mutex, Condvar};
use std::sync::atomic::{AtomicBool, Ordering};
use crate::compute::state::State;

// Using String for now, but could be a more complex type
pub type Value = String;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum ComputationStatus {
    Pending,
    Running,
    Paused,
    Completed,
    Failed(String),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComputationProgress {
    pub percentage: f32,
    pub description: String,
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Computation {
    pub id: String,
    pub expr: Arc<Expr>,
    pub status: ComputationStatus,
    pub progress: ComputationProgress,
    pub result: Option<Value>,
    pub state: State,
    #[serde(skip)]
    pub pause: Arc<(Mutex<bool>, Condvar)>,
    #[serde(skip)]
	pub cancel_signal: Arc<AtomicBool>, 
}
