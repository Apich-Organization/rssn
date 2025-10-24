use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct State {
    // Placeholder for now. This will hold the intermediate state of a computation.
    pub intermediate_value: String,
}
