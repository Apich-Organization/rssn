use serde::{Deserialize, Serialize};

/// Represents the state of a computation.
///
/// This struct holds intermediate values and other context information
/// required during a computation.
#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct State {
    // Placeholder for now. This will hold the intermediate state of a computation.
    /// An intermediate value string.
    pub intermediate_value: String,
}

impl State {
    /// Creates a new, empty `State`.
    #[must_use]

    pub const fn new() -> Self {

        Self {
            intermediate_value: String::new(),
        }
    }
}

impl Default for State {
    fn default() -> Self {

        Self::new()
    }
}
