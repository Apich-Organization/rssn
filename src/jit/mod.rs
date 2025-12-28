#[cfg(feature = "jit")]
pub mod engine;
#[cfg(feature = "jit")]
pub mod instructions;

#[cfg(feature = "jit")]
pub use engine::JitEngine;
#[cfg(feature = "jit")]
pub use instructions::{Instruction, JitType};
