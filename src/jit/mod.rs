//! The JIT compiler for RSSN

// Unavoidable for functionality issues.
#![allow(unsafe_code)]

#[cfg(feature = "jit")]
pub mod engine;
#[cfg(feature = "jit")]
pub mod instructions;

#[cfg(feature = "jit")]
pub use engine::JitEngine;
#[cfg(feature = "jit")]
pub use instructions::Instruction;
#[cfg(feature = "jit")]
pub use instructions::JitType;
