//! WGSL source from a [`KernelIr`]: pure string generation, no GPU needed.
//!
//! This is the WGSL face of the shared emitter in [`super::emit`], which
//! documents the interface of the generated module (bindings 0 to 4), the
//! precision and the deliberate differences from the CPU. The names here
//! predate that module and are kept because the `gpu-wgpu` device and
//! downstream WebGPU embedders use them.
//!
//! One invocation computes one row. The dispatch is two-dimensional so that
//! more than `65535 * workgroup_size` rows fit:
//! `row = gid.y * (num_workgroups.x * workgroup_size) + gid.x`; the host
//! dispatches `ceil(rows / workgroup_size)` workgroups spread over x and y.

use super::DeviceError;
use super::KernelIr;
use super::emit::EmitOptions;
use super::emit::Precision;
use super::emit::Target;
use super::emit::emit;

/// A source emitter for one target language.
///
/// [`super::emit::emit`] covers the languages this crate knows; implement the
/// trait for a target it does not (PTX, SPIR-V assembly, ...).
pub trait KernelEmitter {
    /// Target name (`wgsl`, `ptx`, ...).
    fn target(&self) -> &'static str;

    /// Emits the program.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] for a malformed program,
    /// [`DeviceError::Unsupported`] for something the target cannot express.
    fn emit(
        &self,
        ir: &KernelIr,
    ) -> Result<String, DeviceError>;
}

/// Any target of [`super::emit`] as a [`KernelEmitter`].
#[derive(Clone, Debug)]
pub struct SourceEmitter {
    /// Options.
    pub options: EmitOptions,
}

impl KernelEmitter for SourceEmitter {
    fn target(&self) -> &'static str {
        self.options.target.name()
    }

    fn emit(
        &self,
        ir: &KernelIr,
    ) -> Result<String, DeviceError> {
        emit(ir, &self.options)
    }
}

/// Options of the WGSL emitter.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WgslOptions {
    /// `@workgroup_size(n)`.
    pub workgroup_size: u32,
    /// Entry point name.
    pub entry: String,
}

impl Default for WgslOptions {
    fn default() -> Self {
        Self {
            workgroup_size: 64,
            entry: "main".into(),
        }
    }
}

/// The WGSL emitter.
#[derive(Clone, Debug, Default)]
pub struct WgslEmitter {
    /// Options.
    pub options: WgslOptions,
}

impl KernelEmitter for WgslEmitter {
    fn target(&self) -> &'static str {
        "wgsl"
    }

    fn emit(
        &self,
        ir: &KernelIr,
    ) -> Result<String, DeviceError> {
        emit_wgsl(ir, &self.options)
    }
}

/// A WGSL `f32` expression for a constant.
#[must_use]
pub fn f32_literal(c: f64) -> String {
    #[allow(clippy::cast_possible_truncation)]
    let x = c as f32;
    if x.is_nan() {
        "bitcast<f32>(0x7fc00000u)".into()
    } else if x == f32::INFINITY {
        "bitcast<f32>(0x7f800000u)".into()
    } else if x == f32::NEG_INFINITY {
        "bitcast<f32>(0xff800000u)".into()
    } else {
        format!("{x:?}")
    }
}

/// Emits the WGSL module.
///
/// # Errors
/// [`DeviceError::InvalidKernel`] when `ir` does not validate or the
/// workgroup size is zero.
pub fn emit_wgsl(
    ir: &KernelIr,
    options: &WgslOptions,
) -> Result<String, DeviceError> {
    emit(
        ir,
        &EmitOptions {
            target: Target::Wgsl,
            precision: Precision::F32,
            workgroup_size: options.workgroup_size,
            entry: options.entry.clone(),
        },
    )
}
