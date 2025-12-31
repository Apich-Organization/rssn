//! Bincode-based FFI API for numerical FEA functions.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_fea;

#[derive(Deserialize)]

struct LinearElement1DInput {
    length: f64,
    youngs_modulus: f64,
    area: f64,
}

#[derive(Deserialize)]

struct StressInput {
    sx: f64,
    sy: f64,
    txy: f64,
}

#[derive(Serialize)]

struct PrincipalStressOutput {
    sigma1: f64,
    sigma2: f64,
    angle: f64,
}

/// Computes the axial stiffness of a 1D linear finite element using bincode serialization.
///
/// The axial stiffness represents the force-displacement relationship for a bar element:
/// k = (E × A) / L, where E is Young's modulus, A is cross-sectional area, and L is length.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `LinearElement1DInput` with:
///   - `length`: Element length L (m)
///   - `youngs_modulus`: Young's modulus E (Pa)
///   - `area`: Cross-sectional area A (m²)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The computed axial stiffness k (N/m)
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_linear_element_1d_stiffness_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : LinearElement1DInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let stiffness = input
        .youngs_modulus
        * input.area
        / input.length;

    to_bincode_buffer(&FfiResult {
        ok: Some(stiffness),
        err: None::<String>,
    })
}

/// Computes the von Mises equivalent stress from a 2D stress state using bincode serialization.
///
/// The von Mises stress is a scalar measure of stress intensity used in yield criteria:
/// `σ_vm` = √(`σ_x²` - `σ_xσ_y` + `σ_y²` + `3τ_xy²`)
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `StressInput` with:
///   - `sx`: Normal stress in x-direction `σ_x` (Pa)
///   - `sy`: Normal stress in y-direction `σ_y` (Pa)
///   - `txy`: Shear stress `τ_xy` (Pa)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The von Mises stress `σ_vm` (Pa)
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_von_mises_stress_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : StressInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let vm = physics_fea::TriangleElement2D::von_mises_stress(&[
        input.sx,
        input.sy,
        input.txy,
    ]);

    to_bincode_buffer(&FfiResult {
        ok: Some(vm),
        err: None::<String>,
    })
}

/// Computes the principal stresses and orientation from a 2D stress state using bincode serialization.
///
/// Principal stresses are the eigenvalues of the stress tensor, representing maximum and minimum
/// normal stresses. The angle indicates the orientation of the principal stress axes.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `StressInput` with:
///   - `sx`: Normal stress in x-direction `σ_x` (Pa)
///   - `sy`: Normal stress in y-direction `σ_y` (Pa)
///   - `txy`: Shear stress `τ_xy` (Pa)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<PrincipalStressOutput, String>` with either:
/// - `ok`: Object containing:
///   - `sigma1`: Maximum principal stress σ₁ (Pa)
///   - `sigma2`: Minimum principal stress σ₂ (Pa)
///   - `angle`: Orientation angle θ (radians)
/// - `err`: Error message if deserialization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_principal_stresses_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : StressInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<PrincipalStressOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let (sigma1, sigma2, angle) =
        physics_fea::principal_stresses(
            &[
                input.sx,
                input.sy,
                input.txy,
            ],
        );

    let output =
        PrincipalStressOutput {
            sigma1,
            sigma2,
            angle,
        };

    to_bincode_buffer(&FfiResult {
        ok: Some(output),
        err: None::<String>,
    })
}
