//! JSON-based FFI API for numerical FEA functions.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::physics_fea;

// ============================================================================
// Input/Output structs
// ============================================================================

#[derive(Deserialize)]

struct MaterialInput {
    youngs_modulus: f64,
    poissons_ratio: f64,
    density: f64,
    thermal_conductivity: f64,
    thermal_expansion: f64,
    yield_strength: f64,
}

#[derive(Serialize)]

struct MaterialOutput {
    shear_modulus: f64,
    bulk_modulus: f64,
}

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

#[derive(Deserialize)]

struct SafetyFactorInput {
    sx: f64,
    sy: f64,
    txy: f64,
    yield_strength: f64,
}

#[derive(Deserialize)]

struct BeamElement2DInput {
    length: f64,
    youngs_modulus: f64,
    area: f64,
    moment_of_inertia: f64,
    angle: f64,
}

#[derive(Deserialize)]

struct ThermalElement1DInput {
    length: f64,
    conductivity: f64,
    area: f64,
}

#[derive(Deserialize)]

struct MeshInput {
    width: f64,
    height: f64,
    nx: usize,
    ny: usize,
}

#[derive(Serialize)]

struct MeshOutput {
    num_nodes: usize,
    num_elements: usize,
    nodes: Vec<(f64, f64)>,
    elements: Vec<[usize; 3]>,
}

// ============================================================================
// Material Functions
// ============================================================================

/// Computes derived material properties from fundamental elastic constants using JSON serialization.
///
/// This function calculates shear modulus G = E/(2(1+ν)) and bulk modulus K = E/(3(1-2ν))
/// from Young's modulus and Poisson's ratio.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `youngs_modulus`: Young's modulus E (Pa)
///   - `poissons_ratio`: Poisson's ratio ν (dimensionless)
///   - `density`: Material density ρ (kg/m³)
///   - `thermal_conductivity`: Thermal conductivity k (W/(m·K))
///   - `thermal_expansion`: Thermal expansion coefficient α (1/K)
///   - `yield_strength`: Yield strength `σ_y` (Pa)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with:
/// - `shear_modulus`: G (Pa)
/// - `bulk_modulus`: K (Pa)
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_material_properties_json(
    input: *const c_char
) -> *mut c_char {

    let input : MaterialInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<MaterialOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let mat =
        physics_fea::Material::new(
            input.youngs_modulus,
            input.poissons_ratio,
            input.density,
            input.thermal_conductivity,
            input.thermal_expansion,
            input.yield_strength,
        );

    let output = MaterialOutput {
        shear_modulus: mat
            .shear_modulus(),
        bulk_modulus: mat
            .bulk_modulus(),
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Returns standard material properties for structural steel using JSON serialization.
///
/// Provides reference properties for mild steel (E ≈ 200 `GPa`, ν ≈ 0.3, ρ ≈ 7850 kg/m³).
///
/// # Arguments
///
/// * `_input` - Unused parameter for API consistency (can be null)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with steel properties:
/// - `shear_modulus`: G (Pa)
/// - `bulk_modulus`: K (Pa)
///
/// # Safety
///
/// This function is unsafe because it returns a raw pointer that the caller must free.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_material_steel_json(
    _input: *const c_char
) -> *mut c_char {

    let mat =
        physics_fea::Material::steel();

    let output = MaterialOutput {
        shear_modulus: mat
            .shear_modulus(),
        bulk_modulus: mat
            .bulk_modulus(),
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Element Functions
// ============================================================================

/// Computes the axial stiffness of a 1D linear finite element using JSON serialization.
///
/// The axial stiffness represents the force-displacement relationship for a bar element:
/// k = (E × A) / L, where E is Young's modulus, A is cross-sectional area, and L is length.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `length`: Element length L (m)
///   - `youngs_modulus`: Young's modulus E (Pa)
///   - `area`: Cross-sectional area A (m²)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed axial stiffness k (N/m).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_linear_element_1d_stiffness_json(
    input: *const c_char
) -> *mut c_char {

    let input : LinearElement1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let element =
        physics_fea::LinearElement1D {
            length: input.length,
            youngs_modulus: input
                .youngs_modulus,
            area: input.area,
        };

    let stiffness = element
        .youngs_modulus
        * element.area
        / element.length;

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(stiffness),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the global stiffness matrix for a 2D beam element using JSON serialization.
///
/// The beam element includes both axial and bending behavior. The stiffness matrix
/// is transformed from local to global coordinates using the specified angle.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `length`: Element length L (m)
///   - `youngs_modulus`: Young's modulus E (Pa)
///   - `area`: Cross-sectional area A (m²)
///   - `moment_of_inertia`: Second moment of area I (m⁴)
///   - `angle`: Element orientation angle θ (radians)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the 6×6 global stiffness matrix in row-major format.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_beam_element_2d_stiffness_json(
    input: *const c_char
) -> *mut c_char {

    let input : BeamElement2DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let beam =
        physics_fea::BeamElement2D::new(
            input.length,
            input.youngs_modulus,
            input.area,
            input.moment_of_inertia,
            input.angle,
        );

    let k =
        beam.global_stiffness_matrix();

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(k.data()),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the thermal conductivity coefficient for a 1D thermal element using JSON serialization.
///
/// The thermal conductivity represents the heat flow-temperature relationship:
/// `k_thermal` = (k × A) / L, where k is thermal conductivity, A is area, and L is length.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `length`: Element length L (m)
///   - `conductivity`: Material thermal conductivity k (W/(m·K))
///   - `area`: Cross-sectional area A (m²)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed thermal conductivity coefficient (W/K).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_thermal_element_1d_conductivity_json(
    input: *const c_char
) -> *mut c_char {

    let input : ThermalElement1DInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let conductivity =
        input.conductivity * input.area
            / input.length;

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(conductivity),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Stress Analysis Functions
// ============================================================================

/// Computes the von Mises equivalent stress from a 2D stress state using JSON serialization.
///
/// The von Mises stress is a scalar measure of stress intensity used in yield criteria:
/// `σ_vm` = √(`σ_x²` - `σ_xσ_y` + `σ_y²` + `3τ_xy²`)
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `sx`: Normal stress in x-direction `σ_x` (Pa)
///   - `sy`: Normal stress in y-direction `σ_y` (Pa)
///   - `txy`: Shear stress `τ_xy` (Pa)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the von Mises stress `σ_vm` (Pa).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_von_mises_stress_json(
    input: *const c_char
) -> *mut c_char {

    let input : StressInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let vm = physics_fea::TriangleElement2D::von_mises_stress(&[
        input.sx,
        input.sy,
        input.txy,
    ]);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(vm),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the principal stresses and orientation from a 2D stress state using JSON serialization.
///
/// Principal stresses are the eigenvalues of the stress tensor, representing maximum and minimum
/// normal stresses. The angle indicates the orientation of the principal stress axes.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `sx`: Normal stress in x-direction `σ_x` (Pa)
///   - `sy`: Normal stress in y-direction `σ_y` (Pa)
///   - `txy`: Shear stress `τ_xy` (Pa)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with:
/// - `sigma1`: Maximum principal stress σ₁ (Pa)
/// - `sigma2`: Minimum principal stress σ₂ (Pa)
/// - `angle`: Orientation angle θ (radians)
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_principal_stresses_json(
    input: *const c_char
) -> *mut c_char {

    let input : StressInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<PrincipalStressOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
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

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the factor of safety against yielding using the von Mises criterion and JSON serialization.
///
/// The safety factor is the ratio of yield strength to equivalent stress:
/// SF = `σ_yield` / `σ_vm`. A value > 1 indicates the material will not yield.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `sx`: Normal stress in x-direction `σ_x` (Pa)
///   - `sy`: Normal stress in y-direction `σ_y` (Pa)
///   - `txy`: Shear stress `τ_xy` (Pa)
///   - `yield_strength`: Material yield strength `σ_y` (Pa)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the computed safety factor (dimensionless).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_safety_factor_json(
    input: *const c_char
) -> *mut c_char {

    let input : SafetyFactorInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let sf = physics_fea::safety_factor_von_mises(
        &[
            input.sx,
            input.sy,
            input.txy,
        ],
        input.yield_strength,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(sf),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

// ============================================================================
// Mesh Functions
// ============================================================================

/// Generates a structured triangular mesh for a rectangular domain using JSON serialization.
///
/// Creates a finite element mesh by subdividing a rectangle into triangular elements.
/// The mesh is structured with regular spacing in both x and y directions.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `width`: Domain width (m)
///   - `height`: Domain height (m)
///   - `nx`: Number of divisions in x-direction
///   - `ny`: Number of divisions in y-direction
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult` with:
/// - `num_nodes`: Total number of mesh nodes
/// - `num_elements`: Total number of triangular elements
/// - `nodes`: Array of (x, y) node coordinates
/// - `elements`: Array of triangular element connectivity (3 node indices per element)
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_fea_create_rectangular_mesh_json(
    input: *const c_char
) -> *mut c_char {

    let input : MeshInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<MeshOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let (nodes, elements) = physics_fea::create_rectangular_mesh(
        input.width,
        input.height,
        input.nx,
        input.ny,
    );

    let output = MeshOutput {
        num_nodes: nodes.len(),
        num_elements: elements.len(),
        nodes: nodes
            .iter()
            .map(|n| (n.x, n.y))
            .collect(),
        elements,
    };

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(output),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
