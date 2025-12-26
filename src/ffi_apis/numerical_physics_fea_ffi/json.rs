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
    youngs_modulus : f64,
    poissons_ratio : f64,
    density : f64,
    thermal_conductivity : f64,
    thermal_expansion : f64,
    yield_strength : f64,
}

#[derive(Serialize)]

struct MaterialOutput {
    shear_modulus : f64,
    bulk_modulus : f64,
}

#[derive(Deserialize)]

struct LinearElement1DInput {
    length : f64,
    youngs_modulus : f64,
    area : f64,
}

#[derive(Deserialize)]

struct StressInput {
    sx : f64,
    sy : f64,
    txy : f64,
}

#[derive(Serialize)]

struct PrincipalStressOutput {
    sigma1 : f64,
    sigma2 : f64,
    angle : f64,
}

#[derive(Deserialize)]

struct SafetyFactorInput {
    sx : f64,
    sy : f64,
    txy : f64,
    yield_strength : f64,
}

#[derive(Deserialize)]

struct BeamElement2DInput {
    length : f64,
    youngs_modulus : f64,
    area : f64,
    moment_of_inertia : f64,
    angle : f64,
}

#[derive(Deserialize)]

struct ThermalElement1DInput {
    length : f64,
    conductivity : f64,
    area : f64,
}

#[derive(Deserialize)]

struct MeshInput {
    width : f64,
    height : f64,
    nx : usize,
    ny : usize,
}

#[derive(Serialize)]

struct MeshOutput {
    num_nodes : usize,
    num_elements : usize,
    nodes : Vec<(f64, f64)>,
    elements : Vec<[usize; 3]>,
}

// ============================================================================
// Material Functions
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_material_properties_json(
    input : *const c_char
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

    let mat = physics_fea::Material::new(
        input.youngs_modulus,
        input.poissons_ratio,
        input.density,
        input.thermal_conductivity,
        input.thermal_expansion,
        input.yield_strength,
    );

    let output = MaterialOutput {
        shear_modulus : mat.shear_modulus(),
        bulk_modulus : mat.bulk_modulus(),
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(output),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_material_steel_json(_input : *const c_char) -> *mut c_char {

    let mat = physics_fea::Material::steel();

    let output = MaterialOutput {
        shear_modulus : mat.shear_modulus(),
        bulk_modulus : mat.bulk_modulus(),
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(output),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Element Functions
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_linear_element_1d_stiffness_json(
    input : *const c_char
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

    let element = physics_fea::LinearElement1D {
        length : input.length,
        youngs_modulus : input.youngs_modulus,
        area : input.area,
    };

    let stiffness = element.youngs_modulus * element.area / element.length;

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(stiffness),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_beam_element_2d_stiffness_json(
    input : *const c_char
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

    let beam = physics_fea::BeamElement2D::new(
        input.length,
        input.youngs_modulus,
        input.area,
        input.moment_of_inertia,
        input.angle,
    );

    let k = beam.global_stiffness_matrix();

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(k.data()),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_thermal_element_1d_conductivity_json(
    input : *const c_char
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

    let conductivity = input.conductivity * input.area / input.length;

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(conductivity),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Stress Analysis Functions
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_von_mises_stress_json(input : *const c_char) -> *mut c_char {

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
        serde_json::to_string(&FfiResult {
            ok : Some(vm),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_principal_stresses_json(
    input : *const c_char
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

    let (sigma1, sigma2, angle) = physics_fea::principal_stresses(&[
        input.sx,
        input.sy,
        input.txy,
    ]);

    let output = PrincipalStressOutput {
        sigma1,
        sigma2,
        angle,
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(output),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_safety_factor_json(input : *const c_char) -> *mut c_char {

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
        serde_json::to_string(&FfiResult {
            ok : Some(sf),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// ============================================================================
// Mesh Functions
// ============================================================================

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_create_rectangular_mesh_json(
    input : *const c_char
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
        num_nodes : nodes.len(),
        num_elements : elements.len(),
        nodes : nodes
            .iter()
            .map(|n| (n.x, n.y))
            .collect(),
        elements,
    };

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(output),
            err : None::<String>,
        })
        .unwrap(),
    )
}
