//! Handle-based FFI API for numerical FEA functions.

use crate::numerical::physics_fea;

// ============================================================================
// Material Functions
// ============================================================================

/// Creates steel material and returns shear modulus.
#[no_mangle]

pub extern "C" fn rssn_num_fea_material_steel_shear_modulus(
) -> f64 {

    physics_fea::Material::steel()
        .shear_modulus()
}

/// Creates aluminum material and returns shear modulus.
#[no_mangle]

pub extern "C" fn rssn_num_fea_material_aluminum_shear_modulus(
) -> f64 {

    physics_fea::Material::aluminum()
        .shear_modulus()
}

/// Creates copper material and returns shear modulus.
#[no_mangle]

pub extern "C" fn rssn_num_fea_material_copper_shear_modulus(
) -> f64 {

    physics_fea::Material::copper()
        .shear_modulus()
}

/// Computes shear modulus from Young's modulus and Poisson's ratio.
#[no_mangle]

pub extern "C" fn rssn_num_fea_shear_modulus(
    youngs_modulus : f64,
    poissons_ratio : f64,
) -> f64 {

    youngs_modulus
        / (2.0 * (1.0 + poissons_ratio))
}

/// Computes bulk modulus from Young's modulus and Poisson's ratio.
#[no_mangle]

pub extern "C" fn rssn_num_fea_bulk_modulus(
    youngs_modulus : f64,
    poissons_ratio : f64,
) -> f64 {

    youngs_modulus
        / (3.0
            * (1.0
                - 2.0 * poissons_ratio))
}

// ============================================================================
// 1D Element Functions
// ============================================================================

/// Computes and returns the stiffness value for a 1D linear element.
/// k = E * A / L
#[no_mangle]

pub extern "C" fn rssn_num_fea_linear_element_1d_stiffness(
    length : f64,
    youngs_modulus : f64,
    area : f64,
) -> f64 {

    youngs_modulus * area / length
}

// ============================================================================
// Stress Analysis Functions
// ============================================================================

/// Computes von Mises stress from plane stress components.
#[no_mangle]

pub extern "C" fn rssn_num_fea_von_mises_stress(
    sx : f64,
    sy : f64,
    txy : f64,
) -> f64 {

    physics_fea::TriangleElement2D::von_mises_stress(&[
        sx, sy, txy,
    ])
}

/// Computes maximum shear stress from principal stresses.
#[no_mangle]

pub extern "C" fn rssn_num_fea_max_shear_stress(
    sigma1 : f64,
    sigma2 : f64,
) -> f64 {

    physics_fea::max_shear_stress(
        sigma1,
        sigma2,
    )
}

/// Computes principal stresses from stress components.
/// Returns sigma1 in out_sigma1, sigma2 in out_sigma2, angle in out_angle.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fea_principal_stresses(
    sx : f64,
    sy : f64,
    txy : f64,
    out_sigma1 : *mut f64,
    out_sigma2 : *mut f64,
    out_angle : *mut f64,
) -> i32 {

    if out_sigma1.is_null()
        || out_sigma2.is_null()
        || out_angle.is_null()
    {

        return -1;
    }

    let (s1, s2, angle) =
        physics_fea::principal_stresses(
            &[sx, sy, txy],
        );

    *out_sigma1 = s1;

    *out_sigma2 = s2;

    *out_angle = angle;

    0
}

/// Computes safety factor based on von Mises criterion.
#[no_mangle]

pub extern "C" fn rssn_num_fea_safety_factor_von_mises(
    sx : f64,
    sy : f64,
    txy : f64,
    yield_strength : f64,
) -> f64 {

    physics_fea::safety_factor_von_mises(
        &[sx, sy, txy],
        yield_strength,
    )
}

// ============================================================================
// Thermal Element Functions
// ============================================================================

/// Computes the conductivity value for a 1D thermal element.
/// k = κ * A / L
#[no_mangle]

pub extern "C" fn rssn_num_fea_thermal_element_1d_conductivity(
    length : f64,
    conductivity : f64,
    area : f64,
) -> f64 {

    conductivity * area / length
}
