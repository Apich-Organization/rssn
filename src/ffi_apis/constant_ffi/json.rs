//! JSON-based FFI API for constant module.
//!
//! This provides string-based serialization for easy language interop.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::to_c_string;

/// Build information structure for JSON serialization.
#[derive(
    Debug, Clone, Serialize, Deserialize,
)]

pub struct BuildInfo {
    /// The date when the library was built.
    pub build_date: String,
    /// The commit SHA of the build.
    pub commit_sha: String,
    /// The rustc version used for building.
    pub rustc_version: String,
    /// The cargo target triple.
    pub cargo_target_triple: String,
    /// System information.
    pub system_info: String,
}

/// Returns all build information as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_info_json(
) -> *mut c_char {

    let info = BuildInfo {
        build_date : crate::constant::get_build_date().to_string(),
        commit_sha : crate::constant::get_commit_sha().to_string(),
        rustc_version : crate::constant::get_rustc_version().to_string(),
        cargo_target_triple : crate::constant::get_cargo_target_triple().to_string(),
        system_info : crate::constant::get_system_info().to_string(),
    };

    match serde_json::to_string(&info) {
        | Ok(json) => to_c_string(json),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the build date as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_date_json(
) -> *mut c_char {

    let date =
        crate::constant::get_build_date(
        );

    match serde_json::to_string(&date) {
        | Ok(json) => to_c_string(json),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the commit SHA as a JSON string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha_json(
) -> *mut c_char {

    let sha =
        crate::constant::get_commit_sha(
        );

    match serde_json::to_string(&sha) {
        | Ok(json) => to_c_string(json),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

macro_rules! gen_ffi_json {
    (
        $ffi_name:ident,
        $internal_getter:path
    ) => {
        /// Generates a FFI function that retrieves a constant value,
        /// serializes it to a JSON string, and returns a C string pointer.
        /// The caller must free the returned string using rssn_free_string.
        #[no_mangle]

        pub extern "C" fn $ffi_name(
        ) -> *mut std::os::raw::c_char {

            let value =
                $internal_getter();

            // Serialize the value (f64 or &str) to a JSON string
            match serde_json::to_string(
                &value,
            ) {
                | Ok(json) => {
                    to_c_string(json)
                },
                | Err(_) => {
                    std::ptr::null_mut()
                },
            }
        }
    };
}

// --- FFI Implementations ---

// // --- Build Info (Original Strings) ---
// gen_ffi_json!(rssn_get_build_date_json, crate::constant::get_build_date);
// gen_ffi_json!(rssn_get_commit_sha_json, crate::constant::get_commit_sha);
// gen_ffi_json!(rssn_get_rustc_version_json, crate::constant::get_rustc_version);
// gen_ffi_json!(rssn_get_cargo_target_triple_json, crate::constant::get_cargo_target_triple);
// gen_ffi_json!(rssn_get_system_info_json, crate::constant::get_system_info);
// I will comment out these calls for now but personally I think it's a good idea to have them instead of the original functions.

// --- Fundamental Constants ---
gen_ffi_json!(
    rssn_get_speed_of_light_json,
    crate::constant::get_speed_of_light
);

gen_ffi_json!(rssn_get_planck_constant_json, crate::constant::get_planck_constant);

gen_ffi_json!(rssn_get_elementary_charge_json, crate::constant::get_elementary_charge);

gen_ffi_json!(rssn_get_gravitational_constant_json, crate::constant::get_gravitational_constant);

gen_ffi_json!(
    rssn_get_electron_mass_json,
    crate::constant::get_electron_mass
);

gen_ffi_json!(rssn_get_boltzmann_constant_json, crate::constant::get_boltzmann_constant);

gen_ffi_json!(rssn_get_avogadro_constant_json, crate::constant::get_avogadro_constant);

gen_ffi_json!(rssn_get_fine_structure_constant_json, crate::constant::get_fine_structure_constant);

gen_ffi_json!(rssn_get_rydberg_constant_json, crate::constant::get_rydberg_constant);

// --- Electromagnetic Constants ---
gen_ffi_json!(rssn_get_vacuum_electric_permittivity_json, crate::constant::get_vacuum_electric_permittivity);

gen_ffi_json!(rssn_get_vacuum_magnetic_permeability_json, crate::constant::get_vacuum_magnetic_permeability);

gen_ffi_json!(rssn_get_josephson_constant_json, crate::constant::get_josephson_constant);

gen_ffi_json!(rssn_get_von_klitzing_constant_json, crate::constant::get_von_klitzing_constant);

gen_ffi_json!(rssn_get_magnetic_flux_quantum_json, crate::constant::get_magnetic_flux_quantum);

// --- Atomic & Particle Constants ---
gen_ffi_json!(
    rssn_get_proton_mass_kg_json,
    crate::constant::get_proton_mass_kg
);

gen_ffi_json!(
    rssn_get_neutron_mass_u_json,
    crate::constant::get_neutron_mass_u
);

gen_ffi_json!(rssn_get_atomic_mass_constant_json, crate::constant::get_atomic_mass_constant);

gen_ffi_json!(rssn_get_proton_electron_mass_ratio_json, crate::constant::get_proton_electron_mass_ratio);

// --- Magnetic Moments & Factors ---
gen_ffi_json!(
    rssn_get_bohr_magneton_json,
    crate::constant::get_bohr_magneton
);

gen_ffi_json!(rssn_get_nuclear_magneton_json, crate::constant::get_nuclear_magneton);

gen_ffi_json!(rssn_get_electron_g_factor_json, crate::constant::get_electron_g_factor);

// --- Thermodynamic Constants ---
gen_ffi_json!(rssn_get_molar_gas_constant_json, crate::constant::get_molar_gas_constant);

gen_ffi_json!(rssn_get_faraday_constant_json, crate::constant::get_faraday_constant);

gen_ffi_json!(rssn_get_stefan_boltzmann_constant_json, crate::constant::get_stefan_boltzmann_constant);

// --- Planck Constants ---
gen_ffi_json!(rssn_get_reduced_planck_constant_json, crate::constant::get_reduced_planck_constant);

gen_ffi_json!(rssn_get_inverse_fine_structure_constant_json, crate::constant::get_inverse_fine_structure_constant);

// --- Atomic & Molecular Constants ---
gen_ffi_json!(
    rssn_get_bohr_radius_json,
    crate::constant::get_bohr_radius
);

gen_ffi_json!(
    rssn_get_hartree_energy_json,
    crate::constant::get_hartree_energy
);

gen_ffi_json!(rssn_get_electron_mass_u_json, crate::constant::get_electron_mass_u);

gen_ffi_json!(
    rssn_get_proton_mass_u_json,
    crate::constant::get_proton_mass_u
);

gen_ffi_json!(rssn_get_deuteron_mass_u_json, crate::constant::get_deuteron_mass_u);

gen_ffi_json!(rssn_get_alpha_particle_mass_u_json, crate::constant::get_alpha_particle_mass_u);

gen_ffi_json!(rssn_get_classical_electron_radius_json, crate::constant::get_classical_electron_radius);

gen_ffi_json!(rssn_get_thomson_cross_section_json, crate::constant::get_thomson_cross_section);

gen_ffi_json!(rssn_get_wien_displacement_constant_json, crate::constant::get_wien_displacement_constant);

gen_ffi_json!(rssn_get_first_radiation_constant_json, crate::constant::get_first_radiation_constant);

gen_ffi_json!(rssn_get_second_radiation_constant_json, crate::constant::get_second_radiation_constant);

gen_ffi_json!(
    rssn_get_muon_g_factor_json,
    crate::constant::get_muon_g_factor
);

gen_ffi_json!(
    rssn_get_muon_mass_u_json,
    crate::constant::get_muon_mass_u
);

// --- Ratios & Quotients ---
gen_ffi_json!(rssn_get_muon_electron_mass_ratio_json, crate::constant::get_muon_electron_mass_ratio);

gen_ffi_json!(rssn_get_neutron_proton_mass_ratio_json, crate::constant::get_neutron_proton_mass_ratio);

gen_ffi_json!(rssn_get_electron_muon_mass_ratio_json, crate::constant::get_electron_muon_mass_ratio);

gen_ffi_json!(rssn_get_proton_magnetic_moment_json, crate::constant::get_proton_magnetic_moment);

gen_ffi_json!(rssn_get_neutron_magnetic_moment_json, crate::constant::get_neutron_magnetic_moment);

gen_ffi_json!(rssn_get_muon_magnetic_moment_json, crate::constant::get_muon_magnetic_moment);

gen_ffi_json!(rssn_get_shielded_proton_gyromagnetic_ratio_json, crate::constant::get_shielded_proton_gyromagnetic_ratio);

gen_ffi_json!(rssn_get_molar_volume_ideal_gas_json, crate::constant::get_molar_volume_ideal_gas);

gen_ffi_json!(rssn_get_electron_charge_to_mass_quotient_json, crate::constant::get_electron_charge_to_mass_quotient);
