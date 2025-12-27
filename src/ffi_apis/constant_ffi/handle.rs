//! Handle-based FFI API for constant module.
//!
//! This provides traditional C-style opaque pointer functions.

use std::ffi::CString;
use std::os::raw::c_char;

/// Returns the build date as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_build_date(
) -> *mut c_char {

    let date =
        crate::constant::get_build_date(
        );

    match CString::new(date) {
        | Ok(c_str) => c_str.into_raw(),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the commit SHA as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha(
) -> *mut c_char {

    let sha =
        crate::constant::get_commit_sha(
        );

    match CString::new(sha) {
        | Ok(c_str) => c_str.into_raw(),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the rustc version as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_rustc_version(
) -> *mut c_char {

    let version = crate::constant::get_rustc_version();

    match CString::new(version) {
        | Ok(c_str) => c_str.into_raw(),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the cargo target triple as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_cargo_target_triple(
) -> *mut c_char {

    let triple = crate::constant::get_cargo_target_triple();

    match CString::new(triple) {
        | Ok(c_str) => c_str.into_raw(),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the system info as a C string.
/// The caller must free the returned string using rssn_free_string.
#[no_mangle]

pub extern "C" fn rssn_get_system_info(
) -> *mut c_char {

    let info = crate::constant::get_system_info();

    match CString::new(info) {
        | Ok(c_str) => c_str.into_raw(),
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

macro_rules! gen_ffi_handle {
    ($ffi_name:ident, $internal_getter:path) => {
        #[no_mangle]
        pub extern "C" fn $ffi_name() -> *mut std::os::raw::c_char {
            let value = $internal_getter();

            // Convert to string (works for both f64 and &str)
            let val_str = value.to_string();

            match std::ffi::CString::new(val_str) {
                Ok(c_str) => c_str.into_raw(),
                Err(_) => std::ptr::null_mut(),
            }
        }
    };
}

// --- FFI Implementations ---

// // --- Build Info (Original Strings) ---
// gen_ffi_handle!(rssn_get_build_date, crate::constant::get_build_date);
// gen_ffi_handle!(rssn_get_commit_sha, crate::constant::get_commit_sha);
// gen_ffi_handle!(rssn_get_rustc_version, crate::constant::get_rustc_version);
// gen_ffi_handle!(rssn_get_cargo_target_triple, crate::constant::get_cargo_target_triple);
// gen_ffi_handle!(rssn_get_system_info, crate::constant::get_system_info);
// I will comment out these calls for now but personally I think it's a good idea to have them instead of the original functions.

// --- Fundamental Constants ---
gen_ffi_handle!(
    rssn_get_speed_of_light,
    crate::constant::get_speed_of_light
);

gen_ffi_handle!(rssn_get_planck_constant, crate::constant::get_planck_constant);

gen_ffi_handle!(rssn_get_reduced_planck_constant, crate::constant::get_reduced_planck_constant);

gen_ffi_handle!(rssn_get_elementary_charge, crate::constant::get_elementary_charge);

gen_ffi_handle!(rssn_get_gravitational_constant, crate::constant::get_gravitational_constant);

gen_ffi_handle!(rssn_get_fine_structure_constant, crate::constant::get_fine_structure_constant);

gen_ffi_handle!(rssn_get_inverse_fine_structure_constant, crate::constant::get_inverse_fine_structure_constant);

// --- Electromagnetic Constants ---
gen_ffi_handle!(rssn_get_vacuum_magnetic_permeability, crate::constant::get_vacuum_magnetic_permeability);

gen_ffi_handle!(rssn_get_vacuum_electric_permittivity, crate::constant::get_vacuum_electric_permittivity);

gen_ffi_handle!(rssn_get_josephson_constant, crate::constant::get_josephson_constant);

gen_ffi_handle!(rssn_get_von_klitzing_constant, crate::constant::get_von_klitzing_constant);

gen_ffi_handle!(rssn_get_magnetic_flux_quantum, crate::constant::get_magnetic_flux_quantum);

// --- Atomic & Particle Constants ---
gen_ffi_handle!(
    rssn_get_electron_mass,
    crate::constant::get_electron_mass
);

gen_ffi_handle!(rssn_get_electron_mass_u, crate::constant::get_electron_mass_u);

gen_ffi_handle!(
    rssn_get_proton_mass_kg,
    crate::constant::get_proton_mass_kg
);

gen_ffi_handle!(
    rssn_get_proton_mass_u,
    crate::constant::get_proton_mass_u
);

gen_ffi_handle!(
    rssn_get_neutron_mass_u,
    crate::constant::get_neutron_mass_u
);

gen_ffi_handle!(rssn_get_deuteron_mass_u, crate::constant::get_deuteron_mass_u);

gen_ffi_handle!(rssn_get_alpha_particle_mass_u, crate::constant::get_alpha_particle_mass_u);

gen_ffi_handle!(rssn_get_atomic_mass_constant, crate::constant::get_atomic_mass_constant);

gen_ffi_handle!(rssn_get_rydberg_constant, crate::constant::get_rydberg_constant);

gen_ffi_handle!(
    rssn_get_bohr_radius,
    crate::constant::get_bohr_radius
);

gen_ffi_handle!(
    rssn_get_hartree_energy,
    crate::constant::get_hartree_energy
);

gen_ffi_handle!(rssn_get_classical_electron_radius, crate::constant::get_classical_electron_radius);

gen_ffi_handle!(rssn_get_thomson_cross_section, crate::constant::get_thomson_cross_section);

// --- Magnetic Moments & Factors ---
gen_ffi_handle!(
    rssn_get_bohr_magneton,
    crate::constant::get_bohr_magneton
);

gen_ffi_handle!(rssn_get_nuclear_magneton, crate::constant::get_nuclear_magneton);

gen_ffi_handle!(rssn_get_electron_g_factor, crate::constant::get_electron_g_factor);

gen_ffi_handle!(
    rssn_get_muon_g_factor,
    crate::constant::get_muon_g_factor
);

gen_ffi_handle!(rssn_get_proton_magnetic_moment, crate::constant::get_proton_magnetic_moment);

gen_ffi_handle!(rssn_get_neutron_magnetic_moment, crate::constant::get_neutron_magnetic_moment);

gen_ffi_handle!(rssn_get_shielded_proton_gyromagnetic_ratio, crate::constant::get_shielded_proton_gyromagnetic_ratio);

// --- Physico-Chemical & Thermodynamic ---
gen_ffi_handle!(rssn_get_boltzmann_constant, crate::constant::get_boltzmann_constant);

gen_ffi_handle!(rssn_get_avogadro_constant, crate::constant::get_avogadro_constant);

gen_ffi_handle!(rssn_get_molar_gas_constant, crate::constant::get_molar_gas_constant);

gen_ffi_handle!(rssn_get_faraday_constant, crate::constant::get_faraday_constant);

gen_ffi_handle!(rssn_get_stefan_boltzmann_constant, crate::constant::get_stefan_boltzmann_constant);

gen_ffi_handle!(rssn_get_wien_displacement_constant, crate::constant::get_wien_displacement_constant);

gen_ffi_handle!(rssn_get_molar_volume_ideal_gas, crate::constant::get_molar_volume_ideal_gas);

gen_ffi_handle!(rssn_get_first_radiation_constant, crate::constant::get_first_radiation_constant);

gen_ffi_handle!(rssn_get_second_radiation_constant, crate::constant::get_second_radiation_constant);

// --- Ratios & Quotients ---
gen_ffi_handle!(rssn_get_proton_electron_mass_ratio, crate::constant::get_proton_electron_mass_ratio);

gen_ffi_handle!(rssn_get_muon_electron_mass_ratio, crate::constant::get_muon_electron_mass_ratio);

gen_ffi_handle!(rssn_get_neutron_proton_mass_ratio, crate::constant::get_neutron_proton_mass_ratio);

gen_ffi_handle!(rssn_get_electron_muon_mass_ratio, crate::constant::get_electron_muon_mass_ratio);

gen_ffi_handle!(rssn_get_deuteron_proton_mass_ratio, crate::constant::get_deuteron_proton_mass_ratio);

gen_ffi_handle!(rssn_get_electron_charge_to_mass_quotient, crate::constant::get_electron_charge_to_mass_quotient);

gen_ffi_handle!(rssn_get_muon_magnetic_moment, crate::constant::get_muon_magnetic_moment);

#[inline(always)]
#[no_mangle]

pub unsafe extern "C" fn rssn_free_string_constant(
    ptr: *mut c_char
) {

    if !ptr.is_null() {

        drop(
            std::ffi::CString::from_raw(
                ptr,
            ),
        );
    }
}
