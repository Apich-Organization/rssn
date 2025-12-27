//! Bincode-based FFI API for constant module.
//!
//! This provides binary serialization for high-performance interop.

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::constant_ffi::json::BuildInfo;

/// Returns all build information as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_build_info_bincode(
) -> BincodeBuffer {

    let info = BuildInfo {
        build_date : crate::constant::get_build_date().to_string(),
        commit_sha : crate::constant::get_commit_sha().to_string(),
        rustc_version : crate::constant::get_rustc_version().to_string(),
        cargo_target_triple : crate::constant::get_cargo_target_triple().to_string(),
        system_info : crate::constant::get_system_info().to_string(),
    };

    match bincode_next::serde::encode_to_vec(
        &info,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Returns the build date as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_build_date_bincode(
) -> BincodeBuffer {

    let date =
        crate::constant::get_build_date(
        );

    match bincode_next::serde::encode_to_vec(
        &date,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Returns the commit SHA as a bincode_next buffer.
/// The caller must free the returned buffer using rssn_free_bincode_buffer.
#[no_mangle]

pub extern "C" fn rssn_get_commit_sha_bincode(
) -> BincodeBuffer {

    let sha =
        crate::constant::get_commit_sha(
        );

    match bincode_next::serde::encode_to_vec(
        &sha,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

macro_rules! gen_ffi_bincode {
    ($ffi_name:ident, $internal_getter:path) => {
        #[no_mangle]
        pub extern "C" fn $ffi_name() -> BincodeBuffer {
            let value = $internal_getter();

            match bincode_next::serde::encode_to_vec(
                &value,
                bincode_next::config::standard(),
            ) {
                Ok(bytes) => BincodeBuffer::from_vec(bytes),
                Err(_) => BincodeBuffer::empty(),
            }
        }
    };
}

// --- FFI Implementations ---

// // --- Build Info ---
// gen_ffi_bincode!(rssn_get_build_date_bincode, crate::constant::get_build_date);
// gen_ffi_bincode!(rssn_get_commit_sha_bincode, crate::constant::get_commit_sha);
// gen_ffi_bincode!(rssn_get_rustc_version_bincode, crate::constant::get_rustc_version);
// gen_ffi_bincode!(rssn_get_cargo_target_triple_bincode, crate::constant::get_cargo_target_triple);
// gen_ffi_bincode!(rssn_get_system_info_bincode, crate::constant::get_system_info);
// I will comment out these calls for now but personally I think it's a good idea to have them instead of the original functions.

// --- Fundamental & Electromagnetic ---
gen_ffi_bincode!(rssn_get_speed_of_light_bincode, crate::constant::get_speed_of_light);
gen_ffi_bincode!(rssn_get_planck_constant_bincode, crate::constant::get_planck_constant);
gen_ffi_bincode!(rssn_get_reduced_planck_constant_bincode, crate::constant::get_reduced_planck_constant);
gen_ffi_bincode!(rssn_get_elementary_charge_bincode, crate::constant::get_elementary_charge);
gen_ffi_bincode!(rssn_get_gravitational_constant_bincode, crate::constant::get_gravitational_constant);
gen_ffi_bincode!(rssn_get_fine_structure_constant_bincode, crate::constant::get_fine_structure_constant);
gen_ffi_bincode!(rssn_get_inverse_fine_structure_constant_bincode, crate::constant::get_inverse_fine_structure_constant);
gen_ffi_bincode!(rssn_get_vacuum_magnetic_permeability_bincode, crate::constant::get_vacuum_magnetic_permeability);
gen_ffi_bincode!(rssn_get_vacuum_electric_permittivity_bincode, crate::constant::get_vacuum_electric_permittivity);
gen_ffi_bincode!(rssn_get_josephson_constant_bincode, crate::constant::get_josephson_constant);
gen_ffi_bincode!(rssn_get_von_klitzing_constant_bincode, crate::constant::get_von_klitzing_constant);
gen_ffi_bincode!(rssn_get_magnetic_flux_quantum_bincode, crate::constant::get_magnetic_flux_quantum);

// --- Atomic & Particle masses ---
gen_ffi_bincode!(rssn_get_electron_mass_bincode, crate::constant::get_electron_mass);
gen_ffi_bincode!(rssn_get_electron_mass_u_bincode, crate::constant::get_electron_mass_u);
gen_ffi_bincode!(rssn_get_proton_mass_kg_bincode, crate::constant::get_proton_mass_kg);
gen_ffi_bincode!(rssn_get_proton_mass_u_bincode, crate::constant::get_proton_mass_u);
gen_ffi_bincode!(rssn_get_neutron_mass_u_bincode, crate::constant::get_neutron_mass_u);
gen_ffi_bincode!(rssn_get_deuteron_mass_u_bincode, crate::constant::get_deuteron_mass_u);
gen_ffi_bincode!(rssn_get_alpha_particle_mass_u_bincode, crate::constant::get_alpha_particle_mass_u);
gen_ffi_bincode!(rssn_get_atomic_mass_constant_bincode, crate::constant::get_atomic_mass_constant);
gen_ffi_bincode!(rssn_get_rydberg_constant_bincode, crate::constant::get_rydberg_constant);
gen_ffi_bincode!(rssn_get_bohr_radius_bincode, crate::constant::get_bohr_radius);
gen_ffi_bincode!(rssn_get_hartree_energy_bincode, crate::constant::get_hartree_energy);
gen_ffi_bincode!(rssn_get_classical_electron_radius_bincode, crate::constant::get_classical_electron_radius);
gen_ffi_bincode!(rssn_get_thomson_cross_section_bincode, crate::constant::get_thomson_cross_section);

// --- Magnetic & G-factors ---
gen_ffi_bincode!(rssn_get_bohr_magneton_bincode, crate::constant::get_bohr_magneton);
gen_ffi_bincode!(rssn_get_nuclear_magneton_bincode, crate::constant::get_nuclear_magneton);
gen_ffi_bincode!(rssn_get_electron_g_factor_bincode, crate::constant::get_electron_g_factor);
gen_ffi_bincode!(rssn_get_muon_g_factor_bincode, crate::constant::get_muon_g_factor);
gen_ffi_bincode!(rssn_get_proton_magnetic_moment_bincode, crate::constant::get_proton_magnetic_moment);
gen_ffi_bincode!(rssn_get_neutron_magnetic_moment_bincode, crate::constant::get_neutron_magnetic_moment);
gen_ffi_bincode!(rssn_get_shielded_proton_gyromagnetic_ratio_bincode, crate::constant::get_shielded_proton_gyromagnetic_ratio);

// --- Thermodynamic & Physicochemical ---
gen_ffi_bincode!(rssn_get_boltzmann_constant_bincode, crate::constant::get_boltzmann_constant);
gen_ffi_bincode!(rssn_get_avogadro_constant_bincode, crate::constant::get_avogadro_constant);
gen_ffi_bincode!(rssn_get_molar_gas_constant_bincode, crate::constant::get_molar_gas_constant);
gen_ffi_bincode!(rssn_get_faraday_constant_bincode, crate::constant::get_faraday_constant);
gen_ffi_bincode!(rssn_get_stefan_boltzmann_constant_bincode, crate::constant::get_stefan_boltzmann_constant);
gen_ffi_bincode!(rssn_get_wien_displacement_constant_bincode, crate::constant::get_wien_displacement_constant);
gen_ffi_bincode!(rssn_get_molar_volume_ideal_gas_bincode, crate::constant::get_molar_volume_ideal_gas);
gen_ffi_bincode!(rssn_get_first_radiation_constant_bincode, crate::constant::get_first_radiation_constant);
gen_ffi_bincode!(rssn_get_second_radiation_constant_bincode, crate::constant::get_second_radiation_constant);

// --- Ratios ---
gen_ffi_bincode!(rssn_get_proton_electron_mass_ratio_bincode, crate::constant::get_proton_electron_mass_ratio);
gen_ffi_bincode!(rssn_get_muon_electron_mass_ratio_bincode, crate::constant::get_muon_electron_mass_ratio);
gen_ffi_bincode!(rssn_get_neutron_proton_mass_ratio_bincode, crate::constant::get_neutron_proton_mass_ratio);
gen_ffi_bincode!(rssn_get_electron_muon_mass_ratio_bincode, crate::constant::get_electron_muon_mass_ratio);
gen_ffi_bincode!(rssn_get_deuteron_proton_mass_ratio_bincode, crate::constant::get_deuteron_proton_mass_ratio);
gen_ffi_bincode!(rssn_get_electron_charge_to_mass_quotient_bincode, crate::constant::get_electron_charge_to_mass_quotient);
gen_ffi_bincode!(rssn_get_muon_magnetic_moment_bincode, crate::constant::get_muon_magnetic_moment);