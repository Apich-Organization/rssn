//! FFI Tests for numerical physics module.

use std::ffi::CString;

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_speed_of_light_handle() {

    let c = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_speed_of_light();

    assert!(
        (c - 299_792_458.0).abs() < 1.0
    );
}

#[test]

fn test_planck_constant_handle() {

    let h = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_planck_constant();

    assert!(
        h > 6.62e-34 && h < 6.63e-34
    );
}

#[test]

fn test_gravitational_constant_handle()
{

    let g =
        rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_gravitational_constant();

    assert!(
        g > 6.67e-11 && g < 6.68e-11
    );
}

#[test]

fn test_simple_harmonic_oscillator_handle()
 {

    let x =
        rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_simple_harmonic_oscillator(
            2.0, 1.0, 0.0, 0.0,
        );

    assert!((x - 2.0).abs() < 1e-10);
}

#[test]

fn test_coulomb_force_handle() {

    let f = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_coulomb_force(
        1.0, 1.0, 1.0,
    );

    assert!(f > 8e9 && f < 9e9);
}

#[test]

fn test_ideal_gas_pressure_handle() {

    let p = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_ideal_gas_pressure(
        1.0, 300.0, 1.0,
    );

    assert!(p > 2400.0 && p < 2600.0);
}

#[test]

fn test_lorentz_factor_handle() {

    let gamma =
        rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_lorentz_factor(1000.0);

    assert!(
        (gamma - 1.0).abs() < 1e-10
    );
}

#[test]

fn test_mass_energy_handle() {

    let e = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_mass_energy(1.0);

    assert!(e > 8e16 && e < 1e17);
}

#[test]

fn test_quantum_harmonic_oscillator_energy_handle()
 {

    let e = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_quantum_harmonic_oscillator_energy(0, 1.0);

    let hbar = 1.054_571_817e-34;

    assert!(
        (e - hbar * 0.5).abs() < 1e-40
    );
}

#[test]

fn test_hydrogen_energy_level_handle() {

    let e =
        rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_hydrogen_energy_level(1);

    assert!(e < 0.0);

    assert!(
        e > -2.2e-18 && e < -2.1e-18
    );
}

#[test]

fn test_de_broglie_wavelength_handle() {

    let lambda =
        rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_de_broglie_wavelength(
            1e-24,
        );

    assert!(
        lambda > 6e-10
            && lambda < 7e-10
    );
}

#[test]

fn test_photon_energy_handle() {

    let e = rssn::ffi_apis::numerical_physics_ffi::handle::rssn_num_physics_photon_energy(500e-9);

    assert!(e > 3e-19 && e < 5e-19);
}

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_simple_harmonic_oscillator_json()
 {

    let input = r#"{"amplitude": 2.0, "omega": 1.0, "phase": 0.0, "time": 0.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_simple_harmonic_oscillator_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let x = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(
            (x - 2.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_coulomb_force_json() {

    let input = r#"{"q1": 1.0, "q2": 1.0, "r": 1.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_coulomb_force_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let f = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(f > 8e9 && f < 9e9);
    }
}

#[test]

fn test_ideal_gas_pressure_json() {

    let input = r#"{"n": 1.0, "t": 300.0, "v": 1.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_ideal_gas_pressure_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let p = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(
            p > 2400.0 && p < 2600.0
        );
    }
}

#[test]

fn test_lorentz_factor_json() {

    let input =
        r#"{"velocity": 1000.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_lorentz_factor_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let gamma = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(
            (gamma - 1.0).abs() < 1e-10
        );
    }
}

#[test]

fn test_time_dilation_json() {

    let v = 0.9 * 299_792_458.0; // 0.9c
    let input = format!(
        r#"{{"proper_time": 1.0, "velocity": {}}}"#,
        v
    );

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_time_dilation_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let dilated = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(
            dilated > 2.2
                && dilated < 2.4
        );
    }
}

#[test]

fn test_mass_energy_json() {

    let input = r#"{"mass": 1.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_mass_energy_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let e = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(e > 8e16 && e < 1e17);
    }
}

#[test]

fn test_quantum_harmonic_oscillator_energy_json()
 {

    let input =
        r#"{"n": 0, "omega": 1.0}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_quantum_harmonic_oscillator_energy_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let e = parsed["ok"]
            .as_f64()
            .unwrap();

        let hbar = 1.054_571_817e-34;

        assert!(
            (e - hbar * 0.5).abs()
                < 1e-40
        );
    }
}

#[test]

fn test_hydrogen_energy_level_json() {

    let input = r#"{"n": 1}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_hydrogen_energy_level_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let e = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(e < 0.0);

        assert!(
            e > -2.2e-18
                && e < -2.1e-18
        );
    }
}

#[test]

fn test_photon_energy_json() {

    let input =
        r#"{"wavelength": 5e-7}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_ffi::json::rssn_num_physics_photon_energy_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str =
            std::ffi::CStr::from_ptr(
                result,
            )
            .to_string_lossy();

        let parsed: serde_json::Value =
            serde_json::from_str(
                &result_str,
            )
            .unwrap();

        let e = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(e > 3e-19 && e < 5e-19);
    }
}
