//! FFI Tests for numerical MD module.

use std::ffi::CString;

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_boltzmann_constant_handle() {

    let kb = rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_boltzmann_constant_si();

    assert!(kb > 1.38e-23 && kb < 1.39e-23);
}

#[test]

fn test_avogadro_number_handle() {

    let na = rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_avogadro_number();

    assert!(na > 6e23 && na < 6.1e23);
}

#[test]

fn test_temperature_unit_argon_handle() {

    let t = rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_temperature_unit_argon();

    assert!((t - 119.8).abs() < 1.0);
}

#[test]

fn test_minimum_image_1d_handle() {

    let d =
        rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_minimum_image_1d(8.0, 10.0);

    assert!((d + 2.0).abs() < 1e-10);
}

#[test]

fn test_apply_pbc_1d_handle() {

    let x = rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_apply_pbc_1d(11.0, 10.0);

    assert!((x - 1.0).abs() < 1e-10);

    let x2 = rssn::ffi_apis::numerical_physics_md_ffi::handle::rssn_num_md_apply_pbc_1d(-1.0, 10.0);

    assert!((x2 - 9.0).abs() < 1e-10);
}

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_lennard_jones_json() {

    let r_eq = 2.0_f64.powf(1.0 / 6.0);

    let input = format!(
        r#"{{"p1_position": [0.0, 0.0, 0.0], "p2_position": [{}, 0.0, 0.0], "epsilon": 1.0, "sigma": 1.0}}"#,
        r_eq
    );

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_lennard_jones_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let potential = parsed["ok"]["potential"]
            .as_f64()
            .unwrap();

        // At equilibrium, potential = -ε = -1.0
        assert!((potential + 1.0).abs() < 1e-10);
    }
}

#[test]

fn test_morse_json() {

    let input = r#"{"p1_position": [0.0, 0.0, 0.0], "p2_position": [1.0, 0.0, 0.0], "de": 1.0, "a": 1.0, "re": 1.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_morse_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let potential = parsed["ok"]["potential"]
            .as_f64()
            .unwrap();

        // At equilibrium, potential = 0
        assert!(potential.abs() < 1e-10);
    }
}

#[test]

fn test_harmonic_json() {

    let input = r#"{"p1_position": [0.0, 0.0, 0.0], "p2_position": [1.5, 0.0, 0.0], "k": 100.0, "r0": 1.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_harmonic_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let potential = parsed["ok"]["potential"]
            .as_f64()
            .unwrap();

        // V = 0.5 * 100 * 0.5² = 12.5
        assert!((potential - 12.5).abs() < 1e-10);
    }
}

#[test]

fn test_system_properties_json() {

    let input = r#"{"particles": [{"id": 0, "mass": 1.0, "position": [0.0, 0.0, 0.0], "velocity": [1.0, 0.0, 0.0]}, {"id": 1, "mass": 1.0, "position": [1.0, 0.0, 0.0], "velocity": [-1.0, 0.0, 0.0]}]}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_system_properties_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let ke = parsed["ok"]["kinetic_energy"]
            .as_f64()
            .unwrap();

        // KE = 0.5 + 0.5 = 1.0
        assert!((ke - 1.0).abs() < 1e-10);
    }
}

#[test]

fn test_create_cubic_lattice_json() {

    let input = r#"{"n_per_side": 2, "lattice_constant": 1.0, "mass": 1.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_create_cubic_lattice_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let particles = parsed["ok"]
            .as_array()
            .unwrap();

        assert_eq!(particles.len(), 8); // 2³ = 8
    }
}

#[test]

fn test_apply_pbc_json() {

    let input = r#"{"position": [11.0, -1.0, 5.0], "box_size": [10.0, 10.0, 10.0]}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_apply_pbc_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let wrapped = parsed["ok"]
            .as_array()
            .unwrap();

        assert!(
            (wrapped[0]
                .as_f64()
                .unwrap()
                - 1.0)
                .abs()
                < 1e-10
        );

        assert!(
            (wrapped[1]
                .as_f64()
                .unwrap()
                - 9.0)
                .abs()
                < 1e-10
        );
    }
}

#[test]

fn test_minimum_image_json() {

    let input = r#"{"position": [8.0, -8.0, 0.0], "box_size": [10.0, 10.0, 10.0]}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_md_ffi::json::rssn_num_md_minimum_image_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let r_mic = parsed["ok"]
            .as_array()
            .unwrap();

        assert!(
            (r_mic[0]
                .as_f64()
                .unwrap()
                + 2.0)
                .abs()
                < 1e-10
        );

        assert!(
            (r_mic[1]
                .as_f64()
                .unwrap()
                - 2.0)
                .abs()
                < 1e-10
        );
    }
}
