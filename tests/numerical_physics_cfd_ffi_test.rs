//! FFI Tests for numerical CFD module.

use std::ffi::CString;

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_air_kinematic_viscosity_handle() {

    let nu =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_air_kinematic_viscosity();

    // Air kinematic viscosity ≈ 1.5e-5 m²/s
    assert!(nu > 1e-5 && nu < 2e-5);
}

#[test]

fn test_water_kinematic_viscosity_handle() {

    let nu =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_water_kinematic_viscosity();

    // Water kinematic viscosity ≈ 1e-6 m²/s
    assert!(nu > 9e-7 && nu < 1.1e-6);
}

#[test]

fn test_air_prandtl_number_handle() {

    let pr = rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_air_prandtl_number();

    // Air Prandtl number ≈ 0.71
    assert!(pr > 0.7 && pr < 0.75);
}

#[test]

fn test_water_prandtl_number_handle() {

    let pr = rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_water_prandtl_number();

    // Water Prandtl number ≈ 7
    assert!(pr > 6.0 && pr < 8.0);
}

#[test]

fn test_reynolds_number_handle() {

    let re = rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_reynolds_number(
        1.0, 1.0, 1e-6,
    );

    assert!((re - 1e6).abs() < 1.0);
}

#[test]

fn test_mach_number_handle() {

    let ma =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_mach_number(250.0, 340.0);

    assert!((ma - 0.735).abs() < 0.01);
}

#[test]

fn test_froude_number_handle() {

    let fr = rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_froude_number(
        5.0, 10.0, 9.81,
    );

    assert!(fr > 0.4 && fr < 0.6);
}

#[test]

fn test_cfl_number_handle() {

    let cfl =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_cfl_number(1.0, 0.01, 0.1);

    assert!((cfl - 0.1).abs() < 1e-10);
}

#[test]

fn test_check_cfl_stability_handle() {

    let stable =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_check_cfl_stability(
            1.0, 0.01, 0.1, 1.0,
        );

    assert!(stable);

    let unstable =
        rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_check_cfl_stability(
            10.0, 0.1, 0.1, 1.0,
        );

    assert!(!unstable);
}

#[test]

fn test_diffusion_number_handle() {

    let r = rssn::ffi_apis::numerical_physics_cfd_ffi::handle::rssn_num_cfd_diffusion_number(
        0.01, 0.001, 0.01,
    );

    assert!((r - 0.1).abs() < 1e-10);
}

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_air_properties_json() {

    let input = r#"{}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_air_properties_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let pr = parsed["ok"]["prandtl_number"]
            .as_f64()
            .unwrap();

        assert!(pr > 0.7 && pr < 0.75);
    }
}

#[test]

fn test_water_properties_json() {

    let input = r#"{}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_water_properties_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let pr = parsed["ok"]["prandtl_number"]
            .as_f64()
            .unwrap();

        assert!(pr > 6.0 && pr < 8.0);
    }
}

#[test]

fn test_reynolds_number_json() {

    let input = r#"{"velocity": 1.0, "length": 1.0, "kinematic_viscosity": 1e-6}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_reynolds_number_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let re = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!((re - 1e6).abs() < 1.0);
    }
}

#[test]

fn test_cfl_number_json() {

    let input = r#"{"velocity": 1.0, "dt": 0.01, "dx": 0.1}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_cfl_number_json(
            c_input.as_ptr(),
        );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let cfl = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!((cfl - 0.1).abs() < 1e-10);
    }
}

#[test]

fn test_solve_advection_1d_json() {

    let input = r#"{"u0": [0.0, 0.0, 1.0, 1.0, 0.0, 0.0], "c": 1.0, "dx": 0.1, "dt": 0.01, "num_steps": 5}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_solve_advection_1d_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let data = parsed["ok"]
            .as_array()
            .unwrap();

        assert_eq!(data.len(), 6); // 5 steps + initial
    }
}

#[test]

fn test_solve_diffusion_1d_json() {

    let input = r#"{"u0": [0.0, 0.0, 1.0, 0.0, 0.0], "alpha": 0.01, "dx": 0.1, "dt": 0.01, "num_steps": 5}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_solve_diffusion_1d_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let data = parsed["ok"]
            .as_array()
            .unwrap();

        assert_eq!(data.len(), 6);
    }
}

#[test]

fn test_solve_burgers_1d_json() {

    let input =
        r#"{"u0": [0.0, 0.5, 1.0, 0.5, 0.0], "nu": 0.01, "dx": 0.1, "dt": 0.001, "num_steps": 5}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_physics_cfd_ffi::json::rssn_num_cfd_solve_burgers_1d_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed : serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let data = parsed["ok"]
            .as_array()
            .unwrap();

        assert_eq!(data.len(), 6);
    }
}
