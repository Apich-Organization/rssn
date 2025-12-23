//! FFI Tests for numerical FEA module.

use std::ffi::CString;

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]
fn test_material_steel_shear_modulus_handle() {
    let g = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_material_steel_shear_modulus();
    assert!(g > 76e9 && g < 78e9);
}

#[test]
fn test_material_aluminum_shear_modulus_handle() {
    let g = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_material_aluminum_shear_modulus();
    assert!(g > 26e9 && g < 27e9);
}

#[test]
fn test_shear_modulus_handle() {
    let g = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_shear_modulus(200e9, 0.3);
    assert!((g - 200e9 / 2.6).abs() < 1e6);
}

#[test]
fn test_bulk_modulus_handle() {
    let k = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_bulk_modulus(200e9, 0.3);
    assert!(k > 165e9 && k < 168e9);
}

#[test]
fn test_linear_element_1d_stiffness_handle() {
    let k = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_linear_element_1d_stiffness(
        1.0, 200e9, 0.001,
    );
    assert!((k - 200e6).abs() < 1e-6);
}

#[test]
fn test_von_mises_stress_handle() {
    let vm = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_von_mises_stress(
        100e6, 0.0, 0.0,
    );
    assert!((vm - 100e6).abs() < 1e-6);
}

#[test]
fn test_max_shear_stress_handle() {
    let tau = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_max_shear_stress(
        100e6, -100e6,
    );
    assert!((tau - 100e6).abs() < 1e-6);
}

#[test]
fn test_principal_stresses_handle() {
    let mut sigma1 = 0.0;
    let mut sigma2 = 0.0;
    let mut angle = 0.0;
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_principal_stresses(
            100e6, 0.0, 0.0,
            &mut sigma1, &mut sigma2, &mut angle,
        );
        assert_eq!(result, 0);
        assert!((sigma1 - 100e6).abs() < 1e-6);
        assert!(sigma2.abs() < 1e-6);
    }
}

#[test]
fn test_safety_factor_handle() {
    let sf = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_safety_factor_von_mises(
        100e6, 0.0, 0.0, 250e6,
    );
    assert!((sf - 2.5).abs() < 1e-6);
}

#[test]
fn test_thermal_element_1d_conductivity_handle() {
    let k = rssn::ffi_apis::numerical_physics_fea_ffi::handle::rssn_num_fea_thermal_element_1d_conductivity(
        1.0, 50.0, 0.001,
    );
    assert!((k - 0.05).abs() < 1e-10);
}

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]
fn test_material_steel_json() {
    let input = r#"{}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_material_steel_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let g = parsed["ok"]["shear_modulus"].as_f64().unwrap();
        assert!(g > 76e9 && g < 78e9);
    }
}

#[test]
fn test_linear_element_1d_stiffness_json() {
    let input = r#"{"length": 1.0, "youngs_modulus": 200e9, "area": 0.001}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_linear_element_1d_stiffness_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let k = parsed["ok"].as_f64().unwrap();
        assert!((k - 200e6).abs() < 1e-6);
    }
}

#[test]
fn test_von_mises_stress_json() {
    let input = r#"{"sx": 100e6, "sy": 0.0, "txy": 0.0}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_von_mises_stress_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let vm = parsed["ok"].as_f64().unwrap();
        assert!((vm - 100e6).abs() < 1e-6);
    }
}

#[test]
fn test_principal_stresses_json() {
    let input = r#"{"sx": 100e6, "sy": 0.0, "txy": 0.0}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_principal_stresses_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let sigma1 = parsed["ok"]["sigma1"].as_f64().unwrap();
        assert!((sigma1 - 100e6).abs() < 1e-6);
    }
}

#[test]
fn test_safety_factor_json() {
    let input = r#"{"sx": 100e6, "sy": 0.0, "txy": 0.0, "yield_strength": 250e6}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_safety_factor_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let sf = parsed["ok"].as_f64().unwrap();
        assert!((sf - 2.5).abs() < 1e-6);
    }
}

#[test]
fn test_create_rectangular_mesh_json() {
    let input = r#"{"width": 1.0, "height": 1.0, "nx": 2, "ny": 2}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_create_rectangular_mesh_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let num_nodes = parsed["ok"]["num_nodes"].as_u64().unwrap();
        let num_elements = parsed["ok"]["num_elements"].as_u64().unwrap();
        assert_eq!(num_nodes, 9);
        assert_eq!(num_elements, 8);
    }
}

#[test]
fn test_beam_element_2d_stiffness_json() {
    let input = r#"{"length": 1.0, "youngs_modulus": 200e9, "area": 0.001, "moment_of_inertia": 1e-6, "angle": 0.0}"#;
    let c_input = CString::new(input).unwrap();
    
    unsafe {
        let result = rssn::ffi_apis::numerical_physics_fea_ffi::json::rssn_num_fea_beam_element_2d_stiffness_json(c_input.as_ptr());
        assert!(!result.is_null());
        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();
        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();
        let data = parsed["ok"].as_array().unwrap();
        assert_eq!(data.len(), 36); // 6x6 matrix
    }
}
