//! FFI Tests for numerical computer graphics module.

use std::ffi::CString;

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_dot_product_json() {

    let input = r#"{"v1": {"x": 1.0, "y": 0.0, "z": 0.0}, "v2": {"x": 0.0, "y": 1.0, "z": 0.0}}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_dot_product_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        assert_eq!(parsed["ok"], 0.0);
    }
}

#[test]

fn test_cross_product_json() {

    let input = r#"{"v1": {"x": 1.0, "y": 0.0, "z": 0.0}, "v2": {"x": 0.0, "y": 1.0, "z": 0.0}}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_cross_product_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        assert_eq!(parsed["ok"]["z"], 1.0);
    }
}

#[test]

fn test_normalize_json() {

    let input = r#"{"x": 3.0, "y": 4.0, "z": 0.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_normalize_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let x = parsed["ok"]["x"]
            .as_f64()
            .unwrap();

        let y = parsed["ok"]["y"]
            .as_f64()
            .unwrap();

        assert!((x - 0.6).abs() < 1e-10);

        assert!((y - 0.8).abs() < 1e-10);
    }
}

#[test]

fn test_magnitude_json() {

    let input = r#"{"x": 3.0, "y": 4.0, "z": 0.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result =
            rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_magnitude_json(
                c_input.as_ptr(),
            );

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        assert_eq!(parsed["ok"], 5.0);
    }
}

#[test]

fn test_translation_matrix_json() {

    let input = r#"{"dx": 1.0, "dy": 2.0, "dz": 3.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_translation_matrix_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let data = parsed["ok"]
            .as_array()
            .unwrap();

        assert_eq!(data.len(), 16);

        // Check translation components (indices 3, 7, 11 in row-major order)
        assert_eq!(data[3], 1.0);

        assert_eq!(data[7], 2.0);

        assert_eq!(data[11], 3.0);
    }
}

#[test]

fn test_rotation_matrix_x_json() {

    let input = r#"{"angle": 0.0}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_rotation_matrix_x_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let data = parsed["ok"]
            .as_array()
            .unwrap();

        // Identity matrix for angle = 0
        assert_eq!(data[0], 1.0);

        assert_eq!(data[5], 1.0);

        assert_eq!(data[10], 1.0);

        assert_eq!(data[15], 1.0);
    }
}

#[test]

fn test_quaternion_multiply_json() {

    let input = r#"{"q1": {"w": 1.0, "x": 0.0, "y": 0.0, "z": 0.0}, "q2": {"w": 1.0, "x": 0.0, "y": 0.0, "z": 0.0}}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_quaternion_multiply_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        // Identity * Identity = Identity
        assert_eq!(parsed["ok"]["w"], 1.0);

        assert_eq!(parsed["ok"]["x"], 0.0);
    }
}

#[test]

fn test_ray_sphere_intersection_json() {

    let input = r#"{
        "ray_origin": {"x": 0.0, "y": 0.0, "z": -5.0},
        "ray_direction": {"x": 0.0, "y": 0.0, "z": 1.0},
        "sphere_center": {"x": 0.0, "y": 0.0, "z": 0.0},
        "sphere_radius": 1.0
    }"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_ray_sphere_intersection_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        // Should hit at t=4
        assert!(parsed["ok"].is_object());

        let t = parsed["ok"]["t"]
            .as_f64()
            .unwrap();

        assert!((t - 4.0).abs() < 1e-10);
    }
}

#[test]

fn test_bezier_cubic_json() {

    let input = r#"{
        "p0": {"x": 0.0, "y": 0.0, "z": 0.0},
        "p1": {"x": 0.25, "y": 1.0, "z": 0.0},
        "p2": {"x": 0.75, "y": 1.0, "z": 0.0},
        "p3": {"x": 1.0, "y": 0.0, "z": 0.0},
        "t": 0.0
    }"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_bezier_cubic_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        // At t=0, should be at p0
        assert_eq!(parsed["ok"]["x"], 0.0);

        assert_eq!(parsed["ok"]["y"], 0.0);
    }
}

#[test]

fn test_angle_between_json() {

    let input = r#"{"v1": {"x": 1.0, "y": 0.0, "z": 0.0}, "v2": {"x": 0.0, "y": 1.0, "z": 0.0}}"#;

    let c_input = CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::json::rssn_num_graphics_angle_between_json(c_input.as_ptr());

        assert!(!result.is_null());

        let result_str = std::ffi::CStr::from_ptr(result).to_string_lossy();

        let parsed: serde_json::Value = serde_json::from_str(&result_str).unwrap();

        let angle = parsed["ok"]
            .as_f64()
            .unwrap();

        // Should be pi/2
        assert!((angle - std::f64::consts::PI / 2.0).abs() < 1e-10);
    }
}

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_dot_product_handle() {

    let result =
        rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_dot_product(
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        );

    assert_eq!(result, 0.0);
}

#[test]

fn test_magnitude_handle() {

    let result =
        rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_magnitude(
            3.0, 4.0, 0.0,
        );

    assert!((result - 5.0).abs() < 1e-10);
}

#[test]

fn test_cross_product_handle() {

    let mut out_x = 0.0;

    let mut out_y = 0.0;

    let mut out_z = 0.0;

    unsafe {

        let result = rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_cross_product(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            &mut out_x, &mut out_y, &mut out_z,
        );

        assert_eq!(result, 0);

        assert_eq!(out_z, 1.0);
    }
}

#[test]

fn test_normalize_handle() {

    let mut out_x = 0.0;

    let mut out_y = 0.0;

    let mut out_z = 0.0;

    unsafe {

        let result =
            rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_normalize(
                3.0, 4.0, 0.0, &mut out_x, &mut out_y, &mut out_z,
            );

        assert_eq!(result, 0);

        assert!((out_x - 0.6).abs() < 1e-10);

        assert!((out_y - 0.8).abs() < 1e-10);
    }
}

#[test]

fn test_angle_between_handle() {

    let result =
        rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_angle_between(
            1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
        );

    assert!((result - std::f64::consts::PI / 2.0).abs() < 1e-10);
}

#[test]

fn test_degrees_to_radians_handle() {

    let result = rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_degrees_to_radians(180.0);

    assert!((result - std::f64::consts::PI).abs() < 1e-10);
}

#[test]

fn test_radians_to_degrees_handle() {

    let result = rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_radians_to_degrees(std::f64::consts::PI);

    assert!((result - 180.0).abs() < 1e-10);
}

#[test]

fn test_ray_sphere_intersection_handle() {

    let result = rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_ray_sphere_intersection(
        0.0, 0.0, -5.0,  // ray origin
        0.0, 0.0, 1.0,   // ray direction
        0.0, 0.0, 0.0,   // sphere center
        1.0,             // sphere radius
    );

    assert!((result - 4.0).abs() < 1e-10);
}

#[test]

fn test_ray_sphere_intersection_miss_handle() {

    let result = rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_ray_sphere_intersection(
        0.0, 10.0, -5.0, // ray origin (misses)
        0.0, 0.0, 1.0,   // ray direction
        0.0, 0.0, 0.0,   // sphere center
        1.0,             // sphere radius
    );

    assert_eq!(result, -1.0);
}

#[test]

fn test_bezier_cubic_handle() {

    let mut out_x = 0.0;

    let mut out_y = 0.0;

    let mut out_z = 0.0;

    unsafe {

        let result =
            rssn::ffi_apis::numerical_computer_graphics_ffi::handle::rssn_num_graphics_bezier_cubic(
                0.0, 0.0, 0.0, // p0
                0.25, 1.0, 0.0, // p1
                0.75, 1.0, 0.0, // p2
                1.0, 0.0, 0.0, // p3
                0.0, // t
                &mut out_x, &mut out_y, &mut out_z,
            );

        assert_eq!(result, 0);

        assert_eq!(out_x, 0.0);

        assert_eq!(out_y, 0.0);
    }
}
