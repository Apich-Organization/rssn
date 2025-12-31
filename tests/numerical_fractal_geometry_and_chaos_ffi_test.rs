//! FFI Tests for numerical fractal geometry and chaos module.

use std::ffi::CString;

// ============================================================================
// JSON FFI Tests
// ============================================================================

#[test]

fn test_mandelbrot_escape_time_json() {

    let input = r#"{"c_real": 0.0, "c_imag": 0.0, "max_iter": 100}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_mandelbrot_escape_time_json(c_input.as_ptr());

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

        assert_eq!(parsed["ok"], 100);
    }
}

#[test]

fn test_julia_escape_time_json() {

    let input = r#"{"z_real": 0.0, "z_imag": 0.0, "c_real": 0.0, "c_imag": 0.0, "max_iter": 100}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_julia_escape_time_json(c_input.as_ptr());

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

        assert_eq!(parsed["ok"], 100);
    }
}

#[test]

fn test_lorenz_attractor_json() {

    let input = r#"{"start_point": [1.0, 1.0, 1.0], "dt": 0.01, "num_steps": 10}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_lorenz_attractor_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            10
        );
    }
}

#[test]

fn test_henon_map_json() {

    let input = r#"{"start_point": [0.0, 0.0], "num_steps": 10, "a": 1.4, "b": 0.3}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_henon_map_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            10
        );
    }
}

#[test]

fn test_logistic_map_json() {

    let input = r#"{"x0": 0.5, "r": 3.5, "num_steps": 10}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_logistic_map_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            11
        ); // x0 + 10 iterations
    }
}

#[test]

fn test_lyapunov_logistic_json() {

    let input = r#"{"r": 4.0, "x0": 0.5, "transient": 100, "num_iterations": 500}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_lyapunov_logistic_json(c_input.as_ptr());

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

        // r=4 should give positive Lyapunov exponent
        let lyap = parsed["ok"]
            .as_f64()
            .unwrap();

        assert!(lyap > 0.0);
    }
}

#[test]

fn test_mandelbrot_set_json() {

    let input = r#"{"width": 5, "height": 5, "x_range": [-2.0, 1.0], "y_range": [-1.5, 1.5], "max_iter": 20}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_mandelbrot_set_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            5
        );
    }
}

#[test]

fn test_julia_set_json() {

    let input = r#"{"width": 5, "height": 5, "x_range": [-2.0, 2.0], "y_range": [-2.0, 2.0], "c": [-0.4, 0.6], "max_iter": 20}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_julia_set_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            5
        );
    }
}

#[test]

fn test_rossler_attractor_json() {

    let input = r#"{"start_point": [1.0, 1.0, 1.0], "dt": 0.01, "num_steps": 10, "a": 0.2, "b": 0.2, "c": 5.7}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_rossler_attractor_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            10
        );
    }
}

#[test]

fn test_tinkerbell_map_json() {

    let input = r#"{"start_point": [-0.72, -0.64], "num_steps": 10, "a": 0.9, "b": -0.6013, "c": 2.0, "d": 0.5}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_tinkerbell_map_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            10
        );
    }
}

#[test]

fn test_bifurcation_json() {

    let input = r#"{"r_range": [2.5, 4.0], "num_r_values": 5, "transient": 50, "num_points": 3, "x0": 0.5}"#;

    let c_input =
        CString::new(input).unwrap();

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::json::rssn_num_fractal_bifurcation_json(c_input.as_ptr());

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

        assert!(
            parsed["ok"].is_array()
        );

        assert_eq!(
            parsed["ok"]
                .as_array()
                .unwrap()
                .len(),
            15
        ); // 5 * 3
    }
}

// ============================================================================
// Handle-based FFI Tests
// ============================================================================

#[test]

fn test_mandelbrot_escape_time_handle()
{

    let escape = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_mandelbrot_escape_time(
        0.0,
        0.0,
        100,
    );

    assert_eq!(escape, 100);
}

#[test]

fn test_julia_escape_time_handle() {

    let escape = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_julia_escape_time(
        0.0,
        0.0,
        0.0,
        0.0,
        100,
    );

    assert_eq!(escape, 100);
}

#[test]

fn test_lyapunov_logistic_handle() {

    let lyap = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_lyapunov_logistic(
        4.0,
        0.5,
        100,
        500,
    );

    assert!(lyap > 0.0);
}

#[test]

fn test_lorenz_attractor_handle() {

    let mut output = vec![0.0f64; 30]; // 10 points * 3

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_lorenz_attractor(
            1.0,
            1.0,
            1.0,
            0.01,
            10,
            output.as_mut_ptr(),
        );

        assert_eq!(result, 0);
    }
}

#[test]

fn test_henon_map_handle() {

    let mut output = vec![0.0f64; 20]; // 10 points * 2

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_henon_map(
            0.0,
            0.0,
            10,
            1.4,
            0.3,
            output.as_mut_ptr(),
        );

        assert_eq!(result, 0);
    }
}

#[test]

fn test_logistic_map_handle() {

    let mut output = vec![0.0f64; 11]; // x0 + 10 iterations

    unsafe {

        let result = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_logistic_map(
            0.5,
            3.5,
            10,
            output.as_mut_ptr(),
        );

        assert_eq!(result, 0);

        assert!(
            (output[0] - 0.5).abs()
                < 1e-10
        );
    }
}

#[test]

fn test_box_counting_dim_handle() {

    let points: Vec<f64> = (0 .. 100)
        .flat_map(|i| {

            vec![
                i as f64 / 100.0,
                0.0,
            ]
        })
        .collect();

    unsafe {

        let dim = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_box_counting_dim(
            points.as_ptr(),
            50,
            8,
        );

        assert!(
            dim >= 0.5 && dim <= 1.5
        );
    }
}

#[test]

fn test_correlation_dim_handle() {

    let points: Vec<f64> = (0 .. 100)
        .flat_map(|i| {

            vec![
                i as f64 / 100.0,
                0.0,
            ]
        })
        .collect();

    unsafe {

        let dim = rssn::ffi_apis::numerical_fractal_geometry_and_chaos_ffi::handle::rssn_num_fractal_correlation_dim(
            points.as_ptr(),
            50,
            8,
        );

        assert!(dim.is_finite());
    }
}
