use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_interpolate_ffi::*;
use rssn::numerical::polynomial::Polynomial;
use std::ffi::{CStr, CString};

#[test]
fn test_numerical_interpolate_handle_ffi() {
    unsafe {
        // Lagrange
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, 4.0];
        let poly_ptr = handle::rssn_num_lagrange_interpolation(x.as_ptr(), y.as_ptr(), 3);
        assert!(!poly_ptr.is_null());
        let poly = &*poly_ptr;
        assert_approx_eq!(poly.eval(1.5), 2.25, 1e-9);
        let _ = Box::from_raw(poly_ptr);

        // Cubic Spline
        let x2 = vec![0.0, 1.0, 2.0];
        let y2 = vec![0.0, 1.0, 0.0];
        let handle = handle::rssn_num_cubic_spline_interpolation(x2.as_ptr(), y2.as_ptr(), 3);
        assert!(!handle.is_null());
        let val = handle::rssn_num_cubic_spline_evaluate(handle, 0.5);
        assert_approx_eq!(val, 0.6875, 1e-9);
        handle::rssn_num_cubic_spline_free(handle);

        // Bezier
        let cp = vec![
            0.0, 0.0, 1.0, 2.0, 2.0, 0.0,
        ];
        let mut out = vec![0.0, 0.0];
        let status = handle::rssn_num_bezier_curve(cp.as_ptr(), 3, 2, 0.5, out.as_mut_ptr());
        assert_eq!(status, 0);
        assert_approx_eq!(out[0], 1.0, 1e-9);
        assert_approx_eq!(out[1], 1.0, 1e-9);
    }
}

#[test]
fn test_numerical_interpolate_json_ffi() {
    unsafe {
        // Lagrange
        let json_input = r#"{"points": [[0.0, 0.0], [1.0, 1.0], [2.0, 4.0]]}"#;
        let c_json = CString::new(json_input).unwrap();
        let res_ptr = json::rssn_num_lagrange_interpolation_json(c_json.as_ptr());
        assert!(!res_ptr.is_null());
        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();
        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();
        assert!(v["ok"].is_object());
        rssn_free_string(res_ptr);

        // Cubic Spline
        let json_input2 = r#"{"points": [[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], "x_eval": 0.5}"#;
        let c_json2 = CString::new(json_input2).unwrap();
        let res_ptr2 = json::rssn_num_cubic_spline_interpolation_json(c_json2.as_ptr());
        assert_approx_eq!(
            serde_json::from_str::<serde_json::Value>(CStr::from_ptr(res_ptr2).to_str().unwrap())
                .unwrap()["ok"]
                .as_f64()
                .unwrap(),
            0.6875,
            1e-9
        );
        rssn_free_string(res_ptr2);
    }
}

#[test]
fn test_numerical_interpolate_bincode_ffi() {
    unsafe {
        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]
        struct LagrangeInput {
            points: Vec<(f64, f64)>,
        }

        let input = LagrangeInput {
            points: vec![
                (0.0, 0.0),
                (1.0, 1.0),
                (2.0, 4.0),
            ],
        };
        let buffer = to_bincode_buffer(&input);
        let res_buffer = bincode_api::rssn_num_lagrange_interpolation_bincode(buffer);
        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]
        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }
        let res: FfiResult<Polynomial, String> = from_bincode_buffer(&res_buffer).unwrap();
        assert!(res.ok.is_some());
        assert_approx_eq!(res.ok.unwrap().eval(1.5), 2.25, 1e-9);

        rssn_free_bincode_buffer(res_buffer);
        rssn_free_bincode_buffer(buffer);
    }
}
