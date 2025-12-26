use assert_approx_eq::assert_approx_eq;
use num_complex::Complex;
use rssn::ffi_apis::common::{rssn_free_bincode_buffer, rssn_free_string};
use rssn::ffi_apis::numerical_transforms_ffi::{bincode_api, handle, json};
use std::ffi::{CStr, CString};

#[test]

fn test_numerical_transforms_handle_ffi() {

    unsafe {

        let mut re = vec![1.0, 1.0, 0.0, 0.0];

        let mut im = vec![0.0, 0.0, 0.0, 0.0];

        // FFT Inplace
        let status = handle::rssn_num_fft_inplace(re.as_mut_ptr(), im.as_mut_ptr(), 4);

        assert_eq!(status, 0);

        assert_approx_eq!(re[0], 2.0f64, 1e-9);

        assert_approx_eq!(re[1], 1.0f64, 1e-9);

        assert_approx_eq!(im[1], -1.0f64, 1e-9);

        // IFFT Inplace
        let status2 = handle::rssn_num_ifft_inplace(re.as_mut_ptr(), im.as_mut_ptr(), 4);

        assert_eq!(status2, 0);

        assert_approx_eq!(re[0], 1.0f64, 1e-9);

        assert_approx_eq!(re[1], 1.0f64, 1e-9);

        assert_approx_eq!(re[2], 0.0f64, 1e-9);

        assert_approx_eq!(re[3], 0.0f64, 1e-9);
    }
}

#[test]

fn test_numerical_transforms_json_ffi() {

    unsafe {

        let input_json = r#"{"data": [[1.0, 0.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0]]}"#;

        let c_json = CString::new(input_json).unwrap();

        let res_ptr = json::rssn_num_fft_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        println!("JSON Response: {}", res_str);

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        let results = v["ok"]
            .as_array()
            .unwrap();

        assert_approx_eq!(
            results[0][0]
                .as_f64()
                .unwrap(),
            2.0f64,
            1e-9
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_transforms_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer};
        use serde::{Deserialize, Serialize};

        #[derive(Serialize)]

        struct TransformInput {
            data: Vec<Complex<f64>>,
        }

        let input = TransformInput {
            data: vec![
                Complex::new(1.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_fft_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<Vec<Complex<f64>>, String> = from_bincode_buffer(&res_buffer).unwrap();

        let ok_res = res
            .ok
            .as_ref()
            .unwrap();

        assert_approx_eq!(ok_res[0].re, 2.0, 1e-9);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
