use std::ffi::CStr;
use std::ffi::CString;

use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::numerical_signal_ffi::bincode_api;
use rssn::ffi_apis::numerical_signal_ffi::handle;
use rssn::ffi_apis::numerical_signal_ffi::json;
use rssn::numerical::matrix::Matrix;
use rustfft::num_complex::Complex;

#[test]

fn test_numerical_signal_handle_ffi() {

    unsafe {

        let re =
            vec![1.0, 1.0, 1.0, 1.0];

        let im =
            vec![0.0, 0.0, 0.0, 0.0];

        // FFT
        let matrix_ptr =
            handle::rssn_num_signal_fft(
                re.as_ptr(),
                im.as_ptr(),
                4,
            );

        assert!(!matrix_ptr.is_null());

        let matrix = &*matrix_ptr;

        assert_approx_eq!(
            matrix.data()[0],
            4.0,
            1e-9
        ); // Real part of DC
        assert_approx_eq!(
            matrix.data()[1],
            0.0,
            1e-9
        ); // Imag part of DC
        let _ =
            Box::from_raw(matrix_ptr);

        // Convolve
        let a = vec![1.0, 2.0, 3.0];

        let v = vec![0.0, 1.0, 0.5];

        let res_ptr = handle::rssn_num_signal_convolve(
            a.as_ptr(),
            3,
            v.as_ptr(),
            3,
        );

        assert!(!res_ptr.is_null());

        let res_matrix = &*res_ptr;

        assert_eq!(
            res_matrix.data(),
            &vec![
                0.0, 1.0, 2.5, 4.0, 1.5
            ]
        );

        let _ = Box::from_raw(res_ptr);
    }
}

#[test]

fn test_numerical_signal_json_ffi() {

    unsafe {

        let input_json = r#"{"a": [1.0, 2.0], "v": [1.0, 0.5]}"#;

        let c_json =
            CString::new(input_json)
                .unwrap();

        let res_ptr = json::rssn_num_signal_convolve_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v: serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        let results = v["ok"]
            .as_array()
            .unwrap();

        assert_approx_eq!(
            results[0]
                .as_f64()
                .unwrap(),
            1.0,
            1e-9
        );

        assert_approx_eq!(
            results[1]
                .as_f64()
                .unwrap(),
            2.5,
            1e-9
        );

        assert_approx_eq!(
            results[2]
                .as_f64()
                .unwrap(),
            1.0,
            1e-9
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_numerical_signal_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::from_bincode_buffer;
        use rssn::ffi_apis::common::to_bincode_buffer;
        use serde::Deserialize;
        use serde::Serialize;

        #[derive(Serialize)]

        struct ConvolveInput {
            a: Vec<f64>,
            v: Vec<f64>,
        }

        let input = ConvolveInput {
            a: vec![1.0, 2.0],
            v: vec![1.0, 0.5],
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_signal_convolve_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<
            Vec<f64>,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        let ok_res = res
            .ok
            .as_ref()
            .unwrap();

        assert_approx_eq!(
            ok_res[0],
            1.0,
            1e-9
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
