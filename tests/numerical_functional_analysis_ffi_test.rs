use std::ffi::CStr;
use std::ffi::CString;

use rssn::ffi_apis::common::from_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_bincode_buffer;
use rssn::ffi_apis::common::rssn_free_string;
use rssn::ffi_apis::common::to_bincode_buffer;
use rssn::ffi_apis::common::BincodeBuffer;
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_functional_analysis_ffi::bincode_api;
use rssn::ffi_apis::numerical_functional_analysis_ffi::handle;
use rssn::ffi_apis::numerical_functional_analysis_ffi::json;
use serde::Deserialize;
use serde::Serialize;

#[test]

fn test_fa_handle_ffi() {

    unsafe {

        let x = [0.0, 1.0, 2.0];

        let y = [1.0, 1.0, 1.0];

        let res =
            handle::rssn_num_fa_l2_norm(
                x.as_ptr(),
                y.as_ptr(),
                3,
            );

        assert!(
            (res - 2.0f64.sqrt()).abs()
                < 1e-9
        );
    }
}

#[test]

fn test_fa_json_ffi() {

    unsafe {

        let json_input = r#"{
            "points": [[0.0, 1.0], [1.0, 1.0], [2.0, 1.0]]
        }"#;

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_num_fa_l2_norm_json(c_json.as_ptr());

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

        let res = v["ok"]
            .as_f64()
            .unwrap();

        assert!(
            (res - 2.0f64.sqrt()).abs()
                < 1e-9
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_fa_bincode_ffi() {

    unsafe {

        #[derive(Serialize)]

        struct PointsInput {
            points: Vec<(f64, f64)>,
        }

        let input = PointsInput {
            points: vec![
                (0.0, 1.0),
                (1.0, 1.0),
                (2.0, 1.0),
            ],
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_fa_l2_norm_bincode(buffer);

        assert!(!res_buffer.is_null());

        let res: FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert!(
            (res.ok.unwrap()
                - 2.0f64.sqrt())
            .abs()
                < 1e-9
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
