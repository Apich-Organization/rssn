use assert_approx_eq::assert_approx_eq;
use rssn::ffi_apis::common::{
    rssn_free_bincode_buffer,
    rssn_free_string,
};
use rssn::ffi_apis::numerical_differential_geometry_ffi::{
    bincode_api,
    handle,
    json,
};
use rssn::symbolic::coordinates::CoordinateSystem;
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_dg_handle_ffi() {

    unsafe {

        let point = vec![1.0, 2.0, 3.0];

        let mut result = 0.0;

        let status = handle::rssn_num_dg_ricci_scalar(
            CoordinateSystem::Cartesian,
            point.as_ptr(),
            3,
            &mut result,
        );

        assert_eq!(status, 0);

        assert_approx_eq!(result, 0.0, 1e-9);

        let g_ptr = handle::rssn_num_dg_metric_tensor(
            CoordinateSystem::Cartesian,
            point.as_ptr(),
            3,
        );

        assert!(!g_ptr.is_null());

        let g = &*g_ptr;

        assert_eq!(g.data()[0], 1.0);

        let _ = Box::from_raw(g_ptr);
    }
}

#[test]

fn test_dg_json_ffi() {

    unsafe {

        let _point = vec![
            2.0,
            0.0,
            std::f64::consts::PI / 2.0,
        ];

        let json_input = format!(
            r#"{{"system": "Spherical", "point": [2.0, 0.0, {}]}}"#,
            std::f64::consts::PI / 2.0
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = json::rssn_num_dg_ricci_scalar_json(c_json.as_ptr());

        assert!(!res_ptr.is_null());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert_approx_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            0.0,
            1e-9
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_dg_bincode_ffi() {

    unsafe {

        use rssn::ffi_apis::common::{
            from_bincode_buffer,
            to_bincode_buffer,
        };
        use serde::{
            Deserialize,
            Serialize,
        };

        #[derive(Serialize)]

        struct DgInput {
            system: CoordinateSystem,
            point: Vec<f64>,
        }

        let input = DgInput {
            system: CoordinateSystem::Spherical,
            point: vec![
                1.0,
                0.0,
                std::f64::consts::PI / 2.0,
            ],
        };

        let buffer = to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_dg_ricci_scalar_bincode(buffer);

        assert!(!res_buffer.is_null());

        #[derive(Deserialize)]

        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }

        let res: FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert_approx_eq!(
            res.ok.unwrap(),
            0.0,
            1e-9
        );

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);
    }
}
