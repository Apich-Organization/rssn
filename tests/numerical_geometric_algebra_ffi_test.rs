use rssn::ffi_apis::numerical_geometric_algebra_ffi::handle::*;
use rssn::ffi_apis::numerical_geometric_algebra_ffi::json::*;
use rssn::prelude::numerical::numerical_Multivector3D as Multivector3D;
use std::ffi::{CStr, CString};

#[test]

fn test_ga_handle_ffi() {

    unsafe {

        let mv1 = rssn_num_ga_create(
            1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0,
        );

        let mv2 = rssn_num_ga_create(
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        );

        let mv_sum = rssn_num_ga_add(mv1, mv2);

        let mut s = 0.0;

        let mut v1 = 0.0;

        rssn_num_ga_get_components(
            mv_sum,
            &mut s,
            &mut v1,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );

        assert_eq!(s, 1.0);

        assert_eq!(v1, 3.0);

        let mv_prod = rssn_num_ga_mul(mv1, mv2);

        let mut b12 = 0.0;

        rssn_num_ga_get_components(
            mv_prod,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            &mut b12,
            std::ptr::null_mut(),
            std::ptr::null_mut(),
            std::ptr::null_mut(),
        );

        // e1 * e1 = 1, e1*e2 = e12
        // mv1 = 1 + 2e1 + 3e2 + 4e3
        // mv2 = e1
        // mv1 * mv2 = e1 + 2(e1*e1) + 3(e2*e1) + 4(e3*e1) = e1 + 2 - 3e12 + 4e31
        assert_eq!(b12, -3.0);

        rssn_num_ga_free(mv1);

        rssn_num_ga_free(mv2);

        rssn_num_ga_free(mv_sum);

        rssn_num_ga_free(mv_prod);
    }
}

#[test]

fn test_ga_json_ffi() {

    unsafe {

        let mv1 = Multivector3D::new(
            1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        );

        let mv2 = Multivector3D::new(
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        );

        let json_input = format!(
            "{{\"mv1\": {}, \"mv2\": {}}}",
            serde_json::to_string(&mv1).unwrap(),
            serde_json::to_string(&mv2).unwrap()
        );

        let c_json = CString::new(json_input).unwrap();

        let res_ptr = rssn_num_ga_add_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        // FfiResult { ok: Some(Multivector3D), err: None }
        let res: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert_eq!(res["ok"]["s"], 1.0);

        assert_eq!(res["ok"]["v1"], 3.0);

        // Use CString::from_raw to free the pointer allocated by CString::into_raw in the FFI
        let _ = CString::from_raw(res_ptr);
    }
}
