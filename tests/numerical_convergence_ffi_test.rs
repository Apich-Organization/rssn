use rssn::ffi_apis::numerical_convergence_ffi::*;
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_convergence_handle_ffi() {

    unsafe {

        // Aitken
        let seq = vec![
            1.5,
            1.25,
            1.125,
            1.0625,
        ]; // 1 + 0.5^i for i=1..4 (converges to 1)
        let ptr = handle::rssn_convergence_aitken(
            seq.as_ptr(),
            seq.len(),
        );

        assert!(!ptr.is_null());

        let len = handle::rssn_convergence_get_vec_len(ptr);

        let mut out = vec![0.0; len];

        handle::rssn_convergence_get_vec_data(
            ptr,
            out.as_mut_ptr(),
        );

        handle::rssn_convergence_free_vec(ptr);

        assert!(len > 0);

        assert!(
            (out.last().unwrap() - 1.0)
                .abs()
                < 1e-10
        );

        // Richardson
        let rich_seq =
            vec![1.08, 1.02, 1.005]; // Dummy sequence converging to 1
        let r_ptr = handle::rssn_convergence_richardson(
            rich_seq.as_ptr(),
            rich_seq.len(),
        );

        handle::rssn_convergence_free_vec(r_ptr);
    }
}

#[test]

fn test_convergence_json_ffi() {

    unsafe {

        // Wynn Epsilon
        let seq = vec![
            1.0,
            0.66666,
            0.86666,
            0.7238,
        ]; // Alternating series partial sums
        let json_input = format!(
            r#"{{"sequence": {:?}}}"#,
            seq
        );

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_convergence_wynn_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        // Check valid JSON
        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert!(v
            .get("ok")
            .is_some());

        // Free
        let _ =
            CString::from_raw(res_ptr);
    }
}
