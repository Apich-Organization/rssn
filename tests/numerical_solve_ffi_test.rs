use rssn::ffi_apis::numerical_solve_ffi::*;
use rssn::prelude::numerical::*;
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_solve_handle_ffi() {

    unsafe {

        // 3x = 9 => x = 3
        let a = numerical_Matrix::new(
            1,
            1,
            vec![3.0],
        );

        let b = vec![9.0];

        let sol_ptr = handle::rssn_num_solve_linear_system_handle(
            &a,
            b.as_ptr(),
            b.len(),
        );

        assert!(!sol_ptr.is_null());

        assert!(handle::rssn_num_solve_is_unique(sol_ptr));

        let len = handle::rssn_num_solve_get_unique_solution_len(sol_ptr);

        assert_eq!(len, 1);

        let mut out = vec![0.0; len];

        handle::rssn_num_solve_get_unique_solution(
            sol_ptr,
            out.as_mut_ptr(),
        );

        assert!(
            (out[0] - 3.0).abs() < 1e-9
        );

        handle::rssn_num_solve_free_solution(sol_ptr);
    }
}

#[test]

fn test_solve_json_ffi() {

    unsafe {

        // x + y = 2
        // x - y = 0
        // Solution: x=1, y=1
        let a = numerical_Matrix::new(
            2,
            2,
            vec![1.0, 1.0, 1.0, -1.0],
        );

        let b = vec![2.0, 0.0];

        let json_input = format!(
            r#"{{"matrix": {}, "vector": {:?}}}"#,
            serde_json::to_string(&a)
                .unwrap(),
            b
        );

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_solve_linear_system_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        // Parse result
        // Expected: Ok(Unique([1.0, 1.0]))
        let v: serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        // Structure is FfiResult { ok: Some(LinearSolution), ... }
        // LinearSolution::Unique is { "Unique": [...] }
        let ok = &v["ok"];

        let sol = &ok["Unique"];

        assert_eq!(sol[0], 1.0);

        assert_eq!(sol[1], 1.0);

        let _ =
            CString::from_raw(res_ptr);
    }
}
