use rssn::ffi_apis::common::{
    from_bincode_buffer,
    rssn_free_bincode_buffer,
    rssn_free_string,
    to_bincode_buffer,
    BincodeBuffer,
};
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_combinatorics_ffi::{
    bincode_api,
    handle,
    json,
};
use serde::{
    Deserialize,
    Serialize,
};
use std::ffi::{
    CStr,
    CString,
};

#[test]

fn test_comb_handle_ffi() {

    unsafe {

        let mut res = 0.0;

        // Factorial
        handle::rssn_num_comb_factorial(
            5,
            &mut res,
        );

        assert_eq!(res, 120.0);

        // Permutations
        handle::rssn_num_comb_permutations(5, 2, &mut res);

        assert_eq!(res, 20.0);

        // Combinations
        handle::rssn_num_comb_combinations(5, 2, &mut res);

        assert_eq!(res, 10.0);

        // Recurrence (Fibonacci)
        let coeffs = [1.0, 1.0];

        let initial = [0.0, 1.0];

        handle::rssn_num_comb_solve_recurrence(
            coeffs.as_ptr(),
            coeffs.len(),
            initial.as_ptr(),
            initial.len(),
            5,
            &mut res,
        );

        assert_eq!(res, 5.0);
    }
}

#[test]

fn test_comb_json_ffi() {

    unsafe {

        // Factorial
        let json_input = r#"{"n": 5}"#;

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_num_comb_factorial_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            120.0
        );

        rssn_free_string(res_ptr);

        // Recurrence
        let json_input = r#"{
            "coeffs": [1.0, 1.0],
            "initial_conditions": [0.0, 1.0],
            "target_n": 5
        }"#;

        let c_json =
            CString::new(json_input)
                .unwrap();

        let res_ptr = json::rssn_num_comb_solve_recurrence_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            5.0
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_comb_bincode_ffi() {

    unsafe {

        #[derive(Serialize)]

        struct NInput {
            n : u64,
        }

        let input = NInput {
            n : 5,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_comb_factorial_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert_eq!(
            res.ok.unwrap(),
            120.0
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );

        #[derive(Serialize)]

        struct RecurrenceInput {
            coeffs : Vec<f64>,
            initial_conditions :
                Vec<f64>,
            target_n : usize,
        }

        let input = RecurrenceInput {
            coeffs : vec![1.0, 1.0],
            initial_conditions : vec![
                0.0, 1.0,
            ],
            target_n : 5,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_comb_solve_recurrence_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert_eq!(
            res.ok.unwrap(),
            5.0
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}

#[test]

fn test_comb_handle_others() {

    unsafe {

        let mut res = 0.0;

        // Stirling second: S(3, 2) = 3
        handle::rssn_num_comb_stirling_second(3, 2, &mut res);

        assert_eq!(res, 3.0);

        // Bell: B(3) = 5
        handle::rssn_num_comb_bell(
            3,
            &mut res,
        );

        assert_eq!(res, 5.0);

        // Catalan: C_3 = 5
        handle::rssn_num_comb_catalan(
            3,
            &mut res,
        );

        assert_eq!(res, 5.0);

        // Rising factorial: 2^(3) = 24
        handle::rssn_num_comb_rising_factorial(2.0, 3, &mut res);

        assert_eq!(res, 24.0);

        // Falling factorial: 4_2 = 12
        handle::rssn_num_comb_falling_factorial(4.0, 2, &mut res);

        assert_eq!(res, 12.0);
    }
}

#[test]

fn test_comb_json_others() {

    unsafe {

        #[derive(Serialize)]

        struct NInput {
            n : u64,
        }

        #[derive(Serialize)]

        struct NKInput {
            n : u64,
            k : u64,
        }

        #[derive(Serialize)]

        struct XNInput {
            x : f64,
            n : u64,
        }

        // Stirling S(3, 2) = 3
        let input = NKInput {
            n : 3,
            k : 2,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_comb_stirling_second_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            3.0
        );

        rssn_free_string(res_ptr);

        // Rising factorial 2^(3) = 24
        let input = XNInput {
            x : 2.0,
            n : 3,
        };

        let json_str =
            serde_json::to_string(
                &input,
            )
            .unwrap();

        let c_json =
            CString::new(json_str)
                .unwrap();

        let res_ptr = json::rssn_num_comb_rising_factorial_json(c_json.as_ptr());

        let res_str =
            CStr::from_ptr(res_ptr)
                .to_str()
                .unwrap();

        let v : serde_json::Value =
            serde_json::from_str(
                res_str,
            )
            .unwrap();

        assert_eq!(
            v["ok"]
                .as_f64()
                .unwrap(),
            24.0
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_comb_bincode_others() {

    unsafe {

        #[derive(Serialize)]

        struct NInput {
            n : u64,
        }

        #[derive(Serialize)]

        struct NKInput {
            n : u64,
            k : u64,
        }

        // Bell B(3) = 5
        let input = NInput {
            n : 3,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_comb_bell_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert_eq!(
            res.ok.unwrap(),
            5.0
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );

        // Catalan C_3 = 5
        let input = NInput {
            n : 3,
        };

        let buffer =
            to_bincode_buffer(&input);

        let res_buffer = bincode_api::rssn_num_comb_catalan_bincode(buffer);

        let res : FfiResult<
            f64,
            String,
        > = from_bincode_buffer(
            &res_buffer,
        )
        .unwrap();

        assert_eq!(
            res.ok.unwrap(),
            5.0
        );

        rssn_free_bincode_buffer(
            res_buffer,
        );

        rssn_free_bincode_buffer(
            buffer,
        );
    }
}
