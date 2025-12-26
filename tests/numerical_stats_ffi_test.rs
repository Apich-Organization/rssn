use rssn::ffi_apis::common::{
    from_bincode_buffer,
    rssn_free_bincode_buffer,
    rssn_free_string,
    to_bincode_buffer,
};
use rssn::ffi_apis::ffi_api::FfiResult;
use rssn::ffi_apis::numerical_stats_ffi::{
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

fn test_stats_handle_ffi() {

    unsafe {

        let data = vec![
            1.0, 2.0, 3.0, 4.0, 5.0,
        ];

        // Mean
        let m = handle::rssn_num_stats_mean(
            data.as_ptr(),
            data.len(),
        );

        assert!((m - 3.0).abs() < 1e-10);

        // Variance
        let v = handle::rssn_num_stats_variance(
            data.as_ptr(),
            data.len(),
        );

        assert!(v > 0.0);

        // Std Dev
        let s = handle::rssn_num_stats_std_dev(
            data.as_ptr(),
            data.len(),
        );

        assert!(s > 0.0);

        // Geometric Mean
        let gm = handle::rssn_num_stats_geometric_mean(
            data.as_ptr(),
            data.len(),
        );

        assert!(gm > 0.0);

        // Harmonic Mean
        let hm = handle::rssn_num_stats_harmonic_mean(
            data.as_ptr(),
            data.len(),
        );

        assert!(hm > 0.0);

        // Range
        let r = handle::rssn_num_stats_range(
            data.as_ptr(),
            data.len(),
        );

        assert!((r - 4.0).abs() < 1e-10);

        // CV
        let cv = handle::rssn_num_stats_cv(
            data.as_ptr(),
            data.len(),
        );

        assert!(cv > 0.0);

        // Standard Error
        let se = handle::rssn_num_stats_standard_error(
            data.as_ptr(),
            data.len(),
        );

        assert!(se > 0.0);

        // Covariance
        let data2 = vec![
            2.0, 4.0, 6.0, 8.0, 10.0,
        ];

        let cov = handle::rssn_num_stats_covariance(
            data.as_ptr(),
            data.len(),
            data2.as_ptr(),
            data2.len(),
        );

        assert!(cov > 0.0);

        // Correlation
        let corr = handle::rssn_num_stats_correlation(
            data.as_ptr(),
            data.len(),
            data2.as_ptr(),
            data2.len(),
        );

        assert!((corr - 1.0).abs() < 1e-10);

        // Two-sample t-test
        let mut t = 0.0;

        let mut p = 0.0;

        handle::rssn_num_stats_two_sample_t_test(
            data.as_ptr(),
            data.len(),
            data.as_ptr(),
            data.len(),
            &mut t,
            &mut p,
        );

        assert!(t.abs() < 1e-10);

        // Welch t-test
        handle::rssn_num_stats_welch_t_test(
            data.as_ptr(),
            data.len(),
            data.as_ptr(),
            data.len(),
            &mut t,
            &mut p,
        );

        assert!(t.abs() < 1e-10);

        // Chi-squared test
        let obs = vec![10.0, 20.0, 30.0];

        let exp = vec![10.0, 20.0, 30.0];

        let mut chi = 0.0;

        handle::rssn_num_stats_chi_squared_test(
            obs.as_ptr(),
            exp.as_ptr(),
            obs.len(),
            &mut chi,
            &mut p,
        );

        assert!(chi.abs() < 1e-10);

        // Linear regression
        let x = vec![1.0, 2.0, 3.0, 4.0];

        let y = vec![3.0, 5.0, 7.0, 9.0]; // y = 2x + 1
        let mut slope = 0.0;

        let mut intercept = 0.0;

        handle::rssn_num_stats_linear_regression(
            x.as_ptr(),
            y.as_ptr(),
            x.len(),
            &mut slope,
            &mut intercept,
        );

        assert!(
            (slope - 2.0).abs() < 1e-10,
            "slope was {}",
            slope
        );

        assert!(
            (intercept - 1.0).abs() < 1e-10,
            "intercept was {}",
            intercept
        );
    }
}

#[derive(Serialize)]

struct DataInput {
    data: Vec<f64>,
}

#[derive(Serialize)]

struct TwoDataInput {
    data1: Vec<f64>,
    data2: Vec<f64>,
}

#[derive(Serialize)]

struct RegressionInput {
    x: Vec<f64>,
    y: Vec<f64>,
}

#[test]

fn test_stats_json_ffi() {

    unsafe {

        let input = DataInput {
            data: vec![
                1.0, 2.0, 3.0, 4.0, 5.0,
            ],
        };

        let json_str = serde_json::to_string(&input).unwrap();

        let c_json = CString::new(json_str).unwrap();

        // Mean
        let res_ptr = json::rssn_num_stats_mean_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - 3.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Variance
        let res_ptr = json::rssn_num_stats_variance_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            v["ok"]
                .as_f64()
                .unwrap()
                > 0.0
        );

        rssn_free_string(res_ptr);

        // Std Dev
        let res_ptr = json::rssn_num_stats_std_dev_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            v["ok"]
                .as_f64()
                .unwrap()
                > 0.0
        );

        rssn_free_string(res_ptr);

        // Z-scores
        let res_ptr = json::rssn_num_stats_z_scores_json(c_json.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        let z = v["ok"]
            .as_array()
            .unwrap();

        assert_eq!(z.len(), 5);

        rssn_free_string(res_ptr);

        // Covariance
        let two_input = TwoDataInput {
            data1: vec![
                1.0, 2.0, 3.0, 4.0, 5.0,
            ],
            data2: vec![
                2.0, 4.0, 6.0, 8.0, 10.0,
            ],
        };

        let json_str = serde_json::to_string(&two_input).unwrap();

        let c_json2 = CString::new(json_str).unwrap();

        let res_ptr = json::rssn_num_stats_covariance_json(c_json2.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            v["ok"]
                .as_f64()
                .unwrap()
                > 0.0
        );

        rssn_free_string(res_ptr);

        // Correlation
        let res_ptr = json::rssn_num_stats_correlation_json(c_json2.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            (v["ok"]
                .as_f64()
                .unwrap()
                - 1.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);

        // Linear regression
        let reg_input = RegressionInput {
            x: vec![1.0, 2.0, 3.0, 4.0],
            y: vec![3.0, 5.0, 7.0, 9.0],
        };

        let json_str = serde_json::to_string(&reg_input).unwrap();

        let c_json3 = CString::new(json_str).unwrap();

        let res_ptr = json::rssn_num_stats_linear_regression_json(c_json3.as_ptr());

        let res_str = CStr::from_ptr(res_ptr)
            .to_str()
            .unwrap();

        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();

        assert!(
            (v["ok"]["slope"]
                .as_f64()
                .unwrap()
                - 2.0)
                .abs()
                < 1e-10
        );

        assert!(
            (v["ok"]["intercept"]
                .as_f64()
                .unwrap()
                - 1.0)
                .abs()
                < 1e-10
        );

        rssn_free_string(res_ptr);
    }
}

#[test]

fn test_stats_bincode_ffi() {

    unsafe {

        #[derive(Serialize)]

        struct DataInputB {
            data: Vec<f64>,
        }

        let input = DataInputB {
            data: vec![
                1.0, 2.0, 3.0, 4.0, 5.0,
            ],
        };

        let buffer = to_bincode_buffer(&input);

        // Mean
        let res_buffer = bincode_api::rssn_num_stats_mean_bincode(buffer);

        let res: FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert!((res.ok.unwrap() - 3.0).abs() < 1e-10);

        rssn_free_bincode_buffer(res_buffer);

        // Variance
        let res_buffer = bincode_api::rssn_num_stats_variance_bincode(buffer);

        let res: FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert!(res.ok.unwrap() > 0.0);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer);

        // Two data input
        #[derive(Serialize)]

        struct TwoDataB {
            data1: Vec<f64>,
            data2: Vec<f64>,
        }

        let two_input = TwoDataB {
            data1: vec![
                1.0, 2.0, 3.0, 4.0, 5.0,
            ],
            data2: vec![
                2.0, 4.0, 6.0, 8.0, 10.0,
            ],
        };

        let buffer2 = to_bincode_buffer(&two_input);

        // Correlation
        let res_buffer = bincode_api::rssn_num_stats_correlation_bincode(buffer2);

        let res: FfiResult<f64, String> = from_bincode_buffer(&res_buffer).unwrap();

        assert!((res.ok.unwrap() - 1.0).abs() < 1e-10);

        rssn_free_bincode_buffer(res_buffer);

        rssn_free_bincode_buffer(buffer2);
    }
}
