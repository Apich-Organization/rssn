use rssn::ffi_apis::numerical_ode_ffi::*;
use rssn::symbolic::core::Expr;
use std::ffi::{CStr, CString};
use rssn::ffi_apis::common::{rssn_free_string, rssn_free_bincode_buffer};
use assert_approx_eq::assert_approx_eq;

#[test]
fn test_numerical_ode_handle_ffi() {
    unsafe {
        let y0_expr = Expr::new_variable("y0");
        let funcs = vec![&y0_expr as *const Expr];
        let y_init = vec![1.0];
        
        let matrix_ptr = handle::rssn_num_ode_solve(
            funcs.as_ptr(),
            1,
            y_init.as_ptr(),
            1,
            0.0,
            1.0,
            100,
            2, // RungeKutta4
        );
        
        assert!(!matrix_ptr.is_null());
        let matrix = &*matrix_ptr;
        // Last element in matrix data should be approx e^1
        let last_val = *matrix.data().last().unwrap();
        assert_approx_eq!(last_val, 2.71828, 1e-5);
        
        let _ = Box::from_raw(matrix_ptr);
    }
}

#[test]
fn test_numerical_ode_json_ffi() {
    unsafe {
        let f = Expr::new_variable("y0");
        let f_json = serde_json::to_string(&f).unwrap();
        let json_input = format!(
            r#"{{"funcs": [{}], "y0": [1.0], "x_range": [0.0, 1.0], "num_steps": 100, "method": "RungeKutta4"}}"#,
            f_json
        );
        let c_json = CString::new(json_input).unwrap();
        
        let res_ptr = json::rssn_num_ode_solve_json(c_json.as_ptr());
        assert!(!res_ptr.is_null());
        
        let res_str = CStr::from_ptr(res_ptr).to_str().unwrap();
        let v: serde_json::Value = serde_json::from_str(res_str).unwrap();
        
        let results = v["ok"].as_array().unwrap();
        let last_y = results.last().unwrap().as_array().unwrap()[0].as_f64().unwrap();
        assert_approx_eq!(last_y, 2.71828, 1e-5);
        
        rssn_free_string(res_ptr);
    }
}

#[test]
fn test_numerical_ode_bincode_ffi() {
    unsafe {
        use rssn::ffi_apis::common::{to_bincode_buffer, from_bincode_buffer};
        use serde::{Serialize, Deserialize};
        use rssn::numerical::ode::OdeSolverMethod;

        #[derive(Serialize)]
        struct OdeInput {
            funcs: Vec<Expr>,
            y0: Vec<f64>,
            x_range: (f64, f64),
            num_steps: usize,
            method: OdeSolverMethod,
        }

        let input = OdeInput {
            funcs: vec![Expr::new_variable("y0")],
            y0: vec![1.0],
            x_range: (0.0, 1.0),
            num_steps: 100,
            method: OdeSolverMethod::RungeKutta4,
        };
        
        let buffer = to_bincode_buffer(&input);
        let res_buffer = bincode_api::rssn_num_ode_solve_bincode(buffer);
        assert!(!res_buffer.is_null());
        
        #[derive(Deserialize)]
        struct FfiResult<T, E> {
            ok: Option<T>,
            #[allow(dead_code)]
            err: Option<E>,
        }
        let res: FfiResult<Vec<Vec<f64>>, String> = from_bincode_buffer(&res_buffer).unwrap();
        let ok_res = res.ok.as_ref().unwrap();
        assert_approx_eq!(ok_res.last().unwrap()[0], 2.71828, 1e-5);
        
        rssn_free_bincode_buffer(res_buffer);
        rssn_free_bincode_buffer(buffer);
    }
}
