use rssn::symbolic::cad::*;
use rssn::symbolic::core::{Expr, Monomial, SparsePolynomial};
use std::collections::BTreeMap;

#[test]
fn test_cad_simple_1d() {
    // p(x) = x^2 - 1
    let mut terms = BTreeMap::new();
    let mut vars_map = BTreeMap::new();
    vars_map.insert("x".to_string(), 2);
    terms.insert(Monomial(vars_map), Expr::Constant(1.0));
    let vars_map_0 = BTreeMap::new();
    terms.insert(Monomial(vars_map_0), Expr::Constant(-1.0));
    let p = SparsePolynomial { terms };

    let result = cad(&[p], &["x"]).unwrap();
    // Roots are -1 and 1.
    // Intervals: (-inf, -1), {-1}, (-1, 1), {1}, (1, inf)
    assert_eq!(result.cells.len(), 5);

    // Check sample points
    assert!(result.cells[0].sample_point[0] < -1.0);
    assert!((result.cells[1].sample_point[0] - (-1.0)).abs() < 1e-7);
    assert!(result.cells[2].sample_point[0] > -1.0 && result.cells[2].sample_point[0] < 1.0);
    assert!((result.cells[3].sample_point[0] - 1.0).abs() < 1e-7);
    assert!(result.cells[4].sample_point[0] > 1.0);
}

#[test]
fn test_cad_simple_2d() {
    // p(x, y) = x^2 + y^2 - 1 (Unit circle)
    // vars = [x, y]
    let mut terms = BTreeMap::new();

    let mut vars_x2 = BTreeMap::new();
    vars_x2.insert("x".to_string(), 2);
    terms.insert(Monomial(vars_x2), Expr::Constant(1.0));

    let mut vars_y2 = BTreeMap::new();
    vars_y2.insert("y".to_string(), 2);
    terms.insert(Monomial(vars_y2), Expr::Constant(1.0));

    let vars_0 = BTreeMap::new();
    terms.insert(Monomial(vars_0), Expr::Constant(-1.0));

    let p = SparsePolynomial { terms };

    let result = cad(&[p], &["x", "y"]).unwrap();

    // Projection phase:
    // proj_var = y
    // resultant(p, p_prime, y) where p = y^2 + (x^2 - 1)
    // p_prime = 2y
    // Sylvia matrix: [[1, 0, x^2-1], [2, 0, 0], [0, 2, 0]]
    // Wait, deg(p)=2, deg(p_prime)=1. Sylvia is 3x3.
    // Resultant = Discriminant of y^2 + (x^2-1) = -4(x^2-1)
    // Roots of -4(x^2-1) are x = -1, 1.

    // First level (x): 5 cells (-inf, -1), {-1}, (-1, 1), {1}, (1, inf)
    // For x in (-inf, -1) or (1, inf), x^2-1 > 0, so y^2 + (x^2-1) has no real roots. 1 cell above each.
    // For x = -1 or 1, x^2-1 = 0, so y^2 = 0 has 1 root y=0. 3 cells above each.
    // For x in (-1, 1), x^2-1 < 0, so y^2 + (x^2-1) has 2 roots. 5 cells above.

    // Total cells: 1 (for (-inf, -1)) + 3 (for {-1}) + 5 (for (-1, 1)) + 3 (for {1}) + 1 (for (1, inf)) = 13
    assert_eq!(result.cells.len(), 13);
}

#[test]
fn test_cad_serialization() {
    // p(x) = x
    let mut terms = BTreeMap::new();
    let mut vars_x = BTreeMap::new();
    vars_x.insert("x".to_string(), 1);
    terms.insert(Monomial(vars_x), Expr::Constant(1.0));
    let p = SparsePolynomial { terms };

    let result = cad(&[p], &["x"]).unwrap();
    let json = serde_json::to_string(&result).unwrap();
    let deserialized: Cad = serde_json::from_str(&json).unwrap();
    assert_eq!(result.cells.len(), deserialized.cells.len());
}

#[test]
#[cfg(feature = "ffi_api")]
fn test_cad_ffi_handle() {
    use rssn::ffi_apis::symbolic_cad_ffi::handle::*;
    use std::ffi::CString;

    // p(x) = x^2 - 1
    let p = Expr::new_sub(
        Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
        Expr::new_constant(1.0),
    );
    let sp = rssn::symbolic::polynomial::expr_to_sparse_poly(&p, &["x"]);
    let p_ptr = Box::into_raw(Box::new(Expr::SparsePolynomial(sp)));
    let polys = [p_ptr as *const Expr];
    let x_str = CString::new("x").unwrap();
    let vars = [x_str.as_ptr()];

    let handle = rssn_cad_handle(polys.as_ptr(), 1, vars.as_ptr(), 1);
    assert!(!handle.is_null());

    let cell_count = rssn_cad_get_cell_count(handle);
    assert_eq!(cell_count, 5);

    rssn_free_cad_handle(handle);
    unsafe {
        let _ = Box::from_raw(p_ptr);
    }
}

#[test]
#[cfg(feature = "ffi_api")]
fn test_cad_ffi_json() {
    use rssn::ffi_apis::common::rssn_free_string;
    use rssn::ffi_apis::symbolic_cad_ffi::json::*;
    use std::ffi::{CStr, CString};

    #[derive(serde::Serialize)]
    struct CadInput {
        polys: Vec<Expr>,
        vars: Vec<String>,
    }

    // p(x) = x^2 - 1
    let p = Expr::new_sub(
        Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
        Expr::new_constant(1.0),
    );
    let input = CadInput {
        polys: vec![p],
        vars: vec!["x".to_string()],
    };
    let input_json = serde_json::to_string(&input).unwrap();
    let input_c_str = CString::new(input_json).unwrap();

    let result_ptr = unsafe { rssn_json_cad(input_c_str.as_ptr()) };
    assert!(!result_ptr.is_null());

    let result_json = unsafe { CStr::from_ptr(result_ptr) }.to_str().unwrap();
    assert!(result_json.contains("\"cells\""));
    assert!(result_json.contains("5"));

    rssn_free_string(result_ptr);
}

#[test]
#[cfg(feature = "ffi_api")]
fn test_cad_ffi_bincode() {
    use rssn::ffi_apis::common::{
        from_bincode_buffer, rssn_free_bincode_buffer, to_bincode_buffer, BincodeBuffer,
    };
    use rssn::ffi_apis::symbolic_cad_ffi::bincode_api::*;
    use serde::Serialize;

    #[derive(Serialize)]
    struct CadInput {
        polys: Vec<Expr>,
        vars: Vec<String>,
    }

    // p(x) = x^2 - 1
    let p = Expr::new_sub(
        Expr::new_pow(Expr::new_variable("x"), Expr::new_constant(2.0)),
        Expr::new_constant(1.0),
    );
    let input = CadInput {
        polys: vec![p],
        vars: vec!["x".to_string()],
    };

    let input_buf = to_bincode_buffer(&input);
    assert!(!input_buf.is_null());

    let result_buf = rssn_bincode_cad(input_buf);
    assert!(!result_buf.is_null());

    let result: Option<Cad> = from_bincode_buffer(&result_buf);
    assert!(result.is_some());
    assert_eq!(result.unwrap().cells.len(), 5);

    rssn_free_bincode_buffer(input_buf);
    rssn_free_bincode_buffer(result_buf);
}
