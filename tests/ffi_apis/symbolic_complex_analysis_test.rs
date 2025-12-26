use std::ffi::CString;
use std::ptr;

use rssn::ffi_apis::symbolic_complex_analysis_ffi::handle::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_path_continuation_new() {

    let func = Expr::new_sin(
        Expr::Variable("z".to_string()),
    );

    let var =
        CString::new("z").unwrap();

    let start_point =
        Expr::Constant(0.0);

    let order = 5;

    let pc = unsafe {

        path_continuation_new(
            &func as *const Expr,
            var.as_ptr(),
            &start_point as *const Expr,
            order,
        )
    };

    assert!(!pc.is_null());

    let pc_ref = unsafe {

        &*pc
    };

    assert_eq!(pc_ref.var, "z");

    assert_eq!(pc_ref.order, 5);

    assert_eq!(
        pc_ref.pieces.len(),
        1
    );

    unsafe {

        Box::from_raw(pc);
    }
}
