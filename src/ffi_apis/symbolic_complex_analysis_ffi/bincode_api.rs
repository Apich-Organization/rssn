use crate::ffi_apis::common::*;
use crate::symbolic::complex_analysis::*;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

#[no_mangle]

pub unsafe extern "C" fn path_continuation_new_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    start_point_bincode: BincodeBuffer,
    order: usize,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let start_point: Expr = match from_bincode_buffer(&start_point_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let path_continuation = PathContinuation::new(&func, var_str, &start_point, order);

    to_bincode_buffer(&path_continuation)
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_continue_along_path_bincode(
    pc_bincode: BincodeBuffer,
    path_points_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let mut pc: PathContinuation = match from_bincode_buffer(&pc_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let path_points: Vec<Expr> = match from_bincode_buffer(&path_points_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    match pc.continue_along_path(&path_points) {
        Ok(_) => to_bincode_buffer(&"OK".to_string()),
        Err(e) => to_bincode_buffer(&e),
    }
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_get_final_expression_bincode(
    pc_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let pc: PathContinuation = match from_bincode_buffer(&pc_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&pc.get_final_expression())
}

#[no_mangle]

pub unsafe extern "C" fn estimate_radius_of_convergence_bincode(
    series_expr_bincode: BincodeBuffer,
    var: *const c_char,
    center_bincode: BincodeBuffer,
    order: usize,
) -> f64 {

    let series_expr: Expr = match from_bincode_buffer(&series_expr_bincode) {
        Some(e) => e,
        None => return 0.0,
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let center: Expr = match from_bincode_buffer(&center_bincode) {
        Some(e) => e,
        None => return 0.0,
    };

    crate::symbolic::complex_analysis::estimate_radius_of_convergence(
        &series_expr,
        var_str,
        &center,
        order,
    )
    .unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn complex_distance_bincode(
    p1_bincode: BincodeBuffer,
    p2_bincode: BincodeBuffer,
) -> f64 {

    let p1: Expr = match from_bincode_buffer(&p1_bincode) {
        Some(e) => e,
        None => return 0.0,
    };

    let p2: Expr = match from_bincode_buffer(&p2_bincode) {
        Some(e) => e,
        None => return 0.0,
    };

    crate::symbolic::complex_analysis::complex_distance(&p1, &p2).unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn classify_singularity_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    singularity_bincode: BincodeBuffer,
    order: usize,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let singularity: Expr = match from_bincode_buffer(&singularity_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let singularity_type = crate::symbolic::complex_analysis::classify_singularity(
        &func,
        var_str,
        &singularity,
        order,
    );

    to_bincode_buffer(&singularity_type)
}

#[no_mangle]

pub unsafe extern "C" fn laurent_series_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    center_bincode: BincodeBuffer,
    order: usize,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let center: Expr = match from_bincode_buffer(&center_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let series = crate::symbolic::complex_analysis::laurent_series(&func, var_str, &center, order);

    to_bincode_buffer(&series)
}

#[no_mangle]

pub unsafe extern "C" fn calculate_residue_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    singularity_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let singularity: Expr = match from_bincode_buffer(&singularity_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let residue =
        crate::symbolic::complex_analysis::calculate_residue(&func, var_str, &singularity);

    to_bincode_buffer(&residue)
}

#[no_mangle]

pub unsafe extern "C" fn contour_integral_residue_theorem_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    singularities_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let singularities: Vec<Expr> = match from_bincode_buffer(&singularities_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::contour_integral_residue_theorem(
        &func,
        var_str,
        &singularities,
    );

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_new_bincode(
    a_bincode: BincodeBuffer,
    b_bincode: BincodeBuffer,
    c_bincode: BincodeBuffer,
    d_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let a: Expr = match from_bincode_buffer(&a_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let b: Expr = match from_bincode_buffer(&b_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let c: Expr = match from_bincode_buffer(&c_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let d: Expr = match from_bincode_buffer(&d_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let mobius = MobiusTransformation::new(a, b, c, d);

    to_bincode_buffer(&mobius)
}

#[no_mangle]

pub extern "C" fn mobius_transformation_identity_bincode() -> BincodeBuffer {

    let mobius = MobiusTransformation::identity();

    to_bincode_buffer(&mobius)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_apply_bincode(
    mobius_bincode: BincodeBuffer,
    z_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let mobius: MobiusTransformation = match from_bincode_buffer(&mobius_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let z: Expr = match from_bincode_buffer(&z_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = mobius.apply(&z);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_compose_bincode(
    mobius1_bincode: BincodeBuffer,
    mobius2_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let mobius1: MobiusTransformation = match from_bincode_buffer(&mobius1_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let mobius2: MobiusTransformation = match from_bincode_buffer(&mobius2_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = mobius1.compose(&mobius2);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_inverse_bincode(
    mobius_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let mobius: MobiusTransformation = match from_bincode_buffer(&mobius_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = mobius.inverse();

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_integral_formula_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    z0_bincode: BincodeBuffer,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let z0: Expr = match from_bincode_buffer(&z0_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::cauchy_integral_formula(&func, var_str, &z0);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_derivative_formula_bincode(
    func_bincode: BincodeBuffer,
    var: *const c_char,
    z0_bincode: BincodeBuffer,
    n: usize,
) -> BincodeBuffer {

    let func: Expr = match from_bincode_buffer(&func_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap();

    let z0: Expr = match from_bincode_buffer(&z0_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result =
        crate::symbolic::complex_analysis::cauchy_derivative_formula(&func, var_str, &z0, n);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_exp_bincode(z_bincode: BincodeBuffer) -> BincodeBuffer {

    let z: Expr = match from_bincode_buffer(&z_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::complex_exp(&z);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_log_bincode(z_bincode: BincodeBuffer) -> BincodeBuffer {

    let z: Expr = match from_bincode_buffer(&z_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::complex_log(&z);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_arg_bincode(z_bincode: BincodeBuffer) -> BincodeBuffer {

    let z: Expr = match from_bincode_buffer(&z_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::complex_arg(&z);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn complex_modulus_bincode(z_bincode: BincodeBuffer) -> BincodeBuffer {

    let z: Expr = match from_bincode_buffer(&z_bincode) {
        Some(e) => e,
        None => return BincodeBuffer::empty(),
    };

    let result = crate::symbolic::complex_analysis::complex_modulus(&z);

    to_bincode_buffer(&result)
}
