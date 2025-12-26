use crate::ffi_apis::common::*;
use crate::symbolic::complex_analysis::*;
use crate::symbolic::core::Expr;
use std::os::raw::c_char;

#[no_mangle]

pub unsafe extern "C" fn path_continuation_new(
    func: *const Expr,
    var: *const c_char,
    start_point: *const Expr,
    order: usize,
) -> *mut PathContinuation {

    if func.is_null() || start_point.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let start_point_ref = &*start_point;

    let path_continuation = PathContinuation::new(func_ref, var_str, start_point_ref, order);

    Box::into_raw(Box::new(path_continuation))
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_continue_along_path(
    pc: *mut PathContinuation,
    path_points: *const *const Expr,
    path_points_len: usize,
) -> *mut c_char {

    if pc.is_null() || path_points.is_null() {

        return std::ptr::null_mut();
    }

    let pc_ref = &mut *pc;

    let path_points_slice = std::slice::from_raw_parts(path_points, path_points_len);

    let path_points_vec: Vec<Expr> = path_points_slice
        .iter()
        .map(|&ptr| (*ptr).clone())
        .collect();

    match pc_ref.continue_along_path(&path_points_vec) {
        Ok(_) => to_c_string("OK".to_string()),
        Err(e) => to_c_string(e),
    }
}

#[no_mangle]

pub unsafe extern "C" fn path_continuation_get_final_expression(
    pc: *const PathContinuation
) -> *mut Expr {

    if pc.is_null() {

        return std::ptr::null_mut();
    }

    let pc_ref = &*pc;

    match pc_ref.get_final_expression() {
        Some(expr) => Box::into_raw(Box::new(expr.clone())),
        None => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn estimate_radius_of_convergence(
    series_expr: *const Expr,
    var: *const c_char,
    center: *const Expr,
    order: usize,
) -> f64 {

    if series_expr.is_null() || center.is_null() {

        return 0.0;
    }

    let series_expr_ref = &*series_expr;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let center_ref = &*center;

    crate::symbolic::complex_analysis::estimate_radius_of_convergence(
        series_expr_ref,
        var_str,
        center_ref,
        order,
    )
    .unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn complex_distance(
    p1: *const Expr,
    p2: *const Expr,
) -> f64 {

    if p1.is_null() || p2.is_null() {

        return 0.0;
    }

    let p1_ref = &*p1;

    let p2_ref = &*p2;

    crate::symbolic::complex_analysis::complex_distance(p1_ref, p2_ref).unwrap_or(0.0)
}

#[no_mangle]

pub unsafe extern "C" fn classify_singularity(
    func: *const Expr,
    var: *const c_char,
    singularity: *const Expr,
    order: usize,
) -> *mut SingularityType {

    if func.is_null() || singularity.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let singularity_ref = &*singularity;

    let singularity_type = crate::symbolic::complex_analysis::classify_singularity(
        func_ref,
        var_str,
        singularity_ref,
        order,
    );

    Box::into_raw(Box::new(singularity_type))
}

#[no_mangle]

pub unsafe extern "C" fn laurent_series(
    func: *const Expr,
    var: *const c_char,
    center: *const Expr,
    order: usize,
) -> *mut Expr {

    if func.is_null() || center.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let center_ref = &*center;

    let series =
        crate::symbolic::complex_analysis::laurent_series(func_ref, var_str, center_ref, order);

    Box::into_raw(Box::new(series))
}

#[no_mangle]

pub unsafe extern "C" fn calculate_residue(
    func: *const Expr,
    var: *const c_char,
    singularity: *const Expr,
) -> *mut Expr {

    if func.is_null() || singularity.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let singularity_ref = &*singularity;

    let residue =
        crate::symbolic::complex_analysis::calculate_residue(func_ref, var_str, singularity_ref);

    Box::into_raw(Box::new(residue))
}

#[no_mangle]

pub unsafe extern "C" fn contour_integral_residue_theorem(
    func: *const Expr,
    var: *const c_char,
    singularities: *const *const Expr,
    singularities_len: usize,
) -> *mut Expr {

    if func.is_null() || singularities.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let singularities_slice = std::slice::from_raw_parts(singularities, singularities_len);

    let singularities_vec: Vec<Expr> = singularities_slice
        .iter()
        .map(|&ptr| (*ptr).clone())
        .collect();

    let result = crate::symbolic::complex_analysis::contour_integral_residue_theorem(
        func_ref,
        var_str,
        &singularities_vec,
    );

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_new(
    a: *const Expr,
    b: *const Expr,
    c: *const Expr,
    d: *const Expr,
) -> *mut MobiusTransformation {

    if a.is_null() || b.is_null() || c.is_null() || d.is_null() {

        return std::ptr::null_mut();
    }

    let a_ref = &*a;

    let b_ref = &*b;

    let c_ref = &*c;

    let d_ref = &*d;

    let mobius =
        MobiusTransformation::new(a_ref.clone(), b_ref.clone(), c_ref.clone(), d_ref.clone());

    Box::into_raw(Box::new(mobius))
}

#[no_mangle]

pub extern "C" fn mobius_transformation_identity() -> *mut MobiusTransformation {

    let mobius = MobiusTransformation::identity();

    Box::into_raw(Box::new(mobius))
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_apply(
    mobius: *const MobiusTransformation,
    z: *const Expr,
) -> *mut Expr {

    if mobius.is_null() || z.is_null() {

        return std::ptr::null_mut();
    }

    let mobius_ref = &*mobius;

    let z_ref = &*z;

    let result = mobius_ref.apply(z_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_compose(
    mobius1: *const MobiusTransformation,
    mobius2: *const MobiusTransformation,
) -> *mut MobiusTransformation {

    if mobius1.is_null() || mobius2.is_null() {

        return std::ptr::null_mut();
    }

    let mobius1_ref = &*mobius1;

    let mobius2_ref = &*mobius2;

    let result = mobius1_ref.compose(mobius2_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn mobius_transformation_inverse(
    mobius: *const MobiusTransformation
) -> *mut MobiusTransformation {

    if mobius.is_null() {

        return std::ptr::null_mut();
    }

    let mobius_ref = &*mobius;

    let result = mobius_ref.inverse();

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_integral_formula(
    func: *const Expr,
    var: *const c_char,
    z0: *const Expr,
) -> *mut Expr {

    if func.is_null() || z0.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let z0_ref = &*z0;

    let result =
        crate::symbolic::complex_analysis::cauchy_integral_formula(func_ref, var_str, z0_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn cauchy_derivative_formula(
    func: *const Expr,
    var: *const c_char,
    z0: *const Expr,
    n: usize,
) -> *mut Expr {

    if func.is_null() || z0.is_null() {

        return std::ptr::null_mut();
    }

    let func_ref = &*func;

    let var_str = std::ffi::CStr::from_ptr(var)
        .to_str()
        .unwrap();

    let z0_ref = &*z0;

    let result =
        crate::symbolic::complex_analysis::cauchy_derivative_formula(func_ref, var_str, z0_ref, n);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn complex_exp(z: *const Expr) -> *mut Expr {

    if z.is_null() {

        return std::ptr::null_mut();
    }

    let z_ref = &*z;

    let result = crate::symbolic::complex_analysis::complex_exp(z_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn complex_log(z: *const Expr) -> *mut Expr {

    if z.is_null() {

        return std::ptr::null_mut();
    }

    let z_ref = &*z;

    let result = crate::symbolic::complex_analysis::complex_log(z_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn complex_arg(z: *const Expr) -> *mut Expr {

    if z.is_null() {

        return std::ptr::null_mut();
    }

    let z_ref = &*z;

    let result = crate::symbolic::complex_analysis::complex_arg(z_ref);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn complex_modulus(z: *const Expr) -> *mut Expr {

    if z.is_null() {

        return std::ptr::null_mut();
    }

    let z_ref = &*z;

    let result = crate::symbolic::complex_analysis::complex_modulus(z_ref);

    Box::into_raw(Box::new(result))
}
