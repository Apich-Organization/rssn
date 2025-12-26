use crate::symbolic::core::Expr;
use crate::symbolic::fractal_geometry_and_chaos::*;

// --- IteratedFunctionSystem ---

/// Creates a new IteratedFunctionSystem (Handle)
#[no_mangle]

pub extern "C" fn rssn_ifs_create(
    functions_ptr : *const *mut Expr,
    functions_len : usize,
    probabilities_ptr : *const *mut Expr,
    probabilities_len : usize,
    variables_ptr : *const *const std::os::raw::c_char,
    variables_len : usize,
) -> *mut IteratedFunctionSystem {

    if functions_ptr.is_null()
        || probabilities_ptr.is_null()
        || variables_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let functions_slice =
            std::slice::from_raw_parts(
                functions_ptr,
                functions_len,
            );

        let functions : Vec<Expr> =
            functions_slice
                .iter()
                .map(|&p| (*p).clone())
                .collect();

        let probabilities_slice =
            std::slice::from_raw_parts(
                probabilities_ptr,
                probabilities_len,
            );

        let probabilities : Vec<Expr> =
            probabilities_slice
                .iter()
                .map(|&p| (*p).clone())
                .collect();

        let variables_slice =
            std::slice::from_raw_parts(
                variables_ptr,
                variables_len,
            );

        let variables : Vec<String> = variables_slice
            .iter()
            .filter_map(|&p| {
                if p.is_null() {

                    None
                } else {

                    std::ffi::CStr::from_ptr(p)
                        .to_str()
                        .ok()
                        .map(|s| s.to_string())
                }
            })
            .collect();

        let ifs =
            IteratedFunctionSystem::new(
                functions,
                probabilities,
                variables,
            );

        Box::into_raw(Box::new(ifs))
    }
}

/// Frees an IteratedFunctionSystem handle
#[no_mangle]

pub extern "C" fn rssn_ifs_free(
    ptr : *mut IteratedFunctionSystem
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Calculates similarity dimension (Handle)
#[no_mangle]

pub extern "C" fn rssn_ifs_similarity_dimension(
    scaling_factors_ptr : *const *mut Expr,
    len : usize,
) -> *mut Expr {

    if scaling_factors_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let slice =
            std::slice::from_raw_parts(
                scaling_factors_ptr,
                len,
            );

        let factors : Vec<Expr> = slice
            .iter()
            .map(|&p| (*p).clone())
            .collect();

        let result = IteratedFunctionSystem::similarity_dimension(&factors);

        Box::into_raw(Box::new(result))
    }
}

// --- ComplexDynamicalSystem ---

/// Creates a new Mandelbrot family system (Handle)
#[no_mangle]

pub extern "C" fn rssn_complex_system_new_mandelbrot(
    c_ptr : *const Expr
) -> *mut ComplexDynamicalSystem {

    if c_ptr.is_null() {

        return std::ptr::null_mut();
    }

    unsafe {

        let c = (*c_ptr).clone();

        let system = ComplexDynamicalSystem::new_mandelbrot_family(c);

        Box::into_raw(Box::new(system))
    }
}

/// Frees a ComplexDynamicalSystem handle
#[no_mangle]

pub extern "C" fn rssn_complex_system_free(
    ptr : *mut ComplexDynamicalSystem
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Iterates the system once (Handle)
#[no_mangle]

pub extern "C" fn rssn_complex_system_iterate(
    system_ptr : *const ComplexDynamicalSystem,
    z_ptr : *const Expr,
) -> *mut Expr {

    if system_ptr.is_null()
        || z_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let system = &*system_ptr;

        let z = &*z_ptr;

        let result = system.iterate(z);

        Box::into_raw(Box::new(result))
    }
}

/// Finds fixed points (Handle)
#[no_mangle]

pub extern "C" fn rssn_complex_system_fixed_points(
    system_ptr : *const ComplexDynamicalSystem,
    out_len : *mut usize,
) -> *mut *mut Expr {

    if system_ptr.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let system = &*system_ptr;

        let points =
            system.fixed_points();

        *out_len = points.len();

        let boxed_points : Vec<
            *mut Expr,
        > = points
            .into_iter()
            .map(|p| {
                Box::into_raw(Box::new(
                    p,
                ))
            })
            .collect();

        let ptr = boxed_points.as_ptr()
            as *mut *mut Expr;

        std::mem::forget(boxed_points);

        ptr
    }
}

// --- Chaos Theory Functions ---

/// Finds fixed points of a 1D map (Handle)
#[no_mangle]

pub extern "C" fn rssn_find_fixed_points(
    map_ptr : *const Expr,
    var : *const std::os::raw::c_char,
    out_len : *mut usize,
) -> *mut *mut Expr {

    if map_ptr.is_null()
        || var.is_null()
        || out_len.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let map = &*map_ptr;

        let var_str =
            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .unwrap_or("x");

        let points = find_fixed_points(
            map,
            var_str,
        );

        *out_len = points.len();

        let boxed_points : Vec<
            *mut Expr,
        > = points
            .into_iter()
            .map(|p| {
                Box::into_raw(Box::new(
                    p,
                ))
            })
            .collect();

        let ptr = boxed_points.as_ptr()
            as *mut *mut Expr;

        std::mem::forget(boxed_points);

        ptr
    }
}

/// Analyzes stability of a fixed point (Handle)
#[no_mangle]

pub extern "C" fn rssn_analyze_stability(
    map_ptr : *const Expr,
    var : *const std::os::raw::c_char,
    fixed_point_ptr : *const Expr,
) -> *mut Expr {

    if map_ptr.is_null()
        || var.is_null()
        || fixed_point_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let map = &*map_ptr;

        let var_str =
            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .unwrap_or("x");

        let fixed_point =
            &*fixed_point_ptr;

        let result = analyze_stability(
            map,
            var_str,
            fixed_point,
        );

        Box::into_raw(Box::new(result))
    }
}

/// Calculates Lyapunov exponent (Handle)
#[no_mangle]

pub extern "C" fn rssn_lyapunov_exponent(
    map_ptr : *const Expr,
    var : *const std::os::raw::c_char,
    initial_x_ptr : *const Expr,
    n_iterations : usize,
) -> *mut Expr {

    if map_ptr.is_null()
        || var.is_null()
        || initial_x_ptr.is_null()
    {

        return std::ptr::null_mut();
    }

    unsafe {

        let map = &*map_ptr;

        let var_str =
            std::ffi::CStr::from_ptr(
                var,
            )
            .to_str()
            .unwrap_or("x");

        let initial_x = &*initial_x_ptr;

        let result = lyapunov_exponent(
            map,
            var_str,
            initial_x,
            n_iterations,
        );

        Box::into_raw(Box::new(result))
    }
}

/// Returns Lorenz system equations (Handle)
#[no_mangle]

pub extern "C" fn rssn_lorenz_system(
    dx_out : *mut *mut Expr,
    dy_out : *mut *mut Expr,
    dz_out : *mut *mut Expr,
) -> bool {

    if dx_out.is_null()
        || dy_out.is_null()
        || dz_out.is_null()
    {

        return false;
    }

    unsafe {

        let (dx, dy, dz) =
            lorenz_system();

        *dx_out =
            Box::into_raw(Box::new(dx));

        *dy_out =
            Box::into_raw(Box::new(dy));

        *dz_out =
            Box::into_raw(Box::new(dz));

        true
    }
}
