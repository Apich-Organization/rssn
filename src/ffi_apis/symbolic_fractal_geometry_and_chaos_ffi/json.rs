use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::fractal_geometry_and_chaos::*;

/// Creates a new IteratedFunctionSystem (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_ifs_create(
    functions_json: *const std::os::raw::c_char,
    probabilities_json: *const std::os::raw::c_char,
    variables_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let functions: Option<Vec<Expr>> = from_json_string(functions_json);
    let probabilities: Option<Vec<Expr>> = from_json_string(probabilities_json);
    let variables: Option<Vec<String>> = from_json_string(variables_json);

    if let (Some(f), Some(p), Some(v)) = (functions, probabilities, variables) {
        let ifs = IteratedFunctionSystem::new(f, p, v);
        to_json_string(&ifs)
    } else {
        std::ptr::null_mut()
    }
}

/// Calculates similarity dimension (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_ifs_similarity_dimension(
    scaling_factors_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let factors: Option<Vec<Expr>> = from_json_string(scaling_factors_json);

    if let Some(f) = factors {
        let result = IteratedFunctionSystem::similarity_dimension(&f);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Creates a new Mandelbrot family system (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_complex_system_new_mandelbrot(
    c_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let c: Option<Expr> = from_json_string(c_json);

    if let Some(c_val) = c {
        let system = ComplexDynamicalSystem::new_mandelbrot_family(c_val);
        to_json_string(&system)
    } else {
        std::ptr::null_mut()
    }
}

/// Iterates the system once (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_complex_system_iterate(
    system_json: *const std::os::raw::c_char,
    z_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let system: Option<ComplexDynamicalSystem> = from_json_string(system_json);
    let z: Option<Expr> = from_json_string(z_json);

    if let (Some(sys), Some(z_val)) = (system, z) {
        let result = sys.iterate(&z_val);
        to_json_string(&result)
    } else {
        std::ptr::null_mut()
    }
}

/// Finds fixed points (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_complex_system_fixed_points(
    system_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let system: Option<ComplexDynamicalSystem> = from_json_string(system_json);

    if let Some(sys) = system {
        let points = sys.fixed_points();
        to_json_string(&points)
    } else {
        std::ptr::null_mut()
    }
}

/// Finds fixed points of a 1D map (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_find_fixed_points(
    map_json: *const std::os::raw::c_char,
    var: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let map: Option<Expr> = from_json_string(map_json);

    if map.is_none() || var.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap_or("x");
        let points = find_fixed_points(&map.unwrap(), var_str);
        to_json_string(&points)
    }
}

/// Analyzes stability of a fixed point (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_analyze_stability(
    map_json: *const std::os::raw::c_char,
    var: *const std::os::raw::c_char,
    fixed_point_json: *const std::os::raw::c_char,
) -> *mut std::os::raw::c_char {
    let map: Option<Expr> = from_json_string(map_json);
    let fixed_point: Option<Expr> = from_json_string(fixed_point_json);

    if map.is_none() || var.is_null() || fixed_point.is_none() {
        return std::ptr::null_mut();
    }

    unsafe {
        let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap_or("x");
        let result = analyze_stability(&map.unwrap(), var_str, &fixed_point.unwrap());
        to_json_string(&result)
    }
}

/// Calculates Lyapunov exponent (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_lyapunov_exponent(
    map_json: *const std::os::raw::c_char,
    var: *const std::os::raw::c_char,
    initial_x_json: *const std::os::raw::c_char,
    n_iterations: usize,
) -> *mut std::os::raw::c_char {
    let map: Option<Expr> = from_json_string(map_json);
    let initial_x: Option<Expr> = from_json_string(initial_x_json);

    if map.is_none() || var.is_null() || initial_x.is_none() {
        return std::ptr::null_mut();
    }

    unsafe {
        let var_str = std::ffi::CStr::from_ptr(var).to_str().unwrap_or("x");
        let result = lyapunov_exponent(&map.unwrap(), var_str, &initial_x.unwrap(), n_iterations);
        to_json_string(&result)
    }
}

/// Returns Lorenz system equations (JSON)
#[no_mangle]
pub extern "C" fn rssn_json_lorenz_system() -> *mut std::os::raw::c_char {
    let (dx, dy, dz) = lorenz_system();
    let result = vec![dx, dy, dz];
    to_json_string(&result)
}
