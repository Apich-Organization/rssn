use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::fractal_geometry_and_chaos::*;

/// Creates a new IteratedFunctionSystem (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_ifs_create(
    functions_buf: BincodeBuffer,
    probabilities_buf: BincodeBuffer,
    variables_buf: BincodeBuffer,
) -> BincodeBuffer {

    let functions: Option<Vec<Expr>> = from_bincode_buffer(&functions_buf);

    let probabilities: Option<Vec<Expr>> = from_bincode_buffer(&probabilities_buf);

    let variables: Option<Vec<String>> = from_bincode_buffer(&variables_buf);

    if let (Some(f), Some(p), Some(v)) = (functions, probabilities, variables) {

        let ifs = IteratedFunctionSystem::new(f, p, v);

        to_bincode_buffer(&ifs)
    } else {

        BincodeBuffer::empty()
    }
}

/// Calculates similarity dimension (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_ifs_similarity_dimension(
    scaling_factors_buf: BincodeBuffer,
) -> BincodeBuffer {

    let factors: Option<Vec<Expr>> = from_bincode_buffer(&scaling_factors_buf);

    if let Some(f) = factors {

        let result = IteratedFunctionSystem::similarity_dimension(&f);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Creates a new Mandelbrot family system (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_complex_system_new_mandelbrot(
    c_buf: BincodeBuffer,
) -> BincodeBuffer {

    let c: Option<Expr> = from_bincode_buffer(&c_buf);

    if let Some(c_val) = c {

        let system = ComplexDynamicalSystem::new_mandelbrot_family(c_val);

        to_bincode_buffer(&system)
    } else {

        BincodeBuffer::empty()
    }
}

/// Iterates the system once (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_complex_system_iterate(
    system_buf: BincodeBuffer,
    z_buf: BincodeBuffer,
) -> BincodeBuffer {

    let system: Option<ComplexDynamicalSystem> = from_bincode_buffer(&system_buf);

    let z: Option<Expr> = from_bincode_buffer(&z_buf);

    if let (Some(sys), Some(z_val)) = (system, z) {

        let result = sys.iterate(&z_val);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Finds fixed points (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_complex_system_fixed_points(
    system_buf: BincodeBuffer,
) -> BincodeBuffer {

    let system: Option<ComplexDynamicalSystem> = from_bincode_buffer(&system_buf);

    if let Some(sys) = system {

        let points = sys.fixed_points();

        to_bincode_buffer(&points)
    } else {

        BincodeBuffer::empty()
    }
}

/// Finds fixed points of a 1D map (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_find_fixed_points(
    map_buf: BincodeBuffer,
    var: *const std::os::raw::c_char,
) -> BincodeBuffer {

    let map: Option<Expr> = from_bincode_buffer(&map_buf);

    if map.is_none() || var.is_null() {

        return BincodeBuffer::empty();
    }

    unsafe {

        let var_str = std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap_or("x");

        let points = find_fixed_points(&map.unwrap(), var_str);

        to_bincode_buffer(&points)
    }
}

/// Analyzes stability of a fixed point (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_analyze_stability(
    map_buf: BincodeBuffer,
    var: *const std::os::raw::c_char,
    fixed_point_buf: BincodeBuffer,
) -> BincodeBuffer {

    let map: Option<Expr> = from_bincode_buffer(&map_buf);

    let fixed_point: Option<Expr> = from_bincode_buffer(&fixed_point_buf);

    if map.is_none() || var.is_null() || fixed_point.is_none() {

        return BincodeBuffer::empty();
    }

    unsafe {

        let var_str = std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap_or("x");

        let result = analyze_stability(&map.unwrap(), var_str, &fixed_point.unwrap());

        to_bincode_buffer(&result)
    }
}

/// Calculates Lyapunov exponent (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_lyapunov_exponent(
    map_buf: BincodeBuffer,
    var: *const std::os::raw::c_char,
    initial_x_buf: BincodeBuffer,
    n_iterations: usize,
) -> BincodeBuffer {

    let map: Option<Expr> = from_bincode_buffer(&map_buf);

    let initial_x: Option<Expr> = from_bincode_buffer(&initial_x_buf);

    if map.is_none() || var.is_null() || initial_x.is_none() {

        return BincodeBuffer::empty();
    }

    unsafe {

        let var_str = std::ffi::CStr::from_ptr(var)
            .to_str()
            .unwrap_or("x");

        let result = lyapunov_exponent(&map.unwrap(), var_str, &initial_x.unwrap(), n_iterations);

        to_bincode_buffer(&result)
    }
}

/// Returns Lorenz system equations (Bincode)
#[no_mangle]

pub extern "C" fn rssn_bincode_lorenz_system() -> BincodeBuffer {

    let (dx, dy, dz) = lorenz_system();

    let result = vec![dx, dy, dz];

    to_bincode_buffer(&result)
}
