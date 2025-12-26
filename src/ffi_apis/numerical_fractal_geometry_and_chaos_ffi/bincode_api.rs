//! Bincode-based FFI API for numerical fractal geometry and chaos functions.

use serde::Deserialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::fractal_geometry_and_chaos;

#[derive(Deserialize)]

struct MandelbrotSetInput {
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    max_iter: u32,
}

#[derive(Deserialize)]

struct MandelbrotPointInput {
    c_real: f64,
    c_imag: f64,
    max_iter: u32,
}

#[derive(Deserialize)]

struct JuliaSetInput {
    width: usize,
    height: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    c: (f64, f64),
    max_iter: u32,
}

#[derive(Deserialize)]

struct LorenzInput {
    start_point: (f64, f64, f64),
    dt: f64,
    num_steps: usize,
}

#[derive(Deserialize)]

struct HenonInput {
    start_point: (f64, f64),
    num_steps: usize,
    a: f64,
    b: f64,
}

#[derive(Deserialize)]

struct LogisticMapInput {
    x0: f64,
    r: f64,
    num_steps: usize,
}

#[derive(Deserialize)]

struct LyapunovLogisticInput {
    r: f64,
    x0: f64,
    transient: usize,
    num_iterations: usize,
}

#[derive(Deserialize)]

struct DimensionInput {
    points: Vec<(f64, f64)>,
    num_scales: usize,
}

// Mandelbrot set
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_mandelbrot_set_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: MandelbrotSetInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<Vec<u32>>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::generate_mandelbrot_set(
            input.width,
            input.height,
            input.x_range,
            input.y_range,
            input.max_iter,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_mandelbrot_escape_time_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: MandelbrotPointInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    u32,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::mandelbrot_escape_time(
            input.c_real,
            input.c_imag,
            input.max_iter,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Julia set
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_julia_set_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: JuliaSetInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<Vec<u32>>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::generate_julia_set(
            input.width,
            input.height,
            input.x_range,
            input.y_range,
            input.c,
            input.max_iter,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Lorenz attractor
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lorenz_attractor_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LorenzInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<(f64, f64, f64)>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result = fractal_geometry_and_chaos::generate_lorenz_attractor(
        input.start_point,
        input.dt,
        input.num_steps,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Henon map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_henon_map_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: HenonInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<(f64, f64)>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::generate_henon_map(
            input.start_point,
            input.num_steps,
            input.a,
            input.b,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Logistic map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_logistic_map_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LogisticMapInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::logistic_map_iterate(
            input.x0,
            input.r,
            input.num_steps,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Lyapunov exponent for logistic map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lyapunov_logistic_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LyapunovLogisticInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result = fractal_geometry_and_chaos::lyapunov_exponent_logistic(
        input.r,
        input.x0,
        input.transient,
        input.num_iterations,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Box-counting dimension
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_box_counting_dim_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DimensionInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::box_counting_dimension(
            &input.points,
            input.num_scales,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

// Correlation dimension
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_correlation_dim_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DimensionInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        fractal_geometry_and_chaos::correlation_dimension(
            &input.points,
            input.num_scales,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}
