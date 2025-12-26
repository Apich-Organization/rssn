//! JSON-based FFI API for numerical fractal geometry and chaos functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::fractal_geometry_and_chaos;

#[derive(Deserialize)]

struct MandelbrotSetInput {
    width : usize,
    height : usize,
    x_range : (f64, f64),
    y_range : (f64, f64),
    max_iter : u32,
}

#[derive(Deserialize)]

struct MandelbrotPointInput {
    c_real : f64,
    c_imag : f64,
    max_iter : u32,
}

#[derive(Deserialize)]

struct JuliaSetInput {
    width : usize,
    height : usize,
    x_range : (f64, f64),
    y_range : (f64, f64),
    c : (f64, f64),
    max_iter : u32,
}

#[derive(Deserialize)]

struct JuliaPointInput {
    z_real : f64,
    z_imag : f64,
    c_real : f64,
    c_imag : f64,
    max_iter : u32,
}

#[derive(Deserialize)]

struct LorenzInput {
    start_point : (f64, f64, f64),
    dt : f64,
    num_steps : usize,
}

#[derive(Deserialize)]

struct LorenzCustomInput {
    start_point : (f64, f64, f64),
    dt : f64,
    num_steps : usize,
    sigma : f64,
    rho : f64,
    beta : f64,
}

#[derive(Deserialize)]

struct RosslerInput {
    start_point : (f64, f64, f64),
    dt : f64,
    num_steps : usize,
    a : f64,
    b : f64,
    c : f64,
}

#[derive(Deserialize)]

struct HenonInput {
    start_point : (f64, f64),
    num_steps : usize,
    a : f64,
    b : f64,
}

#[derive(Deserialize)]

struct TinkerbellInput {
    start_point : (f64, f64),
    num_steps : usize,
    a : f64,
    b : f64,
    c : f64,
    d : f64,
}

#[derive(Deserialize)]

struct LogisticMapInput {
    x0 : f64,
    r : f64,
    num_steps : usize,
}

#[derive(Deserialize)]

struct BifurcationInput {
    r_range : (f64, f64),
    num_r_values : usize,
    transient : usize,
    num_points : usize,
    x0 : f64,
}

#[derive(Deserialize)]

struct LyapunovLogisticInput {
    r : f64,
    x0 : f64,
    transient : usize,
    num_iterations : usize,
}

#[derive(Deserialize)]

struct LyapunovLorenzInput {
    start_point : (f64, f64, f64),
    dt : f64,
    num_steps : usize,
    sigma : f64,
    rho : f64,
    beta : f64,
}

#[derive(Deserialize)]

struct DimensionInput {
    points : Vec<(f64, f64)>,
    num_scales : usize,
}

// Mandelbrot set
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_mandelbrot_set_json(
    input : *const c_char
) -> *mut c_char {

    let input: MandelbrotSetInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Vec<u32>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
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

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_mandelbrot_escape_time_json(
    input : *const c_char
) -> *mut c_char {

    let input: MandelbrotPointInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        u32,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::mandelbrot_escape_time(
            input.c_real,
            input.c_imag,
            input.max_iter,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Julia set
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_julia_set_json(
    input : *const c_char
) -> *mut c_char {

    let input: JuliaSetInput = match from_json_string(input)
    {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<Vec<u32>>,
                    String,
                > {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
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

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_julia_escape_time_json(
    input : *const c_char
) -> *mut c_char {

    let input: JuliaPointInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        u32,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::julia_escape_time(
            input.z_real,
            input.z_imag,
            input.c_real,
            input.c_imag,
            input.max_iter,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Lorenz attractor
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lorenz_attractor_json(
    input : *const c_char
) -> *mut c_char {

    let input: LorenzInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, f64, f64)>,
                    String,
                > {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        },
    };

    let result = fractal_geometry_and_chaos::generate_lorenz_attractor(
        input.start_point,
        input.dt,
        input.num_steps,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lorenz_attractor_custom_json(
    input : *const c_char
) -> *mut c_char {

    let input: LorenzCustomInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<(f64, f64, f64)>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result = fractal_geometry_and_chaos::generate_lorenz_attractor_custom(
        input.start_point,
        input.dt,
        input.num_steps,
        input.sigma,
        input.rho,
        input.beta,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Rossler attractor
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_rossler_attractor_json(
    input : *const c_char
) -> *mut c_char {

    let input: RosslerInput = match from_json_string(input)
    {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, f64, f64)>,
                    String,
                > {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        },
    };

    let result = fractal_geometry_and_chaos::generate_rossler_attractor(
        input.start_point,
        input.dt,
        input.num_steps,
        input.a,
        input.b,
        input.c,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Henon map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_henon_map_json(
    input : *const c_char
) -> *mut c_char {

    let input: HenonInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<
                    Vec<(f64, f64)>,
                    String,
                > {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        },
    };

    let result =
        fractal_geometry_and_chaos::generate_henon_map(
            input.start_point,
            input.num_steps,
            input.a,
            input.b,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Tinkerbell map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_tinkerbell_map_json(
    input : *const c_char
) -> *mut c_char {

    let input: TinkerbellInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<(f64, f64)>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::generate_tinkerbell_map(
            input.start_point,
            input.num_steps,
            input.a,
            input.b,
            input.c,
            input.d,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Logistic map
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_logistic_map_json(
    input : *const c_char
) -> *mut c_char {

    let input: LogisticMapInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<f64>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::logistic_map_iterate(
            input.x0,
            input.r,
            input.num_steps,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Bifurcation diagram
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_bifurcation_json(
    input : *const c_char
) -> *mut c_char {

    let input: BifurcationInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<(f64, f64)>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::logistic_bifurcation(
            input.r_range,
            input.num_r_values,
            input.transient,
            input.num_points,
            input.x0,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Lyapunov exponents
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lyapunov_logistic_json(
    input : *const c_char
) -> *mut c_char {

    let input: LyapunovLogisticInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        f64,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result = fractal_geometry_and_chaos::lyapunov_exponent_logistic(
        input.r,
        input.x0,
        input.transient,
        input.num_iterations,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_lyapunov_lorenz_json(
    input : *const c_char
) -> *mut c_char {

    let input: LyapunovLorenzInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        f64,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result = fractal_geometry_and_chaos::lyapunov_exponent_lorenz(
        input.start_point,
        input.dt,
        input.num_steps,
        input.sigma,
        input.rho,
        input.beta,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

// Dimension estimation
#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_box_counting_dim_json(
    input : *const c_char
) -> *mut c_char {

    let input: DimensionInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        f64,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::box_counting_dimension(
            &input.points,
            input.num_scales,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fractal_correlation_dim_json(
    input : *const c_char
) -> *mut c_char {

    let input: DimensionInput =
        match from_json_string(input) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        f64,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON".to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    let result =
        fractal_geometry_and_chaos::correlation_dimension(
            &input.points,
            input.num_scales,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok : Some(result),
                err : None::<String>,
            },
        )
        .unwrap(),
    )
}
