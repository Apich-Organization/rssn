//! JSON-based FFI API for numerical computer graphics functions.

use crate::ffi_apis::common::{from_json_string, to_c_string};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::computer_graphics;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

#[derive(Deserialize)]
struct Vector3DInput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]
struct TwoVectors3DInput {
    v1: Vector3DInput,
    v2: Vector3DInput,
}

#[derive(Serialize)]
struct Vector3DOutput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]
struct ReflectInput {
    incident: Vector3DInput,
    normal: Vector3DInput,
}

#[derive(Deserialize)]
struct TransformInput {
    dx: f64,
    dy: f64,
    dz: f64,
}

#[derive(Deserialize)]
struct ScaleInput {
    sx: f64,
    sy: f64,
    sz: f64,
}

#[derive(Deserialize)]
struct RotationAxisInput {
    axis: Vector3DInput,
    angle: f64,
}

#[derive(Deserialize)]
struct AngleInput {
    angle: f64,
}

#[derive(Deserialize)]
struct QuaternionInput {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Serialize)]
struct QuaternionOutput {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]
struct TwoQuaternionsInput {
    q1: QuaternionInput,
    q2: QuaternionInput,
}

#[derive(Deserialize)]
struct QuaternionSlerpInput {
    q1: QuaternionInput,
    q2: QuaternionInput,
    t: f64,
}

#[derive(Deserialize)]
struct RaySphereInput {
    ray_origin: Vector3DInput,
    ray_direction: Vector3DInput,
    sphere_center: Vector3DInput,
    sphere_radius: f64,
}

#[derive(Serialize)]
struct IntersectionOutput {
    t: f64,
    point: Vector3DOutput,
    normal: Vector3DOutput,
}

#[derive(Deserialize)]
struct Point3DInput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Serialize)]
struct Point3DOutput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]
struct BezierCubicInput {
    p0: Point3DInput,
    p1: Point3DInput,
    p2: Point3DInput,
    p3: Point3DInput,
    t: f64,
}

#[derive(Deserialize)]
struct LookAtInput {
    eye: Vector3DInput,
    center: Vector3DInput,
    up: Vector3DInput,
}

#[derive(Deserialize)]
struct PerspectiveInput {
    fov_y: f64,
    aspect: f64,
    near: f64,
    far: f64,
}

// Vector operations
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_dot_product_json(input: *const c_char) -> *mut c_char {
    let input: TwoVectors3DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let v1 = computer_graphics::Vector3D::new(input.v1.x, input.v1.y, input.v1.z);
    let v2 = computer_graphics::Vector3D::new(input.v2.x, input.v2.y, input.v2.z);
    let result = computer_graphics::dot_product(&v1, &v2);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_cross_product_json(input: *const c_char) -> *mut c_char {
    let input: TwoVectors3DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vector3DOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let v1 = computer_graphics::Vector3D::new(input.v1.x, input.v1.y, input.v1.z);
    let v2 = computer_graphics::Vector3D::new(input.v2.x, input.v2.y, input.v2.z);
    let result = computer_graphics::cross_product(&v1, &v2);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(Vector3DOutput {
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_normalize_json(input: *const c_char) -> *mut c_char {
    let input: Vector3DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vector3DOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let v = computer_graphics::Vector3D::new(input.x, input.y, input.z);
    let result = v.normalize();
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(Vector3DOutput {
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_magnitude_json(input: *const c_char) -> *mut c_char {
    let input: Vector3DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let v = computer_graphics::Vector3D::new(input.x, input.y, input.z);
    let result = v.magnitude();
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_reflect_json(input: *const c_char) -> *mut c_char {
    let input: ReflectInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vector3DOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let incident =
        computer_graphics::Vector3D::new(input.incident.x, input.incident.y, input.incident.z);
    let normal = computer_graphics::Vector3D::new(input.normal.x, input.normal.y, input.normal.z);
    let result = computer_graphics::reflect(&incident, &normal);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(Vector3DOutput {
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_angle_between_json(input: *const c_char) -> *mut c_char {
    let input: TwoVectors3DInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let v1 = computer_graphics::Vector3D::new(input.v1.x, input.v1.y, input.v1.z);
    let v2 = computer_graphics::Vector3D::new(input.v2.x, input.v2.y, input.v2.z);
    let result = computer_graphics::angle_between(&v1, &v2);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(result),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// Transformation matrices
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_translation_matrix_json(
    input: *const c_char,
) -> *mut c_char {
    let input: TransformInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix = computer_graphics::translation_matrix(input.dx, input.dy, input.dz);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_scaling_matrix_json(
    input: *const c_char,
) -> *mut c_char {
    let input: ScaleInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix = computer_graphics::scaling_matrix(input.sx, input.sy, input.sz);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_x_json(
    input: *const c_char,
) -> *mut c_char {
    let input: AngleInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix = computer_graphics::rotation_matrix_x(input.angle);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_y_json(
    input: *const c_char,
) -> *mut c_char {
    let input: AngleInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix = computer_graphics::rotation_matrix_y(input.angle);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_z_json(
    input: *const c_char,
) -> *mut c_char {
    let input: AngleInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix = computer_graphics::rotation_matrix_z(input.angle);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_axis_json(
    input: *const c_char,
) -> *mut c_char {
    let input: RotationAxisInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let axis = computer_graphics::Vector3D::new(input.axis.x, input.axis.y, input.axis.z);
    let matrix = computer_graphics::rotation_matrix_axis(&axis, input.angle);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// Quaternion operations
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_quaternion_multiply_json(
    input: *const c_char,
) -> *mut c_char {
    let input: TwoQuaternionsInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<QuaternionOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let q1 = computer_graphics::Quaternion::new(input.q1.w, input.q1.x, input.q1.y, input.q1.z);
    let q2 = computer_graphics::Quaternion::new(input.q2.w, input.q2.x, input.q2.y, input.q2.z);
    let result = q1.multiply(&q2);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(QuaternionOutput {
                w: result.w,
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_quaternion_slerp_json(
    input: *const c_char,
) -> *mut c_char {
    let input: QuaternionSlerpInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<QuaternionOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let q1 = computer_graphics::Quaternion::new(input.q1.w, input.q1.x, input.q1.y, input.q1.z);
    let q2 = computer_graphics::Quaternion::new(input.q2.w, input.q2.x, input.q2.y, input.q2.z);
    let result = q1.slerp(&q2, input.t);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(QuaternionOutput {
                w: result.w,
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// Ray tracing
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_ray_sphere_intersection_json(
    input: *const c_char,
) -> *mut c_char {
    let input: RaySphereInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Option<IntersectionOutput>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let ray = computer_graphics::Ray::new(
        computer_graphics::Point3D::new(input.ray_origin.x, input.ray_origin.y, input.ray_origin.z),
        computer_graphics::Vector3D::new(
            input.ray_direction.x,
            input.ray_direction.y,
            input.ray_direction.z,
        ),
    );
    let sphere = computer_graphics::Sphere::new(
        computer_graphics::Point3D::new(
            input.sphere_center.x,
            input.sphere_center.y,
            input.sphere_center.z,
        ),
        input.sphere_radius,
    );
    let result = computer_graphics::ray_sphere_intersection(&ray, &sphere);
    let output = result.map(|i| IntersectionOutput {
        t: i.t,
        point: Vector3DOutput {
            x: i.point.x,
            y: i.point.y,
            z: i.point.z,
        },
        normal: Vector3DOutput {
            x: i.normal.x,
            y: i.normal.y,
            z: i.normal.z,
        },
    });
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(output),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// Curves
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_bezier_cubic_json(input: *const c_char) -> *mut c_char {
    let input: BezierCubicInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Point3DOutput, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let p0 = computer_graphics::Point3D::new(input.p0.x, input.p0.y, input.p0.z);
    let p1 = computer_graphics::Point3D::new(input.p1.x, input.p1.y, input.p1.z);
    let p2 = computer_graphics::Point3D::new(input.p2.x, input.p2.y, input.p2.z);
    let p3 = computer_graphics::Point3D::new(input.p3.x, input.p3.y, input.p3.z);
    let result = computer_graphics::bezier_cubic(&p0, &p1, &p2, &p3, input.t);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(Point3DOutput {
                x: result.x,
                y: result.y,
                z: result.z,
            }),
            err: None::<String>,
        })
        .unwrap(),
    )
}

// View matrices
#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_look_at_matrix_json(
    input: *const c_char,
) -> *mut c_char {
    let input: LookAtInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let eye = computer_graphics::Vector3D::new(input.eye.x, input.eye.y, input.eye.z);
    let center = computer_graphics::Vector3D::new(input.center.x, input.center.y, input.center.z);
    let up = computer_graphics::Vector3D::new(input.up.x, input.up.y, input.up.z);
    let matrix = computer_graphics::look_at_matrix(&eye, &center, &up);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]
pub unsafe extern "C" fn rssn_num_graphics_perspective_matrix_json(
    input: *const c_char,
) -> *mut c_char {
    let input: PerspectiveInput = match from_json_string(input) {
        Some(i) => i,
        None => {
            return to_c_string(
                serde_json::to_string(&FfiResult::<Vec<f64>, String> {
                    ok: None,
                    err: Some("Invalid JSON".to_string()),
                })
                .unwrap(),
            )
        }
    };
    let matrix =
        computer_graphics::perspective_matrix(input.fov_y, input.aspect, input.near, input.far);
    to_c_string(
        serde_json::to_string(&FfiResult {
            ok: Some(matrix.data()),
            err: None::<String>,
        })
        .unwrap(),
    )
}
