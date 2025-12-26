//! Handle-based FFI API for numerical computer graphics functions.

use crate::numerical::computer_graphics;

/// Computes the dot product of two 3D vectors.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_dot_product(
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
) -> f64 {

    let v1 = computer_graphics::Vector3D::new(x1, y1, z1);

    let v2 = computer_graphics::Vector3D::new(x2, y2, z2);

    computer_graphics::dot_product(
        &v1, &v2,
    )
}

/// Computes the cross product of two 3D vectors.
/// Result is stored in out_x, out_y, out_z.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_cross_product(
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> i32 {

    if out_x.is_null()
        || out_y.is_null()
        || out_z.is_null()
    {

        return -1;
    }

    let v1 = computer_graphics::Vector3D::new(x1, y1, z1);

    let v2 = computer_graphics::Vector3D::new(x2, y2, z2);

    let result = computer_graphics::cross_product(&v1, &v2);

    *out_x = result.x;

    *out_y = result.y;

    *out_z = result.z;

    0
}

/// Normalizes a 3D vector.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_normalize(
    x: f64,
    y: f64,
    z: f64,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> i32 {

    if out_x.is_null()
        || out_y.is_null()
        || out_z.is_null()
    {

        return -1;
    }

    let v = computer_graphics::Vector3D::new(x, y, z);

    let result = v.normalize();

    *out_x = result.x;

    *out_y = result.y;

    *out_z = result.z;

    0
}

/// Computes the magnitude of a 3D vector.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_magnitude(
    x: f64,
    y: f64,
    z: f64,
) -> f64 {

    let v = computer_graphics::Vector3D::new(x, y, z);

    v.magnitude()
}

/// Computes the reflection vector.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_reflect(
    ix: f64,
    iy: f64,
    iz: f64,
    nx: f64,
    ny: f64,
    nz: f64,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> i32 {

    if out_x.is_null()
        || out_y.is_null()
        || out_z.is_null()
    {

        return -1;
    }

    let incident =
        computer_graphics::Vector3D::new(ix, iy, iz);

    let normal =
        computer_graphics::Vector3D::new(nx, ny, nz);

    let result =
        computer_graphics::reflect(
            &incident,
            &normal,
        );

    *out_x = result.x;

    *out_y = result.y;

    *out_z = result.z;

    0
}

/// Computes the angle between two 3D vectors in radians.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_angle_between(
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
) -> f64 {

    let v1 = computer_graphics::Vector3D::new(x1, y1, z1);

    let v2 = computer_graphics::Vector3D::new(x2, y2, z2);

    computer_graphics::angle_between(
        &v1, &v2,
    )
}

/// Converts degrees to radians.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_degrees_to_radians(
    degrees: f64
) -> f64 {

    computer_graphics::degrees_to_radians(degrees)
}

/// Converts radians to degrees.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_radians_to_degrees(
    radians: f64
) -> f64 {

    computer_graphics::radians_to_degrees(radians)
}

/// Quaternion multiply.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_quaternion_multiply(
    w1: f64,
    x1: f64,
    y1: f64,
    z1: f64,
    w2: f64,
    x2: f64,
    y2: f64,
    z2: f64,
    out_w: *mut f64,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> i32 {

    if out_w.is_null()
        || out_x.is_null()
        || out_y.is_null()
        || out_z.is_null()
    {

        return -1;
    }

    let q1 =
        computer_graphics::Quaternion::new(w1, x1, y1, z1);

    let q2 =
        computer_graphics::Quaternion::new(w2, x2, y2, z2);

    let result = q1.multiply(&q2);

    *out_w = result.w;

    *out_x = result.x;

    *out_y = result.y;

    *out_z = result.z;

    0
}

/// Rotation matrix around X axis.
/// Output: 16 f64 values in row-major order.
///
/// # Safety
/// `out_ptr` must point to at least 16 f64 values.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_x(
    angle_rad: f64,
    out_ptr: *mut f64,
) -> i32 {

    if out_ptr.is_null() {

        return -1;
    }

    let matrix =
        computer_graphics::rotation_matrix_x(angle_rad);

    for i in 0..16 {

        *out_ptr.add(i) =
            matrix.data()[i];
    }

    0
}

/// Ray-sphere intersection.
/// Returns t value or -1 if no intersection.
#[no_mangle]

pub extern "C" fn rssn_num_graphics_ray_sphere_intersection(
    ray_ox: f64,
    ray_oy: f64,
    ray_oz: f64,
    ray_dx: f64,
    ray_dy: f64,
    ray_dz: f64,
    sphere_cx: f64,
    sphere_cy: f64,
    sphere_cz: f64,
    sphere_r: f64,
) -> f64 {

    let ray = computer_graphics::Ray::new(
        computer_graphics::Point3D::new(
            ray_ox, ray_oy, ray_oz,
        ),
        computer_graphics::Vector3D::new(
            ray_dx, ray_dy, ray_dz,
        ),
    );

    let sphere = computer_graphics::Sphere::new(
        computer_graphics::Point3D::new(
            sphere_cx, sphere_cy, sphere_cz,
        ),
        sphere_r,
    );

    match computer_graphics::ray_sphere_intersection(
        &ray, &sphere,
    ) {
        | Some(intersection) => intersection.t,
        | None => -1.0,
    }
}

/// Bezier cubic curve evaluation.
///
/// # Safety
/// Pointers must be valid.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_bezier_cubic(
    p0x: f64,
    p0y: f64,
    p0z: f64,
    p1x: f64,
    p1y: f64,
    p1z: f64,
    p2x: f64,
    p2y: f64,
    p2z: f64,
    p3x: f64,
    p3y: f64,
    p3z: f64,
    t: f64,
    out_x: *mut f64,
    out_y: *mut f64,
    out_z: *mut f64,
) -> i32 {

    if out_x.is_null()
        || out_y.is_null()
        || out_z.is_null()
    {

        return -1;
    }

    let p0 =
        computer_graphics::Point3D::new(
            p0x, p0y, p0z,
        );

    let p1 =
        computer_graphics::Point3D::new(
            p1x, p1y, p1z,
        );

    let p2 =
        computer_graphics::Point3D::new(
            p2x, p2y, p2z,
        );

    let p3 =
        computer_graphics::Point3D::new(
            p3x, p3y, p3z,
        );

    let result =
        computer_graphics::bezier_cubic(
            &p0, &p1, &p2, &p3, t,
        );

    *out_x = result.x;

    *out_y = result.y;

    *out_z = result.z;

    0
}
