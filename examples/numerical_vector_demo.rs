use rssn::prelude::numerical::*;
use rssn::prelude::*;

fn main() {

    println!(
        "--- Numerical Vector \
         Operations Demo ---"
    );

    // Define some vectors
    let v1 = vec![1.0, 2.0, 3.0];

    let v2 = vec![4.0, 5.0, 6.0];

    println!("v1: {:?}", v1);

    println!("v2: {:?}", v2);

    // Basic Arithmetic
    let sum =
        numerical_vec_add(&v1, &v2)
            .unwrap();

    let diff =
        numerical_vec_sub(&v2, &v1)
            .unwrap();

    let scaled =
        numerical_scalar_mul(&v1, 1.5);

    println!("v1 + v2: {:?}", sum);

    println!(
        "v2 - v1: {:?}",
        diff
    );

    println!(
        "1.5 * v1: {:?}",
        scaled
    );

    // Products and Norms
    let dot =
        numerical_dot_product(&v1, &v2)
            .unwrap();

    let norm_v1 = numerical_norm(&v1);

    let dist =
        numerical_distance(&v1, &v2)
            .unwrap();

    println!("v1 . v2: {}", dot);

    println!(
        "||v1||: {}",
        norm_v1
    );

    println!(
        "dist(v1, v2): {}",
        dist
    );

    // 3D Specific
    if let Ok(cross) =
        numerical_cross_product(
            &v1, &v2,
        )
    {

        println!(
            "v1 x v2: {:?}",
            cross
        );
    }

    // Geometry
    let unit_v1 =
        numerical_normalize(&v1)
            .unwrap();

    let angle_rad =
        numerical_angle(&v1, &v2)
            .unwrap();

    let angle_deg =
        angle_rad.to_degrees();

    println!(
        "Normalized v1: {:?}",
        unit_v1
    );

    println!(
        "Angle between v1 and v2: {} \
         rad ({} deg)",
        angle_rad, angle_deg
    );

    // Projections and Reflections
    let v_proj = vec![1.0, 1.0, 0.0];

    let v_onto = vec![1.0, 0.0, 0.0];

    let proj = numerical_project(
        &v_proj, &v_onto,
    )
    .unwrap();

    println!(
        "Project {:?} onto {:?}: {:?}",
        v_proj, v_onto, proj
    );

    let v_refl = vec![1.0, -1.0, 0.0];

    let normal = vec![0.0, 1.0, 0.0];

    let refl = numerical_reflect(
        &v_refl, &normal,
    )
    .unwrap();

    println!(
        "Reflect {:?} about normal \
         {:?}: {:?}",
        v_refl, normal, refl
    );

    // Interpolation
    let t = 0.5;

    let midpoint =
        numerical_lerp(&v1, &v2, t)
            .unwrap();

    println!(
        "Lerp between v1 and v2 at \
         t={}: {:?}",
        t, midpoint
    );
}
