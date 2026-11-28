//! Area invariants for finite binary64 polygons across coordinate scales.

use rssn::kernels::geometry::polygon_area;

fn rectangle(
    x: f64,
    y: f64,
    width: f64,
    height: f64,
) -> [[f64; 2]; 4] {
    [
        [x, y],
        [x + width, y],
        [x + width, y + height],
        [x, y + height],
    ]
}

#[test]
fn polygon_area_survives_translation_and_extreme_products() {
    for (mut points, expected) in [
        (rectangle(2f64.powi(54), 2f64.powi(54), 4.0, 4.0), 16.0),
        (
            rectangle(
                2f64.powi(1020),
                2f64.powi(40),
                2f64.powi(970),
                2f64.powi(-10),
            ),
            2f64.powi(960),
        ),
        (
            rectangle(0.0, 0.0, 2f64.powi(512), 2f64.powi(511)),
            2f64.powi(1023),
        ),
        (
            rectangle(0.0, 0.0, 2f64.powi(-537), 2f64.powi(-537)),
            f64::from_bits(1),
        ),
        (
            [
                [-1e308, 0.0],
                [1e308, 0.0],
                [1e308, 1e-300],
                [-1e308, 1e-300],
            ],
            2e8,
        ),
    ] {
        for _ in 0..4 {
            assert_eq!(polygon_area(points.as_flattened()), expected);
            points.reverse();
            assert_eq!(polygon_area(points.as_flattened()), expected);
            points.reverse();
            points.rotate_left(1);
        }
    }
}

#[test]
fn polygon_area_preserves_thin_areas_and_flat_input_conventions() {
    let m = 134_217_728.0;
    assert_eq!(polygon_area(&[0.0, 0.0, m, m - 1.0, m + 1.0, m]), 0.5);
    assert_eq!(polygon_area(&[0.0, 0.0, 2.0, 2.0, 3.0, 3.0]), 0.0);
    assert_eq!(polygon_area(&[0.0, 0.0, 2.0, 0.0, 0.0, 3.0, f64::NAN]), 3.0);
    assert_eq!(polygon_area(&[0.0, 0.0, 2.0, 3.0]), 0.0);
    assert_eq!(polygon_area(&[]), 0.0);
}
