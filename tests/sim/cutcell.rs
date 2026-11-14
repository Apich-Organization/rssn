//! `sim::cutcell`: open fractions of segments, faces and cells cut by a
//! level-set body. Planar interfaces are reproduced to round-off (the clipped
//! piecewise-linear interpolant is exact for linear level sets); curved
//! bodies converge to the analytic area and volume as the refinement depth
//! grows.

use rssn::sim::cutcell::default_depth;
use rssn::sim::cutcell::open_fraction;
use rssn::sim::cutcell::open_length_fraction;
use rssn::sim::cutcell::simplex_positive_fraction;
use std::f64::consts::PI;

const UNIT: [f64; 3] = [1.0, 1.0, 1.0];
const ORIGIN: [f64; 3] = [0.0, 0.0, 0.0];

#[test]
fn simplex_fractions_have_the_closed_forms() {
    // Segment [f0, f1]: the positive part is f0/(f0 - f1) for f0 > 0 > f1.
    assert!((simplex_positive_fraction(&[1.0, -1.0]) - 0.5).abs() < 1e-15);
    assert!((simplex_positive_fraction(&[2.0, -1.0]) - 2.0 / 3.0).abs() < 1e-15);
    assert!((simplex_positive_fraction(&[-1.0, 3.0]) - 0.75).abs() < 1e-15);
    // Triangle: one positive vertex cuts off (f/(f - f_j))(f/(f - f_k)), two
    // positive vertices leave the complement.
    assert!((simplex_positive_fraction(&[1.0, -1.0, -1.0]) - 0.25).abs() < 1e-15);
    assert!((simplex_positive_fraction(&[1.0, 1.0, -1.0]) - 0.75).abs() < 1e-15);
    assert!(
        (simplex_positive_fraction(&[3.0, -1.0, -2.0]) - (3.0 / 4.0) * (3.0 / 5.0)).abs() < 1e-15
    );
    // Tetrahedron: one positive vertex (3/4)³, three positive vertices the
    // complement, two and two by the symmetry f -> -f.
    assert!((simplex_positive_fraction(&[3.0, -1.0, -1.0, -1.0]) - 27.0 / 64.0).abs() < 1e-14);
    assert!(
        (simplex_positive_fraction(&[1.0, 1.0, 1.0, -3.0]) - (1.0 - 27.0 / 64.0)).abs() < 1e-14
    );
    assert!((simplex_positive_fraction(&[1.0, 1.0, -1.0, -1.0]) - 0.5).abs() < 1e-14);
    let skew = simplex_positive_fraction(&[2.0, 1.0, -1.0, -3.0]);
    let mirror = simplex_positive_fraction(&[-2.0, -1.0, 1.0, 3.0]);
    assert!(
        (skew + mirror - 1.0).abs() < 1e-13,
        "complementary: {skew} + {mirror}"
    );
    assert!(skew > 0.0 && skew < 1.0);
    // Entirely open and entirely closed.
    assert_eq!(simplex_positive_fraction(&[0.5, 2.0, 1.0]), 1.0);
    assert_eq!(simplex_positive_fraction(&[-0.5, -2.0, 0.0]), 0.0);
    // A vertex at zero counts as closed.
    assert!((simplex_positive_fraction(&[1.0, 0.0, 0.0]) - 1.0).abs() < 1e-15);
}

#[test]
fn planar_interfaces_are_exact_in_one_two_and_three_dimensions() {
    let x_cut = |p: [f64; 3]| p[0] - 0.3;
    // Segment, face, cell: the open part is x > 0.3.
    assert!((open_fraction(&x_cut, ORIGIN, UNIT, &[0], 0) - 0.7).abs() < 1e-12);
    assert!((open_fraction(&x_cut, ORIGIN, UNIT, &[0, 1], 3) - 0.7).abs() < 1e-12);
    assert!((open_fraction(&x_cut, ORIGIN, UNIT, &[0, 1, 2], 2) - 0.7).abs() < 1e-12);
    // Oblique planes: the triangle x + y > 0.5 of the square and the corner
    // simplex cut from the cube.
    let diagonal = |p: [f64; 3]| p[0] + p[1] - 0.5;
    assert!((open_fraction(&diagonal, ORIGIN, UNIT, &[0, 1], 3) - 0.875).abs() < 1e-12);
    for (c, want) in [
        (1.0_f64, 5.0 / 6.0),
        (1.5, 0.5),
        (0.5, 47.0 / 48.0),
        (2.5, 1.0 / 48.0),
    ] {
        let plane = move |p: [f64; 3]| p[0] + p[1] + p[2] - c;
        let got = open_fraction(&plane, ORIGIN, UNIT, &[0, 1, 2], 1);
        assert!(
            (got - want).abs() < 1e-12,
            "x + y + z > {c}: {got} vs {want}"
        );
    }
    // An anisotropic box: lengths scale the sub-box but not the fraction.
    let skew = |p: [f64; 3]| 2.0 * p[0] - p[1] - 0.6;
    let wide = open_fraction(&skew, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], &[0, 1], 3);
    // Area of {2x - y > 0.6} in the unit square: for x in [0.3, 0.8] the
    // open height is 2x - 0.6 (integral 0.25), for x in [0.8, 1] it is 1.
    let reference = 0.25 + 0.2;
    assert!((wide - reference).abs() < 1e-12, "{wide} vs {reference}");
    // A face of a 3D cell spans two of the three axes and ignores the third.
    let xz = |p: [f64; 3]| p[0] + p[2] - 1.0;
    let face = open_fraction(&xz, [0.0, 5.0, 0.0], UNIT, &[0, 2], 3);
    assert!((face - 0.5).abs() < 1e-12, "{face}");
    // Shifted and scaled boxes.
    let shifted = open_fraction(&x_cut, [0.2, 0.0, 0.0], [0.4, 2.0, 1.0], &[0, 1], 3);
    // x in [0.2, 0.6], open for x > 0.3: 0.3/0.4.
    assert!((shifted - 0.75).abs() < 1e-12, "{shifted}");
}

#[test]
fn trivial_and_degenerate_boxes() {
    let open = |_: [f64; 3]| 1.0;
    let closed = |_: [f64; 3]| -1.0;
    for axes in [vec![0], vec![0, 1], vec![0, 1, 2]] {
        assert_eq!(open_fraction(&open, ORIGIN, UNIT, &axes, 3), 1.0);
        assert_eq!(open_fraction(&closed, ORIGIN, UNIT, &axes, 3), 0.0);
    }
    // No axes: a point, open or closed.
    assert_eq!(open_fraction(&open, [0.3, 0.4, 0.5], UNIT, &[], 0), 1.0);
    assert_eq!(open_fraction(&closed, [0.3, 0.4, 0.5], UNIT, &[], 0), 0.0);
    // The level set is positive in the open region and closed where it
    // vanishes: a function that is zero everywhere is a body.
    assert_eq!(
        open_fraction(&|_: [f64; 3]| 0.0, ORIGIN, UNIT, &[0, 1], 2),
        0.0
    );
    // Defaults: more refinement in 1D/2D than in 3D.
    assert_eq!(default_depth(1), 4);
    assert_eq!(default_depth(2), 4);
    assert_eq!(default_depth(3), 3);
}

#[test]
fn segments_cut_by_a_circle_use_the_exact_chord() {
    let r = 0.3;
    let disk = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - r;
    // The line y = 0.6 meets the circle in a chord of length 2 sqrt(r² - 0.1²).
    let chord = 2.0 * (r * r - 0.01_f64).sqrt();
    let got = open_length_fraction(&disk, [0.0, 0.6, 0.0], 0, 1.0);
    assert!(
        (got - (1.0 - chord)).abs() < 1e-12,
        "{got} vs {}",
        1.0 - chord
    );
    // Along y the same chord by symmetry, and a tangent-free miss is fully open.
    let got_y = open_length_fraction(&disk, [0.6, 0.0, 0.0], 1, 1.0);
    assert!((got_y - (1.0 - chord)).abs() < 1e-12);
    assert_eq!(open_length_fraction(&disk, [0.0, 0.95, 0.0], 0, 1.0), 1.0);
    // A segment through the centre loses the diameter; a segment that ends
    // inside the body loses the part inside.
    let diameter = open_length_fraction(&disk, [0.0, 0.5, 0.0], 0, 1.0);
    assert!((diameter - (1.0 - 2.0 * r)).abs() < 1e-12);
    let half = open_length_fraction(&disk, [0.0, 0.5, 0.0], 0, 0.5);
    // x in [0, 0.5]: the body occupies [0.2, 0.5].
    assert!((half - 0.4).abs() < 1e-12, "{half}");
    // Two separate bodies on one segment (sign changes at four places).
    let pair = |p: [f64; 3]| ((p[0] - 0.25).abs() - 0.05).min((p[0] - 0.75).abs() - 0.1);
    let got = open_length_fraction(&pair, ORIGIN, 0, 1.0);
    assert!((got - (1.0 - 0.1 - 0.2)).abs() < 1e-12, "{got}");
}

#[test]
fn curved_bodies_converge_to_the_analytic_area_and_volume() {
    let r = 0.3;
    let disk = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - r;
    let exact = 1.0 - PI * r * r;
    let errors: Vec<f64> = [1_usize, 3, 5]
        .iter()
        .map(|&depth| (open_fraction(&disk, ORIGIN, UNIT, &[0, 1], depth) - exact).abs())
        .collect();
    assert!(
        errors[2] < errors[1] && errors[1] < errors[0],
        "errors {errors:?} do not shrink"
    );
    assert!(errors[2] < 3e-3, "depth 5 error {}", errors[2]);
    assert!((open_fraction(&disk, ORIGIN, UNIT, &[0, 1], default_depth(2)) - exact).abs() < 6e-3);
    // The body can sit strictly inside a box whose corners are all open: the
    // centre probe must still refine.
    let small = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - 0.1;
    let got = open_fraction(&small, ORIGIN, UNIT, &[0, 1], 4);
    assert!((got - (1.0 - PI * 0.01)).abs() < 3e-3, "{got}");
    // A sphere in a cell: 1 - 4/3 π r³.
    let ball = |p: [f64; 3]| {
        ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2) + (p[2] - 0.5).powi(2)).sqrt() - 0.35
    };
    let exact = 1.0 - 4.0 / 3.0 * PI * 0.35_f64.powi(3);
    let coarse = open_fraction(&ball, ORIGIN, UNIT, &[0, 1, 2], 2);
    let fine = open_fraction(&ball, ORIGIN, UNIT, &[0, 1, 2], 4);
    assert!(
        (fine - exact).abs() < (coarse - exact).abs(),
        "{fine} vs {coarse} (exact {exact})"
    );
    assert!((fine - exact).abs() < 1.2e-2, "sphere: {fine} vs {exact}");
    // A circle cutting through one face of a 3D cell: a cylinder along z,
    // whose open face in the xz-plane at y = 0.5 loses the strip of the diameter.
    let cylinder = |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - 0.25;
    let face = open_fraction(&cylinder, [0.0, 0.5, 0.0], UNIT, &[0, 2], 3);
    assert!((face - 0.5).abs() < 1e-3, "{face}");
}

#[test]
fn open_and_closed_fractions_are_complementary_and_additive() {
    let r = 0.31;
    let body = |p: [f64; 3]| ((p[0] - 0.47).powi(2) + (p[1] - 0.52).powi(2)).sqrt() - r;
    let anti = |p: [f64; 3]| -body(p);
    let (open, closed) = (
        open_fraction(&body, ORIGIN, UNIT, &[0, 1], 4),
        open_fraction(&anti, ORIGIN, UNIT, &[0, 1], 4),
    );
    assert!((open + closed - 1.0).abs() < 1e-12, "{open} + {closed}");
    // The fraction of a box is the mean of the fractions of its quarters when
    // the interface is planar (the refinement tree is then irrelevant).
    let plane = |p: [f64; 3]| 3.0 * p[0] - 2.0 * p[1] - 0.7;
    let whole = open_fraction(&plane, ORIGIN, UNIT, &[0, 1], 2);
    let mut quarters = 0.0;
    for (i, j) in [(0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5)] {
        quarters += open_fraction(&plane, [i, j, 0.0], [0.5, 0.5, 1.0], &[0, 1], 2);
    }
    assert!(
        (whole - quarters / 4.0).abs() < 1e-12,
        "{whole} vs {}",
        quarters / 4.0
    );
    // Monotone in the level set: enlarging the body shrinks the open part.
    let mut previous = 1.0;
    for radius in [0.1, 0.2, 0.3, 0.4] {
        let phi = move |p: [f64; 3]| ((p[0] - 0.5).powi(2) + (p[1] - 0.5).powi(2)).sqrt() - radius;
        let f = open_fraction(&phi, ORIGIN, UNIT, &[0, 1], 4);
        assert!(f < previous, "radius {radius}: {f} !< {previous}");
        previous = f;
    }
}
