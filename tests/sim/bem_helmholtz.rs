//! `sim::bem_helmholtz`: Burton–Miller exterior scattering. The Mie
//! comparisons live in `bem_research.rs`; these tests check properties that
//! hold independently of the series: input validation, invariance under
//! translation of the body (a pure phase of the incident wave) and under the
//! scaling of the direction, the extinction theorem inside the body, the
//! Rayleigh (low-frequency) limits of the soft and hard sphere, and the
//! agreement of the Burton–Miller and conventional formulations away from
//! the irregular frequencies.

use num_complex::Complex64;
use rssn::sim::bem_helmholtz::PlaneWave;
use rssn::sim::bem_helmholtz::Scatterer;
use rssn::sim::bem_helmholtz::solve_helmholtz_exterior;
use rssn::sim::physics_bem::SurfaceMesh3D;

const Z: PlaneWave = PlaneWave {
    direction: [0.0, 0.0, 1.0],
};

fn unit_sphere(sub: usize) -> SurfaceMesh3D {
    SurfaceMesh3D::icosphere(1.0, [0.0; 3], sub)
}

#[test]
fn invalid_input_is_rejected() {
    let mesh = unit_sphere(1);
    for k in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            solve_helmholtz_exterior(&mesh, k, Scatterer::SoundSoft, Z, None).is_err(),
            "k = {k}"
        );
    }
    let still = PlaneWave {
        direction: [0.0, 0.0, 0.0],
    };
    assert!(solve_helmholtz_exterior(&mesh, 1.0, Scatterer::SoundHard, still, None).is_err());
    // One unknown per triangle.
    let sol = solve_helmholtz_exterior(&mesh, 1.0, Scatterer::SoundHard, Z, None).unwrap();
    assert_eq!(sol.density.len(), mesh.triangles.len());
    assert_eq!(sol.k, 1.0);
}

#[test]
fn the_incident_field_is_the_plane_wave() {
    let mesh = unit_sphere(1);
    let k = 2.0;
    let sol = solve_helmholtz_exterior(
        &mesh,
        k,
        Scatterer::SoundSoft,
        PlaneWave {
            direction: [0.0, 3.0, 4.0],
        },
        None,
    )
    .unwrap();
    // e^{i k d.x} with d normalised: d = (0, 0.6, 0.8).
    for x in [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [-0.5, 0.7, 0.2]] {
        let phase = k * (0.6 * x[1] + 0.8 * x[2]);
        let want = Complex64::new(phase.cos(), phase.sin());
        assert!((sol.incident_field(x) - want).norm() < 1e-14, "{x:?}");
    }
}

#[test]
fn the_direction_is_normalised_and_translation_is_a_phase() {
    let k = 1.5;
    let mesh = unit_sphere(2);
    let a = solve_helmholtz_exterior(&mesh, k, Scatterer::SoundHard, Z, None).unwrap();
    let b = solve_helmholtz_exterior(
        &mesh,
        k,
        Scatterer::SoundHard,
        PlaneWave {
            direction: [0.0, 0.0, 7.0],
        },
        None,
    )
    .unwrap();
    let diff = a
        .density
        .iter()
        .zip(&b.density)
        .map(|(x, y)| (x - y).norm())
        .fold(0.0, f64::max);
    assert!(
        diff < 1e-12,
        "direction scaling changes the solution by {diff}"
    );
    // Moving the body by c multiplies the incident wave at the body by
    // e^{i k d.c}; the scattered field at x + c is that phase times the
    // original scattered field at x.
    let shift = [0.4, -0.3, 0.6];
    let moved = SurfaceMesh3D::icosphere(1.0, shift, 2);
    let c = solve_helmholtz_exterior(&moved, k, Scatterer::SoundHard, Z, None).unwrap();
    let phase = Complex64::new(0.0, k * shift[2]).exp();
    for x in [[2.0, 0.0, 0.0], [0.0, 1.5, 1.5], [-1.0, -1.0, -2.0]] {
        let original = a.scattered_field(&mesh, x);
        let translated =
            c.scattered_field(&moved, [x[0] + shift[0], x[1] + shift[1], x[2] + shift[2]]);
        assert!(
            (translated - phase * original).norm() < 1e-9 * original.norm().max(1e-3),
            "{x:?}: {translated} vs {}",
            phase * original
        );
    }
}

#[test]
fn the_scattered_field_extinguishes_the_incident_wave_inside_the_body() {
    // Inside the scatterer the representation formula gives -p_inc (the
    // extinction theorem), for both boundary conditions.
    let k = 1.0;
    let mesh = unit_sphere(2);
    for kind in [Scatterer::SoundSoft, Scatterer::SoundHard] {
        let sol = solve_helmholtz_exterior(&mesh, k, kind, Z, None).unwrap();
        for x in [[0.0, 0.0, 0.0], [0.3, 0.0, 0.2], [-0.2, 0.3, -0.4]] {
            let scattered = sol.scattered_field(&mesh, x);
            let incident = sol.incident_field(x);
            assert!(
                (scattered + incident).norm() < 0.12,
                "{kind:?} at {x:?}: {scattered} + {incident}"
            );
        }
    }
}

#[test]
fn the_soft_sphere_scatters_a_monopole_at_low_frequency() {
    // For ka << 1 a sound-soft sphere radiates the monopole
    // p_s = -a e^{ika} / (1 + ika) e^{ik(r - a)}/r ... the leading term is
    // -a/(1 + ika) e^{ikr}/r.
    let k = 0.1;
    let mesh = unit_sphere(3);
    // The conventional equation is regular at low frequency (there is no
    // irregular frequency below ka = pi); the Burton-Miller coupling i/k is
    // large there and amplifies the discretisation error of its double-layer
    // terms into a phase error of order 2 ka.
    let sol = solve_helmholtz_exterior(
        &mesh,
        k,
        Scatterer::SoundSoft,
        Z,
        Some(Complex64::new(0.0, 0.0)),
    )
    .unwrap();
    for x in [
        [3.0_f64, 0.0, 0.0],
        [0.0, 0.0, 4.0],
        [0.0, -3.5, 0.0],
        [0.0, 0.0, -3.0],
    ] {
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let want = -Complex64::new(1.0, 0.0) / Complex64::new(1.0, k)
            * Complex64::new(0.0, k * r).exp()
            / r;
        let got = sol.scattered_field(&mesh, x);
        assert!(
            (got - want).norm() < 0.10 * want.norm(),
            "{x:?}: {got} vs {want}"
        );
    }
}

#[test]
fn the_hard_sphere_surface_pressure_matches_the_rayleigh_expansion() {
    // Total pressure on a rigid sphere for ka << 1:
    // p = 1 + (3/2) i ka cos(theta) + O((ka)^2), theta from the incidence axis.
    let k = 0.1;
    let mesh = unit_sphere(3);
    let sol = solve_helmholtz_exterior(&mesh, k, Scatterer::SoundHard, Z, None).unwrap();
    // Least-squares slope of Im p against cos(theta) and mean of Re p.
    let (mut sxx, mut sxy, mut re_sum) = (0.0, 0.0, 0.0);
    let nt = mesh.triangles.len();
    for t in 0..nt {
        let c = mesh.centroid(t).unwrap();
        let cos_theta = c[2] / (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        sxx += cos_theta * cos_theta;
        sxy += cos_theta * sol.density[t].im;
        re_sum += sol.density[t].re;
    }
    let slope = sxy / sxx;
    assert!(
        (slope - 1.5 * k).abs() < 0.2 * 1.5 * k,
        "slope {slope} vs {}",
        1.5 * k
    );
    assert!(
        (re_sum / nt as f64 - 1.0).abs() < 0.05,
        "mean Re p {}",
        re_sum / nt as f64
    );
}

#[test]
fn burton_miller_and_the_conventional_equation_agree_away_from_irregular_frequencies() {
    // At ka = 1 (no interior eigenfrequency of the unit sphere near it) the
    // coupling term only perturbs the discretisation error: the scattered
    // fields of alpha = i/k and of the conventional equation (alpha = 0) agree
    // to a few per cent, while the density unknowns are of order one.
    let k = 1.0;
    let mesh = unit_sphere(2);
    for kind in [Scatterer::SoundSoft, Scatterer::SoundHard] {
        let burton_miller = solve_helmholtz_exterior(&mesh, k, kind, Z, None).unwrap();
        let conventional =
            solve_helmholtz_exterior(&mesh, k, kind, Z, Some(Complex64::new(0.0, 0.0))).unwrap();
        let mut scale: f64 = 0.0;
        let mut worst: f64 = 0.0;
        for x in [
            [2.5, 0.0, 0.0],
            [0.0, 0.0, 2.5],
            [0.0, 0.0, -2.5],
            [1.5, 1.5, 1.0],
        ] {
            let (a, b) = (
                burton_miller.scattered_field(&mesh, x),
                conventional.scattered_field(&mesh, x),
            );
            scale = scale.max(a.norm());
            worst = worst.max((a - b).norm());
        }
        assert!(scale > 0.05, "{kind:?}: the body scatters, |p_s| = {scale}");
        assert!(
            worst < 0.10 * scale,
            "{kind:?}: formulations differ by {worst} of {scale}"
        );
    }
}
