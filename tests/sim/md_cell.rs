//! `sim::md_cell`: the small-matrix helpers behind the flexible-cell
//! barostat (product, determinant, inverse, `e^A` and `phi1(A)`) against
//! exact identities, and the barostat's constructors, degrees of freedom and
//! conserved-quantity contribution. The coupled dynamics (NPH/NPT, SLLOD,
//! Lees-Edwards) are in `md_engine.rs`.

use rssn::kernels::random::Rng;
// `sim::md_cell` is a private submodule of `md_engine`, which re-exports it.
use rssn::sim::md_engine::CellMode;
use rssn::sim::md_engine::FlexibleBarostat;
use rssn::sim::md_engine::mat_det;
use rssn::sim::md_engine::mat_exp_phi1;
use rssn::sim::md_engine::mat_inverse;
use rssn::sim::md_engine::mat_mul;
use rssn::sim::md_engine::mat_vec;

fn random_matrix<const D: usize>(
    rng: &mut Rng,
    scale: f64,
) -> [[f64; D]; D] {
    let mut m = [[0.0; D]; D];
    for row in &mut m {
        for x in row.iter_mut() {
            *x = scale * rng.normal();
        }
    }
    m
}

fn max_abs_diff<const D: usize>(
    a: &[[f64; D]; D],
    b: &[[f64; D]; D],
) -> f64 {
    let mut worst: f64 = 0.0;
    for i in 0..D {
        for j in 0..D {
            worst = worst.max((a[i][j] - b[i][j]).abs());
        }
    }
    worst
}

fn identity<const D: usize>() -> [[f64; D]; D] {
    let mut m = [[0.0; D]; D];
    for i in 0..D {
        m[i][i] = 1.0;
    }
    m
}

#[test]
fn products_and_determinants_follow_the_definitions() {
    let a = [[1.0, 2.0, 0.5], [-1.0, 0.0, 3.0], [2.0, 1.0, 1.0]];
    let b = [[0.0, 1.0, 0.0], [1.0, 0.0, 2.0], [0.5, -1.0, 1.0]];
    let c = mat_mul(&a, &b);
    for i in 0..3 {
        for j in 0..3 {
            let want: f64 = (0..3).map(|k| a[i][k] * b[k][j]).sum();
            assert!((c[i][j] - want).abs() < 1e-14);
        }
    }
    assert_eq!(mat_mul(&a, &identity()), a);
    let y = mat_vec(&a, &[1.0, -2.0, 4.0]);
    assert_eq!(y, [1.0 - 4.0 + 2.0, -1.0 + 12.0, 2.0 - 2.0 + 4.0]);
    // Determinants: the rule of Sarrus, the permutation sign, singular matrices
    // and the product rule.
    let sarrus = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    assert!(
        (mat_det(&a) - sarrus).abs() < 1e-13,
        "{} vs {sarrus}",
        mat_det(&a)
    );
    assert!(
        (mat_det(&[[0.0, 1.0], [1.0, 0.0]]) + 1.0).abs() < 1e-15,
        "a swap has determinant -1"
    );
    assert_eq!(mat_det(&[[1.0, 2.0], [2.0, 4.0]]), 0.0);
    assert_eq!(mat_det(&[[0.0, 0.0], [1.0, 3.0]]), 0.0);
    assert!((mat_det(&identity::<4>()) - 1.0).abs() < 1e-15);
    assert!((mat_det(&c) - mat_det(&a) * mat_det(&b)).abs() < 1e-12);
    // A triclinic cell: lower-triangular rows, determinant the product of the diagonal.
    let cell = [[3.0, 0.0, 0.0], [1.0, 4.0, 0.0], [0.5, -2.0, 5.0]];
    assert!((mat_det(&cell) - 60.0).abs() < 1e-12);
    assert!((mat_det(&[[2.5]]) - 2.5).abs() < 1e-15);
}

#[test]
fn inverses_invert_and_singular_matrices_are_refused() {
    let mut rng = Rng::new(3);
    for _ in 0..20 {
        let a: [[f64; 3]; 3] = random_matrix(&mut rng, 1.0);
        if mat_det(&a).abs() < 1e-2 {
            continue;
        }
        let inv = mat_inverse(&a).unwrap();
        assert!(max_abs_diff(&mat_mul(&a, &inv), &identity()) < 1e-10);
        assert!(max_abs_diff(&mat_mul(&inv, &a), &identity()) < 1e-10);
        assert!((mat_det(&inv) * mat_det(&a) - 1.0).abs() < 1e-10);
    }
    // A matrix that needs a row exchange.
    let p = [[0.0, 2.0], [3.0, 1.0]];
    let inv = mat_inverse(&p).unwrap();
    assert!(
        max_abs_diff(&inv, &[[-1.0 / 6.0, 1.0 / 3.0], [0.5, 0.0]]) < 1e-14,
        "{inv:?}"
    );
    assert!(mat_inverse(&[[1.0, 2.0], [2.0, 4.0]]).is_none());
    assert!(mat_inverse(&[[0.0, 0.0], [0.0, 0.0]]).is_none());
    // Four dimensions as well.
    let a: [[f64; 4]; 4] = random_matrix(&mut rng, 1.0);
    let inv = mat_inverse(&a).unwrap();
    assert!(max_abs_diff(&mat_mul(&a, &inv), &identity()) < 1e-8);
}

#[test]
fn exponential_and_phi1_satisfy_their_defining_identities() {
    // Nilpotent: e^A = I + A, phi1(A) = I + A/2.
    let n = [[0.0, 1.0], [0.0, 0.0]];
    let (e, phi) = mat_exp_phi1(&n);
    assert!(max_abs_diff(&e, &[[1.0, 1.0], [0.0, 1.0]]) < 1e-15);
    assert!(max_abs_diff(&phi, &[[1.0, 0.5], [0.0, 1.0]]) < 1e-15);
    // Zero: both are the identity.
    let (e, phi) = mat_exp_phi1(&[[0.0; 3]; 3]);
    assert!(max_abs_diff(&e, &identity()) < 1e-15 && max_abs_diff(&phi, &identity()) < 1e-15);
    // Scalar blocks, including arguments that need many squarings:
    // e^a and (e^a - 1) / a.
    for a in [0.3_f64, -0.7, 3.0, -7.0, 25.0] {
        let (e, phi) = mat_exp_phi1(&[[a]]);
        let want_e = a.exp();
        let want_phi = (a.exp() - 1.0) / a;
        assert!(
            (e[0][0] / want_e - 1.0).abs() < 1e-12,
            "e^{a}: {} vs {want_e}",
            e[0][0]
        );
        assert!(
            (phi[0][0] / want_phi - 1.0).abs() < 1e-12,
            "phi1({a}): {} vs {want_phi}",
            phi[0][0]
        );
    }
    // A rotation generator: e^A is the rotation by the angle.
    let w = 1.3_f64;
    let (e, _) = mat_exp_phi1(&[[0.0, -w], [w, 0.0]]);
    assert!(max_abs_diff(&e, &[[w.cos(), -w.sin()], [w.sin(), w.cos()]]) < 1e-14);
    // For random matrices: A phi1(A) = e^A - I, det e^A = e^{tr A}, and
    // e^A e^{-A} = I.
    let mut rng = Rng::new(9);
    for scale in [0.1, 0.8, 2.5] {
        for _ in 0..5 {
            let a: [[f64; 3]; 3] = random_matrix(&mut rng, scale);
            let (e, phi) = mat_exp_phi1(&a);
            let mut shifted = e;
            for i in 0..3 {
                shifted[i][i] -= 1.0;
            }
            let scale_ref = e.iter().flatten().fold(1.0_f64, |m, x| m.max(x.abs()));
            assert!(
                max_abs_diff(&mat_mul(&a, &phi), &shifted) < 1e-11 * scale_ref,
                "A phi1(A) = e^A - I at scale {scale}"
            );
            let trace = a[0][0] + a[1][1] + a[2][2];
            assert!(
                (mat_det(&e) / trace.exp() - 1.0).abs() < 1e-10,
                "det e^A = e^tr A"
            );
            let negated = a.map(|r| r.map(|x| -x));
            let (e_inverse, _) = mat_exp_phi1(&negated);
            assert!(
                max_abs_diff(&mat_mul(&e, &e_inverse), &identity()) < 1e-9 * scale_ref * scale_ref
            );
        }
    }
    // The solution of x' = A x + b over a time t: x(t) = e^{At} x0 + t phi1(At) b,
    // checked for the scalar case x' = -2 x + 3 (x0 = 1): 1.5 - 0.5 e^{-2 t}.
    let t = 0.8;
    let (e, phi) = mat_exp_phi1(&[[-2.0 * t]]);
    let x = e[0][0] * 1.0 + t * phi[0][0] * 3.0;
    assert!(
        (x - (1.5 - 0.5 * (-2.0_f64 * t).exp())).abs() < 1e-14,
        "{x}"
    );
}

#[test]
fn barostat_constructors_degrees_of_freedom_and_energy() {
    // Degrees of freedom: D(D+1)/2 for the full cell, D for the box, 1 for the volume.
    for (d, full) in [(2_usize, 3_usize), (3, 6)] {
        assert_eq!(
            FlexibleBarostat::hydrostatic(d, 1.0, 10.0, CellMode::Full).cell_dof(),
            full
        );
        assert_eq!(
            FlexibleBarostat::hydrostatic(d, 1.0, 10.0, CellMode::Orthorhombic).cell_dof(),
            d
        );
        assert_eq!(
            FlexibleBarostat::hydrostatic(d, 1.0, 10.0, CellMode::Isotropic).cell_dof(),
            1
        );
    }
    // The hydrostatic target is p0 I.
    let b = FlexibleBarostat::hydrostatic(3, 2.5, 40.0, CellMode::Full);
    assert_eq!(b.d, 3);
    assert_eq!(b.w, 40.0);
    assert_eq!(b.p_ext, vec![2.5, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 2.5]);
    assert!(
        b.pg.iter().all(|&x| x == 0.0)
            && b.chain.is_none()
            && b.prescribed.is_none()
            && b.cell_friction == 0.0
    );
    // A general stress tensor is kept as given, the energy uses its mean diagonal.
    let tensor = vec![1.0, 0.1, 0.0, 0.1, 2.0, 0.0, 0.0, 0.0, 3.0];
    let t = FlexibleBarostat::with_tensor(3, tensor.clone(), 5.0, CellMode::Orthorhombic);
    assert_eq!(t.p_ext, tensor);
    assert!(
        (t.energy(10.0) - 2.0 * 10.0).abs() < 1e-14,
        "P V with the mean pressure 2: {}",
        t.energy(10.0)
    );
    // Energy = sum p_g^2 / (2 W) + p0 V.
    let mut moving = FlexibleBarostat::hydrostatic(2, 1.5, 4.0, CellMode::Full);
    moving.pg = vec![1.0, 2.0, -1.0, 0.5];
    let kinetic = (1.0 + 4.0 + 1.0 + 0.25) / 8.0;
    assert!(
        (moving.energy(6.0) - (kinetic + 9.0)).abs() < 1e-14,
        "{}",
        moving.energy(6.0)
    );
    // A Nose-Hoover chain on the cell adds its energy and sizes itself to the cell DOF.
    let chained =
        FlexibleBarostat::hydrostatic(3, 1.0, 10.0, CellMode::Full).with_chain(1.0, 1.0, 3);
    assert!(chained.chain.is_some());
    assert!(chained.energy(5.0) >= 5.0 - 1e-12);
    // The Langevin piston and the shear flow.
    let langevin =
        FlexibleBarostat::hydrostatic(3, 1.0, 10.0, CellMode::Isotropic).with_langevin_cell(0.7);
    assert_eq!(langevin.cell_friction, 0.7);
    let shear = FlexibleBarostat::shear_flow(3, 0.25);
    let g = shear.prescribed.as_ref().unwrap();
    assert_eq!(g.len(), 9);
    assert_eq!(g[1], 0.25, "u_x = rate * y: G_xy");
    assert!(g.iter().enumerate().all(|(i, &x)| i == 1 || x == 0.0));
    let flat = FlexibleBarostat::shear_flow(1, 0.25);
    assert!(
        flat.prescribed.unwrap().iter().all(|&x| x == 0.0),
        "no shear in one dimension"
    );
}
