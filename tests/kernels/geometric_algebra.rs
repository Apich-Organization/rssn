//! Cl(3,0) multivectors (ported from `numerical_geometric_algebra_test.rs`).

use rssn::kernels::geometric_algebra::Multivector3D;

#[test]

fn test_multivector_addition() {
    let a = Multivector3D::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);

    let b = Multivector3D::new(8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0);

    let c = a + b;

    assert_eq!(c.s, 9.0);

    assert_eq!(c.v1, 9.0);

    assert_eq!(c.v2, 9.0);

    assert_eq!(c.v3, 9.0);

    assert_eq!(c.b12, 9.0);

    assert_eq!(c.b23, 9.0);

    assert_eq!(c.b31, 9.0);

    assert_eq!(c.pss, 9.0);
}

#[test]

fn test_geometric_product_basis() {
    let e1 = Multivector3D::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    let e2 = Multivector3D::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    let e3 = Multivector3D::new(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);

    // e1 * e1 = 1
    assert_eq!((e1 * e1).s, 1.0);

    // e1 * e2 = e12
    let e12 = e1 * e2;

    assert_eq!(e12.b12, 1.0);

    assert_eq!(e12.s, 0.0);

    // e2 * e1 = -e12
    let e21 = e2 * e1;

    assert_eq!(e21.b12, -1.0);

    // e1 * e2 * e3 = pss
    let pss = e1 * e2 * e3;

    assert_eq!(pss.pss, 1.0);

    // pss * pss = (e1 e2 e3)(e1 e2 e3) = -1
    // (e1 e2 e3)(e1 e2 e3) = e1 e2 (e3 e1) e2 e3 = e1 e2 (-e1 e3) e2 e3 = -e1^2 e2^2 e3^2 = -1
    assert_eq!((pss * pss).s, -1.0);
}

#[test]

fn test_inner_outer_product() {
    let e1 = Multivector3D::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    let e2 = Multivector3D::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // e1 . e2 = 0
    assert_eq!(e1.dot(e2).s, 0.0);

    // e1 ^ e2 = e12
    assert_eq!(e1.wedge(e2).b12, 1.0);

    // v . v = |v|^2
    let v = Multivector3D::new(0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    assert_eq!(v.dot(v).s, 25.0);
}

#[test]

fn test_reverse_conjugate() {
    let a = Multivector3D::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0);

    let rev = a.reverse();

    assert_eq!(rev.s, 1.0);

    assert_eq!(rev.v1, 2.0);

    assert_eq!(rev.b12, -5.0);

    assert_eq!(rev.pss, -8.0);

    let conj = a.conjugate();

    assert_eq!(conj.s, 1.0);

    assert_eq!(conj.v1, -2.0);

    assert_eq!(conj.b12, -5.0);

    assert_eq!(conj.pss, 8.0);
}

#[test]

fn test_inverse() {
    let v = Multivector3D::new(0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0);

    let v_inv = v.inv().unwrap_or_else(|| panic!("vector is invertible"));

    let res = v * v_inv;

    assert!((res.s - 1.0).abs() < 1e-12);

    assert!(res.v1.abs() < 1e-12);
}

#[cfg(test)]
mod proptests {

    use assert_approx_eq::assert_approx_eq;
    use proptest::prelude::*;
    use proptest::test_runner::RngSeed;

    use super::*;

    prop_compose! {
        fn arb_multivector()(
            s in -100.0..100.0f64,
            v1 in -100.0..100.0f64,
            v2 in -100.0..100.0f64,
            v3 in -100.0..100.0f64,
            b12 in -100.0..100.0f64,
            b23 in -100.0..100.0f64,
            b31 in -100.0..100.0f64,
            pss in -100.0..100.0f64
        ) -> Multivector3D {
            Multivector3D::new(s, v1, v2, v3, b12, b23, b31, pss)
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig { rng_seed: RngSeed::Fixed(0x5EED), failure_persistence: None, ..ProptestConfig::default() })]

        #[test]
        fn test_addition_commutativity(a in arb_multivector(), b in arb_multivector()) {
            let res1 = a + b;
            let res2 = b + a;

            assert_approx_eq!(res1.s, res2.s);
            assert_approx_eq!(res1.v1, res2.v1);
            assert_approx_eq!(res1.v2, res2.v2);
            assert_approx_eq!(res1.v3, res2.v3);
            assert_approx_eq!(res1.b12, res2.b12);
            assert_approx_eq!(res1.b23, res2.b23);
            assert_approx_eq!(res1.b31, res2.b31);
            assert_approx_eq!(res1.pss, res2.pss);
        }

        #[test]
        fn test_addition_associativity(a in arb_multivector(), b in arb_multivector(), c in arb_multivector()) {
            let res1 = (a + b) + c;
            let res2 = a + (b + c);

            assert_approx_eq!(res1.s, res2.s);
            assert_approx_eq!(res1.v1, res2.v1);
            assert_approx_eq!(res1.v2, res2.v2);
            assert_approx_eq!(res1.v3, res2.v3);
            assert_approx_eq!(res1.b12, res2.b12);
            assert_approx_eq!(res1.b23, res2.b23);
            assert_approx_eq!(res1.b31, res2.b31);
            assert_approx_eq!(res1.pss, res2.pss);
        }

        #[test]
        fn test_geometric_product_associativity(a in arb_multivector(), b in arb_multivector(), c in arb_multivector()) {
            // (ab)c = a(bc)
            let ab_c = (a * b) * c;
            let a_bc = a * (b * c);

            // Due to floating point errors, we might need a looser tolerance or check differences
            let diff = ab_c - a_bc;
            assert!(diff.norm_sq() < 1e-6, "Associativity failed: norm_sq(diff) = {}", diff.norm_sq());
        }

        #[test]
        fn test_distributivity(a in arb_multivector(), b in arb_multivector(), c in arb_multivector()) {
            // a(b + c) = ab + ac
            let lhs = a * (b + c);
            let rhs = (a * b) + (a * c);

            let diff = lhs - rhs;
            assert!(diff.norm_sq() < 1e-6, "Distributivity failed: norm_sq(diff) = {}", diff.norm_sq());
        }

        #[test]
        fn test_reverse_property(a in arb_multivector(), b in arb_multivector()) {
            // reverse(ab) = reverse(b) * reverse(a)
            let lhs = (a * b).reverse();
            let rhs = b.reverse() * a.reverse();

            let diff = lhs - rhs;
            assert!(diff.norm_sq() < 1e-6, "Reverse property failed: norm_sq(diff) = {}", diff.norm_sq());
        }
    }
}

// ============================================================================
// Added: algebraic identities with independently known results
// ============================================================================

fn vec3(x: f64, y: f64, z: f64) -> Multivector3D {
    Multivector3D::new(0.0, x, y, z, 0.0, 0.0, 0.0, 0.0)
}

fn close(a: Multivector3D, b: Multivector3D, tol: f64) -> bool {
    (a - b).norm() < tol
}

#[test]
fn basis_bivectors_square_to_minus_one_and_pss_is_central() {
    let e1 = vec3(1.0, 0.0, 0.0);
    let e2 = vec3(0.0, 1.0, 0.0);
    let e3 = vec3(0.0, 0.0, 1.0);
    for (a, b) in [(e1, e2), (e2, e3), (e3, e1)] {
        let biv = a * b;
        assert_eq!((biv * biv).s, -1.0);
    }
    let i = e1 * e2 * e3;
    assert_eq!(i.pss, 1.0);
    for x in [e1, e2, e3, e1 * e2, Multivector3D::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)] {
        assert!(close(i * x, x * i, 1e-12), "pseudoscalar must commute with everything in G3");
    }
    // e2 e3 = I e1, e3 e1 = I e2 (duality)
    assert!(close(e2 * e3, i * e1, 1e-12));
    assert!(close(e3 * e1, i * e2, 1e-12));
}

#[test]
fn wedge_of_vectors_matches_cross_product_dual() {
    let a = vec3(1.0, 2.0, 3.0);
    let b = vec3(4.0, 5.0, 6.0);
    let w = a.wedge(b);
    // a x b = (-3, 6, -3) ; a ^ b = (a x b)_1 e23 + (a x b)_2 e31 + (a x b)_3 e12
    assert_eq!((w.b23, w.b31, w.b12), (-3.0, 6.0, -3.0));
    assert_eq!(w.s, 0.0);
    assert_eq!(a.wedge(a).norm(), 0.0);
    assert!(close(b.wedge(a), -w, 1e-12));
}

#[test]
fn geometric_product_of_vectors_is_dot_plus_wedge() {
    let a = vec3(1.0, -2.0, 0.5);
    let b = vec3(3.0, 1.0, -4.0);
    let prod = a * b;
    assert!(close(prod, a.dot(b) + a.wedge(b), 1e-12));
    assert!((prod.s - (3.0 - 2.0 - 2.0)).abs() < 1e-12);
}

#[test]
fn vector_bivector_product_splits_into_contraction_and_wedge() {
    let a = vec3(1.0, 2.0, 3.0);
    let b = Multivector3D::new(0.0, 0.0, 0.0, 0.0, 1.5, -0.5, 2.0, 0.0);
    assert!(close(a * b, a.dot(b) + a.wedge(b), 1e-12));
    assert!(close(b * a, b.dot(a) + b.wedge(a), 1e-12));
}

#[test]
fn rotor_rotates_a_vector_in_the_e12_plane() {
    let theta = 0.7f64;
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    let r = Multivector3D::new(c, 0.0, 0.0, 0.0, -s, 0.0, 0.0, 0.0);
    let e1 = vec3(1.0, 0.0, 0.0);
    let rotated = r * e1 * r.reverse();
    assert!(close(rotated, vec3(theta.cos(), theta.sin(), 0.0), 1e-12), "{rotated:?}");
    // The rotor has unit norm and R * reverse(R) = 1.
    assert!((r.norm() - 1.0).abs() < 1e-12);
    assert!(close(r * r.reverse(), Multivector3D::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 1e-12));
}

#[test]
fn norm_negation_subtraction_and_zero_inverse() {
    let a = Multivector3D::new(1.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(a.norm_sq(), 9.0);
    assert_eq!(a.norm(), 3.0);
    assert!(close(-a + a, Multivector3D::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 1e-15));
    assert!(close(a - a, Multivector3D::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 1e-15));
    assert!(Multivector3D::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0).inv().is_none());
}

#[test]
fn pure_pseudoscalar_and_vector_inverses() {
    let i = Multivector3D::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    let inv = i.inv().unwrap_or_else(|| panic!("I is invertible"));
    assert_eq!(inv.pss, -1.0); // I^-1 = -I since I^2 = -1
    let v = vec3(2.0, 0.0, 0.0);
    let vi = v.inv().unwrap_or_else(|| panic!("v is invertible"));
    assert_eq!(vi.v1, 0.5);
}

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig { rng_seed: proptest::test_runner::RngSeed::Fixed(0x5EED), failure_persistence: None, ..proptest::prelude::ProptestConfig::default() })]

    /// A vector times its inverse is the scalar 1.
    #[test]
    fn prop_vector_inverse(x in -10.0..10.0f64, y in -10.0..10.0f64, z in -10.0..10.0f64) {
        let v = vec3(x, y, z);
        proptest::prop_assume!(v.norm() > 1e-3);
        let inv = v.inv().ok_or_else(|| proptest::test_runner::TestCaseError::fail("no inverse"))?;
        proptest::prop_assert!(close(v * inv, Multivector3D::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 1e-9));
    }

    /// v * v = |v|^2 for vectors.
    #[test]
    fn prop_vector_squares_to_norm_squared(x in -10.0..10.0f64, y in -10.0..10.0f64, z in -10.0..10.0f64) {
        let v = vec3(x, y, z);
        let sq = v * v;
        proptest::prop_assert!((sq.s - (x * x + y * y + z * z)).abs() < 1e-9);
        proptest::prop_assert!(sq.v1.abs() + sq.v2.abs() + sq.v3.abs() + sq.b12.abs() + sq.b23.abs() + sq.b31.abs() < 1e-9);
    }

    /// Wedge and dot decompose the product of two vectors.
    #[test]
    fn prop_vector_product_decomposition(
        a in proptest::collection::vec(-10.0..10.0f64, 3), b in proptest::collection::vec(-10.0..10.0f64, 3),
    ) {
        let (u, w) = (vec3(a[0], a[1], a[2]), vec3(b[0], b[1], b[2]));
        proptest::prop_assert!(close(u * w, u.dot(w) + u.wedge(w), 1e-9));
    }
}
