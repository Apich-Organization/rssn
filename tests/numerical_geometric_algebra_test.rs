use rssn::prelude::numerical::numerical_Multivector3D as Multivector3D;

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

    let v_inv = v.inv().unwrap();

    let res = v * v_inv;

    assert!((res.s - 1.0).abs() < 1e-12);

    assert!(res.v1.abs() < 1e-12);
}

#[cfg(test)]

mod proptests {

    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use proptest::prelude::*;

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
