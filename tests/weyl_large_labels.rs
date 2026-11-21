//! Exact dimension and Casimir formulas at the full Dynkin label range.

use num_bigint::BigInt;
use rssn::kernels::lie_structure::RootSystem;
use rssn::kernels::lie_structure::cartan_matrix;
use rssn::kernels::qlinalg::Q;
use rssn::kernels::qlinalg::q;

fn root_system(
    letter: char,
    rank: usize,
) -> RootSystem {
    RootSystem::new(cartan_matrix(letter, rank).unwrap()).unwrap()
}

#[test]
fn a1_exact_formulas_cover_the_maximum_label_and_shift_boundaries() {
    let root = root_system('A', 1);
    for label in [0, 1, i64::MAX - 2, i64::MAX - 1, i64::MAX] {
        let l = BigInt::from(label);
        let dimension = &l + 1_i64;
        let casimir = Q::new(&l * (&l + 2_i64), 2.into());
        let index = Q::from_integer(dimension.clone()) * &casimir / q(3);
        assert_eq!(root.weyl_dim(&[label]), Some(dimension));
        assert_eq!(root.casimir(&[label]), Some(casimir));
        assert_eq!(root.dynkin_index(&[label]), Some(index));
    }
}

#[test]
fn signed_casimir_does_not_treat_negative_labels_as_highest_weights() {
    let root = root_system('A', 1);
    for label in [i64::MIN, -2, -1] {
        let l = BigInt::from(label);
        let expected = Q::new(&l * (&l + 2_i64), 2.into());
        assert_eq!(root.casimir(&[label]), Some(expected));
        assert!(root.weyl_dim(&[label]).is_none());
        assert!(root.dynkin_index(&[label]).is_none());
    }
    assert!(root.weyl_dim(&[]).is_none());
    assert!(root.casimir(&[0, 0]).is_none());
    assert!(root.dynkin_index(&[0, 0]).is_none());
}

#[test]
fn a2_dimension_and_weight_inner_product_match_closed_forms() {
    let root = root_system('A', 2);
    for labels in [[i64::MAX, i64::MAX], [i64::MAX, 0], [10, 19], [0, 0]] {
        let [a, b] = labels.map(BigInt::from);
        let dimension = (&a + 1_i64) * (&b + 1_i64) * (&a + &b + 2_i64) / 2_i64;
        let quadratic = &a * &a + &a * &b + &b * &b + 3_i64 * (&a + &b);
        let casimir = Q::new(2_i64 * quadratic, 3.into());
        assert_eq!(root.weyl_dim(&labels), Some(dimension));
        assert_eq!(root.casimir(&labels), Some(casimir));
    }
    let lhs = [i64::MIN, i64::MAX];
    let rhs = [i64::MAX - 1, i64::MIN + 2];
    let [a, b] = lhs.map(BigInt::from);
    let [c, d] = rhs.map(BigInt::from);
    let expected = Q::new(2_i64 * (&a * &c + &b * &d) + &a * &d + &b * &c, 3.into());
    assert_eq!(root.inner(&lhs, &rhs), expected);
    assert_eq!(root.inner(&rhs, &lhs), expected);
}

fn g2_dimension(
    a: &BigInt,
    b: &BigInt,
) -> BigInt {
    // The first fundamental representation has dimension seven and the
    // second fourteen with this Cartan matrix's numbering.
    (a + 1_i64)
        * (b + 1_i64)
        * (a + b + 2_i64)
        * (a + 2_i64 * b + 3_i64)
        * (a + 3_i64 * b + 4_i64)
        * (2_i64 * a + 3_i64 * b + 5_i64)
        / 120_i64
}

#[test]
fn g2_and_a_reducible_system_preserve_short_root_factors() {
    let g2 = root_system('G', 2);
    for labels in [[i64::MAX, i64::MAX - 1], [1, 0], [0, 1], [i64::MAX, 0]] {
        let [a, b] = labels.map(BigInt::from);
        let expected = g2_dimension(&a, &b);
        let casimir = Q::new(
            2_i64 * &a * &a + 6_i64 * &a * &b + 6_i64 * &b * &b + 10_i64 * &a + 18_i64 * &b,
            3.into(),
        );
        assert_eq!(g2.weyl_dim(&labels), Some(expected));
        assert_eq!(g2.casimir(&labels), Some(casimir));
    }

    let separate = RootSystem::new(vec![vec![2, 0, 0], vec![0, 2, -1], vec![0, -3, 2]]).unwrap();
    let labels = [i64::MAX, i64::MAX, 2];
    let [a, b, c] = labels.map(BigInt::from);
    assert_eq!(
        separate.weyl_dim(&labels),
        Some((a + 1_i64) * g2_dimension(&b, &c))
    );
    assert!(separate.casimir(&labels).is_none());
    assert!(separate.dynkin_index(&labels).is_none());
}
