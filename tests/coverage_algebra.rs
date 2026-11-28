//! Coverage tests: Lie algebra classification / structure theory and exact
//! homology (Smith normal form, standard spaces, persistence, cubical).
//! Public API only; known tables and independent identities are used.

#![allow(clippy::needless_range_loop)]

use std::collections::HashMap;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::homology::*;
use rssn::kernels::lie_structure::*;
use rssn::kernels::qlinalg::Q;
use rssn::kernels::qlinalg::QMat;
use rssn::kernels::qlinalg::q;
use rssn::kernels::qlinalg::{
    self,
};
use rssn::kernels::topology::Simplex;
use rssn::kernels::topology::euler_characteristic;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    fn int(
        &mut self,
        lo: i64,
        hi: i64,
    ) -> i64 {
        lo + (self.next() % (hi - lo + 1) as u64) as i64
    }
}

fn big(n: u64) -> BigInt {
    BigInt::from(n)
}

// ------------------------------------------------------- root systems

/// (letter, rank, #positive roots, |W|, Coxeter number, dual Coxeter number, det of Cartan matrix)
const TABLE: &[(char, usize, usize, u64, usize, i64, i64)] = &[
    ('A', 1, 1, 2, 2, 2, 2),
    ('A', 2, 3, 6, 3, 3, 3),
    ('A', 3, 6, 24, 4, 4, 4),
    ('A', 5, 15, 720, 6, 6, 6),
    ('B', 2, 4, 8, 4, 3, 2),
    ('B', 3, 9, 48, 6, 5, 2),
    ('B', 4, 16, 384, 8, 7, 2),
    ('C', 3, 9, 48, 6, 4, 2),
    ('C', 4, 16, 384, 8, 5, 2),
    ('D', 4, 12, 192, 6, 6, 4),
    ('D', 5, 20, 1920, 8, 8, 4),
    ('E', 6, 36, 51840, 12, 12, 3),
    ('E', 7, 63, 2_903_040, 18, 18, 2),
    ('E', 8, 120, 696_729_600, 30, 30, 1),
    ('F', 4, 24, 1152, 12, 9, 1),
    ('G', 2, 6, 12, 6, 4, 1),
];

fn root_system(
    letter: char,
    n: usize,
) -> RootSystem {
    RootSystem::new(cartan_matrix(letter, n).unwrap()).unwrap()
}

/// Dynkin labels of the highest root: `<theta, alpha_j^vee> = sum_i theta_i a_ij`.
fn theta_labels(rs: &RootSystem) -> Vec<i64> {
    let theta = rs.highest_root().unwrap();
    (0..rs.n)
        .map(|j| (0..rs.n).map(|i| theta[i] * rs.cartan[i][j]).sum())
        .collect()
}

#[test]
fn root_system_invariants_match_the_classification_tables() {
    for &(letter, n, npos, order, h, hv, det) in TABLE {
        let name = format!("{letter}{n}");
        let a = cartan_matrix(letter, n).unwrap();
        assert_eq!(parse_type(&name).as_ref(), Some(&a), "{name}");
        // Cartan matrices: diagonal 2, nonpositive off-diagonal, determinant table
        for i in 0..n {
            assert_eq!(a[i][i], 2);
            for j in 0..n {
                if i != j {
                    assert!(a[i][j] <= 0 && (a[i][j] == 0) == (a[j][i] == 0), "{name}");
                }
            }
        }
        let aq: QMat = a
            .iter()
            .map(|r| r.iter().map(|&v| q(v)).collect())
            .collect();
        assert_eq!(qlinalg::det(&aq), q(det), "det {name}");
        assert_eq!(recognize(&a).map(|v| v.len()), Some(1), "{name}");
        assert_eq!(positive_roots(&a).len(), npos, "{name}");
        let rs = RootSystem::new(a.clone()).unwrap();
        assert!(rs.is_irreducible());
        assert_eq!(rs.positive.len(), npos);
        assert_eq!(rs.weyl_order(), big(order), "|W| {name}");
        assert_eq!(rs.coxeter_number(), Some(h), "h {name}");
        assert_eq!(rs.dual_coxeter_number(), Some(q(hv)), "h^v {name}");
        // exponents sum to the number of positive roots; largest is h - 1
        let ex = rs.exponents();
        assert_eq!(ex.len(), n);
        assert_eq!(ex.iter().sum::<usize>(), npos, "{name}");
        assert_eq!(*ex.last().unwrap(), h - 1, "{name}");
        // positive roots are distinct, ordered by height
        let heights: Vec<i64> = rs.positive.iter().map(|r| r.iter().sum()).collect();
        assert!(heights.windows(2).all(|w| w[0] <= w[1]), "{name}");
        let mut uniq = rs.positive.clone();
        uniq.sort();
        uniq.dedup();
        assert_eq!(uniq.len(), npos);
        assert!(rs.positive.iter().all(|r| r.iter().all(|&c| c >= 0)));
        // Gram matrix symmetric
        let g = rs.gram();
        for i in 0..n {
            for j in 0..n {
                assert_eq!(g[i][j], g[j][i], "gram {name}");
            }
        }
        // adjoint representation: dimension n + 2N, Casimir and Dynkin index 2 h^vee
        let labels = theta_labels(&rs);
        assert!(labels.iter().all(|&l| l >= 0), "{name}");
        assert_eq!(
            rs.weyl_dim(&labels),
            Some(big((n + 2 * npos) as u64)),
            "dim adjoint {name}"
        );
        assert_eq!(rs.casimir(&labels), Some(q(2 * hv)), "casimir {name}");
        assert_eq!(rs.dynkin_index(&labels), Some(q(2 * hv)), "index {name}");
        // trivial representation
        let zero = vec![0; n];
        assert_eq!(rs.weyl_dim(&zero), Some(BigInt::one()));
        assert_eq!(rs.casimir(&zero), Some(q(0)));
    }
}

#[test]
fn well_known_representation_dimensions() {
    let dim = |l: char, n: usize, labels: &[i64]| root_system(l, n).weyl_dim(labels).unwrap();
    let e = |n: usize, k: usize| {
        let mut v = vec![0; n];
        v[k] = 1;
        v
    };
    assert_eq!(dim('E', 6, &e(6, 0)), big(27));
    assert_eq!(dim('E', 7, &e(7, 6)), big(56));
    assert_eq!(dim('E', 8, &e(8, 7)), big(248));
    assert_eq!(dim('F', 4, &e(4, 3)), big(26));
    assert_eq!(dim('F', 4, &e(4, 0)), big(52));
    let mut g2 = [dim('G', 2, &e(2, 0)), dim('G', 2, &e(2, 1))];
    g2.sort();
    assert_eq!(g2, [big(7), big(14)]);
    // B3: vector 7 and spinor 8; C3: 6, 14, 14; D4: all three 8-dimensional
    assert_eq!(dim('B', 3, &e(3, 0)), big(7));
    assert_eq!(dim('B', 3, &e(3, 2)), big(8));
    assert_eq!(dim('C', 3, &e(3, 0)), big(6));
    assert_eq!(dim('C', 3, &e(3, 1)), big(14));
    assert_eq!(dim('D', 4, &e(4, 0)), big(8));
    assert_eq!(dim('D', 4, &e(4, 2)), big(8));
    assert_eq!(dim('D', 4, &e(4, 3)), big(8));
    // A_n: exterior powers C(n+1, k); sl2: k + 1; sl3: (a+1)(b+1)(a+b+2)/2
    for (k, want) in [5u64, 10, 10, 5].into_iter().enumerate() {
        assert_eq!(dim('A', 4, &e(4, k)), big(want));
    }
    for k in 0..30 {
        assert_eq!(dim('A', 1, &[k]), big(k as u64 + 1));
    }
    for a in 0..6i64 {
        for b in 0..6i64 {
            assert_eq!(
                dim('A', 2, &[a, b]),
                big(((a + 1) * (b + 1) * (a + b + 2) / 2) as u64)
            );
        }
    }
    // invalid labels
    let a2 = root_system('A', 2);
    assert!(a2.weyl_dim(&[1]).is_none());
    assert!(a2.weyl_dim(&[1, -1]).is_none());
    assert!(a2.casimir(&[1]).is_none());
    // fundamental of su(n): Casimir (n^2-1)/n and Dynkin index 1
    for n in 2..6usize {
        let rs = root_system('A', n);
        // A_n is su(n+1)
        let m = n + 1;
        assert_eq!(
            rs.casimir(&e(n, 0)),
            Some(Q::new(((m * m - 1) as i64).into(), (m as i64).into()))
        );
        assert_eq!(rs.dynkin_index(&e(n, 0)), Some(q(1)), "A{n}");
    }
}

#[test]
fn cartan_matrix_parsing_and_recognition() {
    for bad in ["", "A", "A0", "X3", "A65", "3A", "E9x", "A-1"] {
        assert!(parse_type(bad).is_none(), "{bad:?}");
    }
    assert!(cartan_matrix('A', 0).is_none());
    assert!(cartan_matrix('A', 65).is_none());
    assert!(cartan_matrix('Z', 3).is_none());
    assert_eq!(parse_type(" e8 ").unwrap(), cartan_matrix('E', 8).unwrap());
    // reducible: A2 + G2 in block form
    let (a, g) = (
        cartan_matrix('A', 2).unwrap(),
        cartan_matrix('G', 2).unwrap(),
    );
    let mut m = vec![vec![0; 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            m[i][j] = a[i][j];
            m[2 + i][2 + j] = g[i][j];
        }
    }
    let kinds = recognize(&m).unwrap();
    assert_eq!(kinds.len(), 2);
    assert!(kinds.contains(&('A', 2)) && kinds.contains(&('G', 2)));
    assert_eq!(components(&m).len(), 2);
    let rs = RootSystem::new(m).unwrap();
    assert!(!rs.is_irreducible());
    assert!(rs.highest_root().is_none());
    assert!(rs.coxeter_number().is_none());
    assert!(rs.casimir(&[1, 0, 0, 0]).is_none());
    assert!(rs.dominant_multiplicities(&[1, 0, 0, 0], 100).is_none());
    assert_eq!(rs.weyl_order(), big(6 * 12));
    // not Cartan matrices
    assert!(recognize(&vec![vec![2, 1], vec![1, 2]]).is_none());
    assert!(recognize(&vec![vec![2, -2], vec![-2, 2]]).is_none());
    assert!(RootSystem::new(vec![vec![2, -3], vec![-3, 2]]).is_none());
    // type-recognition is invariant under relabelling of the nodes (D4, swap 0 and 3)
    let d4 = cartan_matrix('D', 4).unwrap();
    let perm = [3, 1, 2, 0];
    let relabelled: Vec<Vec<i64>> = (0..4)
        .map(|i| (0..4).map(|j| d4[perm[i]][perm[j]]).collect())
        .collect();
    assert_eq!(recognize(&relabelled), Some(vec![('D', 4)]));
    // Euclidean simple roots reproduce the Cartan matrix: a_ij = 2 (a_i, a_j) / (a_j, a_j)
    for (l, n) in [
        ('A', 3),
        ('B', 3),
        ('C', 3),
        ('D', 4),
        ('F', 4),
        ('G', 2),
        ('E', 6),
    ] {
        if let Some(roots) = euclidean_simple_roots(l, n) {
            let a = cartan_matrix(l, n).unwrap();
            for i in 0..n {
                for j in 0..n {
                    let dot = roots[i]
                        .iter()
                        .zip(&roots[j])
                        .fold(q(0), |acc, (x, y)| acc + x.clone() * y.clone());
                    let nj = roots[j]
                        .iter()
                        .fold(q(0), |acc, y| acc + y.clone() * y.clone());
                    assert_eq!(q(2) * dot / nj, q(a[i][j]), "{l}{n} ({i},{j})");
                }
            }
        }
    }
}

#[test]
fn weights_multiplicities_and_tensor_products() {
    // adjoint of su(3): zero weight has multiplicity 2, six roots with multiplicity 1
    let a2 = root_system('A', 2);
    let w = a2.all_weights(&[1, 1], 1000).unwrap();
    assert_eq!(w.len(), 7);
    assert_eq!(w[&vec![0, 0]], big(2));
    assert_eq!(w.values().sum::<BigInt>(), big(8));
    assert_eq!(a2.weight_multiplicity(&[1, 1], &[0, 0]), Some(big(2)));
    assert_eq!(a2.weight_multiplicity(&[1, 1], &[1, 1]), Some(big(1)));
    assert_eq!(
        a2.weight_multiplicity(&[1, 1], &[5, 5]),
        Some(BigInt::zero())
    );
    assert!(a2.weight_multiplicity(&[1, 1], &[0]).is_none());
    // cap respected, wrong sizes rejected
    assert!(a2.all_weights(&[3, 3], 3).is_none());
    assert!(a2.dominant_multiplicities(&[1], 100).is_none());
    assert!(a2.dominant_multiplicities(&[-1, 0], 100).is_none());
    // the sum of all weight multiplicities is the dimension (several types)
    for (l, n, top) in [
        ('A', 3, vec![1, 0, 1]),
        ('B', 2, vec![1, 1]),
        ('G', 2, vec![0, 1]),
        ('G', 2, vec![2, 0]),
        ('C', 3, vec![0, 1, 0]),
        ('D', 4, vec![1, 0, 0, 0]),
    ] {
        let rs = root_system(l, n);
        let weights = rs.all_weights(&top, 100_000).unwrap();
        assert_eq!(
            weights.values().sum::<BigInt>(),
            rs.weyl_dim(&top).unwrap(),
            "{l}{n} {top:?}"
        );
        // multiplicities are invariant under the Weyl group: reflect each weight in simple root i
        for (wt, m) in &weights {
            for i in 0..n {
                let refl: Vec<i64> = wt
                    .iter()
                    .zip(&rs.cartan[i])
                    .map(|(x, a)| x - wt[i] * a)
                    .collect();
                assert_eq!(weights.get(&refl), Some(m), "{l}{n}");
            }
        }
    }
    // su(2): 1 x 1 = 2 + 0;   su(3): 3 x 3 = 6 + 3bar
    let a1 = root_system('A', 1);
    assert_eq!(
        a1.tensor(&[1], &[1]).unwrap(),
        vec![(vec![0], big(1)), (vec![2], big(1))]
    );
    let t = a2.tensor(&[1, 0], &[1, 0]).unwrap();
    assert_eq!(t, vec![(vec![0, 1], big(1)), (vec![2, 0], big(1))]);
    // dimensions multiply, in several algebras and symmetric in the factors
    for (l, n) in [('A', 2), ('A', 3), ('B', 2), ('G', 2), ('C', 3)] {
        let rs = root_system(l, n);
        let mut rng = Lcg((n * 17) as u64);
        for _ in 0..4 {
            let a: Vec<i64> = (0..n).map(|_| rng.int(0, 1)).collect();
            let b: Vec<i64> = (0..n).map(|_| rng.int(0, 1)).collect();
            let prod = rs.weyl_dim(&a).unwrap() * rs.weyl_dim(&b).unwrap();
            let t = rs.tensor(&a, &b).unwrap();
            let total: BigInt = t.iter().map(|(hw, m)| rs.weyl_dim(hw).unwrap() * m).sum();
            assert_eq!(total, prod, "{l}{n}: {a:?} x {b:?}");
            assert!(t.iter().all(|(_, m)| m > &BigInt::zero()));
            assert_eq!(rs.tensor(&b, &a).unwrap(), t);
        }
    }
}

// ------------------------------------------------ Lie algebra structure

fn real_basis(
    family: Family,
    n: usize,
) -> Vec<QMat> {
    standard_basis(family, n)
        .unwrap()
        .iter()
        .map(|m| {
            m.iter()
                .map(|r| {
                    r.iter()
                        .map(|&(a, b)| {
                            assert_eq!(b, 0);
                            q(a)
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

fn trace_of_product(
    a: &QMat,
    b: &QMat,
) -> Q {
    let p = qlinalg::matmul(a, b);
    (0..p.len()).fold(q(0), |s, i| s + p[i][i].clone())
}

#[test]
fn classical_algebras_have_the_expected_killing_forms() {
    // (family, n, dimension, Killing = c * trace form)
    let cases: Vec<(Family, usize, usize, i64)> = vec![
        (Family::Sl, 2, 3, 4),
        (Family::Sl, 3, 8, 6),
        (Family::Sl, 4, 15, 8),
        (Family::So, 3, 3, 1),
        (Family::So, 4, 6, 2),
        (Family::So, 5, 10, 3),
        (Family::Sp, 1, 3, 4),
        (Family::Sp, 2, 10, 6),
    ];
    for (family, n, dim, c) in cases {
        let basis = real_basis(family, n);
        assert_eq!(basis.len(), dim, "{family:?}({n})");
        let lie =
            Lie::from_matrices(&basis).unwrap_or_else(|| panic!("{family:?}({n}) not closed"));
        assert_eq!(lie.n, dim);
        let k = lie.killing();
        for a in 0..dim {
            for b in 0..dim {
                assert_eq!(
                    k[a][b],
                    q(c) * trace_of_product(&basis[a], &basis[b]),
                    "{family:?}({n}) K[{a}][{b}]"
                );
            }
        }
        assert!(lie.is_semisimple(), "{family:?}({n})");
        assert!(!lie.is_solvable());
        assert!(lie.center().is_empty());
        assert!(lie.radical().is_empty());
        assert!(!lie.cartan_solvable_test());
        assert!(
            lie.derived_series().last().is_some_and(|d| d.len() == dim),
            "perfect"
        );
        let (s, r) = lie.levi().unwrap();
        assert_eq!(s.len(), dim);
        assert!(r.is_empty());
    }
    // so(2) is abelian, not semisimple
    let so2 = Lie::from_matrices(&real_basis(Family::So, 2)).unwrap();
    assert!(!so2.is_semisimple());
    assert!(so2.is_solvable() && so2.is_nilpotent());
    // out-of-range requests
    assert!(standard_basis(Family::Sl, 1).is_none());
    assert!(standard_basis(Family::Gl, 0).is_none());
    assert!(standard_basis(Family::So, 33).is_none());
}

#[test]
fn gl_is_reductive_with_a_one_dimensional_center() {
    for n in 2..4usize {
        let lie = Lie::from_matrices(&real_basis(Family::Gl, n)).unwrap();
        assert_eq!(lie.n, n * n);
        assert!(!lie.is_semisimple());
        assert!(!lie.is_solvable());
        assert_eq!(lie.center().len(), 1);
        assert_eq!(lie.radical().len(), 1);
        let (s, r) = lie.levi().unwrap();
        assert_eq!((s.len(), r.len()), (n * n - 1, 1));
        assert!(lie.is_subalgebra(&s));
        assert!(lie.is_subalgebra(&r));
    }
    // a non-closed set of matrices is not a Lie algebra; dependent ones are rejected
    let e12 = vec![vec![q(0), q(1)], vec![q(0), q(0)]];
    let e21 = vec![vec![q(0), q(0)], vec![q(1), q(0)]];
    assert!(Lie::from_matrices(&[e12.clone(), e21]).is_none());
    assert!(Lie::from_matrices(&[e12.clone(), e12]).is_none());
    assert!(Lie::from_matrices(&[]).is_none());
}

fn constants(
    n: usize,
    entries: &[(usize, usize, usize, i64)],
) -> Vec<Vec<Vec<Q>>> {
    let mut c = vec![vec![vec![q(0); n]; n]; n];
    for &(i, j, k, v) in entries {
        c[i][j][k] = q(v);
        c[j][i][k] = q(-v);
    }
    c
}

#[test]
fn solvable_nilpotent_and_abelian_algebras() {
    // Heisenberg [x, y] = z
    let h = Lie::from_constants(constants(3, &[(0, 1, 2, 1)])).unwrap();
    assert!(h.is_nilpotent() && h.is_solvable());
    assert!(!h.is_semisimple());
    assert!(h.cartan_solvable_test());
    assert_eq!(h.center().len(), 1);
    assert_eq!(h.upper_central_series().last().unwrap().len(), 3);
    assert!(h.lower_central_series().last().unwrap().is_empty());
    assert_eq!(h.radical().len(), 3);
    assert_eq!(h.derived_series()[1].len(), 1);
    // ax + b: [x, y] = y is solvable but not nilpotent
    let aff = Lie::from_constants(constants(2, &[(0, 1, 1, 1)])).unwrap();
    assert!(aff.is_solvable());
    assert!(!aff.is_nilpotent());
    assert!(aff.center().is_empty());
    assert!(aff.cartan_solvable_test());
    // abelian
    let ab = Lie::from_constants(constants(3, &[])).unwrap();
    assert!(ab.is_nilpotent());
    assert_eq!(ab.center().len(), 3);
    assert_eq!(
        ab.bracket(&[q(1), q(2), q(3)], &[q(4), q(5), q(6)]),
        vec![q(0); 3]
    );
    // sl2 from structure constants: [h,e]=2e, [h,f]=-2f, [e,f]=h with K(h,h)=8, K(e,f)=4
    let sl2 =
        Lie::from_constants(constants(3, &[(0, 1, 1, 2), (0, 2, 2, -2), (1, 2, 0, 1)])).unwrap();
    let k = sl2.killing();
    assert_eq!(
        (k[0][0].clone(), k[1][2].clone(), k[0][1].clone()),
        (q(8), q(4), q(0))
    );
    assert!(sl2.is_semisimple());
    // ad is a homomorphism: ad[x,y] = [ad x, ad y]
    let mut rng = Lcg(5);
    for _ in 0..10 {
        let x: Vec<Q> = (0..3).map(|_| q(rng.int(-3, 3))).collect();
        let y: Vec<Q> = (0..3).map(|_| q(rng.int(-3, 3))).collect();
        let (ax, ay) = (sl2.ad(&x), sl2.ad(&y));
        let comm: QMat = {
            let (p, r) = (qlinalg::matmul(&ax, &ay), qlinalg::matmul(&ay, &ax));
            (0..3)
                .map(|i| (0..3).map(|j| p[i][j].clone() - r[i][j].clone()).collect())
                .collect()
        };
        assert_eq!(sl2.ad(&sl2.bracket(&x, &y)), comm);
        // antisymmetry
        let (xy, yx) = (sl2.bracket(&x, &y), sl2.bracket(&y, &x));
        assert!(
            xy.iter()
                .zip(&yx)
                .all(|(a, b)| a.clone() + b.clone() == q(0))
        );
    }
    // invalid constants
    assert!(Lie::from_constants(vec![]).is_none());
    let mut skew_broken = constants(2, &[(0, 1, 1, 1)]);
    skew_broken[1][0][1] = q(5);
    assert!(Lie::from_constants(skew_broken).is_none());
    assert!(
        Lie::from_constants(constants(3, &[(0, 1, 1, 1), (1, 2, 2, 1), (2, 0, 0, 1)])).is_none(),
        "Jacobi fails"
    );
    let mut ragged = constants(2, &[]);
    ragged[0].pop();
    assert!(Lie::from_constants(ragged).is_none());
}

#[test]
fn cartan_subalgebra_and_root_decomposition_of_sl3() {
    let basis = real_basis(Family::Sl, 3);
    let lie = Lie::from_matrices(&basis).unwrap();
    let h = lie.cartan_subalgebra();
    assert_eq!(h.len(), 2, "rank of sl3");
    assert!(lie.is_subalgebra(&h));
    for x in &h {
        for y in &h {
            assert!(qlinalg::is_zero_vec(&lie.bracket(x, y)));
        }
    }
    let (spaces, leftover) = lie.root_decomposition(&h).unwrap();
    assert_eq!(leftover, 0);
    assert_eq!(spaces.iter().map(|s| s.space.len()).sum::<usize>(), 8);
    let nonzero: Vec<_> = spaces
        .iter()
        .filter(|s| s.root.iter().any(|v| !v.is_zero()))
        .collect();
    assert_eq!(nonzero.len(), 6);
    assert!(nonzero.iter().all(|s| s.space.len() == 1));
    // roots come in +/- pairs
    for s in &nonzero {
        let neg: Vec<Q> = s.root.iter().map(|v| -v.clone()).collect();
        assert!(nonzero.iter().any(|t| t.root == neg));
    }
    // non-abelian input is refused
    let all = qlinalg::identity(8);
    assert!(lie.root_decomposition(&all).is_none());
    // sl2: roots +-2 on the Chevalley element
    let sl2 = Lie::from_matrices(&real_basis(Family::Sl, 2)).unwrap();
    let h2 = sl2.cartan_subalgebra();
    assert_eq!(h2.len(), 1);
    let (sp, left) = sl2.root_decomposition(&h2).unwrap();
    assert_eq!((sp.len(), left), (3, 0));
}

#[test]
fn casimir_of_the_defining_representation_is_a_scalar() {
    for n in [2usize, 3] {
        let basis = real_basis(Family::Sl, n);
        let lie = Lie::from_matrices(&basis).unwrap();
        let c = lie.casimir(&basis).unwrap();
        let want = Q::new(((n * n - 1) as i64).into(), ((2 * n * n) as i64).into());
        for i in 0..n {
            for j in 0..n {
                assert_eq!(c[i][j], if i == j { want.clone() } else { q(0) }, "sl{n}");
            }
        }
    }
    // wrong number of matrices, non-semisimple algebra
    let lie = Lie::from_matrices(&real_basis(Family::Sl, 2)).unwrap();
    assert!(lie.casimir(&real_basis(Family::Sl, 2)[..2]).is_none());
    let ab = Lie::from_constants(constants(2, &[])).unwrap();
    assert!(
        ab.casimir(&[qlinalg::identity(2), qlinalg::identity(2)])
            .is_none()
    );
}

#[test]
fn charpoly_and_rational_roots() {
    let m: QMat = vec![vec![q(2), q(1)], vec![q(1), q(2)]];
    assert_eq!(charpoly(&m), vec![q(3), q(-4), q(1)]);
    assert_eq!(charpoly(&[]), vec![q(1)]);
    // Cayley-Hamilton on random matrices
    let mut rng = Lcg(3);
    for n in 1..6usize {
        let a: QMat = (0..n)
            .map(|_| (0..n).map(|_| q(rng.int(-3, 3))).collect())
            .collect();
        let p = charpoly(&a);
        assert_eq!(p.len(), n + 1);
        assert_eq!(p[n], q(1));
        let mut acc: QMat = vec![vec![q(0); n]; n];
        let mut power = qlinalg::identity(n);
        for coeff in &p {
            for i in 0..n {
                for j in 0..n {
                    acc[i][j] = acc[i][j].clone() + coeff.clone() * power[i][j].clone();
                }
            }
            power = qlinalg::matmul(&a, &power);
        }
        assert!(
            acc.iter().all(|r| qlinalg::is_zero_vec(r)),
            "p(A) = 0, n = {n}"
        );
        // constant term is (-1)^n det
        let sign = if n % 2 == 0 { 1 } else { -1 };
        assert_eq!(p[0], q(sign) * qlinalg::det(&a));
    }
    let mut roots = rational_roots(&[q(3), q(-4), q(1)]).unwrap();
    roots.sort();
    assert_eq!(roots, vec![q(1), q(3)]);
    let mut roots = rational_roots(&[q(-2), q(3), q(2)]).unwrap();
    roots.sort();
    assert_eq!(roots, vec![q(-2), Q::new(1.into(), 2.into())]);
    assert!(rational_roots(&[q(1), q(0), q(1)]).unwrap().is_empty());
    assert!(
        rational_roots(&[q(0), q(-1), q(1)])
            .unwrap()
            .contains(&q(0))
    );
    assert!(rational_roots(&[q(0), q(0)]).is_none());
}

// ------------------------------------------------------------- homology

fn bi(v: &[u64]) -> Vec<BigInt> {
    v.iter().map(|&x| big(x)).collect()
}

#[test]
fn smith_normal_form_properties() {
    assert_eq!(smith_invariants(&[vec![2, 0], vec![0, 3]]), bi(&[1, 6]));
    assert_eq!(smith_invariants(&[vec![2, 4], vec![6, 8]]), bi(&[2, 4]));
    assert_eq!(
        smith_invariants(&[vec![0, 0], vec![0, 0]]),
        Vec::<BigInt>::new()
    );
    assert_eq!(smith_invariants(&[]), Vec::<BigInt>::new());
    assert_eq!(smith_invariants(&[vec![6]]), bi(&[6]));
    assert_eq!(smith_invariants(&[vec![-6]]), bi(&[6]));
    let mut rng = Lcg(12);
    for _ in 0..40 {
        let (m, n) = (rng.int(1, 5) as usize, rng.int(1, 5) as usize);
        let a: Vec<Vec<i64>> = (0..m)
            .map(|_| (0..n).map(|_| rng.int(-4, 4)).collect())
            .collect();
        let inv = smith_invariants(&a);
        // rank agrees with rational rank and with rank mod a large prime
        let aq: QMat = a
            .iter()
            .map(|r| r.iter().map(|&v| q(v)).collect())
            .collect();
        assert_eq!(inv.len(), qlinalg::rank(&aq, n));
        assert_eq!(rank_mod_p(&a, 1_000_003), inv.len());
        // divisibility chain
        for w in inv.windows(2) {
            assert!((&w[1] % &w[0]).is_zero(), "{inv:?}");
        }
        // product of invariants = |det| for square nonsingular matrices
        if m == n && inv.len() == n {
            let d = qlinalg::det(&aq).to_integer().magnitude().clone();
            assert_eq!(inv.iter().product::<BigInt>().magnitude().clone(), d);
        }
        // invariants are invariant under transposition
        let t: Vec<Vec<i64>> = (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect();
        assert_eq!(smith_invariants(&t), inv);
    }
    assert_eq!(rank_mod_p(&[vec![2, 0], vec![0, 3]], 2), 1);
    assert_eq!(rank_mod_p(&[vec![2, 0], vec![0, 3]], 3), 1);
    assert_eq!(rank_mod_p(&[vec![2, 0], vec![0, 3]], 5), 2);
    assert_eq!(rank_mod_p(&[], 5), 0);
    let sieve = |n: u64| (2..n).all(|d| !n.is_multiple_of(d)) && n >= 2;
    for n in 0..300 {
        assert_eq!(is_prime(n), sieve(n), "{n}");
    }
}

fn circle() -> Vec<Simplex> {
    sphere(1)
}

#[test]
fn homology_of_standard_spaces() {
    for n in 1..5usize {
        let s = sphere(n);
        let mut want = vec![0; n + 1];
        want[0] = 1;
        want[n] = 1;
        assert_eq!(betti_list(&s), want, "S^{n}");
        let mut contractible = vec![0; n + 1];
        contractible[0] = 1;
        assert_eq!(betti_list(&full_simplex(n)), contractible, "simplex {n}");
        // Euler characteristic from simplex counts equals the alternating Betti sum
        let chi: i64 = betti_list(&s)
            .iter()
            .enumerate()
            .map(|(k, &b)| {
                if k % 2 == 0 {
                    b as i64
                } else {
                    -(b as i64)
                }
            })
            .sum();
        assert_eq!(euler_characteristic(&s), chi);
    }
    assert_eq!(betti_list(&sphere(0)), vec![2]);
    // cone is contractible whatever it is built on; suspension raises the sphere dimension
    for base in [circle(), sphere(2), projective_plane(0)] {
        let c = complex_chains(&cone(&base));
        for k in 0..4 {
            let g = reduced_homology(&c, k);
            assert!(g.rank == 0 && g.torsion.is_empty(), "cone, degree {k}");
        }
    }
    assert_eq!(betti_list(&suspension(&sphere(2))), betti_list(&sphere(3)));
    // join of spheres: S^a * S^b = S^(a+b+1)
    assert_eq!(
        betti_list(&join(&sphere(0), &sphere(0))),
        betti_list(&sphere(1))
    );
    assert_eq!(
        betti_list(&join(&circle(), &circle())),
        betti_list(&sphere(3))
    );
    // wedge of two circles
    assert_eq!(betti_list(&wedge(&circle(), &circle(), 0, 0)), vec![1, 2]);
    // torus as a product of circles
    let torus = product(&circle(), &circle(), 100_000).unwrap();
    assert_eq!(betti_list(&torus), vec![1, 2, 1]);
    assert_eq!(euler_characteristic(&torus), 0);
    assert!(product(&circle(), &circle(), 5).is_none());
    // subcomplex helpers: |A u B| + |A n B| = |A| + |B|
    let a = full_simplex(2);
    let b = sphere(2)
        .into_iter()
        .filter(|s| s.iter().all(|&v| v != 0) || s.len() == 1)
        .collect::<Vec<_>>();
    let b = union(&b, &[vec![3]]);
    let (u, i) = (union(&a, &b), intersection(&a, &b));
    assert_eq!(u.len() + i.len(), a.len() + b.len());
    assert!(
        is_subcomplex(&a, &u)
            && is_subcomplex(&b, &u)
            && is_subcomplex(&i, &a)
            && is_subcomplex(&i, &b)
    );
    assert!(!is_subcomplex(&u, &a) || a.len() == u.len());
}

#[test]
fn torsion_in_surfaces() {
    // RP^2: H_1 = Z/2, H_2 = 0 over Z; over Z/2 the Betti numbers are 1,1,1; over Z/3 they are 1,0,0
    for plane in [projective_plane(0), projective_plane(3)] {
        let cc = complex_chains(&plane);
        assert_eq!(
            cc.homology(0),
            Group {
                rank: 1,
                torsion: vec![]
            }
        );
        assert_eq!(
            cc.homology(1),
            Group {
                rank: 0,
                torsion: bi(&[2])
            }
        );
        assert_eq!(
            cc.homology(2),
            Group {
                rank: 0,
                torsion: vec![]
            }
        );
        assert_eq!(
            [
                cc.homology_mod_p(0, 2),
                cc.homology_mod_p(1, 2),
                cc.homology_mod_p(2, 2)
            ],
            [1, 1, 1]
        );
        assert_eq!(
            [
                cc.homology_mod_p(0, 3),
                cc.homology_mod_p(1, 3),
                cc.homology_mod_p(2, 3)
            ],
            [1, 0, 0]
        );
        // universal coefficients: H^2 carries the torsion of H_1
        assert_eq!(
            cc.cohomology(1),
            Group {
                rank: 0,
                torsion: vec![]
            }
        );
        assert_eq!(
            cc.cohomology(2),
            Group {
                rank: 0,
                torsion: bi(&[2])
            }
        );
        assert_eq!(cc.top_degree(), Some(2));
    }
    // Klein bottle: H_1 = Z + Z/2
    let kb = complex_chains(&klein_bottle(4, 4));
    assert_eq!(kb.homology(0).rank, 1);
    assert_eq!(
        kb.homology(1),
        Group {
            rank: 1,
            torsion: bi(&[2])
        }
    );
    assert_eq!(
        kb.homology(2),
        Group {
            rank: 0,
            torsion: vec![]
        }
    );
    // dunce cap is acyclic
    assert_eq!(betti_list(&dunce_cap(3)), vec![1, 0, 0]);
    let dc = complex_chains(&dunce_cap(3));
    assert!(dc.homology(1).torsion.is_empty() && dc.homology(2).torsion.is_empty());
    // relative homology of (disc, boundary) is that of S^2 shifted: H_2 = Z
    let disc = full_simplex(2);
    let rel = relative_chains(&disc, &circle());
    assert_eq!(
        [
            rel.homology(0).rank,
            rel.homology(1).rank,
            rel.homology(2).rank
        ],
        [0, 0, 1]
    );
    // empty complex
    let empty = complex_chains(&[]);
    assert_eq!(empty.top_degree(), None);
    assert_eq!(
        empty.homology(0),
        Group {
            rank: 0,
            torsion: vec![]
        }
    );
}

#[test]
fn fundamental_groups() {
    let abel =
        |c: &[Simplex]| abelianization(&simplify_presentation(&pi1_presentation(c).unwrap()));
    assert_eq!(
        abel(&sphere(2)),
        Group {
            rank: 0,
            torsion: vec![]
        }
    );
    assert_eq!(
        abel(&circle()),
        Group {
            rank: 1,
            torsion: vec![]
        }
    );
    assert_eq!(
        abel(&product(&circle(), &circle(), 100_000).unwrap()),
        Group {
            rank: 2,
            torsion: vec![]
        }
    );
    assert_eq!(
        abel(&projective_plane(0)),
        Group {
            rank: 0,
            torsion: bi(&[2])
        }
    );
    assert_eq!(
        abel(&klein_bottle(4, 4)),
        Group {
            rank: 1,
            torsion: bi(&[2])
        }
    );
    assert_eq!(
        abel(&wedge(&circle(), &circle(), 0, 0)),
        Group {
            rank: 2,
            torsion: vec![]
        }
    );
    assert_eq!(
        abel(&dunce_cap(3)),
        Group {
            rank: 0,
            torsion: vec![]
        }
    );
    assert!(pi1_presentation(&[]).is_none());
    // H_1 is the abelianisation of pi_1 (Hurewicz) for every connected space above
    for c in [sphere(2), projective_plane(0), klein_bottle(4, 4), circle()] {
        let h1 = complex_chains(&c).homology(1);
        assert_eq!(abel(&c), h1);
    }
    // a presentation with a trivial-word relator and a power relator
    let p = Presentation {
        count: 2,
        relators: vec![vec![1, 1, 1], vec![2, -2]],
        edges: vec![],
    };
    assert_eq!(
        abelianization(&p),
        Group {
            rank: 1,
            torsion: bi(&[3])
        }
    );
    let s = simplify_presentation(&p);
    assert_eq!(
        abelianization(&s),
        Group {
            rank: 1,
            torsion: bi(&[3])
        }
    );
    assert_eq!(
        abelianization(&Presentation {
            count: 0,
            relators: vec![],
            edges: vec![]
        }),
        Group {
            rank: 0,
            torsion: vec![]
        }
    );
}

#[test]
fn simplicial_maps_and_induced_homology() {
    let disc = full_simplex(2);
    let id: HashMap<usize, usize> = (0..3).map(|v| (v, v)).collect();
    assert!(is_simplicial_map(&disc, &disc, &id));
    let collapse: HashMap<usize, usize> = [(0, 0), (1, 1), (2, 1)].into_iter().collect();
    assert!(is_simplicial_map(&disc, &disc, &collapse));
    let missing: HashMap<usize, usize> = [(0, 0), (1, 1), (2, 9)].into_iter().collect();
    assert!(!is_simplicial_map(&disc, &disc, &missing));
    assert!(induced_homology(&disc, &disc, &missing, 0).is_none());
    // identity on the torus induces the identity on H_1 and H_2
    let torus = product(&circle(), &circle(), 100_000).unwrap();
    let vertices: Vec<usize> = torus
        .iter()
        .filter(|s| s.len() == 1)
        .map(|s| s[0])
        .collect();
    let tid: HashMap<usize, usize> = vertices.iter().map(|&v| (v, v)).collect();
    for (dim, size) in [(1usize, 2usize), (2, 1)] {
        let m = induced_homology(&torus, &torus, &tid, dim).unwrap();
        assert_eq!(m, qlinalg::identity(size), "H_{dim}");
    }
    // the chain map of the identity is the identity matrix
    let cm = induced_chain_map(&disc, &disc, &id, 1).unwrap();
    for (i, row) in cm.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            assert_eq!(v, i64::from(i == j));
        }
    }
}

#[test]
fn mayer_vietoris_on_the_two_sphere() {
    // S^2 = boundary of the tetrahedron, as two pairs of faces meeting in a 4-cycle
    let a = close(&[vec![0, 1, 2], vec![0, 1, 3]]);
    let b = close(&[vec![0, 2, 3], vec![1, 2, 3]]);
    assert!(mayer_vietoris(&a, &b).is_some());
    assert_eq!(betti_list(&union(&a, &b)), vec![1, 0, 1]);
    assert_eq!(betti_list(&intersection(&a, &b)), vec![1, 1]);
    // the circle as two arcs meeting in two points
    let c1 = close(&[vec![0, 1], vec![1, 2]]);
    let c2 = close(&[vec![0, 2]]);
    assert!(mayer_vietoris(&c1, &c2).is_some());
}

fn close(s: &[Simplex]) -> Vec<Simplex> {
    rssn::kernels::topology::close_complex(s)
}

#[test]
fn persistent_homology_of_a_circle_sample() {
    let n = 8;
    let pts: Vec<Vec<f64>> = (0..n)
        .map(|k| {
            let t = std::f64::consts::TAU * k as f64 / n as f64;
            vec![t.cos(), t.sin()]
        })
        .collect();
    let d = distance_matrix(&pts);
    for i in 0..n {
        assert_eq!(d[i][i], 0.0);
        for j in 0..n {
            assert_eq!(d[i][j], d[j][i]);
        }
    }
    let edge = 2.0 * (std::f64::consts::PI / n as f64).sin();
    assert!((d[0][1] - edge).abs() < 1e-12);
    let filt = rips_filtration(&d, 3.0, 3, 1_000_000).unwrap();
    let bars = persistence_bars(&filt, 1).unwrap();
    let dim0: Vec<_> = bars.iter().filter(|b| b.dim == 0).collect();
    let dim1: Vec<_> = bars.iter().filter(|b| b.dim == 1).collect();
    assert_eq!(dim0.iter().filter(|b| b.death.is_infinite()).count(), 1);
    assert_eq!(
        dim0.len(),
        n - 1 + 1 - dim0.iter().filter(|b| b.birth == b.death).count()
    );
    for b in dim0.iter().filter(|b| b.death.is_finite()) {
        assert!((b.death - edge).abs() < 1e-12 && b.birth == 0.0);
    }
    assert_eq!(dim1.len(), 1, "one loop: {dim1:?}");
    assert!((dim1[0].birth - edge).abs() < 1e-12);
    assert!(dim1[0].death > dim1[0].birth + 0.3);
    // bars are sorted by (dim, birth, death) and have positive length
    assert!(bars.iter().all(|b| b.death > b.birth));
    // limit and empty input
    assert!(rips_filtration(&d, 3.0, 3, 10).is_none());
    assert_eq!(persistence_bars(&[], 2), Some(vec![]));
    // a filtration with a missing face is rejected rather than panicking
    let _ = persistence_bars(&[(0.0, vec![0, 1])], 1);
}

#[test]
fn bottleneck_distance_is_a_metric() {
    let a = [(0.0, 2.0), (1.0, 4.0)];
    let b = [(0.0, 2.5)];
    assert_eq!(bottleneck_distance(&a, &a), 0.0);
    assert_eq!(bottleneck_distance(&[], &[]), 0.0);
    assert!(
        (bottleneck_distance(&[(0.0, 2.0)], &[]) - 1.0).abs() < 1e-12,
        "matched to the diagonal"
    );
    assert!((bottleneck_distance(&[(0.0, 1.0)], &[(0.0, 1.5)]) - 0.5).abs() < 1e-12);
    assert!((bottleneck_distance(&a, &b) - bottleneck_distance(&b, &a)).abs() < 1e-12);
    // infinite bars: counts must match, otherwise the distance is infinite
    assert_eq!(
        bottleneck_distance(&[(0.0, f64::INFINITY)], &[]),
        f64::INFINITY
    );
    assert!(
        (bottleneck_distance(&[(0.0, f64::INFINITY)], &[(0.3, f64::INFINITY)]) - 0.3).abs() < 1e-12
    );
    // triangle inequality on random small diagrams
    let mut rng = Lcg(21);
    let mut diagram = |k: usize| -> Vec<(f64, f64)> {
        (0..k)
            .map(|_| {
                let birth = rng.int(0, 100) as f64 / 20.0;
                (birth, birth + rng.int(0, 100) as f64 / 20.0)
            })
            .collect()
    };
    for _ in 0..30 {
        let (x, y, z) = (diagram(3), diagram(4), diagram(2));
        let (dxy, dyz, dxz) = (
            bottleneck_distance(&x, &y),
            bottleneck_distance(&y, &z),
            bottleneck_distance(&x, &z),
        );
        assert!(dxz <= dxy + dyz + 1e-12, "{dxz} > {dxy} + {dyz}");
        assert!((dxy - bottleneck_distance(&y, &x)).abs() < 1e-12);
        assert!(dxy >= 0.0);
    }
}

#[test]
fn cubical_betti_numbers_agree_with_smith_normal_form() {
    let ring = vec![
        vec![true, true, true],
        vec![true, false, true],
        vec![true, true, true],
    ];
    assert_eq!(cubical_betti(&ring, false), (1, 1));
    assert_eq!(cubical_betti(&[vec![true]], false), (1, 0));
    assert_eq!(cubical_betti(&[], false), (0, 0));
    assert_eq!(cubical_betti(&[vec![false, false]], true), (0, 0));
    // two pixels touching only at a corner
    let diag = vec![vec![true, false], vec![false, true]];
    assert_eq!(cubical_betti(&diag, true).0, 1);
    assert_eq!(cubical_betti(&diag, false).0, 2);
    // random images: union-find / Euler formula vs Smith normal form of the chain complex
    let mut rng = Lcg(99);
    for trial in 0..40 {
        let (h, w) = (rng.int(1, 6) as usize, rng.int(1, 6) as usize);
        let img: Vec<Vec<bool>> = (0..h)
            .map(|_| (0..w).map(|_| rng.next() % 100 < 60).collect())
            .collect();
        for eight in [false, true] {
            let (b0, b1) = cubical_betti(&img, eight);
            let cc = cubical_chains(&img, eight);
            assert_eq!(
                cc.homology(0).rank,
                b0,
                "trial {trial} eight={eight} {img:?}"
            );
            assert_eq!(
                cc.homology(1).rank,
                b1,
                "trial {trial} eight={eight} {img:?}"
            );
            assert_eq!(cc.homology(2).rank, 0, "planar");
            // chain complex axiom d_1 d_2 = 0
            if cc.dims[2] > 0 && cc.dims[1] > 0 {
                let (d1, d2) = (&cc.d[1], &cc.d[2]);
                for i in 0..cc.dims[0] {
                    for j in 0..cc.dims[2] {
                        let s: i64 = (0..cc.dims[1]).map(|k| d1[i][k] * d2[k][j]).sum();
                        assert_eq!(s, 0);
                    }
                }
            }
        }
    }
}
