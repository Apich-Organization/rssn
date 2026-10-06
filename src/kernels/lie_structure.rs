//! # Lie algebras: classification data and structure theory
//!
//! Two exact (rational and `BigInt`) toolkits.
//!
//! **Semisimple Lie algebras.** Cartan matrices of the types `A_n`, `B_n`,
//! `C_n`, `D_n`, `E_6`, `E_7`, `E_8`, `F_4`, `G_2` (Bourbaki numbering, with
//! `a_ij = <alpha_i, alpha_j^vee>`, so `B_n` has `-2` in position
//! `(n-1, n)` and `C_n` in `(n, n-1)`), the recognition of the type of a
//! Cartan matrix, positive roots in simple-root coordinates, the Weyl group
//! order and Coxeter numbers from the heights of the roots, fundamental
//! weights, the Weyl dimension formula, Freudenthal's multiplicity formula,
//! the Brauer-Klimyk tensor product, Casimir eigenvalues and Dynkin indices.
//! Weights are written by their Dynkin labels `<lambda, alpha_i^vee>`.
//!
//! **Structure theory.** A finite dimensional Lie algebra over `Q` given by
//! structure constants: derived, lower and upper central series, Killing
//! form, Cartan's criteria, radical, Levi decomposition, Cartan
//! subalgebras, root decomposition and Casimir operators.

use std::collections::BTreeMap;
use std::collections::HashSet;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use super::qlinalg as ql;
use super::qlinalg::Q;
use super::qlinalg::QMat;

// ----------------------------------------------------------------------
// Cartan matrices
// ----------------------------------------------------------------------

/// A Cartan matrix, `a[i][j] = <alpha_i, alpha_j^vee>`.
pub type Cartan = Vec<Vec<i64>>;

/// Largest rank accepted.
pub const MAX_RANK: usize = 64;

/// The Cartan matrix of type `letter` and rank `n`.
#[must_use]
pub fn cartan_matrix(
    letter: char,
    n: usize,
) -> Option<Cartan> {
    if n == 0 || n > MAX_RANK {
        return None;
    }
    let mut a = vec![vec![0_i64; n]; n];
    for (i, row) in a.iter_mut().enumerate() {
        row[i] = 2;
    }
    match letter {
        | 'A' => (0..n - 1).for_each(|i| link(&mut a, i, i + 1)),
        | 'B' | 'C' if n >= 2 => {
            (0..n - 1).for_each(|i| link(&mut a, i, i + 1));
            if letter == 'B' {
                a[n - 2][n - 1] = -2;
            } else {
                a[n - 1][n - 2] = -2;
            }
        },
        | 'D' if n >= 4 => {
            (0..n - 2).for_each(|i| link(&mut a, i, i + 1));
            link(&mut a, n - 3, n - 1);
        },
        | 'E' if (6..=8).contains(&n) => {
            link(&mut a, 0, 2);
            (2..n - 1).for_each(|i| link(&mut a, i, i + 1));
            link(&mut a, 1, 3);
        },
        | 'F' if n == 4 => {
            link(&mut a, 0, 1);
            a[1][2] = -2;
            a[2][1] = -1;
            link(&mut a, 2, 3);
        },
        | 'G' if n == 2 => {
            a[0][1] = -1;
            a[1][0] = -3;
        },
        | _ => return None,
    }
    Some(a)
}

fn link(
    a: &mut Cartan,
    i: usize,
    j: usize,
) {
    a[i][j] = -1;
    a[j][i] = -1;
}

/// Parses a type name such as `A2`, `E8`, `G2`.
#[must_use]
pub fn parse_type(name: &str) -> Option<Cartan> {
    let mut chars = name.trim().chars();
    let letter = chars.next()?.to_ascii_uppercase();
    let n: usize = chars.as_str().trim().parse().ok()?;
    cartan_matrix(letter, n)
}

/// The standard simple roots of a type in Euclidean coordinates.
///
/// Those of Bourbaki, not normalised to a common length: `A_n` lives in `R^{n+1}`,
/// the exceptional `E`-types in `R^8`, `F_4` in `R^4`, `G_2` in `R^3`.
#[must_use]
pub fn euclidean_simple_roots(
    letter: char,
    n: usize,
) -> Option<Vec<Vec<Q>>> {
    let half = Q::new(BigInt::one(), BigInt::from(2));
    let vec_of = |dim: usize, entries: &[(usize, Q)]| -> Vec<Q> {
        let mut v = ql::zeros(dim);
        for (i, x) in entries {
            v[*i] = x.clone();
        }
        v
    };
    let one = ql::q(1);
    let chain = |dim: usize, count: usize| -> Vec<Vec<Q>> {
        (0..count).map(|i| vec_of(dim, &[(i, ql::q(1)), (i + 1, ql::q(-1))])).collect()
    };
    match letter.to_ascii_uppercase() {
        | 'A' if n >= 1 => Some(chain(n + 1, n)),
        | 'B' if n >= 2 => {
            let mut r = chain(n, n - 1);
            r.push(vec_of(n, &[(n - 1, one)]));
            Some(r)
        },
        | 'C' if n >= 2 => {
            let mut r = chain(n, n - 1);
            r.push(vec_of(n, &[(n - 1, ql::q(2))]));
            Some(r)
        },
        | 'D' if n >= 4 => {
            let mut r = chain(n, n - 1);
            r.push(vec_of(n, &[(n - 2, one.clone()), (n - 1, one)]));
            Some(r)
        },
        | 'E' if (6..=8).contains(&n) => {
            let mut all = vec![
                vec![half.clone(), -half.clone(), -half.clone(), -half.clone(), -half.clone(), -half.clone(), -half.clone(), half],
                vec_of(8, &[(0, ql::q(1)), (1, ql::q(1))]),
            ];
            all.push(vec_of(8, &[(1, ql::q(1)), (0, ql::q(-1))]));
            for i in 2..7 {
                all.push(vec_of(8, &[(i, ql::q(1)), (i - 1, ql::q(-1))]));
            }
            all.truncate(n);
            Some(all)
        },
        | 'F' if n == 4 => Some(vec![
            vec_of(4, &[(1, ql::q(1)), (2, ql::q(-1))]),
            vec_of(4, &[(2, ql::q(1)), (3, ql::q(-1))]),
            vec_of(4, &[(3, ql::q(1))]),
            vec_of(4, &[(0, half.clone()), (1, -half.clone()), (2, -half.clone()), (3, -half)]),
        ]),
        | 'G' if n == 2 => Some(vec![vec_of(3, &[(0, ql::q(1)), (1, ql::q(-1))]), vec_of(3, &[(0, ql::q(-2)), (1, ql::q(1)), (2, ql::q(1))])]),
        | _ => None,
    }
}

fn is_cartan_shape(a: &Cartan) -> bool {
    let n = a.len();
    n > 0
        && n <= MAX_RANK
        && a.iter().all(|r| r.len() == n)
        && (0..n).all(|i| {
            a[i][i] == 2
                && (0..n).all(|j| i == j || (a[i][j] <= 0 && (a[i][j] == 0) == (a[j][i] == 0) && a[i][j] * a[j][i] <= 3))
        })
}

/// The connected components (as index lists) of the Dynkin diagram.
#[must_use]
pub fn components(a: &Cartan) -> Vec<Vec<usize>> {
    let n = a.len();
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut comp = vec![start];
        seen[start] = true;
        let mut at = 0;
        while at < comp.len() {
            let u = comp[at];
            at += 1;
            for v in 0..n {
                if !seen[v] && a[u][v] != 0 {
                    seen[v] = true;
                    comp.push(v);
                }
            }
        }
        comp.sort_unstable();
        out.push(comp);
    }
    out
}

struct Search<'a> {
    a: &'a Cartan,
    b: &'a Cartan,
    order: Vec<usize>,
    parent: Vec<usize>,
    image: Vec<usize>,
    used: Vec<bool>,
}
fn go(
    s: &mut Search<'_>,
    k: usize,
) -> bool {
    let n = s.a.len();
    if k == n {
        return true;
    }
    let v = s.order[k];
    let candidates: Vec<usize> = if k == 0 {
        (0..n).collect()
    } else {
        let anchor = s.image[s.parent[v]];
        (0..n).filter(|&w| !s.used[w] && s.b[anchor][w] != 0).collect()
    };
    for w in candidates {
        let consistent = s.order[..k].iter().all(|&u| s.a[v][u] == s.b[w][s.image[u]] && s.a[u][v] == s.b[s.image[u]][w]);
        if consistent {
            s.image[v] = w;
            s.used[w] = true;
            if go(s, k + 1) {
                return true;
            }
            s.used[w] = false;
        }
    }
    false
}

fn isomorphic(
    a: &Cartan,
    b: &Cartan,
) -> bool {
    let n = a.len();
    if n != b.len() {
        return false;
    }
    let mut order = vec![0_usize];
    let mut parent = vec![usize::MAX; n];
    let mut seen = vec![false; n];
    seen[0] = true;
    let mut at = 0;
    while at < order.len() {
        let u = order[at];
        at += 1;
        for v in 0..n {
            if !seen[v] && a[u][v] != 0 {
                seen[v] = true;
                parent[v] = u;
                order.push(v);
            }
        }
    }
    if order.len() != n {
        return false;
    }
    let mut s = Search { a, b, order, parent, image: vec![usize::MAX; n], used: vec![false; n] };
    go(&mut s, 0)
}

fn sub_matrix(
    a: &Cartan,
    idx: &[usize],
) -> Cartan {
    idx.iter().map(|&i| idx.iter().map(|&j| a[i][j]).collect()).collect()
}

/// The simple factors of the algebra with Cartan matrix `a` as
/// `(letter, rank)` pairs, or `None` if `a` is not the Cartan matrix of a
/// semisimple Lie algebra.
#[must_use]
pub fn recognize(a: &Cartan) -> Option<Vec<(char, usize)>> {
    if !is_cartan_shape(a) {
        return None;
    }
    let mut out = Vec::new();
    for comp in components(a) {
        let sub = sub_matrix(a, &comp);
        let n = comp.len();
        let letters: &[char] = &['A', 'B', 'C', 'D', 'E', 'F', 'G'];
        let found = letters.iter().find(|&&l| cartan_matrix(l, n).is_some_and(|c| isomorphic(&sub, &c)))?;
        out.push((*found, n));
    }
    Some(out)
}

// ----------------------------------------------------------------------
// Root systems and representations
// ----------------------------------------------------------------------

/// The root system of a semisimple Lie algebra.
#[derive(Debug, Clone)]
pub struct RootSystem {
    /// The rank.
    pub n: usize,
    /// The Cartan matrix.
    pub cartan: Cartan,
    /// The positive roots in simple-root coordinates, by height.
    pub positive: Vec<Vec<i64>>,
    /// The squared lengths of the simple roots (the longest root of every
    /// simple factor has squared length 2).
    pub lengths: Vec<Q>,
    inverse: QMat,
}

impl RootSystem {
    /// The root system of a Cartan matrix of finite type.
    #[must_use]
    pub fn new(cartan: Cartan) -> Option<Self> {
        recognize(&cartan)?;
        let n = cartan.len();
        let positive = positive_roots(&cartan);
        let inverse = ql::inverse(&cartan.iter().map(|r| r.iter().map(|&x| ql::q(x)).collect()).collect::<Vec<Vec<Q>>>())?;
        let lengths = root_lengths(&cartan);
        Some(Self { n, cartan, positive, lengths, inverse })
    }

    /// The number of positive roots.
    #[must_use]
    pub const fn num_positive(&self) -> usize {
        self.positive.len()
    }

    /// The dimension `n + 2 N` of the Lie algebra.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.n + 2 * self.positive.len()
    }

    /// Whether the Dynkin diagram is connected.
    #[must_use]
    pub fn is_irreducible(&self) -> bool {
        components(&self.cartan).len() == 1
    }

    /// The Gram matrix `(alpha_i, alpha_j)` of the simple roots.
    #[must_use]
    pub fn gram(&self) -> QMat {
        (0..self.n)
            .map(|i| (0..self.n).map(|j| ql::q(self.cartan[i][j]) * &self.lengths[j] / ql::q(2)).collect())
            .collect()
    }

    /// The fundamental weights in simple-root coordinates (rows of the
    /// inverse Cartan matrix).
    #[must_use]
    pub const fn fundamental_weights(&self) -> &QMat {
        &self.inverse
    }

    /// The exponents: the dual partition of the numbers of roots of each height.
    #[must_use]
    pub fn exponents(&self) -> Vec<usize> {
        let heights: Vec<i64> = self.positive.iter().map(|r| r.iter().sum()).collect();
        let top = heights.iter().copied().max().unwrap_or(0);
        let counts: Vec<usize> = (1..=top).map(|h| heights.iter().filter(|&&x| x == h).count()).collect();
        let mut out: Vec<usize> = (1..=self.n).map(|j| counts.iter().filter(|&&c| c >= j).count()).collect();
        out.sort_unstable();
        out
    }

    /// The order of the Weyl group, `prod (e_i + 1)`.
    #[must_use]
    pub fn weyl_order(&self) -> BigInt {
        self.exponents().iter().fold(BigInt::one(), |acc, &e| acc * BigInt::from(e + 1))
    }

    /// The highest root (irreducible systems only).
    #[must_use]
    pub fn highest_root(&self) -> Option<Vec<i64>> {
        self.is_irreducible().then(|| self.positive.iter().max_by_key(|r| r.iter().sum::<i64>()).cloned()).flatten()
    }

    /// The Coxeter number `2 N / n` (irreducible systems only).
    #[must_use]
    pub fn coxeter_number(&self) -> Option<usize> {
        self.is_irreducible().then(|| 2 * self.positive.len() / self.n)
    }

    /// The dual Coxeter number `1 + (theta, rho)` (irreducible systems only).
    #[must_use]
    pub fn dual_coxeter_number(&self) -> Option<Q> {
        let theta = self.highest_root()?;
        let sum = theta.iter().zip(&self.lengths).fold(Q::zero(), |acc, (&c, s)| acc + ql::q(c) * s / ql::q(2));
        Some(sum + ql::q(1))
    }

    /// The Weyl dimension formula for the highest weight with Dynkin
    /// labels `labels` (nonnegative integers).
    #[must_use]
    pub fn weyl_dim(
        &self,
        labels: &[i64],
    ) -> Option<BigInt> {
        if labels.len() != self.n || labels.iter().any(|&l| l < 0) {
            return None;
        }
        let mut value = Q::one();
        for root in &self.positive {
            let (mut num, mut den) = (Q::zero(), Q::zero());
            for ((&c, s), &l) in root.iter().zip(&self.lengths).zip(labels) {
                num += ql::q(c) * s * ql::q(l + 1);
                den += ql::q(c) * s;
            }
            value = value * num / den;
        }
        value.is_integer().then(|| value.to_integer())
    }

    /// The inner product of two weights given by Dynkin labels.
    #[must_use]
    pub fn inner(
        &self,
        a: &[i64],
        b: &[i64],
    ) -> Q {
        let mut total = Q::zero();
        for (k, (&ak, length)) in a.iter().zip(&self.lengths).enumerate() {
            let x: Q = (0..self.n).fold(Q::zero(), |acc, i| acc + ql::q(b[i]) * &self.inverse[i][k]);
            total += x * ql::q(ak) * length / ql::q(2);
        }
        total
    }

    /// The Casimir eigenvalue `(L, L + 2 rho)` of the irreducible
    /// representation of highest weight `L` (long roots of squared length 2).
    #[must_use]
    pub fn casimir(
        &self,
        labels: &[i64],
    ) -> Option<Q> {
        if labels.len() != self.n || !self.is_irreducible() {
            return None;
        }
        let shifted: Vec<i64> = labels.iter().map(|&l| l + 2).collect();
        Some(self.inner(labels, &shifted))
    }

    /// The Dynkin index `dim(V) C(V) / dim(g)` (the fundamental of `su(n)`
    /// has index 1, the adjoint `2 h^vee`).
    #[must_use]
    pub fn dynkin_index(
        &self,
        labels: &[i64],
    ) -> Option<Q> {
        let dim = Q::from_integer(self.weyl_dim(labels)?);
        Some(dim * self.casimir(labels)? / ql::q(i64::try_from(self.dimension()).ok()?))
    }

    fn reflect_dominant(
        &self,
        w: &[i64],
    ) -> (Vec<i64>, bool) {
        let mut w = w.to_vec();
        let mut odd = false;
        while let Some(i) = w.iter().position(|&x| x < 0) {
            let f = w[i];
            for (x, a) in w.iter_mut().zip(&self.cartan[i]) {
                *x -= f * a;
            }
            odd = !odd;
        }
        (w, odd)
    }

    fn root_data(&self) -> Vec<(Vec<i64>, Vec<i64>, Q)> {
        let gram = self.gram();
        self.positive
            .iter()
            .map(|c| {
                let labels: Vec<i64> = (0..self.n).map(|i| (0..self.n).map(|j| c[j] * self.cartan[j][i]).sum()).collect();
                let mut norm = Q::zero();
                for i in 0..self.n {
                    for j in 0..self.n {
                        norm += ql::q(c[i] * c[j]) * &gram[i][j];
                    }
                }
                (c.clone(), labels, norm)
            })
            .collect()
    }

    fn height_below(
        &self,
        top: &[i64],
        w: &[i64],
    ) -> i64 {
        let diff: Vec<i64> = top.iter().zip(w).map(|(a, b)| a - b).collect();
        let mut total = Q::zero();
        for k in 0..self.n {
            for (i, d) in diff.iter().enumerate() {
                total += ql::q(*d) * &self.inverse[i][k];
            }
        }
        total.to_integer().to_i64().unwrap_or(0)
    }

    /// The multiplicities of the dominant weights of the irreducible
    /// representation of highest weight `top`, by Freudenthal's formula;
    /// `None` for reducible systems or when more than `cap` dominant
    /// weights occur.
    #[must_use]
    pub fn dominant_multiplicities(
        &self,
        top: &[i64],
        cap: usize,
    ) -> Option<BTreeMap<Vec<i64>, BigInt>> {
        if top.len() != self.n || top.iter().any(|&l| l < 0) || !self.is_irreducible() {
            return None;
        }
        let data = self.root_data();
        let mut seen: HashSet<Vec<i64>> = HashSet::from([top.to_vec()]);
        let mut frontier = vec![top.to_vec()];
        while let Some(mu) = frontier.pop() {
            for (_, labels, _) in &data {
                let next: Vec<i64> = mu.iter().zip(labels).map(|(a, b)| a - b).collect();
                if next.iter().all(|&x| x >= 0) && seen.insert(next.clone()) {
                    if seen.len() > cap {
                        return None;
                    }
                    frontier.push(next);
                }
            }
        }
        let mut order: Vec<(i64, Vec<i64>)> = seen.into_iter().map(|w| (self.height_below(top, &w), w)).collect();
        order.sort();
        let two_rho: Vec<i64> = vec![2; self.n];
        let casimir = |w: &[i64]| -> Q {
            let shifted: Vec<i64> = w.iter().zip(&two_rho).map(|(a, b)| a + b).collect();
            self.inner(w, &shifted)
        };
        let top_casimir = casimir(top);
        let mut mult: BTreeMap<Vec<i64>, BigInt> = BTreeMap::new();
        for (_, mu) in order {
            if mu == top {
                mult.insert(mu, BigInt::one());
                continue;
            }
            let mut num = Q::zero();
            for (c, labels, norm) in &data {
                let pairing: Q = c.iter().zip(&mu).zip(&self.lengths).fold(Q::zero(), |acc, ((&ci, &li), s)| acc + ql::q(ci * li) * s / ql::q(2));
                for k in 1_i64.. {
                    let moved: Vec<i64> = mu.iter().zip(labels).map(|(a, b)| a + k * b).collect();
                    let (dominant, _) = self.reflect_dominant(&moved);
                    let Some(m) = mult.get(&dominant) else { break };
                    num += Q::from_integer(m.clone()) * (&pairing + ql::q(k) * norm);
                }
            }
            let value = num * ql::q(2) / (&top_casimir - casimir(&mu));
            if !value.is_integer() || value.is_negative() {
                return None;
            }
            mult.insert(mu, value.to_integer());
        }
        Some(mult)
    }

    /// All weights with multiplicity (dominant weights and their Weyl
    /// orbits); `None` if there are more than `cap` of them.
    #[must_use]
    pub fn all_weights(
        &self,
        top: &[i64],
        cap: usize,
    ) -> Option<BTreeMap<Vec<i64>, BigInt>> {
        let dominant = self.dominant_multiplicities(top, cap)?;
        let mut out = BTreeMap::new();
        for (d, m) in dominant {
            let mut orbit: HashSet<Vec<i64>> = HashSet::from([d.clone()]);
            let mut stack = vec![d];
            while let Some(w) = stack.pop() {
                for i in 0..self.n {
                    let f = w[i];
                    if f == 0 {
                        continue;
                    }
                    let next: Vec<i64> = w.iter().zip(&self.cartan[i]).map(|(x, a)| x - f * a).collect();
                    if orbit.insert(next.clone()) {
                        if out.len() + orbit.len() > cap {
                            return None;
                        }
                        stack.push(next);
                    }
                }
            }
            for w in orbit {
                out.insert(w, m.clone());
            }
        }
        Some(out)
    }

    /// The multiplicity of the weight `w` (any weight) in the irreducible
    /// representation of highest weight `top`.
    #[must_use]
    pub fn weight_multiplicity(
        &self,
        top: &[i64],
        w: &[i64],
    ) -> Option<BigInt> {
        if w.len() != self.n {
            return None;
        }
        let (dominant, _) = self.reflect_dominant(w);
        let table = self.dominant_multiplicities(top, 100_000)?;
        Some(table.get(&dominant).cloned().unwrap_or_else(BigInt::zero))
    }

    /// The decomposition of `V(a) ⊗ V(b)` into irreducibles (Brauer-Klimyk):
    /// `(highest weight, multiplicity)` pairs.
    #[must_use]
    pub fn tensor(
        &self,
        a: &[i64],
        b: &[i64],
    ) -> Option<Vec<(Vec<i64>, BigInt)>> {
        let (da, db) = (self.weyl_dim(a)?, self.weyl_dim(b)?);
        let (big, small) = if da >= db { (a, b) } else { (b, a) };
        let weights = self.all_weights(small, 200_000)?;
        let mut out: BTreeMap<Vec<i64>, BigInt> = BTreeMap::new();
        for (nu, m) in weights {
            let shifted: Vec<i64> = big.iter().zip(&nu).map(|(x, y)| x + y + 1).collect();
            let (dominant, odd) = self.reflect_dominant(&shifted);
            if dominant.contains(&0) {
                continue;
            }
            let highest: Vec<i64> = dominant.iter().map(|x| x - 1).collect();
            let entry = out.entry(highest).or_insert_with(BigInt::zero);
            if odd {
                *entry -= m;
            } else {
                *entry += m;
            }
        }
        Some(out.into_iter().filter(|(_, m)| !m.is_zero()).collect())
    }
}

/// The positive roots of the root system with Cartan matrix `a`, in
/// simple-root coordinates, ordered by height (the root strings through the
/// simple roots are followed level by level).
#[must_use]
pub fn positive_roots(a: &Cartan) -> Vec<Vec<i64>> {
    let n = a.len();
    let simple: Vec<Vec<i64>> = (0..n).map(|i| (0..n).map(|j| i64::from(i == j)).collect()).collect();
    let mut all: HashSet<Vec<i64>> = simple.iter().cloned().collect();
    let mut level = simple;
    let mut out: Vec<Vec<i64>> = Vec::new();
    while !level.is_empty() {
        let mut next: Vec<Vec<i64>> = Vec::new();
        for beta in &level {
            for i in 0..n {
                let pairing: i64 = (0..n).map(|j| beta[j] * a[j][i]).sum();
                let mut p = 0_i64;
                let mut probe = beta.clone();
                loop {
                    probe[i] -= 1;
                    if all.contains(&probe) {
                        p += 1;
                    } else {
                        break;
                    }
                }
                if p - pairing > 0 {
                    let mut up = beta.clone();
                    up[i] += 1;
                    if all.insert(up.clone()) {
                        next.push(up);
                    }
                }
            }
        }
        out.extend(level);
        next.sort();
        level = next;
    }
    out
}

fn root_lengths(a: &Cartan) -> Vec<Q> {
    let n = a.len();
    let mut s = vec![Q::zero(); n];
    for comp in components(a) {
        s[comp[0]] = Q::one();
        let mut queue = vec![comp[0]];
        let mut at = 0;
        while at < queue.len() {
            let i = queue[at];
            at += 1;
            for &j in &comp {
                if a[i][j] != 0 && s[j].is_zero() && i != j {
                    s[j] = s[i].clone() * ql::q(a[j][i]) / ql::q(a[i][j]);
                    queue.push(j);
                }
            }
        }
        let longest = comp.iter().map(|&i| s[i].clone()).max().unwrap_or_else(Q::one);
        for &i in &comp {
            s[i] = s[i].clone() * ql::q(2) / &longest;
        }
    }
    s
}

// ----------------------------------------------------------------------
// Standard bases
// ----------------------------------------------------------------------

/// A Gaussian-integer matrix: entries `(re, im)`.
pub type GMat = Vec<Vec<(i64, i64)>>;

fn zero_g(n: usize) -> GMat {
    vec![vec![(0, 0); n]; n]
}

/// The elementary matrix with a one at `(i, j)`.
fn unit(
    n: usize,
    i: usize,
    j: usize,
    value: (i64, i64),
) -> GMat {
    let mut m = zero_g(n);
    m[i][j] = value;
    m
}

fn add_g(
    a: &GMat,
    b: &GMat,
) -> GMat {
    a.iter().zip(b).map(|(r, s)| r.iter().zip(s).map(|(x, y)| (x.0 + y.0, x.1 + y.1)).collect()).collect()
}

/// The matrix algebras whose standard bases are available.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// `gl(n)`.
    Gl,
    /// `sl(n)`.
    Sl,
    /// `so(n)`.
    So,
    /// `sp(2n)` (the argument is `n`).
    Sp,
    /// `su(n)`.
    Su,
}

/// The standard basis of a classical matrix Lie algebra.
///
/// Matrices of
/// Gaussian integers (`gl`, `sl`: elementary matrices and the Chevalley
/// `H_k = E_kk - E_{k+1,k+1}`; `so(n)`: `E_ij - E_ji`; `sp(2n)`:
/// `[[A, B], [C, -A^T]]` with `B`, `C` symmetric; `su(n)`: `E_ij - E_ji`,
/// `i (E_ij + E_ji)`, `i H_k`).
#[must_use]
pub fn standard_basis(
    family: Family,
    n: usize,
) -> Option<Vec<GMat>> {
    if n == 0 || n > 32 {
        return None;
    }
    let one = (1, 0);
    let mut out = Vec::new();
    match family {
        | Family::Gl => {
            for i in 0..n {
                for j in 0..n {
                    out.push(unit(n, i, j, one));
                }
            }
        },
        | Family::Sl | Family::Su => {
            if n < 2 {
                return None;
            }
            let su = family == Family::Su;
            for i in 0..n {
                for j in i + 1..n {
                    if su {
                        out.push(add_g(&unit(n, i, j, one), &unit(n, j, i, (-1, 0))));
                        out.push(add_g(&unit(n, i, j, (0, 1)), &unit(n, j, i, (0, 1))));
                    } else {
                        out.push(unit(n, i, j, one));
                    }
                }
            }
            if !su {
                for i in 0..n {
                    for j in 0..i {
                        out.push(unit(n, i, j, one));
                    }
                }
            }
            for k in 0..n - 1 {
                let h = add_g(&unit(n, k, k, one), &unit(n, k + 1, k + 1, (-1, 0)));
                out.push(if su { h.iter().map(|r| r.iter().map(|x| (-x.1, x.0)).collect()).collect() } else { h });
            }
        },
        | Family::So => {
            if n < 2 {
                return None;
            }
            for i in 0..n {
                for j in i + 1..n {
                    out.push(add_g(&unit(n, i, j, one), &unit(n, j, i, (-1, 0))));
                }
            }
        },
        | Family::Sp => {
            let m = 2 * n;
            for i in 0..n {
                for j in 0..n {
                    out.push(add_g(&unit(m, i, j, one), &unit(m, n + j, n + i, (-1, 0))));
                }
            }
            for i in 0..n {
                for j in i..n {
                    let mut b = unit(m, i, n + j, one);
                    if i != j {
                        b = add_g(&b, &unit(m, j, n + i, one));
                    }
                    out.push(b);
                }
            }
            for i in 0..n {
                for j in i..n {
                    let mut c = unit(m, n + i, j, one);
                    if i != j {
                        c = add_g(&c, &unit(m, n + j, i, one));
                    }
                    out.push(c);
                }
            }
        },
    }
    Some(out)
}

// ----------------------------------------------------------------------
// Structure theory
// ----------------------------------------------------------------------

/// A finite dimensional Lie algebra over `Q` with structure constants
/// `[e_i, e_j] = sum_k c[i][j][k] e_k`.
#[derive(Debug, Clone)]
pub struct Lie {
    /// The dimension.
    pub n: usize,
    /// The structure constants.
    pub c: Vec<Vec<Vec<Q>>>,
}

/// A root space: the values of the root on the Cartan basis and a basis
/// of the space (coordinate vectors).
#[derive(Debug, Clone)]
pub struct RootSpace {
    /// The root, as its values on the basis of the Cartan subalgebra.
    pub root: Vec<Q>,
    /// A basis of the root space.
    pub space: QMat,
}

fn flatten(m: &[Vec<Q>]) -> Vec<Q> {
    m.iter().flatten().cloned().collect()
}

impl Lie {
    /// The Lie algebra spanned by a basis of rational matrices; `None`
    /// unless they are independent and closed under the commutator.
    #[must_use]
    pub fn from_matrices(basis: &[QMat]) -> Option<Self> {
        let n = basis.len();
        let flat: Vec<Vec<Q>> = basis.iter().map(|m| flatten(m)).collect();
        let len = flat.first()?.len();
        if ql::rank(&flat, len) != n {
            return None;
        }
        let system = ql::transpose(&flat, len);
        let mut c = vec![vec![Vec::new(); n]; n];
        for i in 0..n {
            for j in 0..n {
                let ab = ql::matmul(&basis[i], &basis[j]);
                let ba = ql::matmul(&basis[j], &basis[i]);
                let diff: Vec<Q> = flatten(&ab).into_iter().zip(flatten(&ba)).map(|(x, y)| x - y).collect();
                c[i][j] = ql::solve(&system, n, &diff)?;
                // the solution must reproduce the commutator exactly
                let back = ql::matvec(&system, &c[i][j]);
                if back != diff {
                    return None;
                }
            }
        }
        Some(Self { n, c })
    }

    /// A Lie algebra from structure constants; `None` unless the bracket
    /// is antisymmetric and satisfies the Jacobi identity.
    #[must_use]
    pub fn from_constants(c: Vec<Vec<Vec<Q>>>) -> Option<Self> {
        let n = c.len();
        if n == 0 || c.iter().any(|p| p.len() != n || p.iter().any(|r| r.len() != n)) {
            return None;
        }
        let lie = Self { n, c };
        let unit = ql::identity(n);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    if lie.c[i][j][k] != -lie.c[j][i][k].clone() {
                        return None;
                    }
                }
            }
        }
        for i in 0..n {
            for j in i + 1..n {
                for k in j + 1..n {
                    let mut total = ql::zeros(n);
                    for (a, b, d) in [(i, j, k), (j, k, i), (k, i, j)] {
                        let inner = lie.bracket(&unit[b], &unit[d]);
                        let outer = lie.bracket(&unit[a], &inner);
                        for (t, o) in total.iter_mut().zip(outer) {
                            *t += o;
                        }
                    }
                    if !ql::is_zero_vec(&total) {
                        return None;
                    }
                }
            }
        }
        Some(lie)
    }

    /// `[x, y]` in coordinates.
    #[must_use]
    pub fn bracket(
        &self,
        x: &[Q],
        y: &[Q],
    ) -> Vec<Q> {
        let mut out = ql::zeros(self.n);
        for (i, xi) in x.iter().enumerate() {
            if xi.is_zero() {
                continue;
            }
            for (j, yj) in y.iter().enumerate() {
                if yj.is_zero() {
                    continue;
                }
                let f = xi * yj;
                for (o, cijk) in out.iter_mut().zip(&self.c[i][j]) {
                    *o += &f * cijk;
                }
            }
        }
        out
    }

    /// The matrix of `ad x`: column `j` holds the coordinates of `[x, e_j]`.
    #[must_use]
    pub fn ad(
        &self,
        x: &[Q],
    ) -> QMat {
        let unit = ql::identity(self.n);
        let columns: Vec<Vec<Q>> = unit.iter().map(|e| self.bracket(x, e)).collect();
        ql::transpose(&columns, self.n)
    }

    /// The Killing form `K_ab = tr(ad e_a ad e_b)`.
    #[must_use]
    pub fn killing(&self) -> QMat {
        let n = self.n;
        (0..n)
            .map(|a| {
                (0..n)
                    .map(|b| {
                        let mut s = Q::zero();
                        for j in 0..n {
                            for k in 0..n {
                                s += &self.c[a][j][k] * &self.c[b][k][j];
                            }
                        }
                        s
                    })
                    .collect()
            })
            .collect()
    }

    /// The span of all `[a_i, b_j]` (echelon basis).
    #[must_use]
    pub fn bracket_spaces(
        &self,
        a: &[Vec<Q>],
        b: &[Vec<Q>],
    ) -> QMat {
        let all: Vec<Vec<Q>> = a.iter().flat_map(|x| b.iter().map(move |y| (x, y))).map(|(x, y)| self.bracket(x, y)).collect();
        ql::span_basis(&all, self.n)
    }

    /// The derived series `g ⊇ [g, g] ⊇ ...` up to the point it stabilises.
    #[must_use]
    pub fn derived_series(&self) -> Vec<QMat> {
        let mut cur = ql::identity(self.n);
        let mut out = vec![cur.clone()];
        loop {
            let next = self.bracket_spaces(&cur, &cur);
            if next.len() == cur.len() {
                break;
            }
            out.push(next.clone());
            if next.is_empty() {
                break;
            }
            cur = next;
        }
        out
    }

    /// The lower central series `g ⊇ [g, g] ⊇ [g, [g, g]] ⊇ ...`.
    #[must_use]
    pub fn lower_central_series(&self) -> Vec<QMat> {
        let all = ql::identity(self.n);
        let mut cur = all.clone();
        let mut out = vec![cur.clone()];
        loop {
            let next = self.bracket_spaces(&all, &cur);
            if next.len() == cur.len() {
                break;
            }
            out.push(next.clone());
            if next.is_empty() {
                break;
            }
            cur = next;
        }
        out
    }

    /// The upper central series `0 ⊆ Z_1 ⊆ Z_2 ⊆ ...` (without the leading 0).
    #[must_use]
    pub fn upper_central_series(&self) -> Vec<QMat> {
        let n = self.n;
        let mut out: Vec<QMat> = Vec::new();
        let mut cur: QMat = Vec::new();
        loop {
            let mut echelon = cur.clone();
            let pivots = ql::rref(&mut echelon, n);
            // x -> ([x, e_j] mod Z)_j, as a matrix with columns indexed by x_i
            let mut columns: Vec<Vec<Q>> = Vec::with_capacity(n);
            for i in 0..n {
                let mut col = Vec::with_capacity(n * n);
                for j in 0..n {
                    col.extend(ql::reduce_mod(&self.c[i][j], &echelon, &pivots));
                }
                columns.push(col);
            }
            let system = ql::transpose(&columns, n * n);
            let next = ql::span_basis(&ql::nullspace(&system, n), n);
            if next.len() == cur.len() {
                break;
            }
            out.push(next.clone());
            cur = next;
        }
        out
    }

    /// Whether the derived series reaches zero.
    #[must_use]
    pub fn is_solvable(&self) -> bool {
        self.derived_series().last().is_some_and(Vec::is_empty)
    }

    /// Whether the lower central series reaches zero.
    #[must_use]
    pub fn is_nilpotent(&self) -> bool {
        self.lower_central_series().last().is_some_and(Vec::is_empty)
    }

    /// Cartan's first criterion: `g` is solvable iff `K(g, [g, g]) = 0`.
    #[must_use]
    pub fn cartan_solvable_test(&self) -> bool {
        let k = self.killing();
        let derived = self.bracket_spaces(&ql::identity(self.n), &ql::identity(self.n));
        derived.iter().all(|y| ql::is_zero_vec(&ql::matvec(&k, y)))
    }

    /// Cartan's second criterion: `g` is semisimple iff `K` is nondegenerate.
    #[must_use]
    pub fn is_semisimple(&self) -> bool {
        !ql::det(&self.killing()).is_zero()
    }

    /// The center.
    #[must_use]
    pub fn center(&self) -> QMat {
        self.upper_central_series().into_iter().next().unwrap_or_default()
    }

    /// The radical (largest solvable ideal) `= [g, g]^⊥` for the Killing form.
    #[must_use]
    pub fn radical(&self) -> QMat {
        let k = self.killing();
        let derived = self.bracket_spaces(&ql::identity(self.n), &ql::identity(self.n));
        let rows: Vec<Vec<Q>> = derived.iter().map(|y| ql::matvec(&k, y)).collect();
        ql::span_basis(&ql::nullspace(&rows, self.n), self.n)
    }

    /// Whether the span of `vectors` is closed under the bracket.
    #[must_use]
    pub fn is_subalgebra(
        &self,
        vectors: &[Vec<Q>],
    ) -> bool {
        let mut echelon = vectors.to_vec();
        let pivots = ql::rref(&mut echelon, self.n);
        vectors.iter().all(|x| vectors.iter().all(|y| ql::is_zero_vec(&ql::reduce_mod(&self.bracket(x, y), &echelon, &pivots))))
    }

    /// A Levi decomposition `g = s ⋉ r`: a semisimple subalgebra `s`
    /// complementary to the radical `r` (returned as `(s, r)`), lifted
    /// through the derived series of `r` (Whitehead's lemma makes each lift
    /// a linear problem).
    #[must_use]
    #[allow(clippy::needless_range_loop)]
    pub fn levi(&self) -> Option<(QMat, QMat)> {
        let n = self.n;
        let radical = self.radical();
        let mut echelon = radical.clone();
        let pivots = ql::rref(&mut echelon, n);
        let unit = ql::identity(n);
        let mut s: Vec<Vec<Q>> = (0..n).filter(|c| !pivots.contains(c)).map(|c| unit[c].clone()).collect();
        let p = s.len();
        if p == 0 {
            return Some((Vec::new(), radical));
        }
        // the structure constants of g / r in the basis s
        let complement_system = ql::transpose(&s.iter().chain(&radical).cloned().collect::<Vec<_>>(), n);
        let mut quotient = vec![vec![Vec::new(); p]; p];
        for i in 0..p {
            for j in 0..p {
                let x = ql::solve(&complement_system, p + radical.len(), &self.bracket(&s[i], &s[j]))?;
                quotient[i][j] = x[..p].to_vec();
            }
        }
        // the derived series of the radical
        let mut chain: Vec<QMat> = vec![radical];
        while let Some(last) = chain.last() {
            if last.is_empty() {
                break;
            }
            let next = self.bracket_spaces(last, last);
            chain.push(next);
        }
        for level in chain.windows(2) {
            let (upper, lower) = (&level[0], &level[1]);
            let mut lower_ech = lower.clone();
            let lower_piv = ql::rref(&mut lower_ech, n);
            let mut reduced: Vec<Vec<Q>> = upper.iter().map(|v| ql::reduce_mod(v, &lower_ech, &lower_piv)).collect();
            ql::rref(&mut reduced, n);
            let w = reduced;
            let d = w.len();
            if d == 0 {
                continue;
            }
            // right-hand sides: m_ij = [s_i, s_j] - sum_l c_ij^l s_l
            let mut rhs: Vec<Vec<Q>> = Vec::new();
            for i in 0..p {
                for j in i + 1..p {
                    let mut m = self.bracket(&s[i], &s[j]);
                    for (l, coefficient) in quotient[i][j].iter().enumerate() {
                        for (mk, sk) in m.iter_mut().zip(&s[l]) {
                            *mk -= coefficient * sk;
                        }
                    }
                    rhs.push(ql::reduce_mod(&m, &lower_ech, &lower_piv));
                }
            }
            // unknowns x[i][a]: u_i = sum_a x[i][a] w_a
            let variables = p * d;
            let mut columns: Vec<Vec<Q>> = Vec::with_capacity(variables);
            for i in 0..p {
                for wa in &w {
                    let mut col: Vec<Q> = Vec::new();
                    for i0 in 0..p {
                        for j0 in i0 + 1..p {
                            let mut t = ql::zeros(n);
                            if j0 == i {
                                for (tk, bk) in t.iter_mut().zip(self.bracket(&s[i0], wa)) {
                                    *tk += bk;
                                }
                            }
                            if i0 == i {
                                for (tk, bk) in t.iter_mut().zip(self.bracket(&s[j0], wa)) {
                                    *tk -= bk;
                                }
                            }
                            let coefficient = &quotient[i0][j0][i];
                            for (tk, wk) in t.iter_mut().zip(wa) {
                                *tk -= coefficient * wk;
                            }
                            col.extend(ql::reduce_mod(&t, &lower_ech, &lower_piv));
                        }
                    }
                    columns.push(col);
                }
            }
            let target: Vec<Q> = rhs.iter().flatten().map(|x| -x.clone()).collect();
            let system = ql::transpose(&columns, target.len());
            let x = ql::solve(&system, variables, &target)?;
            for (i, si) in s.iter_mut().enumerate() {
                for (a, wa) in w.iter().enumerate() {
                    for (sk, wk) in si.iter_mut().zip(wa) {
                        *sk += &x[i * d + a] * wk;
                    }
                }
            }
        }
        // verify: [s_i, s_j] = sum c_ij^l s_l exactly
        let mut echelon_s = s.clone();
        let piv_s = ql::rref(&mut echelon_s, n);
        let closed = s.iter().all(|x| s.iter().all(|y| ql::is_zero_vec(&ql::reduce_mod(&self.bracket(x, y), &echelon_s, &piv_s))));
        let radical = self.radical();
        closed.then_some((s, radical))
    }

    /// The generalised zero eigenspace of `ad x` (a Cartan subalgebra when
    /// `x` is regular).
    fn zero_eigenspace(
        &self,
        x: &[Q],
    ) -> QMat {
        let m = self.ad(x);
        let mut power = ql::identity(self.n);
        for _ in 0..self.n {
            power = ql::matmul(&m, &power);
        }
        ql::span_basis(&ql::nullspace(&power, self.n), self.n)
    }

    /// A Cartan subalgebra: the generalised zero eigenspace of a regular
    /// element, chosen among dense and sparse pseudo-random candidates for
    /// the smallest dimension and then for the most rational roots (sparse
    /// candidates find the split tori of matrix algebras with a diagonal
    /// part).
    #[must_use]
    pub fn cartan_subalgebra(&self) -> QMat {
        let mut state: u64 = 0x2545_f491_4f6c_dd1d;
        let mut next = move |bound: u64| {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % bound
        };
        let mut best: Option<(QMat, usize)> = None;
        for attempt in 0..40 {
            let mut x = ql::zeros(self.n);
            if attempt % 2 == 0 {
                for c in &mut x {
                    *c = ql::q(i64::try_from(next(11)).unwrap_or(0) - 5);
                }
            } else {
                for _ in 0..=next(3) {
                    let at = usize::try_from(next(self.n as u64)).unwrap_or(0);
                    x[at] = ql::q(i64::try_from(next(9)).unwrap_or(0) + 1);
                }
            }
            let h = self.zero_eigenspace(&x);
            if best.as_ref().is_some_and(|(b, _)| h.len() > b.len()) {
                continue;
            }
            let missing = self.root_decomposition(&h).map_or(self.n, |(_, m)| m);
            if best.as_ref().is_none_or(|(b, m)| h.len() < b.len() || missing < *m) {
                best = Some((h, missing));
            }
        }
        best.map(|(h, _)| h).unwrap_or_default()
    }

    /// The decomposition of `g` into simultaneous generalised eigenspaces
    /// of `ad h` for the abelian subalgebra spanned by `h`. Only roots with
    /// rational values are found; the second component counts the
    /// dimensions that belong to other (irrational or complex) roots.
    #[must_use]
    pub fn root_decomposition(
        &self,
        h: &[Vec<Q>],
    ) -> Option<(Vec<RootSpace>, usize)> {
        let n = self.n;
        if !h.iter().all(|x| h.iter().all(|y| ql::is_zero_vec(&self.bracket(x, y)))) {
            return None;
        }
        let mut parts: Vec<(Vec<Q>, QMat)> = vec![(Vec::new(), ql::identity(n))];
        for hi in h {
            let mut refined = Vec::new();
            for (root, space) in parts {
                let d = space.len();
                let system = ql::transpose(&space, n);
                let columns: Vec<Vec<Q>> = space.iter().map(|e| ql::solve(&system, d, &self.bracket(hi, e))).collect::<Option<_>>()?;
                let m = ql::transpose(&columns, d);
                for lambda in rational_roots(&charpoly(&m)).unwrap_or_default() {
                    let mut shifted = m.clone();
                    for (i, row) in shifted.iter_mut().enumerate() {
                        row[i] -= &lambda;
                    }
                    let mut power = ql::identity(d);
                    for _ in 0..d {
                        power = ql::matmul(&shifted, &power);
                    }
                    let kernel = ql::nullspace(&power, d);
                    let vectors: QMat = kernel
                        .iter()
                        .map(|y| {
                            let mut v = ql::zeros(n);
                            for (ya, ea) in y.iter().zip(&space) {
                                for (vk, ek) in v.iter_mut().zip(ea) {
                                    *vk += ya * ek;
                                }
                            }
                            v
                        })
                        .collect();
                    if !vectors.is_empty() {
                        let mut r = root.clone();
                        r.push(lambda);
                        refined.push((r, ql::span_basis(&vectors, n)));
                    }
                }
            }
            parts = refined;
        }
        let found: usize = parts.iter().map(|(_, s)| s.len()).sum();
        let spaces = parts.into_iter().map(|(root, space)| RootSpace { root, space }).collect();
        Some((spaces, n - found))
    }

    /// The quadratic Casimir operator `sum K^{ab} rho(e_a) rho(e_b)` of a
    /// representation `rho` (one matrix per basis element) for the Killing
    /// form; the algebra must be semisimple.
    #[must_use]
    pub fn casimir(
        &self,
        rho: &[QMat],
    ) -> Option<QMat> {
        if rho.len() != self.n {
            return None;
        }
        let inverse = ql::inverse(&self.killing())?;
        let size = rho.first()?.len();
        let mut total = vec![ql::zeros(size); size];
        for a in 0..self.n {
            for b in 0..self.n {
                if inverse[a][b].is_zero() {
                    continue;
                }
                let product = ql::matmul(&rho[a], &rho[b]);
                for (trow, prow) in total.iter_mut().zip(&product) {
                    for (t, p) in trow.iter_mut().zip(prow) {
                        *t += &inverse[a][b] * p;
                    }
                }
            }
        }
        Some(total)
    }
}

/// The characteristic polynomial `det(x I - m)` as coefficients of
/// `x^0 ..= x^n` (Faddeev-LeVerrier).
#[must_use]
pub fn charpoly(m: &[Vec<Q>]) -> Vec<Q> {
    let n = m.len();
    let mut coefficients = vec![Q::zero(); n + 1];
    coefficients[n] = Q::one();
    let mut aux = vec![ql::zeros(n); n];
    for k in 1..=n {
        let mut next = ql::matmul(m, &aux);
        for (i, row) in next.iter_mut().enumerate() {
            row[i] += &coefficients[n - k + 1];
        }
        let product = ql::matmul(m, &next);
        let trace = (0..n).fold(Q::zero(), |acc, i| acc + &product[i][i]);
        coefficients[n - k] = -trace / ql::q(i64::try_from(k).unwrap_or(1));
        aux = next;
    }
    coefficients
}

fn divisors(n: &BigInt) -> Option<Vec<BigInt>> {
    let v = n.abs().to_u64()?;
    if v == 0 || v > 1_000_000_000_000 {
        return None;
    }
    let mut out = Vec::new();
    for d in (1..).take_while(|d| d * d <= v) {
        if v % d == 0 {
            out.push(BigInt::from(d));
            if d * d != v {
                out.push(BigInt::from(v / d));
            }
        }
    }
    Some(out)
}

/// The distinct rational roots of the polynomial with coefficients
/// `p[0] ..= p[n]` (rational root theorem); `None` if the coefficients are
/// too large to factor.
#[must_use]
pub fn rational_roots(p: &[Q]) -> Option<Vec<Q>> {
    let denominator = p.iter().fold(BigInt::one(), |acc, c| num_integer::lcm(acc, c.denom().clone()));
    let ints: Vec<BigInt> = p.iter().map(|c| (c * Q::from_integer(denominator.clone())).to_integer()).collect();
    let low = ints.iter().position(|c| !c.is_zero())?;
    let mut roots = Vec::new();
    if low > 0 {
        roots.push(Q::zero());
    }
    let trimmed = &ints[low..];
    if trimmed.len() > 1 {
        let (a0, an) = (trimmed.first()?, trimmed.last()?);
        let (num_d, den_d) = (divisors(a0)?, divisors(an)?);
        let mut seen: HashSet<Q> = HashSet::new();
        for a in &num_d {
            for b in &den_d {
                for sign in [1, -1] {
                    let candidate = Q::new(a.clone() * BigInt::from(sign), b.clone());
                    if !seen.insert(candidate.clone()) {
                        continue;
                    }
                    let value = trimmed.iter().rev().fold(Q::zero(), |acc, c| acc * &candidate + Q::from_integer(c.clone()));
                    if value.is_zero() {
                        roots.push(candidate);
                    }
                }
            }
        }
    }
    Some(roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rs(name: &str) -> RootSystem {
        RootSystem::new(parse_type(name).expect("type")).expect("finite type")
    }

    fn dim(
        name: &str,
        labels: &[i64],
    ) -> i64 {
        rs(name).weyl_dim(labels).and_then(|d| d.to_i64()).expect("dimension")
    }

    #[test]
    fn root_system_data() {
        let table = [
            ("A2", 3, 6, 3),
            ("A3", 6, 24, 4),
            ("B2", 4, 8, 4),
            ("B3", 9, 48, 6),
            ("C3", 9, 48, 6),
            ("D4", 12, 192, 6),
            ("E6", 36, 51_840, 12),
            ("E7", 63, 2_903_040, 18),
            ("E8", 120, 696_729_600, 30),
            ("F4", 24, 1152, 12),
            ("G2", 6, 12, 6),
        ];
        for (name, positive, order, h) in table {
            let r = rs(name);
            assert_eq!(r.num_positive(), positive, "{name}");
            assert_eq!(r.weyl_order(), BigInt::from(order), "{name}");
            assert_eq!(r.coxeter_number(), Some(h), "{name}");
        }
        assert_eq!(rs("E8").dimension(), 248);
        assert_eq!(rs("G2").highest_root(), Some(vec![3, 2]));
        assert_eq!(rs("A3").highest_root(), Some(vec![1, 1, 1]));
        assert_eq!(rs("E8").dual_coxeter_number(), Some(ql::q(30)));
        assert_eq!(rs("G2").dual_coxeter_number(), Some(ql::q(4)));
        assert_eq!(rs("B3").dual_coxeter_number(), Some(ql::q(5)));
        assert_eq!(rs("C3").dual_coxeter_number(), Some(ql::q(4)));
        assert_eq!(rs("F4").dual_coxeter_number(), Some(ql::q(9)));
        assert_eq!(rs("A2").exponents(), vec![1, 2]);
        assert_eq!(rs("E8").exponents(), vec![1, 7, 11, 13, 17, 19, 23, 29]);
    }

    #[test]
    fn recognition() {
        for (letter, n) in [('A', 1), ('A', 5), ('B', 2), ('B', 4), ('C', 3), ('C', 5), ('D', 4), ('D', 6), ('E', 6), ('E', 7), ('E', 8), ('F', 4), ('G', 2)] {
            let c = cartan_matrix(letter, n).expect("cartan");
            let expected = if (letter, n) == ('C', 2) { ('B', 2) } else { (letter, n) };
            assert_eq!(recognize(&c), Some(vec![expected]));
        }
        // permuted and reducible
        let c = vec![vec![2, 0, -1], vec![0, 2, 0], vec![-1, 0, 2]];
        assert_eq!(recognize(&c), Some(vec![('A', 2), ('A', 1)]));
        let c = vec![vec![2, -1], vec![-2, 2]];
        assert_eq!(recognize(&c), Some(vec![('B', 2)]));
        let b3 = cartan_matrix('B', 3).expect("cartan");
        let reversed: Cartan = (0..3).rev().map(|i| (0..3).rev().map(|j| b3[i][j]).collect()).collect();
        assert_eq!(recognize(&reversed), Some(vec![('B', 3)]));
        // affine A1 and a non-Cartan matrix are rejected
        assert!(recognize(&vec![vec![2, -2], vec![-2, 2]]).is_none());
        assert!(recognize(&vec![vec![2, -1], vec![0, 2]]).is_none());
        // A4 relabelled
        let a4 = vec![vec![2, 0, 0, -1], vec![0, 2, -1, 0], vec![0, -1, 2, -1], vec![-1, 0, -1, 2]];
        assert_eq!(recognize(&a4), Some(vec![('A', 4)]));
    }

    #[test]
    fn euclidean_roots_reproduce_the_cartan_matrices() {
        for (letter, n) in [('A', 3), ('B', 3), ('C', 3), ('D', 5), ('E', 6), ('E', 7), ('E', 8), ('F', 4), ('G', 2)] {
            let roots = euclidean_simple_roots(letter, n).expect("roots");
            let cartan = cartan_matrix(letter, n).expect("cartan");
            let dot = |a: &[Q], b: &[Q]| a.iter().zip(b).fold(Q::zero(), |acc, (x, y)| acc + x * y);
            for i in 0..n {
                for j in 0..n {
                    let value = ql::q(2) * dot(&roots[i], &roots[j]) / dot(&roots[j], &roots[j]);
                    assert_eq!(value, ql::q(cartan[i][j]), "{letter}{n} ({i},{j})");
                }
            }
        }
    }

    #[test]
    fn weyl_dimensions() {
        assert_eq!(dim("A1", &[1]), 2);
        assert_eq!(dim("A2", &[1, 0]), 3);
        assert_eq!(dim("A2", &[0, 1]), 3);
        assert_eq!(dim("A2", &[2, 0]), 6);
        assert_eq!(dim("A2", &[1, 1]), 8);
        assert_eq!(dim("A2", &[3, 0]), 10);
        assert_eq!(dim("A2", &[2, 2]), 27);
        assert_eq!(dim("A3", &[0, 1, 0]), 6);
        assert_eq!(dim("B2", &[0, 1]), 4);
        assert_eq!(dim("B2", &[1, 0]), 5);
        assert_eq!(dim("G2", &[1, 0]), 7);
        assert_eq!(dim("G2", &[0, 1]), 14);
        assert_eq!(dim("F4", &[0, 0, 0, 1]), 26);
        assert_eq!(dim("E6", &[1, 0, 0, 0, 0, 0]), 27);
        assert_eq!(dim("E7", &[0, 0, 0, 0, 0, 0, 1]), 56);
        assert_eq!(dim("E8", &[0, 0, 0, 0, 0, 0, 0, 1]), 248);
        assert_eq!(dim("E8", &[1, 0, 0, 0, 0, 0, 0, 0]), 3875);
        assert_eq!(dim("D4", &[1, 0, 0, 0]), 8);
        assert_eq!(dim("D4", &[0, 1, 0, 0]), 28);
        // the adjoint representation has the dimension of the algebra
        for name in ["A4", "B3", "C4", "D5", "E6", "E7", "E8", "F4", "G2"] {
            let r = rs(name);
            let theta = r.highest_root().expect("irreducible");
            let labels: Vec<i64> = (0..r.n).map(|i| (0..r.n).map(|j| theta[j] * r.cartan[j][i]).sum()).collect();
            assert_eq!(r.weyl_dim(&labels), Some(BigInt::from(r.dimension())), "{name}");
        }
    }

    #[test]
    fn casimir_and_index() {
        let su3 = rs("A2");
        assert_eq!(su3.casimir(&[1, 0]), Some(Q::new(BigInt::from(8), BigInt::from(3))));
        assert_eq!(su3.casimir(&[1, 1]), Some(ql::q(6)));
        assert_eq!(su3.dynkin_index(&[1, 0]), Some(ql::q(1)));
        assert_eq!(su3.dynkin_index(&[1, 1]), Some(ql::q(6)));
        assert_eq!(su3.dynkin_index(&[2, 0]), Some(ql::q(5)));
        let e8 = rs("E8");
        assert_eq!(e8.casimir(&[0, 0, 0, 0, 0, 0, 0, 1]), Some(ql::q(60)));
        assert_eq!(rs("G2").dynkin_index(&[1, 0]), Some(ql::q(2)));
    }

    #[test]
    fn freudenthal() {
        let su3 = rs("A2");
        assert_eq!(su3.weight_multiplicity(&[1, 1], &[0, 0]), Some(BigInt::from(2)));
        assert_eq!(su3.weight_multiplicity(&[1, 1], &[1, 1]), Some(BigInt::one()));
        assert_eq!(su3.weight_multiplicity(&[1, 1], &[-1, 2]), Some(BigInt::one()));
        assert_eq!(su3.weight_multiplicity(&[2, 2], &[0, 0]), Some(BigInt::from(3)));
        assert_eq!(su3.weight_multiplicity(&[1, 0], &[0, 0]), Some(BigInt::zero()));
        for (name, top) in [("A2", vec![2, 1]), ("A3", vec![1, 1, 1]), ("B3", vec![1, 0, 1]), ("G2", vec![1, 1]), ("D4", vec![0, 1, 1, 0]), ("F4", vec![1, 0, 0, 0]), ("E6", vec![0, 0, 0, 0, 0, 1])] {
            let r = rs(name);
            let total: BigInt = r.all_weights(&top, 300_000).expect("weights").values().sum();
            assert_eq!(Some(total), r.weyl_dim(&top), "{name}");
        }
    }

    #[test]
    fn tensor_products() {
        let su3 = rs("A2");
        let as_pairs = |v: Vec<(Vec<i64>, BigInt)>| -> Vec<(Vec<i64>, i64)> { v.into_iter().map(|(w, m)| (w, m.to_i64().expect("small"))).collect() };
        assert_eq!(as_pairs(su3.tensor(&[1, 0], &[0, 1]).expect("tensor")), vec![(vec![0, 0], 1), (vec![1, 1], 1)]);
        assert_eq!(as_pairs(su3.tensor(&[1, 0], &[1, 0]).expect("tensor")), vec![(vec![0, 1], 1), (vec![2, 0], 1)]);
        assert_eq!(
            as_pairs(su3.tensor(&[1, 1], &[1, 1]).expect("tensor")),
            vec![(vec![0, 0], 1), (vec![0, 3], 1), (vec![1, 1], 2), (vec![2, 2], 1), (vec![3, 0], 1)]
        );
        for (name, a, b) in [("A2", vec![2, 1], vec![1, 1]), ("B2", vec![1, 1], vec![0, 1]), ("G2", vec![1, 0], vec![1, 0]), ("A3", vec![1, 0, 1], vec![0, 1, 0]), ("D4", vec![1, 0, 0, 0], vec![0, 0, 1, 0])] {
            let r = rs(name);
            let parts = r.tensor(&a, &b).expect("tensor");
            let total: BigInt = parts.iter().map(|(w, m)| r.weyl_dim(w).expect("dim") * m).sum();
            assert_eq!(total, r.weyl_dim(&a).expect("dim") * r.weyl_dim(&b).expect("dim"), "{name}");
        }
        // 7 x 7 of G2 = 1 + 7 + 14 + 27
        let g2 = rs("G2");
        let parts = g2.tensor(&[1, 0], &[1, 0]).expect("tensor");
        let dims: Vec<i64> = parts.iter().map(|(w, _)| g2.weyl_dim(w).and_then(|d| d.to_i64()).expect("dim")).collect();
        let mut dims = dims;
        dims.sort_unstable();
        assert_eq!(dims, vec![1, 7, 14, 27]);
    }

    fn matrices(basis: &[GMat]) -> Vec<QMat> {
        basis.iter().map(|m| m.iter().map(|r| r.iter().map(|x| ql::q(x.0)).collect()).collect()).collect()
    }

    fn real_algebra(
        family: Family,
        n: usize,
    ) -> Lie {
        Lie::from_matrices(&matrices(&standard_basis(family, n).expect("basis"))).expect("a Lie algebra")
    }

    #[test]
    fn standard_bases_are_lie_algebras() {
        for (family, n, dim) in [(Family::Sl, 3, 8), (Family::So, 4, 6), (Family::So, 5, 10), (Family::Sp, 2, 10), (Family::Sp, 3, 21), (Family::Gl, 3, 9)] {
            assert_eq!(real_algebra(family, n).n, dim);
        }
        assert_eq!(standard_basis(Family::Su, 3).expect("su3").len(), 8);
    }

    #[test]
    fn solvable_nilpotent_semisimple() {
        let upper: Vec<QMat> = {
            let mut v = Vec::new();
            for i in 0..3 {
                for j in i..3 {
                    let mut m = vec![ql::zeros(3); 3];
                    m[i][j] = ql::q(1);
                    v.push(m);
                }
            }
            v
        };
        let b = Lie::from_matrices(&upper).expect("closed");
        assert!(b.is_solvable() && !b.is_nilpotent() && !b.is_semisimple() && b.cartan_solvable_test());
        assert_eq!(b.derived_series().iter().map(Vec::len).collect::<Vec<_>>(), vec![6, 3, 1, 0]);
        assert_eq!(b.radical().len(), 6);
        let heisenberg: Vec<QMat> = [(0, 1), (1, 2), (0, 2)]
            .iter()
            .map(|&(i, j)| {
                let mut m = vec![ql::zeros(3); 3];
                m[i][j] = ql::q(1);
                m
            })
            .collect();
        let h = Lie::from_matrices(&heisenberg).expect("closed");
        assert!(h.is_nilpotent() && h.is_solvable());
        assert_eq!(h.lower_central_series().iter().map(Vec::len).collect::<Vec<_>>(), vec![3, 1, 0]);
        assert_eq!(h.upper_central_series().iter().map(Vec::len).collect::<Vec<_>>(), vec![1, 3]);
        assert_eq!(h.center().len(), 1);
        for (family, n) in [(Family::Sl, 2), (Family::Sl, 3), (Family::So, 4), (Family::Sp, 2)] {
            let g = real_algebra(family, n);
            assert!(g.is_semisimple() && !g.is_solvable() && !g.cartan_solvable_test() && g.radical().is_empty(), "{family:?}{n}");
        }
        assert!(!real_algebra(Family::So, 3).is_solvable());
        let gl2 = real_algebra(Family::Gl, 2);
        assert!(!gl2.is_semisimple() && !gl2.is_solvable());
        assert_eq!(gl2.radical().len(), 1);
        assert_eq!(gl2.center().len(), 1);
    }

    #[test]
    fn levi_decomposition() {
        // sl2 acting on its standard representation: g = sl2 + C^2 inside 3x3 matrices
        let e = |i: usize, j: usize| -> QMat {
            let mut m = vec![ql::zeros(3); 3];
            m[i][j] = ql::q(1);
            m
        };
        let mut h = e(0, 0);
        h[1][1] = ql::q(-1);
        let basis = vec![h, e(0, 1), e(1, 0), e(0, 2), e(1, 2)];
        let g = Lie::from_matrices(&basis).expect("closed");
        assert!(!g.is_semisimple() && !g.is_solvable());
        let (levi, radical) = g.levi().expect("levi decomposition");
        assert_eq!((levi.len(), radical.len()), (3, 2));
        assert!(g.is_subalgebra(&levi));
        // the Levi factor is semisimple: its Killing form is nondegenerate
        let structure: Vec<Vec<Vec<Q>>> = {
            let system = ql::transpose(&levi, 5);
            levi.iter().map(|x| levi.iter().map(|y| ql::solve(&system, 3, &g.bracket(x, y)).expect("closed")).collect()).collect()
        };
        assert!(Lie::from_constants(structure).expect("Lie algebra").is_semisimple());
        // a solvable algebra has zero Levi factor, a semisimple one is its own
        let (l, r) = real_algebra(Family::Sl, 2).levi().expect("levi");
        assert_eq!((l.len(), r.len()), (3, 0));
        let (l, r) = real_algebra(Family::Gl, 2).levi().expect("levi");
        assert_eq!((l.len(), r.len()), (3, 1));
    }

    #[test]
    fn cartan_subalgebras_and_roots() {
        let sl3 = real_algebra(Family::Sl, 3);
        let h = sl3.cartan_subalgebra();
        assert_eq!(h.len(), 2);
        let (spaces, missing) = sl3.root_decomposition(&h).expect("abelian");
        assert_eq!(missing, 0);
        assert_eq!(spaces.len(), 7);
        assert_eq!(spaces.iter().filter(|s| s.root.iter().all(Zero::is_zero)).map(|s| s.space.len()).sum::<usize>(), 2);
        assert!(spaces.iter().filter(|s| s.root.iter().any(|r| !r.is_zero())).all(|s| s.space.len() == 1));
        let sp4 = real_algebra(Family::Sp, 2);
        let h = sp4.cartan_subalgebra();
        assert_eq!(h.len(), 2);
        let (spaces, missing) = sp4.root_decomposition(&h).expect("abelian");
        assert_eq!((missing, spaces.len()), (0, 9));
        // so(3) is not split over Q: only the zero weight space is found
        let so3 = real_algebra(Family::So, 3);
        let h = so3.cartan_subalgebra();
        assert_eq!(h.len(), 1);
        let (spaces, missing) = so3.root_decomposition(&h).expect("abelian");
        assert_eq!((spaces.len(), missing), (1, 2));
        // a nilpotent algebra is its own Cartan subalgebra
        let upper_unit = Lie::from_constants(vec![
            vec![vec![ql::q(0), ql::q(0), ql::q(0)], vec![ql::q(0), ql::q(0), ql::q(1)], vec![ql::q(0), ql::q(0), ql::q(0)]],
            vec![vec![ql::q(0), ql::q(0), ql::q(-1)], vec![ql::q(0), ql::q(0), ql::q(0)], vec![ql::q(0), ql::q(0), ql::q(0)]],
            vec![vec![ql::q(0); 3]; 3],
        ]);
        assert_eq!(upper_unit.expect("heisenberg").cartan_subalgebra().len(), 3);
    }

    #[test]
    fn casimir_operators() {
        let sl2 = real_algebra(Family::Sl, 2);
        let rho = matrices(&standard_basis(Family::Sl, 2).expect("basis"));
        let c = sl2.casimir(&rho).expect("semisimple");
        // the defining representation of sl2: Casimir for the Killing form is (3/8) I
        assert_eq!(c[0][0], Q::new(BigInt::from(3), BigInt::from(8)));
        assert_eq!(c[0][1], Q::zero());
        assert_eq!(c[1][1], c[0][0]);
        // on the adjoint representation it is the identity
        let ad: Vec<QMat> = ql::identity(3).iter().map(|e| sl2.ad(e)).collect();
        let c = sl2.casimir(&ad).expect("semisimple");
        assert_eq!(c, ql::identity(3));
    }

    #[test]
    fn polynomial_roots() {
        let m: QMat = vec![vec![ql::q(1), ql::q(2)], vec![ql::q(3), ql::q(0)]];
        assert_eq!(charpoly(&m), vec![ql::q(-6), ql::q(-1), ql::q(1)]);
        let mut roots = rational_roots(&charpoly(&m)).expect("small");
        roots.sort();
        assert_eq!(roots, vec![ql::q(-2), ql::q(3)]);
        assert_eq!(rational_roots(&[ql::q(1), ql::q(0), ql::q(1)]), Some(Vec::new()));
    }
}
