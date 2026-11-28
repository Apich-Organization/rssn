//! Permutation groups: Schreier–Sims, Todd–Coxeter and derived structure.
//!
//! A permutation of `n` points is the list of the images of `1, ..., n`
//! (see [`groups`](super::groups)); a permutation group is given by a list
//! of generators of equal degree. The product is `(p q)(i) = p(q(i))`.
//!
//! The group is represented by a *stabiliser chain* (a base and a strong
//! generating set) computed by the deterministic Schreier–Sims algorithm, so
//! order, membership and orbit–stabiliser data are exact for groups of any
//! order that fits the degree (the order is an arbitrary-precision integer).
//!
//! | operator | value |
//! |---|---|
//! | `perm_group_order(gens)` | the order, by Schreier–Sims |
//! | `perm_group_contains(gens, p)` | membership test (sifting) |
//! | `perm_group_base(gens)` | a base: points whose pointwise stabiliser is trivial |
//! | `perm_group_strong_generators(gens)` | a strong generating set relative to the base |
//! | `perm_group_basic_orbits(gens)` | the lengths of the basic orbits (their product is the order) |
//! | `perm_group_orbits(gens)` | the orbits on the points |
//! | `perm_group_is_transitive(gens)` | one orbit |
//! | `perm_group_stabilizer_order(gens, i)` | order of the stabiliser of point `i` |
//! | `perm_group_elements(gens)` | all elements (order at most 5040) |
//! | `group_from_perms(gens)` | the `group(elements, table)` term (order at most 720) |
//! | `perm_group_is_abelian(gens)` | commutativity |
//! | `perm_group_derived_series(gens)` | orders of the derived series `G = G0 > G1 > ...` (stops at the perfect core) |
//! | `perm_group_is_solvable(gens)` | the derived series reaches the trivial group |
//! | `perm_group_lower_central_series(gens)` | orders of `G = g1 > g2 > ...` with `g(i+1) = [gi, G]` |
//! | `perm_group_is_nilpotent(gens)` | the lower central series reaches the trivial group |
//! | `perm_group_derived_subgroup(gens)` | generators of `G'` |
//! | `todd_coxeter(gens, relators, subgroup_gens)` | `list(index, list(perm per generator))`: the coset enumeration of `<gens \| relators>` over the subgroup, by the HLT strategy with coincidences |
//! | `todd_coxeter_index(gens, relators, subgroup_gens)` | only the index |
//! | `fp_group_order(gens, relators)` | the order of a presented group (index of the trivial subgroup) |
//! | `fp_group(gens, relators)` | the presented group as a `group` term (order at most 720) |
//! | `fp_group_abelian_invariants(gens, relators)` | the invariant factors `> 1` of the abelianisation, then `0` per free factor (Smith normal form of the exponent sums) |
//!
//! Presentations: `gens` is a list of generator symbols; a word is a list of
//! letters, each a generator `a`, a power `a^k` (negative for the inverse),
//! or a non-zero integer `+-i` for the `i`-th generator or its inverse. The
//! relators are words; `list(a, a, a)`, `list(a^3)` and `list(1, 1, 1)` all
//! stand for `a^3`.
//! Index-only requests verify live transitions without constructing action
//! arrays. `fp_group` checks freely reduced relators and the generated
//! permutation group's order before building a Cayley table.

// Index loops over several parallel tables and float comparisons with explicit
// tolerances read better than the iterator forms; `as` casts convert small,
// non-negative, bounded values.
#![allow(
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::cast_sign_loss,
    clippy::manual_midpoint
)]

use std::borrow::Cow;
use std::collections::HashMap;

use num_bigint::BigInt;

use super::V;
use super::def_doc;
use super::groups::build_group;
use super::idx;
use super::items;
use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::RuleError;
use crate::graph::op::core;
use crate::graph::rule::Installer;

/// A permutation as 0-based images.
pub(crate) type Perm = Vec<usize>;

/// The largest order for which all elements are listed.
const ELEMENT_LIMIT: usize = 5040;
/// The largest order for which a table is built.
pub(crate) const TABLE_LIMIT: usize = 720;
/// The largest number of cosets of a Todd–Coxeter enumeration.
const COSET_LIMIT: usize = 400_000;

pub(crate) fn identity(n: usize) -> Perm {
    (0..n).collect()
}

/// `a ∘ b`: apply `b` first.
pub(crate) fn compose(
    a: &[usize],
    b: &[usize],
) -> Perm {
    b.iter().map(|&x| a[x]).collect()
}

pub(crate) fn inverse(p: &[usize]) -> Perm {
    let mut inv = vec![0; p.len()];
    for (i, &x) in p.iter().enumerate() {
        inv[x] = i;
    }
    inv
}

#[cfg(test)]
pub(crate) fn is_identity(p: &[usize]) -> bool {
    p.iter().enumerate().all(|(i, &x)| i == x)
}

/// `a^-1 b^-1 a b`.
fn commutator(
    a: &[usize],
    b: &[usize],
) -> Perm {
    compose(&compose(&inverse(a), &inverse(b)), &compose(a, b))
}

/// A stabiliser chain: the kernel's [`Group`](crate::kernels::permgroup::Group)
/// (incremental deterministic Schreier–Sims) behind this module's `usize`
/// permutations.
#[derive(Clone, Debug)]
pub(crate) struct Chain {
    group: crate::kernels::permgroup::Group,
    base: Vec<usize>,
}

fn to_kernel(p: &[usize]) -> crate::kernels::permgroup::Perm {
    p.iter()
        .map(|&x| u32::try_from(x).unwrap_or(u32::MAX))
        .collect()
}

fn from_kernel(p: &[u32]) -> Perm {
    p.iter().map(|&x| x as usize).collect()
}

impl Chain {
    /// The chain of the group generated by `gens` (all of degree `n`).
    pub(crate) fn new(
        n: usize,
        gens: &[Perm],
    ) -> Self {
        let gens: Vec<_> = gens.iter().map(|g| to_kernel(g)).collect();
        let group = crate::kernels::permgroup::Group::new(n, &gens);
        let base = group.base().iter().map(|&b| b as usize).collect();
        Self { group, base }
    }

    /// Adds `g` to the group.
    pub(crate) fn insert(
        &mut self,
        g: &Perm,
    ) {
        self.group.extend(to_kernel(g));
        self.base = self.group.base().iter().map(|&b| b as usize).collect();
    }

    /// Whether `g` lies in the group.
    pub(crate) fn contains(
        &self,
        g: &[usize],
    ) -> bool {
        self.group.contains(&to_kernel(g))
    }

    /// The order.
    pub(crate) fn order(&self) -> BigInt {
        BigInt::from(self.group.order())
    }

    /// The base.
    pub(crate) fn base(&self) -> &[usize] {
        &self.base
    }

    /// The lengths of the basic orbits.
    pub(crate) fn basic_orbit_lengths(&self) -> Vec<usize> {
        self.group.basic_orbit_lengths()
    }

    /// A strong generating set.
    pub(crate) fn strong_generators(&self) -> Vec<Perm> {
        self.group
            .strong_generators()
            .iter()
            .map(|p| from_kernel(p))
            .collect()
    }

    /// The order of the stabiliser of `point`: `|G| / |point^G|`.
    pub(crate) fn stabilizer_order(
        &self,
        point: usize,
    ) -> BigInt {
        let orbit = self
            .orbits()
            .into_iter()
            .find(|o| o.contains(&point))
            .map_or(1, |o| o.len());
        self.order() / BigInt::from(orbit)
    }

    /// All elements, sorted, if there are at most `limit`.
    pub(crate) fn elements(
        &self,
        limit: usize,
    ) -> Option<Vec<Perm>> {
        let mut all: Vec<Perm> = self
            .group
            .elements(limit)?
            .iter()
            .map(|p| from_kernel(p))
            .collect();
        all.sort();
        Some(all)
    }

    /// The orbits on the points, each sorted.
    pub(crate) fn orbits(&self) -> Vec<Vec<usize>> {
        let gens = self.strong_generators();
        orbits_of(self.group.degree(), &gens)
    }
}

pub(crate) fn orbits_of(
    n: usize,
    gens: &[Perm],
) -> Vec<Vec<usize>> {
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut orbit = vec![start];
        seen[start] = true;
        let mut head = 0;
        while head < orbit.len() {
            let pt = orbit[head];
            head += 1;
            for g in gens {
                if !seen[g[pt]] {
                    seen[g[pt]] = true;
                    orbit.push(g[pt]);
                }
            }
        }
        orbit.sort_unstable();
        out.push(orbit);
    }
    out
}

/// The normal closure in `<group_gens>` of `subset`.
fn normal_closure(
    n: usize,
    group_gens: &[Perm],
    subset: &[Perm],
) -> Chain {
    let mut chain = Chain::new(n, subset);
    loop {
        let mut changed = false;
        for h in chain.strong_generators() {
            for g in group_gens {
                let c = compose(&compose(g, &h), &inverse(g));
                if !chain.contains(&c) {
                    chain.insert(&c);
                    changed = true;
                }
            }
        }
        if !changed {
            return chain;
        }
    }
}

/// The derived subgroup of `<gens>`.
pub(crate) fn derived_subgroup(
    n: usize,
    gens: &[Perm],
) -> Chain {
    let mut commutators = Vec::new();
    for (i, a) in gens.iter().enumerate() {
        for b in &gens[i + 1..] {
            commutators.push(commutator(a, b));
        }
    }
    normal_closure(n, gens, &commutators)
}

/// `[<sub>, <gens>]` for a normal subgroup `<sub>`.
fn commutator_with(
    n: usize,
    gens: &[Perm],
    sub: &[Perm],
) -> Chain {
    let mut commutators = Vec::new();
    for a in sub {
        for g in gens {
            commutators.push(commutator(a, g));
        }
    }
    normal_closure(n, gens, &commutators)
}

/// The orders of the derived series, ending at the first repeated term.
pub(crate) fn derived_series(
    n: usize,
    gens: &[Perm],
) -> Vec<BigInt> {
    let mut chain = Chain::new(n, gens);
    let mut out = vec![chain.order()];
    loop {
        let g = chain.strong_generators();
        let next = derived_subgroup(n, &g);
        let o = next.order();
        let stop = o == *out.last().unwrap_or(&o);
        if !stop {
            out.push(o);
        }
        if stop || next.order() == BigInt::from(1) {
            return out;
        }
        chain = next;
    }
}

/// The orders of the lower central series, ending at the first repeat.
pub(crate) fn lower_central_series(
    n: usize,
    gens: &[Perm],
) -> Vec<BigInt> {
    let all = Chain::new(n, gens).strong_generators();
    let mut cur = all.clone();
    let mut out = vec![Chain::new(n, &all).order()];
    loop {
        let next = commutator_with(n, &all, &cur);
        let o = next.order();
        if o == *out.last().unwrap_or(&o) {
            return out;
        }
        let trivial = o == BigInt::from(1);
        out.push(o);
        if trivial {
            return out;
        }
        cur = next.strong_generators();
    }
}

// ----------------------------------------------------------------------
// Todd–Coxeter
// ----------------------------------------------------------------------

const NONE: usize = usize::MAX;

/// A coset table under construction (HLT strategy).
struct Cosets {
    ncols: usize,
    table: Vec<usize>,
    parent: Vec<usize>,
    pending: Vec<usize>,
    live: usize,
}

impl Cosets {
    fn new(ncols: usize) -> Self {
        let mut c = Self {
            ncols,
            table: Vec::new(),
            parent: Vec::new(),
            pending: Vec::new(),
            live: 0,
        };
        c.new_coset();
        c
    }

    fn new_coset(&mut self) -> usize {
        let id = self.parent.len();
        self.parent.push(id);
        self.table.extend(std::iter::repeat_n(NONE, self.ncols));
        self.live += 1;
        id
    }

    fn rep(
        &mut self,
        mut k: usize,
    ) -> usize {
        let mut root = k;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        while self.parent[k] != root {
            let next = self.parent[k];
            self.parent[k] = root;
            k = next;
        }
        root
    }

    fn get(
        &self,
        c: usize,
        x: usize,
    ) -> usize {
        self.table[c * self.ncols + x]
    }

    fn set(
        &mut self,
        c: usize,
        x: usize,
        v: usize,
    ) {
        self.table[c * self.ncols + x] = v;
    }

    fn merge(
        &mut self,
        k: usize,
        l: usize,
        queue: &mut Vec<usize>,
    ) {
        let (k, l) = (self.rep(k), self.rep(l));
        if k != l {
            let (mu, nu) = (k.min(l), k.max(l));
            self.parent[nu] = mu;
            self.live -= 1;
            queue.push(nu);
        }
    }

    fn coincidence(
        &mut self,
        a: usize,
        b: usize,
    ) {
        let mut queue = std::mem::take(&mut self.pending);
        self.merge(a, b, &mut queue);
        let mut i = 0;
        while i < queue.len() {
            let gamma = queue[i];
            i += 1;
            for x in 0..self.ncols {
                let delta = self.get(gamma, x);
                if delta == NONE {
                    continue;
                }
                self.set(delta, x ^ 1, NONE);
                let mu = self.rep(gamma);
                let nu = self.rep(delta);
                let t = self.get(mu, x);
                if t == NONE {
                    let t2 = self.get(nu, x ^ 1);
                    if t2 == NONE {
                        self.set(mu, x, nu);
                        self.set(nu, x ^ 1, mu);
                    } else {
                        self.merge(mu, t2, &mut queue);
                    }
                } else {
                    self.merge(nu, t, &mut queue);
                }
            }
        }
        queue.clear();
        self.pending = queue;
    }

    /// Scans the word `w` (columns) from coset `c`, defining cosets as
    /// needed.
    fn scan_and_fill(
        &mut self,
        c: usize,
        w: &[usize],
    ) -> bool {
        let (mut f, mut b) = (c, c);
        let (mut i, mut j) = (0_usize, w.len());
        // `j` is one past the last unscanned letter.
        loop {
            while i < j {
                let next = self.get(f, w[i]);
                if next == NONE {
                    break;
                }
                f = next;
                i += 1;
            }
            if i == j {
                if f != b {
                    self.coincidence(f, b);
                }
                return true;
            }
            while j > i {
                let next = self.get(b, w[j - 1] ^ 1);
                if next == NONE {
                    break;
                }
                b = next;
                j -= 1;
            }
            if j == i {
                self.coincidence(f, b);
                return true;
            }
            if j == i + 1 {
                self.set(f, w[i], b);
                self.set(b, w[i] ^ 1, f);
                return true;
            }
            if self.parent.len() >= COSET_LIMIT {
                return false;
            }
            let n = self.new_coset();
            self.set(f, w[i], n);
            self.set(n, w[i] ^ 1, f);
        }
    }

    fn live_rows(&self) -> impl Iterator<Item = usize> + '_ {
        self.parent
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(i, root)| (i == root).then_some(i))
    }

    /// Flatten representatives and check that every live transition has its
    /// inverse before accepting the index or constructing an action.
    fn seal(&mut self) -> Option<()> {
        for row in 0..self.parent.len() {
            self.rep(row);
        }
        for row in self.live_rows() {
            for column in 0..self.ncols {
                let next = *self.parent.get(self.get(row, column))?;
                let previous = self.parent.get(self.get(next, column ^ 1))?;
                if *previous != row {
                    return None;
                }
            }
        }
        Some(())
    }

    fn permutations(&self) -> Vec<Perm> {
        let mut labels = vec![NONE; self.parent.len()];
        for (index, row) in self.live_rows().enumerate() {
            labels[row] = index;
        }
        let mut perms: Vec<_> = (0..self.ncols / 2)
            .map(|_| Vec::with_capacity(self.live))
            .collect();
        for row in self.live_rows() {
            for (g, perm) in perms.iter_mut().enumerate() {
                perm.push(labels[self.parent[self.get(row, 2 * g)]]);
            }
        }
        perms
    }

    fn word_image(
        &self,
        start: usize,
        word: &[usize],
    ) -> usize {
        word.iter()
            .fold(start, |row, &column| self.parent[self.get(row, column)])
    }

    fn verifies(
        &self,
        relators: &[Vec<usize>],
        subgens: &[Vec<usize>],
    ) -> bool {
        let relators: Vec<_> = relators.iter().map(|word| reduced_word(word)).collect();
        self.live_rows().all(|row| {
            relators
                .iter()
                .all(|word| self.word_image(row, word) == row)
        }) && subgens
            .iter()
            .all(|word| self.word_image(0, &reduced_word(word)) == 0)
    }
}

fn reduced_word(word: &[usize]) -> Cow<'_, [usize]> {
    if word.windows(2).all(|pair| pair[0] != pair[1] ^ 1) {
        return Cow::Borrowed(word);
    }
    let mut reduced = Vec::with_capacity(word.len());
    for &letter in word {
        if reduced.last() == Some(&(letter ^ 1)) {
            reduced.pop();
        } else {
            reduced.push(letter);
        }
    }
    Cow::Owned(reduced)
}

fn enumerate_cosets(
    ngens: usize,
    relators: &[Vec<usize>],
    subgens: &[Vec<usize>],
) -> Option<Cosets> {
    let ncols = ngens.checked_mul(2)?;
    if relators
        .iter()
        .chain(subgens)
        .flatten()
        .any(|&letter| letter >= ncols)
    {
        return None;
    }
    let mut t = Cosets::new(ncols);
    for w in subgens {
        if !t.scan_and_fill(0, w) {
            return None;
        }
    }
    let mut c = 0;
    while c < t.parent.len() {
        if t.parent[c] == c {
            for w in relators {
                if t.parent[c] != c {
                    break;
                }
                if !t.scan_and_fill(c, w) {
                    return None;
                }
            }
            if t.parent[c] == c {
                for x in 0..t.ncols {
                    if t.parent[c] == c && t.get(c, x) == NONE {
                        if t.parent.len() >= COSET_LIMIT {
                            return None;
                        }
                        let n = t.new_coset();
                        t.set(c, x, n);
                        t.set(n, x ^ 1, c);
                    }
                }
            }
        }
        c += 1;
    }
    t.seal()?;
    Some(t)
}

/// The coset table of `<gens | relators>` over the subgroup `subgens`:
/// the number of cosets and the generators' 0-based actions.
pub(crate) fn todd_coxeter_table(
    ngens: usize,
    relators: &[Vec<usize>],
    subgens: &[Vec<usize>],
) -> Option<(usize, Vec<Perm>)> {
    let table = enumerate_cosets(ngens, relators, subgens)?;
    Some((table.live, table.permutations()))
}

/// Reads a word over the generators `gens` into columns (`2i` for the
/// generator `i`, `2i + 1` for its inverse).
fn read_word(
    g: &Graph,
    gens: &[NodeId],
    word: NodeId,
) -> Option<Vec<usize>> {
    let letters = items(g, word)?;
    let mut out = Vec::new();
    for l in letters {
        read_letter(g, gens, l, &mut out)?;
    }
    Some(out)
}

fn read_letter(
    g: &Graph,
    gens: &[NodeId],
    l: NodeId,
    out: &mut Vec<usize>,
) -> Option<()> {
    let pos = |n: NodeId| gens.iter().position(|&s| g.same(s, n));
    if let Some(i) = pos(l) {
        out.push(2 * i);
        return Some(());
    }
    if let Some(v) = super::small(g, l) {
        let i = usize::try_from(v.unsigned_abs())
            .ok()?
            .checked_sub(1)
            .filter(|&i| i < gens.len())?;
        out.push(2 * i + usize::from(v < 0));
        return Some(());
    }
    if g.op(l) == core::POW {
        let [b, e] = g.children(l) else { return None };
        let i = pos(*b)?;
        let k = super::small(g, *e)?;
        if k.unsigned_abs() > 10_000 {
            return None;
        }
        for _ in 0..k.unsigned_abs() {
            out.push(2 * i + usize::from(k < 0));
        }
        return Some(());
    }
    None
}

/// The number of generators, the relators and the subgroup generators.
type Presentation = (usize, Vec<Vec<usize>>, Vec<Vec<usize>>);

/// Reads a presentation: generators, relators and subgroup generators.
fn read_presentation(
    g: &Graph,
    a: &[NodeId],
) -> Option<Presentation> {
    let gens = items(g, *a.first()?)?;
    if gens.is_empty() {
        return None;
    }
    let rels: Vec<Vec<usize>> = items(g, *a.get(1)?)?
        .into_iter()
        .map(|w| read_word(g, &gens, w))
        .collect::<Option<_>>()?;
    let subs: Vec<Vec<usize>> = match a.get(2) {
        | Some(&s) => {
            items(g, s)?
                .into_iter()
                .map(|w| read_word(g, &gens, w))
                .collect::<Option<_>>()?
        },
        | None => Vec::new(),
    };
    Some((gens.len(), rels, subs))
}

fn perm_v(p: &[usize]) -> V {
    V::ints(p.iter().map(|&x| i64::try_from(x + 1).unwrap_or(0)))
}

fn todd_coxeter(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, rels, subs) = read_presentation(cx.graph, a)?;
    let (index, perms) = todd_coxeter_table(n, &rels, &subs)?;
    Some(V::List(vec![
        V::uint(index),
        V::List(perms.iter().map(|p| perm_v(p)).collect()),
    ]))
}

fn todd_coxeter_index(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, rels, subs) = read_presentation(cx.graph, a)?;
    Some(V::uint(enumerate_cosets(n, &rels, &subs)?.live))
}

fn fp_group_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, rels, _) = read_presentation(cx.graph, a)?;
    Some(V::uint(enumerate_cosets(n, &rels, &[])?.live))
}

fn fp_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, rels, _) = read_presentation(cx.graph, a)?;
    let table = enumerate_cosets(n, &rels, &[])?;
    if table.live > TABLE_LIMIT || !table.verifies(&rels, &[]) {
        return None;
    }
    let chain = Chain::new(table.live, &table.permutations());
    if chain.order() != BigInt::from(table.live) {
        return None;
    }
    group_term_from_chain(cx.graph, &chain)
}

// ----------------------------------------------------------------------
// Terms
// ----------------------------------------------------------------------

/// The generators of a permutation group term: degree and permutations.
fn read_gens(
    g: &Graph,
    n: NodeId,
) -> Option<(usize, Vec<Perm>)> {
    let list = items(g, n)?;
    let mut perms = Vec::new();
    for p in list {
        perms.push(read_perm(g, p)?);
    }
    let degree = perms.first()?.len();
    if perms.iter().any(|p| p.len() != degree) || degree == 0 || degree > 1024 {
        return None;
    }
    Some((degree, perms))
}

/// A permutation as 0-based images.
pub(crate) fn read_perm(
    g: &Graph,
    n: NodeId,
) -> Option<Perm> {
    let images = super::idxs(g, n)?;
    let m = images.len();
    let mut seen = vec![false; m];
    let mut out = Vec::with_capacity(m);
    for v in images {
        let v = v.checked_sub(1).filter(|&v| v < m)?;
        if std::mem::replace(&mut seen[v], true) {
            return None;
        }
        out.push(v);
    }
    Some(out)
}

fn perm_node(
    g: &mut Graph,
    p: &[usize],
) -> NodeId {
    let nodes: Vec<NodeId> = p
        .iter()
        .map(|&x| g.int(i64::try_from(x + 1).unwrap_or(0)))
        .collect();
    super::list(g, &nodes)
}

/// The group term of the permutation group generated by `gens` of degree
/// `n`, if its order is at most [`TABLE_LIMIT`].
pub(crate) fn perm_group_term(
    g: &mut Graph,
    n: usize,
    gens: &[Perm],
) -> Option<V> {
    let chain = Chain::new(n, gens);
    group_term_from_chain(g, &chain)
}

fn group_term_from_chain(
    g: &mut Graph,
    chain: &Chain,
) -> Option<V> {
    let elems = chain.elements(TABLE_LIMIT)?;
    let table = cayley_rows(&elems)?;
    let nodes: Vec<NodeId> = elems.iter().map(|p| perm_node(g, p)).collect();
    build_group(g, &nodes, &table)
}

/// The product table in the supplied element order. Each row reuses its
/// permutation buffer, and independent rows can use the active CPU pool.
fn cayley_rows(elems: &[Perm]) -> Option<Vec<Vec<usize>>> {
    let degree = elems.first()?.len();
    let positions: HashMap<&[usize], usize> = elems
        .iter()
        .enumerate()
        .map(|(i, p)| (p.as_slice(), i))
        .collect();
    let row = |p: &Perm, scratch: &mut [usize]| -> Option<Vec<usize>> {
        let mut row = Vec::with_capacity(elems.len());
        for q in elems {
            for (slot, &point) in scratch.iter_mut().zip(q) {
                *slot = p[point];
            }
            row.push(*positions.get(&*scratch)?);
        }
        Some(row)
    };
    let work = elems
        .len()
        .saturating_mul(elems.len())
        .saturating_mul(degree);
    if work >= 32_768 && crate::execution::should_split(elems.len(), 8) {
        crate::kernels::par::map_slice(elems, 8, |p| row(p, &mut vec![0; degree]))
            .into_iter()
            .collect()
    } else {
        let mut scratch = vec![0; degree];
        elems.iter().map(|p| row(p, &mut scratch)).collect()
    }
}

fn chain_of(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<Chain> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(Chain::new(n, &gens))
}

fn perm_group_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Int(chain_of(cx, a)?.order()))
}

fn perm_group_contains(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [gens, p] = a else { return None };
    let (n, gens) = read_gens(cx.graph, *gens)?;
    let p = read_perm(cx.graph, *p)?;
    if p.len() != n {
        return Some(V::Bool(false));
    }
    Some(V::Bool(Chain::new(n, &gens).contains(&p)))
}

fn perm_group_base(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = chain_of(cx, a)?;
    Some(V::ints(
        c.base().iter().map(|&b| i64::try_from(b + 1).unwrap_or(0)),
    ))
}

fn perm_group_strong_generators(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = chain_of(cx, a)?;
    Some(V::List(
        c.strong_generators().iter().map(|p| perm_v(p)).collect(),
    ))
}

fn perm_group_basic_orbits(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = chain_of(cx, a)?;
    Some(V::ints(
        c.basic_orbit_lengths()
            .into_iter()
            .map(|l| i64::try_from(l).unwrap_or(0)),
    ))
}

fn perm_group_orbits(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = chain_of(cx, a)?;
    Some(V::List(
        c.orbits()
            .iter()
            .map(|o| V::ints(o.iter().map(|&x| i64::try_from(x + 1).unwrap_or(0))))
            .collect(),
    ))
}

fn perm_group_is_transitive(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(chain_of(cx, a)?.orbits().len() == 1))
}

fn perm_group_stabilizer_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [gens, pt] = a else { return None };
    let (n, gens) = read_gens(cx.graph, *gens)?;
    let pt = idx(cx.graph, *pt)?.checked_sub(1).filter(|&p| p < n)?;
    Some(V::Int(Chain::new(n, &gens).stabilizer_order(pt)))
}

fn perm_group_elements(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let c = chain_of(cx, a)?;
    Some(V::List(
        c.elements(ELEMENT_LIMIT)?
            .iter()
            .map(|p| perm_v(p))
            .collect(),
    ))
}

fn group_from_perms(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    perm_group_term(cx.graph, n, &gens)
}

fn perm_group_is_abelian(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (_, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(V::Bool(gens.iter().all(|x| {
        gens.iter().all(|y| compose(x, y) == compose(y, x))
    })))
}

fn orders_v(orders: Vec<BigInt>) -> V {
    V::List(orders.into_iter().map(V::Int).collect())
}

fn perm_group_derived_series(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(orders_v(derived_series(n, &gens)))
}

fn perm_group_is_solvable(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(V::Bool(
        derived_series(n, &gens)
            .last()
            .is_some_and(|o| *o == BigInt::from(1)),
    ))
}

fn perm_group_lower_central_series(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(orders_v(lower_central_series(n, &gens)))
}

fn perm_group_is_nilpotent(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(V::Bool(
        lower_central_series(n, &gens)
            .last()
            .is_some_and(|o| *o == BigInt::from(1)),
    ))
}

fn perm_group_derived_subgroup(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    let d = derived_subgroup(n, &gens);
    let g = d.strong_generators();
    if g.is_empty() {
        return Some(V::List(vec![perm_v(&identity(n))]));
    }
    Some(V::List(g.iter().map(|p| perm_v(p)).collect()))
}

/// The abelian invariants of a presented group: the invariant factors
/// `d > 1` of its abelianisation (`d_1 | d_2 | ...`), then a `0` for every
/// infinite cyclic factor. The abelianisation is `Z^n` modulo the exponent
/// sums of the relators, so the invariants are the Smith normal form of that
/// integer matrix.
fn fp_group_abelian_invariants(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, rels, _) = read_presentation(cx.graph, a)?;
    let matrix: Vec<Vec<i64>> = rels
        .iter()
        .map(|word| {
            let mut row = vec![0_i64; n];
            for &letter in word {
                if let Some(entry) = row.get_mut(letter / 2) {
                    *entry += if letter % 2 == 0 { 1 } else { -1 };
                }
            }
            row
        })
        .collect();
    let factors = if matrix.is_empty() {
        Vec::new()
    } else {
        crate::kernels::homology::smith_invariants(&matrix)
    };
    let rank = factors.len();
    let mut out: Vec<V> = factors
        .into_iter()
        .map(|d| num_traits::Signed::abs(&d))
        .filter(|d| !num_traits::One::is_one(d))
        .map(V::Int)
        .collect();
    out.extend(std::iter::repeat_with(|| V::int(0)).take(n.saturating_sub(rank)));
    Some(V::List(out))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def_doc(
        i,
        "perm_group_order",
        Arity::Fixed(1),
        perm_group_order,
        "Computes the exact order of the permutation group generated by `gens` using the Schreier-Sims algorithm.",
    )?;
    def_doc(
        i,
        "perm_group_contains",
        Arity::Fixed(2),
        perm_group_contains,
        "Tests whether permutation `p` belongs to the permutation group generated by `gens` by sifting through the stabilizer chain.",
    )?;
    def_doc(
        i,
        "perm_group_base",
        Arity::Fixed(1),
        perm_group_base,
        "Returns a Schreier-Sims base (points whose pointwise stabilizer is trivial) for the group generated by `gens`.",
    )?;
    def_doc(
        i,
        "perm_group_strong_generators",
        Arity::Fixed(1),
        perm_group_strong_generators,
        "Returns a strong generating set relative to the Schreier-Sims base for the group generated by `gens`.",
    )?;
    def_doc(
        i,
        "perm_group_basic_orbits",
        Arity::Fixed(1),
        perm_group_basic_orbits,
        "Returns the lengths of the basic orbits along the stabilizer chain (whose product equals the group order).",
    )?;
    def_doc(
        i,
        "perm_group_orbits",
        Arity::Fixed(1),
        perm_group_orbits,
        "Returns the orbits of the permutation group generated by `gens` on `1..n`.",
    )?;
    def_doc(
        i,
        "perm_group_is_transitive",
        Arity::Fixed(1),
        perm_group_is_transitive,
        "Returns `true` if the permutation group generated by `gens` acts transitively on `1..n`.",
    )?;
    def_doc(
        i,
        "perm_group_stabilizer_order",
        Arity::Fixed(2),
        perm_group_stabilizer_order,
        "Computes the order of the point stabilizer of `i` in the permutation group generated by `gens`.",
    )?;
    def_doc(
        i,
        "perm_group_elements",
        Arity::Fixed(1),
        perm_group_elements,
        "Enumerates all elements of the permutation group generated by `gens` (for group orders up to 5040).",
    )?;
    def_doc(
        i,
        "group_from_perms",
        Arity::Fixed(1),
        group_from_perms,
        "Builds a `group(elements, table)` term from permutation generators `gens` (for group orders up to 720).",
    )?;
    def_doc(
        i,
        "perm_group_is_abelian",
        Arity::Fixed(1),
        perm_group_is_abelian,
        "Returns `true` if the permutation group generated by `gens` is abelian.",
    )?;
    def_doc(
        i,
        "perm_group_derived_series",
        Arity::Fixed(1),
        perm_group_derived_series,
        "Computes the orders of the derived series `G = G_0 > G_1 > ...` of the group generated by `gens`, stopping at the perfect core.",
    )?;
    def_doc(
        i,
        "perm_group_is_solvable",
        Arity::Fixed(1),
        perm_group_is_solvable,
        "Returns `true` if the permutation group generated by `gens` is solvable.",
    )?;
    def_doc(
        i,
        "perm_group_lower_central_series",
        Arity::Fixed(1),
        perm_group_lower_central_series,
        "Computes the orders of the lower central series `G = gamma_1 > gamma_2 > ...` of the group generated by `gens`.",
    )?;
    def_doc(
        i,
        "perm_group_is_nilpotent",
        Arity::Fixed(1),
        perm_group_is_nilpotent,
        "Returns `true` if the permutation group generated by `gens` is nilpotent.",
    )?;
    def_doc(
        i,
        "perm_group_derived_subgroup",
        Arity::Fixed(1),
        perm_group_derived_subgroup,
        "Returns generators of the commutator (derived) subgroup `G'` of the group generated by `gens`.",
    )?;
    def_doc(
        i,
        "todd_coxeter",
        Arity::Fixed(3),
        todd_coxeter,
        "Runs Todd-Coxeter coset enumeration on `<gens | relators>` over `subgroup_gens`, returning `list(index, list(perm_per_gen))`.",
    )?;
    def_doc(
        i,
        "todd_coxeter_index",
        Arity::Fixed(3),
        todd_coxeter_index,
        "Computes the index of `subgroup_gens` in the finitely presented group `<gens | relators>` by Todd–Coxeter enumeration. The result counts live cosets with complete inverse transitions.",
    )?;
    def_doc(
        i,
        "fp_group_order",
        Arity::Fixed(2),
        fp_group_order,
        "Computes the order of the finitely presented group `<gens | relators>` via Todd-Coxeter.",
    )?;
    def_doc(
        i,
        "fp_group",
        Arity::Fixed(2),
        fp_group,
        "Constructs a `group(elements, table)` term for the finitely presented group `<gens | relators>` of order at most 720. The generator action satisfies the relators and has the same order as the enumerated coset count.",
    )?;
    def_doc(
        i,
        "fp_group_abelian_invariants",
        Arity::Fixed(2),
        fp_group_abelian_invariants,
        "Computes the abelianization invariant factors (torsion factors `> 1` followed by `0` per free factor) of `<gens | relators>` via Smith normal form.",
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;
    use super::*;

    fn cyc(n: usize) -> Perm {
        (0..n).map(|i| (i + 1) % n).collect()
    }

    fn swap01(n: usize) -> Perm {
        let mut p = identity(n);
        p.swap(0, 1);
        p
    }

    #[test]
    fn cayley_rows_match_independent_composition_in_every_pool() {
        let chain = Chain::new(5, &[cyc(5), swap01(5)]);
        let elements = chain.elements(TABLE_LIMIT).unwrap();
        assert_eq!(elements.len(), 120);
        let expected: Vec<Vec<_>> = elements
            .iter()
            .map(|p| {
                elements
                    .iter()
                    .map(|q| {
                        let product = compose(p, q);
                        elements.iter().position(|value| *value == product).unwrap()
                    })
                    .collect()
            })
            .collect();
        for workers in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            let actual = pool.install(|| cayley_rows(&elements)).unwrap();
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn abelian_invariants_of_presentations() {
        // C2 x C3 = C6; the Klein four-group; Z^2; Z x C4.
        assert_eq!(
            s("fp_group_abelian_invariants(list(a, b), list(list(a^2), list(b^3)))"),
            "list(6)"
        );
        assert_eq!(
            s(
                "fp_group_abelian_invariants(list(a, b), list(list(a^2), list(b^2), list(a, b, a^-1, b^-1)))"
            ),
            "list(2, 2)"
        );
        assert_eq!(
            s("fp_group_abelian_invariants(list(a, b), list())"),
            "list(0, 0)"
        );
        assert_eq!(
            s("fp_group_abelian_invariants(list(a, b), list(list(b^4)))"),
            "list(4, 0)"
        );
        // The binary icosahedral group <a, b | (ab)^2 = a^3 = b^5> is perfect.
        assert_eq!(
            s(
                "fp_group_abelian_invariants(list(a, b), list(list(a, b, a, b, a^-3), list(a, b, a, b, b^-5)))"
            ),
            "list()"
        );
    }

    #[test]
    fn symmetric_groups_have_factorial_order() {
        for n in 2..=9 {
            let c = Chain::new(n, &[cyc(n), swap01(n)]);
            let fact: usize = (1..=n).product();
            assert_eq!(c.order(), BigInt::from(fact), "S{n}");
        }
        assert_eq!(
            Chain::new(12, &[cyc(12), swap01(12)]).order(),
            BigInt::from(479_001_600_u64)
        );
    }

    #[test]
    fn mathieu_groups() {
        assert_eq!(
            s(
                "perm_group_order(list(list(2,3,4,5,6,7,8,9,10,11,1), list(1,2,7,10,6,4,11,3,9,5,8)))"
            ),
            "7920"
        );
        // M12 on twelve points
        let m12 = "list(list(2,3,4,5,6,7,8,9,10,11,1,12), list(1,2,7,10,6,4,11,3,9,5,8,12), list(12,11,6,8,9,3,10,4,5,7,2,1))";
        assert_eq!(s(&format!("perm_group_order({m12})")), "95040");
        assert_eq!(s(&format!("perm_group_is_transitive({m12})")), "true");
        // the point stabiliser of M12 is M11
        assert_eq!(
            s(&format!("perm_group_stabilizer_order({m12}, 12)")),
            "7920"
        );
        assert_eq!(s(&format!("perm_group_stabilizer_order({m12}, 1)")), "7920");
        assert_eq!(s(&format!("perm_group_is_solvable({m12})")), "false");
        assert_eq!(
            s(&format!("perm_group_derived_series({m12})")),
            "list(95040)"
        );
    }

    #[test]
    fn membership_and_base() {
        let a4 = "list(list(2,3,1,4), list(1,3,4,2))";
        assert_eq!(s(&format!("perm_group_order({a4})")), "12");
        assert_eq!(
            s(&format!("perm_group_contains({a4}, list(2,1,4,3))")),
            "true"
        );
        assert_eq!(
            s(&format!("perm_group_contains({a4}, list(2,1,3,4))")),
            "false"
        );
        assert_eq!(
            s(&format!("perm_group_contains({a4}, list(1,2,3))")),
            "false"
        );
        let c = Chain::new(4, &[vec![1, 2, 0, 3], vec![0, 2, 3, 1]]);
        // the product of the basic orbit lengths is the order
        assert_eq!(c.basic_orbit_lengths().iter().product::<usize>(), 12);
        assert!(c.base().len() >= 2);
        // the base has trivial pointwise stabiliser
        for g in c.elements(100).unwrap_or_default() {
            if c.base().iter().all(|&b| g[b] == b) {
                assert!(is_identity(&g));
            }
        }
        assert_eq!(s(&format!("perm_group_base({a4})")), "list(1, 2)");
    }

    #[test]
    fn elements_and_tables() {
        assert_eq!(s("perm_group_order(list(list(2,3,1), list(2,1,3)))"), "6");
        let e = s("perm_group_elements(list(list(2,3,1)))");
        assert_eq!(e, "list(list(1, 2, 3), list(2, 3, 1), list(3, 1, 2))");
        assert_eq!(
            s("group_order(group_from_perms(list(list(2,3,4,1), list(1,4,3,2))))"),
            "8"
        );
        assert_eq!(
            s("group_is_valid(group_from_perms(list(list(2,3,4,1), list(1,4,3,2))))"),
            "true"
        );
        assert_eq!(
            s("perm_group_orbits(list(list(2,1,3,5,4)))"),
            "list(list(1, 2), list(3), list(4, 5))"
        );
    }

    #[test]
    fn solvability_and_nilpotency() {
        // S4: derived series 24 > 12 > 4 > 1, lower central series 24 > 12
        let s4 = "list(list(2,3,4,1), list(2,1,3,4))";
        assert_eq!(
            s(&format!("perm_group_derived_series({s4})")),
            "list(24, 12, 4, 1)"
        );
        assert_eq!(s(&format!("perm_group_is_solvable({s4})")), "true");
        assert_eq!(
            s(&format!("perm_group_lower_central_series({s4})")),
            "list(24, 12)"
        );
        assert_eq!(s(&format!("perm_group_is_nilpotent({s4})")), "false");
        // S5 is not solvable; A5 is perfect
        let s5 = "list(list(2,3,4,5,1), list(2,1,3,4,5))";
        assert_eq!(
            s(&format!("perm_group_derived_series({s5})")),
            "list(120, 60)"
        );
        assert_eq!(s(&format!("perm_group_is_solvable({s5})")), "false");
        // D4 is nilpotent of class 2
        let d4 = "list(list(2,3,4,1), list(1,4,3,2))";
        assert_eq!(
            s(&format!("perm_group_lower_central_series({d4})")),
            "list(8, 2, 1)"
        );
        assert_eq!(s(&format!("perm_group_is_nilpotent({d4})")), "true");
        assert_eq!(
            s(&format!("perm_group_derived_series({d4})")),
            "list(8, 2, 1)"
        );
        assert_eq!(s(&format!("perm_group_is_abelian({d4})")), "false");
        assert_eq!(
            s("perm_group_is_abelian(list(list(2,1,3,4), list(1,2,4,3)))"),
            "true"
        );
        assert!(
            s(&format!("perm_group_derived_subgroup({s4})"))
                .matches("list(")
                .count()
                >= 2
        );
    }

    #[test]
    fn todd_coxeter_basics() {
        // <a, b | a^3, b^2, (ab)^2> = S3
        assert_eq!(
            s(
                "todd_coxeter_index(list(a, b), list(list(a, a, a), list(b, b), list(a, b, a, b)), list())"
            ),
            "6"
        );
        assert_eq!(
            s("fp_group_order(list(a, b), list(list(a^3), list(b^2), list(a, b, a, b)))"),
            "6"
        );
        // over the subgroup <b> the index is 3
        assert_eq!(
            s(
                "todd_coxeter_index(list(a, b), list(list(a, a, a), list(b, b), list(a, b, a, b)), list(list(b)))"
            ),
            "3"
        );
        // the action is a permutation representation of degree 3
        let act = s(
            "todd_coxeter(list(a, b), list(list(a, a, a), list(b, b), list(a, b, a, b)), list(list(b)))",
        );
        assert!(act.starts_with("list(3, list(list("), "{act}");
        // Z4 x Z2 and the quaternion group Q8
        assert_eq!(
            s("fp_group_order(list(a, b), list(list(a^4), list(b^2), list(a, b, -1, -2)))"),
            "8"
        );
        assert_eq!(
            s("fp_group_order(list(a, b), list(list(a^4), list(a, a, -2, -2), list(b, a, -2, a)))"),
            "8"
        );
        // the trivial and the free abelian quotient
        assert_eq!(s("fp_group_order(list(a), list(list(a)))"), "1");
        assert_eq!(s("fp_group_order(list(a), list(list(a^12)))"), "12");
        // an infinite group exceeds the coset limit and stays unreduced
        assert_eq!(
            s("fp_group_order(list(a), list())"),
            "fp_group_order(list(a), list())"
        );
    }

    #[test]
    fn coset_seal_rejects_missing_or_inconsistent_inverse_edges() {
        let mut table = Cosets::new(2);
        assert_eq!(table.seal(), None);
        let second = table.new_coset();
        table.set(0, 0, second);
        table.set(0, 1, 0);
        table.set(second, 0, 0);
        table.set(second, 1, second);
        assert_eq!(table.seal(), None);
        table.set(0, 1, second);
        table.set(second, 1, 0);
        assert_eq!(table.seal(), Some(()));
        assert_eq!(table.permutations(), [vec![1, 0]]);
    }

    #[test]
    fn coset_action_checks_relators_and_the_subgroup_at_the_base() {
        let relators = [vec![0, 0], vec![2, 2, 2], vec![0, 2, 0, 2]];
        let subgroup = [vec![0]];
        let table = enumerate_cosets(2, &relators, &subgroup).unwrap();
        assert_eq!(table.live, 3);
        assert!(table.verifies(&relators, &subgroup));
        assert!(!table.verifies(&[vec![2]], &[]));
        assert!(!table.verifies(&[], &[vec![2]]));
        assert!(reduced_word(&[0, 2, 3, 1]).is_empty());
        assert!(enumerate_cosets(2, &[vec![4]], &[]).is_none());
        assert!(enumerate_cosets(usize::MAX, &[], &[]).is_none());
    }

    #[test]
    fn fp_group_is_a_group_term() {
        let g = "fp_group(list(a, b), list(list(a^3), list(b^2), list(a, b, a, b)))";
        assert_eq!(s(&format!("group_order({g})")), "6");
        assert_eq!(s(&format!("group_is_abelian({g})")), "false");
        assert_eq!(s(&format!("group_is_valid({g})")), "true");
    }

    #[test]
    fn psl27() {
        // PSL(2,7) = <a, b | a^2, b^3, (ab)^7, [a,b]^4> of order 168
        let words: Vec<Vec<usize>> = vec![
            vec![0, 0],
            vec![2, 2, 2],
            std::iter::repeat_n([0, 2], 7).flatten().collect(),
            std::iter::repeat_n([0, 2, 1, 3], 4).flatten().collect(),
        ];
        // here columns: a = 0, a^-1 = 1, b = 2, b^-1 = 3; b^3 = (2,2,2)
        let (index, perms) = todd_coxeter_table(2, &words, &[]).unwrap_or((0, Vec::new()));
        assert_eq!(index, 168);
        assert_eq!(perms.len(), 2);
    }
}
