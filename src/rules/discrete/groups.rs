//! Finite groups, their representations and permutations.
//!
//! A finite group is the inert term `group(list(e1, ..., en), table)` where
//! `table` is the Cayley table `list(list(...), ...)` whose entry `(i, j)` is
//! the element `ei * ej`. Elements are arbitrary terms (integers, symbols,
//! permutation lists); two elements are equal when the engine knows them to
//! be equal. Table entries that are not syntactically an element are
//! simplified once before they are looked up.
//!
//! A *representation* is a list of square matrices (`list(list(..), ..)`)
//! aligned with the element list of the group. A *permutation* of `n` points
//! is the list of the images of `1, ..., n`: `list(1, 2, 3)` is the identity
//! and the product is `(p q)(i) = p(q(i))`.
//!
//! | operator | value |
//! |---|---|
//! | `cyclic_group(n)` | `Z_n`: elements `0..n-1`, addition mod `n` |
//! | `dihedral_group(n)` | `D_n` of order `2n`: `r0..r(n-1)`, `s0..s(n-1)` with `r_a r_b = r_(a+b)`, `r_a s_b = s_(b-a)`, `s_a r_b = s_(a+b)`, `s_a s_b = r_(b-a)` |
//! | `symmetric_group(n)` | `S_n` (`n <= 6`), permutations in lexicographic order |
//! | `klein_four_group()` | elements `e, a, b, c` |
//! | `group_from_table(elements, table)` | the group term, if the table is a valid group |
//! | `group_elements(G)`, `group_order(G)`, `group_identity(G)` | the element list, its length, the identity |
//! | `group_mul(G, a, b)`, `group_inverse(G, a)`, `group_element_order(G, a)` | product, inverse, order of an element |
//! | `group_is_abelian(G)`, `group_is_valid(G)` | truth values (validity: closure, associativity, identity, inverses) |
//! | `group_conjugacy_classes(G)`, `group_center(G)` | list of classes, centre |
//! | `group_subgroups(G)` | all subgroups as element lists (orders up to 48) |
//! | `group_cosets(G, H)`, `group_is_normal(G, H)` | left cosets of a subgroup, normality test |
//! | `representation_is_valid(G, mats)` | `M(ab) = M(a) M(b)` for all pairs |
//! | `group_character(mats)` | list of the traces |
//! | `perm_compose(p, q)`, `perm_inverse(p)`, `perm_order(p)`, `perm_cycles(p)`, `perm_sign(p)` | permutation arithmetic; `perm_cycles` lists the non-trivial cycles |

use std::collections::BTreeSet;
use std::collections::HashMap;

use num_bigint::BigInt;
use num_integer::Integer;

use super::apply;
use super::def;
use super::def_inert;
use super::idx;
use super::idxs;
use super::items;
use super::matmul;
use super::neg;
use super::rows;
use super::sum;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::ClassId;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::RuleError;

/// The largest group order for which all subgroups are enumerated.
const SUBGROUP_LIMIT: usize = 48;

/// A group read from a term, before validation.
struct Raw {
    elems: Vec<NodeId>,
    index: HashMap<ClassId, usize>,
    distinct: bool,
    table: Vec<Vec<Option<usize>>>,
}

/// A group whose table is square and made of elements.
struct Grp {
    elems: Vec<NodeId>,
    index: HashMap<ClassId, usize>,
    table: Vec<Vec<usize>>,
}

/// The position of the element `n`, simplifying it if necessary.
fn lookup(
    cx: &mut Cx<'_>,
    index: &HashMap<ClassId, usize>,
    n: NodeId,
) -> Option<usize> {
    if let Some(&i) = index.get(&cx.graph.find(n)) {
        return Some(i);
    }
    let s = cx.simplify(n);
    index.get(&cx.graph.find(s)).copied()
}

fn raw_from(
    cx: &mut Cx<'_>,
    elems: Vec<NodeId>,
    table: &[Vec<NodeId>],
) -> Raw {
    let mut index = HashMap::new();
    let mut distinct = true;
    for (i, &e) in elems.iter().enumerate() {
        if index.insert(cx.graph.find(e), i).is_some() {
            distinct = false;
        }
    }
    let table = table
        .iter()
        .map(|row| row.iter().map(|&e| lookup(cx, &index, e)).collect())
        .collect();
    Raw { elems, index, distinct, table }
}

fn find_group(
    g: &Graph,
    n: NodeId,
) -> Option<NodeId> {
    let op = g.ops().lookup("group")?;
    if g.op(n) == op {
        return Some(n);
    }
    g.members(g.find(n)).find(|&e| g.op(e) == op)
}

fn read_raw(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<Raw> {
    let gn = find_group(cx.graph, node)?;
    let children = cx.graph.children(gn).to_vec();
    let [el, tb] = children.as_slice() else { return None };
    let elems = items(cx.graph, *el)?;
    let table = rows(cx.graph, *tb)?;
    Some(raw_from(cx, elems, &table))
}

impl Raw {
    /// The group, if the elements are distinct and the table is a square
    /// array of elements.
    fn checked(self) -> Option<Grp> {
        let n = self.elems.len();
        if !self.distinct || self.table.len() != n {
            return None;
        }
        let mut table = Vec::with_capacity(n);
        for row in self.table {
            if row.len() != n {
                return None;
            }
            table.push(row.into_iter().collect::<Option<Vec<usize>>>()?);
        }
        Some(Grp { elems: self.elems, index: self.index, table })
    }
}

fn read_group(
    cx: &mut Cx<'_>,
    node: NodeId,
) -> Option<Grp> {
    read_raw(cx, node)?.checked()
}

impl Grp {
    const fn n(&self) -> usize {
        self.elems.len()
    }

    fn identity(&self) -> Option<usize> {
        (0..self.n()).find(|&e| (0..self.n()).all(|j| self.table[e][j] == j && self.table[j][e] == j))
    }

    fn inverse(
        &self,
        a: usize,
    ) -> Option<usize> {
        let e = self.identity()?;
        (0..self.n()).find(|&j| self.table[a][j] == e && self.table[j][a] == e)
    }

    fn is_abelian(&self) -> bool {
        (0..self.n()).all(|a| (0..self.n()).all(|b| self.table[a][b] == self.table[b][a]))
    }

    fn is_associative(&self) -> bool {
        let n = self.n();
        let t = &self.table;
        (0..n).all(|a| (0..n).all(|b| (0..n).all(|c| t[t[a][b]][c] == t[a][t[b][c]])))
    }

    fn is_valid(&self) -> bool {
        self.identity().is_some() && (0..self.n()).all(|a| self.inverse(a).is_some()) && self.is_associative()
    }

    /// The order of `a` (smallest `k >= 1` with `a^k = e`).
    fn element_order(
        &self,
        a: usize,
    ) -> Option<usize> {
        let e = self.identity()?;
        let mut cur = a;
        for k in 1..=self.n() {
            if cur == e {
                return Some(k);
            }
            cur = self.table[cur][a];
        }
        None
    }

    /// The subgroup generated by `gens` (sorted), by closure under
    /// right multiplication by the generators.
    fn closure(
        &self,
        gens: &[usize],
    ) -> Option<Vec<usize>> {
        let e = self.identity()?;
        let mut seen = vec![false; self.n()];
        seen[e] = true;
        let mut stack = vec![e];
        while let Some(x) = stack.pop() {
            for &g in gens {
                let y = self.table[x][g];
                if !seen[y] {
                    seen[y] = true;
                    stack.push(y);
                }
            }
        }
        Some((0..self.n()).filter(|&i| seen[i]).collect())
    }

    fn is_subgroup(
        &self,
        h: &[usize],
    ) -> bool {
        let Some(e) = self.identity() else { return false };
        h.contains(&e) && h.iter().all(|&a| h.iter().all(|&b| h.contains(&self.table[a][b])))
    }

    fn elems_v(
        &self,
        ids: &[usize],
    ) -> V {
        V::List(ids.iter().map(|&i| V::Node(self.elems[i])).collect())
    }

    /// Reads a list of elements as sorted, distinct positions.
    fn read_subset(
        &self,
        cx: &mut Cx<'_>,
        n: NodeId,
    ) -> Option<Vec<usize>> {
        let mut out = BTreeSet::new();
        for e in items(cx.graph, n)? {
            out.insert(lookup(cx, &self.index, e)?);
        }
        Some(out.into_iter().collect())
    }
}

/// The group term for `elems` and the table of positions.
fn build_group(
    g: &mut Graph,
    elems: &[NodeId],
    table: &[Vec<usize>],
) -> Option<V> {
    let el = super::list(g, elems);
    let rows: Vec<Vec<NodeId>> = table.iter().map(|r| r.iter().map(|&i| elems[i]).collect()).collect();
    let tb = super::matrix(g, &rows);
    apply(g, "group", &[el, tb]).map(V::Node)
}

// ----------------------------------------------------------------------
// Constructors
// ----------------------------------------------------------------------

/// A group size argument within `1..=limit`.
fn size(
    cx: &Cx<'_>,
    a: &[NodeId],
    limit: usize,
) -> Option<usize> {
    let n = idx(cx.graph, *a.first()?)?;
    (1..=limit).contains(&n).then_some(n)
}

fn cyclic_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = size(cx, a, 256)?;
    let elems: Vec<NodeId> = (0..n).map(|i| cx.graph.int(i64::try_from(i).unwrap_or(0))).collect();
    let table: Vec<Vec<usize>> = (0..n).map(|i| (0..n).map(|j| (i + j) % n).collect()).collect();
    build_group(cx.graph, &elems, &table)
}

fn dihedral_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = size(cx, a, 128)?;
    let mut elems: Vec<NodeId> = (0..n).map(|i| cx.graph.sym(&format!("r{i}"))).collect();
    elems.extend((0..n).map(|i| cx.graph.sym(&format!("s{i}"))));
    let sub = |x: usize, y: usize| (y + n - x) % n;
    let mut table = vec![vec![0; 2 * n]; 2 * n];
    for i in 0..n {
        for j in 0..n {
            table[i][j] = (i + j) % n;
            table[i][n + j] = n + sub(i, j);
            table[n + i][j] = n + (i + j) % n;
            table[n + i][n + j] = sub(i, j);
        }
    }
    build_group(cx.graph, &elems, &table)
}

/// All permutations of `0..n` in lexicographic order.
fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut p: Vec<usize> = (0..n).collect();
    loop {
        out.push(p.clone());
        let Some(i) = (1..n).rev().find(|&i| p[i - 1] < p[i]) else { break };
        let j = (i..n).rev().find(|&j| p[j] > p[i - 1]).unwrap_or(i);
        p.swap(i - 1, j);
        p[i..].reverse();
    }
    out
}

fn symmetric_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = size(cx, a, 6)?;
    let perms = permutations(n);
    let pos: HashMap<&Vec<usize>, usize> = perms.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let mut elems = Vec::with_capacity(perms.len());
    for p in &perms {
        let images: Vec<NodeId> = p.iter().map(|&x| cx.graph.int(i64::try_from(x + 1).unwrap_or(0))).collect();
        elems.push(super::list(cx.graph, &images));
    }
    let mut table = Vec::with_capacity(perms.len());
    for p in &perms {
        let mut row = Vec::with_capacity(perms.len());
        for q in &perms {
            let r: Vec<usize> = q.iter().map(|&x| p[x]).collect();
            row.push(*pos.get(&r)?);
        }
        table.push(row);
    }
    build_group(cx.graph, &elems, &table)
}

fn klein_four_group(
    cx: &mut Cx<'_>,
    _a: &[NodeId],
) -> Option<V> {
    let elems: Vec<NodeId> = ["e", "a", "b", "c"].iter().map(|s| cx.graph.sym(s)).collect();
    // The Klein group is Z2 x Z2: position bits xor.
    let table: Vec<Vec<usize>> = (0..4).map(|i| (0..4).map(|j| i ^ j).collect()).collect();
    build_group(cx.graph, &elems, &table)
}

fn group_from_table(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [el, tb] = a else { return None };
    let elems = items(cx.graph, *el)?;
    let elems = super::simplified(cx, &elems);
    let table = rows(cx.graph, *tb)?;
    let grp = raw_from(cx, elems, &table).checked()?;
    if !grp.is_valid() {
        return None;
    }
    build_group(cx.graph, &grp.elems, &grp.table)
}

// ----------------------------------------------------------------------
// Queries
// ----------------------------------------------------------------------

fn group_elements(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let grp = read_group(cx, *a.first()?)?;
    Some(V::nodes(&grp.elems))
}

fn group_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(read_group(cx, *a.first()?)?.n()))
}

fn group_identity(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let grp = read_group(cx, *a.first()?)?;
    Some(V::Node(grp.elems[grp.identity()?]))
}

fn group_mul(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, x, y] = a else { return None };
    let grp = read_group(cx, *g)?;
    let (x, y) = (lookup(cx, &grp.index, *x)?, lookup(cx, &grp.index, *y)?);
    Some(V::Node(grp.elems[grp.table[x][y]]))
}

fn group_inverse(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, x] = a else { return None };
    let grp = read_group(cx, *g)?;
    let x = lookup(cx, &grp.index, *x)?;
    Some(V::Node(grp.elems[grp.inverse(x)?]))
}

fn group_is_abelian(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::Bool(read_group(cx, *a.first()?)?.is_abelian()))
}

fn group_element_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, x] = a else { return None };
    let grp = read_group(cx, *g)?;
    let x = lookup(cx, &grp.index, *x)?;
    Some(V::uint(grp.element_order(x)?))
}

fn group_conjugacy_classes(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let grp = read_group(cx, *a.first()?)?;
    let inv: Vec<usize> = (0..grp.n()).map(|x| grp.inverse(x)).collect::<Option<_>>()?;
    let mut done = vec![false; grp.n()];
    let mut classes = Vec::new();
    for x in 0..grp.n() {
        if done[x] {
            continue;
        }
        let class: BTreeSet<usize> = (0..grp.n()).map(|g| grp.table[grp.table[g][x]][inv[g]]).collect();
        for &c in &class {
            done[c] = true;
        }
        let class: Vec<usize> = class.into_iter().collect();
        classes.push(grp.elems_v(&class));
    }
    Some(V::List(classes))
}

fn group_center(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let grp = read_group(cx, *a.first()?)?;
    let center: Vec<usize> =
        (0..grp.n()).filter(|&z| (0..grp.n()).all(|g| grp.table[z][g] == grp.table[g][z])).collect();
    Some(grp.elems_v(&center))
}

fn group_is_valid(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let raw = read_raw(cx, *a.first()?)?;
    Some(V::Bool(raw.checked().is_some_and(|g| g.is_valid())))
}

fn group_subgroups(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let grp = read_group(cx, *a.first()?)?;
    if grp.n() > SUBGROUP_LIMIT {
        return None;
    }
    let trivial = grp.closure(&[])?;
    let mut found: BTreeSet<(usize, Vec<usize>)> = BTreeSet::new();
    found.insert((trivial.len(), trivial.clone()));
    let mut queue = vec![trivial];
    while let Some(h) = queue.pop() {
        for g in 0..grp.n() {
            if h.contains(&g) {
                continue;
            }
            let mut gens = h.clone();
            gens.push(g);
            let k = grp.closure(&gens)?;
            if found.insert((k.len(), k.clone())) {
                queue.push(k);
            }
        }
    }
    Some(V::List(found.iter().map(|(_, h)| grp.elems_v(h)).collect()))
}

fn group_cosets(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, h] = a else { return None };
    let grp = read_group(cx, *g)?;
    let h = grp.read_subset(cx, *h)?;
    if !grp.is_subgroup(&h) {
        return None;
    }
    let mut done = vec![false; grp.n()];
    let mut cosets = Vec::new();
    for x in 0..grp.n() {
        if done[x] {
            continue;
        }
        let coset: BTreeSet<usize> = h.iter().map(|&y| grp.table[x][y]).collect();
        for &c in &coset {
            done[c] = true;
        }
        let coset: Vec<usize> = coset.into_iter().collect();
        cosets.push(grp.elems_v(&coset));
    }
    Some(V::List(cosets))
}

fn group_is_normal(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, h] = a else { return None };
    let grp = read_group(cx, *g)?;
    let h = grp.read_subset(cx, *h)?;
    if !grp.is_subgroup(&h) {
        return None;
    }
    let inv: Vec<usize> = (0..grp.n()).map(|x| grp.inverse(x)).collect::<Option<_>>()?;
    let normal = (0..grp.n()).all(|x| h.iter().all(|&y| h.contains(&grp.table[grp.table[x][y]][inv[x]])));
    Some(V::Bool(normal))
}

// ----------------------------------------------------------------------
// Representations
// ----------------------------------------------------------------------

fn read_matrices(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<Vec<Vec<NodeId>>>> {
    items(g, n)?.into_iter().map(|m| rows(g, m)).collect()
}

/// Whether `a` and `b` are equal terms or differ by zero.
fn same_entry(
    cx: &mut Cx<'_>,
    a: NodeId,
    b: NodeId,
) -> bool {
    if cx.graph.same(a, b) {
        return true;
    }
    let nb = neg(cx.graph, b);
    let diff = sum(cx.graph, &[a, nb]);
    cx.is_zero(diff)
}

fn representation_is_valid(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [g, m] = a else { return None };
    let grp = read_group(cx, *g)?;
    let mats = read_matrices(cx.graph, *m)?;
    if mats.len() != grp.n() {
        return None;
    }
    for i in 0..grp.n() {
        for j in 0..grp.n() {
            let target = &mats[grp.table[i][j]];
            let Some(product) = matmul(cx, &mats[i], &mats[j]) else { return Some(V::Bool(false)) };
            if product.len() != target.len() {
                return Some(V::Bool(false));
            }
            for (pr, tr) in product.iter().zip(target) {
                if pr.len() != tr.len() {
                    return Some(V::Bool(false));
                }
                for (&x, &y) in pr.iter().zip(tr) {
                    if !same_entry(cx, x, y) {
                        return Some(V::Bool(false));
                    }
                }
            }
        }
    }
    Some(V::Bool(true))
}

fn group_character(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let mats = read_matrices(cx.graph, *a.first()?)?;
    let mut traces = Vec::with_capacity(mats.len());
    for m in &mats {
        if m.iter().any(|r| r.len() != m.len()) {
            return None;
        }
        let diag: Vec<NodeId> = m.iter().enumerate().map(|(i, r)| r[i]).collect();
        let t = sum(cx.graph, &diag);
        traces.push(V::Node(cx.simplify(t)));
    }
    Some(V::List(traces))
}

// ----------------------------------------------------------------------
// Permutations
// ----------------------------------------------------------------------

/// A permutation as 0-based images; `None` unless it is a bijection of
/// `1..=n`.
fn read_perm(
    g: &Graph,
    n: NodeId,
) -> Option<Vec<usize>> {
    let images = idxs(g, n)?;
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

fn perm_v(p: &[usize]) -> V {
    V::ints(p.iter().map(|&x| i64::try_from(x + 1).unwrap_or(0)))
}

/// The cycles of `p`, fixed points included, by smallest element.
fn cycles(p: &[usize]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; p.len()];
    let mut out = Vec::new();
    for start in 0..p.len() {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut x = start;
        while !seen[x] {
            seen[x] = true;
            cycle.push(x);
            x = p[x];
        }
        out.push(cycle);
    }
    out
}

fn perm_compose(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let [p, q] = a else { return None };
    let (p, q) = (read_perm(cx.graph, *p)?, read_perm(cx.graph, *q)?);
    if p.len() != q.len() {
        return None;
    }
    let r: Vec<usize> = q.iter().map(|&x| p[x]).collect();
    Some(perm_v(&r))
}

fn perm_inverse(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let p = read_perm(cx.graph, *a.first()?)?;
    let mut inv = vec![0; p.len()];
    for (i, &x) in p.iter().enumerate() {
        inv[x] = i;
    }
    Some(perm_v(&inv))
}

fn perm_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let p = read_perm(cx.graph, *a.first()?)?;
    let order = cycles(&p).iter().fold(BigInt::from(1), |acc, c| acc.lcm(&BigInt::from(c.len())));
    Some(V::Int(order))
}

fn perm_cycles(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let p = read_perm(cx.graph, *a.first()?)?;
    Some(V::List(cycles(&p).iter().filter(|c| c.len() > 1).map(|c| perm_v(c)).collect()))
}

fn perm_sign(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let p = read_perm(cx.graph, *a.first()?)?;
    Some(V::int(if (p.len() - cycles(&p).len()).is_multiple_of(2) { 1 } else { -1 }))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def_inert(i, "group", Arity::Fixed(2))?;
    def(i, "cyclic_group", Arity::Fixed(1), cyclic_group)?;
    def(i, "dihedral_group", Arity::Fixed(1), dihedral_group)?;
    def(i, "symmetric_group", Arity::Fixed(1), symmetric_group)?;
    def(i, "klein_four_group", Arity::Fixed(0), klein_four_group)?;
    def(i, "group_from_table", Arity::Fixed(2), group_from_table)?;
    def(i, "group_elements", Arity::Fixed(1), group_elements)?;
    def(i, "group_order", Arity::Fixed(1), group_order)?;
    def(i, "group_identity", Arity::Fixed(1), group_identity)?;
    def(i, "group_mul", Arity::Fixed(3), group_mul)?;
    def(i, "group_inverse", Arity::Fixed(2), group_inverse)?;
    def(i, "group_is_abelian", Arity::Fixed(1), group_is_abelian)?;
    def(i, "group_element_order", Arity::Fixed(2), group_element_order)?;
    def(i, "group_conjugacy_classes", Arity::Fixed(1), group_conjugacy_classes)?;
    def(i, "group_center", Arity::Fixed(1), group_center)?;
    def(i, "group_is_valid", Arity::Fixed(1), group_is_valid)?;
    def(i, "group_subgroups", Arity::Fixed(1), group_subgroups)?;
    def(i, "group_cosets", Arity::Fixed(2), group_cosets)?;
    def(i, "group_is_normal", Arity::Fixed(2), group_is_normal)?;
    def(i, "representation_is_valid", Arity::Fixed(2), representation_is_valid)?;
    def(i, "group_character", Arity::Fixed(1), group_character)?;
    def(i, "perm_compose", Arity::Fixed(2), perm_compose)?;
    def(i, "perm_inverse", Arity::Fixed(1), perm_inverse)?;
    def(i, "perm_order", Arity::Fixed(1), perm_order)?;
    def(i, "perm_cycles", Arity::Fixed(1), perm_cycles)?;
    def(i, "perm_sign", Arity::Fixed(1), perm_sign)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;

    #[test]
    fn constructors() {
        assert_eq!(
            s("cyclic_group(3)"),
            "group(list(0, 1, 2), list(list(0, 1, 2), list(1, 2, 0), list(2, 0, 1)))"
        );
        assert_eq!(s("group_order(dihedral_group(4))"), "8");
        assert_eq!(s("group_order(symmetric_group(4))"), "24");
        assert_eq!(s("group_elements(klein_four_group())"), "list(e, a, b, c)");
        assert_eq!(s("group_mul(klein_four_group(), a, b)"), "c");
        assert_eq!(s("group_mul(dihedral_group(5), r2, s1)"), "s4");
        assert_eq!(s("group_mul(dihedral_group(5), s1, r2)"), "s3");
        assert_eq!(s("group_mul(dihedral_group(5), s1, s3)"), "r2");
        assert_eq!(s("group_identity(symmetric_group(3))"), "list(1, 2, 3)");
        // refused sizes stay unreduced
        assert_eq!(s("group_order(symmetric_group(7))"), "group_order(symmetric_group(7))");
        assert_eq!(s("cyclic_group(0)"), "cyclic_group(0)");
    }

    #[test]
    fn from_table_validates() {
        assert_eq!(
            s("group_is_valid(group_from_table(list(e, x), list(list(e, x), list(x, e))))"),
            "true"
        );
        assert_eq!(
            s("group_order(group_from_table(list(e, x), list(list(e, x), list(x, e))))"),
            "2"
        );
        // not associative / no identity: stays unreduced
        let bad = "group_from_table(list(a, b), list(list(b, a), list(a, a)))";
        assert_eq!(s(bad), bad);
        assert_eq!(s(&format!("group_is_valid(group({}))", &bad["group_from_table(".len()..bad.len() - 1])), "false");
        // entries are simplified
        assert_eq!(
            s("group_order(group_from_table(list(0, 1), list(list(0, 1), list(1, 1 + 1 - 2))))"),
            "2"
        );
    }

    #[test]
    fn element_queries() {
        assert_eq!(s("group_inverse(cyclic_group(5), 2)"), "3");
        assert_eq!(s("group_inverse(dihedral_group(3), r1)"), "r2");
        assert_eq!(s("group_inverse(dihedral_group(3), s1)"), "s1");
        assert_eq!(s("group_element_order(cyclic_group(6), 2)"), "3");
        assert_eq!(s("group_element_order(dihedral_group(4), r1)"), "4");
        assert_eq!(s("group_element_order(dihedral_group(4), s2)"), "2");
        assert_eq!(s("group_element_order(symmetric_group(3), list(2, 3, 1))"), "3");
        assert_eq!(s("group_is_abelian(cyclic_group(6))"), "true");
        assert_eq!(s("group_is_abelian(klein_four_group())"), "true");
        assert_eq!(s("group_is_abelian(dihedral_group(3))"), "false");
        assert_eq!(s("group_is_abelian(symmetric_group(3))"), "false");
        assert_eq!(s("group_is_valid(symmetric_group(3))"), "true");
        assert_eq!(s("group_is_valid(dihedral_group(5))"), "true");
    }

    #[test]
    fn classes_and_centre() {
        assert_eq!(
            s("group_conjugacy_classes(dihedral_group(3))"),
            "list(list(r0), list(r1, r2), list(s0, s1, s2))"
        );
        assert_eq!(s("group_conjugacy_classes(klein_four_group())"), "list(list(e), list(a), list(b), list(c))");
        assert_eq!(s("group_center(dihedral_group(4))"), "list(r0, r2)");
        assert_eq!(s("group_center(dihedral_group(3))"), "list(r0)");
        assert_eq!(s("group_center(cyclic_group(3))"), "list(0, 1, 2)");
    }

    #[test]
    fn subgroups_and_cosets() {
        // D3 has 6 subgroups: 1, three of order 2, A3, D3
        let subs = s("group_subgroups(dihedral_group(3))");
        assert_eq!(subs.matches("list(").count(), 7);
        assert!(subs.starts_with("list(list(r0), "));
        // S4 has 30 subgroups
        let subs = s("group_subgroups(symmetric_group(4))");
        assert!(subs.starts_with("list(list(list(1, 2, 3, 4))"), "{subs}");
        assert_eq!(s("group_subgroups(cyclic_group(6))"), "list(list(0), list(0, 3), list(0, 2, 4), list(0, 1, 2, 3, 4, 5))");
        assert_eq!(s("group_cosets(cyclic_group(6), list(0, 3))"), "list(list(0, 3), list(1, 4), list(2, 5))");
        assert_eq!(s("group_cosets(dihedral_group(3), list(r0, s0))"), "list(list(r0, s0), list(r1, s2), list(r2, s1))");
        // not a subgroup: unreduced
        assert!(s("group_cosets(cyclic_group(6), list(0, 1))").starts_with("group_cosets(group("));
    }

    #[test]
    fn normality() {
        assert_eq!(s("group_is_normal(dihedral_group(3), list(r0, r1, r2))"), "true");
        assert_eq!(s("group_is_normal(dihedral_group(3), list(r0, s0))"), "false");
        assert_eq!(s("group_is_normal(cyclic_group(6), list(0, 2, 4))"), "true");
        assert_eq!(s("group_is_normal(dihedral_group(4), list(r0, r2))"), "true");
    }

    #[test]
    fn representations() {
        // the sign representation of S3 on the 1x1 matrices
        let sign = "list(list(list(1)), list(list(-1)), list(list(-1)), list(list(-1)), list(list(1)), list(list(1)))";
        // lexicographic order: 123, 132, 213, 231, 312, 321
        let sign = sign.replace("list(list(1)), list(list(-1)), list(list(-1)), list(list(-1)), list(list(1)), list(list(1))",
            "list(list(1)), list(list(-1)), list(list(-1)), list(list(1)), list(list(1)), list(list(-1))");
        assert_eq!(s(&format!("representation_is_valid(symmetric_group(3), {sign})")), "true");
        assert_eq!(s(&format!("group_character({sign})")), "list(1, -1, -1, 1, 1, -1)");
        let wrong = "list(list(list(1)), list(list(-1)), list(list(-1)), list(list(-1)), list(list(1)), list(list(1)))";
        assert_eq!(s(&format!("representation_is_valid(symmetric_group(3), {wrong})")), "false");
        // the regular representation of Z2 by 2x2 matrices, with a symbolic entry
        let rep = "list(list(list(1, 0), list(0, 1)), list(list(0, 1), list(1, 0)))";
        assert_eq!(s(&format!("representation_is_valid(cyclic_group(2), {rep})")), "true");
        assert_eq!(s(&format!("group_character({rep})")), "list(2, 0)");
        let sym = "list(list(list(x, 0), list(0, x)), list(list(1, 0), list(0, 1)))";
        assert_eq!(s(&format!("group_character({sym})")), "list(2*x, 2)");
    }

    #[test]
    fn permutations() {
        assert_eq!(s("perm_compose(list(2, 3, 1), list(1, 3, 2))"), "list(2, 1, 3)");
        assert_eq!(s("perm_inverse(list(2, 3, 1))"), "list(3, 1, 2)");
        assert_eq!(s("perm_order(list(2, 1, 4, 5, 3))"), "6");
        assert_eq!(s("perm_cycles(list(2, 1, 4, 5, 3, 6))"), "list(list(1, 2), list(3, 4, 5))");
        assert_eq!(s("perm_cycles(list(1, 2))"), "list()");
        assert_eq!(s("perm_sign(list(2, 3, 1))"), "1");
        assert_eq!(s("perm_sign(list(2, 1, 3))"), "-1");
        // not permutations
        assert_eq!(s("perm_order(list(1, 1))"), "perm_order(list(1, 1))");
        assert_eq!(s("perm_compose(list(1, 2), list(1, 2, 3))"), "perm_compose(list(1, 2), list(1, 2, 3))");
    }
}
