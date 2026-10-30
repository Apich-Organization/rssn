//! # Canonical labelling, partition refinement and automorphism groups
//!
//! The layer has three parts, each usable on its own:
//!
//! 1. [`ColourGraph`] — a vertex-coloured graph with `u32` edge labels, the
//!    common input. Parallel edges of a term-level graph are recorded as
//!    labels, so nothing is lost by dropping repeated arcs.
//! 2. [`refine`] — 1-dimensional Weisfeiler–Leman / equitable-partition
//!    refinement, with a [`RefinementCertificate`] that is *equal for
//!    isomorphic graphs* (a necessary condition, never a sufficient one).
//! 3. [`canonical_labelling`] — a nauty-style individualisation–refinement
//!    search returning the canonical adjacency matrix, the automorphism group
//!    generators, their orbits and the group order.
//!
//! # Partition refinement
//!
//! The initial partition is given by the vertex colours (sorted by colour
//! value, so the partition is a canonical function of the input). A round
//! replaces the colour of every vertex by
//! `(cell, multiset of (cell, label) over out-neighbours, multiset over
//! in-neighbours)` and sorts the vertices of each cell by that signature; a
//! cell splits where its signatures differ. Names are assigned in
//! lexicographic order of the signatures, which makes the colouring of every
//! round invariant under a renaming of the vertices, hence the certificate.
//! The process stops at an *equitable* partition: every two vertices of a
//! cell have the same number of neighbours (per label) in every cell.
//!
//! # Individualisation–refinement
//!
//! The search keeps an ordered partition. At a node it refines the partition;
//! if the partition is discrete the flattened order is a candidate labelling,
//! whose *code* — the vertex colours followed by the adjacency matrix in that
//! order — is compared with the best code seen so far. Otherwise the first
//! non-singleton cell is the target and every vertex of the target is
//! individualised in turn (its cell is replaced by the singleton and the
//! rest), which is the branching rule of `nauty`.
//!
//! Two prunings keep the tree small:
//!
//! * **Automorphism pruning.** Every leaf whose code equals the best code
//!   yields an automorphism `best_order[i] -> order[i]`; those permutations
//!   generate a subgroup of `Aut(G)` (a stabiliser chain is maintained so
//!   redundant generators are dropped). At a node, a target vertex that lies
//!   in the orbit — under the found generators that fix the individualised
//!   prefix pointwise — of an *already branched* vertex is skipped: its
//!   subtree is the image of the branched one and cannot lower the code nor
//!   produce a new generator. The search therefore explores one representative
//!   per orbit instead of one per vertex.
//! * **Leaf pruning.** A candidate code is built and compared row by row and
//!   abandoned at the first row that is already larger than the best code.
//!
//! Because the tree and the codes are functions of the isomorphism class
//! only, the minimum code — the *canonical form* — is the same for every
//! input ordering of an isomorphic graph, and two graphs are isomorphic if
//! and only if their canonical codes agree. (The canonical *order* is unique
//! only up to an automorphism: when the group is non-trivial several leaves
//! realise the minimum, exactly as in nauty.)
//!
//! # Budget and fallback
//!
//! The search is bounded by a [`SearchBudget`]. When the budget runs out the
//! search unwinds, [`CanonicalLabelling::is_exact`] is `false`, the returned
//! order and code are only the best found (a valid labelling, *not*
//! guaranteed to be canonical) and [`CanonicalLabelling::group_order`] is the
//! order of the subgroup found so far (a divisor of the true order). Callers
//! that need a guaranteed answer must test `is_exact()`; the term-level
//! operators in `rules::discrete::graphs` fall back to the exhaustive
//! matcher for graphs of at most 64 vertices and otherwise leave the request
//! unreduced.

use std::cmp::Ordering;

use num_bigint::BigInt;

/// A vertex-coloured graph with integer edge labels.
///
/// Vertices are `0..n`, every vertex carries a `u32` colour, and an arc
/// carries a non-zero `u32` label. An undirected graph stores each edge in
/// both directions, so `label(u, v) == label(v, u)` holds for every pair.
/// A self-loop is stored the same way as any other arc.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColourGraph {
    n: usize,
    directed: bool,
    /// `matrix[u][v]`, the label of the arc `u -> v`; `0` means "no arc".
    matrix: Vec<Vec<u32>>,
    /// Out-arcs of each vertex, as `(head, label)`, sorted by head.
    out: Vec<Vec<(usize, u32)>>,
    /// In-arcs of each vertex, as `(tail, label)`; empty when undirected.
    incoming: Vec<Vec<(usize, u32)>>,
    colours: Vec<u32>,
}

impl ColourGraph {
    /// A graph on `n` vertices with no edges and the single colour class `0`.
    #[must_use]
    pub fn new(
        n: usize,
        directed: bool,
    ) -> Self {
        Self {
            n,
            directed,
            matrix: vec![vec![0; n]; n],
            out: vec![Vec::new(); n],
            incoming: vec![Vec::new(); n],
            colours: vec![0; n],
        }
    }

    /// A graph on `n` vertices with the given colours.
    ///
    /// Returns `None` when `colours` does not have one entry per vertex.
    #[must_use]
    pub fn with_colours(
        n: usize,
        directed: bool,
        colours: &[u32],
    ) -> Option<Self> {
        let mut graph = Self::new(n, directed);
        graph.set_colours(colours).then_some(graph)
    }

    /// Replaces the vertex colours.
    ///
    /// Returns `false` when `colours` does not have one entry per vertex.
    pub fn set_colours(
        &mut self,
        colours: &[u32],
    ) -> bool {
        if colours.len() != self.n {
            return false;
        }
        self.colours.copy_from_slice(colours);
        true
    }

    /// Sets the label of the edge `u - v` (`u -> v` for a directed graph),
    /// replacing any previous label; a label of `0` removes the edge.
    ///
    /// For an undirected graph the arc `v -> u` is set as well. Returns
    /// `false` when an endpoint is out of range.
    pub fn set_edge(
        &mut self,
        u: usize,
        v: usize,
        label: u32,
    ) -> bool {
        if u >= self.n || v >= self.n {
            return false;
        }
        self.set_arc(u, v, label);
        if !self.directed && u != v {
            self.set_arc(v, u, label);
        }
        true
    }

    /// Sets the label of the arc `u -> v` only, in either kind of graph.
    fn set_arc(
        &mut self,
        u: usize,
        v: usize,
        label: u32,
    ) {
        let row = &mut self.matrix[u];
        let old = row[v];
        if old == label {
            return;
        }
        row[v] = label;
        update_list(&mut self.out[u], v, label);
        if self.directed {
            update_list(&mut self.incoming[v], u, label);
        } else {
            update_list(&mut self.out[v], u, label);
        }
    }

    /// The number of vertices.
    #[must_use]
    pub const fn num_nodes(&self) -> usize {
        self.n
    }

    /// Whether the arcs are directed.
    #[must_use]
    pub const fn is_directed(&self) -> bool {
        self.directed
    }

    /// The vertex colours.
    #[must_use]
    pub fn colours(&self) -> &[u32] {
        &self.colours
    }

    /// The label of the arc `u -> v`; `0` when there is none (or when an
    /// endpoint is out of range).
    #[must_use]
    pub fn edge_label(
        &self,
        u: usize,
        v: usize,
    ) -> u32 {
        self.matrix
            .get(u)
            .and_then(|row| row.get(v))
            .copied()
            .unwrap_or(0)
    }

    /// Whether the arc `u -> v` exists.
    #[must_use]
    pub fn has_edge(
        &self,
        u: usize,
        v: usize,
    ) -> bool {
        self.edge_label(u, v) != 0
    }

    /// The out-arcs of `v` as `(head, label)`, sorted by head.
    #[must_use]
    pub fn out_arcs(
        &self,
        v: usize,
    ) -> &[(usize, u32)] {
        self.out.get(v).map_or(&[], Vec::as_slice)
    }

    /// The in-arcs of `v` as `(tail, label)`; empty for an undirected graph.
    #[must_use]
    pub fn in_arcs(
        &self,
        v: usize,
    ) -> &[(usize, u32)] {
        self.incoming.get(v).map_or(&[], Vec::as_slice)
    }

    /// The number of vertices that are joined to `v` by an arc (in either
    /// direction for an undirected graph), self-loops included once.
    #[must_use]
    pub fn degree(
        &self,
        v: usize,
    ) -> usize {
        self.out.get(v).map_or(0, Vec::len)
    }

    /// The adjacency matrix, row by row.
    #[must_use]
    pub fn matrix(&self) -> &[Vec<u32>] {
        &self.matrix
    }

    /// The initial partition of the refinement: one cell per colour value, in
    /// increasing order of the value, the vertices of a cell in increasing
    /// order.
    fn initial_cells(&self) -> Vec<Vec<usize>> {
        let mut order: Vec<usize> = (0..self.n).collect();
        order.sort_by_key(|&v| (self.colours[v], v));
        let mut cells: Vec<Vec<usize>> = Vec::new();
        for v in order {
            match cells.last_mut() {
                | Some(cell) if self.colours[cell[0]] == self.colours[v] => cell.push(v),
                | _ => cells.push(vec![v]),
            }
        }
        cells
    }
}

/// Inserts, updates or removes `(head, label)` in a head-sorted arc list.
fn update_list(
    list: &mut Vec<(usize, u32)>,
    head: usize,
    label: u32,
) {
    match list.binary_search_by_key(&head, |&(h, _)| h) {
        | Ok(i) => {
            if label == 0 {
                list.remove(i);
            } else {
                list[i].1 = label;
            }
        },
        | Err(i) => {
            if label != 0 {
                list.insert(i, (head, label));
            }
        },
    }
}

// ---------------------------------------------------------------------
// Partition refinement
// ---------------------------------------------------------------------

/// The neighbour profile of one vertex: how many of its arcs go to each cell,
/// per label, as `(cell, label, count)` sorted by `(cell, label)`.
type Profile = Vec<(usize, u32, usize)>;

/// The refinement signature of a vertex: its cell and its two profiles.
type Signature = (usize, Profile, Profile);

/// Counts the arcs of `arcs` by `(cell of the other end, label)`.
fn profile(
    arcs: &[(usize, u32)],
    cell_of: &[usize],
) -> Profile {
    let mut keys: Vec<(usize, u32)> = arcs
        .iter()
        .map(|&(u, label)| (cell_of.get(u).copied().unwrap_or(0), label))
        .collect();
    keys.sort_unstable();
    let mut out: Profile = Vec::new();
    for (cell, label) in keys {
        match out.last_mut() {
            | Some(last) if (last.0, last.1) == (cell, label) => last.2 += 1,
            | _ => out.push((cell, label, 1)),
        }
    }
    out
}

/// The signature of `v` in the ordered partition given by `cell_of`.
fn signature(
    graph: &ColourGraph,
    v: usize,
    cell_of: &[usize],
) -> Signature {
    let cell = cell_of.get(v).copied().unwrap_or(0);
    let out = profile(graph.out_arcs(v), cell_of);
    let incoming = if graph.is_directed() {
        profile(graph.in_arcs(v), cell_of)
    } else {
        Vec::new()
    };
    (cell, out, incoming)
}

/// The cell index of every vertex.
fn cell_indices(
    n: usize,
    cells: &[Vec<usize>],
) -> Vec<usize> {
    let mut cell_of = vec![0; n];
    for (i, cell) in cells.iter().enumerate() {
        for &v in cell {
            if let Some(slot) = cell_of.get_mut(v) {
                *slot = i;
            }
        }
    }
    cell_of
}

/// Refines an ordered partition in place until it is equitable, returning the
/// number of rounds performed (at least one).
///
/// Each cell is replaced, in place, by its sub-cells ordered by signature, so
/// the result is a canonical function of the input partition and of the
/// coloured graph. Vertices of a cell that never split keep their relative
/// order; that order is irrelevant to the codes the search builds, because a
/// cell of more than one vertex is never part of a leaf.
fn refine_cells(
    graph: &ColourGraph,
    cells: &mut Vec<Vec<usize>>,
) -> usize {
    let mut rounds = 0;
    loop {
        rounds += 1;
        let cell_of = cell_indices(graph.num_nodes(), cells);
        let mut next: Vec<Vec<usize>> = Vec::with_capacity(cells.len());
        let mut split = false;
        for cell in cells.iter() {
            if cell.len() < 2 {
                next.push(cell.clone());
                continue;
            }
            let mut keyed: Vec<(Signature, usize)> = cell
                .iter()
                .map(|&v| (signature(graph, v, &cell_of), v))
                .collect();
            keyed.sort();
            let mut start = 0;
            while start < keyed.len() {
                let mut end = start + 1;
                while end < keyed.len() && keyed[end].0 == keyed[start].0 {
                    end += 1;
                }
                next.push(keyed[start..end].iter().map(|&(_, v)| v).collect());
                if end - start < cell.len() {
                    split = true;
                }
                start = end;
            }
        }
        *cells = next;
        if !split {
            return rounds;
        }
    }
}

/// One cell of a stable partition, with the data that identifies it.
///
/// Two cells are equal when they have the same initial colour, the same size
/// and the same neighbour counts into the cells of their partition in the same
/// *canonical cell order* — so a cell certificate is comparable between two
/// graphs.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CellCertificate {
    colour: u32,
    size: usize,
    out: Profile,
    incoming: Profile,
}

impl CellCertificate {
    /// The initial colour shared by the vertices of the cell.
    #[must_use]
    pub const fn colour(&self) -> u32 {
        self.colour
    }

    /// The number of vertices in the cell.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// How many out-neighbours a vertex of the cell has in each cell, per
    /// label: `(cell index, label, count)`.
    #[must_use]
    pub fn out_counts(&self) -> &[(usize, u32, usize)] {
        &self.out
    }

    /// The same for in-neighbours; empty for an undirected graph.
    #[must_use]
    pub fn in_counts(&self) -> &[(usize, u32, usize)] {
        &self.incoming
    }
}

/// An isomorphism-invariant fingerprint of a stable partition.
///
/// The certificate of a graph is equal to the certificate of every graph
/// isomorphic to it. The converse fails: all strongly regular pairs (the 4x4
/// rook graph and the Shrikhande graph, for instance) share a certificate,
/// which is exactly why the canonical labelling below is needed.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RefinementCertificate {
    n: usize,
    directed: bool,
    cells: Vec<CellCertificate>,
}

impl RefinementCertificate {
    /// The number of vertices.
    #[must_use]
    pub const fn num_nodes(&self) -> usize {
        self.n
    }

    /// Whether the graph is directed.
    #[must_use]
    pub const fn is_directed(&self) -> bool {
        self.directed
    }

    /// The number of cells, that is, of distinct stable colours.
    #[must_use]
    pub const fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Whether the colouring is discrete (every cell is a singleton).
    #[must_use]
    pub fn is_discrete(&self) -> bool {
        self.cells.iter().all(|cell| cell.size == 1)
    }

    /// The sizes of the cells, in canonical cell order.
    #[must_use]
    pub fn cell_sizes(&self) -> Vec<usize> {
        self.cells.iter().map(|cell| cell.size).collect()
    }

    /// The cells, in canonical cell order.
    #[must_use]
    pub fn cells(&self) -> &[CellCertificate] {
        &self.cells
    }
}

/// The stable (equitable) partition of a coloured graph.
///
/// Produced by [`refine`]. The cells are in a canonical order: by increasing
/// initial colour, and inside a colour class by increasing signature, so the
/// whole partition — not just its certificate — is independent of the vertex
/// numbering up to the renumbering itself.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refinement {
    n: usize,
    directed: bool,
    cells: Vec<Vec<usize>>,
    cell_of: Vec<usize>,
    rounds: usize,
    certificate: RefinementCertificate,
}

impl Refinement {
    /// The cells, in canonical order.
    #[must_use]
    pub fn cells(&self) -> &[Vec<usize>] {
        &self.cells
    }

    /// The cell index of every vertex.
    #[must_use]
    pub fn cell_of(&self) -> &[usize] {
        &self.cell_of
    }

    /// The number of cells.
    #[must_use]
    pub const fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Whether the graph the partition came from is directed.
    #[must_use]
    pub const fn is_directed(&self) -> bool {
        self.directed
    }

    /// The sizes of the cells, in canonical cell order.
    #[must_use]
    pub fn cell_sizes(&self) -> Vec<usize> {
        self.cells.iter().map(Vec::len).collect()
    }

    /// Whether every cell is a singleton.
    #[must_use]
    pub fn is_discrete(&self) -> bool {
        self.cells.iter().all(|cell| cell.len() == 1)
    }

    /// The number of refinement rounds performed.
    #[must_use]
    pub const fn rounds(&self) -> usize {
        self.rounds
    }

    /// The isomorphism-invariant certificate of the partition.
    #[must_use]
    pub const fn certificate(&self) -> &RefinementCertificate {
        &self.certificate
    }

    /// Checks the result against a graph: the cells partition the vertices,
    /// each cell has one colour, and each cell is equitable — every vertex of
    /// a cell has the same neighbour counts (per label) into every cell.
    ///
    /// Returns `false` when `graph` is not the graph the partition came from
    /// (a different number of vertices, or a vertex out of range).
    #[must_use]
    pub fn verify(
        &self,
        graph: &ColourGraph,
    ) -> bool {
        if self.n != graph.num_nodes() {
            return false;
        }
        let mut seen = vec![false; self.n];
        for cell in &self.cells {
            if cell.is_empty() {
                return false;
            }
            let colour = graph.colours[cell[0]];
            let mut first: Option<Signature> = None;
            for &v in cell {
                if v >= self.n || seen[v] || graph.colours[v] != colour {
                    return false;
                }
                seen[v] = true;
                let sig = signature(graph, v, &self.cell_of);
                match &first {
                    | None => first = Some(sig),
                    | Some(want) if *want != sig => return false,
                    | Some(_) => {},
                }
            }
        }
        seen.iter().all(|&s| s)
    }
}

/// Refines a coloured graph to its coarsest equitable partition: the
/// 1-dimensional Weisfeiler–Leman colouring.
///
/// The result depends on the isomorphism class only, so two graphs whose
/// certificates differ are certainly not isomorphic. Equal certificates prove
/// nothing beyond the first-order neighbourhood structure of the graph.
#[must_use]
pub fn refine(graph: &ColourGraph) -> Refinement {
    let n = graph.num_nodes();
    let mut cells = graph.initial_cells();
    let rounds = refine_cells(graph, &mut cells);
    let cell_of = cell_indices(n, &cells);
    let certificate = RefinementCertificate {
        n,
        directed: graph.is_directed(),
        cells: cells
            .iter()
            .map(|cell| {
                let v = cell[0];
                CellCertificate {
                    colour: graph.colours[v],
                    size: cell.len(),
                    out: profile(graph.out_arcs(v), &cell_of),
                    incoming: if graph.is_directed() {
                        profile(graph.in_arcs(v), &cell_of)
                    } else {
                        Vec::new()
                    },
                }
            })
            .collect(),
    };
    Refinement {
        n,
        directed: graph.is_directed(),
        cells,
        cell_of,
        rounds,
        certificate,
    }
}

/// Whether the 1-dimensional Weisfeiler–Leman colourings of two graphs agree.
///
/// A necessary condition for isomorphism, and an exact test on trees; on
/// random and on most structured inputs it decides as well, but it is fooled
/// by every pair of strongly regular graphs with the same parameters (the 4x4
/// rook graph and the Shrikhande graph, for instance), by which time only
/// [`canonical_labelling`] can tell them apart.
#[must_use]
pub fn wl_equivalent(
    a: &ColourGraph,
    b: &ColourGraph,
) -> bool {
    a.num_nodes() == b.num_nodes()
        && a.is_directed() == b.is_directed()
        && refine(a).certificate == refine(b).certificate
}

// ---------------------------------------------------------------------
// Permutations and their groups
// ---------------------------------------------------------------------

/// A permutation of `0..n`, stored as the list of images.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Permutation {
    image: Vec<usize>,
}

impl Permutation {
    /// The identity of `n` points.
    #[must_use]
    pub fn identity(n: usize) -> Self {
        Self {
            image: (0..n).collect(),
        }
    }

    /// A permutation from its image list.
    ///
    /// Returns `None` unless `image` is a permutation of `0..image.len()`.
    #[must_use]
    pub fn from_images(image: Vec<usize>) -> Option<Self> {
        let n = image.len();
        let mut seen = vec![false; n];
        for &x in &image {
            if x >= n || seen[x] {
                return None;
            }
            seen[x] = true;
        }
        Some(Self { image })
    }

    /// A permutation known to be valid, from its image list.
    const fn exact(image: Vec<usize>) -> Self {
        Self { image }
    }

    /// The image list: entry `i` is the image of `i`.
    #[must_use]
    pub fn image(&self) -> &[usize] {
        &self.image
    }

    /// The number of points moved by the permutation.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.image.len()
    }

    /// Whether the permutation acts on no point at all.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.image.is_empty()
    }

    /// The image of `v`, or `None` when `v` is not a point.
    #[must_use]
    pub fn apply(
        &self,
        v: usize,
    ) -> Option<usize> {
        self.image.get(v).copied()
    }

    /// `self ∘ other` (apply `other` first), or `None` on a length mismatch.
    #[must_use]
    pub fn compose(
        &self,
        other: &Self,
    ) -> Option<Self> {
        if self.len() != other.len() {
            return None;
        }
        Some(Self::exact(
            other.image.iter().map(|&x| self.image[x]).collect(),
        ))
    }

    /// The inverse permutation.
    #[must_use]
    pub fn inverse(&self) -> Self {
        let mut image = vec![0; self.image.len()];
        for (i, &x) in self.image.iter().enumerate() {
            image[x] = i;
        }
        Self::exact(image)
    }

    /// Whether every point is fixed.
    #[must_use]
    pub fn is_identity(&self) -> bool {
        is_identity(&self.image)
    }

    /// The non-trivial cycles, each starting at its smallest point.
    #[must_use]
    pub fn cycles(&self) -> Vec<Vec<usize>> {
        let mut seen = vec![false; self.image.len()];
        let mut out = Vec::new();
        for start in 0..self.image.len() {
            if seen[start] || self.image[start] == start {
                continue;
            }
            let mut cycle = vec![start];
            seen[start] = true;
            let mut v = self.image[start];
            while v != start {
                seen[v] = true;
                cycle.push(v);
                v = self.image[v];
            }
            out.push(cycle);
        }
        out
    }
}

/// Whether `p` fixes every point.
fn is_identity(p: &[usize]) -> bool {
    p.iter().enumerate().all(|(i, &x)| i == x)
}

/// `a ∘ b`: apply `b` first.
fn compose(
    a: &[usize],
    b: &[usize],
) -> Vec<usize> {
    b.iter().map(|&x| a[x]).collect()
}

/// The inverse of a permutation given by its images.
fn inverse(p: &[usize]) -> Vec<usize> {
    let mut inv = vec![0; p.len()];
    for (i, &x) in p.iter().enumerate() {
        inv[x] = i;
    }
    inv
}

/// A stabiliser chain (base plus strong generating set) of a permutation
/// group, maintained incrementally with the Schreier–Sims algorithm.
///
/// Only what the canonical labelling needs is here: [`insert`] a generator
/// (redundant ones are sifted away), [`contains`] a permutation, and the
/// [`order`] of the group as a product of basic orbit lengths.
///
/// [`insert`]: Self::insert
/// [`contains`]: Self::contains
/// [`order`]: Self::order
#[derive(Clone, Debug)]
struct StabiliserChain {
    n: usize,
    base: Vec<usize>,
    /// `gens[level]`: generators fixing `base[..level]` pointwise.
    gens: Vec<Vec<Vec<usize>>>,
    /// `trans[level][point]`: a group element taking `base[level]` to
    /// `point`.
    trans: Vec<Vec<Option<Vec<usize>>>>,
    /// `orbit[level]`: the orbit of `base[level]` under `gens[level]`.
    orbit: Vec<Vec<usize>>,
}

impl StabiliserChain {
    /// The trivial group on `n` points.
    const fn new(n: usize) -> Self {
        Self {
            n,
            base: Vec::new(),
            gens: Vec::new(),
            trans: Vec::new(),
            orbit: Vec::new(),
        }
    }

    /// Recomputes the orbit and the transversal of one level.
    fn recompute_orbit(
        &mut self,
        level: usize,
    ) {
        let b = self.base[level];
        let mut trans: Vec<Option<Vec<usize>>> = vec![None; self.n];
        trans[b] = Some((0..self.n).collect());
        let mut orbit = vec![b];
        let mut head = 0;
        while head < orbit.len() {
            let point = orbit[head];
            head += 1;
            let u = trans[point]
                .clone()
                .unwrap_or_else(|| (0..self.n).collect());
            for s in &self.gens[level] {
                let image = s[point];
                if trans[image].is_none() {
                    trans[image] = Some(compose(s, &u));
                    orbit.push(image);
                }
            }
        }
        self.trans[level] = trans;
        self.orbit[level] = orbit;
    }

    /// Sifts `g` through the levels from `from` on, returning the residue and
    /// the level at which it got stuck (the depth when it went through).
    fn sift(
        &self,
        mut g: Vec<usize>,
        from: usize,
    ) -> (Vec<usize>, usize) {
        for level in from..self.base.len() {
            let point = g[self.base[level]];
            match &self.trans[level][point] {
                | Some(u) => g = compose(&inverse(u), &g),
                | None => return (g, level),
            }
        }
        (g, self.base.len())
    }

    /// Adds the non-identity residue `r` as a strong generator of the levels
    /// `from..=to`, creating a level when `to` is the depth.
    fn add_strong(
        &mut self,
        r: &[usize],
        from: usize,
        to: usize,
    ) {
        if to == self.base.len() {
            let b = (0..self.n).find(|&x| r[x] != x).unwrap_or(0);
            self.base.push(b);
            self.gens.push(Vec::new());
            self.trans.push(Vec::new());
            self.orbit.push(Vec::new());
        }
        for level in from..=to {
            self.gens[level].push(r.to_vec());
            self.recompute_orbit(level);
        }
    }

    /// Verifies the Schreier generators from level `start - 1` upwards and
    /// repairs the chain.
    fn complete(
        &mut self,
        start: usize,
    ) {
        let mut i = start;
        while i > 0 {
            let level = i - 1;
            let mut added = None;
            'outer: for &point in &self.orbit[level].clone() {
                for s in &self.gens[level].clone() {
                    let image = s[point];
                    let (Some(up), Some(ui)) =
                        (&self.trans[level][point], &self.trans[level][image])
                    else {
                        continue;
                    };
                    let h = compose(&inverse(ui), &compose(s, up));
                    let (r, j) = self.sift(h, level + 1);
                    if !is_identity(&r) {
                        self.add_strong(&r, level + 1, j);
                        added = Some(j + 1);
                        break 'outer;
                    }
                }
            }
            match added {
                | Some(next) => i = next,
                | None => i -= 1,
            }
        }
    }

    /// Adds a generator to the group.
    fn insert(
        &mut self,
        g: &[usize],
    ) {
        let (r, j) = self.sift(g.to_vec(), 0);
        if !is_identity(&r) {
            self.add_strong(&r, 0, j);
            let depth = self.base.len();
            self.complete(depth);
        }
    }

    /// Whether `g` lies in the group.
    fn contains(
        &self,
        g: &[usize],
    ) -> bool {
        g.len() == self.n && {
            let (r, level) = self.sift(g.to_vec(), 0);
            level == self.base.len() && is_identity(&r)
        }
    }

    /// The order of the group.
    fn order(&self) -> BigInt {
        self.orbit.iter().fold(BigInt::from(1), |acc, orbit| {
            acc * BigInt::from(orbit.len())
        })
    }
}

/// The orbits of the group generated by `gens`, each sorted, ordered by their
/// smallest point.
fn orbits_of(
    n: usize,
    gens: &[Permutation],
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
            let point = orbit[head];
            head += 1;
            for g in gens {
                let image = g.image[point];
                if !seen[image] {
                    seen[image] = true;
                    orbit.push(image);
                }
            }
        }
        orbit.sort_unstable();
        out.push(orbit);
    }
    out
}

/// A disjoint-set forest used for the orbit test of the automorphism pruning.
struct Orbits {
    parent: Vec<usize>,
}

impl Orbits {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(
        &mut self,
        mut v: usize,
    ) -> usize {
        while self.parent[v] != v {
            let grand = self.parent[self.parent[v]];
            self.parent[v] = grand;
            v = grand;
        }
        v
    }

    fn union(
        &mut self,
        a: usize,
        b: usize,
    ) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

// ---------------------------------------------------------------------
// Individualisation–refinement
// ---------------------------------------------------------------------

/// The bound on the individualisation–refinement search.
///
/// The budget counts *search nodes*, not leaves, so it bounds the work
/// regardless of the size of the automorphism group. When it runs out the
/// search stops and the result is marked inexact (see the module docs for the
/// fallback contract).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SearchBudget {
    max_nodes: usize,
}

impl SearchBudget {
    /// A budget of `max_nodes` search nodes.
    #[must_use]
    pub const fn nodes(max_nodes: usize) -> Self {
        Self { max_nodes }
    }

    /// The default budget for a graph on `n` vertices: `2_000 + 20 n` nodes.
    ///
    /// The constant term covers small graphs with large automorphism groups
    /// (where the pruning makes the tree quadratic), the linear term covers
    /// graphs whose colouring stays coarse.
    #[must_use]
    pub const fn for_size(n: usize) -> Self {
        Self {
            max_nodes: 2_000usize.saturating_add(20usize.saturating_mul(n)),
        }
    }

    /// The maximum number of search nodes.
    #[must_use]
    pub const fn max_nodes(&self) -> usize {
        self.max_nodes
    }
}

/// A canonical labelling of a coloured graph, with its automorphism group.
///
/// Produced by [`canonical_labelling`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CanonicalLabelling {
    n: usize,
    order: Vec<usize>,
    position: Vec<usize>,
    code: Vec<u32>,
    generators: Vec<Permutation>,
    orbits: Vec<Vec<usize>>,
    group_order: BigInt,
    exact: bool,
    nodes: usize,
}

impl CanonicalLabelling {
    /// The canonical vertex order: `order()[i]` is the original vertex that
    /// becomes vertex `i` of the canonical form.
    ///
    /// When the automorphism group is non-trivial the minimum code is reached
    /// by several leaves; any of their orders is returned, and all of them
    /// give the same [`code`](Self::code).
    #[must_use]
    pub fn order(&self) -> &[usize] {
        &self.order
    }

    /// The canonical position of an original vertex, or `None` when out of
    /// range.
    #[must_use]
    pub fn position(
        &self,
        v: usize,
    ) -> Option<usize> {
        self.position.get(v).copied()
    }

    /// The canonical code: the `n` vertex colours in canonical order,
    /// followed by the `n * n` adjacency matrix rows in canonical order.
    ///
    /// Two graphs are isomorphic if and only if their canonical codes are
    /// equal, provided both searches were exact.
    #[must_use]
    pub fn code(&self) -> &[u32] {
        &self.code
    }

    /// The canonical adjacency matrix, rows in canonical order.
    #[must_use]
    pub fn matrix(&self) -> Vec<Vec<u32>> {
        let n = self.n;
        (0..n)
            .map(|i| (0..n).map(|j| self.code[n + i * n + j]).collect())
            .collect()
    }

    /// The generators of the automorphism group found by the search.
    ///
    /// They generate the full group when [`is_exact`](Self::is_exact) holds,
    /// and a subgroup of it otherwise.
    #[must_use]
    pub fn generators(&self) -> &[Permutation] {
        &self.generators
    }

    /// The orbits of the automorphism group, each sorted, ordered by their
    /// smallest vertex.
    #[must_use]
    pub fn orbits(&self) -> &[Vec<usize>] {
        &self.orbits
    }

    /// The order of the automorphism group: `n!` for the complete graph, `2n`
    /// for the cycle, `120` for the Petersen graph.
    #[must_use]
    pub const fn group_order(&self) -> &BigInt {
        &self.group_order
    }

    /// The number of vertices.
    #[must_use]
    pub const fn num_nodes(&self) -> usize {
        self.n
    }

    /// Whether the search finished within its budget, that is, whether the
    /// code is guaranteed to be the canonical one and the group the full
    /// automorphism group.
    #[must_use]
    pub const fn is_exact(&self) -> bool {
        self.exact
    }

    /// The number of search nodes visited.
    #[must_use]
    pub const fn num_nodes_visited(&self) -> usize {
        self.nodes
    }
}

/// The mutable state of one search.
struct Search<'a> {
    graph: &'a ColourGraph,
    n: usize,
    best: Option<Vec<u32>>,
    best_order: Vec<usize>,
    candidate: Vec<u32>,
    chain: StabiliserChain,
    generators: Vec<Permutation>,
    nodes: usize,
    max_nodes: usize,
    aborted: bool,
}

impl Search<'_> {
    /// Refines and branches, or records a leaf.
    fn run(
        &mut self,
        mut cells: Vec<Vec<usize>>,
        prefix: &mut Vec<usize>,
    ) {
        self.nodes += 1;
        if self.nodes > self.max_nodes {
            self.aborted = true;
            return;
        }
        refine_cells(self.graph, &mut cells);
        let target = cells.iter().position(|cell| cell.len() > 1);
        let Some(target) = target else {
            let order: Vec<usize> = cells.iter().flatten().copied().collect();
            self.leaf(&order);
            return;
        };
        let targets = cells[target].clone();
        let mut branched: Vec<usize> = Vec::new();
        // The orbits of the found automorphisms that fix the prefix are only
        // invalidated by a new generator, so one disjoint-set forest serves
        // every vertex of the cell until then.
        let mut orbit_cache: Option<(usize, Orbits)> = None;
        for v in targets {
            if self.aborted {
                return;
            }
            if !branched.is_empty() {
                if orbit_cache
                    .as_ref()
                    .is_none_or(|(count, _)| *count != self.generators.len())
                {
                    orbit_cache = Some((self.generators.len(), self.orbits_fixing(prefix)));
                }
                if let Some((_, orbits)) = &mut orbit_cache {
                    let root = orbits.find(v);
                    if branched.iter().any(|&u| orbits.find(u) == root) {
                        continue;
                    }
                }
            }
            let mut child: Vec<Vec<usize>> = Vec::with_capacity(cells.len() + 1);
            for (i, cell) in cells.iter().enumerate() {
                if i == target {
                    child.push(vec![v]);
                    child.push(cell.iter().copied().filter(|&u| u != v).collect());
                } else {
                    child.push(cell.clone());
                }
            }
            prefix.push(v);
            self.run(child, prefix);
            prefix.pop();
            branched.push(v);
        }
    }

    /// The orbits of the found automorphisms that fix `prefix` pointwise, as a
    /// disjoint-set forest over all vertices.
    fn orbits_fixing(
        &self,
        prefix: &[usize],
    ) -> Orbits {
        let mut orbits = Orbits::new(self.n);
        for g in &self.generators {
            if prefix.iter().all(|&p| g.image[p] == p) {
                for x in 0..self.n {
                    orbits.union(x, g.image[x]);
                }
            }
        }
        orbits
    }

    /// Builds the code of a leaf order and compares it with the best code.
    fn leaf(
        &mut self,
        order: &[usize],
    ) {
        self.candidate.clear();
        self.candidate.reserve(self.n * (self.n + 1));
        for &v in order {
            self.candidate.push(self.graph.colours[v]);
        }
        for &u in order {
            for &v in order {
                self.candidate.push(self.graph.matrix[u][v]);
            }
        }
        let ordering = match &self.best {
            | Some(best) => self.candidate.cmp(best),
            | None => Ordering::Less,
        };
        match ordering {
            | Ordering::Less => {
                self.best = Some(self.candidate.clone());
                self.best_order.clear();
                self.best_order.extend_from_slice(order);
            },
            | Ordering::Equal => {
                let mut sigma = vec![0; self.n];
                for (i, &v) in order.iter().enumerate() {
                    sigma[self.best_order[i]] = v;
                }
                if !is_identity(&sigma) && !self.chain.contains(&sigma) {
                    self.chain.insert(&sigma);
                    self.generators.push(Permutation::exact(sigma));
                }
            },
            | Ordering::Greater => {},
        }
    }
}

/// Searches for the canonical labelling of `graph` and its automorphism
/// group.
///
/// Returns the labelling with the smallest code found within `budget` nodes,
/// the generators and orbits of the automorphism group found on the way, and
/// whether the search was exact (see the module docs for the fallback).
///
/// # Complexity
///
/// A node costs one equitable refinement, `O(n + m)` (times the number of
/// rounds, at most `n`), and a leaf costs one code, `O(n^2)`. The number of
/// nodes is that of the individualisation tree after automorphism pruning;
/// for a graph whose colouring is discrete it is one.
#[must_use]
pub fn canonical_labelling(
    graph: &ColourGraph,
    budget: SearchBudget,
) -> CanonicalLabelling {
    let n = graph.num_nodes();
    let mut search = Search {
        graph,
        n,
        best: None,
        best_order: Vec::new(),
        candidate: Vec::new(),
        chain: StabiliserChain::new(n),
        generators: Vec::new(),
        nodes: 0,
        max_nodes: budget.max_nodes,
        aborted: false,
    };
    let mut prefix = Vec::new();
    search.run(graph.initial_cells(), &mut prefix);
    // A leaf always has `n` vertices in its order; the fallback keeps the
    // function total when the budget is exhausted at the root.
    let (order, code, exact) = if let Some(best) = search.best.take()
        && search.best_order.len() == n
    {
        (
            std::mem::take(&mut search.best_order),
            best,
            !search.aborted,
        )
    } else {
        let order: Vec<usize> = (0..n).collect();
        let code = identity_code(graph, &order);
        (order, code, false)
    };
    let mut position = vec![0; n];
    for (i, &v) in order.iter().enumerate() {
        position[v] = i;
    }
    let orbits = orbits_of(n, &search.generators);
    let group_order = search.chain.order();
    CanonicalLabelling {
        n,
        order,
        position,
        code,
        generators: search.generators,
        orbits,
        group_order,
        exact,
        nodes: search.nodes,
    }
}

/// The code of `graph` in the given vertex order.
fn identity_code(
    graph: &ColourGraph,
    order: &[usize],
) -> Vec<u32> {
    let n = graph.num_nodes();
    let mut code = Vec::with_capacity(n * (n + 1));
    for &v in order {
        code.push(graph.colours[v]);
    }
    for &u in order {
        for &v in order {
            code.push(graph.matrix[u][v]);
        }
    }
    code
}

/// Whether two coloured graphs are isomorphic.
///
/// Returns `None` when the budget ran out before either search finished: the
/// graphs then remain undecided. A `Some(false)` answer may also be reached
/// from the refinement certificate alone, which is cheaper.
#[must_use]
pub fn isomorphic(
    a: &ColourGraph,
    b: &ColourGraph,
    budget: SearchBudget,
) -> Option<bool> {
    if a.num_nodes() != b.num_nodes() || a.is_directed() != b.is_directed() {
        return Some(false);
    }
    if !wl_equivalent(a, b) {
        return Some(false);
    }
    let (fa, fb) = (
        canonical_labelling(a, budget),
        canonical_labelling(b, budget),
    );
    (fa.exact && fb.exact).then(|| fa.code == fb.code)
}
