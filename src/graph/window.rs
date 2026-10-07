//! Tree windows: mutable linked-list trees carved out of the DAG for cheap
//! local optimisation.
//!
//! Equality saturation pays for every intermediate term with a permanent
//! e-node. For the bread-and-butter work of algebra — reassociating sums,
//! collecting like terms, folding constants — that price is absurd: the
//! intermediates are worthless and associativity/commutativity alone blow
//! the graph up exponentially.
//!
//! A [`TreeWindow`] is the alternative. Starting at a concrete node it
//! copies the *unshared* region below it into a slab of doubly linked cells
//! (parent, first/last child, previous/next sibling). Shared subterms and
//! leaves stay behind as opaque *atoms* that point back into the graph.
//! Inside the window, rewriting is destructive and O(1) per splice, search
//! can be greedy or a beam, and nothing touches the graph until
//! [`TreeWindow::commit`] interns the winner — a single new term that the
//! caller unions with the original.
//!
//! So the DAG owns sharing and equivalence, the window owns local search,
//! and carve/commit is the bridge between them.

use std::collections::HashSet;

use super::id::NodeId;
use super::id::OpId;
use super::id::to_u32;
use super::number::Number;
use super::op::OpFlags;
use super::pattern::Pat;
use super::rule::Rewrite;
use super::store::Graph;

/// Index of a cell inside a [`TreeWindow`].
pub type CellId = u32;

const NIL: CellId = u32::MAX;

#[derive(Clone, Debug)]
struct Cell {
    op: OpId,
    /// Graph node this cell stands for, or `NONE` for an interior cell.
    atom: NodeId,
    parent: CellId,
    first: CellId,
    last: CellId,
    prev: CellId,
    next: CellId,
    arity: u32,
    live: bool,
}

/// A custom transformation over a whole window (constant folding, like-term
/// collection, ...). Must only perform identities.
pub trait WindowPass: Send + Sync {
    /// Runs the pass; returns `true` if the window changed.
    fn run(
        &self,
        graph: &mut Graph,
        window: &mut TreeWindow,
    ) -> bool;
}

/// A mutable tree view of the unshared region below a DAG node.
#[derive(Clone, Debug)]
pub struct TreeWindow {
    cells: Vec<Cell>,
    free: Vec<CellId>,
    root: CellId,
    live: usize,
}

impl TreeWindow {
    /// Carves the window rooted at `root`.
    ///
    /// A child is copied into the window when it is an operator application
    /// that is either unshared (its class has a single parent) or
    /// [`TRANSPARENT`](OpFlags::TRANSPARENT), and the window holds fewer
    /// than `max_cells` cells; otherwise it becomes an atom.
    #[must_use]
    pub fn carve(
        graph: &Graph,
        root: NodeId,
        max_cells: usize,
    ) -> Self {
        let mut window = Self {
            cells: Vec::new(),
            free: Vec::new(),
            root: NIL,
            live: 0,
        };
        // (node, parent cell); children are pushed in reverse to keep order.
        let mut stack = vec![(root, NIL)];
        while let Some((node, parent)) = stack.pop() {
            let children = graph.children(node);
            let transparent = graph.ops().get(graph.op(node)).flags.has(OpFlags::TRANSPARENT);
            let inline = !children.is_empty()
                && (parent == NIL
                    || (window.live < max_cells
                        && (transparent || graph.parents(graph.find(node)).len() <= 1)));
            let cell = if inline {
                window.alloc(graph.op(node), NodeId::NONE)
            } else {
                window.new_atom(graph, node)
            };
            if parent == NIL {
                window.root = cell;
            } else {
                window.append_child(parent, cell);
            }
            if inline {
                stack.extend(children.iter().rev().map(|&c| (c, cell)));
            }
        }
        // The stack visits children in order, but siblings' subtrees finish
        // before later siblings are attached, so order is preserved.
        window
    }

    fn alloc(
        &mut self,
        op: OpId,
        atom: NodeId,
    ) -> CellId {
        let cell = Cell {
            op,
            atom,
            parent: NIL,
            first: NIL,
            last: NIL,
            prev: NIL,
            next: NIL,
            arity: 0,
            live: true,
        };
        self.live = self.live.saturating_add(1);
        if let Some(id) = self.free.pop() {
            if let Some(slot) = self.cells.get_mut(id as usize) {
                *slot = cell;
            }
            id
        } else {
            self.cells.push(cell);
            to_u32(self.cells.len().saturating_sub(1))
        }
    }

    #[inline]
    fn cell(
        &self,
        id: CellId,
    ) -> &Cell {
        match self.cells.get(id as usize) {
            | Some(c) => c,
            | None => bad_cell(id),
        }
    }

    #[inline]
    fn cell_mut(
        &mut self,
        id: CellId,
    ) -> &mut Cell {
        match self.cells.get_mut(id as usize) {
            | Some(c) => c,
            | None => bad_cell(id),
        }
    }

    // ------------------------------------------------------------------
    // Reading
    // ------------------------------------------------------------------

    /// The root cell.
    #[must_use]
    pub const fn root(&self) -> CellId {
        self.root
    }

    /// Number of live cells.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.live
    }

    /// Whether the window holds no cells.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.live == 0
    }

    /// Whether `cell` is still part of the window. Cells freed by an edit
    /// stay addressable until their slot is reused.
    #[must_use]
    pub fn is_live(
        &self,
        cell: CellId,
    ) -> bool {
        self.cells.get(cell as usize).is_some_and(|c| c.live)
    }

    /// Operator of a cell. For an atom this is the operator of the graph
    /// node it stands for.
    #[must_use]
    pub fn op(
        &self,
        cell: CellId,
    ) -> OpId {
        self.cell(cell).op
    }

    /// The graph node an atom stands for; `None` for interior cells.
    #[must_use]
    pub fn atom(
        &self,
        cell: CellId,
    ) -> Option<NodeId> {
        let atom = self.cell(cell).atom;
        (!atom.is_none()).then_some(atom)
    }

    /// Parent of a cell, `None` at the root.
    #[must_use]
    pub fn parent(
        &self,
        cell: CellId,
    ) -> Option<CellId> {
        let parent = self.cell(cell).parent;
        (parent != NIL).then_some(parent)
    }

    /// Number of children.
    #[must_use]
    pub fn arity(
        &self,
        cell: CellId,
    ) -> usize {
        self.cell(cell).arity as usize
    }

    /// Children of a cell, in order.
    pub fn children(
        &self,
        cell: CellId,
    ) -> impl Iterator<Item = CellId> + '_ {
        let mut cur = self.cell(cell).first;
        std::iter::from_fn(move || {
            if cur == NIL {
                return None;
            }
            let id = cur;
            cur = self.cell(id).next;
            Some(id)
        })
    }

    /// The literal number an atom stands for.
    #[must_use]
    pub fn number<'g>(
        &self,
        graph: &'g Graph,
        cell: CellId,
    ) -> Option<&'g Number> {
        self.atom(cell).and_then(|n| graph.number_of(n))
    }

    /// All live cells in post-order (children before parents).
    #[must_use]
    pub fn post_order(&self) -> Vec<CellId> {
        let mut out = Vec::with_capacity(self.live);
        if self.root == NIL {
            return out;
        }
        let mut stack = vec![(self.root, false)];
        while let Some((cell, expanded)) = stack.pop() {
            if expanded {
                out.push(cell);
            } else {
                stack.push((cell, true));
                let mut child = self.cell(cell).last;
                while child != NIL {
                    stack.push((child, false));
                    child = self.cell(child).prev;
                }
            }
        }
        out
    }

    /// Sum of operator costs of interior cells plus one per atom.
    #[must_use]
    pub fn cost(
        &self,
        graph: &Graph,
    ) -> u64 {
        self.post_order()
            .into_iter()
            .map(|c| {
                if self.atom(c).is_some() {
                    1
                } else {
                    u64::from(graph.ops().get(self.op(c)).cost)
                }
            })
            .sum()
    }

    /// Order-insensitive (for commutative operators) structural hash of the
    /// subtree at `cell`, modulo the equalities the graph knows for atoms.
    #[must_use]
    pub fn fingerprint(
        &self,
        graph: &Graph,
        cell: CellId,
    ) -> u64 {
        const K: u64 = 0x9e37_79b9_7f4a_7c15;
        if let Some(atom) = self.atom(cell) {
            return u64::from(graph.find(atom).raw()).wrapping_mul(K) ^ 0x5555;
        }
        let op = self.op(cell);
        let mut hashes: Vec<u64> = self
            .children(cell)
            .map(|c| self.fingerprint(graph, c))
            .collect();
        if graph.ops().get(op).flags.has(OpFlags::COMMUTATIVE) {
            hashes.sort_unstable();
        }
        hashes.into_iter().fold(
            u64::from(op.raw()).wrapping_add(1).wrapping_mul(K),
            |h, c| (h.rotate_left(7) ^ c).wrapping_mul(K),
        )
    }

    /// Structural equality of two subtrees (commutative children compared
    /// as multisets, atoms modulo graph equalities).
    #[must_use]
    pub fn same(
        &self,
        graph: &Graph,
        a: CellId,
        b: CellId,
    ) -> bool {
        if a == b {
            return true;
        }
        match (self.atom(a), self.atom(b)) {
            | (Some(x), Some(y)) => return graph.same(x, y),
            | (None, None) => {},
            | _ => return false,
        }
        if self.op(a) != self.op(b) || self.arity(a) != self.arity(b) {
            return false;
        }
        let mut ca: Vec<CellId> = self.children(a).collect();
        let mut cb: Vec<CellId> = self.children(b).collect();
        if graph.ops().get(self.op(a)).flags.has(OpFlags::COMMUTATIVE) {
            ca.sort_by_key(|&c| self.fingerprint(graph, c));
            cb.sort_by_key(|&c| self.fingerprint(graph, c));
        }
        ca.iter().zip(&cb).all(|(&x, &y)| self.same(graph, x, y))
    }

    // ------------------------------------------------------------------
    // Editing
    // ------------------------------------------------------------------

    /// Creates a detached atom standing for `node`.
    pub fn new_atom(
        &mut self,
        graph: &Graph,
        node: NodeId,
    ) -> CellId {
        self.alloc(graph.op(node), node)
    }

    /// Creates a detached interior cell with the given (detached) children.
    pub fn new_node(
        &mut self,
        op: OpId,
        children: &[CellId],
    ) -> CellId {
        let cell = self.alloc(op, NodeId::NONE);
        for &child in children {
            self.append_child(cell, child);
        }
        cell
    }

    /// Appends the detached cell `child` as the last child of `parent`.
    pub fn append_child(
        &mut self,
        parent: CellId,
        child: CellId,
    ) {
        let last = self.cell(parent).last;
        {
            let c = self.cell_mut(child);
            c.parent = parent;
            c.prev = last;
            c.next = NIL;
        }
        if last == NIL {
            self.cell_mut(parent).first = child;
        } else {
            self.cell_mut(last).next = child;
        }
        let p = self.cell_mut(parent);
        p.last = child;
        p.arity = p.arity.saturating_add(1);
    }

    /// Unlinks `cell` from its parent, leaving it a detached subtree.
    pub fn detach(
        &mut self,
        cell: CellId,
    ) {
        let Cell {
            parent, prev, next, ..
        } = *self.cell(cell);
        if parent == NIL {
            if self.root == cell {
                self.root = NIL;
            }
            return;
        }
        if prev == NIL {
            self.cell_mut(parent).first = next;
        } else {
            self.cell_mut(prev).next = next;
        }
        if next == NIL {
            self.cell_mut(parent).last = prev;
        } else {
            self.cell_mut(next).prev = prev;
        }
        let p = self.cell_mut(parent);
        p.arity = p.arity.saturating_sub(1);
        let c = self.cell_mut(cell);
        c.parent = NIL;
        c.prev = NIL;
        c.next = NIL;
    }

    /// Frees the detached subtree at `cell`.
    pub fn discard(
        &mut self,
        cell: CellId,
    ) {
        let mut stack = vec![cell];
        while let Some(id) = stack.pop() {
            stack.extend(self.children(id));
            self.cell_mut(id).live = false;
            self.live = self.live.saturating_sub(1);
            self.free.push(id);
        }
    }

    /// Puts the detached subtree `new` where `old` is and frees `old`.
    pub fn replace(
        &mut self,
        old: CellId,
        new: CellId,
    ) {
        if old == new {
            return;
        }
        let Cell {
            parent, prev, next, ..
        } = *self.cell(old);
        {
            let n = self.cell_mut(new);
            n.parent = parent;
            n.prev = prev;
            n.next = next;
        }
        if parent == NIL {
            if self.root == old {
                self.root = new;
            }
        } else {
            if prev == NIL {
                self.cell_mut(parent).first = new;
            } else {
                self.cell_mut(prev).next = new;
            }
            if next == NIL {
                self.cell_mut(parent).last = new;
            } else {
                self.cell_mut(next).prev = new;
            }
        }
        let o = self.cell_mut(old);
        o.parent = NIL;
        o.prev = NIL;
        o.next = NIL;
        self.discard(old);
    }

    /// Replaces `child` by its own children, in place: the O(1) list splice
    /// that flattens `f(.., f(a, b), ..)` into `f(.., a, b, ..)`.
    pub fn splice_up(
        &mut self,
        child: CellId,
    ) {
        let Cell {
            parent,
            prev,
            next,
            first,
            last,
            arity,
            ..
        } = *self.cell(child);
        if parent == NIL {
            return;
        }
        if first == NIL {
            self.detach(child);
            self.discard(child);
            return;
        }
        let mut cur = first;
        while cur != NIL {
            self.cell_mut(cur).parent = parent;
            cur = self.cell(cur).next;
        }
        self.cell_mut(first).prev = prev;
        self.cell_mut(last).next = next;
        if prev == NIL {
            self.cell_mut(parent).first = first;
        } else {
            self.cell_mut(prev).next = first;
        }
        if next == NIL {
            self.cell_mut(parent).last = last;
        } else {
            self.cell_mut(next).prev = last;
        }
        let p = self.cell_mut(parent);
        p.arity = p.arity.saturating_add(arity).saturating_sub(1);
        let c = self.cell_mut(child);
        c.first = NIL;
        c.last = NIL;
        c.parent = NIL;
        c.prev = NIL;
        c.next = NIL;
        c.live = false;
        self.live = self.live.saturating_sub(1);
        self.free.push(child);
    }

    /// Deep-copies the subtree at `cell` into a new detached subtree.
    pub fn clone_subtree(
        &mut self,
        cell: CellId,
    ) -> CellId {
        let Cell { op, atom, .. } = *self.cell(cell);
        let copy = self.alloc(op, atom);
        let children: Vec<CellId> = self.children(cell).collect();
        for child in children {
            let child_copy = self.clone_subtree(child);
            self.append_child(copy, child_copy);
        }
        copy
    }

    /// Interns the subtree at `cell` and returns its node.
    pub fn node_of(
        &self,
        graph: &mut Graph,
        cell: CellId,
    ) -> NodeId {
        if let Some(atom) = self.atom(cell) {
            return atom;
        }
        let children: Vec<NodeId> = self
            .children(cell)
            .collect::<Vec<_>>()
            .into_iter()
            .map(|c| self.node_of(graph, c))
            .collect();
        graph.node(self.op(cell), &children)
    }

    /// Interns the whole window. The result equals the node the window was
    /// carved from provided only identities were applied.
    pub fn commit(
        &self,
        graph: &mut Graph,
    ) -> NodeId {
        self.node_of(graph, self.root)
    }

    // ------------------------------------------------------------------
    // Generic normalisation
    // ------------------------------------------------------------------

    /// Flattens nested associative operators and collapses associative
    /// cells with a single child. Returns `true` if anything changed.
    pub fn flatten(
        &mut self,
        graph: &Graph,
    ) -> bool {
        let mut changed = false;
        for cell in self.post_order() {
            if !self.cell(cell).live || self.atom(cell).is_some() {
                continue;
            }
            let op = self.op(cell);
            if !graph.ops().get(op).flags.has(OpFlags::ASSOCIATIVE) {
                continue;
            }
            let nested: Vec<CellId> = self
                .children(cell)
                .filter(|&c| self.atom(c).is_none() && self.op(c) == op)
                .collect();
            for child in nested {
                self.splice_up(child);
                changed = true;
            }
            if self.arity(cell) == 1 {
                let only = self.cell(cell).first;
                self.detach(only);
                self.replace(cell, only);
                changed = true;
            }
        }
        changed
    }

    // ------------------------------------------------------------------
    // Rewriting
    // ------------------------------------------------------------------

    /// Finds the first way `rewrite` matches at `cell`, checking guards.
    pub fn find_match(
        &self,
        graph: &mut Graph,
        rewrite: &Rewrite,
        cell: CellId,
    ) -> Option<WindowMatch> {
        let Pat::Node(op, pats) = &rewrite.lhs else {
            return None;
        };
        if self.atom(cell).is_some() || self.op(cell) != *op {
            return None;
        }
        let children: Vec<CellId> = self.children(cell).collect();
        let flags = graph.ops().get(*op).flags;
        let mut found = WindowMatch {
            subst: vec![Vec::new(); rewrite.nvars],
            ops: vec![OpId::NONE; rewrite.nvars],
            rest: Vec::new(),
        };
        let mut search = Search {
            window: self,
            graph,
            rewrite,
            rest: Vec::new(),
            ops: found.ops.clone(),
        };
        let matched = if flags.has(OpFlags::COMMUTATIVE) {
            let allow_rest = flags.has(OpFlags::ASSOCIATIVE);
            (pats.len() == children.len() || (allow_rest && pats.len() < children.len())) && {
                let mut used = vec![false; children.len()];
                search.assign(
                    pats,
                    &children,
                    0,
                    &mut used,
                    &mut Vec::new(),
                    &mut found.subst,
                    Mode::Root,
                )
            }
        } else {
            pats.len() == children.len() && {
                let mut goals: Vec<(&Pat, CellId)> =
                    pats.iter().zip(children.iter().copied()).rev().collect();
                search.solve(&mut goals, &mut found.subst)
            }
        };
        found.rest = std::mem::take(&mut search.rest);
        found.ops = std::mem::take(&mut search.ops);
        matched.then_some(found)
    }

    /// Replaces `cell` by the right-hand side of `rewrite` under a match
    /// returned by [`TreeWindow::find_match`]. Returns the new cell.
    pub fn apply(
        &mut self,
        graph: &mut Graph,
        rewrite: &Rewrite,
        cell: CellId,
        found: &WindowMatch,
    ) -> Option<CellId> {
        let mut built = self.instantiate(graph, &rewrite.rhs, found)?;
        if !found.rest.is_empty() {
            let mut children = vec![built];
            for &r in &found.rest {
                children.push(self.clone_subtree(r));
            }
            built = self.new_node(self.op(cell), &children);
        }
        self.replace(cell, built);
        Some(built)
    }

    fn instantiate(
        &mut self,
        graph: &mut Graph,
        pat: &Pat,
        found: &WindowMatch,
    ) -> Option<CellId> {
        match pat {
            | Pat::Var(v) => match found.subst.get(*v as usize)?.as_slice() {
                | [] => None,
                | [only] => Some(self.clone_subtree(*only)),
                | several => {
                    let copies: Vec<CellId> =
                        several.iter().map(|&c| self.clone_subtree(c)).collect();
                    Some(self.new_node(*found.ops.get(*v as usize)?, &copies))
                },
            },
            | Pat::Num(n) => {
                let node = graph.num(n.clone());
                Some(self.new_atom(graph, node))
            },
            | Pat::Exact(node) => Some(self.new_atom(graph, *node)),
            | Pat::Node(op, pats) => {
                if !graph.ops().get(*op).arity.accepts(pats.len()) {
                    return None;
                }
                if pats.is_empty() {
                    let node = graph.try_node(*op, &[])?;
                    return Some(self.new_atom(graph, node));
                }
                let mut children = Vec::with_capacity(pats.len());
                for p in pats {
                    children.push(self.instantiate(graph, p, found)?);
                }
                Some(self.new_node(*op, &children))
            },
        }
    }

    /// Whether two lists of cells are equal as multisets of subtrees.
    fn same_multiset(
        &self,
        graph: &Graph,
        a: &[CellId],
        b: &[CellId],
    ) -> bool {
        if a.len() != b.len() {
            return false;
        }
        let sorted = |cells: &[CellId]| {
            let mut v = cells.to_vec();
            v.sort_by_key(|&c| self.fingerprint(graph, c));
            v
        };
        sorted(a)
            .iter()
            .zip(sorted(b))
            .all(|(&x, y)| self.same(graph, x, y))
    }

    /// Greedy local optimisation: flattens, runs `passes` and applies
    /// `rewrites` bottom-up until nothing changes or `max_steps` changes
    /// were made. Returns the number of changes (rewrites applied plus
    /// rounds in which flattening or a pass changed the window).
    ///
    /// The rewrites should be simplifying (tier
    /// [`Normalize`](super::rule::Tier::Normalize)); nothing here protects
    /// against a rule set that loops other than the step budget.
    pub fn optimize(
        &mut self,
        graph: &mut Graph,
        rewrites: &[&Rewrite],
        passes: &[&dyn WindowPass],
        max_steps: usize,
    ) -> usize {
        let mut steps = 0_usize;
        let mut rounds = 0_usize;
        loop {
            rounds = rounds.saturating_add(1);
            let mut changed = self.flatten(graph);
            for pass in passes {
                changed |= pass.run(graph, self);
            }
            if changed {
                steps = steps.saturating_add(1);
            }
            'cells: for cell in self.post_order() {
                for rewrite in rewrites {
                    if !self.cell(cell).live {
                        continue 'cells;
                    }
                    if let Some(found) = self.find_match(graph, rewrite, cell)
                        && self.apply(graph, rewrite, cell, &found).is_some() {
                            changed = true;
                            steps = steps.saturating_add(1);
                            // The cell list is stale now; restart the sweep.
                            break 'cells;
                        }
                }
            }
            // `rounds` also bounds passes that keep reporting changes.
            if !changed || steps >= max_steps || rounds > max_steps {
                return steps;
            }
        }
    }

    /// Beam search over structure-changing rewrites.
    ///
    /// Every state in the beam is expanded by applying each of `moves` at
    /// each cell; successors are normalised with [`TreeWindow::optimize`],
    /// deduplicated by fingerprint and the `width` cheapest survive. After
    /// `depth` rounds (or when no new state appears) the cheapest window
    /// ever seen is returned — possibly the start.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn beam(
        &self,
        graph: &mut Graph,
        moves: &[&Rewrite],
        rewrites: &[&Rewrite],
        passes: &[&dyn WindowPass],
        width: usize,
        depth: usize,
        max_steps: usize,
    ) -> Self {
        let mut best = self.clone();
        let mut best_cost = best.cost(graph);
        let mut seen: HashSet<u64> = HashSet::new();
        seen.insert(self.fingerprint(graph, self.root));
        let mut frontier = vec![self.clone()];
        for _ in 0..depth {
            let mut next: Vec<(u64, Self)> = Vec::new();
            for state in &frontier {
                for cell in state.post_order() {
                    for rewrite in moves {
                        let Some(found) = state.find_match(graph, rewrite, cell) else {
                            continue;
                        };
                        let mut succ = state.clone();
                        if succ.apply(graph, rewrite, cell, &found).is_none() {
                            continue;
                        }
                        succ.optimize(graph, rewrites, passes, max_steps);
                        if seen.insert(succ.fingerprint(graph, succ.root)) {
                            next.push((succ.cost(graph), succ));
                        }
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            next.sort_by_key(|(cost, _)| *cost);
            next.truncate(width.max(1));
            if let Some((cost, state)) = next.first()
                && *cost < best_cost {
                    best_cost = *cost;
                    best = state.clone();
                }
            frontier = next.into_iter().map(|(_, s)| s).collect();
        }
        best
    }
}

/// A match of a rewrite inside a [`TreeWindow`].
#[derive(Clone, Debug, Default)]
pub struct WindowMatch {
    /// Cells bound to each pattern variable: none if unbound, one for an
    /// ordinary binding, several when the variable absorbed the leftover
    /// operands of a nested associative operator.
    pub subst: Vec<Vec<CellId>>,
    /// For variables bound to several cells, the operator joining them.
    pub ops: Vec<OpId>,
    /// Unmatched children of an associative-commutative root.
    pub rest: Vec<CellId>,
}

/// Where a commutative assignment takes place.
#[derive(Copy, Clone)]
enum Mode {
    /// At the pattern root: leftover children become the match's rest.
    Root,
    /// Below the root: every child must be consumed.
    Nested,
    /// Below the root, with a variable absorbing the leftover children.
    Absorb(u32, OpId),
}

/// Backtracking matcher of a pattern against window cells.
struct Search<'a> {
    window: &'a TreeWindow,
    graph: &'a mut Graph,
    rewrite: &'a Rewrite,
    rest: Vec<CellId>,
    ops: Vec<OpId>,
}

impl Search<'_> {
    #[allow(clippy::too_many_arguments)]
    fn assign<'p>(
        &mut self,
        pats: &'p [Pat],
        children: &[CellId],
        index: usize,
        used: &mut Vec<bool>,
        goals: &mut Vec<(&'p Pat, CellId)>,
        subst: &mut Vec<Vec<CellId>>,
        mode: Mode,
    ) -> bool {
        let Some(pat) = pats.get(index) else {
            let left: Vec<CellId> = children
                .iter()
                .zip(used.iter())
                .filter(|(_, u)| !**u)
                .map(|(c, _)| *c)
                .collect();
            let mut pending = goals.clone();
            return match mode {
                | Mode::Root => {
                    self.rest = left;
                    self.solve(&mut pending, subst)
                },
                | Mode::Nested => self.solve(&mut pending, subst),
                | Mode::Absorb(var, op) => {
                    let index = var as usize;
                    if let Some(slot) = subst.get_mut(index) {
                        *slot = left;
                    }
                    if let Some(slot) = self.ops.get_mut(index) {
                        *slot = op;
                    }
                    let ok = self.solve(&mut pending, subst);
                    if !ok
                        && let Some(slot) = subst.get_mut(index) {
                            slot.clear();
                        }
                    ok
                },
            };
        };
        for (i, &child) in children.iter().enumerate() {
            if used.get(i).copied().unwrap_or(true) {
                continue;
            }
            if let Some(slot) = used.get_mut(i) {
                *slot = true;
            }
            goals.push((pat, child));
            let saved = subst.clone();
            if self.assign(
                pats,
                children,
                index.saturating_add(1),
                used,
                goals,
                subst,
                mode,
            ) {
                return true;
            }
            subst.clone_from(&saved);
            goals.pop();
            if let Some(slot) = used.get_mut(i) {
                *slot = false;
            }
        }
        false
    }

    /// Discharges all goals; on success `subst` holds the bindings and the
    /// guards have been checked.
    fn solve(
        &mut self,
        goals: &mut Vec<(&Pat, CellId)>,
        subst: &mut Vec<Vec<CellId>>,
    ) -> bool {
        let Some((pat, cell)) = goals.pop() else {
            return self.guards_hold(subst);
        };
        let window = self.window;
        match pat {
            | Pat::Var(v) => {
                let index = *v as usize;
                let bound = subst.get(index).cloned().unwrap_or_default();
                match bound.as_slice() {
                    | [] => {
                        if let Some(slot) = subst.get_mut(index) {
                            *slot = vec![cell];
                        }
                        if self.solve(goals, subst) {
                            return true;
                        }
                        if let Some(slot) = subst.get_mut(index) {
                            slot.clear();
                        }
                        false
                    },
                    | [only] => window.same(self.graph, *only, cell) && self.solve(goals, subst),
                    | several => {
                        let operands: Vec<CellId> = window.children(cell).collect();
                        window.atom(cell).is_none()
                            && self.ops.get(index).is_some_and(|&op| op == window.op(cell))
                            && window.same_multiset(self.graph, several, &operands)
                            && self.solve(goals, subst)
                    },
                }
            },
            | Pat::Num(n) => window.number(self.graph, cell) == Some(n) && self.solve(goals, subst),
            | Pat::Exact(node) => {
                window.atom(cell).is_some_and(|a| self.graph.same(a, *node))
                    && self.solve(goals, subst)
            },
            | Pat::Node(op, pats) => {
                if window.op(cell) != *op {
                    return false;
                }
                if window.atom(cell).is_some() {
                    // Only a nullary operator can match an atom structurally.
                    return pats.is_empty() && self.solve(goals, subst);
                }
                let children: Vec<CellId> = window.children(cell).collect();
                let flags = self.graph.ops().get(*op).flags;
                let mut inner = goals.clone();
                if flags.has(OpFlags::COMMUTATIVE) && children.len() > pats.len() {
                    let Some((Pat::Var(v), front)) = pats.split_last() else {
                        return false;
                    };
                    let free = subst.get(*v as usize).is_some_and(Vec::is_empty);
                    if front.is_empty() || !free || !flags.has(OpFlags::ASSOCIATIVE) {
                        return false;
                    }
                    let mut used = vec![false; children.len()];
                    return self.assign(
                        front,
                        &children,
                        0,
                        &mut used,
                        &mut inner,
                        subst,
                        Mode::Absorb(*v, *op),
                    );
                }
                if children.len() != pats.len() {
                    return false;
                }
                if flags.has(OpFlags::COMMUTATIVE) {
                    let mut used = vec![false; children.len()];
                    self.assign(
                        pats,
                        &children,
                        0,
                        &mut used,
                        &mut inner,
                        subst,
                        Mode::Nested,
                    )
                } else {
                    inner.extend(pats.iter().zip(children.iter().copied()).rev());
                    self.solve(&mut inner, subst)
                }
            },
        }
    }

    fn guards_hold(
        &mut self,
        subst: &[Vec<CellId>],
    ) -> bool {
        if self.rewrite.guards.is_empty() {
            return true;
        }
        let mut nodes = Vec::with_capacity(subst.len());
        for (index, cells) in subst.iter().enumerate() {
            let node = match cells.as_slice() {
                | [] => NodeId::NONE,
                | [only] => self.window.node_of(self.graph, *only),
                | several => {
                    let operands: Vec<NodeId> = several
                        .iter()
                        .map(|&c| self.window.node_of(self.graph, c))
                        .collect();
                    let op = self.ops.get(index).copied().unwrap_or(OpId::NONE);
                    self.graph.try_node(op, &operands).unwrap_or(NodeId::NONE)
                },
            };
            nodes.push(node);
        }
        self.rewrite.admits(self.graph, &nodes)
    }
}

#[cold]
fn bad_cell(id: CellId) -> ! {
    panic!("rssn graph kernel: window cell {id} does not exist")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;
    use crate::graph::op::core;
    use crate::graph::rule::Action;
    use crate::graph::rule::Program;
    use crate::graph::rule::RuleSet;
    use crate::graph::rule::Tier;

    fn graph() -> Graph {
        let mut g = Graph::new();
        for name in ["sin", "cos", "f"] {
            assert!(
                g.ops_mut()
                    .register(OpDescriptor::new(name, Arity::Fixed(1)))
                    .is_ok()
            );
        }
        g
    }

    fn rewrites(
        g: &mut Graph,
        texts: &'static [&'static str],
    ) -> Vec<Rewrite> {
        let set = RuleSet::new("t", move |i| i.rewrites(Tier::Normalize, texts));
        let mut program = Program::default();
        assert_eq!(set.install_into(g, &mut program), Ok(()));
        program
            .rules
            .into_iter()
            .filter_map(|r| match r.action {
                | Action::Rewrite(rw) => Some(rw),
                | Action::Kernel(_) => None,
            })
            .collect()
    }

    fn shape(
        g: &Graph,
        w: &TreeWindow,
        cell: CellId,
    ) -> String {
        if let Some(atom) = w.atom(cell) {
            return format!("[{}]", g.display(atom));
        }
        let kids: Vec<String> = w.children(cell).map(|c| shape(g, w, c)).collect();
        format!("{}({})", g.ops().get(w.op(cell)).name, kids.join(" "))
    }

    #[test]
    fn carve_commit_round_trip() {
        let mut g = graph();
        let node = g.parse("sin(a + b) * (c + f(d))^2").unwrap_or(NodeId::NONE);
        let w = TreeWindow::carve(&g, node, 64);
        assert_eq!(w.commit(&mut g), node);
        assert_eq!(w.post_order().len(), w.len());
    }

    #[test]
    fn shared_subterms_become_atoms() {
        let mut g = graph();
        // `f(a)` is used by both sin and cos: it has two parents.
        let node = g.parse("sin(f(a)) * cos(f(a))").unwrap_or(NodeId::NONE);
        let w = TreeWindow::carve(&g, node, 64);
        assert_eq!(shape(&g, &w, w.root()), "mul(sin([f(a)]) cos([f(a)]))");
        // Arithmetic is transparent: shared or not, the window sees it.
        let node = g.parse("sin(a + b) * cos(a + b)").unwrap_or(NodeId::NONE);
        let w = TreeWindow::carve(&g, node, 64);
        assert_eq!(shape(&g, &w, w.root()), "mul(sin(add([a] [b])) cos(add([a] [b])))");
        // Used once: inlined.
        let node = g.parse("sin(p + q)").unwrap_or(NodeId::NONE);
        let w = TreeWindow::carve(&g, node, 64);
        assert_eq!(shape(&g, &w, w.root()), "sin(add([p] [q]))");
    }

    #[test]
    fn cell_budget_is_respected() {
        let mut g = graph();
        let node = g.parse("f(f(f(f(f(f(x))))))").unwrap_or(NodeId::NONE);
        let w = TreeWindow::carve(&g, node, 3);
        assert!(w.len() <= 5, "window has {} cells", w.len());
        assert_eq!(w.commit(&mut g), node);
    }

    #[test]
    fn splice_up_flattens_in_place() {
        let mut g = graph();
        let (a, b, c, d) = (g.sym("a"), g.sym("b"), g.sym("c"), g.sym("d"));
        let mut w = TreeWindow::carve(&g, a, 8);
        let cells: Vec<CellId> = [a, b, c, d].iter().map(|&n| w.new_atom(&g, n)).collect();
        let [ca, cb, cc, cd] = cells[..] else {
            panic!("four cells")
        };
        let inner = w.new_node(core::ADD, &[cb, cc]);
        let outer = w.new_node(core::ADD, &[ca, inner, cd]);
        let old_root = w.root();
        w.replace(old_root, outer);
        assert_eq!(shape(&g, &w, w.root()), "add([a] add([b] [c]) [d])");
        assert!(w.flatten(&g));
        assert_eq!(shape(&g, &w, w.root()), "add([a] [b] [c] [d])");
        assert_eq!(w.arity(w.root()), 4);
        assert_eq!(w.len(), 5);
        assert!(!w.flatten(&g));
    }

    #[test]
    fn greedy_rewriting_with_rest_and_nonlinear_vars() {
        let mut g = graph();
        let rws = rewrites(
            &mut g,
            &[
                "pyth: sin(?x)^2 + cos(?x)^2 => 1",
                "add0: ?a + 0 => ?a",
                "mul1: ?a * 1 => ?a",
            ],
        );
        let refs: Vec<&Rewrite> = rws.iter().collect();
        let node = g
            .parse("k * (cos(u + v)^2 + 0 + sin(v + u)^2) * 1")
            .unwrap_or(NodeId::NONE);
        let before = g.len();
        let mut w = TreeWindow::carve(&g, node, 64);
        let steps = w.optimize(&mut g, &refs, &[], 100);
        assert!(steps >= 2);
        let out = w.commit(&mut g);
        assert_eq!(g.display(out), "k");
        // Only the literals introduced by rules and guard-free commits were
        // interned: the intermediates never reached the graph.
        assert!(g.len() <= before + 1, "graph grew by {}", g.len() - before);
    }

    #[test]
    fn nested_ac_patterns_absorb_in_windows() {
        let mut g = graph();
        let rws = rewrites(&mut g, &["pull: f(?c * ?r) => ?c * f(?r) if number(?c)"]);
        let refs: Vec<&Rewrite> = rws.iter().collect();
        let node = g.parse("f(x * 3 * y)").unwrap_or(NodeId::NONE);
        let mut w = TreeWindow::carve(&g, node, 64);
        assert_eq!(w.optimize(&mut g, &refs, &[], 10), 1);
        let out = w.commit(&mut g);
        assert_eq!(g.display(out), "3*f(x*y)");
        // Nothing numeric to pull out: the guard rejects every choice.
        let node = g.parse("f(x * y * z)").unwrap_or(NodeId::NONE);
        let mut w = TreeWindow::carve(&g, node, 64);
        assert_eq!(w.optimize(&mut g, &refs, &[], 10), 0);
    }

    #[test]
    fn guards_are_checked_in_windows() {
        let mut g = graph();
        let rws = rewrites(&mut g, &["cancel: ?a * ?a^(-1) => 1 if nonzero(?a)"]);
        let refs: Vec<&Rewrite> = rws.iter().collect();
        let unknown = g.parse("f(x * x^(-1))").unwrap_or(NodeId::NONE);
        let mut w = TreeWindow::carve(&g, unknown, 64);
        assert_eq!(w.optimize(&mut g, &refs, &[], 10), 0, "x may be zero");
        // Built by hand: the parser would fold `3 * 3^(-1)` to `1` itself.
        let f = g.ops().lookup("f").unwrap_or(OpId::NONE);
        let (three, minus_one) = (g.int(3), g.int(-1));
        let inverse = g.node(core::POW, &[three, minus_one]);
        let product = g.node(core::MUL, &[three, inverse]);
        let known = g.node(f, &[product]);
        let mut w = TreeWindow::carve(&g, known, 64);
        assert_eq!(w.optimize(&mut g, &refs, &[], 10), 1);
        let out = w.commit(&mut g);
        assert_eq!(g.display(out), "f(1)");
    }

    #[test]
    #[allow(clippy::tuple_array_conversions)] // false positive: the tuple is a destructuring of separate values, not a conversion
    fn fingerprint_ignores_commutative_order() {
        let mut g = graph();
        let (a, b) = (g.sym("a"), g.sym("b"));
        let mut w = TreeWindow::carve(&g, a, 8);
        let (a1, b1, a2, b2) = (
            w.new_atom(&g, a),
            w.new_atom(&g, b),
            w.new_atom(&g, a),
            w.new_atom(&g, b),
        );
        let x = w.new_node(core::ADD, &[a1, b1]);
        let y = w.new_node(core::ADD, &[b2, a2]);
        assert_eq!(w.fingerprint(&g, x), w.fingerprint(&g, y));
        assert!(w.same(&g, x, y));
        let (a3, b3) = (w.new_atom(&g, a), w.new_atom(&g, b));
        let p = w.new_node(core::POW, &[a3, b3]);
        let (a4, b4) = (w.new_atom(&g, a), w.new_atom(&g, b));
        let q = w.new_node(core::POW, &[b4, a4]);
        assert!(!w.same(&g, p, q));
    }

    #[test]
    fn beam_escapes_a_local_minimum() {
        let mut g = graph();
        let simplify = rewrites(
            &mut g,
            &["pyth: sin(?x)^2 + cos(?x)^2 => 1", "mul1: ?a * 1 => ?a"],
        );
        let moves = rewrites(&mut g, &["factor: ?a*?b + ?a*?c => ?a*(?b + ?c)"]);
        let simplify: Vec<&Rewrite> = simplify.iter().collect();
        let moves: Vec<&Rewrite> = moves.iter().collect();
        // Greedy simplification is stuck; factoring k out exposes the identity.
        let node = g.parse("k*sin(t)^2 + k*cos(t)^2").unwrap_or(NodeId::NONE);
        let mut w = TreeWindow::carve(&g, node, 64);
        assert_eq!(w.optimize(&mut g, &simplify, &[], 100), 0);
        let best = w.beam(&mut g, &moves, &simplify, &[], 4, 3, 100);
        let out = best.commit(&mut g);
        assert_eq!(g.display(out), "k");
    }
}
