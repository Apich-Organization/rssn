//! The graph store: a hash-consed DAG whose nodes are threaded into
//! equivalence classes by intrusive circular lists.
//!
//! # Two layers over one arena
//!
//! 1. **The DAG (structural layer).** Every node is a concrete term. Its
//!    children are concrete nodes. Building the same term twice returns the
//!    same [`NodeId`], so memory is shared maximally and construction is
//!    O(arity).
//! 2. **The e-graph (congruence layer).** Nodes proven equal are linked into
//!    a *ring*: each node stores the id of the next member of its class, and
//!    a union splices two rings in O(1) without copying anything. A
//!    union-find gives the canonical [`ClassId`]; a second hash table keyed
//!    by *canonical* children maintains congruence closure (`a = b` implies
//!    `f(a) = f(b)`).
//!
//! Parent pointers are intrusive singly linked lists as well ("use lists"),
//! concatenated in O(1) on union and walked by [`Graph::rebuild`] to repair
//! congruence.
//!
//! Nodes whose canonical form coincides with an older node of the same class
//! are flagged *redundant*: they stay valid concrete terms but are skipped
//! by matching, so the e-graph never processes the same e-node twice.

use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use std::hash::Hasher;

use super::facts::Facts;
use super::id::ClassId;
use super::id::NodeId;
use super::id::OpId;
use super::id::PayloadId;
use super::id::SymbolId;
use super::id::to_u32;
use super::number::Number;
use super::op::OpFlags;
use super::op::OpTable;
use super::op::core;
use super::payload::Interner;
use super::payload::Payload;

/// A hasher for keys that already are well-mixed 64-bit hashes.
#[derive(Default)]
pub(crate) struct PassThrough(u64);

impl Hasher for PassThrough {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(
        &mut self,
        bytes: &[u8],
    ) {
        for &b in bytes {
            self.0 = mix(self.0, u64::from(b));
        }
    }

    fn write_u64(
        &mut self,
        value: u64,
    ) {
        self.0 = value;
    }
}

type HashTable = HashMap<u64, NodeId, BuildHasherDefault<PassThrough>>;

#[inline]
const fn mix(
    h: u64,
    x: u64,
) -> u64 {
    (h.rotate_left(5) ^ x).wrapping_mul(0x517c_c1b7_2722_0a95)
}

fn hash_key(
    op: OpId,
    payload: PayloadId,
    children: &[u32],
) -> u64 {
    let mut h = mix(0x9e37_79b9_7f4a_7c15, u64::from(op.0));
    h = mix(h, u64::from(payload.0));
    for &c in children {
        h = mix(h, u64::from(c));
    }
    // Final avalanche so that the low bits used by the table are good.
    h ^= h >> 32;
    h.wrapping_mul(0xd6e8_feb8_6659_fd93)
}

/// A numeric enclosure `mid ± rad` attached to a class as a *witness*.
///
/// Approximations are deliberately not e-nodes: `∫₀¹ x² dx` and `0.3333333`
/// are not identical, so merging them would poison congruence closure. A
/// witness records what numeric kernels found out about a class while
/// leaving its exact members untouched.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Ball {
    /// Centre of the enclosure.
    pub mid: f64,
    /// Radius: an error estimate, not a rigorous bound unless a kernel says
    /// so.
    pub rad: f64,
}

impl Ball {
    /// An enclosure with zero radius.
    #[must_use]
    pub const fn exact(mid: f64) -> Self {
        Self { mid, rad: 0.0 }
    }

    /// Whether zero is certainly outside the enclosure.
    #[must_use]
    pub fn excludes_zero(&self) -> bool {
        self.mid.abs() > self.rad
    }
}

/// Largest operand count produced by flattening at construction.
const MAX_FLAT_ARITY: usize = 64;
const REDUNDANT: u8 = 1;
const NO_USE: u32 = u32::MAX;

#[derive(Clone, Debug)]
struct Node {
    op: OpId,
    payload: PayloadId,
    child_start: u32,
    child_len: u32,
    /// Union-find parent.
    parent: NodeId,
    /// Next member of the class ring.
    ring_next: NodeId,
    /// Next node in the structural hash bucket.
    struct_next: NodeId,
    /// Next node in the congruence hash bucket.
    memo_next: NodeId,
    /// Hash under which this node currently sits in the congruence table.
    chash: u64,
    /// Visit stamp used to deduplicate use-list walks.
    stamp: u32,
    flags: u8,
}

#[derive(Copy, Clone, Debug)]
struct UseCell {
    node: NodeId,
    next: u32,
}

/// Per-class analysis data. Only meaningful at class representatives.
#[derive(Clone, Debug)]
struct ClassData {
    size: u32,
    uses_head: u32,
    uses_tail: u32,
    /// Payload of the literal number this class is equal to, if known.
    konst: PayloadId,
    /// Sorted symbols the class may depend on.
    free: Vec<SymbolId>,
    approx: Option<Ball>,
    /// Concrete term that must be returned for this class, if any.
    pin: NodeId,
    /// A childless member (literal, symbol, nullary operator), if any.
    leaf: NodeId,
}

/// The expression graph. See the [module documentation](self).
#[derive(Clone, Debug)]
pub struct Graph {
    ops: OpTable,
    interner: Interner,
    nodes: Vec<Node>,
    classes: Vec<ClassData>,
    child_pool: Vec<NodeId>,
    uses: Vec<UseCell>,
    structural: HashTable,
    memo: HashTable,
    by_op: Vec<Vec<NodeId>>,
    pending: Vec<NodeId>,
    analysis_pending: Vec<NodeId>,
    stamp: u32,
    class_count: usize,
    unions: u64,
    conflicts: Vec<(Number, Number)>,
    scratch: Vec<u32>,
    assumptions: HashMap<SymbolId, Facts>,
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}

impl Graph {
    /// Creates an empty graph with the core operators.
    #[must_use]
    pub fn new() -> Self {
        Self::with_ops(OpTable::new())
    }

    /// Creates an empty graph using a prepared operator table.
    #[must_use]
    pub fn with_ops(ops: OpTable) -> Self {
        Self {
            ops,
            interner: Interner::default(),
            nodes: Vec::new(),
            classes: Vec::new(),
            child_pool: Vec::new(),
            uses: Vec::new(),
            structural: HashTable::default(),
            memo: HashTable::default(),
            by_op: Vec::new(),
            pending: Vec::new(),
            analysis_pending: Vec::new(),
            stamp: 0,
            class_count: 0,
            unions: 0,
            conflicts: Vec::new(),
            scratch: Vec::new(),
            assumptions: HashMap::new(),
        }
    }

    // ------------------------------------------------------------------
    // Registries
    // ------------------------------------------------------------------

    /// The operator table.
    #[must_use]
    pub const fn ops(&self) -> &OpTable {
        &self.ops
    }

    /// Mutable access to the operator table, for registering operators.
    pub const fn ops_mut(&mut self) -> &mut OpTable {
        &mut self.ops
    }

    /// The payload and symbol interner.
    #[must_use]
    pub const fn interner(&self) -> &Interner {
        &self.interner
    }

    /// Mutable access to the interner.
    pub const fn interner_mut(&mut self) -> &mut Interner {
        &mut self.interner
    }

    // ------------------------------------------------------------------
    // Node accessors
    // ------------------------------------------------------------------

    #[inline]
    fn n(
        &self,
        id: NodeId,
    ) -> &Node {
        match self.nodes.get(id.index()) {
            | Some(n) => n,
            | None => foreign_node(id),
        }
    }

    #[inline]
    fn n_mut(
        &mut self,
        id: NodeId,
    ) -> &mut Node {
        match self.nodes.get_mut(id.index()) {
            | Some(n) => n,
            | None => foreign_node(id),
        }
    }

    #[inline]
    fn c(
        &self,
        id: ClassId,
    ) -> &ClassData {
        match self.classes.get(id.index()) {
            | Some(c) => c,
            | None => foreign_node(NodeId(id.0)),
        }
    }

    #[inline]
    fn c_mut(
        &mut self,
        id: ClassId,
    ) -> &mut ClassData {
        match self.classes.get_mut(id.index()) {
            | Some(c) => c,
            | None => foreign_node(NodeId(id.0)),
        }
    }

    /// Operator of a node.
    #[must_use]
    pub fn op(
        &self,
        id: NodeId,
    ) -> OpId {
        self.n(id).op
    }

    /// Concrete children of a node.
    #[must_use]
    pub fn children(
        &self,
        id: NodeId,
    ) -> &[NodeId] {
        let n = self.n(id);
        let start = n.child_start as usize;
        self.child_pool
            .get(start..start.saturating_add(n.child_len as usize))
            .unwrap_or(&[])
    }

    /// Payload of a leaf node.
    #[must_use]
    pub fn payload(
        &self,
        id: NodeId,
    ) -> Option<&Payload> {
        self.interner.get(self.n(id).payload)
    }

    /// The literal number held by this very node, if it is a numeric leaf.
    #[must_use]
    pub fn as_number(
        &self,
        id: NodeId,
    ) -> Option<&Number> {
        match self.payload(id) {
            | Some(Payload::Num(n)) if self.n(id).op == core::LIT => Some(n),
            | _ => None,
        }
    }

    /// The symbol held by this very node, if it is a symbol leaf.
    #[must_use]
    pub fn as_symbol(
        &self,
        id: NodeId,
    ) -> Option<SymbolId> {
        match self.payload(id) {
            | Some(Payload::Sym(s)) if self.n(id).op == core::SYM => Some(*s),
            | _ => None,
        }
    }

    /// Whether the node duplicates an older e-node of its class. Redundant
    /// nodes are skipped by [`Graph::enodes`].
    #[must_use]
    pub fn is_redundant(
        &self,
        id: NodeId,
    ) -> bool {
        self.n(id).flags & REDUNDANT != 0
    }

    /// Number of concrete nodes.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Whether the graph holds no nodes.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Number of equivalence classes.
    #[must_use]
    pub const fn class_count(&self) -> usize {
        self.class_count
    }

    /// Number of successful unions performed so far. Schedulers use the
    /// difference between two readings to detect a fixpoint.
    #[must_use]
    pub const fn union_count(&self) -> u64 {
        self.unions
    }

    /// Pairs of distinct literal numbers that were proven equal. A non-empty
    /// result means an unsound rule fired; the graph must not be trusted.
    #[must_use]
    pub fn conflicts(&self) -> &[(Number, Number)] {
        &self.conflicts
    }

    /// All nodes built with `op`, in creation order (including redundant
    /// ones). This is the index rules use to find their candidates.
    #[must_use]
    pub fn nodes_with_op(
        &self,
        op: OpId,
    ) -> &[NodeId] {
        self.by_op.get(op.index()).map_or(&[], Vec::as_slice)
    }

    // ------------------------------------------------------------------
    // Classes
    // ------------------------------------------------------------------

    /// Canonical class of a node.
    #[must_use]
    pub fn find(
        &self,
        id: NodeId,
    ) -> ClassId {
        let mut cur = id;
        loop {
            let parent = self.n(cur).parent;
            if parent == cur {
                return ClassId(cur.0);
            }
            cur = parent;
        }
    }

    fn find_compress(
        &mut self,
        id: NodeId,
    ) -> NodeId {
        let mut cur = id;
        loop {
            let parent = self.n(cur).parent;
            if parent == cur {
                return cur;
            }
            let grand = self.n(parent).parent;
            self.n_mut(cur).parent = grand;
            cur = grand;
        }
    }

    /// Whether two nodes are known to be equal.
    #[must_use]
    pub fn same(
        &self,
        a: NodeId,
        b: NodeId,
    ) -> bool {
        self.find(a) == self.find(b)
    }

    /// All concrete members of a class, redundant ones included.
    pub fn members(
        &self,
        class: ClassId,
    ) -> impl Iterator<Item = NodeId> + '_ {
        let start = NodeId(self.find(NodeId(class.0)).0);
        let mut cur = Some(start);
        std::iter::from_fn(move || {
            let id = cur?;
            let next = self.n(id).ring_next;
            cur = (next != start).then_some(next);
            Some(id)
        })
    }

    /// The e-nodes of a class: what rules match against and extraction
    /// chooses from.
    ///
    /// Redundant members are skipped. A class that contains a *leaf* (a
    /// literal, a symbol, a nullary operator) exposes only that leaf: once
    /// something is known to equal `x` or `2`, its other spellings can
    /// never be a better answer, and matching through them (`x = x^1 =
    /// (x^1)^1 = ...`) is what makes naive saturation explode.
    pub fn enodes(
        &self,
        class: ClassId,
    ) -> impl Iterator<Item = NodeId> + '_ {
        let root = self.find(NodeId(class.0));
        let leaf = self.c(root).leaf;
        let (collapsed, ring) = if leaf.is_none() {
            (None, Some(self.members(root)))
        } else {
            (Some(leaf), None)
        };
        collapsed.into_iter().chain(
            ring.into_iter()
                .flatten()
                .filter(|&id| !self.is_redundant(id)),
        )
    }

    /// Whether `id` is currently an e-node of its class (see
    /// [`Graph::enodes`]).
    #[must_use]
    pub fn is_enode(
        &self,
        id: NodeId,
    ) -> bool {
        let leaf = self.c(self.find(id)).leaf;
        if leaf.is_none() {
            !self.is_redundant(id)
        } else {
            leaf == id
        }
    }

    /// Number of concrete members of a class.
    #[must_use]
    pub fn class_size(
        &self,
        class: ClassId,
    ) -> usize {
        self.c(self.find(NodeId(class.0))).size as usize
    }

    /// The literal number this class is known to equal.
    #[must_use]
    pub fn class_number(
        &self,
        class: ClassId,
    ) -> Option<&Number> {
        match self.interner.get(self.c(self.find(NodeId(class.0))).konst) {
            | Some(Payload::Num(n)) => Some(n),
            | _ => None,
        }
    }

    /// Shorthand for [`Graph::class_number`] of a node's class.
    #[must_use]
    pub fn number_of(
        &self,
        id: NodeId,
    ) -> Option<&Number> {
        self.class_number(self.find(id))
    }

    /// Sorted set of symbols the class may depend on.
    ///
    /// This is an over-approximation that only shrinks: when two classes
    /// are merged the intersection is kept, because a value equal to an
    /// expression without `x` does not depend on `x`.
    #[must_use]
    pub fn free_symbols(
        &self,
        class: ClassId,
    ) -> &[SymbolId] {
        &self.c(self.find(NodeId(class.0))).free
    }

    /// Whether the class may depend on `symbol`.
    #[must_use]
    pub fn depends_on(
        &self,
        class: ClassId,
        symbol: SymbolId,
    ) -> bool {
        self.free_symbols(class).binary_search(&symbol).is_ok()
    }

    /// The numeric witness of a class.
    #[must_use]
    pub fn approx(
        &self,
        class: ClassId,
    ) -> Option<Ball> {
        let root = self.find(NodeId(class.0));
        if let Some(n) = self.class_number(root) {
            return Some(Ball::exact(n.to_f64()));
        }
        self.c(root).approx
    }

    /// Records a numeric witness, keeping the tighter one if the class
    /// already has one. Returns `true` when the stored witness changed.
    pub fn set_approx(
        &mut self,
        class: ClassId,
        ball: Ball,
    ) -> bool {
        let root = self.find(NodeId(class.0));
        let slot = &mut self.c_mut(root).approx;
        match slot {
            | Some(old) if old.rad <= ball.rad => false,
            | _ => {
                *slot = Some(ball);
                true
            },
        }
    }

    /// Pins `term` as *the* answer of its class.
    ///
    /// Some identity transformations exist only to choose a form: `expand`,
    /// `factor`, `collect`. Their result is equal to many cheaper terms, so
    /// cost-based extraction would undo them. A pin tells extraction to
    /// return this concrete term for the class, verbatim. The first pin of
    /// a class wins.
    pub fn pin(
        &mut self,
        term: NodeId,
    ) {
        let root = self.find(term);
        let slot = &mut self.c_mut(root).pin;
        if slot.is_none() {
            *slot = term;
        }
    }

    /// The pinned term of a class, if any.
    #[must_use]
    pub fn pinned(
        &self,
        class: ClassId,
    ) -> Option<NodeId> {
        let pin = self.c(self.find(NodeId(class.0))).pin;
        (!pin.is_none()).then_some(pin)
    }

    /// What has been assumed about a symbol.
    #[must_use]
    pub fn assumption(
        &self,
        symbol: SymbolId,
    ) -> Facts {
        self.assumptions
            .get(&symbol)
            .copied()
            .unwrap_or(Facts::NONE)
    }

    pub(super) fn set_assumption(
        &mut self,
        symbol: SymbolId,
        facts: Facts,
    ) {
        self.assumptions.insert(symbol, facts);
    }

    /// Forgets every numeric witness. Witnesses depend on variable bindings,
    /// so a run with different bindings must start clean.
    pub fn clear_witnesses(&mut self) {
        for class in &mut self.classes {
            class.approx = None;
        }
    }

    /// The e-nodes that have `class` among their children.
    #[must_use]
    pub fn parents(
        &self,
        class: ClassId,
    ) -> Vec<NodeId> {
        let mut out = self.users(class);
        out.retain(|&n| self.is_enode(n));
        out
    }

    /// Every non-redundant node with `class` among its children, whether
    /// or not its own class has collapsed to a leaf. Congruence repair and
    /// analysis propagation need all of them.
    fn users(
        &self,
        class: ClassId,
    ) -> Vec<NodeId> {
        let root = self.find(NodeId(class.0));
        let mut out = Vec::new();
        let mut cell = self.c(root).uses_head;
        while let Some(u) = self.uses.get(cell as usize) {
            if !self.is_redundant(u.node) && !out.contains(&u.node) {
                out.push(u.node);
            }
            cell = u.next;
        }
        out
    }

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Interns a literal leaf.
    pub fn lit(
        &mut self,
        payload: Payload,
    ) -> NodeId {
        let (op, payload) = match payload {
            | Payload::Sym(_) => (core::SYM, payload),
            | other => (core::LIT, other),
        };
        let pid = self.interner.payload(payload);
        self.intern(op, pid, &[])
    }

    /// Interns a literal number.
    pub fn num(
        &mut self,
        n: Number,
    ) -> NodeId {
        self.lit(Payload::Num(n))
    }

    /// Interns an exact integer.
    pub fn int(
        &mut self,
        v: i64,
    ) -> NodeId {
        self.num(Number::from(v))
    }

    /// Interns a float.
    pub fn float(
        &mut self,
        v: f64,
    ) -> NodeId {
        self.num(Number::from(v))
    }

    /// Interns the symbol named `name`.
    pub fn sym(
        &mut self,
        name: &str,
    ) -> NodeId {
        let s = self.interner.symbol(name);
        self.symbol_node(s)
    }

    /// Interns the leaf for an already interned symbol.
    pub fn symbol_node(
        &mut self,
        symbol: SymbolId,
    ) -> NodeId {
        self.lit(Payload::Sym(symbol))
    }

    /// Builds `op(children...)`.
    ///
    /// Associative operators are flattened one level (`f(f(a,b),c)` becomes
    /// `f(a,b,c)`, up to a fixed operand count so that shared subterms are
    /// not duplicated without bound) and collapse to their only child when
    /// given exactly one;
    /// commutative operators get their children sorted. Both are identities
    /// of the operator, so the returned node is always equal to the request.
    ///
    /// # Panics
    /// Panics when the number of children does not match the operator's
    /// arity. Use [`Graph::try_node`] for untrusted input.
    pub fn node(
        &mut self,
        op: OpId,
        children: &[NodeId],
    ) -> NodeId {
        match self.try_node(op, children) {
            | Some(id) => id,
            | None => arity_mismatch(&self.ops.get(op).name, children.len()),
        }
    }

    /// Like [`Graph::node`] but returns `None` on an arity mismatch.
    pub fn try_node(
        &mut self,
        op: OpId,
        children: &[NodeId],
    ) -> Option<NodeId> {
        let desc = self.ops.get(op);
        if !desc.arity.accepts(children.len()) || desc.flags.has(OpFlags::LEAF) {
            return None;
        }
        let flags = desc.flags;
        let mut buf: Vec<NodeId> = Vec::with_capacity(children.len());
        if flags.has(OpFlags::ASSOCIATIVE) {
            for &c in children {
                let inner = self.children(c);
                // Flattening copies the child's operands. Without a cap,
                // `s = s + s` repeated n times would build 2^n operands out
                // of n shared nodes.
                if self.n(c).op == op && buf.len().saturating_add(inner.len()) <= MAX_FLAT_ARITY {
                    buf.extend_from_slice(inner);
                } else {
                    buf.push(c);
                }
            }
            if op == core::ADD || op == core::MUL {
                self.merge_literals(op, &mut buf);
            }
            if let [only] = buf.as_slice() {
                return Some(*only);
            }
        } else {
            buf.extend_from_slice(children);
        }
        if flags.has(OpFlags::COMMUTATIVE) {
            buf.sort_unstable();
        }
        Some(self.intern(op, PayloadId::NONE, &buf))
    }

    /// Combines the literal operands of a sum or product into one and
    /// drops it when it is the identity element: `2 * x * 3` is built as
    /// `6 * x`, `(-1) * (-1) * x` as `x`. An identity of the operator, and
    /// without it repeated negation piles up ever longer products of `-1`
    /// in one class.
    fn merge_literals(
        &mut self,
        op: OpId,
        buf: &mut Vec<NodeId>,
    ) {
        let literals = buf.iter().filter(|&&c| self.as_number(c).is_some()).count();
        let identity = Number::from(i64::from(op == core::MUL));
        let has_identity = buf.iter().any(|&c| self.as_number(c) == Some(&identity));
        if literals < 2 && !has_identity {
            return;
        }
        let mut value = identity.clone();
        let mut rest = Vec::with_capacity(buf.len());
        for &c in buf.iter() {
            match self.as_number(c) {
                | Some(n) => value = if op == core::ADD { value.add(n) } else { value.mul(n) },
                | None => rest.push(c),
            }
        }
        if value != identity || rest.is_empty() {
            let literal = self.num(value);
            rest.insert(0, literal);
        }
        *buf = rest;
    }

    fn intern(
        &mut self,
        op: OpId,
        payload: PayloadId,
        children: &[NodeId],
    ) -> NodeId {
        // --- structural layer -------------------------------------------
        let mut key = std::mem::take(&mut self.scratch);
        key.clear();
        key.extend(children.iter().map(|c| c.0));
        let shash = hash_key(op, payload, &key);
        let mut cur = self.structural.get(&shash).copied().unwrap_or(NodeId::NONE);
        while !cur.is_none() {
            let n = self.n(cur);
            if n.op == op && n.payload == payload && self.children(cur) == children {
                self.scratch = key;
                return cur;
            }
            cur = n.struct_next;
        }

        let id = NodeId(to_u32(self.nodes.len()));
        let child_start = to_u32(self.child_pool.len());
        self.child_pool.extend_from_slice(children);
        let bucket = self.structural.insert(shash, id).unwrap_or(NodeId::NONE);
        self.nodes.push(Node {
            op,
            payload,
            child_start,
            child_len: to_u32(children.len()),
            parent: id,
            ring_next: id,
            struct_next: bucket,
            memo_next: NodeId::NONE,
            chash: 0,
            stamp: 0,
            flags: 0,
        });
        let data = self.make(id);
        self.classes.push(data);
        if self.by_op.len() <= op.index() {
            self.by_op
                .resize_with(op.index().saturating_add(1), Vec::new);
        }
        if let Some(list) = self.by_op.get_mut(op.index()) {
            list.push(id);
        }

        // --- congruence layer -------------------------------------------
        self.canonical_key(id, &mut key);
        let chash = hash_key(op, payload, &key);
        if let Some(rep) = self.memo_lookup(chash, op, payload, &key) {
            // Congruent to an existing e-node: join its class as a redundant
            // member. The fresh singleton has no parents, so nothing needs
            // repairing.
            let root = self.find_compress(rep);
            let node = self.n_mut(id);
            node.flags |= REDUNDANT;
            node.parent = root;
            self.splice(root, id);
            let class = self.c_mut(ClassId(root.0));
            class.size = class.size.saturating_add(1);
        } else {
            self.class_count = self.class_count.saturating_add(1);
            self.memo_insert(chash, id);
            key.sort_unstable();
            key.dedup();
            for &child in &key {
                self.push_use(ClassId(child), id);
            }
        }
        self.scratch = key;
        id
    }

    /// Writes the canonical child classes of `id` into `out`.
    fn canonical_key(
        &self,
        id: NodeId,
        out: &mut Vec<u32>,
    ) {
        out.clear();
        out.extend(self.children(id).iter().map(|&c| self.find(c).0));
        if self.ops.get(self.n(id).op).flags.has(OpFlags::COMMUTATIVE) {
            out.sort_unstable();
        }
    }

    fn memo_lookup(
        &self,
        chash: u64,
        op: OpId,
        payload: PayloadId,
        key: &[u32],
    ) -> Option<NodeId> {
        let mut cur = self.memo.get(&chash).copied().unwrap_or(NodeId::NONE);
        let mut other = Vec::new();
        while !cur.is_none() {
            let n = self.n(cur);
            if n.op == op && n.payload == payload && n.child_len as usize == key.len() {
                self.canonical_key(cur, &mut other);
                if other == key {
                    return Some(cur);
                }
            }
            cur = n.memo_next;
        }
        None
    }

    fn memo_insert(
        &mut self,
        chash: u64,
        id: NodeId,
    ) {
        let head = self.memo.insert(chash, id).unwrap_or(NodeId::NONE);
        let node = self.n_mut(id);
        node.chash = chash;
        node.memo_next = head;
    }

    fn memo_remove(
        &mut self,
        id: NodeId,
    ) {
        let chash = self.n(id).chash;
        let next = self.n(id).memo_next;
        let Some(&head) = self.memo.get(&chash) else {
            return;
        };
        if head == id {
            if next.is_none() {
                self.memo.remove(&chash);
            } else {
                self.memo.insert(chash, next);
            }
        } else {
            let mut cur = head;
            while !cur.is_none() {
                let after = self.n(cur).memo_next;
                if after == id {
                    self.n_mut(cur).memo_next = next;
                    break;
                }
                cur = after;
            }
        }
        self.n_mut(id).memo_next = NodeId::NONE;
    }

    fn push_use(
        &mut self,
        class: ClassId,
        user: NodeId,
    ) {
        let cell = to_u32(self.uses.len());
        self.uses.push(UseCell {
            node: user,
            next: NO_USE,
        });
        let tail = self.c(class).uses_tail;
        if let Some(prev) = self.uses.get_mut(tail as usize) {
            prev.next = cell;
        } else {
            self.c_mut(class).uses_head = cell;
        }
        self.c_mut(class).uses_tail = cell;
    }

    /// Splices the rings containing `a` and `b`, which must be disjoint.
    fn splice(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) {
        let an = self.n(a).ring_next;
        let bn = self.n(b).ring_next;
        self.n_mut(a).ring_next = bn;
        self.n_mut(b).ring_next = an;
    }

    /// Computes the analysis data of a single node from its children.
    fn make(
        &self,
        id: NodeId,
    ) -> ClassData {
        let node = self.n(id);
        let mut data = ClassData {
            size: 1,
            uses_head: NO_USE,
            uses_tail: NO_USE,
            konst: PayloadId::NONE,
            free: Vec::new(),
            approx: None,
            pin: NodeId::NONE,
            // A nullary *request* (a defined constant such as a matrix)
            // is not a final answer and must not hide what it reduces to.
            leaf: if node.child_len == 0 && !self.ops.get(node.op).flags.has(OpFlags::HEAVY) {
                id
            } else {
                NodeId::NONE
            },
        };
        match self.interner.get(node.payload) {
            | Some(Payload::Num(_)) if node.op == core::LIT => data.konst = node.payload,
            | Some(Payload::Sym(s)) if node.op == core::SYM => data.free.push(*s),
            | _ => {},
        }
        let children = self.children(id);
        let binder = self.ops.get(node.op).binder.and_then(|b| {
            let var = children
                .get(usize::from(b.var))
                .and_then(|&v| self.as_symbol(v))?;
            Some((usize::from(b.var), b.scope, var))
        });
        for (i, &child) in children.iter().enumerate() {
            let free = self.free_symbols(self.find(child));
            match binder {
                | Some((var_index, _, _)) if var_index == i => {},
                | Some((_, scope, var)) if i < 32 && scope & (1 << i) != 0 => {
                    data.free.extend(free.iter().copied().filter(|&s| s != var));
                },
                | _ => data.free.extend_from_slice(free),
            }
        }
        data.free.sort_unstable();
        data.free.dedup();
        data
    }

    /// Copies the concrete term `node` of another graph into this one.
    ///
    /// Operators are matched by name (and registered here if missing),
    /// symbols by name, other payloads by value. Equalities known to `src`
    /// are not transferred: only the term itself.
    pub fn import(
        &mut self,
        src: &Self,
        node: NodeId,
    ) -> NodeId {
        let mut copied: HashMap<NodeId, NodeId> = HashMap::new();
        let mut stack = vec![(node, false)];
        while let Some((cur, expanded)) = stack.pop() {
            if copied.contains_key(&cur) {
                continue;
            }
            let children = src.children(cur);
            if !expanded && !children.is_empty() {
                stack.push((cur, true));
                stack.extend(
                    children
                        .iter()
                        .filter(|c| !copied.contains_key(c))
                        .map(|&c| (c, false)),
                );
                continue;
            }
            let desc = src.ops.get(src.op(cur));
            let op = match self.ops.lookup(&desc.name) {
                | Some(op) => op,
                | None => self.ops.register(desc.clone()).unwrap_or(OpId::NONE),
            };
            let payload = match src.payload(cur) {
                | Some(Payload::Sym(s)) => {
                    let symbol = self.interner.symbol(src.interner.symbol_name(*s));
                    self.interner.payload(Payload::Sym(symbol))
                },
                | Some(other) => self.interner.payload(other.clone()),
                | None => PayloadId::NONE,
            };
            let kids: Vec<NodeId> = children
                .iter()
                .filter_map(|c| copied.get(c).copied())
                .collect();
            // The source already is in canonical form: intern verbatim.
            let new = self.intern(op, payload, &kids);
            copied.insert(cur, new);
        }
        copied.get(&node).copied().unwrap_or(NodeId::NONE)
    }

    // ------------------------------------------------------------------
    // Equality
    // ------------------------------------------------------------------

    /// Asserts that `a` and `b` denote the same value.
    ///
    /// Returns `true` when two distinct classes were merged. Congruence is
    /// restored lazily: call [`Graph::rebuild`] before matching again.
    pub fn union(
        &mut self,
        a: NodeId,
        b: NodeId,
    ) -> bool {
        let ra = self.find_compress(a);
        let rb = self.find_compress(b);
        if ra == rb {
            return false;
        }
        let (win, lose) = if self.c(ClassId(ra.0)).size >= self.c(ClassId(rb.0)).size {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.n_mut(lose).parent = win;
        self.splice(win, lose);
        self.class_count = self.class_count.saturating_sub(1);
        self.unions = self.unions.saturating_add(1);

        let lost = std::mem::replace(
            self.c_mut(ClassId(lose.0)),
            ClassData {
                size: 0,
                uses_head: NO_USE,
                uses_tail: NO_USE,
                konst: PayloadId::NONE,
                free: Vec::new(),
                approx: None,
                pin: NodeId::NONE,
                leaf: NodeId::NONE,
            },
        );
        // Concatenate use lists.
        if lost.uses_head != NO_USE {
            let tail = self.c(ClassId(win.0)).uses_tail;
            if let Some(prev) = self.uses.get_mut(tail as usize) {
                prev.next = lost.uses_head;
            } else {
                self.c_mut(ClassId(win.0)).uses_head = lost.uses_head;
            }
            self.c_mut(ClassId(win.0)).uses_tail = lost.uses_tail;
        }
        // Merge analysis data.
        // Two literal numbers in one class: a contradiction if both are
        // exact, or if they differ by more than rounding allows. An exact
        // value next to its own floating-point approximation is what the
        // numeric phase produces and is fine; the exact one is kept.
        let (conflict, prefer_lost) = {
            let a = self.interner.get(self.c(ClassId(win.0)).konst);
            let b = self.interner.get(lost.konst);
            match (a, b) {
                | (Some(Payload::Num(x)), Some(Payload::Num(y))) if x != y => {
                    let (fx, fy) = (x.to_f64(), y.to_f64());
                    let close = (fx - fy).abs() <= 1e-6 * (1.0 + fx.abs().max(fy.abs()));
                    let contradiction = (x.is_exact() && y.is_exact()) || !close;
                    (contradiction.then(|| (x.clone(), y.clone())), !x.is_exact() && y.is_exact())
                },
                | (None, Some(_)) => (None, true),
                | _ => (None, false),
            }
        };
        if let Some(pair) = conflict {
            self.conflicts.push(pair);
        }
        let class = self.c_mut(ClassId(win.0));
        class.size = class.size.saturating_add(lost.size);
        let had_konst = !class.konst.is_none();
        if prefer_lost {
            class.konst = lost.konst;
        }
        if class.pin.is_none() {
            class.pin = lost.pin;
        }
        // The class's leaf is its literal number if it has one, else
        // whichever childless member was there first.
        if prefer_lost || class.leaf.is_none() {
            class.leaf = if lost.leaf.is_none() { class.leaf } else { lost.leaf };
        }
        let before = class.free.len();
        class.free.retain(|s| lost.free.binary_search(s).is_ok());
        let shrunk = class.free.len() != before || class.free.len() != lost.free.len();
        class.approx = match (class.approx, lost.approx) {
            | (Some(x), Some(y)) => Some(if x.rad <= y.rad { x } else { y }),
            | (x, y) => x.or(y),
        };
        if shrunk || (!had_konst && !lost.konst.is_none()) {
            self.analysis_pending.push(win);
        }
        self.pending.push(win);
        true
    }

    /// Restores congruence closure and re-propagates analysis data after a
    /// batch of unions. Returns the number of additional unions implied by
    /// congruence.
    pub fn rebuild(&mut self) -> u64 {
        let before = self.unions;
        loop {
            while let Some(class) = self.pending.pop() {
                self.repair(class);
            }
            let Some(class) = self.analysis_pending.pop() else {
                break;
            };
            self.propagate(class);
        }
        self.unions.saturating_sub(before)
    }

    fn next_stamp(&mut self) -> u32 {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            for n in &mut self.nodes {
                n.stamp = 0;
            }
            self.stamp = 1;
        }
        self.stamp
    }

    /// Re-canonicalises the parents of `class` after it absorbed another.
    fn repair(
        &mut self,
        class: NodeId,
    ) {
        let root = self.find_compress(class);
        let stamp = self.next_stamp();
        let mut cell = {
            let data = self.c_mut(ClassId(root.0));
            data.uses_tail = NO_USE;
            std::mem::replace(&mut data.uses_head, NO_USE)
        };
        let mut key = Vec::new();
        while let Some(&UseCell { node: user, next }) = self.uses.get(cell as usize) {
            cell = next;
            if self.n(user).flags & REDUNDANT != 0 || self.n(user).stamp == stamp {
                continue;
            }
            self.n_mut(user).stamp = stamp;
            self.memo_remove(user);
            self.canonical_key(user, &mut key);
            let (op, payload) = (self.n(user).op, self.n(user).payload);
            let chash = hash_key(op, payload, &key);
            if let Some(twin) = self.memo_lookup(chash, op, payload, &key) {
                self.n_mut(user).flags |= REDUNDANT;
                self.union(user, twin);
            } else {
                self.memo_insert(chash, user);
                let now = self.find_compress(root);
                self.push_use(ClassId(now.0), user);
            }
        }
    }

    /// Recomputes the analysis data of the parents of `class`.
    fn propagate(
        &mut self,
        class: NodeId,
    ) {
        let root = self.find(class);
        for user in self.users(root) {
            let fresh = self.make(user);
            let target = self.find_compress(user);
            let data = self.c_mut(ClassId(target.0));
            let before = data.free.len();
            data.free.retain(|s| fresh.free.binary_search(s).is_ok());
            if data.free.len() != before {
                self.analysis_pending.push(target);
            }
        }
    }

    // ------------------------------------------------------------------
    // Diagnostics
    // ------------------------------------------------------------------

    /// Checks the structural invariants of the store.
    ///
    /// Intended for tests and for debugging rule sets; it is O(n²) in the
    /// worst case.
    ///
    /// # Errors
    /// Returns a description of the first violated invariant.
    pub fn validate(&self) -> Result<(), String> {
        if !self.pending.is_empty() || !self.analysis_pending.is_empty() {
            return Err("validate() called with pending unions; call rebuild() first".to_owned());
        }
        let mut roots = 0_usize;
        let mut seen: HashMap<(OpId, PayloadId, Vec<u32>), ClassId> = HashMap::new();
        let mut key = Vec::new();
        for i in 0..self.nodes.len() {
            let id = NodeId(to_u32(i));
            let class = self.find(id);
            if class.0 == id.0 {
                roots = roots.saturating_add(1);
                let ring: Vec<NodeId> = self.members(class).collect();
                if ring.len() != self.c(class).size as usize {
                    return Err(format!(
                        "{class:?}: ring has {} members, size says {}",
                        ring.len(),
                        self.c(class).size
                    ));
                }
                if let Some(stray) = ring.iter().find(|&&m| self.find(m) != class) {
                    return Err(format!(
                        "{class:?}: ring contains {stray:?} of another class"
                    ));
                }
            }
            self.canonical_key(id, &mut key);
            let node = self.n(id);
            let entry = (node.op, node.payload, key.clone());
            match seen.get(&entry) {
                | Some(&other) if other != class => {
                    return Err(format!(
                        "congruence violated: {id:?} in {class:?} matches a node of {other:?}"
                    ));
                },
                | Some(_) => {},
                | None => {
                    seen.insert(entry, class);
                },
            }
            if !self.is_redundant(id) {
                let chash = hash_key(node.op, node.payload, &key);
                if self.memo_lookup(chash, node.op, node.payload, &key) != Some(id) {
                    return Err(format!(
                        "{id:?} is an e-node but the congruence table does not resolve to it"
                    ));
                }
                for &child in self.children(id) {
                    if !self.users(self.find(child)).contains(&id) {
                        return Err(format!(
                            "{id:?} missing from the use list of its child {child:?}"
                        ));
                    }
                }
            }
        }
        if roots != self.class_count {
            return Err(format!(
                "class_count is {} but {roots} roots exist",
                self.class_count
            ));
        }
        Ok(())
    }
}

#[cold]
fn foreign_node(id: NodeId) -> ! {
    panic!("rssn graph kernel: {id:?} does not belong to this graph")
}

#[cold]
fn arity_mismatch(
    name: &str,
    got: usize,
) -> ! {
    panic!("rssn graph kernel: operator `{name}` cannot take {got} children")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::op::Arity;
    use crate::graph::op::OpDescriptor;

    fn f_op(g: &mut Graph) -> OpId {
        g.ops_mut()
            .register(OpDescriptor::new("f", Arity::Fixed(1)))
            .unwrap_or(OpId::NONE)
    }

    #[test]
    fn structural_sharing() {
        let mut g = Graph::new();
        let x = g.sym("x");
        let y = g.sym("y");
        let a = g.node(core::ADD, &[x, y]);
        let b = g.node(core::ADD, &[y, x]);
        assert_eq!(a, b, "commutative children are sorted");
        assert_eq!(g.sym("x"), x);
        assert_eq!(g.len(), 3);
        assert_eq!(g.validate(), Ok(()));
    }

    #[test]
    fn associative_flattening_and_singleton_collapse() {
        let mut g = Graph::new();
        let (x, y, z) = (g.sym("x"), g.sym("y"), g.sym("z"));
        let xy = g.node(core::ADD, &[x, y]);
        let nested = g.node(core::ADD, &[xy, z]);
        let flat = g.node(core::ADD, &[x, y, z]);
        assert_eq!(nested, flat);
        assert_eq!(g.node(core::ADD, &[x]), x);
    }

    #[test]
    fn arity_is_checked() {
        let mut g = Graph::new();
        let x = g.sym("x");
        assert_eq!(g.try_node(core::POW, &[x]), None);
        assert_eq!(g.try_node(core::SYM, &[]), None, "leaves need a payload");
    }

    #[test]
    fn union_propagates_by_congruence() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let (a, b) = (g.sym("a"), g.sym("b"));
        let fa = g.node(f, &[a]);
        let fb = g.node(f, &[b]);
        let ffa = g.node(f, &[fa]);
        let ffb = g.node(f, &[fb]);
        assert!(!g.same(ffa, ffb));
        assert!(g.union(a, b));
        assert_eq!(g.rebuild(), 2, "f(a)=f(b) and f(f(a))=f(f(b)) follow");
        assert!(g.same(fa, fb));
        assert!(g.same(ffa, ffb));
        assert_eq!(g.validate(), Ok(()));
        assert_eq!(g.class_count(), 3);
    }

    #[test]
    fn ring_enumerates_members_and_skips_redundant() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let (a, b) = (g.sym("a"), g.sym("b"));
        let fa = g.node(f, &[a]);
        let fb = g.node(f, &[b]);
        g.union(a, b);
        g.rebuild();
        let class = g.find(fa);
        let mut members: Vec<NodeId> = g.members(class).collect();
        members.sort_unstable();
        assert_eq!(members, vec![fa, fb]);
        assert_eq!(g.enodes(class).count(), 1, "f(a) and f(b) are one e-node");
        assert_eq!(g.class_size(class), 2);
    }

    #[test]
    fn late_construction_joins_existing_class() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let (a, b) = (g.sym("a"), g.sym("b"));
        let fa = g.node(f, &[a]);
        g.union(a, b);
        g.rebuild();
        // Built after the union: a new concrete term, but the same e-node.
        let fb = g.node(f, &[b]);
        assert_ne!(fa, fb);
        assert!(g.same(fa, fb));
        assert!(g.is_redundant(fb));
        assert_eq!(g.validate(), Ok(()));
    }

    #[test]
    fn commutative_congruence() {
        let mut g = Graph::new();
        let (a, b, c) = (g.sym("a"), g.sym("b"), g.sym("c"));
        let ac = g.node(core::MUL, &[a, c]);
        let cb = g.node(core::MUL, &[c, b]);
        g.union(a, b);
        g.rebuild();
        assert!(g.same(ac, cb));
        assert_eq!(g.validate(), Ok(()));
    }

    #[test]
    fn constants_and_conflicts() {
        let mut g = Graph::new();
        let x = g.sym("x");
        let two = g.int(2);
        assert_eq!(g.number_of(x), None);
        g.union(x, two);
        g.rebuild();
        assert_eq!(g.number_of(x), Some(&Number::from(2)));
        assert!(g.conflicts().is_empty());
        let three = g.int(3);
        g.union(x, three);
        g.rebuild();
        assert_eq!(g.conflicts().len(), 1, "2 = 3 must be reported");
    }

    #[test]
    fn a_float_next_to_the_exact_value_it_approximates_is_not_a_conflict() {
        let mut g = Graph::new();
        let third = g.parse("1/3").unwrap_or(NodeId::NONE);
        let approx = g.float(1.0 / 3.0);
        g.union(approx, third);
        g.rebuild();
        assert!(g.conflicts().is_empty());
        assert_eq!(g.number_of(approx).map(Number::is_exact), Some(true), "the exact value is kept");
        assert_eq!(g.enodes(g.find(approx)).collect::<Vec<_>>(), vec![third]);
        let wrong = g.float(0.34);
        g.union(wrong, third);
        g.rebuild();
        assert_eq!(g.conflicts().len(), 1, "0.34 is not 1/3");
    }

    #[test]
    fn free_symbols_shrink_on_union_and_propagate() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let x = g.sym("x");
        let sx = g.interner().find_symbol("x").unwrap_or(SymbolId::NONE);
        let minus_one = g.int(-1);
        let neg_x = g.node(core::MUL, &[minus_one, x]);
        let diff = g.node(core::ADD, &[x, neg_x]);
        let outer = g.node(f, &[diff]);
        assert!(g.depends_on(g.find(outer), sx));
        let zero = g.int(0);
        g.union(diff, zero);
        g.rebuild();
        assert!(!g.depends_on(g.find(diff), sx));
        assert!(
            !g.depends_on(g.find(outer), sx),
            "f(x - x) no longer depends on x"
        );
    }

    #[test]
    fn binders_hide_their_variable() {
        let mut g = Graph::new();
        let int = g
            .ops_mut()
            .register(OpDescriptor::new("integral", Arity::Fixed(4)).binder(1, 0b1))
            .unwrap_or(OpId::NONE);
        let (x, a) = (g.sym("x"), g.sym("a"));
        let body = g.node(core::MUL, &[a, x]);
        let (lo, hi) = (g.int(0), g.sym("b"));
        let node = g.node(int, &[body, x, lo, hi]);
        let names: Vec<&str> = g
            .free_symbols(g.find(node))
            .iter()
            .map(|&s| g.interner().symbol_name(s))
            .collect();
        assert_eq!(names, vec!["a", "b"]);
    }

    #[test]
    fn witnesses_keep_the_tighter_ball() {
        let mut g = Graph::new();
        let x = g.sym("x");
        let c = g.find(x);
        assert!(g.set_approx(c, Ball { mid: 1.0, rad: 0.1 }));
        assert!(!g.set_approx(c, Ball { mid: 1.0, rad: 0.5 }));
        assert!(g.set_approx(c, Ball { mid: 1.01, rad: 0.01 }));
        assert_eq!(g.approx(c).map(|b| b.mid), Some(1.01));
        g.clear_witnesses();
        assert_eq!(g.approx(c), None);
        let two = g.int(2);
        assert_eq!(g.approx(g.find(two)), Some(Ball::exact(2.0)));
    }

    #[test]
    fn classes_with_a_leaf_expose_only_the_leaf() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let x = g.sym("x");
        let one = g.int(1);
        let x_pow_1 = g.node(core::POW, &[x, one]);
        let outer = g.node(f, &[x_pow_1]);
        assert_eq!(g.enodes(g.find(x_pow_1)).count(), 1);
        g.union(x, x_pow_1);
        g.rebuild();
        let class = g.find(x);
        assert_eq!(
            g.enodes(class).collect::<Vec<_>>(),
            vec![x],
            "x^1 is no longer matched against"
        );
        assert_eq!(g.members(class).count(), 2, "but it is still a member");
        assert!(g.is_enode(x) && !g.is_enode(x_pow_1));
        // `x^1` dropped out of the parents of `x` and of `1`.
        assert!(!g.parents(g.find(one)).contains(&x_pow_1));
        // Congruence still sees through it.
        let direct = g.node(f, &[x]);
        assert!(g.same(outer, direct));
        assert_eq!(g.validate(), Ok(()));
        // A number beats a symbol as the leaf.
        let two = g.int(2);
        g.union(x, two);
        g.rebuild();
        assert_eq!(g.enodes(g.find(x)).collect::<Vec<_>>(), vec![two]);
    }

    #[test]
    fn import_copies_terms_between_graphs() {
        let mut a = Graph::new();
        let fa = f_op(&mut a);
        let inner = a.parse("(x + 1/2)^2 * y").unwrap_or(NodeId::NONE);
        let term = a.node(fa, &[inner]);
        let mut b = Graph::new();
        // Different symbol numbering in the destination.
        b.sym("y");
        let copy = b.import(&a, term);
        assert_eq!(b.display(copy), a.display(term));
        assert_eq!(b.import(&a, term), copy, "importing twice shares");
        assert!(
            b.ops().lookup("f").is_some(),
            "unknown operators are registered"
        );
        assert_eq!(b.validate(), Ok(()));
    }

    #[test]
    fn use_lists_survive_chained_unions() {
        let mut g = Graph::new();
        let f = f_op(&mut g);
        let leaves: Vec<NodeId> = (0..8).map(|i| g.sym(&format!("v{i}"))).collect();
        let apps: Vec<NodeId> = leaves.iter().map(|&l| g.node(f, &[l])).collect();
        for pair in leaves.windows(2) {
            if let [a, b] = pair {
                g.union(*a, *b);
            }
        }
        g.rebuild();
        assert_eq!(g.validate(), Ok(()));
        let class = g.find(*apps.first().unwrap_or(&NodeId::NONE));
        assert!(apps.iter().all(|&n| g.find(n) == class));
        let leaf_class = g.find(*leaves.first().unwrap_or(&NodeId::NONE));
        assert_eq!(g.parents(leaf_class).len(), 1);
    }
}
