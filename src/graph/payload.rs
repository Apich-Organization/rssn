//! Interned leaf payloads.
//!
//! Nodes are fixed-size; anything that does not fit in a node (numbers,
//! symbol references, dense numeric data) lives here and is referenced by a
//! [`PayloadId`]. Payloads are deduplicated so that equal leaves hash-cons
//! to the same node.

use std::any::Any;
use std::collections::HashMap;
use std::fmt;
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;

use super::id::PayloadId;
use super::id::SymbolId;
use super::id::to_u32;
use super::number::Number;

/// Opaque, shared, immutable data attached to a leaf.
///
/// Tensors, trajectories, packed polynomials, or anything a domain or plugin
/// wants to move through the graph without the kernel understanding it.
///
/// Two `Blob`s are equal only when they are the *same allocation*; content
/// hashing of megabyte-sized tensors would defeat O(1) node construction.
#[derive(Clone)]
pub struct Blob {
    kind: &'static str,
    data: Arc<dyn Any + Send + Sync>,
}

impl Blob {
    /// Wraps `value`, tagging it with a human readable `kind` such as
    /// `"tensor"` or `"trajectory"`.
    #[must_use]
    pub fn new<T: Any + Send + Sync>(
        kind: &'static str,
        value: T,
    ) -> Self {
        Self {
            kind,
            data: Arc::new(value),
        }
    }

    /// The tag given at construction.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        self.kind
    }

    /// Borrows the content as `T` if that is what it holds.
    #[must_use]
    pub fn downcast_ref<T: Any>(&self) -> Option<&T> {
        self.data.downcast_ref::<T>()
    }

    fn addr(&self) -> usize {
        Arc::as_ptr(&self.data).cast::<()>() as usize
    }
}

impl fmt::Debug for Blob {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        write!(f, "<{}@{:x}>", self.kind, self.addr())
    }
}

impl PartialEq for Blob {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.addr() == other.addr()
    }
}

impl Eq for Blob {}

impl Hash for Blob {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {
        self.addr().hash(state);
    }
}

/// The content of a leaf node.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Payload {
    /// A literal number.
    Num(Number),
    /// A reference to an interned symbol.
    Sym(SymbolId),
    /// A truth value.
    Bool(bool),
    /// A short string attribute (names of options, units, ...).
    Str(Arc<str>),
    /// Opaque shared data.
    Blob(Blob),
}

/// Deduplicating store for [`Payload`]s and symbol names.
#[derive(Clone, Debug, Default)]
pub struct Interner {
    payloads: Vec<Payload>,
    payload_ids: HashMap<Payload, PayloadId>,
    symbols: Vec<Arc<str>>,
    symbol_ids: HashMap<Arc<str>, SymbolId>,
}

impl Interner {
    /// Interns `payload`, returning the id of the unique stored copy.
    pub fn payload(
        &mut self,
        payload: Payload,
    ) -> PayloadId {
        if let Some(&id) = self.payload_ids.get(&payload) {
            return id;
        }
        let id = PayloadId(to_u32(self.payloads.len()));
        self.payloads.push(payload.clone());
        self.payload_ids.insert(payload, id);
        id
    }

    /// Looks up a payload. Returns `None` for [`PayloadId::NONE`].
    #[must_use]
    pub fn get(
        &self,
        id: PayloadId,
    ) -> Option<&Payload> {
        self.payloads.get(id.index())
    }

    /// Interns a symbol name.
    pub fn symbol(
        &mut self,
        name: &str,
    ) -> SymbolId {
        if let Some(&id) = self.symbol_ids.get(name) {
            return id;
        }
        let id = SymbolId(to_u32(self.symbols.len()));
        let name: Arc<str> = Arc::from(name);
        self.symbols.push(Arc::clone(&name));
        self.symbol_ids.insert(name, id);
        id
    }

    /// Creates a symbol that is guaranteed not to clash with any existing
    /// one, for bound variables and integration constants.
    pub fn fresh_symbol(
        &mut self,
        stem: &str,
    ) -> SymbolId {
        let mut n = self.symbols.len();
        loop {
            let candidate = format!("{stem}#{n}");
            if !self.symbol_ids.contains_key(candidate.as_str()) {
                return self.symbol(&candidate);
            }
            n = n.saturating_add(1);
        }
    }

    /// The name of a symbol.
    #[must_use]
    pub fn symbol_name(
        &self,
        id: SymbolId,
    ) -> &str {
        self.symbols.get(id.index()).map_or("?", |s| &**s)
    }

    /// Looks a symbol up by name without creating it.
    #[must_use]
    pub fn find_symbol(
        &self,
        name: &str,
    ) -> Option<SymbolId> {
        self.symbol_ids.get(name).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn payloads_are_deduplicated() {
        let mut i = Interner::default();
        let a = i.payload(Payload::Num(Number::from(3)));
        let b = i.payload(Payload::Num(Number::from(3)));
        let c = i.payload(Payload::Num(Number::from(3.0)));
        assert_eq!(a, b);
        assert_ne!(a, c, "exact 3 and float 3.0 are different literals");
    }

    #[test]
    fn blobs_compare_by_identity() {
        let x = Blob::new("tensor", vec![1.0_f64, 2.0]);
        let y = Blob::new("tensor", vec![1.0_f64, 2.0]);
        assert_eq!(x, x.clone());
        assert_ne!(x, y);
        assert_eq!(x.downcast_ref::<Vec<f64>>().map(Vec::len), Some(2));
    }

    #[test]
    fn fresh_symbols_never_clash() {
        let mut i = Interner::default();
        let x = i.symbol("x");
        let f1 = i.fresh_symbol("x");
        let f2 = i.fresh_symbol("x");
        assert_ne!(x, f1);
        assert_ne!(f1, f2);
        assert_eq!(i.symbol("x"), x);
    }
}
