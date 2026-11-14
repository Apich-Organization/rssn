//! Session-owned, checked opaque native objects.
//!
//! Wire handles contain only a scope token and object number. Native values are
//! retained by their session, never deserialized from untrusted object bytes.
use crate::kernel::call::OpError;
use crate::kernel::signature::DimExpr;
use crate::kernel::value::OpaqueValue;
use crate::kernel::value::Value;
use std::any::Any;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt;
use std::rc::Rc;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;

static NEXT_SCOPE: AtomicU64 = AtomicU64::new(1);
const MAX_OBJECTS: usize = 4096;
const MAX_BYTES: usize = 64 * 1024 * 1024;
struct Entry {
    kind: String,
    dimensions: Vec<DimExpr>,
    retained_bytes: usize,
    value: RefCell<Box<dyn Any>>,
    _parents: Vec<Rc<Self>>,
}
struct Store {
    next: u64,
    bytes: usize,
    entries: BTreeMap<u64, Rc<Entry>>,
}
/// A bounded object store shared only by calls of one session.
#[derive(Clone)]
pub struct NativeObjects {
    scope: u64,
    store: Rc<RefCell<Store>>,
}
impl fmt::Debug for NativeObjects {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.debug_struct("NativeObjects")
            .field("scope", &self.scope)
            .finish_non_exhaustive()
    }
}
impl Default for NativeObjects {
    fn default() -> Self {
        Self::new()
    }
}
impl NativeObjects {
    /// Allocate an empty, unique object scope.
    ///
    /// # Panics
    /// Panics after the process exhausts all 64-bit scope identifiers.
    #[must_use]
    pub fn new() -> Self {
        let scope = loop {
            let current = NEXT_SCOPE.load(Ordering::Relaxed);
            let next = current
                .checked_add(1)
                .expect("native object scope space exhausted");
            if NEXT_SCOPE
                .compare_exchange_weak(current, next, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                break current;
            }
        };
        Self {
            scope,
            store: Rc::new(RefCell::new(Store {
                next: 1,
                bytes: 0,
                entries: BTreeMap::new(),
            })),
        }
    }

    /// Retain a native value with concrete unit domains and an explicit memory charge.
    ///
    /// `retained_bytes` includes dynamic allocations owned by the value. Adapters
    /// validate constructor work and allocation bounds before constructing it.
    ///
    /// # Errors
    /// Rejects unknown dimensions, exhausted object/memory budgets, or a busy store.
    pub fn insert<T: Any>(
        &self,
        kind: &str,
        dimensions: Vec<DimExpr>,
        value: T,
        retained_bytes: usize,
    ) -> Result<Value, OpError> {
        self.insert_dependent(kind, dimensions, value, retained_bytes, &[])
    }

    /// Retain a value whose backing storage depends on other native objects.
    ///
    /// Parent handles cannot be released until all dependent handles are released.
    ///
    /// # Errors
    /// Rejects invalid parents or the same bounds as [`Self::insert`].
    pub fn insert_dependent<T: Any>(
        &self,
        kind: &str,
        dimensions: Vec<DimExpr>,
        value: T,
        retained_bytes: usize,
        parents: &[Value],
    ) -> Result<Value, OpError> {
        if kind.is_empty() || kind.len() > 256 || dimensions.len() > 64 || parents.len() > 64 {
            return Err(OpError::Domain(
                "native object metadata exceeds bounds".into(),
            ));
        }
        let parent_entries = parents
            .iter()
            .map(|p| self.entry(p))
            .collect::<Result<Vec<_>, _>>()?;
        let mut metadata = kind.len()
            + std::mem::size_of::<Entry>()
            + parents.len() * std::mem::size_of::<Rc<Entry>>();
        for dim in &dimensions {
            let mut stack = vec![(dim, 0usize)];
            let mut nodes = 0usize;
            while let Some((part, depth)) = stack.pop() {
                nodes += 1;
                if nodes > 256 || depth > 32 {
                    return Err(OpError::Domain(
                        "native object dimension is too complex".into(),
                    ));
                }
                metadata = metadata
                    .checked_add(std::mem::size_of::<DimExpr>())
                    .ok_or_else(|| OpError::Domain("native metadata overflow".into()))?;
                match part {
                    | DimExpr::Dimless => {},
                    | DimExpr::Base(name) if !name.is_empty() && name.len() <= 128 => {
                        metadata += name.len();
                    },
                    | DimExpr::Mul(parts) => {
                        if parts.len() > 64 {
                            return Err(OpError::Domain(
                                "native dimension has too many factors".into(),
                            ));
                        }
                        stack.extend(parts.iter().map(|p| (p, depth + 1)));
                    },
                    | DimExpr::Inverse(inner) => stack.push((inner, depth + 1)),
                    | DimExpr::Pow(inner, crate::kernel::signature::Exp::Rat(power))
                        if power.denom > 0 =>
                    {
                        stack.push((inner, depth + 1));
                    },
                    | _ => {
                        return Err(OpError::Domain(
                            "native objects require concrete unit domains".into(),
                        ));
                    },
                }
            }
            if crate::kernel::overload::canonical_concrete_dimension(dim).is_none() {
                return Err(OpError::Domain(
                    "native objects require concrete unit domains".into(),
                ));
            }
        }
        let retained_bytes = retained_bytes
            .max(std::mem::size_of::<T>())
            .checked_add(metadata)
            .ok_or_else(|| OpError::Domain("native memory charge overflow".into()))?;
        let mut store = self.store.try_borrow_mut().map_err(|_| busy())?;
        let bytes = store
            .bytes
            .checked_add(retained_bytes)
            .filter(|n| *n <= MAX_BYTES)
            .ok_or_else(|| {
                OpError::Host(
                    "native object memory budget exhausted; release unused objects".into(),
                )
            })?;
        if store.entries.len() >= MAX_OBJECTS {
            return Err(OpError::Host(
                "native object count budget exhausted; release unused objects".into(),
            ));
        }
        let id = store.next;
        store.next = id
            .checked_add(1)
            .ok_or_else(|| OpError::Host("native object handle space exhausted".into()))?;
        store.entries.insert(
            id,
            Rc::new(Entry {
                kind: kind.into(),
                dimensions,
                retained_bytes,
                value: RefCell::new(Box::new(value)),
                _parents: parent_entries,
            }),
        );
        store.bytes = bytes;
        let mut bytes = self.scope.to_le_bytes().to_vec();
        bytes.extend_from_slice(&id.to_le_bytes());
        Ok(Value::Opaque(OpaqueValue {
            kind: kind.into(),
            bytes,
        }))
    }

    fn id(
        &self,
        handle: &Value,
    ) -> Result<u64, OpError> {
        let Value::Opaque(value) = handle else {
            return Err(OpError::Domain("expected a native opaque handle".into()));
        };
        let bytes: &[u8; 16] = value
            .bytes
            .as_slice()
            .try_into()
            .map_err(|_| OpError::Domain("invalid native handle encoding".into()))?;
        let scope = u64::from_le_bytes(
            bytes[..8]
                .try_into()
                .map_err(|_| OpError::Internal("native scope encoding".into()))?,
        );
        if scope != self.scope {
            return Err(OpError::Domain(
                "native object belongs to another session".into(),
            ));
        }
        Ok(u64::from_le_bytes(bytes[8..].try_into().map_err(|_| {
            OpError::Internal("native identifier encoding".into())
        })?))
    }

    fn entry(
        &self,
        handle: &Value,
    ) -> Result<Rc<Entry>, OpError> {
        let id = self.id(handle)?;
        let store = self.store.try_borrow().map_err(|_| busy())?;
        let entry = store
            .entries
            .get(&id)
            .ok_or_else(|| OpError::Domain("native object handle is released or unknown".into()))?;
        let Value::Opaque(value) = handle else {
            return Err(OpError::Internal(
                "native handle changed representation".into(),
            ));
        };
        if value.kind != entry.kind {
            return Err(OpError::Domain(
                "native object kind does not match its handle".into(),
            ));
        }
        Ok(entry.clone())
    }

    /// Read a native object of the expected Rust type and opaque kind.
    ///
    /// # Errors
    /// Rejects foreign/stale/forged handles, wrong types, and conflicting mutable use.
    pub fn with<T: Any, R>(
        &self,
        handle: &Value,
        kind: &str,
        f: impl FnOnce(&T) -> Result<R, OpError>,
    ) -> Result<R, OpError> {
        let entry = self.entry(handle)?;
        if entry.kind != kind {
            return Err(OpError::Domain("wrong native object kind".into()));
        }
        let value = entry.value.try_borrow().map_err(|_| busy())?;
        f(value
            .downcast_ref::<T>()
            .ok_or_else(|| OpError::Domain("wrong native object type".into()))?)
    }

    /// Mutate a native object with checked exclusive access.
    ///
    /// # Errors
    /// Rejects foreign/stale/forged handles, wrong types, and reentrant aliasing.
    pub fn with_mut<T: Any, R>(
        &self,
        handle: &Value,
        kind: &str,
        f: impl FnOnce(&mut T) -> Result<R, OpError>,
    ) -> Result<R, OpError> {
        let entry = self.entry(handle)?;
        if entry.kind != kind {
            return Err(OpError::Domain("wrong native object kind".into()));
        }
        let mut value = entry.value.try_borrow_mut().map_err(|_| busy())?;
        f(value
            .downcast_mut::<T>()
            .ok_or_else(|| OpError::Domain("wrong native object type".into()))?)
    }

    /// Unit domains recorded when the native object was constructed.
    ///
    /// # Errors
    /// Rejects a foreign, stale or forged handle.
    pub fn dimensions(
        &self,
        handle: &Value,
    ) -> Result<Vec<DimExpr>, OpError> {
        Ok(self.entry(handle)?.dimensions.clone())
    }

    /// Release one handle; all copies become invalid.
    ///
    /// # Errors
    /// Rejects foreign/stale handles and objects currently used by a callback.
    pub fn release(
        &self,
        handle: &Value,
    ) -> Result<(), OpError> {
        let entry = self.entry(handle)?;
        let id = self.id(handle)?;
        if Rc::strong_count(&entry) > 2 {
            return Err(OpError::Host(
                "native object has dependent handles; release those first".into(),
            ));
        }
        let _guard = entry.value.try_borrow_mut().map_err(|_| busy())?;
        let mut store = self.store.try_borrow_mut().map_err(|_| busy())?;
        let removed = store
            .entries
            .remove(&id)
            .ok_or_else(|| OpError::Domain("native object is already released".into()))?;
        store.bytes -= removed.retained_bytes;
        Ok(())
    }

    /// Limit newly created handles to a serial callback's dynamic extent.
    ///
    /// # Errors
    /// Rejects an active mutation of the store.
    pub fn scoped(&self) -> Result<NativeObjectScope, OpError> {
        let start = self.store.try_borrow().map_err(|_| busy())?.next;
        Ok(NativeObjectScope {
            objects: self.clone(),
            start,
            committed: false,
        })
    }

    /// Current object count and charged retained bytes.
    ///
    /// # Errors
    /// Rejects a concurrent store mutation.
    pub fn usage(&self) -> Result<(usize, usize), OpError> {
        let store = self.store.try_borrow().map_err(|_| busy())?;
        Ok((store.entries.len(), store.bytes))
    }
}
fn busy() -> OpError {
    OpError::Host("native object is already in conflicting use".into())
}

/// Invalidates native handles created during one serial callback scope.
pub struct NativeObjectScope {
    objects: NativeObjects,
    start: u64,
    committed: bool,
}
impl NativeObjectScope {
    /// Keep newly constructed objects after a successful atomic adapter result.
    /// Uncommitted scopes invalidate all handles created since the checkpoint.
    pub fn commit(mut self) {
        self.committed = true;
    }
}
impl Drop for NativeObjectScope {
    fn drop(&mut self) {
        if self.committed {
            return;
        }
        if let Ok(mut store) = self.objects.store.try_borrow_mut() {
            let removed = store.entries.split_off(&self.start);
            let charge: usize = removed.values().map(|e| e.retained_bytes).sum();
            store.bytes -= charge;
        }
    }
}
