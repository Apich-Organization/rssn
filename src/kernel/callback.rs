//! Checked host function handles: type information comes from their callback
//! bank, never from a fabricated function type on an unscoped integer handle.
use crate::kernel::call::CallCtx;
use crate::kernel::call::OpError;
use crate::kernel::overload::dimensions_equivalent;
use crate::kernel::overload::matches_type;
use crate::kernel::signature::ArityInfo;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::NatExpr;
use crate::kernel::signature::Signature;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::CallbackRef;
use crate::kernel::value::Value;
use std::collections::BTreeSet;

/// An immutable signature snapshot alongside each registered host function.
/// Cloning a bank retains the same entries without running host code.
///
/// Slots are positional: a [`CallbackRef`] is an index into this bank, so a
/// retired slot is left in place as a tombstone instead of shifting the
/// handles of the callbacks after it. [`CallbackBank::add`] reuses the first
/// tombstone, which keeps a long-lived host (a REPL, a server, a loop that
/// builds a fresh closure per iteration) from growing the bank without bound.
#[derive(Clone, Debug, Default)]
pub struct CallbackBank {
    entries: Vec<Option<std::sync::Arc<CallbackEntry>>>,
    live: usize,
}
#[derive(Debug)]
struct CallbackEntry {
    callback: std::sync::Arc<dyn crate::kernel::plugin::HostCallback>,
    signature: Signature,
}
impl CallbackBank {
    /// An empty callback bank.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
            live: 0,
        }
    }

    /// Registered callback count; retired slots are not counted.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.live
    }

    /// Whether no callbacks are registered.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.live == 0
    }

    /// The number of slots, including retired ones. A handle is valid when it
    /// is below this and [`CallbackBank::get`] returns a callback for it.
    #[must_use]
    pub const fn slots(&self) -> usize {
        self.entries.len()
    }

    /// The registered host object, without re-querying its signature.
    #[must_use]
    pub fn get(
        &self,
        index: usize,
    ) -> Option<&std::sync::Arc<dyn crate::kernel::plugin::HostCallback>> {
        self.entries
            .get(index)
            .and_then(Option::as_ref)
            .map(|entry| &entry.callback)
    }

    /// The immutable signature captured at registration.
    #[must_use]
    pub fn signature(
        &self,
        index: usize,
    ) -> Option<&Signature> {
        self.entries
            .get(index)
            .and_then(Option::as_ref)
            .map(|entry| &entry.signature)
    }

    /// Snapshot and validate a host callback before assigning its handle.
    ///
    /// The returned handle names the slot, which is either the first retired
    /// slot or a new one at the end. A retired slot's handle is never handed
    /// out again while it is live, so a stale [`CallbackRef`] resolves to
    /// nothing rather than to an unrelated function.
    ///
    /// # Errors
    /// Returns `Unsupported` for an unrepresentable signature and `Host` if
    /// all handle indices are exhausted. Host signature code runs once here.
    pub fn add(
        &mut self,
        callback: std::sync::Arc<dyn crate::kernel::plugin::HostCallback>,
    ) -> Result<CallbackRef, OpError> {
        let signature = callback.signature();
        validate_signature(&signature)?;
        let entry = std::sync::Arc::new(CallbackEntry { callback, signature });
        if let Some(index) = self.entries.iter().position(Option::is_none) {
            self.entries[index] = Some(entry);
            self.live += 1;
            return Ok(CallbackRef(u32::try_from(index).map_err(|_| {
                OpError::Host("callback handle space exhausted".into())
            })?));
        }
        let index = u32::try_from(self.entries.len())
            .map_err(|_| OpError::Host("callback handle space exhausted".into()))?;
        self.entries.push(Some(entry));
        self.live += 1;
        Ok(CallbackRef(index))
    }

    /// Releases the callback in slot `index`, leaving a reusable tombstone.
    ///
    /// Returns `false` when the slot is out of range or already retired. A
    /// handle that has been retired no longer resolves: [`CallbackBank::get`]
    /// and [`CallbackBank::signature`] return `None` for it.
    pub fn retire(
        &mut self,
        index: usize,
    ) -> bool {
        match self.entries.get_mut(index) {
            | Some(slot @ Some(_)) => {
                *slot = None;
                self.live -= 1;
                true
            },
            | _ => false,
        }
    }
}

fn validate_signature(signature: &Signature) -> Result<(), OpError> {
    if !signature.type_params.is_empty()
        || signature.variadic.is_some()
        || !signature.constraints.is_empty()
        || !signature.params.iter().all(concrete_type)
        || !concrete_type(&signature.ret)
    {
        return Err(OpError::Unsupported {
            backend: "host callback".into(),
            feature: "callback signature must be concrete, fixed-arity and unconstrained".into(),
        });
    }
    Ok(())
}

impl CallCtx<'_> {
    /// The concrete fixed-arity signature registered for a callback handle.
    ///
    /// # Errors
    /// Returns `Host` for an invalid handle and `Unsupported` for an unbound,
    /// constrained, variadic, or unsupported structured callback signature.
    pub fn callback_signature(
        &self,
        handle: CallbackRef,
    ) -> Result<&Signature, OpError> {
        self.callbacks.signature(handle.0 as usize).ok_or_else(|| {
            OpError::Host(format!(
                "callback {} was not supplied in the context",
                handle.0
            ))
        })
    }

    /// The type of a value, resolving function handles through this context.
    ///
    /// # Errors
    /// Returns the errors of [`Self::callback_signature`] for function handles.
    pub fn value_type(
        &self,
        value: &Value,
    ) -> Result<TypeExpr, OpError> {
        match value {
            | Value::Dyn(value) => {
                let inner = self.value_type(&value.value)?;
                Ok(TypeExpr::Dim(Box::new(inner), value.dim.clone()))
            },
            | Value::Function(handle) => {
                let signature = self.callback_signature(*handle)?;
                Ok(TypeExpr::Function(
                    signature.params.clone(),
                    Box::new(signature.ret.clone()),
                ))
            },
            | Value::Record(record) => {
                Ok(TypeExpr::Record(
                    record
                        .fields
                        .iter()
                        .map(|(name, value)| self.value_type(value).map(|ty| (name.clone(), ty)))
                        .collect::<Result<_, _>>()?,
                ))
            },
            | Value::Variant(variant) => {
                Ok(TypeExpr::Variant(vec![(
                    variant.tag.clone(),
                    variant
                        .payload
                        .iter()
                        .map(|value| self.value_type(value))
                        .collect::<Result<_, _>>()?,
                )]))
            },
            | _ => Ok(value.type_expr()),
        }
    }

    /// Types of arguments, with concrete callback signatures in place of handles.
    ///
    /// # Errors
    /// Returns the errors of [`Self::value_type`].
    pub fn argument_types(
        &self,
        args: &[Value],
    ) -> Result<Vec<TypeExpr>, OpError> {
        args.iter().map(|value| self.value_type(value)).collect()
    }

    /// Invoke a callback after checking input arity, types and output evidence.
    ///
    /// User errors propagate unchanged. No callback-bank lock or graph borrow
    /// is held during invocation, so a host may make a nested kernel call.
    /// Tensor shape/data lengths and record/variant contents are checked, not
    /// merely their top-level metadata. Unknown dynamic data cannot prove a
    /// concrete callback signature. Callbacks are not implicitly parallelized.
    ///
    /// # Errors
    /// Returns registration/signature errors, `Arity`, `CallbackType`, or the
    /// callback's original error.
    pub fn invoke_callback(
        &self,
        handle: CallbackRef,
        args: &[Value],
    ) -> Result<Value, OpError> {
        let signature = self.callback_signature(handle)?;
        if args.len() != signature.params.len() {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(u8::try_from(signature.params.len()).map_err(|_| {
                    OpError::Unsupported {
                        backend: "host callback".into(),
                        feature: "callback arity exceeds ABI limit".into(),
                    }
                })?),
                got: args.len(),
            });
        }
        for (index, (value, expected)) in args.iter().zip(&signature.params).enumerate() {
            self.check_callback_value(handle, Some(index), expected, value)?;
        }
        let callback = self
            .callbacks
            .get(handle.0 as usize)
            .cloned()
            .ok_or_else(|| {
                OpError::Host(format!("callback {} disappeared from context", handle.0))
            })?;
        let result = callback.call(args)?;
        self.check_callback_value(handle, None, &signature.ret, &result)?;
        Ok(result)
    }

    fn check_callback_value(
        &self,
        handle: CallbackRef,
        argument: Option<usize>,
        expected: &TypeExpr,
        value: &Value,
    ) -> Result<(), OpError> {
        if self.accepts_callback_value(expected, value)? {
            return Ok(());
        }
        Err(OpError::CallbackType {
            callback: handle,
            argument,
            expected: Box::new(expected.clone()),
            got: Box::new(self.value_type(value)?),
        })
    }

    fn accepts_callback_value(
        &self,
        expected: &TypeExpr,
        value: &Value,
    ) -> Result<bool, OpError> {
        if matches!(value, Value::Term(_)) {
            return Ok(false);
        }
        let actual = self.value_type(value)?;
        if !concrete_type(&actual) {
            return Ok(false);
        }
        let expected = strip_dimensionless(expected);
        if !matches_type(expected, strip_dimensionless(&actual)).accepted() {
            return Ok(false);
        }
        match (expected, value) {
            | (TypeExpr::Tensor { rank, shape, elem }, Value::Tensor(tensor)) => {
                if !matches!(rank,NatExpr::Lit(rank) if *rank == tensor.shape.len() as u64)
                    || tensor
                        .shape
                        .iter()
                        .try_fold(1usize, |n, &d| n.checked_mul(d))
                        != Some(tensor.data.len())
                {
                    return Ok(false);
                }
                if let Some(extents) = shape
                    && (extents.len() != tensor.shape.len()
                        || extents
                            .iter()
                            .zip(&tensor.shape)
                            .any(|(extent, &n)| !matches!(extent,NatExpr::Lit(e) if *e==n as u64)))
                {
                    return Ok(false);
                }
                if !matches_type(strip_dimensionless(elem), strip_dimensionless(&tensor.elem))
                    .accepted()
                {
                    return Ok(false);
                }
                for item in &tensor.data {
                    let item_expected = match (elem.as_ref(), &tensor.elem, item) {
                        | (
                            TypeExpr::Dim(inner, a),
                            TypeExpr::Dim(_, b),
                            Value::Int(_) | Value::Real(_) | Value::Float(_) | Value::Complex(_),
                        ) if dimensions_equivalent(a, b) => inner.as_ref(),
                        | _ => elem,
                    };
                    if !self.accepts_callback_value(item_expected, item)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            },
            | (TypeExpr::Record(fields), Value::Record(record)) => {
                if fields.len() != record.fields.len() {
                    return Ok(false);
                }
                for (name, ty) in fields {
                    let mut found = record.fields.iter().filter(|(label, _)| label == name);
                    let Some((_, item)) = found.next() else {
                        return Ok(false);
                    };
                    if found.next().is_some() || !self.accepts_callback_value(ty, item)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            },
            | (TypeExpr::Variant(variants), Value::Variant(value)) => {
                let Some((_, fields)) = variants.iter().find(|(tag, _)| tag == &value.tag) else {
                    return Ok(false);
                };
                if fields.len() != value.payload.len() {
                    return Ok(false);
                }
                for (ty, value) in fields.iter().zip(&value.payload) {
                    if !self.accepts_callback_value(ty, value)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            },
            | _ => Ok(true),
        }
    }
}

fn strip_dimensionless(ty: &TypeExpr) -> &TypeExpr {
    if let TypeExpr::Dim(inner, dim) = ty
        && dimensions_equivalent(dim, &DimExpr::Dimless)
    {
        strip_dimensionless(inner)
    } else {
        ty
    }
}
fn concrete_dimension(dim: &DimExpr) -> bool {
    match dim {
        | DimExpr::Dimless | DimExpr::Base(_) => true,
        | DimExpr::Mul(parts) => parts.iter().all(concrete_dimension),
        | DimExpr::Inverse(inner) => concrete_dimension(inner),
        | DimExpr::Pow(inner, Exp::Rat(r)) => r.denom != 0 && concrete_dimension(inner),
        | _ => false,
    }
}
fn concrete_type(ty: &TypeExpr) -> bool {
    match ty {
        | TypeExpr::Dim(inner, dim) => {
            matches!(
                inner.as_ref(),
                TypeExpr::Int | TypeExpr::Real | TypeExpr::Float | TypeExpr::Complex
            ) && concrete_dimension(dim)
        },
        | TypeExpr::Unit(dim) => concrete_dimension(dim),
        | TypeExpr::Tensor { rank, shape, elem } => {
            let NatExpr::Lit(rank) = rank else {
                return false;
            };
            let Ok(rank) = usize::try_from(*rank) else {
                return false;
            };
            shape.as_ref().is_none_or(|shape| {
                shape.len() == rank
                    && shape.iter().all(
                        |extent| matches!(extent, NatExpr::Lit(n) if usize::try_from(*n).is_ok()),
                    )
            }) && concrete_type(elem)
        },
        | TypeExpr::Record(fields) => {
            let names: BTreeSet<_> = fields.iter().map(|(name, _)| name).collect();
            names.len() == fields.len() && fields.iter().all(|(_, ty)| concrete_type(ty))
        },
        | TypeExpr::Variant(variants) => {
            let names: BTreeSet<_> = variants.iter().map(|(name, _)| name).collect();
            names.len() == variants.len()
                && variants
                    .iter()
                    .all(|(_, fields)| fields.iter().all(concrete_type))
        },
        | TypeExpr::Function(args, ret) => args.iter().all(concrete_type) && concrete_type(ret),
        | TypeExpr::Dyn(_)
        | TypeExpr::Var(_)
        | TypeExpr::Sparse { .. }
        | TypeExpr::Table(_)
        | TypeExpr::Point(_)
        | TypeExpr::Device => false,
        | _ => true,
    }
}
