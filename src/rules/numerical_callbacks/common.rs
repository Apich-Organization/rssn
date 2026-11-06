//! Shared serial callback execution and registration contracts.
use super::Arc;
use super::Arity;
use super::BackendSupport;
use super::CallCtx;
use super::Determinism;
use super::Effect;
use super::EffectRow;
use super::EvalDyn;
use super::Installer;
use super::OpDescriptor;
use super::OpError;
use super::Purity;
use super::RefCell;
use super::RuleError;
use super::Signature;
use super::Stability;
use super::ThreadSafety;
use super::ToPrimitive;
use super::TypeExpr;
use super::Value;
use crate::kernel::signature::ArityInfo;
use crate::kernel::value::CallbackRef;
use std::cell::Cell;

/// The first callback/conversion failure wins; a failed operation never calls
/// the host again just because its native solver requests another sample.
pub(super) struct CallbackRun<'a, 'ctx> {
    context: &'a CallCtx<'ctx>,
    handle: CallbackRef,
    failure: RefCell<Option<OpError>>,
    evaluations: Cell<usize>,
}
impl<'a, 'ctx> CallbackRun<'a, 'ctx> {
    pub(super) fn new(
        context: &'a CallCtx<'ctx>,
        value: &Value,
        index: usize,
        arity: u8,
        expected: TypeExpr,
    ) -> Result<(Self, Signature), OpError> {
        let Value::Function(handle) = value else {
            return Err(OpError::Type {
                expected: Box::new(expected),
                got: Box::new(context.value_type(value)?),
                arg: index,
            });
        };
        let signature = context.callback_signature(*handle)?;
        if signature.params.len() != usize::from(arity) {
            return Err(OpError::Arity {
                expected: ArityInfo::fixed(arity),
                got: signature.params.len(),
            });
        }
        Ok((
            Self {
                context,
                handle: *handle,
                failure: RefCell::new(None),
                evaluations: Cell::new(0),
            },
            signature.clone(),
        ))
    }

    pub(super) fn evaluate<T>(
        &self,
        arguments: &[Value],
        decode: impl FnOnce(Value) -> Result<T, OpError>,
    ) -> Option<T> {
        if self.failure.borrow().is_some() {
            return None;
        }
        self.evaluations
            .set(self.evaluations.get().saturating_add(1));
        match self
            .context
            .invoke_callback(self.handle, arguments)
            .and_then(decode)
        {
            | Ok(value) => Some(value),
            | Err(error) => {
                *self.failure.borrow_mut() = Some(error);
                None
            },
        }
    }

    pub(super) fn failed(&self) -> bool {
        self.failure.borrow().is_some()
    }

    pub(super) const fn evaluations(&self) -> usize {
        self.evaluations.get()
    }

    pub(super) fn finish(self) -> Result<(), OpError> {
        match self.failure.into_inner() {
            | Some(error) => Err(error),
            | None => Ok(()),
        }
    }
}

pub(super) fn positive_budget(
    value: &Value,
    index: usize,
    label: &str,
) -> Result<usize, OpError> {
    let Value::Int(number) = value else {
        return Err(OpError::Type {
            expected: Box::new(TypeExpr::Int),
            got: Box::new(value.type_expr()),
            arg: index,
        });
    };
    number
        .to_usize()
        .filter(|count| *count > 0)
        .ok_or_else(|| OpError::Domain(format!("{label} must be positive and fit usize")))
}

pub(super) fn register<E: EvalDyn + 'static>(
    installer: &mut Installer<'_>,
    name: &str,
    path: &str,
    doc: &str,
    arity: u8,
    evaluator: E,
) -> Result<(), RuleError> {
    let signature = evaluator.signature();
    installer.op(OpDescriptor::new(name, Arity::Fixed(arity))
        .rust_path(path)
        .category("numerical_callback")
        .doc(doc)
        .stability(Stability::Experimental)
        .purity(Purity::Effectful)
        .determinism(Determinism::EnvDependent)
        .thread_safety(ThreadSafety::SessionBound)
        .backend(BackendSupport {
            interpreter: true,
            jit: false,
            ffi: false,
            gpu: false,
        })
        .higher_order()
        .signature(signature)
        .eval_dyn(Arc::new(evaluator)))?;
    Ok(())
}

pub(super) fn effectful(mut signature: Signature) -> Signature {
    signature.effects = EffectRow {
        effects: vec![Effect::Stateful, Effect::Diverges],
    };
    signature
}
