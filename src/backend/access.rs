//! Scoped restrictions for program execution hosted by a language interpreter.

use std::cell::Cell;
use std::marker::PhantomData;
use std::rc::Rc;

use crate::kernel::OpError;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Policy {
    Host,
    SourceWithGpu,
    SourceCpuOnly,
}

thread_local! {
    static POLICY: Cell<Policy> = const { Cell::new(Policy::Host) };
}

/// Restores the previous source execution policy on this thread.
#[must_use]
pub struct SourceAccess {
    previous: Option<Policy>,
    _thread: PhantomData<Rc<()>>,
}

impl Drop for SourceAccess {
    fn drop(&mut self) {
        if let Some(previous) = self.previous {
            POLICY.set(previous);
        }
    }
}

/// Restricts calls made by embedded source code on this thread.
///
/// Nested scopes only narrow access. Source calls cannot mutate the host's
/// process default. A host installs the scope separately on each worker.
pub fn restrict_to_source(allow_gpu: bool) -> SourceAccess {
    let previous = POLICY.get();
    let effective = if allow_gpu && previous != Policy::SourceCpuOnly {
        Policy::SourceWithGpu
    } else {
        Policy::SourceCpuOnly
    };
    SourceAccess {
        previous: (previous != effective).then(|| POLICY.replace(effective)),
        _thread: PhantomData,
    }
}

/// Captures the source grant for work that will enter a different thread.
pub fn source_gpu_grant() -> Option<bool> {
    match POLICY.get() {
        | Policy::Host => None,
        | Policy::SourceWithGpu => Some(true),
        | Policy::SourceCpuOnly => Some(false),
    }
}

/// Whether the active source scope permits GPU execution.
#[must_use]
pub fn gpu_allowed() -> bool {
    POLICY.get() != Policy::SourceCpuOnly
}

/// Reports a source permission error before acquiring GPU resources.
///
/// # Errors
/// Returns `E0902` when GPU execution is denied by the active source scope.
pub fn require_gpu() -> Result<(), OpError> {
    if gpu_allowed() {
        Ok(())
    } else {
        Err(OpError::HostDiagnostic {
            code: "E0902",
            message: "permission denied: the `gpu` capability is not granted here".into(),
        })
    }
}

/// Checks the GPU grant when a selected or retained object uses a GPU.
///
/// # Errors
/// Returns `E0902` for GPU-backed work that the active source scope denies.
pub fn require_gpu_if(required: bool) -> Result<(), OpError> {
    if required {
        require_gpu()
    } else {
        Ok(())
    }
}

/// Reserves mutation of the process default for the embedding host.
///
/// # Errors
/// Returns `E0902` when called from an active source scope.
pub fn require_process_default() -> Result<(), OpError> {
    if POLICY.get() == Policy::Host {
        Ok(())
    } else {
        Err(OpError::HostDiagnostic {
            code: "E0902",
            message: "the process RSSN backend default is controlled by the host".into(),
        })
    }
}

/// Recognizes built-in selectors that can resolve without probing devices.
#[must_use]
pub fn is_host_selector(selector: &str) -> bool {
    matches!(
        selector,
        "auto" | "interpreter" | "cpu" | "cpu-parallel" | "tiered" | "cranelift" | "jit"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_access_narrows_and_restores_each_thread_scope() {
        assert!(require_process_default().is_ok());
        assert!(require_gpu().is_ok());
        {
            let _outer = restrict_to_source(true);
            assert!(matches!(
                require_process_default(),
                Err(OpError::HostDiagnostic { code: "E0902", .. })
            ));
            assert!(require_gpu().is_ok());
            {
                let _denied = restrict_to_source(false);
                let _nested = restrict_to_source(true);
                assert!(matches!(
                    require_gpu(),
                    Err(OpError::HostDiagnostic { code: "E0902", .. })
                ));
            }
            assert!(require_gpu().is_ok());
        }
        assert!(require_process_default().is_ok());
    }
}
