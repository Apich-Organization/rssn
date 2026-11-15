//! Cooperative cancellation for embedded computations.
//!
//! A host installs an [`Arc<AtomicBool>`] for one evaluation scope and sets it
//! from another thread to request cancellation. Checkpoints return an error;
//! they never panic, kill a worker, or present an unfinished answer as complete.
//! Scope inheritance is explicit when work moves to a Rayon worker thread.

use std::cell::RefCell;
use std::marker::PhantomData;
use std::rc::Rc;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;

thread_local! {
    static CURRENT: RefCell<Option<Arc<AtomicBool>>> = const { RefCell::new(None) };
}

/// A computation was interrupted by its embedding host.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(
        &self,
        formatter: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        formatter.write_str("computation cancelled")
    }
}
impl std::error::Error for Cancelled {}

/// Restores the previous cancellation scope when dropped on its owning thread.
/// The guard cannot be sent to another thread; send [`current()`]'s token and
/// install a fresh scope there instead.
#[must_use]
pub struct Guard {
    previous: Option<Arc<AtomicBool>>,
    owning_thread: PhantomData<Rc<()>>,
}

impl Drop for Guard {
    fn drop(&mut self) {
        CURRENT.with(|current| *current.borrow_mut() = self.previous.take());
    }
}

/// Installs a cancellation token on the calling thread until the guard drops.
/// Nested scopes restore the outer token, including when a computation fails.
#[must_use]
pub fn install(token: Arc<AtomicBool>) -> Guard {
    install_optional(Some(token))
}

/// Installs a captured worker scope, including an explicitly absent token.
/// `None` clears a previous worker token for this scope; this prevents unrelated
/// nested work stolen by Rayon from inheriting another computation's flag.
#[must_use]
pub fn install_optional(token: Option<Arc<AtomicBool>>) -> Guard {
    Guard {
        previous: CURRENT.with(|current| current.replace(token)),
        owning_thread: PhantomData,
    }
}

/// Captures the current token for explicit inheritance by a worker thread.
#[must_use]
pub fn current() -> Option<Arc<AtomicBool>> {
    CURRENT.with(|current| current.borrow().clone())
}

/// Whether the current scope has received a cancellation request.
#[must_use]
pub fn is_cancelled() -> bool {
    CURRENT.with(|current| {
        current
            .borrow()
            .as_ref()
            .is_some_and(|token| token.load(Ordering::Relaxed))
    })
}

/// Stops at a cooperative checkpoint if the current scope was cancelled.
///
/// # Errors
/// Returns [`Cancelled`] after the host sets the scope's flag.
pub fn checkpoint() -> Result<(), Cancelled> {
    if is_cancelled() {
        Err(Cancelled)
    } else {
        Ok(())
    }
}

impl From<Cancelled> for crate::kernel::OpError {
    fn from(error: Cancelled) -> Self {
        Self::HostDiagnostic {
            code: "E0812",
            message: error.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nested_scopes_restore_the_previous_token() {
        let outer = Arc::new(AtomicBool::new(false));
        let _outer = install(Arc::clone(&outer));
        {
            let _inner = install(Arc::new(AtomicBool::new(true)));
            assert_eq!(checkpoint(), Err(Cancelled));
        }
        assert_eq!(checkpoint(), Ok(()));
        outer.store(true, Ordering::Relaxed);
        assert_eq!(checkpoint(), Err(Cancelled));
    }

    #[test]
    fn worker_scope_inherits_without_contaminating_a_later_job() {
        let token = Arc::new(AtomicBool::new(false));
        let _scope = install(Arc::clone(&token));
        let inherited = current();
        let worker = std::thread::spawn(move || {
            assert!(!is_cancelled());
            {
                let _scope = inherited.map(install);
                token.store(true, Ordering::Relaxed);
                assert_eq!(checkpoint(), Err(Cancelled));
            }
            assert!(!is_cancelled());
        });
        assert!(worker.join().is_ok());
        assert!(is_cancelled());
    }
}
