//! Compact 32-bit identifiers used throughout the graph kernel.

use std::fmt;

macro_rules! define_id {
    ($(#[$meta:meta])* $name:ident, $prefix:literal) => {
        $(#[$meta])*
        #[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub(crate) u32);

        impl $name {
            /// Sentinel meaning "no value".
            pub const NONE: Self = Self(u32::MAX);

            /// Builds an identifier from a raw index.
            #[must_use]
            pub const fn from_raw(raw: u32) -> Self {
                Self(raw)
            }

            /// Returns the raw index.
            #[must_use]
            pub const fn raw(self) -> u32 {
                self.0
            }

            /// Returns the identifier as a `usize` index.
            #[must_use]
            pub const fn index(self) -> usize {
                self.0 as usize
            }

            /// Returns `true` when this is the [`Self::NONE`] sentinel.
            #[must_use]
            pub const fn is_none(self) -> bool {
                self.0 == u32::MAX
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                if self.is_none() {
                    write!(f, concat!($prefix, "-"))
                } else {
                    write!(f, concat!($prefix, "{}"), self.0)
                }
            }
        }
    };
}

define_id!(
    /// A concrete term stored in the DAG. Every node denotes exactly one
    /// expression; structurally equal expressions share one `NodeId`.
    NodeId,
    "n"
);
define_id!(
    /// An equivalence class of nodes, named by its current representative.
    /// Only stable between two unions; re-canonicalise with `Graph::find`.
    ClassId,
    "c"
);
define_id!(
    /// An operator registered in the [`OpTable`](super::op::OpTable).
    OpId,
    "op"
);
define_id!(
    /// An interned symbol (variable or function name).
    SymbolId,
    "s"
);
define_id!(
    /// An interned leaf payload.
    PayloadId,
    "p"
);

/// Converts a collection length into a `u32` index.
///
/// The kernel addresses everything with 32-bit indices; exceeding that range
/// is a capacity error that cannot be recovered from locally.
#[must_use]
pub(crate) fn to_u32(len: usize) -> u32 {
    u32::try_from(len).unwrap_or_else(|_| capacity_overflow())
}

#[cold]
fn capacity_overflow() -> u32 {
    panic!("rssn graph kernel: more than u32::MAX entities allocated")
}
