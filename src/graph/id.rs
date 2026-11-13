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

/// A fast, non-cryptographic hasher for the engine's internal tables
/// (FxHash: one rotate, xor and multiply per word).
///
/// The keys are node, class and operator ids the engine itself allocates,
/// never attacker-chosen data, so `SipHash`'s flooding resistance buys
/// nothing there and costs a large share of a small request.
#[derive(Clone, Copy, Default)]
pub struct FxHasher(u64);

const FX_SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

impl FxHasher {
    #[inline]
    const fn add(
        &mut self,
        word: u64,
    ) {
        self.0 = (self.0.rotate_left(5) ^ word).wrapping_mul(FX_SEED);
    }
}

impl std::hash::Hasher for FxHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    #[inline]
    fn write(
        &mut self,
        bytes: &[u8],
    ) {
        for chunk in bytes.chunks(8) {
            let mut word = [0_u8; 8];
            word.iter_mut().zip(chunk).for_each(|(w, b)| *w = *b);
            self.add(u64::from_le_bytes(word));
        }
    }

    #[inline]
    fn write_u8(
        &mut self,
        i: u8,
    ) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_u32(
        &mut self,
        i: u32,
    ) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_u64(
        &mut self,
        i: u64,
    ) {
        self.add(i);
    }

    #[inline]
    fn write_usize(
        &mut self,
        i: usize,
    ) {
        self.add(i as u64);
    }
}

/// A `HashMap` keyed by engine ids, hashed with [`FxHasher`].
pub type FastMap<K, V> = std::collections::HashMap<K, V, std::hash::BuildHasherDefault<FxHasher>>;
