use std::fmt;

/// An opaque identifier representing an equivalence class (E-Class) in an E-Graph.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Id(pub u32);

impl Id {
    /// Creates a new `Id` from a `u32`.
    #[inline(always)]
    #[must_use]
    pub const fn from_u32(val: u32) -> Self {
        Self(val)
    }

    /// Creates a new `Id` from a `usize`.
    #[inline(always)]
    #[must_use]
    pub const fn from_usize(val: usize) -> Self {
        Self(val as u32)
    }

    /// Converts the `Id` to a `usize` index.
    #[inline(always)]
    #[must_use]
    pub const fn as_usize(self) -> usize {
        self.0 as usize
    }

    /// Returns the raw `u32` value.
    #[inline(always)]
    #[must_use]
    pub const fn as_u32(self) -> u32 {
        self.0
    }
}

impl fmt::Debug for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "e{}", self.0)
    }
}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "e{}", self.0)
    }
}
