use super::id::Id;

/// A disjoint-set (Union-Find) data structure with path compression and union-by-rank.
#[derive(Clone, Debug, Default)]
pub struct UnionFind {
    parents: Vec<Id>,
    ranks: Vec<u8>,
}

impl UnionFind {
    /// Creates a new, empty `UnionFind`.
    #[must_use]
    pub fn new() -> Self {
        Self {
            parents: Vec::new(),
            ranks: Vec::new(),
        }
    }

    /// Allocates a new element in its own singleton set and returns its `Id`.
    pub fn make_set(&mut self) -> Id {
        let id = Id::from_usize(self.parents.len());
        self.parents.push(id);
        self.ranks.push(0);
        id
    }

    /// Finds the canonical representative `Id` of the set containing `id`,
    /// applying path compression.
    pub fn find(&mut self, id: Id) -> Id {
        let mut root = id;
        while root != self.parents[root.as_usize()] {
            root = self.parents[root.as_usize()];
        }
        // Path compression
        let mut curr = id;
        while curr != root {
            let next = self.parents[curr.as_usize()];
            self.parents[curr.as_usize()] = root;
            curr = next;
        }
        root
    }

    /// Immutable find (does not compress path). Useful when a mutable borrow is unavailable.
    #[must_use]
    pub fn find_immut(&self, id: Id) -> Id {
        let mut root = id;
        while root != self.parents[root.as_usize()] {
            root = self.parents[root.as_usize()];
        }
        root
    }

    /// Unions the sets containing `id1` and `id2`.
    /// Returns `(canonical_root, changed)`. If `id1` and `id2` were already in the same set,
    /// returns `(root, false)`.
    pub fn union(&mut self, id1: Id, id2: Id) -> (Id, bool) {
        let root1 = self.find(id1);
        let root2 = self.find(id2);
        if root1 == root2 {
            return (root1, false);
        }

        let r1 = self.ranks[root1.as_usize()];
        let r2 = self.ranks[root2.as_usize()];

        let (winner, loser) = if r1 < r2 {
            (root2, root1)
        } else if r1 > r2 {
            (root1, root2)
        } else {
            self.ranks[root1.as_usize()] = r1.saturating_add(1);
            (root1, root2)
        };

        self.parents[loser.as_usize()] = winner;
        (winner, true)
    }

    /// Returns the number of elements allocated in the union-find.
    #[inline(always)]
    #[must_use]
    pub fn len(&self) -> usize {
        self.parents.len()
    }

    /// Returns whether the union-find is empty.
    #[inline(always)]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.parents.is_empty()
    }
}
