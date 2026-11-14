//! The transitive permutation groups of degree 6 and 7, as the possible
//! Galois groups of irreducible polynomials of those degrees.
//!
//! Each group is given by generators acting on the positions `0..n`; the
//! tables list one representative of every conjugacy class in `S_n` (16 of
//! degree 6, 7 of degree 7), ordered by size. [`TransitiveGroup::elements`]
//! enumerates a group by closure, and [`TransitiveGroup::block_systems`]
//! finds its systems of imprimitivity.

/// A transitive subgroup of `S_n`.
#[derive(Debug)]
pub struct TransitiveGroup {
    /// The degree `n`.
    pub degree: usize,
    /// Conventional name.
    pub label: &'static str,
    /// The order of the group.
    pub order: usize,
    /// Whether the group is solvable.
    pub solvable: bool,
    generators: &'static [&'static [u8]],
}

/// A partition of `0..n` into blocks permuted by a group.
pub type BlockSystem = Vec<Vec<usize>>;

const fn g(
    degree: usize,
    label: &'static str,
    order: usize,
    solvable: bool,
    generators: &'static [&'static [u8]],
) -> TransitiveGroup {
    TransitiveGroup {
        degree,
        label,
        order,
        solvable,
        generators,
    }
}

static DEGREE_6: [TransitiveGroup; 16] = [
    g(6, "C6", 6, true, &[&[1, 2, 3, 4, 5, 0]]),
    g(
        6,
        "S3(6)",
        6,
        true,
        &[&[1, 0, 3, 2, 5, 4], &[2, 4, 0, 5, 1, 3]],
    ),
    g(
        6,
        "A4(6)",
        12,
        true,
        &[&[0, 1, 3, 2, 5, 4], &[2, 3, 4, 5, 0, 1]],
    ),
    g(
        6,
        "D6(6)",
        12,
        true,
        &[&[0, 1, 3, 2, 5, 4], &[2, 4, 0, 5, 1, 3]],
    ),
    g(
        6,
        "F18(6)=C3xS3",
        18,
        true,
        &[&[0, 1, 2, 4, 5, 3], &[3, 4, 5, 0, 1, 2]],
    ),
    g(
        6,
        "2A4(6)=C2^3:C3",
        24,
        true,
        &[&[0, 1, 2, 3, 5, 4], &[1, 4, 3, 5, 0, 2]],
    ),
    g(
        6,
        "S4(6d)",
        24,
        true,
        &[&[0, 1, 3, 2, 5, 4], &[2, 4, 1, 5, 0, 3]],
    ),
    g(
        6,
        "S4(6c)",
        24,
        true,
        &[&[0, 1, 3, 4, 5, 2], &[2, 4, 0, 5, 1, 3]],
    ),
    g(
        6,
        "3^2:4",
        36,
        true,
        &[&[0, 1, 2, 4, 5, 3], &[3, 4, 5, 0, 2, 1]],
    ),
    g(
        6,
        "S3xS3",
        36,
        true,
        &[&[0, 1, 3, 2, 5, 4], &[1, 2, 4, 5, 3, 0]],
    ),
    g(
        6,
        "2S4(6)=C2wrS3",
        48,
        true,
        &[&[0, 1, 3, 2, 5, 4], &[2, 4, 1, 3, 0, 5]],
    ),
    g(
        6,
        "A5(6)=PSL(2,5)",
        60,
        false,
        &[&[0, 1, 3, 2, 5, 4], &[1, 2, 0, 4, 5, 3]],
    ),
    g(
        6,
        "3^2:D4=S3wrC2",
        72,
        true,
        &[&[0, 1, 2, 3, 5, 4], &[1, 2, 4, 5, 0, 3]],
    ),
    g(
        6,
        "S5(6)=PGL(2,5)",
        120,
        false,
        &[&[0, 1, 3, 2, 5, 4], &[1, 2, 3, 4, 5, 0]],
    ),
    g(
        6,
        "A6",
        360,
        false,
        &[&[0, 1, 2, 4, 5, 3], &[1, 2, 3, 0, 5, 4]],
    ),
    g(
        6,
        "S6",
        720,
        false,
        &[&[0, 1, 2, 3, 5, 4], &[1, 2, 3, 4, 0, 5]],
    ),
];

static DEGREE_7: [TransitiveGroup; 7] = [
    g(7, "C7", 7, true, &[&[1, 2, 3, 4, 5, 6, 0]]),
    g(
        7,
        "D7",
        14,
        true,
        &[&[1, 2, 3, 4, 5, 6, 0], &[0, 6, 5, 4, 3, 2, 1]],
    ),
    g(
        7,
        "F21",
        21,
        true,
        &[&[1, 2, 3, 4, 5, 6, 0], &[0, 2, 4, 6, 1, 3, 5]],
    ),
    g(
        7,
        "F42",
        42,
        true,
        &[&[1, 2, 3, 4, 5, 6, 0], &[0, 3, 6, 2, 5, 1, 4]],
    ),
    g(
        7,
        "PSL(3,2)",
        168,
        false,
        &[&[0, 1, 4, 3, 2, 6, 5], &[1, 2, 0, 4, 6, 5, 3]],
    ),
    g(
        7,
        "A7",
        2520,
        false,
        &[&[1, 2, 0, 3, 4, 5, 6], &[0, 1, 3, 4, 5, 6, 2]],
    ),
    g(
        7,
        "S7",
        5040,
        false,
        &[&[1, 0, 2, 3, 4, 5, 6], &[1, 2, 3, 4, 5, 6, 0]],
    ),
];

/// The transitive groups of degree `n` (6 or 7), smallest first; empty for
/// other degrees.
#[must_use]
pub fn transitive_groups(n: usize) -> &'static [TransitiveGroup] {
    match n {
        | 6 => &DEGREE_6,
        | 7 => &DEGREE_7,
        | _ => &[],
    }
}

fn compose(
    a: &[u8],
    b: &[u8],
) -> Vec<u8> {
    // (a after b)(i) = a[b[i]]
    b.iter()
        .map(|&i| a.get(usize::from(i)).copied().unwrap_or(0))
        .collect()
}

impl TransitiveGroup {
    /// All elements, as images `[h(0), …, h(n-1)]`.
    #[must_use]
    pub fn elements(&self) -> Vec<Vec<u8>> {
        let identity: Vec<u8> = (0..self.degree)
            .map(|i| u8::try_from(i).unwrap_or(0))
            .collect();
        let mut seen = std::collections::HashSet::new();
        seen.insert(identity.clone());
        let mut frontier = vec![identity.clone()];
        let mut all = vec![identity];
        while !frontier.is_empty() {
            let mut next = Vec::new();
            for element in &frontier {
                for generator in self.generators {
                    let product = compose(generator, element);
                    if seen.insert(product.clone()) {
                        next.push(product.clone());
                        all.push(product);
                    }
                }
            }
            frontier = next;
        }
        all
    }

    /// Whether every element is an even permutation (the discriminant of a
    /// polynomial with this group is then a square).
    #[must_use]
    pub fn is_even(&self) -> bool {
        self.generators.iter().all(|generator| parity(generator))
    }

    /// Whether the group is abelian.
    #[must_use]
    pub fn is_abelian(&self) -> bool {
        self.generators.iter().all(|a| {
            self.generators
                .iter()
                .all(|b| compose(a, b) == compose(b, a))
        })
    }

    /// The non-trivial systems of imprimitivity (partitions of `0..n` into
    /// blocks of equal size permuted by the group), each with its blocks
    /// sorted.
    #[must_use]
    pub fn block_systems(&self) -> Vec<BlockSystem> {
        let n = self.degree;
        let mut systems: Vec<BlockSystem> = Vec::new();
        for j in 1..n {
            // Smallest block system with 0 and j in one block (union-find
            // closed under the generators).
            let mut parent: Vec<usize> = (0..n).collect();
            let mut pending = vec![(0_usize, j)];
            while let Some((a, b)) = pending.pop() {
                let (ra, rb) = (find(&parent, a), find(&parent, b));
                if ra == rb {
                    continue;
                }
                if let Some(slot) = parent.get_mut(ra) {
                    *slot = rb;
                }
                for generator in self.generators {
                    let image = |i: usize| generator.get(i).map_or(0, |&v| usize::from(v));
                    pending.push((image(a), image(b)));
                }
            }
            let roots: Vec<usize> = (0..n).map(|i| find(&parent, i)).collect();
            let mut blocks: Vec<Vec<usize>> = Vec::new();
            let mut keys: Vec<usize> = Vec::new();
            for (i, &root) in roots.iter().enumerate() {
                match keys.iter().position(|&k| k == root) {
                    | Some(at) => {
                        if let Some(block) = blocks.get_mut(at) {
                            block.push(i);
                        }
                    },
                    | None => {
                        keys.push(root);
                        blocks.push(vec![i]);
                    },
                }
            }
            if blocks.len() > 1 && blocks.len() < n && !systems.contains(&blocks) {
                systems.push(blocks);
            }
        }
        systems
    }
}

fn find(
    parent: &[usize],
    i: usize,
) -> usize {
    let mut root = i;
    while parent.get(root).copied().unwrap_or(root) != root {
        root = parent.get(root).copied().unwrap_or(root);
    }
    root
}

fn parity(p: &[u8]) -> bool {
    let mut inversions = 0_usize;
    for (i, a) in p.iter().enumerate() {
        for b in p.iter().skip(i + 1) {
            if a > b {
                inversions += 1;
            }
        }
    }
    inversions.is_multiple_of(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orders_match_the_tables() {
        for n in [6, 7] {
            for group in transitive_groups(n) {
                assert_eq!(group.elements().len(), group.order, "{}", group.label);
            }
        }
        // Solvable exactly when the order is not one of 60, 120, 168, 360, 720, 2520, 5040.
        for group in transitive_groups(6).iter().chain(transitive_groups(7)) {
            let unsolvable = [60, 120, 168, 360, 720, 2520, 5040].contains(&group.order);
            assert_eq!(group.solvable, !unsolvable, "{}", group.label);
        }
    }

    #[test]
    fn block_systems_of_sextic_groups() {
        // C6: blocks of size 2 and of size 3. S6 and A5(6): primitive.
        let c6 = transitive_groups(6)
            .first()
            .map(TransitiveGroup::block_systems)
            .unwrap_or_default();
        assert_eq!(c6.len(), 2);
        let s6 = transitive_groups(6)
            .last()
            .map(TransitiveGroup::block_systems)
            .unwrap_or_default();
        assert!(s6.is_empty(), "not empty");
        for group in transitive_groups(6) {
            assert_eq!(
                group.block_systems().is_empty(),
                !group.solvable,
                "{}",
                group.label
            );
        }
        assert_eq!(
            transitive_groups(6)
                .iter()
                .filter(|g| g.is_abelian())
                .count(),
            1
        );
        assert!(
            transitive_groups(7)
                .iter()
                .all(|g| g.block_systems().is_empty())
        );
    }
}
