//! Permutation groups: base and strong generating sets, backtrack search, and
//! the structural algorithms built on them.
//!
//! A permutation of `0..n` is the vector of its images; the product is
//! composition, `(p q)(i) = p(q(i))`. A [`Group`] is a *stabiliser chain*
//! (a base `b₀ … b_{k-1}` and, per level, the generators, the basic orbit and
//! a transversal of the stabiliser `G⁽ⁱ⁾` of `b₀ … b_{i-1}`), built by the
//! incremental deterministic Schreier–Sims algorithm. Everything else is
//! derived from it:
//!
//! * **order, membership, uniform random elements** come straight from the
//!   chain (an element is a product of one transversal element per level);
//! * **one backtrack engine** ([`Group::subgroup`]) finds the subgroup of all
//!   elements with a property, walking the tree of base images with a
//!   property-specific pruning and testing only one image per orbit of the
//!   part of the subgroup already found. Centralisers, normalisers, set
//!   stabilisers, centres and intersections are all instances; the
//!   conjugacy test ([`Group::conjugating_element`]) is the same walk looking
//!   for one element of a coset;
//! * **conjugacy classes** by random elements and their powers, each new
//!   class verified by a conjugacy test and sized by its centraliser, until
//!   the class equation sums to the order;
//! * **Sylow subgroups** by growing a p-subgroup inside its normaliser;
//! * **block systems** (Atkinson's union–find), primitivity, and the
//!   composition factors of solvable groups.
//!
//! Sizes are exact (`BigUint` orders). The algorithms are the textbook ones
//! (Seress, *Permutation Group Algorithms*; Holt, Eick & O'Brien, *Handbook of
//! Computational Group Theory*), sized for degrees up to a few hundred and
//! orders that a backtrack can still walk.

use num_bigint::BigUint;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;
use rand_v09::Rng;
use rand_v09::SeedableRng;
use rand_v09::rngs::StdRng;

/// A permutation of `0..n`: `p[i]` is the image of `i`.
pub type Perm = Vec<u32>;

/// The identity of degree `n`.
#[must_use]
pub fn identity(n: usize) -> Perm {
    (0..n)
        .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
        .collect()
}

/// `p ∘ q`: first `q`, then `p`.
#[must_use]
pub fn compose(
    p: &[u32],
    q: &[u32],
) -> Perm {
    q.iter()
        .map(|&i| p.get(i as usize).copied().unwrap_or(i))
        .collect()
}

/// The inverse permutation.
#[must_use]
pub fn inverse(p: &[u32]) -> Perm {
    let mut inv = vec![0; p.len()];
    for (i, &image) in p.iter().enumerate() {
        if let Some(slot) = inv.get_mut(image as usize) {
            *slot = u32::try_from(i).unwrap_or(u32::MAX);
        }
    }
    inv
}

/// Whether `p` fixes every point.
#[must_use]
pub fn is_identity(p: &[u32]) -> bool {
    p.iter().enumerate().all(|(i, &x)| x as usize == i)
}

/// Whether `p` is a permutation of `0..p.len()`.
#[must_use]
pub fn is_permutation(p: &[u32]) -> bool {
    let mut seen = vec![false; p.len()];
    p.iter().all(|&x| {
        seen.get_mut(x as usize)
            .is_some_and(|s| !std::mem::replace(s, true))
    })
}

/// The cycles of `p` (fixed points included), each starting at its least point.
#[must_use]
pub fn cycles(p: &[u32]) -> Vec<Vec<u32>> {
    let mut seen = vec![false; p.len()];
    let mut out = Vec::new();
    for start in 0..p.len() {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut x = start;
        while !seen[x] {
            seen[x] = true;
            cycle.push(u32::try_from(x).unwrap_or(u32::MAX));
            x = p[x] as usize;
        }
        out.push(cycle);
    }
    out
}

/// The cycle type: the sorted cycle lengths, fixed points included.
#[must_use]
pub fn cycle_type(p: &[u32]) -> Vec<usize> {
    let mut t: Vec<usize> = cycles(p).iter().map(Vec::len).collect();
    t.sort_unstable();
    t
}

/// For each point, the length of its cycle under `p`.
fn cycle_lengths(p: &[u32]) -> Vec<usize> {
    let mut out = vec![0; p.len()];
    for c in cycles(p) {
        for &x in &c {
            out[x as usize] = c.len();
        }
    }
    out
}

/// The order of `p` (the lcm of its cycle lengths).
#[must_use]
pub fn element_order(p: &[u32]) -> BigUint {
    cycles(p).iter().fold(BigUint::one(), |acc, c| {
        let len = BigUint::from(c.len());
        let g = num_integer::Integer::gcd(&acc, &len);
        acc * len / g
    })
}

/// `p^k` for `k ≥ 0` by repeated squaring.
#[must_use]
pub fn power(
    p: &[u32],
    mut k: u64,
) -> Perm {
    let mut base = p.to_vec();
    let mut acc = identity(p.len());
    while k > 0 {
        if k & 1 == 1 {
            acc = compose(&acc, &base);
        }
        base = compose(&base, &base);
        k >>= 1;
    }
    acc
}

/// `g p g⁻¹`.
#[must_use]
pub fn conjugate(
    p: &[u32],
    g: &[u32],
) -> Perm {
    compose(&compose(g, p), &inverse(g))
}

/// One level of a stabiliser chain.
#[derive(Clone, Debug)]
struct Level {
    point: u32,
    /// Generators of `G⁽ⁱ⁾` (they fix the earlier base points).
    gens: Vec<Perm>,
    /// The basic orbit, in discovery order.
    orbit: Vec<u32>,
    /// `transversal[β]`: an element of `G⁽ⁱ⁾` mapping `point` to `β`.
    transversal: Vec<Option<Perm>>,
}

/// A permutation group as a stabiliser chain.
#[derive(Clone, Debug)]
pub struct Group {
    n: usize,
    levels: Vec<Level>,
    /// The generators it was built from.
    gens: Vec<Perm>,
}

impl Group {
    /// The group generated by `gens`, all of degree `n`.
    ///
    /// Generators of another degree are ignored; an empty list gives the
    /// trivial group.
    #[must_use]
    pub fn new(
        n: usize,
        gens: &[Perm],
    ) -> Self {
        let mut g = Self {
            n,
            levels: Vec::new(),
            gens: Vec::new(),
        };
        for p in gens {
            if p.len() == n && is_permutation(p) && !is_identity(p) {
                g.gens.push(p.clone());
                g.insert(p.clone(), 0);
            }
        }
        g
    }

    /// Adds a generator, keeping the chain (incremental Schreier–Sims).
    pub fn extend(
        &mut self,
        p: Perm,
    ) {
        if p.len() == self.n && is_permutation(&p) && !self.contains(&p) {
            self.gens.push(p.clone());
            self.insert(p, 0);
        }
    }

    /// The trivial group of degree `n`.
    #[must_use]
    pub fn trivial(n: usize) -> Self {
        Self::new(n, &[])
    }

    /// The degree.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.n
    }

    /// The generators the group was built from.
    #[must_use]
    pub fn generators(&self) -> &[Perm] {
        &self.gens
    }

    /// The base.
    #[must_use]
    pub fn base(&self) -> Vec<u32> {
        self.levels.iter().map(|l| l.point).collect()
    }

    /// A strong generating set (the union of the level generators).
    #[must_use]
    pub fn strong_generators(&self) -> Vec<Perm> {
        let mut out: Vec<Perm> = Vec::new();
        for l in &self.levels {
            for g in &l.gens {
                if !out.contains(g) {
                    out.push(g.clone());
                }
            }
        }
        out
    }

    /// The lengths of the basic orbits.
    #[must_use]
    pub fn basic_orbit_lengths(&self) -> Vec<usize> {
        self.levels.iter().map(|l| l.orbit.len()).collect()
    }

    /// The order: the product of the basic orbit lengths.
    #[must_use]
    pub fn order(&self) -> BigUint {
        self.levels
            .iter()
            .fold(BigUint::one(), |acc, l| acc * BigUint::from(l.orbit.len()))
    }

    fn new_level(
        &self,
        point: u32,
    ) -> Level {
        let mut transversal = vec![None; self.n];
        transversal[point as usize] = Some(identity(self.n));
        Level {
            point,
            gens: Vec::new(),
            orbit: vec![point],
            transversal,
        }
    }

    /// Sifts `g` through levels `from..`: the residue and the level at which
    /// it stopped (`levels.len()` when it went through).
    fn strip(
        &self,
        mut g: Perm,
        from: usize,
    ) -> (Perm, usize) {
        for (i, level) in self.levels.iter().enumerate().skip(from) {
            let beta = g[level.point as usize] as usize;
            match &level.transversal[beta] {
                | Some(u) => g = compose(&inverse(u), &g),
                | None => return (g, i),
            }
        }
        let len = self.levels.len();
        (g, len)
    }

    /// Adds `g` (which fixes the base points before `level`) and restores
    /// the chain invariant: every Schreier generator sifts.
    fn insert(
        &mut self,
        g: Perm,
        level: usize,
    ) {
        let (h, j) = self.strip(g, level);
        if is_identity(&h) {
            return;
        }
        if j == self.levels.len() {
            let Some(moved) = h
                .iter()
                .enumerate()
                .find(|&(i, &x)| x as usize != i)
                .map(|(i, _)| i)
            else {
                return;
            };
            let level = self.new_level(u32::try_from(moved).unwrap_or(u32::MAX));
            self.levels.push(level);
        }
        // `h` fixes the base points before `j`, so it belongs to every level
        // from `level` to `j`; add it from the deepest one up.
        for l in (level..=j).rev() {
            self.add_generator(l, h.clone());
        }
    }

    /// Adds a generator to one level, extends the orbit and sifts the new
    /// Schreier generators into the next level.
    fn add_generator(
        &mut self,
        l: usize,
        h: Perm,
    ) {
        let Some(level) = self.levels.get_mut(l) else {
            return;
        };
        level.gens.push(h.clone());
        // Extend the orbit with every generator; remember the new points.
        let old_len = level.orbit.len();
        let mut k = 0;
        while k < level.orbit.len() {
            let delta = level.orbit[k] as usize;
            for s in &level.gens {
                let image = s[delta] as usize;
                if level.transversal[image].is_none() {
                    let u = level.transversal[delta]
                        .clone()
                        .unwrap_or_else(|| identity(s.len()));
                    level.transversal[image] = Some(compose(s, &u));
                    level.orbit.push(u32::try_from(image).unwrap_or(u32::MAX));
                }
            }
            k += 1;
        }
        // New Schreier generators: `h` with every old point, every
        // generator with every new point.
        let mut schreier: Vec<Perm> = Vec::new();
        let gens = level.gens.clone();
        let orbit = level.orbit.clone();
        let transversal = &level.transversal;
        let schreier_gen = |s: &Perm, delta: usize| -> Option<Perm> {
            let u = transversal[delta].as_ref()?;
            let su = compose(s, u);
            let image = su[level.point as usize] as usize;
            let v = transversal[image].as_ref()?;
            let g = compose(&inverse(v), &su);
            (!is_identity(&g)).then_some(g)
        };
        for &delta in orbit.iter().take(old_len) {
            if let Some(g) = schreier_gen(&h, delta as usize) {
                schreier.push(g);
            }
        }
        for &delta in orbit.iter().skip(old_len) {
            for s in &gens {
                if let Some(g) = schreier_gen(s, delta as usize) {
                    schreier.push(g);
                }
            }
        }
        for g in schreier {
            self.insert(g, l + 1);
        }
    }

    /// Whether `p` belongs to the group.
    #[must_use]
    pub fn contains(
        &self,
        p: &[u32],
    ) -> bool {
        if p.len() != self.n || !is_permutation(p) {
            return false;
        }
        let (h, j) = self.strip(p.to_vec(), 0);
        j == self.levels.len() && is_identity(&h)
    }

    /// A uniformly random element.
    pub fn random_element(
        &self,
        rng: &mut StdRng,
    ) -> Perm {
        let mut g = identity(self.n);
        for level in &self.levels {
            let beta = level.orbit[rng.random_range(0..level.orbit.len())] as usize;
            if let Some(u) = &level.transversal[beta] {
                g = compose(&g, u);
            }
        }
        g
    }

    /// Every element (only for small groups; `None` above `limit`).
    #[must_use]
    pub fn elements(
        &self,
        limit: usize,
    ) -> Option<Vec<Perm>> {
        if self.order() > BigUint::from(limit) {
            return None;
        }
        let mut out = vec![identity(self.n)];
        for level in self.levels.iter().rev() {
            let mut next = Vec::with_capacity(out.len() * level.orbit.len());
            for &beta in &level.orbit {
                let u = level.transversal[beta as usize].as_ref()?;
                for g in &out {
                    next.push(compose(u, g));
                }
            }
            out = next;
        }
        Some(out)
    }

    /// The orbits of the group on `0..n`.
    #[must_use]
    pub fn orbits(&self) -> Vec<Vec<u32>> {
        orbits(self.n, &self.gens)
    }

    /// Whether the group acts transitively.
    #[must_use]
    pub fn is_transitive(&self) -> bool {
        self.orbits().len() <= 1
    }

    /// Whether every pair of generators commutes.
    #[must_use]
    pub fn is_abelian(&self) -> bool {
        self.gens.iter().enumerate().all(|(i, a)| {
            self.gens
                .iter()
                .skip(i + 1)
                .all(|b| compose(a, b) == compose(b, a))
        })
    }

    /// Whether `other` is a subgroup.
    #[must_use]
    pub fn contains_group(
        &self,
        other: &Self,
    ) -> bool {
        other.gens.iter().all(|g| self.contains(g))
    }

    // -----------------------------------------------------------------
    // Backtrack search
    // -----------------------------------------------------------------

    /// The subgroup `{g ∈ G : property(g)}`, which must be a subgroup.
    ///
    /// `prune(level, point, image)` may reject a partial element as soon as
    /// it maps the base point `point` to `image`: return `false` when no
    /// element of the subgroup can do that. Pruning must never reject a
    /// subgroup element; it only makes the search faster.
    pub fn subgroup(
        &self,
        property: &dyn Fn(&Perm) -> bool,
        prune: &dyn Fn(u32, u32) -> bool,
    ) -> Self {
        let mut found = Self::trivial(self.n);
        let k = self.levels.len();
        for i in (0..k).rev() {
            let level = &self.levels[i];
            // The part of the subgroup fixing b₀ … b_{i-1}: its orbits on
            // the basic orbit are tested once each.
            let fixed: Vec<u32> = self.levels[..i].iter().map(|l| l.point).collect();
            let stab = found.pointwise_stabilizer(&fixed);
            let mut covered = vec![false; self.n];
            for o in orbits(self.n, &stab.gens) {
                if o.contains(&level.point) {
                    for &x in &o {
                        covered[x as usize] = true;
                    }
                }
            }
            for &beta in &level.orbit {
                if covered[beta as usize] || !prune(level.point, beta) {
                    continue;
                }
                let Some(u) = &level.transversal[beta as usize] else {
                    continue;
                };
                if let Some(g) = self.search(i + 1, u.clone(), property, prune) {
                    found = Self::new(self.n, &[found.gens.clone(), vec![g]].concat());
                    // Recompute the orbits of the grown stabiliser.
                    let stab = found.pointwise_stabilizer(&fixed);
                    for o in orbits(self.n, &stab.gens) {
                        if o.contains(&level.point) {
                            for &x in &o {
                                covered[x as usize] = true;
                            }
                        }
                    }
                }
            }
        }
        found
    }

    /// Depth-first search for an element `prefix ∘ h` (`h ∈ G⁽ʲ⁾`) with the
    /// property.
    fn search(
        &self,
        j: usize,
        prefix: Perm,
        property: &dyn Fn(&Perm) -> bool,
        prune: &dyn Fn(u32, u32) -> bool,
    ) -> Option<Perm> {
        let Some(level) = self.levels.get(j) else {
            return property(&prefix).then_some(prefix);
        };
        for &gamma in &level.orbit {
            let Some(u) = &level.transversal[gamma as usize] else {
                continue;
            };
            let candidate = compose(&prefix, u);
            if !prune(level.point, candidate[level.point as usize]) {
                continue;
            }
            if let Some(g) = self.search(j + 1, candidate, property, prune) {
                return Some(g);
            }
        }
        None
    }

    /// The pointwise stabiliser of `points`.
    #[must_use]
    pub fn pointwise_stabilizer(
        &self,
        points: &[u32],
    ) -> Self {
        if points.is_empty() {
            return self.clone();
        }
        // With the points as a base prefix the stabiliser is a tail of the
        // chain: rebuild with a base starting at `points`.
        let rebased = self.with_base_prefix(points);
        let depth = points.len();
        let gens: Vec<Perm> = rebased
            .levels
            .iter()
            .skip(depth)
            .flat_map(|l| l.gens.iter().cloned())
            .filter(|g| points.iter().all(|&p| g[p as usize] == p))
            .collect();
        Self::new(self.n, &gens)
    }

    /// The same group with a chain whose base starts with `prefix`.
    #[must_use]
    pub fn with_base_prefix(
        &self,
        prefix: &[u32],
    ) -> Self {
        let mut g = Self {
            n: self.n,
            levels: Vec::new(),
            gens: Vec::new(),
        };
        for &p in prefix {
            if (p as usize) < self.n && !g.levels.iter().any(|l| l.point == p) {
                let level = g.new_level(p);
                g.levels.push(level);
            }
        }
        for s in self.strong_generators() {
            g.gens.push(s.clone());
            g.insert(s, 0);
        }
        g.gens = self.gens.clone();
        // Drop trivial levels of the prefix only at the end of the chain.
        while g
            .levels
            .last()
            .is_some_and(|l| l.orbit.len() == 1 && l.gens.is_empty())
        {
            g.levels.pop();
        }
        g
    }

    /// The centraliser of `x` in the group.
    #[must_use]
    pub fn centralizer(
        &self,
        x: &[u32],
    ) -> Self {
        let lengths = cycle_lengths(x);
        let x = x.to_vec();
        // A base through the support of `x` prunes hardest.
        let support: Vec<u32> = (0..self.n)
            .filter(|&i| x[i] as usize != i)
            .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
            .collect();
        let g = self.with_base_prefix(&support);
        g.subgroup(&|p| compose(p, &x) == compose(&x, p), &|point, image| {
            lengths[point as usize] == lengths[image as usize]
        })
    }

    /// The centre.
    #[must_use]
    pub fn center(&self) -> Self {
        let gens = self.gens.clone();
        let lengths: Vec<Vec<usize>> = gens.iter().map(|g| cycle_lengths(g)).collect();
        self.subgroup(
            &|p| gens.iter().all(|g| compose(p, g) == compose(g, p)),
            &|point, image| {
                lengths
                    .iter()
                    .all(|l| l[point as usize] == l[image as usize])
            },
        )
    }

    /// The stabiliser of the set `points`.
    #[must_use]
    pub fn set_stabilizer(
        &self,
        points: &[u32],
    ) -> Self {
        let mut inside = vec![false; self.n];
        for &p in points {
            if let Some(f) = inside.get_mut(p as usize) {
                *f = true;
            }
        }
        let g = self.with_base_prefix(points);
        let check = inside.clone();
        g.subgroup(
            &|p| (0..self.n).all(|i| check[i] == check[p[i] as usize]),
            &|point, image| inside[point as usize] == inside[image as usize],
        )
    }

    /// The normaliser of the subgroup `h` (generated inside this group's
    /// degree; it need not be a subgroup).
    #[must_use]
    pub fn normalizer(
        &self,
        h: &Self,
    ) -> Self {
        let orbit_size = {
            let mut size = vec![1; self.n];
            for o in h.orbits() {
                for &x in &o {
                    size[x as usize] = o.len();
                }
            }
            size
        };
        let hg = h.gens.clone();
        self.subgroup(
            &|p| hg.iter().all(|x| h.contains(&conjugate(x, p))),
            &|point, image| orbit_size[point as usize] == orbit_size[image as usize],
        )
    }

    /// The intersection with another group of the same degree.
    #[must_use]
    pub fn intersection(
        &self,
        other: &Self,
    ) -> Self {
        let other_orbit = {
            let mut label = vec![usize::MAX; self.n];
            for (k, o) in other.orbits().iter().enumerate() {
                for &x in o {
                    label[x as usize] = k;
                }
            }
            label
        };
        self.subgroup(&|p| other.contains(p), &|point, image| {
            other_orbit[point as usize] == other_orbit[image as usize]
        })
    }

    /// An element `g` with `g x g⁻¹ = y`, if `x` and `y` are conjugate.
    #[must_use]
    pub fn conjugating_element(
        &self,
        x: &[u32],
        y: &[u32],
    ) -> Option<Perm> {
        if cycle_type(x) != cycle_type(y) {
            return None;
        }
        let (lx, ly) = (cycle_lengths(x), cycle_lengths(y));
        let (x, y) = (x.to_vec(), y.to_vec());
        let support: Vec<u32> = (0..self.n)
            .filter(|&i| x[i] as usize != i)
            .map(|i| u32::try_from(i).unwrap_or(u32::MAX))
            .collect();
        let g = self.with_base_prefix(&support);
        let property = |p: &Perm| compose(p, &x) == compose(&y, p);
        let prune = |point: u32, image: u32| lx[point as usize] == ly[image as usize];
        g.search(0, identity(self.n), &property, &prune)
    }

    /// Whether `x` and `y` are conjugate in the group.
    #[must_use]
    pub fn are_conjugate(
        &self,
        x: &[u32],
        y: &[u32],
    ) -> bool {
        self.conjugating_element(x, y).is_some()
    }

    /// The normal closure of `x`'s group in this group.
    #[must_use]
    pub fn normal_closure(
        &self,
        gens: &[Perm],
    ) -> Self {
        let mut closure = Self::new(self.n, gens);
        let mut queue: Vec<Perm> = gens.to_vec();
        while let Some(x) = queue.pop() {
            for g in &self.gens {
                let c = conjugate(&x, g);
                if !closure.contains(&c) {
                    closure = Self::new(self.n, &[closure.gens.clone(), vec![c.clone()]].concat());
                    queue.push(c);
                }
            }
        }
        closure
    }

    /// Whether `h` is a normal subgroup.
    #[must_use]
    pub fn is_normal(
        &self,
        h: &Self,
    ) -> bool {
        self.contains_group(h)
            && h.gens
                .iter()
                .all(|x| self.gens.iter().all(|g| h.contains(&conjugate(x, g))))
    }

    /// The derived subgroup `[G, G]`.
    #[must_use]
    pub fn derived_subgroup(&self) -> Self {
        let mut comms = Vec::new();
        for a in &self.gens {
            for b in &self.gens {
                let c = compose(&compose(&inverse(a), &inverse(b)), &compose(a, b));
                if !is_identity(&c) {
                    comms.push(c);
                }
            }
        }
        self.normal_closure(&comms)
    }

    /// Whether the group is solvable.
    #[must_use]
    pub fn is_solvable(&self) -> bool {
        let mut g = self.clone();
        loop {
            let d = g.derived_subgroup();
            if d.order().is_one() {
                return true;
            }
            if d.order() == g.order() {
                return false;
            }
            g = d;
        }
    }

    // -----------------------------------------------------------------
    // Conjugacy classes, Sylow subgroups
    // -----------------------------------------------------------------

    /// The conjugacy classes as `(representative, size)`, sorted by element
    /// order and then size; `None` when random search did not complete the
    /// class equation within `max_tries` samples.
    ///
    /// Random elements and all their powers are sampled; a sample in no known
    /// class (cycle type first, then a conjugacy test) opens a class of size
    /// `|G| / |C(g)|`. The result is exact: it is returned only when the
    /// class sizes sum to the order.
    #[must_use]
    pub fn conjugacy_classes(
        &self,
        seed: u64,
        max_tries: usize,
    ) -> Option<Vec<(Perm, BigUint)>> {
        let order = self.order();
        let mut rng = StdRng::seed_from_u64(seed);
        let mut classes: Vec<(Perm, BigUint, Vec<usize>)> = Vec::new();
        let mut total = BigUint::zero();
        let consider =
            |g: Perm, classes: &mut Vec<(Perm, BigUint, Vec<usize>)>, total: &mut BigUint| {
                let ty = cycle_type(&g);
                let known = classes
                    .iter()
                    .filter(|(_, _, t)| *t == ty)
                    .any(|(rep, _, _)| self.are_conjugate(rep, &g));
                if !known {
                    let size = &order / self.centralizer(&g).order();
                    *total += &size;
                    classes.push((g, size, ty));
                }
            };
        consider(identity(self.n), &mut classes, &mut total);
        for g in self.gens.clone() {
            if total == order {
                break;
            }
            consider(g, &mut classes, &mut total);
        }
        let mut tries = 0;
        while total < order && tries < max_tries {
            tries += 1;
            let g = self.random_element(&mut rng);
            // Every power of a sample: classes of small elements are rare as
            // random samples but common as powers.
            let ord = element_order(&g).to_u64().unwrap_or(1).min(1024);
            let mut x = g.clone();
            for _ in 1..ord {
                if total == order {
                    break;
                }
                consider(x.clone(), &mut classes, &mut total);
                x = compose(&x, &g);
            }
        }
        if total != order {
            return None;
        }
        let mut out: Vec<(Perm, BigUint)> = classes.into_iter().map(|(g, s, _)| (g, s)).collect();
        out.sort_by(|a, b| {
            element_order(&a.0)
                .cmp(&element_order(&b.0))
                .then(a.1.cmp(&b.1))
        });
        Some(out)
    }

    /// A Sylow `p`-subgroup.
    ///
    /// Grows a `p`-subgroup `P` by `p`-elements of its normaliser outside it
    /// until `|P|` is the full power of `p` dividing `|G|` (a `p`-subgroup is
    /// Sylow exactly when `p` does not divide `[N(P) : P]`, so a `p`-element
    /// of `N(P) \ P` exists until then). `None` if random sampling failed
    /// `max_tries` times in a row.
    #[must_use]
    pub fn sylow_subgroup(
        &self,
        p: u64,
        seed: u64,
        max_tries: usize,
    ) -> Option<Self> {
        let order = self.order();
        let pb = BigUint::from(p);
        let mut target = BigUint::one();
        let mut rest = order.clone();
        while (&rest % &pb).is_zero() {
            rest /= &pb;
            target *= &pb;
        }
        let mut sylow = Self::trivial(self.n);
        let mut rng = StdRng::seed_from_u64(seed);
        while sylow.order() < target {
            let normalizer = if sylow.order().is_one() {
                self.clone()
            } else {
                self.normalizer(&sylow)
            };
            let mut grown = false;
            for _ in 0..max_tries {
                let g = normalizer.random_element(&mut rng);
                let ord = element_order(&g);
                if !(&ord % &pb).is_zero() {
                    continue;
                }
                // The p-part of g: g^(ord / p^a).
                let mut cofactor = ord.clone();
                while (&cofactor % &pb).is_zero() {
                    cofactor /= &pb;
                }
                let y = power(&g, cofactor.to_u64().unwrap_or(1));
                if !sylow.contains(&y) {
                    sylow = Self::new(self.n, &[sylow.gens.clone(), vec![y]].concat());
                    grown = true;
                    break;
                }
            }
            if !grown {
                return None;
            }
        }
        Some(sylow)
    }

    /// The exponent: the lcm of the element orders (over class
    /// representatives).
    #[must_use]
    pub fn exponent(
        &self,
        seed: u64,
    ) -> Option<BigUint> {
        let classes = self.conjugacy_classes(seed, 100_000)?;
        Some(classes.iter().fold(BigUint::one(), |acc, (g, _)| {
            let o = element_order(g);
            let d = num_integer::Integer::gcd(&acc, &o);
            acc * o / d
        }))
    }

    /// The composition factors of a solvable group, as the primes `p` of the
    /// cyclic factors `C_p` (with multiplicity, top of the series first);
    /// `None` for a non-solvable group.
    #[must_use]
    pub fn solvable_composition_factors(&self) -> Option<Vec<u64>> {
        let mut out = Vec::new();
        let mut g = self.clone();
        loop {
            let order = g.order();
            if order.is_one() {
                return Some(out);
            }
            let d = g.derived_subgroup();
            if d.order() == order {
                return None;
            }
            let mut index = order / d.order();
            let mut q = 2_u64;
            while !index.is_one() {
                if BigUint::from(q) * BigUint::from(q) > index {
                    out.push(index.to_u64()?);
                    break;
                }
                while (&index % q).is_zero() {
                    out.push(q);
                    index /= q;
                }
                q += 1;
            }
            g = d;
        }
    }
}

/// The orbits of the group generated by `gens` on `0..n`.
#[must_use]
pub fn orbits(
    n: usize,
    gens: &[Perm],
) -> Vec<Vec<u32>> {
    let mut seen = vec![false; n];
    let mut out = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut orbit = vec![u32::try_from(start).unwrap_or(u32::MAX)];
        let mut k = 0;
        while k < orbit.len() {
            let x = orbit[k] as usize;
            for g in gens {
                let y = g.get(x).copied().unwrap_or(0) as usize;
                if y < n && !seen[y] {
                    seen[y] = true;
                    orbit.push(u32::try_from(y).unwrap_or(u32::MAX));
                }
            }
            k += 1;
        }
        out.push(orbit);
    }
    out
}

// ---------------------------------------------------------------------------
// Blocks
// ---------------------------------------------------------------------------

fn find(
    parent: &mut [usize],
    mut x: usize,
) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// The finest block system of a transitive group in which `a` and `b` lie
/// in one block (Atkinson's algorithm): the blocks, each sorted, in order of
/// their least point.
#[must_use]
pub fn minimal_blocks(
    n: usize,
    gens: &[Perm],
    a: u32,
    b: u32,
) -> Vec<Vec<u32>> {
    let mut parent: Vec<usize> = (0..n).collect();
    let mut queue = vec![(a as usize, b as usize)];
    if let (Some(ra), Some(rb)) = (
        Some(find(&mut parent, a as usize)),
        Some(find(&mut parent, b as usize)),
    ) {
        if ra != rb {
            parent[ra.max(rb)] = ra.min(rb);
        }
    }
    while let Some((x, y)) = queue.pop() {
        for g in gens {
            let (gx, gy) = (g[x] as usize, g[y] as usize);
            let (rx, ry) = (find(&mut parent, gx), find(&mut parent, gy));
            if rx != ry {
                parent[rx.max(ry)] = rx.min(ry);
                queue.push((gx, gy));
            }
        }
    }
    let mut blocks: std::collections::BTreeMap<usize, Vec<u32>> = std::collections::BTreeMap::new();
    for x in 0..n {
        let r = find(&mut parent, x);
        blocks
            .entry(r)
            .or_default()
            .push(u32::try_from(x).unwrap_or(u32::MAX));
    }
    blocks.into_values().collect()
}

/// Whether a transitive group is primitive (no block system other than the
/// trivial ones). Intransitive groups are not primitive.
#[must_use]
pub fn is_primitive(
    n: usize,
    gens: &[Perm],
) -> bool {
    if orbits(n, gens).len() > 1 {
        return false;
    }
    (1..n).all(|b| minimal_blocks(n, gens, 0, u32::try_from(b).unwrap_or(u32::MAX)).len() == 1)
}

/// The minimal non-trivial block systems of a transitive group (each a
/// list of blocks), without repetitions.
#[must_use]
pub fn minimal_block_systems(
    n: usize,
    gens: &[Perm],
) -> Vec<Vec<Vec<u32>>> {
    let mut systems: Vec<Vec<Vec<u32>>> = Vec::new();
    for b in 1..n {
        let s = minimal_blocks(n, gens, 0, u32::try_from(b).unwrap_or(u32::MAX));
        if s.len() > 1 && !systems.contains(&s) {
            systems.push(s);
        }
    }
    // Keep the minimal ones: a system whose block of 0 contains another
    // system's block of 0 is coarser.
    let block0 = |s: &Vec<Vec<u32>>| {
        s.iter()
            .find(|b| b.contains(&0))
            .cloned()
            .unwrap_or_default()
    };
    let all = systems.clone();
    systems.retain(|s| {
        let mine = block0(s);
        !all.iter().any(|t| {
            let theirs = block0(t);
            theirs.len() < mine.len() && theirs.iter().all(|x| mine.contains(x))
        })
    });
    systems
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cyc(n: usize) -> Perm {
        (0..n).map(|i| ((i + 1) % n) as u32).collect()
    }

    fn swap01(n: usize) -> Perm {
        let mut p = identity(n);
        p.swap(0, 1);
        p
    }

    fn sym(n: usize) -> Group {
        Group::new(n, &[cyc(n), swap01(n)])
    }

    fn factorial(n: u64) -> BigUint {
        (1..=n).map(BigUint::from).product()
    }

    #[test]
    fn orders_membership_and_mathieu() {
        for n in 2..=9 {
            assert_eq!(sym(n).order(), factorial(n as u64), "S{n}");
        }
        assert_eq!(sym(12).order(), factorial(12));
        // M11 on 11 points (0-based from the classical generators).
        let a: Perm = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0];
        let b: Perm = vec![0, 1, 6, 9, 5, 3, 10, 2, 8, 4, 7];
        let m11 = Group::new(11, &[a.clone(), b]);
        assert_eq!(m11.order(), BigUint::from(7920_u32));
        assert!(m11.contains(&compose(&a, &a)));
        assert!(!m11.contains(&swap01(11)));
        assert!(is_primitive(11, m11.generators()));
    }

    #[test]
    fn centralizers_centres_and_conjugacy() {
        let s5 = sym(5);
        // C(5-cycle) in S5 is the cyclic group of order 5.
        assert_eq!(s5.centralizer(&cyc(5)).order(), BigUint::from(5_u32));
        // C((0 1)) in S5 = <(0 1)> x S3, order 12.
        assert_eq!(s5.centralizer(&swap01(5)).order(), BigUint::from(12_u32));
        assert!(s5.center().order().is_one());
        // Dihedral group of order 8 on 4 points: centre of order 2.
        let d4 = Group::new(4, &[cyc(4), vec![0, 3, 2, 1]]);
        assert_eq!(d4.order(), BigUint::from(8_u32));
        assert_eq!(d4.center().order(), BigUint::from(2_u32));
        let x: Perm = vec![1, 0, 2, 3, 4];
        let y: Perm = vec![0, 1, 2, 4, 3];
        let g = s5.conjugating_element(&x, &y).unwrap_or_default();
        assert_eq!(conjugate(&x, &g), y);
        assert!(!s5.are_conjugate(&x, &cyc(5)));
    }

    #[test]
    fn class_equations() {
        // S5 has 7 classes; A5 has 5 (the 5-cycles split).
        let s5 = sym(5);
        let classes = s5.conjugacy_classes(1, 10_000).unwrap_or_default();
        assert_eq!(classes.len(), 7);
        let a5 = s5.derived_subgroup();
        assert_eq!(a5.order(), BigUint::from(60_u32));
        let classes = a5.conjugacy_classes(2, 10_000).unwrap_or_default();
        assert_eq!(classes.len(), 5);
        let sizes: Vec<u64> = classes
            .iter()
            .map(|(_, s)| s.to_u64().unwrap_or(0))
            .collect();
        assert_eq!(sizes.iter().sum::<u64>(), 60);
        assert!(sizes.contains(&12) && sizes.contains(&15) && sizes.contains(&20));
        assert_eq!(a5.exponent(3), Some(BigUint::from(30_u32)));
    }

    #[test]
    fn sylow_subgroups_normalizers_and_intersections() {
        let s6 = sym(6);
        let p2 = s6
            .sylow_subgroup(2, 7, 1000)
            .unwrap_or_else(|| Group::trivial(6));
        assert_eq!(p2.order(), BigUint::from(16_u32));
        let p3 = s6
            .sylow_subgroup(3, 7, 1000)
            .unwrap_or_else(|| Group::trivial(6));
        assert_eq!(p3.order(), BigUint::from(9_u32));
        let p5 = s6
            .sylow_subgroup(5, 7, 1000)
            .unwrap_or_else(|| Group::trivial(6));
        assert_eq!(p5.order(), BigUint::from(5_u32));
        // N(P5) in S6 has order 20 (6 Sylow 5-subgroups).
        assert_eq!(s6.normalizer(&p5).order(), BigUint::from(20_u32));
        let a6 = s6.derived_subgroup();
        assert!(s6.is_normal(&a6));
        assert_eq!(a6.intersection(&p2).order(), BigUint::from(8_u32));
        // Set stabiliser of {0, 1} in S6: S2 x S4.
        assert_eq!(s6.set_stabilizer(&[0, 1]).order(), BigUint::from(48_u32));
    }

    #[test]
    fn blocks_and_composition_factors() {
        // D4 on the square's vertices: {0, 2} and {1, 3} form a block system.
        let d4 = [cyc(4), vec![0, 3, 2, 1]];
        assert_eq!(minimal_blocks(4, &d4, 0, 2), vec![vec![0, 2], vec![1, 3]]);
        assert!(!is_primitive(4, &d4));
        assert!(is_primitive(5, sym(5).generators()));
        // S4 is solvable: composition factors 2, 3, 2, 2.
        let mut f = sym(4).solvable_composition_factors().unwrap_or_default();
        f.sort_unstable();
        assert_eq!(f, vec![2, 2, 2, 3]);
        assert!(sym(5).solvable_composition_factors().is_none());
        assert!(!sym(5).is_solvable());
    }
}
