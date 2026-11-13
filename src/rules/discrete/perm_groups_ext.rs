//! Structure of permutation groups: subgroup search, classes, Sylow
//! subgroups and blocks, on top of [`kernels::permgroup`](crate::kernels::permgroup).
//!
//! Conventions are those of [`perm_groups`](super::perm_groups): a
//! permutation is the list of the images of `1, ..., n`, a group is a list of
//! generators of equal degree (up to 1024 points here), and a subgroup is
//! returned as a list of generators (the identity alone for the trivial
//! group, so the degree survives).
//!
//! | operator | value |
//! |---|---|
//! | `perm_group_centralizer(gens, x)` | generators of the centraliser of `x` |
//! | `perm_group_center(gens)` | generators of the centre |
//! | `perm_group_normalizer(gens, hgens)` | generators of the normaliser of `<hgens>` |
//! | `perm_group_set_stabilizer(gens, points)` | generators of the stabiliser of a set of points |
//! | `perm_group_intersection(gens, hgens)` | generators of `<gens> ∩ <hgens>` |
//! | `perm_group_normal_closure(gens, xs)` | generators of the normal closure of `<xs>` |
//! | `perm_group_is_normal(gens, hgens)` | whether `<hgens>` is a normal subgroup |
//! | `perm_group_conjugating_element(gens, x, y)` | some `g` with `g x g⁻¹ = y`, or `false` |
//! | `perm_group_are_conjugate(gens, x, y)` | whether `x` and `y` are conjugate |
//! | `perm_group_conjugacy_classes(gens)` | `list(list(representative, size), ...)`, by element order |
//! | `perm_group_class_count(gens)` | the number of conjugacy classes |
//! | `perm_group_exponent(gens)` | the lcm of the element orders |
//! | `perm_group_sylow(gens, p)` | generators of a Sylow `p`-subgroup |
//! | `perm_group_composition_factors(gens)` | primes `p` of the factors `C_p` of a solvable group |
//! | `perm_group_is_primitive(gens)` | transitive with no non-trivial block system |
//! | `perm_group_blocks(gens, a, b)` | the finest block system with `a` and `b` in one block |
//! | `perm_group_minimal_block_systems(gens)` | every minimal non-trivial block system |
//! | `perm_group_random_element(gens, seed)` | a uniformly random element |
//!
//! The randomised algorithms (classes, Sylow subgroups) are seeded and
//! deterministic; their answers are verified (the class equation sums to the
//! order, the Sylow subgroup has the full prime-power order), and an
//! operator that cannot verify within its budget leaves the request
//! unevaluated rather than answer wrongly.

use num_bigint::BigInt;
use rand_v09::SeedableRng;

use super::V;
use super::def;
use super::idx;
use super::idxs;
use super::items;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::op::Arity;
use crate::graph::rule::Cx;
use crate::graph::rule::Installer;
use crate::graph::rule::RuleError;
use crate::kernels::permgroup;
use crate::kernels::permgroup::Group;
use crate::kernels::permgroup::Perm;

/// Largest degree accepted.
const MAX_DEGREE: usize = 1024;
/// Random samples a randomised algorithm may draw before giving up.
const TRIES: usize = 20_000;
/// Seed of the randomised algorithms (the answers do not depend on it).
const SEED: u64 = 0x5EED_6209;

fn read_perm(
    g: &Graph,
    n: NodeId,
) -> Option<Perm> {
    let images = idxs(g, n)?;
    let m = images.len();
    let perm: Option<Perm> = images
        .into_iter()
        .map(|v| v.checked_sub(1).filter(|&v| v < m).and_then(|v| u32::try_from(v).ok()))
        .collect();
    let perm = perm?;
    permgroup::is_permutation(&perm).then_some(perm)
}

fn read_gens(
    g: &Graph,
    n: NodeId,
) -> Option<(usize, Vec<Perm>)> {
    let perms: Option<Vec<Perm>> = items(g, n)?.into_iter().map(|p| read_perm(g, p)).collect();
    let perms = perms?;
    let degree = perms.first()?.len();
    (degree > 0 && degree <= MAX_DEGREE && perms.iter().all(|p| p.len() == degree))
        .then_some((degree, perms))
}

fn read_group(
    g: &Graph,
    n: NodeId,
) -> Option<Group> {
    let (degree, gens) = read_gens(g, n)?;
    Some(Group::new(degree, &gens))
}

fn perm_v(p: &[u32]) -> V {
    V::ints(p.iter().map(|&x| i64::from(x) + 1))
}

fn group_v(group: &Group) -> V {
    let gens = group.generators();
    if gens.is_empty() {
        return V::List(vec![perm_v(&permgroup::identity(group.degree()))]);
    }
    V::List(gens.iter().map(|p| perm_v(p)).collect())
}

fn big(n: num_bigint::BigUint) -> V {
    V::Int(BigInt::from(n))
}

/// A point argument (1-based) as a 0-based point of degree `n`.
fn point(
    g: &Graph,
    node: NodeId,
    n: usize,
) -> Option<u32> {
    let p = idx(g, node)?.checked_sub(1).filter(|&p| p < n)?;
    u32::try_from(p).ok()
}

/// A permutation argument of the group's degree.
fn element(
    g: &Graph,
    node: NodeId,
    group: &Group,
) -> Option<Perm> {
    read_perm(g, node).filter(|p| p.len() == group.degree())
}

fn perm_group_centralizer(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let x = element(cx.graph, *a.get(1)?, &group)?;
    Some(group_v(&group.centralizer(&x)))
}

fn perm_group_center(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    Some(group_v(&group.center()))
}

fn same_degree(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<(Group, Group)> {
    let group = read_group(cx.graph, *a.first()?)?;
    let other = read_group(cx.graph, *a.get(1)?)?;
    (other.degree() == group.degree()).then_some((group, other))
}

fn perm_group_normalizer(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (group, h) = same_degree(cx, a)?;
    Some(group_v(&group.normalizer(&h)))
}

fn perm_group_intersection(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (group, h) = same_degree(cx, a)?;
    Some(group_v(&group.intersection(&h)))
}

fn perm_group_is_normal(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (group, h) = same_degree(cx, a)?;
    Some(V::Bool(group.is_normal(&h)))
}

fn perm_group_normal_closure(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (group, h) = same_degree(cx, a)?;
    Some(group_v(&group.normal_closure(h.generators())))
}

fn perm_group_set_stabilizer(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let points: Option<Vec<u32>> = items(cx.graph, *a.get(1)?)?
        .into_iter()
        .map(|p| point(cx.graph, p, group.degree()))
        .collect();
    Some(group_v(&group.set_stabilizer(&points?)))
}

fn perm_group_conjugating_element(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let x = element(cx.graph, *a.get(1)?, &group)?;
    let y = element(cx.graph, *a.get(2)?, &group)?;
    Some(match group.conjugating_element(&x, &y) {
        | Some(g) => perm_v(&g),
        | None => V::Bool(false),
    })
}

fn perm_group_are_conjugate(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let x = element(cx.graph, *a.get(1)?, &group)?;
    let y = element(cx.graph, *a.get(2)?, &group)?;
    Some(V::Bool(group.are_conjugate(&x, &y)))
}

fn perm_group_conjugacy_classes(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let classes = group.conjugacy_classes(SEED, TRIES)?;
    Some(V::List(
        classes
            .into_iter()
            .map(|(rep, size)| V::List(vec![perm_v(&rep), big(size)]))
            .collect(),
    ))
}

fn perm_group_class_count(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    Some(V::uint(group.conjugacy_classes(SEED, TRIES)?.len()))
}

fn perm_group_exponent(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    Some(big(group.exponent(SEED)?))
}

fn perm_group_sylow(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let p = u64::try_from(idx(cx.graph, *a.get(1)?)?).ok()?;
    if p < 2 || !(2..p).take_while(|d| d * d <= p).all(|d| p % d != 0) {
        return None;
    }
    Some(group_v(&group.sylow_subgroup(p, SEED, TRIES)?))
}

fn perm_group_composition_factors(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    Some(V::ints(group.solvable_composition_factors()?))
}

fn perm_group_is_primitive(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    Some(V::Bool(permgroup::is_primitive(n, &gens)))
}

fn blocks_v(blocks: &[Vec<u32>]) -> V {
    V::List(
        blocks
            .iter()
            .map(|b| V::ints(b.iter().map(|&x| i64::from(x) + 1)))
            .collect(),
    )
}

fn perm_group_blocks(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    if permgroup::orbits(n, &gens).len() != 1 {
        return None;
    }
    let x = point(cx.graph, *a.get(1)?, n)?;
    let y = point(cx.graph, *a.get(2)?, n)?;
    Some(blocks_v(&permgroup::minimal_blocks(n, &gens, x, y)))
}

fn perm_group_minimal_block_systems(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (n, gens) = read_gens(cx.graph, *a.first()?)?;
    if permgroup::orbits(n, &gens).len() != 1 {
        return None;
    }
    Some(V::List(
        permgroup::minimal_block_systems(n, &gens)
            .iter()
            .map(|s| blocks_v(s))
            .collect(),
    ))
}

fn perm_group_random_element(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let group = read_group(cx.graph, *a.first()?)?;
    let seed = u64::try_from(idx(cx.graph, *a.get(1)?)?).ok()?;
    let mut rng = rand_v09::rngs::StdRng::seed_from_u64(seed);
    Some(perm_v(&group.random_element(&mut rng)))
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let table: &[(&str, u8, super::Run)] = &[
        ("perm_group_centralizer", 2, perm_group_centralizer),
        ("perm_group_center", 1, perm_group_center),
        ("perm_group_normalizer", 2, perm_group_normalizer),
        ("perm_group_intersection", 2, perm_group_intersection),
        ("perm_group_is_normal", 2, perm_group_is_normal),
        ("perm_group_normal_closure", 2, perm_group_normal_closure),
        ("perm_group_set_stabilizer", 2, perm_group_set_stabilizer),
        ("perm_group_conjugating_element", 3, perm_group_conjugating_element),
        ("perm_group_are_conjugate", 3, perm_group_are_conjugate),
        ("perm_group_conjugacy_classes", 1, perm_group_conjugacy_classes),
        ("perm_group_class_count", 1, perm_group_class_count),
        ("perm_group_exponent", 1, perm_group_exponent),
        ("perm_group_sylow", 2, perm_group_sylow),
        ("perm_group_composition_factors", 1, perm_group_composition_factors),
        ("perm_group_is_primitive", 1, perm_group_is_primitive),
        ("perm_group_blocks", 3, perm_group_blocks),
        ("perm_group_minimal_block_systems", 1, perm_group_minimal_block_systems),
        ("perm_group_random_element", 2, perm_group_random_element),
    ];
    for &(name, arity, run) in table {
        def(i, name, Arity::Fixed(arity), run)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;

    const S5: &str = "list(list(2,3,4,5,1), list(2,1,3,4,5))";

    #[test]
    fn subgroup_constructions() {
        assert_eq!(
            s(&format!("perm_group_order(perm_group_centralizer({S5}, list(2,3,4,5,1)))")),
            "5"
        );
        assert_eq!(
            s(&format!("perm_group_order(perm_group_set_stabilizer({S5}, list(1, 2)))")),
            "12"
        );
        assert_eq!(s(&format!("perm_group_order(perm_group_center({S5}))")), "1");
        assert_eq!(s(&format!("perm_group_order(perm_group_sylow({S5}, 2))")), "8");
        assert_eq!(s(&format!("perm_group_class_count({S5})")), "7");
        assert_eq!(s(&format!("perm_group_exponent({S5})")), "60");
        assert_eq!(
            s(&format!("perm_group_are_conjugate({S5}, list(2,1,3,4,5), list(1,2,3,5,4))")),
            "true"
        );
        assert_eq!(s(&format!("perm_group_is_primitive({S5})")), "true");
        assert_eq!(
            s("perm_group_blocks(list(list(2,3,4,1), list(1,4,3,2)), 1, 3)"),
            "list(list(1, 3), list(2, 4))"
        );
        assert_eq!(
            s("perm_group_composition_factors(list(list(2,3,4,1), list(2,1,3,4)))"),
            "list(2, 3, 2, 2)"
        );
    }
}
