//! The quantifier-elimination loop on top of the cylindrical algebraic
//! decomposition.
//!
//! Quantifiers are eliminated from the inside out. For one block
//! `Q x_1 ... Q x_k . f` with a quantifier-free body `f`, the atom
//! polynomials are decomposed over the variables `(parameters, block)` with
//! the parameters first, the truth of `f` is evaluated on every cell from
//! its sign vector, and the truth values of all cells over one parameter
//! cell are combined with `or` (existential) or `and` (universal). The
//! resulting quantifier-free formula is a disjunction of sign conditions of
//! the parameter projection polynomials, minimised by dropping literals
//! that do not exclude any false cell.
//!
//! Every block result is verified independently before it is returned: on a
//! deterministic grid of rational parameter assignments the original block
//! is decided by a separate decomposition over the block variables alone,
//! and the returned formula is evaluated exactly; any disagreement makes
//! the whole request unevaluated (`None`). This is a safety net on top of
//! the exact construction, never a source of answers.

#![allow(
    clippy::missing_const_for_fn,
    clippy::must_use_candidate,
    clippy::option_if_let_else,
    clippy::redundant_pub_crate,
    clippy::single_call_fn
)]

use std::collections::BTreeMap;
use std::collections::BTreeSet;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::op::core;
use crate::rules::poly::repr::Gens;
use crate::rules::poly::repr::Limits as ReprLimits;
use crate::rules::poly::repr::from_term;

use super::cad;
use super::cad::CadLimits;
use super::cad::ParamCell;
use super::field::Elem;
use super::field::Field;
use super::formula::Atom;
use super::formula::Cmp;
use super::formula::Formula;
use super::formula::Parsed;
use super::formula::repr_to_mpoly;
use super::poly::MPoly;

/// Most rational assignments used by the independent verification.
const MAX_SAMPLES: usize = 28;
/// Most cells described by the output formula before it is abandoned.
const MAX_CUBES: usize = 512;

/// The distinct non-zero atom polynomials of a formula.
#[must_use]
fn atom_polys(formula: &Formula) -> Vec<MPoly> {
    let mut set = BTreeSet::new();
    collect_atoms(formula, &mut set);
    set.into_iter().collect()
}

fn collect_atoms(
    formula: &Formula,
    out: &mut BTreeSet<MPoly>,
) {
    match formula {
        | Formula::Truth(_) => {},
        | Formula::Atom(atom) => {
            if !atom.poly.is_zero() {
                out.insert(atom.poly.clone());
            }
        },
        | Formula::Not(x) => collect_atoms(x, out),
        | Formula::And(xs) | Formula::Or(xs) => {
            for x in xs {
                collect_atoms(x, out);
            }
        },
        | Formula::Exists(_, body) | Formula::Forall(_, body) => collect_atoms(body, out),
    }
}

/// Renames the variables of the atoms of a quantifier-free formula.
#[must_use]
fn remap(
    formula: &Formula,
    map: &[u32],
) -> Formula {
    match formula {
        | Formula::Truth(value) => Formula::Truth(*value),
        | Formula::Atom(atom) => Formula::Atom(Atom {
            cmp: atom.cmp,
            poly: atom.poly.renamed(map),
        }),
        | Formula::Not(x) => Formula::Not(Box::new(remap(x, map))),
        | Formula::And(xs) => Formula::And(xs.iter().map(|x| remap(x, map)).collect()),
        | Formula::Or(xs) => Formula::Or(xs.iter().map(|x| remap(x, map)).collect()),
        | Formula::Exists(_, _) | Formula::Forall(_, _) => formula.clone(),
    }
}

/// Substitutes rational values for variables, keeping the others.
#[must_use]
fn specialize(
    formula: &Formula,
    values: &BTreeMap<u32, BigRational>,
) -> Formula {
    match formula {
        | Formula::Truth(value) => Formula::Truth(*value),
        | Formula::Atom(atom) => Formula::Atom(Atom {
            cmp: atom.cmp,
            poly: atom.poly.substitute_rational(values),
        }),
        | Formula::Not(x) => Formula::Not(Box::new(specialize(x, values))),
        | Formula::And(xs) => Formula::And(xs.iter().map(|x| specialize(x, values)).collect()),
        | Formula::Or(xs) => Formula::Or(xs.iter().map(|x| specialize(x, values)).collect()),
        | Formula::Exists(_, _) | Formula::Forall(_, _) => formula.clone(),
    }
}

/// The truth function of a quantifier-free formula over a sign vector
/// aligned with `polys`.
fn truth_function<'a>(
    formula: &'a Formula,
    polys: &'a [MPoly],
) -> impl Fn(&[i32]) -> bool + 'a {
    move |signs: &[i32]| {
        formula.truth_with(&|p: &MPoly| match polys.binary_search(p) {
            | Ok(index) => signs.get(index).copied().unwrap_or(0),
            | Err(_) => 0,
        })
    }
}

/// Decides a quantifier-free formula with the block variables quantified,
/// using a decomposition over those variables only.
#[must_use]
fn decide_block(
    body: &Formula,
    variables: usize,
    exists: bool,
    limits: &CadLimits,
) -> Option<bool> {
    let polys = atom_polys(body);
    let truth = truth_function(body, &polys);
    let (cells, _) = cad::quantify(&polys, variables, 0, exists, &truth, limits)?;
    cells.first().map(|cell| cell.truth)
}

/// One elimination step for one block.
#[must_use]
fn eliminate_block(
    body: &Formula,
    block: &[u32],
    exists: bool,
    total_variables: usize,
    limits: &CadLimits,
) -> Option<Formula> {
    if block.is_empty() {
        return Some(body.clone());
    }
    let mut support = BTreeSet::new();
    body.atom_variables(&mut support);
    let mut present: Vec<u32> = Vec::new();
    for &variable in block {
        if support.contains(&variable) && !present.contains(&variable) {
            present.push(variable);
        }
    }
    if present.is_empty() {
        return Some(body.clone());
    }
    let mut block_set = BTreeSet::new();
    block_set.extend(present.iter().copied());
    let params: Vec<u32> = support
        .iter()
        .copied()
        .filter(|v| !block_set.contains(v))
        .collect();
    let mut order = params.clone();
    order.extend(present.iter().copied());
    let mut map = vec![0_u32; total_variables];
    for (index, slot) in map.iter_mut().enumerate() {
        *slot = u32::try_from(index).ok()?;
    }
    for (local, &global) in order.iter().enumerate() {
        if let Some(slot) = map.get_mut(usize::try_from(global).ok()?) {
            *slot = u32::try_from(local).ok()?;
        }
    }
    let local_body = remap(body, &map);
    let polys = atom_polys(&local_body);
    let truth = truth_function(&local_body, &polys);
    let (cells, parameter_polys) =
        cad::quantify(&polys, order.len(), params.len(), exists, &truth, limits)?;
    if params.is_empty() {
        return Some(Formula::Truth(cells.first()?.truth));
    }
    let mut back = vec![0_u32; params.len()];
    for (index, &global) in params.iter().enumerate() {
        if let Some(slot) = back.get_mut(index) {
            *slot = global;
        }
    }
    let global_polys: Vec<MPoly> = parameter_polys.iter().map(|p| p.renamed(&back)).collect();
    let cubes = minimise(
        &cells,
        &global_polys,
        cells.first().map_or(0, |c| c.signs.len()),
    )?;
    let mut parts = Vec::with_capacity(cubes.len());
    for cube in &cubes {
        let mut literals = Vec::new();
        for (index, sign) in cube.iter().enumerate() {
            if let (Some(sign), Some(poly)) = (sign, global_polys.get(index)) {
                literals.push(Formula::Atom(Atom {
                    cmp: Cmp::from_sign(*sign),
                    poly: poly.clone(),
                }));
            }
        }
        parts.push(Formula::and(literals));
    }
    let result = Formula::or(parts);
    verify(
        &local_body,
        &present,
        &params,
        exists,
        &result,
        total_variables,
        limits,
    )?;
    Some(result)
}

/// A cube: for each parameter projection polynomial either a fixed sign or
/// a don't-care.
type Cube = Vec<Option<i32>>;

/// Whether a cell matches a cube.
fn matches(
    cube: &Cube,
    signs: &[i32],
) -> bool {
    cube.iter()
        .enumerate()
        .all(|(index, value)| match (value, signs.get(index)) {
            | (Some(wanted), Some(actual)) => wanted == actual,
            | (None, _) => true,
            | (Some(_), None) => false,
        })
}

/// Greedy minimisation of the true cells' sign conditions.
#[must_use]
fn minimise(
    cells: &[ParamCell],
    _polys: &[MPoly],
    width: usize,
) -> Option<Vec<Cube>> {
    let has_false = cells.iter().any(|c| !c.truth);
    if !has_false {
        return Some(vec![vec![None; width]]);
    }
    let mut cubes: Vec<Cube> = Vec::new();
    for cell in cells.iter().filter(|c| c.truth) {
        let mut cube: Cube = cell.signs.iter().take(width).map(|&s| Some(s)).collect();
        while cube.len() < width {
            cube.push(None);
        }
        for index in 0..width {
            let saved = cube.get(index).copied().unwrap_or(None);
            if saved.is_none() {
                continue;
            }
            if let Some(slot) = cube.get_mut(index) {
                *slot = None;
            }
            let conflict = cells
                .iter()
                .filter(|c| !c.truth)
                .any(|c| matches(&cube, &c.signs));
            if conflict && let Some(slot) = cube.get_mut(index) {
                *slot = saved;
            }
        }
        if !cubes.contains(&cube) {
            cubes.push(cube);
        }
        if cubes.len() > MAX_CUBES {
            return None;
        }
    }
    // A cube that is a specialisation of another one is redundant.
    let mut kept: Vec<Cube> = Vec::with_capacity(cubes.len());
    for (index, cube) in cubes.iter().enumerate() {
        let subsumed = cubes.iter().enumerate().any(|(other, candidate)| {
            other != index && covers(candidate, cube) && (!covers(cube, candidate) || other < index)
        });
        if !subsumed {
            kept.push(cube.clone());
        }
    }
    Some(kept)
}

/// Whether cube `general` contains cube `specific`.
fn covers(
    general: &Cube,
    specific: &Cube,
) -> bool {
    general
        .iter()
        .enumerate()
        .all(|(index, value)| match value {
            | None => true,
            | Some(wanted) => specific.get(index).copied().flatten() == Some(*wanted),
        })
}

/// Verifies a block result on rational parameter assignments.
fn verify(
    body: &Formula,
    block: &[u32],
    params: &[u32],
    exists: bool,
    result: &Formula,
    total_variables: usize,
    limits: &CadLimits,
) -> Option<()> {
    let block_count = block.len();
    let parameter_count = params.len();
    let samples = assignments(parameter_count);
    for assignment in &samples {
        let mut values = vec![BigRational::zero(); total_variables];
        for (index, value) in assignment.iter().enumerate() {
            let global = params.get(index)?;
            if let Some(slot) = values.get_mut(usize::try_from(*global).ok()?) {
                *slot = value.clone();
            }
        }
        let expected = result.eval_rational(&values)?;
        let mut local_values: BTreeMap<u32, BigRational> = BTreeMap::new();
        for (index, value) in assignment.iter().enumerate() {
            local_values.insert(u32::try_from(index).ok()?, value.clone());
        }
        let specialized = specialize(body, &local_values);
        let mut map = vec![0_u32; parameter_count.saturating_add(block_count)];
        for index in 0..block_count {
            if let Some(slot) = map.get_mut(parameter_count.saturating_add(index)) {
                *slot = u32::try_from(index).ok()?;
            }
        }
        let local = remap(&specialized, &map);
        let actual = decide_block(&local, block_count, exists, limits)?;
        if actual != expected {
            return None;
        }
    }
    Some(())
}

/// The rational parameter assignments used by verification.
#[must_use]
fn assignments(parameters: usize) -> Vec<Vec<BigRational>> {
    let palette: Vec<BigRational> = [
        (0, 1),
        (1, 1),
        (-1, 1),
        (1, 2),
        (-1, 2),
        (2, 1),
        (-2, 1),
        (1, 3),
        (-1, 3),
        (3, 1),
        (-3, 1),
        (1, 4),
        (-1, 4),
        (4, 1),
        (-4, 1),
        (5, 1),
        (-5, 1),
        (2, 3),
        (-2, 3),
        (3, 2),
        (-3, 2),
    ]
    .iter()
    .map(|&(n, d)| BigRational::new(BigInt::from(n), BigInt::from(d)))
    .collect();
    if parameters == 0 {
        return vec![Vec::new()];
    }
    if parameters == 1 {
        return palette.into_iter().map(|value| vec![value]).collect();
    }
    let count = if parameters <= 2 {
        MAX_SAMPLES
    } else {
        MAX_SAMPLES / 2
    };
    let mut out: Vec<Vec<BigRational>> = Vec::new();
    let push = |values: Vec<BigRational>, out: &mut Vec<Vec<BigRational>>| {
        if !out.contains(&values) {
            out.push(values);
        }
    };
    push(vec![BigRational::zero(); parameters], &mut out);
    push(vec![BigRational::one(); parameters], &mut out);
    for index in 0..parameters {
        for sign in [1_i64, -1] {
            let mut values = vec![BigRational::zero(); parameters];
            if let Some(slot) = values.get_mut(index) {
                *slot = BigRational::from_integer(BigInt::from(sign));
            }
            push(values, &mut out);
        }
    }
    for step in 0..count {
        let mut values = Vec::with_capacity(parameters);
        for index in 0..parameters {
            let position = step
                .saturating_mul(3)
                .saturating_add(index.saturating_mul(5))
                .saturating_add(step.saturating_mul(index) % 7)
                % palette.len();
            values.push(
                palette
                    .get(position)
                    .cloned()
                    .unwrap_or_else(BigRational::zero),
            );
        }
        push(values, &mut out);
        if out.len() > MAX_SAMPLES.saturating_mul(2) {
            break;
        }
    }
    out
}

/// Eliminates every quantifier of a formula, from the inside out.
#[must_use]
fn eliminate(
    formula: &Formula,
    total_variables: usize,
    limits: &CadLimits,
) -> Option<Formula> {
    match formula {
        | Formula::Truth(_) | Formula::Atom(_) => Some(formula.clone()),
        | Formula::Not(x) => Some(eliminate(x, total_variables, limits)?.not()),
        | Formula::And(xs) => {
            let mut parts = Vec::with_capacity(xs.len());
            for x in xs {
                parts.push(eliminate(x, total_variables, limits)?);
            }
            Some(Formula::and(parts))
        },
        | Formula::Or(xs) => {
            let mut parts = Vec::with_capacity(xs.len());
            for x in xs {
                parts.push(eliminate(x, total_variables, limits)?);
            }
            Some(Formula::or(parts))
        },
        | Formula::Exists(vars, body) => {
            let body = eliminate(body, total_variables, limits)?;
            eliminate_block(&body, vars, true, total_variables, limits)
        },
        | Formula::Forall(vars, body) => {
            let body = eliminate(body, total_variables, limits)?;
            eliminate_block(&body, vars, false, total_variables, limits)
        },
    }
}

/// Runs the whole request: parses `formula`, adds an existential block for
/// the listed free variables, eliminates and builds the answer term.
#[must_use]
pub(crate) fn run(
    graph: &mut Graph,
    formula: NodeId,
    listed: &[NodeId],
) -> Option<NodeId> {
    let parsed = super::formula::parse(graph, formula, listed)?;
    let result = eliminate_all(&parsed, graph, listed)?;
    super::formula::term(graph, &parsed.ops, &result, &parsed.nodes)
}

/// Eliminates all quantifiers; `listed` free variables become existential.
#[must_use]
fn eliminate_all(
    parsed: &Parsed,
    graph: &Graph,
    listed: &[NodeId],
) -> Option<Formula> {
    let total = parsed.variables.len();
    let mut formula = parsed.formula.clone();
    let free = formula.free_variables();
    let mut extra: Vec<u32> = Vec::new();
    for &node in listed {
        let symbol = graph.as_symbol(node)?;
        let index = parsed
            .variables
            .iter()
            .position(|&s| s == symbol)
            .and_then(|i| u32::try_from(i).ok())?;
        if free.contains(&index) && !extra.contains(&index) {
            extra.push(index);
        }
    }
    if !extra.is_empty() {
        formula = Formula::Exists(extra, Box::new(formula));
    }
    eliminate(&formula, total, &CadLimits::default())
}

/// Handles the logic rule set's `exists`/`forall` when the body is a real
/// polynomial formula. `None` leaves the request to the Boolean path.
#[must_use]
pub(crate) fn try_quantifier(
    graph: &mut Graph,
    prefix: NodeId,
    body: NodeId,
    universal: bool,
) -> Option<NodeId> {
    let extras: Vec<NodeId> = if graph.op(prefix) == core::LIST {
        graph.children(prefix).to_vec()
    } else {
        vec![prefix]
    };
    if extras.iter().any(|&node| graph.as_symbol(node).is_none()) {
        return None;
    }
    let parsed = super::formula::parse(graph, body, &extras)?;
    let mut block = Vec::with_capacity(extras.len());
    for &node in &extras {
        let symbol = graph.as_symbol(node)?;
        let index = parsed
            .variables
            .iter()
            .position(|&s| s == symbol)
            .and_then(|index| u32::try_from(index).ok())
            .filter(|index| !block.contains(index));
        if let Some(index) = index {
            block.push(index);
        }
    }
    let mut support = BTreeSet::new();
    parsed.formula.atom_variables(&mut support);
    block.retain(|index| support.contains(index));
    let formula = if block.is_empty() {
        parsed.formula.clone()
    } else if universal {
        Formula::Forall(block, Box::new(parsed.formula.clone()))
    } else {
        Formula::Exists(block, Box::new(parsed.formula.clone()))
    };
    let total = parsed.variables.len();
    let result = eliminate(&formula, total, &CadLimits::default())?;
    super::formula::term(graph, &parsed.ops, &result, &parsed.nodes)
}

/// Reads a list of polynomials in the given variable order.
#[must_use]
pub(crate) fn read_polynomials(
    graph: &mut Graph,
    polys: NodeId,
    variables: &[NodeId],
) -> Option<Vec<MPoly>> {
    if graph.op(polys) != core::LIST {
        return None;
    }
    let terms = graph.children(polys).to_vec();
    let mut gens = Gens::default();
    for &variable in variables {
        gens.index(graph, variable);
    }
    let mut out = Vec::with_capacity(terms.len());
    for term in terms {
        let repr = from_term(graph, &mut gens, term, ReprLimits::default())?;
        if gens.len() != variables.len() {
            return None;
        }
        out.push(repr_to_mpoly(&repr, variables.len())?);
    }
    Some(out)
}

/// The coordinates of a cell sample as exact terms.
#[must_use]
pub(crate) fn sample_terms(
    graph: &mut Graph,
    field: &Field,
    coords: &[Elem],
) -> Option<Vec<NodeId>> {
    let mut out = Vec::with_capacity(coords.len());
    for coordinate in coords {
        out.push(super::formula::elem_term(graph, field, coordinate)?);
    }
    Some(out)
}

/// Builds the `list(sample_point, open_sections, signs)` description of a
/// cell.
#[must_use]
pub(crate) fn cell_term(
    graph: &mut Graph,
    cell: &cad::Cell,
) -> Option<NodeId> {
    let sample = sample_terms(graph, &cell.field, &cell.sample)?;
    let sample = graph.node(core::LIST, &sample);
    let open: Vec<NodeId> = cell
        .open
        .iter()
        .map(|&flag| graph.lit(crate::graph::Payload::Bool(flag)))
        .collect();
    let open = graph.node(core::LIST, &open);
    let signs: Vec<NodeId> = cell
        .signs
        .iter()
        .map(|&sign| graph.int(i64::from(sign)))
        .collect();
    let signs = graph.node(core::LIST, &signs);
    Some(graph.node(core::LIST, &[sample, open, signs]))
}
