//! Crystallographic and molecular point groups.
//!
//! A point group is named by its Schoenflies symbol, a symbol such as `C2v`,
//! `D4h`, `Td`, `Oh`, `Ih`: the families `Cn`, `Cnv`, `Cnh`, `Dn`, `Dnh`,
//! `Dnd` and `Sn` (`n` even) for `n <= 12`, `Ci`, `Cs` (`S6` is `C3i`),
//! and the polyhedral groups `T`, `Td`, `Th`, `O`, `Oh`, `I` (the
//! constant `I` is accepted as the icosahedral rotation group) and `Ih`.
//! Groups are built from exact orthogonal 3x3 matrices in the standard
//! orientation (principal axis `z`, a `C2` axis along `x` for the dihedral
//! groups, cubic axes along the coordinates).
//!
//! The character table is computed by the Dixon algorithm of
//! [`representations`](super::representations) on the multiplication table
//! of the matrices; complex-conjugate pairs of one-dimensional
//! representations are combined into the physically irreducible
//! two-dimensional `E`. The irreducible representations carry **Mulliken
//! labels** found from the geometry of the classes alone: the letter by the
//! dimension (`A`, `B` for 1, `E` 2, `T` 3, `G` 4, `H` 5); `A`/`B` by the
//! character of the principal rotation (or of `S_2n` for `D_nd` and `S_2n`);
//! the subscripts `1`/`2` by the character of the perpendicular `C2` (else
//! of the vertical mirror), `B1, B2, B3` of `D2` by which `C2` axis is
//! preserved, `E_k`/`T_k` by the character of the generator; `g`/`u` by the
//! inversion, `'`/`''` by the horizontal mirror when there is no inversion.
//!
//! Classes are labelled in the usual way: `E`, `8C3`, `3C2`, `6S4`, `6σd`,
//! `i`, `3σh`, `C2'`, `C5^2`, ... and ordered with the identity first, then
//! rotations by decreasing order, inversion, improper rotations, mirrors.
//!
//! A molecule is a list of atoms, each `list(x, y, z)` or `list(label, x, y,
//! z)` (atoms with different labels are different elements). Symmetry
//! operations are found numerically (default tolerance `1e-3`; an optional
//! last argument gives another) about the centroid by testing rotations,
//! improper rotations and mirrors about candidate axes (atoms, pair sums and
//! differences, cross products); the point group is recognised from the
//! class structure.
//!
//! | operator | value |
//! |---|---|
//! | `point_group(name)` | the `group` term whose elements are 3x3 matrices |
//! | `point_group_order(name)`, `point_group_is_crystallographic(name)` | order; all element orders in `{1, 2, 3, 4, 6}` |
//! | `point_group_classes(name)`, `point_group_class_sizes(name)` | class labels (symbols) and sizes |
//! | `point_group_irreps(name)` | the Mulliken labels |
//! | `point_group_irrep_dimensions(name)` | the dimensions of the irreps (2 for a combined complex pair `E`) |
//! | `point_group_character_table(name)` | the rows (irreps) over the classes; exact entries |
//! | `point_group_decompose(name, chi)` | a class function (one value per class, or per element) as a sum of irreps, e.g. `2*A1 + B2` |
//! | `point_group_multiplicities(name, chi)` | the multiplicities aligned with `point_group_irreps` |
//! | `point_group_vector_character(name)`, `point_group_rotation_character(name)` | the characters of `(x, y, z)` and of the rotations `(Rx, Ry, Rz)` |
//! | `point_group_ir_active(name)`, `point_group_raman_active(name)` | irreps of `x, y, z` / of `x^2, xy, ...` |
//! | `point_group_function_irreps(name)` | for each of `x, y, z, x^2, y^2, z^2, x*y, x*z, y*z` the irreps in which it has a component |
//! | `point_group_hm(name)`, `point_group_from_hm(hm)` | Schoenflies <-> Hermann–Mauguin (crystallographic groups) |
//! | `point_group_crystal_system(name)` | the crystal system |
//! | `crystallographic_point_groups()` | the 32 crystal classes |
//! | `crystal_systems()` | the 7 systems: `list(system, holohedry, point groups, lattices)` |
//! | `bravais_lattices()` | the 14 lattices: `list(Pearson symbol, system, centring, holohedry)` |
//! | `crystallographic_restriction(n)` | whether an `n`-fold rotation is compatible with a lattice (`n` in `{1,2,3,4,6}`); optional second argument: the dimension |
//! | `crystallographic_min_dimension(n)` | the least dimension of a lattice with an `n`-fold symmetry |
//! | `molecule_symmetry_operations(atoms)` | the matrices of the symmetry group |
//! | `molecule_point_group(atoms)` | the Schoenflies symbol (`Cinfv`, `Dinfh` for linear molecules) |
//! | `molecule_decomposition(atoms)` | `list(irrep, n_3N, n_trans, n_rot, n_vib)` rows |
//! | `molecule_vibrations(atoms)` | `Gamma_vib = Gamma_3N - Gamma_trans - Gamma_rot` as a sum of irreps |
//! | `molecule_vibrations_table(atoms)` | `list(list(irrep, multiplicity), ...)` of the vibrations |
//! | `molecule_spectroscopy(atoms)` | `list(irrep, multiplicity, ir_active, raman_active)` of the vibrations |
//!
//! Molecule operators accept an optional last argument, the tolerance.
//! Linear molecules have infinite groups: only their name is reported.

// Index loops over several parallel tables and float comparisons with explicit
// tolerances read better than the iterator forms; `as` casts convert small,
// non-negative, bounded values.
#![allow(
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::cast_sign_loss,
    clippy::manual_midpoint
)]

use std::collections::HashMap;

use num_complex::Complex64;

use super::apply;
use super::def;
use super::groups::build_group;
use super::group_theory::Tab;
use super::items;
use super::matrix;
use super::prod;
use super::representations::character_table;
use super::sum;
use super::V;
use crate::graph::rule::Installer;
use crate::graph::Arity;
use crate::graph::ClassId;
use crate::graph::Cx;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::RuleError;

pub(crate) type M3 = [[f64; 3]; 3];

const ID3: M3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const PI: f64 = std::f64::consts::PI;
const ORDER_LIMIT: usize = 240;

fn mm(
    a: &M3,
    b: &M3,
) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

fn mneg(a: &M3) -> M3 {
    let mut r = *a;
    for row in &mut r {
        for x in row {
            *x = -*x;
        }
    }
    r
}

fn trace(a: &M3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

fn det(a: &M3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn apply_m(
    a: &M3,
    v: [f64; 3],
) -> [f64; 3] {
    [0, 1, 2].map(|i| (0..3).map(|k| a[i][k] * v[k]).sum())
}

fn norm(v: [f64; 3]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = norm(v);
    v.map(|x| x / n)
}

fn dot(
    a: [f64; 3],
    b: [f64; 3],
) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn cross(
    a: [f64; 3],
    b: [f64; 3],
) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// The rotation by `theta` about the axis `u` (Rodrigues).
fn rot(
    u: [f64; 3],
    theta: f64,
) -> M3 {
    let u = unit(u);
    let (c, s) = (theta.cos(), theta.sin());
    let mut r = ID3;
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = c * ID3[i][j] + (1.0 - c) * u[i] * u[j];
        }
    }
    r[0][1] -= s * u[2];
    r[0][2] += s * u[1];
    r[1][0] += s * u[2];
    r[1][2] -= s * u[0];
    r[2][0] -= s * u[1];
    r[2][1] += s * u[0];
    r
}

/// The reflection in the plane with normal `n`.
fn refl(n: [f64; 3]) -> M3 {
    let n = unit(n);
    let mut r = ID3;
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] -= 2.0 * n[i] * n[j];
        }
    }
    r
}

fn key(a: &M3) -> [i64; 9] {
    let mut k = [0_i64; 9];
    for i in 0..3 {
        for j in 0..3 {
            let v = (a[i][j] * 1e5).round() as i64;
            k[3 * i + j] = v;
        }
    }
    k
}

/// The group generated by `gens`, or `None` if it exceeds the limit.
fn generate(gens: &[M3]) -> Option<Vec<M3>> {
    let mut index: HashMap<[i64; 9], usize> = HashMap::new();
    let mut all = vec![ID3];
    index.insert(key(&ID3), 0);
    let mut head = 0;
    while head < all.len() {
        let x = all[head];
        head += 1;
        for g in gens {
            let y = mm(g, &x);
            if let std::collections::hash_map::Entry::Vacant(e) = index.entry(key(&y)) {
                e.insert(all.len());
                all.push(y);
                if all.len() > ORDER_LIMIT {
                    return None;
                }
            }
        }
    }
    Some(all)
}

// ----------------------------------------------------------------------
// Named groups
// ----------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Family {
    C(usize),
    Cv(usize),
    Ch(usize),
    D(usize),
    Dh(usize),
    Dd(usize),
    S(usize),
    T,
    Td,
    Th,
    O,
    Oh,
    I,
    Ih,
}

fn parse_family(name: &str) -> Option<Family> {
    match name {
        | "Ci" => return Some(Family::S(2)),
        | "Cs" => return Some(Family::Ch(1)),
        | "C3i" => return Some(Family::S(6)),
        | "T" => return Some(Family::T),
        | "Td" => return Some(Family::Td),
        | "Th" => return Some(Family::Th),
        | "O" => return Some(Family::O),
        | "Oh" => return Some(Family::Oh),
        | "I" => return Some(Family::I),
        | "Ih" => return Some(Family::Ih),
        | _ => {},
    }
    let mut chars = name.chars();
    let letter = chars.next()?;
    let rest: String = chars.collect();
    let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
    let suffix = &rest[digits.len()..];
    let n: usize = digits.parse().ok()?;
    if !(1..=12).contains(&n) {
        return None;
    }
    match (letter, suffix) {
        | ('C', "") => Some(Family::C(n)),
        | ('C', "v") if n >= 2 => Some(Family::Cv(n)),
        | ('C', "h") => Some(Family::Ch(n)),
        | ('D', "") if n >= 2 => Some(Family::D(n)),
        | ('D', "h") if n >= 2 => Some(Family::Dh(n)),
        | ('D', "d") if n >= 2 => Some(Family::Dd(n)),
        | ('S', "") if n.is_multiple_of(2) => Some(Family::S(n)),
        | _ => None,
    }
}

fn family_generators(f: Family) -> Vec<M3> {
    let z = [0.0, 0.0, 1.0];
    let x = [1.0, 0.0, 0.0];
    let sigma_h = refl(z);
    let cn = |n: usize| rot(z, 2.0 * PI / n as f64);
    let c2x = rot(x, PI);
    let c3 = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let c2z = rot(z, PI);
    let inv = mneg(&ID3);
    match f {
        | Family::C(n) => vec![cn(n)],
        | Family::Cv(n) => vec![cn(n), refl([0.0, 1.0, 0.0])],
        | Family::Ch(n) => vec![cn(n), sigma_h],
        | Family::D(n) => vec![cn(n), c2x],
        | Family::Dh(n) => vec![cn(n), c2x, sigma_h],
        | Family::Dd(n) => vec![cn(n), c2x, mm(&sigma_h, &rot(z, PI / n as f64))],
        | Family::S(n) => vec![mm(&sigma_h, &cn(n))],
        | Family::T => vec![c2z, c3],
        | Family::Td => vec![c2z, c3, refl([1.0, -1.0, 0.0])],
        | Family::Th => vec![c2z, c3, inv],
        | Family::O => vec![c2z, c3, rot(z, PI / 2.0)],
        | Family::Oh => vec![c2z, c3, rot(z, PI / 2.0), inv],
        | Family::I | Family::Ih => {
            let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let mut g = vec![c2z, c3, rot([0.0, 1.0, phi], 2.0 * PI / 5.0)];
            if f == Family::Ih {
                g.push(inv);
            }
            g
        },
    }
}

fn family_matrices(f: Family) -> Option<Vec<M3>> {
    generate(&family_generators(f))
}

fn family_name(f: Family) -> String {
    match f {
        | Family::C(n) => format!("C{n}"),
        | Family::Cv(n) => format!("C{n}v"),
        | Family::Ch(1) => "Cs".into(),
        | Family::Ch(n) => format!("C{n}h"),
        | Family::D(n) => format!("D{n}"),
        | Family::Dh(n) => format!("D{n}h"),
        | Family::Dd(n) => format!("D{n}d"),
        | Family::S(2) => "Ci".into(),
        | Family::S(n) => format!("S{n}"),
        | Family::T => "T".into(),
        | Family::Td => "Td".into(),
        | Family::Th => "Th".into(),
        | Family::O => "O".into(),
        | Family::Oh => "Oh".into(),
        | Family::I => "I".into(),
        | Family::Ih => "Ih".into(),
    }
}

// ----------------------------------------------------------------------
// Operations
// ----------------------------------------------------------------------

#[derive(Clone, Debug)]
enum Kind {
    E,
    Inv,
    Rot { axis: [f64; 3], n: usize, k: usize },
    Mirror { normal: [f64; 3] },
    Imp { axis: [f64; 3], n: usize, k: usize },
}

/// The smallest `n` with `angle / 2 pi = k / n`.
fn fraction(angle: f64) -> (usize, usize) {
    let t = angle / (2.0 * PI);
    for n in 1..=120 {
        let k = (t * n as f64).round();
        if (t * n as f64 - k).abs() < 1e-6 {
            return (n, k as usize);
        }
    }
    (120, (t * 120.0).round() as usize)
}

fn canonical_axis(mut u: [f64; 3]) -> [f64; 3] {
    for c in u {
        if c.abs() > 1e-6 {
            if c < 0.0 {
                u = u.map(|x| -x);
            }
            break;
        }
    }
    u
}

fn classify(m: &M3) -> Kind {
    let proper = det(m) > 0.0;
    let r = if proper { *m } else { mneg(m) };
    let c = ((trace(&r) - 1.0) / 2.0).clamp(-1.0, 1.0);
    let theta = c.acos();
    if theta < 1e-6 {
        return if proper { Kind::E } else { Kind::Inv };
    }
    let axis = if PI - theta < 1e-6 {
        // r = 2 u u^T - I
        let diag = [(r[0][0] + 1.0) / 2.0, (r[1][1] + 1.0) / 2.0, (r[2][2] + 1.0) / 2.0];
        let i = (0..3).max_by(|&a, &b| diag[a].total_cmp(&diag[b])).unwrap_or(0);
        let col = [(r[0][i] + ID3[0][i]) / 2.0, (r[1][i] + ID3[1][i]) / 2.0, (r[2][i] + ID3[2][i]) / 2.0];
        canonical_axis(unit(col))
    } else {
        canonical_axis(unit([r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]]))
    };
    if proper {
        let (n, k) = fraction(theta);
        Kind::Rot { axis, n, k }
    } else if PI - theta < 1e-6 {
        Kind::Mirror { normal: axis }
    } else {
        let (n, k) = fraction(PI - theta);
        Kind::Imp { axis, n, k }
    }
}

// ----------------------------------------------------------------------
// The point group
// ----------------------------------------------------------------------

/// An irreducible representation with its Mulliken label.
#[derive(Clone, Debug)]
pub(crate) struct Irrep {
    pub(crate) label: String,
    pub(crate) dim: usize,
    /// Character values over the classes in display order.
    pub(crate) chi: Vec<f64>,
    /// `<chi, chi>`: 1, or 2 for a combined complex-conjugate pair.
    pub(crate) norm: f64,
}

/// A point group with classes and Mulliken-labelled irreducibles.
pub(crate) struct PointGroup {
    pub(crate) name: String,
    pub(crate) mats: Vec<M3>,
    pub(crate) tab: Tab,
    pub(crate) class_elems: Vec<Vec<usize>>,
    pub(crate) class_labels: Vec<String>,
    pub(crate) irreps: Vec<Irrep>,
    pub(crate) class_of: Vec<usize>,
}

/// A merged row with its preliminary label.
struct Pre {
    dim: usize,
    norm: f64,
    chi: Vec<f64>,
    letter: char,
    suffix: String,
    sub: Option<usize>,
}

/// Geometry of the group used for the labels.
struct Geo {
    principal: Option<[f64; 3]>,
    polyhedral: bool,
    inversion: bool,
    nmax: usize,
}

fn aligned(u: [f64; 3]) -> bool {
    u.iter().filter(|c| c.abs() > 1e-6).count() == 1
}

fn parallel(
    a: [f64; 3],
    b: [f64; 3],
) -> bool {
    dot(a, b).abs() > 1.0 - 1e-6
}

impl Geo {
    fn new(kinds: &[Kind]) -> Self {
        let inversion = kinds.iter().any(|k| matches!(k, Kind::Inv));
        let mut nmax = 1;
        for k in kinds {
            if let Kind::Rot { n, .. } = k {
                nmax = nmax.max(*n);
            }
        }
        let mut axes: Vec<[f64; 3]> = Vec::new();
        for k in kinds {
            if let Kind::Rot { axis, n, .. } = k {
                if *n == nmax && !axes.iter().any(|&a| parallel(a, *axis)) {
                    axes.push(*axis);
                }
            }
        }
        let polyhedral = nmax >= 3 && axes.len() > 1;
        let principal = if nmax == 1 || polyhedral {
            None
        } else if axes.len() == 1 {
            Some(axes[0])
        } else {
            // several C2 axes: D2d has an S4 axis, D2 and D2h have none
            kinds.iter().find_map(|k| match k {
                | Kind::Imp { axis, n: 4, .. } => Some(*axis),
                | _ => None,
            })
        };
        Self { principal, polyhedral, inversion, nmax }
    }

    fn perpendicular(
        &self,
        u: [f64; 3],
    ) -> bool {
        self.principal.is_some_and(|p| dot(p, u).abs() < 1e-6)
    }

    fn along(
        &self,
        u: [f64; 3],
    ) -> bool {
        self.principal.is_some_and(|p| parallel(p, u))
    }
}

fn axis_name(u: [f64; 3]) -> Option<&'static str> {
    if !aligned(u) {
        return None;
    }
    Some(if u[0].abs() > 0.5 {
        "x"
    } else if u[1].abs() > 0.5 {
        "y"
    } else {
        "z"
    })
}

fn power_label(
    base: &str,
    n: usize,
    k: usize,
) -> String {
    if k == 1 {
        format!("{base}{n}")
    } else {
        format!("{base}{n}^{k}")
    }
}

/// Merged character rows: real rows, with complex-conjugate pairs added.
fn merged_rows(ct: &super::representations::CharTable) -> Vec<(usize, f64, Vec<f64>)> {
    let m = ct.z.len();
    let mut used = vec![false; m];
    let mut out = Vec::new();
    for r in 0..m {
        if used[r] {
            continue;
        }
        used[r] = true;
        let real = ct.z[r].iter().all(|z| z.im.abs() < 1e-9);
        if real {
            out.push((ct.degrees[r], 1.0, ct.z[r].iter().map(|z| z.re).collect()));
        } else {
            let partner = (0..m).find(|&j| {
                !used[j] && ct.z[j].iter().zip(&ct.z[r]).all(|(a, b)| (a - b.conj()).norm() < 1e-8)
            });
            if let Some(j) = partner {
                used[j] = true;
            }
            out.push((2 * ct.degrees[r], 2.0, ct.z[r].iter().map(|z| 2.0 * z.re).collect()));
        }
    }
    out
}

impl PointGroup {
    /// The named group in its standard orientation.
    pub(crate) fn named(name: &str) -> Option<Self> {
        let f = parse_family(name)?;
        Self::from_matrices(&family_name(f), family_matrices(f)?)
    }

    /// The group formed by `mats` (closed under multiplication).
    #[allow(clippy::too_many_lines)] // class labels, ordering and irreps in one pass
    pub(crate) fn from_matrices(
        name: &str,
        mut mats: Vec<M3>,
    ) -> Option<Self> {
        mats.sort_by(|a, b| {
            let id_a = key(a) == key(&ID3);
            let id_b = key(b) == key(&ID3);
            id_b.cmp(&id_a)
                .then_with(|| det(b).total_cmp(&det(a)))
                .then_with(|| ((trace(b) * 1e4).round()).total_cmp(&(trace(a) * 1e4).round()))
                .then_with(|| key(a).cmp(&key(b)))
        });
        let n = mats.len();
        let find = |m: &M3| {
            mats.iter().position(|x| (0..3).all(|r| (0..3).all(|c| (x[r][c] - m[r][c]).abs() < 2e-2)))
        };
        let mut table = vec![vec![0; n]; n];
        for i in 0..n {
            for j in 0..n {
                table[i][j] = find(&mm(&mats[i], &mats[j]))?;
            }
        }
        let tab = Tab::new(&table)?;
        let ct = character_table(&tab)?;
        let kinds: Vec<Kind> = mats.iter().map(classify).collect();
        let geo = Geo::new(&kinds);
        let nclass = ct.classes.len();
        let rep_kind: Vec<&Kind> = ct.reps.iter().map(|&r| &kinds[r]).collect();

        // Classes of the proper C2 perpendicular to the principal axis, and
        // of the vertical mirrors.
        let is_perp_c2 = |c: usize| matches!(rep_kind[c], Kind::Rot { axis, n: 2, .. } if geo.perpendicular(*axis));
        let is_vert_mirror =
            |c: usize| matches!(rep_kind[c], Kind::Mirror { normal } if geo.perpendicular(*normal));
        let perp_classes: Vec<usize> = (0..nclass).filter(|&c| is_perp_c2(c)).collect();
        let vert_classes: Vec<usize> = (0..nclass).filter(|&c| is_vert_mirror(c)).collect();
        let has_s2n = kinds.iter().any(|k| matches!(k, Kind::Imp { n, .. } if *n == 2 * geo.nmax));
        let has_c4 = kinds.iter().any(|k| matches!(k, Kind::Rot { n: 4, .. }));
        let class_has = |c: usize, pred: &dyn Fn(&Kind) -> bool| ct.classes[c].iter().any(|&e| pred(&kinds[e]));
        let contains_axis = |c: usize, ax: [f64; 3]| {
            class_has(c, &|k| match k {
                | Kind::Rot { axis, .. } => parallel(*axis, ax),
                | Kind::Mirror { normal } => parallel(*normal, ax),
                | _ => false,
            })
        };

        // The power in a class label: the smaller of `k` and `n - k`, unless the
        // class is not closed under inversion and comes second.
        let label_power = |n: usize, k: usize, c: usize| {
            let inv = ct.inv_class[c];
            if inv == c || c < inv { k } else { n - k }
        };
        // Base labels of the classes.
        let mut base = vec![String::new(); nclass];
        for c in 0..nclass {
            base[c] = match rep_kind[c] {
                | Kind::E => "E".into(),
                | Kind::Inv => "i".into(),
                | Kind::Imp { n, k, .. } => power_label("S", *n, label_power(*n, *k, c)),
                | Kind::Rot { axis, n, k } => {
                    if *n == 2 && geo.principal.is_none() && !geo.polyhedral && geo.nmax == 2 {
                        match axis_name(*axis) {
                            | Some(a) => format!("C2({a})"),
                            | None => "C2".into(),
                        }
                    } else if *n == 2 && geo.perpendicular(*axis) {
                        if perp_classes.len() >= 2 {
                            let first = perp_classes
                                .iter()
                                .copied()
                                .find(|&p| contains_axis(p, [1.0, 0.0, 0.0]))
                                .unwrap_or(perp_classes[0]);
                            if c == first { "C2'".into() } else { "C2''".into() }
                        } else if has_s2n {
                            "C2'".into()
                        } else {
                            "C2".into()
                        }
                    } else if geo.polyhedral && *n == 2 && has_c4 && ct.classes[c].len() > 3 {
                        "C2'".into()
                    } else {
                        power_label("C", *n, label_power(*n, *k, c))
                    }
                },
                | Kind::Mirror { normal } => {
                    if geo.polyhedral {
                        let al = ct.classes[c]
                            .iter()
                            .filter(|&&e| matches!(&kinds[e], Kind::Mirror { normal } if aligned(*normal)))
                            .count();
                        if al == ct.classes[c].len() {
                            "σh".into()
                        } else if al == 0 {
                            "σd".into()
                        } else {
                            "σ".into()
                        }
                    } else if geo.principal.is_none() {
                        match (axis_name(*normal), nclass) {
                            | (_, 2) => "σh".into(),
                            | (Some(a), _) => format!("σ({})", match a {
                                | "x" => "yz",
                                | "y" => "xz",
                                | _ => "xy",
                            }),
                            | (None, _) => "σ".into(),
                        }
                    } else if geo.along(*normal) {
                        "σh".into()
                    } else if has_s2n {
                        "σd".into()
                    } else if geo.nmax == 2 && perp_classes.is_empty() && vert_classes.len() == 2 {
                        // C2v
                        let first = vert_classes
                            .iter()
                            .copied()
                            .find(|&v| class_has(v, &|k| matches!(k, Kind::Mirror { normal } if axis_name(*normal) == Some("y"))))
                            .unwrap_or(vert_classes[0]);
                        if c == first { "σv(xz)".into() } else { "σv'(yz)".into() }
                    } else if vert_classes.len() >= 2 {
                        let contains_x = ct.classes[c]
                            .iter()
                            .any(|&e| matches!(&kinds[e], Kind::Mirror { normal } if normal[0].abs() < 1e-6));
                        if contains_x { "σv".into() } else { "σd".into() }
                    } else {
                        "σv".into()
                    }
                },
            };
        }
        // Make the labels unique.
        for c in 0..nclass {
            let dup: Vec<usize> = (0..nclass).filter(|&d| base[d] == base[c]).collect();
            if dup.len() > 1 {
                let pos = dup.iter().position(|&d| d == c).unwrap_or(0);
                base[c] = format!("{}_{}", base[c], pos + 1);
            }
        }

        // Display order of the classes.
        let rank = |c: usize| match rep_kind[c] {
            | Kind::E => 0,
            | Kind::Rot { .. } => 1,
            | Kind::Inv => 2,
            | Kind::Imp { .. } => 3,
            | Kind::Mirror { .. } => 4,
        };
        let sort_n = |c: usize| match rep_kind[c] {
            | Kind::Rot { n, .. } | Kind::Imp { n, .. } => *n,
            | _ => 0,
        };
        let sort_k = |c: usize| match rep_kind[c] {
            | Kind::Rot { k, .. } | Kind::Imp { k, .. } => *k,
            | _ => 0,
        };
        let perp = |c: usize| match rep_kind[c] {
            | Kind::Rot { axis, .. } => !geo.along(*axis),
            | _ => false,
        };
        let mut order: Vec<usize> = (0..nclass).collect();
        order.sort_by(|&a, &b| {
            rank(a)
                .cmp(&rank(b))
                .then_with(|| sort_n(b).cmp(&sort_n(a)))
                .then_with(|| perp(a).cmp(&perp(b)))
                .then_with(|| sort_k(a).cmp(&sort_k(b)))
                .then_with(|| base[a].cmp(&base[b]))
        });
        let class_elems: Vec<Vec<usize>> = order.iter().map(|&c| ct.classes[c].clone()).collect();
        let class_labels: Vec<String> = order
            .iter()
            .map(|&c| {
                let size = ct.classes[c].len();
                if size > 1 { format!("{size}{}", base[c]) } else { base[c].clone() }
            })
            .collect();
        let mut class_of = vec![0; n];
        for (pos, c) in class_elems.iter().enumerate() {
            for &e in c {
                class_of[e] = pos;
            }
        }
        // Reference classes for the labels.
        let inv_class = (0..nclass).find(|&c| matches!(rep_kind[c], Kind::Inv));
        let sigma_h = (0..nclass)
            .find(|&c| matches!(rep_kind[c], Kind::Mirror { normal } if geo.along(*normal)))
            .or_else(|| {
                (geo.principal.is_none() && !geo.polyhedral && nclass == 2)
                    .then(|| (0..nclass).find(|&c| matches!(rep_kind[c], Kind::Mirror { .. })))
                    .flatten()
            });
        // The generator of the principal axis: S_2n for D_nd and S_2n, else
        // the smallest proper rotation.
        let s2n_gen = has_s2n && geo.nmax.is_multiple_of(2);
        let gen_class = geo.principal.and_then(|p| {
            (0..nclass).find(|&c| match rep_kind[c] {
                | Kind::Imp { axis, n, k } if s2n_gen => parallel(*axis, p) && *n == 2 * geo.nmax && *k == 1,
                | Kind::Rot { axis, n, k } if !s2n_gen => parallel(*axis, p) && *n == geo.nmax && *k == 1,
                | _ => false,
            })
        });
        let gen_order = if s2n_gen { 2 * geo.nmax } else { geo.nmax };
        let c4_class = (0..nclass).find(|&c| matches!(rep_kind[c], Kind::Rot { n: 4, k: 1, .. }));
        let s4_class = (0..nclass).find(|&c| matches!(rep_kind[c], Kind::Imp { n: 4, k: 1, .. }));
        let c5_class = (0..nclass).find(|&c| matches!(rep_kind[c], Kind::Rot { n: 5, k: 1, .. }));
        let sigma_d_class = (0..nclass).find(|&c| base[c] == "σd");
        // The reference operation for subscripts 1 and 2.
        let ref_class = if geo.polyhedral {
            c4_class.or(sigma_d_class)
        } else if !perp_classes.is_empty() {
            Some(
                perp_classes
                    .iter()
                    .copied()
                    .find(|&p| contains_axis(p, [1.0, 0.0, 0.0]))
                    .unwrap_or(perp_classes[0]),
            )
        } else if !vert_classes.is_empty() {
            Some(
                vert_classes
                    .iter()
                    .copied()
                    .find(|&v| class_has(v, &|k| matches!(k, Kind::Mirror { normal } if axis_name(*normal) == Some("y"))))
                    .unwrap_or(vert_classes[0]),
            )
        } else {
            None
        };
        let d2_classes: Vec<(usize, &str)> = if geo.principal.is_none() && !geo.polyhedral && geo.nmax == 2 {
            (0..nclass)
                .filter_map(|c| match rep_kind[c] {
                    | Kind::Rot { axis, n: 2, .. } => axis_name(*axis).map(|a| (c, a)),
                    | _ => None,
                })
                .collect()
        } else {
            Vec::new()
        };

        // Preliminary labels of the merged rows.
        let rows = merged_rows(&ct);
        let mut pre: Vec<Pre> = Vec::new();
        for (dim, nm, chi) in rows {
            let d = dim as f64;
            let sign = |c: Option<usize>| c.map(|c| chi[c] / d > 0.0);
            let suffix = if geo.inversion {
                if sign(inv_class).unwrap_or(true) { "g" } else { "u" }.to_string()
            } else if sign(sigma_h).is_some() && !geo.polyhedral && geo.principal.is_some() {
                if sign(sigma_h).unwrap_or(true) { "'" } else { "''" }.to_string()
            } else if sigma_h.is_some() && geo.principal.is_none() && nclass == 2 {
                // C_s
                if sign(sigma_h).unwrap_or(true) { "'" } else { "''" }.to_string()
            } else {
                String::new()
            };
            let letter = match dim {
                | 1 => {
                    let plus = if geo.polyhedral {
                        true
                    } else if let Some(g) = gen_class {
                        chi[g] > 0.0
                    } else if !d2_classes.is_empty() {
                        d2_classes.iter().all(|&(c, _)| chi[c] > 0.0)
                    } else {
                        true
                    };
                    if plus { 'A' } else { 'B' }
                },
                | 2 => 'E',
                | 3 => 'T',
                | 4 => 'G',
                | 5 => 'H',
                | _ => 'I',
            };
            pre.push(Pre { dim, norm: nm, chi, letter, suffix, sub: None });
        }
        // Subscripts where the letter and parity do not identify the row.
        let groups: Vec<(char, String)> = pre.iter().map(|p| (p.letter, p.suffix.clone())).collect();
        for i in 0..pre.len() {
            let same = groups.iter().filter(|g| **g == groups[i]).count();
            if same < 2 {
                continue;
            }
            let p = &pre[i];
            let sub = match p.letter {
                | 'A' | 'B' => {
                    if !d2_classes.is_empty() && p.letter == 'B' {
                        d2_classes.iter().find(|&&(c, _)| p.chi[c] > 0.0).map(|&(_, a)| match a {
                            | "z" => 1,
                            | "y" => 2,
                            | _ => 3,
                        })
                    } else {
                        ref_class.map(|r| if p.chi[r] > 0.0 { 1 } else { 2 })
                    }
                },
                | 'E' => {
                    if geo.polyhedral {
                        None
                    } else {
                        gen_class.and_then(|g| {
                            (1..=gen_order / 2)
                                .find(|&m| (p.chi[g] - 2.0 * (2.0 * PI * m as f64 / gen_order as f64).cos()).abs() < 1e-6)
                        })
                    }
                },
                | 'T' => {
                    let r = c4_class.or(if has_c4 { None } else { s4_class }).or(c5_class);
                    r.map(|r| if p.chi[r] > 0.0 { 1 } else { 2 })
                },
                | _ => None,
            };
            pre[i].sub = sub;
        }
        let mut irreps: Vec<Irrep> = pre
            .iter()
            .map(|p| {
                let sub = p.sub.map_or_else(String::new, |s| s.to_string());
                Irrep {
                    label: format!("{}{}{}", p.letter, sub, p.suffix),
                    dim: p.dim,
                    chi: order.iter().map(|&c| p.chi[c]).collect(),
                    norm: p.norm,
                }
            })
            .collect();
        // Sort: unprimed/g before double primed/u, then A, B, E, T, G, H, then the subscript.
        let sort_key = |ir: &Irrep| {
            let l = ir.label.chars().next().unwrap_or('A');
            let letter = "ABETGHI".find(l).unwrap_or(9);
            let tail: String = ir.label.chars().skip(1).collect();
            let sub: usize = tail.chars().take_while(char::is_ascii_digit).collect::<String>().parse().unwrap_or(0);
            let suffix: String = tail.chars().skip_while(char::is_ascii_digit).collect();
            let parity = match suffix.as_str() {
                | "g" | "'" | "" => 0,
                | _ => 1,
            };
            (parity, letter, sub)
        };
        irreps.sort_by_key(|ir| sort_key(ir));
        // Unique labels.
        for i in 0..irreps.len() {
            let dup: Vec<usize> = (0..irreps.len()).filter(|&j| irreps[j].label == irreps[i].label).collect();
            if dup.len() > 1 {
                let pos = dup.iter().position(|&d| d == i).unwrap_or(0);
                irreps[i].label = format!("{}_{}", irreps[i].label, pos + 1);
            }
        }
        Some(Self { name: name.to_string(), mats, tab, class_elems, class_labels, irreps, class_of })
    }

    pub(crate) fn class_sizes(&self) -> Vec<usize> {
        self.class_elems.iter().map(Vec::len).collect()
    }

    /// The multiplicities of the irreps in a class function (display order).
    pub(crate) fn multiplicities(
        &self,
        chi: &[f64],
    ) -> Option<Vec<i64>> {
        let n = self.mats.len() as f64;
        let sizes = self.class_sizes();
        let mut out = Vec::new();
        for ir in &self.irreps {
            let ip: f64 = (0..sizes.len()).map(|c| sizes[c] as f64 * chi[c] * ir.chi[c]).sum::<f64>() / n;
            let m = ip / ir.norm;
            if (m - m.round()).abs() > 1e-6 {
                return None;
            }
            out.push(m.round() as i64);
        }
        Some(out)
    }

    /// The class function of the element-wise values `f`.
    pub(crate) fn class_function(
        &self,
        f: impl Fn(&M3) -> f64,
    ) -> Vec<f64> {
        self.class_elems.iter().map(|c| f(&self.mats[c[0]])).collect()
    }

    pub(crate) fn vector_character(&self) -> Vec<f64> {
        self.class_function(trace)
    }

    pub(crate) fn rotation_character(&self) -> Vec<f64> {
        self.class_function(|m| det(m) * trace(m))
    }

    pub(crate) fn quadratic_character(&self) -> Vec<f64> {
        self.class_function(|m| {
            let sq = mm(m, m);
            (trace(m) * trace(m) + trace(&sq)) / 2.0
        })
    }

    pub(crate) fn ir_active(&self) -> Vec<bool> {
        let m = self.multiplicities(&self.vector_character()).unwrap_or_default();
        m.iter().map(|&x| x > 0).collect()
    }

    pub(crate) fn raman_active(&self) -> Vec<bool> {
        let m = self.multiplicities(&self.quadratic_character()).unwrap_or_default();
        m.iter().map(|&x| x > 0).collect()
    }
}

// ----------------------------------------------------------------------
// Exact numbers
// ----------------------------------------------------------------------

fn recognise_quadratic(x: f64) -> Option<(i64, i64, i64, i64)> {
    use num_integer::Integer;
    for c in [1_i64, 2, 3, 4, 6, 8, 12] {
        for d in [1_i64, 2, 3, 5, 6, 7, 10, 11, 13, 15] {
            let sq = (d as f64).sqrt();
            let max_b = if d == 1 { 0 } else { 24 };
            for bb in 0..=max_b {
                for b in if bb == 0 { vec![0] } else { vec![bb, -bb] } {
                    let a = (x * c as f64 - b as f64 * sq).round();
                    if (a + b as f64 * sq - x * c as f64).abs() < 1e-9 && a.abs() < 1e6 {
                        let a = a as i64;
                        let g = a.gcd(&b).gcd(&c);
                        let g = if g == 0 { 1 } else { g };
                        return Some((a / g, b / g, d, c / g));
                    }
                }
            }
        }
    }
    None
}

fn rat(
    g: &mut Graph,
    p: i64,
    q: i64,
) -> NodeId {
    match Number::fraction(p, q) {
        | Some(n) => g.num(n),
        | None => g.int(0),
    }
}

/// An exact term for a real number that is an integer, a quadratic
/// irrational, or the sine or cosine of a rational multiple of `pi`.
pub(crate) fn real_node(
    g: &mut Graph,
    x: f64,
) -> NodeId {
    if (x - x.round()).abs() < 1e-9 {
        return g.int(x.round() as i64);
    }
    if let Some((a, b, d, c)) = recognise_quadratic(x) {
        let mut terms = vec![rat(g, a, c)];
        if b != 0 {
            let dn = g.int(d);
            if let Some(sq) = apply(g, "sqrt", &[dn]) {
                let coef = rat(g, b, c);
                terms.push(prod(g, &[coef, sq]));
            }
        }
        return sum(g, &terms);
    }
    for n in 3..=24_i64 {
        for k in 1..n {
            for (name, val) in [("cos", (2.0 * PI * k as f64 / n as f64).cos()), ("sin", (2.0 * PI * k as f64 / n as f64).sin())] {
                for sign in [1.0, -1.0] {
                    if (x - sign * val).abs() < 1e-9 {
                        let q = rat(g, 2 * k, n);
                        if let Some(pi) = apply(g, "pi", &[]) {
                            let arg = prod(g, &[q, pi]);
                            if let Some(t) = apply(g, name, &[arg]) {
                                let s = g.int(if sign > 0.0 { 1 } else { -1 });
                                return prod(g, &[s, t]);
                            }
                        }
                    }
                }
            }
        }
    }
    g.float(x)
}

fn matrix_nodes(
    g: &mut Graph,
    m: &M3,
) -> NodeId {
    let rows: Vec<Vec<NodeId>> = m.iter().map(|r| r.iter().map(|&x| real_node(g, x)).collect()).collect();
    matrix(g, &rows)
}

// ----------------------------------------------------------------------
// Crystallography data
// ----------------------------------------------------------------------

/// Schoenflies, Hermann–Mauguin, crystal system of the 32 crystal classes.
const CRYSTAL_CLASSES: [(&str, &str, &str); 32] = [
    ("C1", "1", "triclinic"),
    ("Ci", "-1", "triclinic"),
    ("C2", "2", "monoclinic"),
    ("Cs", "m", "monoclinic"),
    ("C2h", "2/m", "monoclinic"),
    ("D2", "222", "orthorhombic"),
    ("C2v", "mm2", "orthorhombic"),
    ("D2h", "mmm", "orthorhombic"),
    ("C4", "4", "tetragonal"),
    ("S4", "-4", "tetragonal"),
    ("C4h", "4/m", "tetragonal"),
    ("D4", "422", "tetragonal"),
    ("C4v", "4mm", "tetragonal"),
    ("D2d", "-42m", "tetragonal"),
    ("D4h", "4/mmm", "tetragonal"),
    ("C3", "3", "trigonal"),
    ("S6", "-3", "trigonal"),
    ("D3", "32", "trigonal"),
    ("C3v", "3m", "trigonal"),
    ("D3d", "-3m", "trigonal"),
    ("C6", "6", "hexagonal"),
    ("C3h", "-6", "hexagonal"),
    ("C6h", "6/m", "hexagonal"),
    ("D6", "622", "hexagonal"),
    ("C6v", "6mm", "hexagonal"),
    ("D3h", "-6m2", "hexagonal"),
    ("D6h", "6/mmm", "hexagonal"),
    ("T", "23", "cubic"),
    ("Th", "m-3", "cubic"),
    ("O", "432", "cubic"),
    ("Td", "-43m", "cubic"),
    ("Oh", "m-3m", "cubic"),
];

/// The 14 Bravais lattices: Pearson symbol, system, centring, holohedry.
const BRAVAIS: [(&str, &str, &str, &str); 14] = [
    ("aP", "triclinic", "P", "Ci"),
    ("mP", "monoclinic", "P", "C2h"),
    ("mC", "monoclinic", "C", "C2h"),
    ("oP", "orthorhombic", "P", "D2h"),
    ("oC", "orthorhombic", "C", "D2h"),
    ("oI", "orthorhombic", "I", "D2h"),
    ("oF", "orthorhombic", "F", "D2h"),
    ("tP", "tetragonal", "P", "D4h"),
    ("tI", "tetragonal", "I", "D4h"),
    ("hR", "trigonal", "R", "D3d"),
    ("hP", "hexagonal", "P", "D6h"),
    ("cP", "cubic", "P", "Oh"),
    ("cI", "cubic", "I", "Oh"),
    ("cF", "cubic", "F", "Oh"),
];

const SYSTEMS: [&str; 7] = ["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"];

/// The key of a Hermann–Mauguin name: bars as `b`, no slashes or underscores.
fn hm_key(s: &str) -> String {
    s.trim_start_matches("hm_").chars().filter(|&c| c != '_' && c != '/').map(|c| if c == '-' { 'b' } else { c }).collect()
}

/// `phi`-based least dimension of a lattice with an `n`-fold rotation.
const fn min_dimension(n: u64) -> u64 {
    if n <= 1 {
        return 0;
    }
    if n == 2 {
        return 1;
    }
    let n = if n % 4 == 2 { n / 2 } else { n };
    let mut m = n;
    let mut total = 0;
    let mut p = 2;
    while p * p <= m {
        if m % p == 0 {
            let mut pk = 1;
            while m % p == 0 {
                m /= p;
                pk *= p;
            }
            total += pk / p * (p - 1);
        }
        p += 1;
    }
    if m > 1 {
        total += m - 1;
    }
    total
}

// ----------------------------------------------------------------------
// Molecules
// ----------------------------------------------------------------------

struct Atoms {
    pos: Vec<[f64; 3]>,
    label: Vec<u64>,
}

fn read_atoms(
    g: &Graph,
    n: NodeId,
) -> Option<Atoms> {
    let bindings = HashMap::new();
    let mut rows = Vec::new();
    let mut seen: Vec<ClassId> = Vec::new();
    for atom in items(g, n)? {
        let parts = items(g, atom)?;
        let (label, coords) = match parts.len() {
            | 3 => (0_u64, &parts[..]),
            | 4 => {
                let class = g.find(parts[0]);
                let pos = seen.iter().position(|&c| c == class).unwrap_or_else(|| {
                    seen.push(class);
                    seen.len() - 1
                });
                (pos as u64 + 1, &parts[1..])
            },
            | _ => return None,
        };
        let mut p = [0.0; 3];
        for (i, &c) in coords.iter().enumerate() {
            let z = g.eval_complex(c, &bindings)?;
            if z.im.abs() > 1e-12 {
                return None;
            }
            p[i] = z.re;
        }
        rows.push((label, p));
    }
    if rows.is_empty() || rows.len() > 80 {
        return None;
    }
    let centroid: [f64; 3] = [0, 1, 2].map(|i| rows.iter().map(|r| r.1[i]).sum::<f64>() / rows.len() as f64);
    Some(Atoms {
        pos: rows.iter().map(|r| [0, 1, 2].map(|i| r.1[i] - centroid[i])).collect(),
        label: rows.iter().map(|r| r.0).collect(),
    })
}

impl Atoms {
    fn is_linear(
        &self,
        tol: f64,
    ) -> Option<[f64; 3]> {
        let dir = self.pos.iter().copied().max_by(|a, b| norm(*a).total_cmp(&norm(*b)))?;
        if norm(dir) < tol {
            return None;
        }
        let u = unit(dir);
        self.pos.iter().all(|&p| norm(cross(p, u)) < tol).then_some(u)
    }

    fn preserved_by(
        &self,
        m: &M3,
        tol: f64,
    ) -> bool {
        let mut used = vec![false; self.pos.len()];
        for (i, &p) in self.pos.iter().enumerate() {
            let q = apply_m(m, p);
            let found = (0..self.pos.len()).find(|&j| {
                !used[j] && self.label[j] == self.label[i] && norm([q[0] - self.pos[j][0], q[1] - self.pos[j][1], q[2] - self.pos[j][2]]) < tol
            });
            match found {
                | Some(j) => used[j] = true,
                | None => return false,
            }
        }
        true
    }

    fn fixed_atoms(
        &self,
        m: &M3,
        tol: f64,
    ) -> usize {
        self.pos
            .iter()
            .filter(|&&p| {
                let q = apply_m(m, p);
                norm([q[0] - p[0], q[1] - p[1], q[2] - p[2]]) < tol
            })
            .count()
    }

    /// The symmetry operations (a group of matrices).
    fn operations(
        &self,
        tol: f64,
    ) -> Option<Vec<M3>> {
        let mut cand: HashMap<[i64; 3], [f64; 3]> = HashMap::new();
        let mut add = |v: [f64; 3]| {
            if norm(v) > 1e-6 {
                let u = canonical_axis(unit(v));
                cand.entry(u.map(|x| (x * 1e4).round() as i64)).or_insert(u);
            }
        };
        for u in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            add(u);
        }
        for (i, &p) in self.pos.iter().enumerate() {
            add(p);
            for &q in &self.pos[i + 1..] {
                add([p[0] + q[0], p[1] + q[1], p[2] + q[2]]);
                add([p[0] - q[0], p[1] - q[1], p[2] - q[2]]);
                add(cross(p, q));
            }
        }
        let mut found: Vec<M3> = vec![ID3];
        let mut push = |m: M3| {
            let same = |x: &M3| (0..3).all(|r| (0..3).all(|c| (x[r][c] - m[r][c]).abs() < 1e-2));
            if !found.iter().any(same) {
                found.push(m);
            }
        };
        let inv = mneg(&ID3);
        if self.preserved_by(&inv, tol) {
            push(inv);
        }
        for u in cand.values() {
            let sigma = refl(*u);
            if self.preserved_by(&sigma, tol) {
                push(sigma);
            }
            for n in 2..=12 {
                for k in 1..n {
                    let angle = 2.0 * PI * f64::from(k) / f64::from(n);
                    let r = rot(*u, angle);
                    if self.preserved_by(&r, tol) {
                        push(r);
                    }
                    let s = mm(&sigma, &r);
                    if self.preserved_by(&s, tol) {
                        push(s);
                    }
                }
            }
        }
        if found.len() > ORDER_LIMIT {
            return None;
        }
        Some(found.into_iter().map(|m| m.map(|row| row.map(snap))).collect())
    }
}

/// Replaces `x` by a nearby simple exact value (a rational with small
/// denominator, or a quadratic irrational) if there is one within `2e-3`.
fn snap(x: f64) -> f64 {
    for (ds, bs) in [(&[1_i64][..], 0_i64), (&[2, 3, 5, 6, 7, 10][..], 12)] {
        for c in [1_i64, 2, 3, 4, 6, 8, 12] {
            for &d in ds {
                let sq = (d as f64).sqrt();
                for b in -bs..=bs {
                    let a = (x * c as f64 - b as f64 * sq).round();
                    let v = (a + b as f64 * sq) / c as f64;
                    if (v - x).abs() < 2e-3 {
                        return v;
                    }
                }
            }
        }
    }
    x
}

/// Fingerprint of a matrix group: sorted `(det sign, trace)` pairs.
fn fingerprint(mats: &[M3]) -> Vec<(i32, i64)> {
    let mut f: Vec<(i32, i64)> =
        mats.iter().map(|m| (if det(m) > 0.0 { 1 } else { -1 }, (trace(m) * 1e3).round() as i64)).collect();
    f.sort_unstable();
    f
}

/// The Schoenflies name of a group of matrices.
fn identify(mats: &[M3]) -> Option<String> {
    let fp = fingerprint(mats);
    let mut families = vec![Family::T, Family::Td, Family::Th, Family::O, Family::Oh, Family::I, Family::Ih];
    for n in 1..=12 {
        families.push(Family::C(n));
        families.push(Family::Ch(n));
        if n >= 2 {
            families.push(Family::Cv(n));
            families.push(Family::D(n));
            families.push(Family::Dh(n));
            families.push(Family::Dd(n));
        }
        if n.is_multiple_of(2) {
            families.push(Family::S(n));
        }
    }
    families.into_iter().find_map(|f| {
        let m = family_matrices(f)?;
        (m.len() == mats.len() && fingerprint(&m) == fp).then(|| family_name(f))
    })
}

enum MoleculeGroup {
    Linear(bool),
    Group(Box<PointGroup>),
}

fn molecule_group(
    atoms: &Atoms,
    tol: f64,
) -> Option<MoleculeGroup> {
    if atoms.pos.len() == 1 {
        return None;
    }
    if atoms.is_linear(tol).is_some() {
        let inversion = atoms.preserved_by(&mneg(&ID3), tol);
        return Some(MoleculeGroup::Linear(inversion));
    }
    let ops = atoms.operations(tol)?;
    let name = identify(&ops)?;
    Some(MoleculeGroup::Group(Box::new(PointGroup::from_matrices(&name, ops)?)))
}

fn tolerance(
    g: &Graph,
    a: &[NodeId],
    pos: usize,
) -> f64 {
    a.get(pos).and_then(|&t| super::float(g, t)).unwrap_or(1e-3)
}

fn molecule_of(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<(Atoms, PointGroup, f64)> {
    let tol = tolerance(cx.graph, a, 1);
    let atoms = read_atoms(cx.graph, *a.first()?)?;
    match molecule_group(&atoms, tol)? {
        | MoleculeGroup::Group(pg) => Some((atoms, *pg, tol)),
        | MoleculeGroup::Linear(_) => None,
    }
}

/// The multiplicities of `Gamma_3N`, `Gamma_trans`, `Gamma_rot`, `Gamma_vib`.
fn vibrational(
    atoms: &Atoms,
    pg: &PointGroup,
    tol: f64,
) -> Option<[Vec<i64>; 4]> {
    let g3n = pg.class_function(|m| atoms.fixed_atoms(m, tol) as f64 * trace(m));
    let trans = pg.vector_character();
    let rot = pg.rotation_character();
    let vib: Vec<f64> = (0..g3n.len()).map(|i| g3n[i] - trans[i] - rot[i]).collect();
    Some([pg.multiplicities(&g3n)?, pg.multiplicities(&trans)?, pg.multiplicities(&rot)?, pg.multiplicities(&vib)?])
}

fn molecule_decomposition(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (atoms, pg, tol) = molecule_of(cx, a)?;
    let [g3n, trans, rot, vib] = vibrational(&atoms, &pg, tol)?;
    let mut rows = Vec::new();
    for (i, ir) in pg.irreps.iter().enumerate() {
        let label = cx.graph.sym(&ir.label);
        rows.push(V::List(vec![
            V::Node(label),
            V::int(g3n[i]),
            V::int(trans[i]),
            V::int(rot[i]),
            V::int(vib[i]),
        ]));
    }
    Some(V::List(rows))
}

fn sum_of_irreps(
    g: &mut Graph,
    pg: &PointGroup,
    mult: &[i64],
) -> NodeId {
    let mut terms = Vec::new();
    for (ir, &m) in pg.irreps.iter().zip(mult) {
        if m != 0 {
            let s = g.sym(&ir.label);
            let c = g.int(m);
            terms.push(prod(g, &[c, s]));
        }
    }
    sum(g, &terms)
}

fn molecule_vibrations(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (atoms, pg, tol) = molecule_of(cx, a)?;
    let [_, _, _, vib] = vibrational(&atoms, &pg, tol)?;
    Some(V::Node(sum_of_irreps(cx.graph, &pg, &vib)))
}

fn molecule_vibrations_table(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (atoms, pg, tol) = molecule_of(cx, a)?;
    let [_, _, _, vib] = vibrational(&atoms, &pg, tol)?;
    let mut rows = Vec::new();
    for (ir, &m) in pg.irreps.iter().zip(&vib) {
        if m != 0 {
            let s = cx.graph.sym(&ir.label);
            rows.push(V::List(vec![V::Node(s), V::int(m)]));
        }
    }
    Some(V::List(rows))
}

fn molecule_spectroscopy(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let (atoms, pg, tol) = molecule_of(cx, a)?;
    let [_, _, _, vib] = vibrational(&atoms, &pg, tol)?;
    let (ir, raman) = (pg.ir_active(), pg.raman_active());
    let mut rows = Vec::new();
    for (i, irrep) in pg.irreps.iter().enumerate() {
        if vib[i] != 0 {
            let s = cx.graph.sym(&irrep.label);
            rows.push(V::List(vec![V::Node(s), V::int(vib[i]), V::Bool(ir[i]), V::Bool(raman[i])]));
        }
    }
    Some(V::List(rows))
}

fn molecule_symmetry_operations(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let tol = tolerance(cx.graph, a, 1);
    let atoms = read_atoms(cx.graph, *a.first()?)?;
    let ops = atoms.operations(tol)?;
    let pg = PointGroup::from_matrices("", ops).map(|p| p.mats)?;
    Some(V::List(pg.iter().map(|m| V::Node(matrix_nodes(cx.graph, m))).collect()))
}

fn molecule_point_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let tol = tolerance(cx.graph, a, 1);
    let atoms = read_atoms(cx.graph, *a.first()?)?;
    let name = match molecule_group(&atoms, tol)? {
        | MoleculeGroup::Linear(true) => "Dinfh".to_string(),
        | MoleculeGroup::Linear(false) => "Cinfv".to_string(),
        | MoleculeGroup::Group(pg) => pg.name.clone(),
    };
    Some(V::Node(cx.graph.sym(&name)))
}

// ----------------------------------------------------------------------
// Operators on named groups
// ----------------------------------------------------------------------

fn read_name(
    g: &Graph,
    n: NodeId,
) -> Option<String> {
    let name_of = |e: NodeId| -> Option<String> {
        if let Some(s) = g.as_symbol(e) {
            return Some(g.interner().symbol_name(s).to_string());
        }
        if g.children(e).is_empty() && &*g.ops().get(g.op(e)).name == "I" {
            return Some("I".into());
        }
        None
    };
    name_of(n).or_else(|| g.enodes(g.find(n)).find_map(name_of))
}

fn group_arg(
    cx: &Cx<'_>,
    a: &[NodeId],
) -> Option<PointGroup> {
    PointGroup::named(&read_name(cx.graph, *a.first()?)?)
}

fn chi_arg(
    cx: &Cx<'_>,
    pg: &PointGroup,
    n: NodeId,
) -> Option<Vec<f64>> {
    let bindings = HashMap::new();
    let vals: Vec<f64> = items(cx.graph, n)?
        .into_iter()
        .map(|x| {
            cx.graph.eval_complex(x, &bindings).and_then(|z: Complex64| (z.im.abs() < 1e-9).then_some(z.re))
        })
        .collect::<Option<_>>()?;
    if vals.len() == pg.class_elems.len() {
        Some(vals)
    } else if vals.len() == pg.mats.len() {
        Some(pg.class_elems.iter().map(|c| vals[c[0]]).collect())
    } else {
        None
    }
}

fn point_group(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let n = pg.mats.len();
    let elems: Vec<NodeId> = pg.mats.iter().map(|m| matrix_nodes(cx.graph, m)).collect();
    let table: Vec<Vec<usize>> = (0..n).map(|i| (0..n).map(|j| pg.tab.m(i, j)).collect()).collect();
    build_group(cx.graph, &elems, &table)
}

fn point_group_order(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    Some(V::uint(group_arg(cx, a)?.mats.len()))
}

fn point_group_is_crystallographic(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(V::Bool((0..pg.mats.len()).all(|x| matches!(pg.tab.order(x), 1 | 2 | 3 | 4 | 6))))
}

fn point_group_classes(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(V::List(pg.class_labels.iter().map(|l| V::Node(cx.graph.sym(l))).collect()))
}

fn point_group_class_sizes(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(V::List(pg.class_sizes().into_iter().map(V::uint).collect()))
}

fn point_group_irreps(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(V::List(pg.irreps.iter().map(|i| V::Node(cx.graph.sym(&i.label))).collect()))
}

fn point_group_irrep_dimensions(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(V::List(pg.irreps.iter().map(|i| V::uint(i.dim)).collect()))
}

fn character_rows(
    g: &mut Graph,
    rows: &[Vec<f64>],
) -> V {
    V::List(rows.iter().map(|r| V::List(r.iter().map(|&x| V::Node(real_node(g, x))).collect())).collect())
}

fn point_group_character_table(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let rows: Vec<Vec<f64>> = pg.irreps.iter().map(|i| i.chi.clone()).collect();
    Some(character_rows(cx.graph, &rows))
}

fn point_group_multiplicities(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let chi = chi_arg(cx, &pg, *a.get(1)?)?;
    Some(V::ints(pg.multiplicities(&chi)?))
}

fn point_group_decompose(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let chi = chi_arg(cx, &pg, *a.get(1)?)?;
    let m = pg.multiplicities(&chi)?;
    Some(V::Node(sum_of_irreps(cx.graph, &pg, &m)))
}

fn point_group_vector_character(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(character_rows(cx.graph, &[pg.vector_character()]).first_row())
}

fn point_group_rotation_character(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    Some(character_rows(cx.graph, &[pg.rotation_character()]).first_row())
}

impl V {
    fn first_row(self) -> Self {
        match self {
            | Self::List(mut rows) if !rows.is_empty() => rows.swap_remove(0),
            | other => other,
        }
    }
}

fn labels_of(
    cx: &mut Cx<'_>,
    pg: &PointGroup,
    mask: &[bool],
) -> V {
    V::List(
        pg.irreps
            .iter()
            .zip(mask)
            .filter(|(_, m)| **m)
            .map(|(i, _)| V::Node(cx.graph.sym(&i.label)))
            .collect(),
    )
}

fn point_group_ir_active(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let mask = pg.ir_active();
    Some(labels_of(cx, &pg, &mask))
}

fn point_group_raman_active(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let mask = pg.raman_active();
    Some(labels_of(cx, &pg, &mask))
}

/// The action of `m` on the coefficient vector of a polynomial of degree
/// `deg` in `x, y, z` (`deg` 1: `x, y, z`; `deg` 2: `x2, y2, z2, xy, xz, yz`).
fn function_action(
    m: &M3,
    deg: usize,
) -> Vec<Vec<f64>> {
    if deg == 1 {
        return m.iter().map(|r| r.to_vec()).collect();
    }
    let monos: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let mut out = vec![vec![0.0; 6]; 6];
    for (col, &(i, j)) in monos.iter().enumerate() {
        // m.f where f = x_i x_j gives sum over k, l of m[k][i] m[l][j] x_k x_l
        for k in 0..3 {
            for l in 0..3 {
                let c = m[k][i] * m[l][j];
                let (a, b) = (k.min(l), k.max(l));
                let row = monos.iter().position(|&p| p == (a, b)).unwrap_or(0);
                out[row][col] += c;
            }
        }
    }
    out
}

fn point_group_function_irreps(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let pg = group_arg(cx, a)?;
    let funcs: [(&str, usize, usize); 9] = [
        ("x", 1, 0),
        ("y", 1, 1),
        ("z", 1, 2),
        ("x^2", 2, 0),
        ("y^2", 2, 1),
        ("z^2", 2, 2),
        ("x*y", 2, 3),
        ("x*z", 2, 4),
        ("y*z", 2, 5),
    ];
    let mut rows = Vec::new();
    for (text, deg, col) in funcs {
        let node = cx.graph.parse(text).ok()?;
        let dim = if deg == 1 { 3 } else { 6 };
        let mut irreps = Vec::new();
        for ir in &pg.irreps {
            let mut v = vec![0.0; dim];
            for (e, m) in pg.mats.iter().enumerate() {
                let rho = function_action(m, deg);
                let chi = ir.chi[pg.class_of[e]];
                for r in 0..dim {
                    v[r] += chi * rho[r][col];
                }
            }
            if v.iter().any(|x| x.abs() > 1e-8) {
                irreps.push(V::Node(cx.graph.sym(&ir.label)));
            }
        }
        rows.push(V::List(vec![V::Node(node), V::List(irreps)]));
    }
    Some(V::List(rows))
}

fn point_group_hm(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let name = read_name(cx.graph, *a.first()?)?;
    let f = parse_family(&name).map(family_name)?;
    let entry = CRYSTAL_CLASSES.iter().find(|c| c.0 == f)?;
    Some(V::Node(cx.graph.sym(entry.1)))
}

fn point_group_from_hm(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = *a.first()?;
    let key = if let Some(v) = super::small(cx.graph, n) {
        if v < 0 { format!("b{}", -v) } else { v.to_string() }
    } else {
        hm_key(&read_name(cx.graph, n)?)
    };
    let entry = CRYSTAL_CLASSES.iter().find(|c| hm_key(c.1) == key)?;
    Some(V::Node(cx.graph.sym(entry.0)))
}

fn point_group_crystal_system(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let name = read_name(cx.graph, *a.first()?)?;
    let f = parse_family(&name).map(family_name)?;
    let entry = CRYSTAL_CLASSES.iter().find(|c| c.0 == f)?;
    Some(V::Node(cx.graph.sym(entry.2)))
}

#[allow(clippy::unnecessary_wraps)] // kernel signature
fn crystallographic_point_groups(
    cx: &mut Cx<'_>,
    _a: &[NodeId],
) -> Option<V> {
    Some(V::List(CRYSTAL_CLASSES.iter().map(|c| V::Node(cx.graph.sym(c.0))).collect()))
}

#[allow(clippy::unnecessary_wraps)] // kernel signature
fn crystal_systems(
    cx: &mut Cx<'_>,
    _a: &[NodeId],
) -> Option<V> {
    let mut rows = Vec::new();
    for system in SYSTEMS {
        let holohedry = BRAVAIS.iter().find(|b| b.1 == system).map_or("Ci", |b| b.3);
        let groups: Vec<V> =
            CRYSTAL_CLASSES.iter().filter(|c| c.2 == system).map(|c| V::Node(cx.graph.sym(c.0))).collect();
        let lattices: Vec<V> =
            BRAVAIS.iter().filter(|b| b.1 == system).map(|b| V::Node(cx.graph.sym(b.0))).collect();
        rows.push(V::List(vec![
            V::Node(cx.graph.sym(system)),
            V::Node(cx.graph.sym(holohedry)),
            V::List(groups),
            V::List(lattices),
        ]));
    }
    Some(V::List(rows))
}

#[allow(clippy::unnecessary_wraps)] // kernel signature
fn bravais_lattices(
    cx: &mut Cx<'_>,
    _a: &[NodeId],
) -> Option<V> {
    Some(V::List(
        BRAVAIS
            .iter()
            .map(|b| V::List([b.0, b.1, b.2, b.3].iter().map(|s| V::Node(cx.graph.sym(s))).collect()))
            .collect(),
    ))
}

fn crystallographic_restriction(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = super::big(cx.graph, *a.first()?)?;
    let n = u64::try_from(n).ok().filter(|&n| n >= 1)?;
    let dim = match a.get(1) {
        | Some(&d) => super::idx(cx.graph, d)? as u64,
        | None => 3,
    };
    Some(V::Bool(min_dimension(n) <= dim))
}

fn crystallographic_min_dimension(
    cx: &mut Cx<'_>,
    a: &[NodeId],
) -> Option<V> {
    let n = super::big(cx.graph, *a.first()?)?;
    let n = u64::try_from(n).ok().filter(|&n| n >= 1)?;
    Some(V::uint(min_dimension(n) as usize))
}

/// Registers an operator taking a list of atoms and an optional tolerance.
fn def_molecule(
    i: &mut Installer<'_>,
    name: &str,
    run: super::Run,
) -> Result<(), RuleError> {
    def(i, name, Arity::Variadic, run)?;
    Ok(())
}

pub(crate) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    def(i, "point_group", Arity::Fixed(1), point_group)?;
    def(i, "point_group_order", Arity::Fixed(1), point_group_order)?;
    def(i, "point_group_is_crystallographic", Arity::Fixed(1), point_group_is_crystallographic)?;
    def(i, "point_group_classes", Arity::Fixed(1), point_group_classes)?;
    def(i, "point_group_class_sizes", Arity::Fixed(1), point_group_class_sizes)?;
    def(i, "point_group_irreps", Arity::Fixed(1), point_group_irreps)?;
    def(i, "point_group_irrep_dimensions", Arity::Fixed(1), point_group_irrep_dimensions)?;
    def(i, "point_group_character_table", Arity::Fixed(1), point_group_character_table)?;
    def(i, "point_group_decompose", Arity::Fixed(2), point_group_decompose)?;
    def(i, "point_group_multiplicities", Arity::Fixed(2), point_group_multiplicities)?;
    def(i, "point_group_vector_character", Arity::Fixed(1), point_group_vector_character)?;
    def(i, "point_group_rotation_character", Arity::Fixed(1), point_group_rotation_character)?;
    def(i, "point_group_ir_active", Arity::Fixed(1), point_group_ir_active)?;
    def(i, "point_group_raman_active", Arity::Fixed(1), point_group_raman_active)?;
    def(i, "point_group_function_irreps", Arity::Fixed(1), point_group_function_irreps)?;
    def(i, "point_group_hm", Arity::Fixed(1), point_group_hm)?;
    def(i, "point_group_from_hm", Arity::Fixed(1), point_group_from_hm)?;
    def(i, "point_group_crystal_system", Arity::Fixed(1), point_group_crystal_system)?;
    def(i, "crystallographic_point_groups", Arity::Fixed(0), crystallographic_point_groups)?;
    def(i, "crystal_systems", Arity::Fixed(0), crystal_systems)?;
    def(i, "bravais_lattices", Arity::Fixed(0), bravais_lattices)?;
    def(i, "crystallographic_restriction", Arity::Variadic, crystallographic_restriction)?;
    def(i, "crystallographic_min_dimension", Arity::Fixed(1), crystallographic_min_dimension)?;
    def_molecule(i, "molecule_symmetry_operations", molecule_symmetry_operations)?;
    def_molecule(i, "molecule_point_group", molecule_point_group)?;
    def_molecule(i, "molecule_decomposition", molecule_decomposition)?;
    def_molecule(i, "molecule_vibrations", molecule_vibrations)?;
    def_molecule(i, "molecule_vibrations_table", molecule_vibrations_table)?;
    def_molecule(i, "molecule_spectroscopy", molecule_spectroscopy)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::test_util::s;
    use super::*;

    /// The character of `irrep` on the class labelled `class`.
    fn chi(
        pg: &PointGroup,
        irrep: &str,
        class: &str,
    ) -> f64 {
        let ir = pg.irreps.iter().find(|i| i.label == irrep).unwrap_or_else(|| panic!("no irrep {irrep}: {:?}", labels(pg)));
        let c = pg.class_labels.iter().position(|l| l == class).unwrap_or_else(|| panic!("no class {class}: {:?}", pg.class_labels));
        ir.chi[c]
    }

    fn labels(pg: &PointGroup) -> Vec<&str> {
        pg.irreps.iter().map(|i| i.label.as_str()).collect()
    }

    fn check_table(
        pg: &PointGroup,
        classes: &[&str],
        rows: &[(&str, &[f64])],
    ) {
        assert_eq!(pg.class_labels, classes, "classes of {}", pg.name);
        assert_eq!(pg.irreps.len(), rows.len(), "{:?}", labels(pg));
        for (label, vals) in rows {
            for (c, v) in classes.iter().zip(*vals) {
                assert!((chi(pg, label, c) - v).abs() < 1e-9, "{} {label} {c}: {} vs {v}", pg.name, chi(pg, label, c));
            }
        }
    }

    #[test]
    fn td_matches_the_textbook() {
        let pg = PointGroup::named("Td").expect("Td");
        assert_eq!(pg.mats.len(), 24);
        assert_eq!(labels(&pg), ["A1", "A2", "E", "T1", "T2"]);
        check_table(
            &pg,
            &["E", "8C3", "3C2", "6S4", "6σd"],
            &[
                ("A1", &[1.0, 1.0, 1.0, 1.0, 1.0]),
                ("A2", &[1.0, 1.0, 1.0, -1.0, -1.0]),
                ("E", &[2.0, -1.0, 2.0, 0.0, 0.0]),
                ("T1", &[3.0, 0.0, -1.0, 1.0, -1.0]),
                ("T2", &[3.0, 0.0, -1.0, -1.0, 1.0]),
            ],
        );
    }

    #[test]
    fn oh_matches_the_textbook() {
        let pg = PointGroup::named("Oh").expect("Oh");
        assert_eq!(pg.mats.len(), 48);
        assert_eq!(labels(&pg), ["A1g", "A2g", "Eg", "T1g", "T2g", "A1u", "A2u", "Eu", "T1u", "T2u"]);
        let classes = ["E", "8C3", "6C4", "3C2", "6C2'", "i", "6S4", "8S6", "3σh", "6σd"];
        let rows: [(&str, [f64; 10]); 10] = [
            ("A1g", [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            ("A2g", [1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0]),
            ("Eg", [2.0, -1.0, 0.0, 2.0, 0.0, 2.0, 0.0, -1.0, 2.0, 0.0]),
            ("T1g", [3.0, 0.0, 1.0, -1.0, -1.0, 3.0, 1.0, 0.0, -1.0, -1.0]),
            ("T2g", [3.0, 0.0, -1.0, -1.0, 1.0, 3.0, -1.0, 0.0, -1.0, 1.0]),
            ("A1u", [1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0]),
            ("A2u", [1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0]),
            ("Eu", [2.0, -1.0, 0.0, 2.0, 0.0, -2.0, 0.0, 1.0, -2.0, 0.0]),
            ("T1u", [3.0, 0.0, 1.0, -1.0, -1.0, -3.0, -1.0, 0.0, 1.0, 1.0]),
            ("T2u", [3.0, 0.0, -1.0, -1.0, 1.0, -3.0, 1.0, 0.0, 1.0, -1.0]),
        ];
        let rows: Vec<(&str, &[f64])> = rows.iter().map(|(l, v)| (*l, &v[..])).collect();
        let mut sizes = pg.class_sizes();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1, 1, 3, 3, 6, 6, 6, 6, 8, 8]);
        let reordered: Vec<&str> = pg.class_labels.iter().map(String::as_str).collect();
        for (label, vals) in &rows {
            for (c, v) in classes.iter().zip(*vals) {
                assert!(reordered.contains(c), "{reordered:?}");
                assert!((chi(&pg, label, c) - v).abs() < 1e-9, "{label} {c}");
            }
        }
    }

    #[test]
    fn axial_groups_are_labelled_correctly() {
        let c2v = PointGroup::named("C2v").expect("C2v");
        assert_eq!(labels(&c2v), ["A1", "A2", "B1", "B2"]);
        assert!((chi(&c2v, "B1", "σv(xz)") - 1.0).abs() < 1e-9);
        assert!((chi(&c2v, "B2", "σv'(yz)") - 1.0).abs() < 1e-9);
        assert!((chi(&c2v, "A2", "C2") - 1.0).abs() < 1e-9 && (chi(&c2v, "A2", "σv(xz)") + 1.0).abs() < 1e-9);

        let c3v = PointGroup::named("C3v").expect("C3v");
        check_table(
            &c3v,
            &["E", "2C3", "3σv"],
            &[("A1", &[1.0, 1.0, 1.0]), ("A2", &[1.0, 1.0, -1.0]), ("E", &[2.0, -1.0, 0.0])],
        );
        let d3h = PointGroup::named("D3h").expect("D3h");
        assert_eq!(labels(&d3h), ["A1'", "A2'", "E'", "A1''", "A2''", "E''"]);
        assert!((chi(&d3h, "A2'", "3C2") + 1.0).abs() < 1e-9);
        assert!((chi(&d3h, "A1''", "σh") + 1.0).abs() < 1e-9);

        let d4h = PointGroup::named("D4h").expect("D4h");
        assert_eq!(labels(&d4h), ["A1g", "A2g", "B1g", "B2g", "Eg", "A1u", "A2u", "B1u", "B2u", "Eu"]);
        assert_eq!(pg_len(&d4h), 16);
        assert!((chi(&d4h, "B1g", "2C2'") - 1.0).abs() < 1e-9 && (chi(&d4h, "B1g", "2C2''") + 1.0).abs() < 1e-9);
        assert!((chi(&d4h, "B2g", "2C2'") + 1.0).abs() < 1e-9);

        let d2h = PointGroup::named("D2h").expect("D2h");
        assert_eq!(labels(&d2h), ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"]);
        assert!((chi(&d2h, "B1g", "C2(z)") - 1.0).abs() < 1e-9);
        let d2d = PointGroup::named("D2d").expect("D2d");
        assert_eq!(labels(&d2d), ["A1", "A2", "B1", "B2", "E"]);
        assert_eq!(d2d.class_labels, ["E", "C2", "2C2'", "2S4", "2σd"]);
        let d3d = PointGroup::named("D3d").expect("D3d");
        assert_eq!(labels(&d3d), ["A1g", "A2g", "Eg", "A1u", "A2u", "Eu"]);
        let d6h = PointGroup::named("D6h").expect("D6h");
        assert_eq!(d6h.irreps.len(), 12);
        assert!(labels(&d6h).contains(&"E2g") && labels(&d6h).contains(&"E1u"));
        let d4d = PointGroup::named("D4d").expect("D4d");
        assert_eq!(labels(&d4d), ["A1", "A2", "B1", "B2", "E1", "E2", "E3"]);
    }

    fn pg_len(pg: &PointGroup) -> usize {
        pg.mats.len()
    }

    #[test]
    fn cyclic_groups_combine_conjugate_pairs() {
        assert_eq!(labels(&PointGroup::named("C2").expect("C2")), ["A", "B"]);
        assert_eq!(labels(&PointGroup::named("C3").expect("C3")), ["A", "E"]);
        assert_eq!(labels(&PointGroup::named("C4").expect("C4")), ["A", "B", "E"]);
        assert_eq!(labels(&PointGroup::named("C6").expect("C6")), ["A", "B", "E1", "E2"]);
        assert_eq!(labels(&PointGroup::named("C5").expect("C5")), ["A", "E1", "E2"]);
        assert_eq!(labels(&PointGroup::named("S4").expect("S4")), ["A", "B", "E"]);
        assert_eq!(labels(&PointGroup::named("S6").expect("S6")), ["Ag", "Eg", "Au", "Eu"]);
        assert_eq!(labels(&PointGroup::named("Cs").expect("Cs")), ["A'", "A''"]);
        assert_eq!(labels(&PointGroup::named("Ci").expect("Ci")), ["Ag", "Au"]);
        assert_eq!(labels(&PointGroup::named("C2h").expect("C2h")), ["Ag", "Bg", "Au", "Bu"]);
        assert_eq!(labels(&PointGroup::named("C3h").expect("C3h")), ["A'", "E'", "A''", "E''"]);
        assert_eq!(labels(&PointGroup::named("T").expect("T")), ["A", "E", "T"]);
        assert_eq!(labels(&PointGroup::named("Th").expect("Th")), ["Ag", "Eg", "Tg", "Au", "Eu", "Tu"]);
        assert_eq!(labels(&PointGroup::named("O").expect("O")), ["A1", "A2", "E", "T1", "T2"]);
        assert_eq!(labels(&PointGroup::named("C1").expect("C1")), ["A"]);
        assert_eq!(labels(&PointGroup::named("D2").expect("D2")), ["A", "B1", "B2", "B3"]);
        assert_eq!(labels(&PointGroup::named("D3").expect("D3")), ["A1", "A2", "E"]);
    }

    #[test]
    fn icosahedral_groups() {
        let i = PointGroup::named("I").expect("I");
        assert_eq!(i.mats.len(), 60);
        assert_eq!(labels(&i), ["A", "T1", "T2", "G", "H"]);
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((chi(&i, "T1", "12C5") - phi).abs() < 1e-9, "{:?}", i.class_labels);
        assert!((chi(&i, "T2", "12C5") - (1.0 - phi)).abs() < 1e-9);
        assert!((chi(&i, "H", "20C3")).abs() < 1e-9 || (chi(&i, "H", "20C3") + 1.0).abs() < 1e-9);
        let ih = PointGroup::named("Ih").expect("Ih");
        assert_eq!(ih.mats.len(), 120);
        assert_eq!(ih.irreps.len(), 10);
        assert_eq!(labels(&ih), ["Ag", "T1g", "T2g", "Gg", "Hg", "Au", "T1u", "T2u", "Gu", "Hu"]);
    }

    #[test]
    fn crystallographic_groups_have_the_right_orders() {
        let orders = [
            ("C1", 1), ("Ci", 2), ("C2", 2), ("Cs", 2), ("C2h", 4), ("D2", 4), ("C2v", 4), ("D2h", 8), ("C4", 4), ("S4", 4),
            ("C4h", 8), ("D4", 8), ("C4v", 8), ("D2d", 8), ("D4h", 16), ("C3", 3), ("S6", 6), ("D3", 6), ("C3v", 6),
            ("D3d", 12), ("C6", 6), ("C3h", 6), ("C6h", 12), ("D6", 12), ("C6v", 12), ("D3h", 12), ("D6h", 24), ("T", 12),
            ("Th", 24), ("O", 24), ("Td", 24), ("Oh", 48),
        ];
        assert_eq!(orders.len(), 32);
        for (name, order) in orders {
            let pg = PointGroup::named(name).unwrap_or_else(|| panic!("{name}"));
            assert_eq!(pg.mats.len(), order, "{name}");
            // sum of squares of the dimensions (counting combined pairs once)
            let total: f64 = pg.irreps.iter().map(|i| (i.dim * i.dim) as f64 / i.norm).sum();
            assert!((total - order as f64).abs() < 1e-9, "{name}");
            assert_eq!(s(&format!("point_group_is_crystallographic({name})")), "true");
        }
        assert_eq!(s("point_group_is_crystallographic(C5)"), "false");
        assert_eq!(s("point_group_is_crystallographic(Ih)"), "false");
        assert_eq!(s("point_group_order(D5h)"), "20");
        assert_eq!(s("point_group_order(I)"), "60");
        assert_eq!(s("point_group_order(C8v)"), "16");
        assert_eq!(s("point_group_order(S8)"), "8");
    }

    #[test]
    fn hermann_mauguin_and_systems() {
        assert_eq!(s("point_group_hm(D4h)"), "4/mmm");
        assert_eq!(s("point_group_hm(C2v)"), "mm2");
        assert_eq!(s("point_group_hm(Oh)"), "m-3m");
        assert_eq!(s("point_group_hm(Td)"), "-43m");
        assert_eq!(s("point_group_from_hm(mmm)"), "D2h");
        assert_eq!(s("point_group_from_hm(hm_4_mmm)"), "D4h");
        assert_eq!(s("point_group_from_hm(hm_b43m)"), "Td");
        assert_eq!(s("point_group_from_hm(mb3m)"), "Oh");
        assert_eq!(s("point_group_from_hm(222)"), "D2");
        assert_eq!(s("point_group_from_hm(-3)"), "S6");
        assert_eq!(s("point_group_from_hm(m)"), "Cs");
        assert_eq!(s("point_group_from_hm(point_group_hm(D3d))"), "D3d");
        assert_eq!(s("point_group_crystal_system(D6h)"), "hexagonal");
        assert_eq!(s("point_group_crystal_system(Th)"), "cubic");
        // every HM name round-trips
        for (sch, hm, _) in CRYSTAL_CLASSES {
            let key = hm_key(hm);
            let back = CRYSTAL_CLASSES.iter().find(|c| hm_key(c.1) == key).map(|c| c.0);
            assert_eq!(back, Some(sch));
        }
        assert_eq!(s("crystallographic_point_groups()").matches(',').count(), 31);
        let sys = s("crystal_systems()");
        assert!(sys.contains("cubic") && sys.contains("hP"), "{sys}");
        assert_eq!(s("bravais_lattices()").matches("list(").count(), 15);
        // 7 systems and 14 lattices
        assert_eq!(super::SYSTEMS.len(), 7);
        assert_eq!(BRAVAIS.len(), 14);
        for system in super::SYSTEMS {
            assert!(CRYSTAL_CLASSES.iter().any(|c| c.2 == system));
            assert!(BRAVAIS.iter().any(|b| b.1 == system));
        }
        // the holohedry of each lattice is the largest group of its system
        for (_, system, _, holo) in BRAVAIS {
            let max = CRYSTAL_CLASSES
                .iter()
                .filter(|c| c.2 == system)
                .map(|c| PointGroup::named(c.0).map_or(0, |p| p.mats.len()))
                .max();
            assert_eq!(PointGroup::named(holo).map(|p| p.mats.len()), max, "{system}");
        }
    }

    #[test]
    fn crystallographic_restriction_theorem() {
        for n in 1..=12 {
            let allowed = matches!(n, 1 | 2 | 3 | 4 | 6);
            assert_eq!(s(&format!("crystallographic_restriction({n})")), allowed.to_string(), "n = {n}");
            assert_eq!(s(&format!("crystallographic_restriction({n}, 2)")), allowed.to_string(), "n = {n} in the plane");
        }
        assert_eq!(s("crystallographic_min_dimension(5)"), "4");
        assert_eq!(s("crystallographic_min_dimension(8)"), "4");
        assert_eq!(s("crystallographic_min_dimension(12)"), "4");
        assert_eq!(s("crystallographic_min_dimension(6)"), "2");
        assert_eq!(s("crystallographic_restriction(5, 4)"), "true");
        assert_eq!(s("crystallographic_restriction(7, 4)"), "false");
        assert_eq!(s("crystallographic_restriction(7, 6)"), "true");
    }

    #[test]
    fn spectroscopic_activity() {
        let pg = PointGroup::named("Td").expect("Td");
        let ir: Vec<&str> = pg.irreps.iter().zip(pg.ir_active()).filter(|(_, a)| *a).map(|(i, _)| i.label.as_str()).collect();
        assert_eq!(ir, ["T2"]);
        let raman: Vec<&str> =
            pg.irreps.iter().zip(pg.raman_active()).filter(|(_, a)| *a).map(|(i, _)| i.label.as_str()).collect();
        assert_eq!(raman, ["A1", "E", "T2"]);
        assert_eq!(s("point_group_ir_active(C2v)"), "list(A1, B1, B2)");
        assert_eq!(s("point_group_raman_active(C2v)"), "list(A1, A2, B1, B2)");
        // centrosymmetric groups: mutual exclusion of IR and Raman
        assert_eq!(s("point_group_ir_active(Oh)"), "list(T1u)");
        assert_eq!(s("point_group_raman_active(Oh)"), "list(A1g, Eg, T2g)");
        assert_eq!(s("point_group_ir_active(D6h)"), "list(A2u, E1u)");
        let funcs = s("point_group_function_irreps(C2v)");
        assert!(funcs.contains("list(z, list(A1))"), "{funcs}");
        assert!(funcs.contains("list(x*y, list(A2))"), "{funcs}");
        let td = s("point_group_function_irreps(Td)");
        assert!(td.contains("list(x, list(T2))"), "{td}");
        assert!(td.contains("list(x*y, list(T2))"), "{td}");
        let c2v = PointGroup::named("C2v").expect("C2v");
        // the regular representation of C2v contains each irrep once
        let regular: Vec<f64> = (0..c2v.class_elems.len()).map(|c| if c == 0 { 4.0 } else { 0.0 }).collect();
        assert_eq!(c2v.multiplicities(&regular), Some(vec![1, 1, 1, 1]));
    }

    #[test]
    fn point_group_terms() {
        assert_eq!(s("group_order(point_group(C2v))"), "4");
        assert_eq!(s("group_is_valid(point_group(Td))"), "true");
        assert_eq!(s("group_order(point_group(Ih))"), "120");
        assert_eq!(s("group_is_isomorphic(point_group(D3), symmetric_group(3))"), "true");
        assert_eq!(s("group_is_isomorphic(point_group(Td), symmetric_group(4))"), "true");
        assert_eq!(s("group_is_isomorphic(point_group(O), symmetric_group(4))"), "true");
        assert_eq!(s("group_is_isomorphic(point_group(Oh), group_direct_product(symmetric_group(4), cyclic_group(2)))"), "true");
        // the icosahedral rotation group is perfect
        assert_eq!(s("group_order(group_subgroup(point_group(I), group_derived_subgroup(point_group(I))))"), "60");
        assert_eq!(s("point_group_irreps(Td)"), "list(A1, A2, E, T1, T2)");
        let table = s("point_group_character_table(C3v)");
        assert_eq!(table, "list(list(1, 1, 1), list(1, 1, -1), list(2, -1, 0))");
        let ico = s("point_group_character_table(I)");
        assert!(ico.contains("sqrt(5)"), "{ico}");
        // the matrices are exact
        assert!(!s("group_elements(point_group(C4))").contains("1.0"));
        assert_eq!(s("point_group_class_sizes(Td)"), "list(1, 8, 3, 6, 6)");
    }

    fn water() -> &'static str {
        "list(list(O, 0, 0, 0.1173), list(H, 0, 0.7572, -0.4692), list(H, 0, -0.7572, -0.4692))"
    }

    fn methane() -> &'static str {
        "list(list(C, 0, 0, 0), list(H, 0.629, 0.629, 0.629), list(H, 0.629, -0.629, -0.629), list(H, -0.629, 0.629, -0.629), list(H, -0.629, -0.629, 0.629))"
    }

    #[test]
    fn vibrations_of_water_and_methane() {
        // H2O in the yz plane: Gamma_vib = 2 A1 + B2
        assert_eq!(s(&format!("molecule_point_group({})", water())), "C2v");
        let table = s(&format!("molecule_vibrations_table({})", water()));
        assert_eq!(table, "list(list(A1, 2), list(B2, 1))");
        let vib = s(&format!("molecule_vibrations({})", water()));
        assert!(vib.contains("2*A1") && vib.contains("B2"), "{vib}");
        let spec = s(&format!("molecule_spectroscopy({})", water()));
        assert_eq!(spec, "list(list(A1, 2, true, true), list(B2, 1, true, true))");
        let dec = s(&format!("molecule_decomposition({})", water()));
        assert!(dec.starts_with("list(list(A1, 3, 1, 0, 2)"), "{dec}");
        // the same molecule in the xz plane gives B1
        let xz = "list(list(O, 0, 0, 0.1173), list(H, 0.7572, 0, -0.4692), list(H, -0.7572, 0, -0.4692))";
        assert_eq!(s(&format!("molecule_vibrations_table({xz})")), "list(list(A1, 2), list(B1, 1))");
        // CH4: A1 + E + 2 T2
        assert_eq!(s(&format!("molecule_point_group({})", methane())), "Td");
        assert_eq!(s(&format!("molecule_vibrations_table({})", methane())), "list(list(A1, 1), list(E, 1), list(T2, 2))");
        let spec = s(&format!("molecule_spectroscopy({})", methane()));
        assert_eq!(spec, "list(list(A1, 1, false, true), list(E, 1, false, true), list(T2, 2, true, true))");
        assert_eq!(s(&format!("group_order(point_group(molecule_point_group({})))", methane())), "24");
    }

    #[test]
    fn point_groups_of_molecules() {
        // NH3 (C3v), BF3 (D3h), benzene (D6h), CO2 and HCl (linear), ethylene (D2h)
        let nh3 = "list(list(N, 0, 0, 0.1), list(H, 0.94, 0, -0.3), list(H, -0.47, 0.8139, -0.3), list(H, -0.47, -0.8139, -0.3))";
        assert_eq!(s(&format!("molecule_point_group({nh3})")), "C3v");
        let bf3 = "list(list(B, 0, 0, 0), list(F, 1.3, 0, 0), list(F, -0.65, 1.1258, 0), list(F, -0.65, -1.1258, 0))";
        assert_eq!(s(&format!("molecule_point_group({bf3})")), "D3h");
        let h = 3.0_f64.sqrt() / 2.0;
        let ring = |r: f64, e: &str| -> String {
            (0..6)
                .map(|k| {
                    let a = std::f64::consts::PI / 3.0 * f64::from(k);
                    format!("list({e}, {}, {}, 0)", r * a.cos(), r * a.sin())
                })
                .collect::<Vec<_>>()
                .join(", ")
        };
        let benzene = format!("list({}, {})", ring(1.4, "C"), ring(2.48, "H"));
        assert_eq!(s(&format!("molecule_point_group({benzene})")), "D6h");
        let _ = h;
        assert_eq!(s("molecule_point_group(list(list(O, 0, 0, 1.16), list(C, 0, 0, 0), list(O, 0, 0, -1.16)))"), "Dinfh");
        assert_eq!(s("molecule_point_group(list(list(H, 0, 0, 0), list(Cl, 0, 0, 1.27)))"), "Cinfv");
        let ethylene = "list(list(C, 0, 0, 0.667), list(C, 0, 0, -0.667), list(H, 0, 0.92, 1.24), list(H, 0, -0.92, 1.24), list(H, 0, 0.92, -1.24), list(H, 0, -0.92, -1.24))";
        assert_eq!(s(&format!("molecule_point_group({ethylene})")), "D2h");
        // a perturbed water molecule is only Cs within a small tolerance
        let skew = "list(list(O, 0, 0, 0.1173), list(H, 0, 0.7572, -0.4692), list(H, 0, -0.7572, -0.4692))";
        assert_eq!(s(&format!("molecule_point_group({skew}, 0.0001)")), "C2v");
        let bent = "list(list(O, 0, 0, 0.1173), list(H, 0, 0.7572, -0.4692), list(H, 0, -0.7, -0.4692))";
        assert_eq!(s(&format!("molecule_point_group({bent}, 0.0001)")), "Cs");
        // symmetry operations of water: E, C2, two mirrors
        let ops = s(&format!("molecule_symmetry_operations({})", water()));
        assert_eq!(ops.matches("list(list(").count(), 4, "{ops}");
    }

    #[test]
    fn decompose_reducible_representations() {
        // the permutation representation of the four H atoms of CH4 in Td
        assert_eq!(s("point_group_decompose(Td, list(4, 1, 0, 0, 2))"), s("A1 + T2"));
        assert_eq!(s("point_group_multiplicities(Td, list(4, 1, 0, 0, 2))"), "list(1, 0, 0, 0, 1)");
        // a non-character stays unreduced
        assert_eq!(s("point_group_decompose(Td, list(1, 0, 0, 0, 0))"), "point_group_decompose(Td, list(1, 0, 0, 0, 0))");
        // the regular representation of C3v, one value per element
        assert_eq!(s("point_group_multiplicities(C3v, list(1, 1, 1))"), "list(1, 0, 0)");
        assert_eq!(s("point_group_decompose(Oh, point_group_vector_character(Oh))"), "T1u");
        assert_eq!(s("point_group_decompose(Oh, point_group_rotation_character(Oh))"), "T1g");
    }
}

