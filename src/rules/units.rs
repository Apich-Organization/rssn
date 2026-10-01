//! Units of measurement and dimensional analysis.
//!
//! A physical quantity is the term `quantity(value, unit)`. The value is
//! any scalar term (exact, float or symbolic); the unit is a product of
//! unit symbols and rational powers of them: `m`, `km/h`, `kg*m/s^2`,
//! `m^2`, `1` (dimensionless). Unit symbols are case sensitive, as in the
//! SI: `m` is the metre, `s` the second, `N` the newton.
//!
//! | operator | value |
//! |---|---|
//! | `quantity(v, u)` | the quantity itself (an inert term) |
//! | `unify_expression(e)` | evaluates an expression of quantities: sums need equal dimensions and are expressed in the first operand's unit, products and powers combine units, a dimensionless result becomes a plain number; a **dimension mismatch leaves the request unreduced** |
//! | `convert(q, u)` | the quantity `q` expressed in the unit `u` (same dimension required) |
//! | `dimension(q)` | the dimension as a product of `Length`, `Mass`, `Time`, `Current`, `Temperature`, `Amount`, `Luminosity` (`1` if dimensionless) |
//! | `simplify_units(e)` | like `unify_expression`, with the result in coherent SI units: a named derived unit (`N`, `Pa`, `J`, `W`, `Hz`, `C`, `V`, `F`, `ohm`, `S`, `Wb`, `T`, `H`) when the dimension has one, otherwise a product of `m`, `kg`, `s`, `A`, `K`, `mol`, `cd` |
//!
//! # Units
//!
//! * SI base units `m`, `g` (so `kg` is `k` + `g`), `s`, `A`, `K`, `mol`,
//!   `cd`; derived `rad`, `sr`, `Hz`, `N`, `Pa`, `J`, `W`, `C`, `V`, `F`,
//!   `ohm`, `S`, `Wb`, `T`, `H`, `lm`, `lx`, `Bq`, `Gy`, `Sv`, `kat`.
//! * Every prefix from `Y` (10^24) to `y` (10^-24), including `da`, `u` or
//!   `µ` for micro, on the SI units and on `L`, `l`, `t`, `bar`, `eV`, `cal`,
//!   `Wh`, `rad`.
//! * Other units: `min`, `h`, `d`, `week`, `yr`, `deg`, `arcmin`, `arcsec`,
//!   `in`, `ft`, `yd`, `mi`, `nmi`, `angstrom`, `lb`, `oz`, `L`, `gal`,
//!   `atm`, `mmHg`, `psi`, `lbf`, `hp`, and the affine temperature units
//!   `degC` and `degF`.
//! * Long names: `meter`, `metre`, `centimeter`, `kilometer`, `kilogram`,
//!   `gram`, `second`, `minute`, `hour`, `day`, `liter`, `litre`, `kelvin`,
//!   `celsius`, `fahrenheit`, and the legacy abbreviations `sqm` (`m^2`) and
//!   `mps` (`m/s`).
//!
//! All scale factors are exact rationals (degrees carry an exact factor of
//! `pi`), so exact values convert exactly: `convert(quantity(5, km), m)` is
//! `quantity(5000, m)`. `degC` and `degF` are affine: they can only be used
//! alone, a lone `quantity(t, degC)` is converted through kelvin when it
//! takes part in arithmetic, and `convert` goes between them exactly.
//! Values with no known unit symbol (or an unknown one inside a unit) leave
//! every request unreduced.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::sync::OnceLock;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_rational::Ratio;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::arith::arith;
use crate::rules::elementary::elementary;
use crate::rules::poly::best;

/// The units rule set.
#[must_use]
pub fn units() -> RuleSet {
    RuleSet::new("units", install).needs(arith()).needs(elementary())
}

type Exp = Ratio<i32>;

/// A unit as a product of symbols with rational exponents.
type UnitMap = BTreeMap<String, Exp>;

/// Names of the seven base dimensions, in the order of [`Dim`].
const DIMENSIONS: [&str; 7] = ["Length", "Mass", "Time", "Current", "Temperature", "Amount", "Luminosity"];

/// Exponents of the seven base dimensions.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Dim([Exp; 7]);

impl Dim {
    fn zero() -> Self {
        Self([Exp::zero(); 7])
    }

    fn of(e: [i32; 7]) -> Self {
        Self(e.map(Exp::from_integer))
    }

    fn is_zero(&self) -> bool {
        self.0.iter().all(Zero::is_zero)
    }

    fn add_scaled(
        &mut self,
        other: &Self,
        k: Exp,
    ) {
        for (a, b) in self.0.iter_mut().zip(&other.0) {
            *a += *b * k;
        }
    }
}

/// An exact scale factor `rat * pi^pi`.
#[derive(Clone, Debug)]
struct Scale {
    rat: BigRational,
    pi: i32,
}

impl Scale {
    fn one() -> Self {
        Self { rat: BigRational::one(), pi: 0 }
    }

    fn is_one(&self) -> bool {
        self.rat.is_one() && self.pi == 0
    }

    fn times(
        &self,
        other: &Self,
    ) -> Self {
        Self { rat: &self.rat * &other.rat, pi: self.pi + other.pi }
    }

    fn powi(
        &self,
        k: i32,
    ) -> Self {
        Self { rat: self.rat.pow(k), pi: self.pi * k }
    }

    fn over(
        &self,
        other: &Self,
    ) -> Self {
        self.times(&other.powi(-1))
    }
}

/// What the table knows about one unit symbol.
#[derive(Clone, Debug)]
struct Resolved {
    dim: Dim,
    scale: Scale,
    /// `SI value = value * scale + offset` for affine units.
    offset: Option<BigRational>,
}

struct Def {
    dim: Dim,
    scale: Scale,
    offset: Option<BigRational>,
    prefixable: bool,
}

/// An exact decimal or fraction (`0.3048`, `1.602e-19`, `5/9`).
fn dec(text: &str) -> BigRational {
    if let Some((n, d)) = text.split_once('/') {
        return dec(n) / dec(d);
    }
    let (mantissa, exponent) = text.split_once('e').map_or((text, 0_i32), |(m, e)| (m, e.parse().unwrap_or(0)));
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits: BigInt = format!("{whole}{fraction}").parse().unwrap_or_else(|_| BigInt::zero());
    let places = i32::try_from(fraction.len()).unwrap_or(0) - exponent;
    let ten = BigInt::from(10);
    if places >= 0 {
        BigRational::new(digits, ten.pow(places.unsigned_abs()))
    } else {
        BigRational::from_integer(digits * ten.pow(places.unsigned_abs()))
    }
}

fn table() -> &'static HashMap<&'static str, Def> {
    static TABLE: OnceLock<HashMap<&'static str, Def>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = HashMap::new();
        let mut add = |name: &'static str, dim: [i32; 7], scale: &str, prefixable: bool| {
            t.insert(name, Def { dim: Dim::of(dim), scale: Scale { rat: dec(scale), pi: 0 }, offset: None, prefixable });
        };
        // (L, M, T, I, Theta, N, J)
        add("m", [1, 0, 0, 0, 0, 0, 0], "1", true);
        add("g", [0, 1, 0, 0, 0, 0, 0], "1/1000", true);
        add("s", [0, 0, 1, 0, 0, 0, 0], "1", true);
        add("A", [0, 0, 0, 1, 0, 0, 0], "1", true);
        add("K", [0, 0, 0, 0, 1, 0, 0], "1", true);
        add("mol", [0, 0, 0, 0, 0, 1, 0], "1", true);
        add("cd", [0, 0, 0, 0, 0, 0, 1], "1", true);
        add("rad", [0; 7], "1", true);
        add("sr", [0; 7], "1", true);
        add("Hz", [0, 0, -1, 0, 0, 0, 0], "1", true);
        add("N", [1, 1, -2, 0, 0, 0, 0], "1", true);
        add("Pa", [-1, 1, -2, 0, 0, 0, 0], "1", true);
        add("J", [2, 1, -2, 0, 0, 0, 0], "1", true);
        add("W", [2, 1, -3, 0, 0, 0, 0], "1", true);
        add("C", [0, 0, 1, 1, 0, 0, 0], "1", true);
        add("V", [2, 1, -3, -1, 0, 0, 0], "1", true);
        add("F", [-2, -1, 4, 2, 0, 0, 0], "1", true);
        add("ohm", [2, 1, -3, -2, 0, 0, 0], "1", true);
        add("S", [-2, -1, 3, 2, 0, 0, 0], "1", true);
        add("Wb", [2, 1, -2, -1, 0, 0, 0], "1", true);
        add("T", [0, 1, -2, -1, 0, 0, 0], "1", true);
        add("H", [2, 1, -2, -2, 0, 0, 0], "1", true);
        add("lm", [0, 0, 0, 0, 0, 0, 1], "1", true);
        add("lx", [-2, 0, 0, 0, 0, 0, 1], "1", true);
        add("Bq", [0, 0, -1, 0, 0, 0, 0], "1", true);
        add("Gy", [2, 0, -2, 0, 0, 0, 0], "1", true);
        add("Sv", [2, 0, -2, 0, 0, 0, 0], "1", true);
        add("kat", [0, 0, -1, 0, 0, 1, 0], "1", true);
        // Non-SI units accepted with the SI.
        add("min", [0, 0, 1, 0, 0, 0, 0], "60", false);
        add("h", [0, 0, 1, 0, 0, 0, 0], "3600", false);
        add("d", [0, 0, 1, 0, 0, 0, 0], "86400", false);
        add("week", [0, 0, 1, 0, 0, 0, 0], "604800", false);
        add("yr", [0, 0, 1, 0, 0, 0, 0], "31557600", false);
        add("in", [1, 0, 0, 0, 0, 0, 0], "0.0254", false);
        add("ft", [1, 0, 0, 0, 0, 0, 0], "0.3048", false);
        add("yd", [1, 0, 0, 0, 0, 0, 0], "0.9144", false);
        add("mi", [1, 0, 0, 0, 0, 0, 0], "1609.344", false);
        add("nmi", [1, 0, 0, 0, 0, 0, 0], "1852", false);
        add("angstrom", [1, 0, 0, 0, 0, 0, 0], "1e-10", false);
        add("lb", [0, 1, 0, 0, 0, 0, 0], "0.45359237", false);
        add("oz", [0, 1, 0, 0, 0, 0, 0], "0.028349523125", false);
        add("t", [0, 1, 0, 0, 0, 0, 0], "1000", true);
        add("L", [3, 0, 0, 0, 0, 0, 0], "1/1000", true);
        add("l", [3, 0, 0, 0, 0, 0, 0], "1/1000", true);
        add("gal", [3, 0, 0, 0, 0, 0, 0], "0.003785411784", false);
        add("bar", [-1, 1, -2, 0, 0, 0, 0], "100000", true);
        add("atm", [-1, 1, -2, 0, 0, 0, 0], "101325", false);
        add("mmHg", [-1, 1, -2, 0, 0, 0, 0], "133.322387415", false);
        add("psi", [-1, 1, -2, 0, 0, 0, 0], "6894.757293168", false);
        add("eV", [2, 1, -2, 0, 0, 0, 0], "1.602176634e-19", true);
        add("cal", [2, 1, -2, 0, 0, 0, 0], "4.184", true);
        add("Wh", [2, 1, -2, 0, 0, 0, 0], "3600", true);
        add("lbf", [1, 1, -2, 0, 0, 0, 0], "4.4482216152605", false);
        // hp = 550 ft lbf / s
        add("hp", [2, 1, -3, 0, 0, 0, 0], "745.69987158227022", false);
        let angle = |name: &'static str, t: &mut HashMap<&'static str, Def>, denominator: &str| {
            t.insert(
                name,
                Def { dim: Dim::zero(), scale: Scale { rat: dec(&format!("1/{denominator}")), pi: 1 }, offset: None, prefixable: false },
            );
        };
        angle("deg", &mut t, "180");
        angle("arcmin", &mut t, "10800");
        angle("arcsec", &mut t, "648000");
        t.insert(
            "degC",
            Def { dim: Dim::of([0, 0, 0, 0, 1, 0, 0]), scale: Scale::one(), offset: Some(dec("273.15")), prefixable: false },
        );
        t.insert(
            "degF",
            Def {
                dim: Dim::of([0, 0, 0, 0, 1, 0, 0]),
                scale: Scale { rat: dec("5/9"), pi: 0 },
                offset: Some(dec("45967/180")),
                prefixable: false,
            },
        );
        t
    })
}

/// SI prefixes and their powers of ten.
const PREFIXES: [(&str, i32); 22] = [
    ("Y", 24),
    ("Z", 21),
    ("E", 18),
    ("P", 15),
    ("T", 12),
    ("G", 9),
    ("M", 6),
    ("k", 3),
    ("h", 2),
    ("da", 1),
    ("d", -1),
    ("c", -2),
    ("m", -3),
    ("u", -6),
    ("\u{b5}", -6),
    ("\u{3bc}", -6),
    ("n", -9),
    ("p", -12),
    ("f", -15),
    ("a", -18),
    ("z", -21),
    ("y", -24),
];

/// Names that stand for a product of other units.
const ALIASES: [(&str, &[(&str, i32)]); 19] = [
    ("sqm", &[("m", 2)]),
    ("mps", &[("m", 1), ("s", -1)]),
    ("meter", &[("m", 1)]),
    ("metre", &[("m", 1)]),
    ("centimeter", &[("cm", 1)]),
    ("kilometer", &[("km", 1)]),
    ("kilogram", &[("kg", 1)]),
    ("gram", &[("g", 1)]),
    ("second", &[("s", 1)]),
    ("minute", &[("min", 1)]),
    ("hour", &[("h", 1)]),
    ("day", &[("d", 1)]),
    ("liter", &[("L", 1)]),
    ("litre", &[("L", 1)]),
    ("kelvin", &[("K", 1)]),
    ("celsius", &[("degC", 1)]),
    ("fahrenheit", &[("degF", 1)]),
    ("inch", &[("in", 1)]),
    ("foot", &[("ft", 1)]),
];

fn resolve_symbol(name: &str) -> Option<Resolved> {
    if let Some((_, parts)) = ALIASES.iter().find(|(alias, _)| *alias == name) {
        let map: UnitMap = parts.iter().map(|&(n, e)| (n.to_owned(), Exp::from_integer(e))).collect();
        return resolve_map(&map);
    }
    let table = table();
    if let Some(d) = table.get(name) {
        return Some(Resolved { dim: d.dim.clone(), scale: d.scale.clone(), offset: d.offset.clone() });
    }
    for (prefix, power) in PREFIXES {
        let Some(rest) = name.strip_prefix(prefix) else {
            continue;
        };
        if let Some(d) = table.get(rest).filter(|d| d.prefixable) {
            let factor = Scale { rat: BigRational::from_integer(BigInt::from(10)).pow(power), pi: 0 };
            return Some(Resolved { dim: d.dim.clone(), scale: d.scale.times(&factor), offset: None });
        }
    }
    None
}

/// Dimension, scale and affine offset of a unit product.
fn resolve_map(map: &UnitMap) -> Option<Resolved> {
    let mut dim = Dim::zero();
    let mut scale = Scale::one();
    let mut offset = None;
    for (name, &e) in map {
        let r = resolve_symbol(name)?;
        dim.add_scaled(&r.dim, e);
        if e.is_integer() {
            scale = scale.times(&r.scale.powi(*e.numer()));
        } else if !r.scale.is_one() {
            return None;
        }
        if r.offset.is_some() {
            if map.len() != 1 || !e.is_one() {
                return None;
            }
            offset = r.offset;
        }
    }
    Some(Resolved { dim, scale, offset })
}

/// The unit product a term denotes.
fn parse_unit(
    cx: &Cx<'_>,
    node: NodeId,
    depth: usize,
) -> Option<UnitMap> {
    if depth > 16 {
        return None;
    }
    let graph = &*cx.graph;
    if let Some(n) = graph.number_of(node) {
        return n.is_one().then(UnitMap::new);
    }
    if let Some(symbol) = graph.symbol_of(node) {
        let mut map = UnitMap::new();
        map.insert(graph.interner().symbol_name(symbol).to_owned(), Exp::one());
        return Some(map);
    }
    // Several forms may share a class: take the first that reads as a unit.
    graph.enodes(graph.find(node)).find_map(|e| {
        let children = graph.children(e).to_vec();
        match graph.op(e) {
            | core::MUL => {
                let mut map = UnitMap::new();
                for c in children {
                    for (name, k) in parse_unit(cx, c, depth + 1)? {
                        *map.entry(name).or_insert_with(Exp::zero) += k;
                    }
                }
                Some(clean(map))
            },
            | core::POW => {
                let &[base, exponent] = children.as_slice() else {
                    return None;
                };
                let r = graph.number_of(exponent)?.to_rational()?;
                let e = Exp::new(r.numer().to_i32()?, r.denom().to_i32()?);
                let map = parse_unit(cx, base, depth + 1)?;
                Some(clean(map.into_iter().map(|(n, k)| (n, k * e)).collect()))
            },
            | _ => None,
        }
    })
}

fn clean(mut map: UnitMap) -> UnitMap {
    map.retain(|_, e| !e.is_zero());
    map
}

/// A scalar or a quantity with a unit.
#[derive(Clone, Debug)]
enum Val {
    Scalar(NodeId),
    Quant { value: NodeId, unit: UnitMap },
}

struct Units<'c, 'a> {
    cx: &'c mut Cx<'a>,
    quantity: OpId,
}

impl Units<'_, '_> {
    fn scale_term(
        &mut self,
        scale: &Scale,
    ) -> NodeId {
        let rat = self.cx.graph.num(Number::rat(scale.rat.clone()));
        if scale.pi == 0 {
            return rat;
        }
        let pi = self.cx.graph.ops().lookup("pi").map(|p| self.cx.graph.node(p, &[]));
        match pi {
            | Some(pi) => {
                let e = self.cx.graph.int(i64::from(scale.pi));
                let power = self.cx.graph.node(core::POW, &[pi, e]);
                self.cx.graph.node(core::MUL, &[rat, power])
            },
            | None => rat,
        }
    }

    /// `value * scale`, simplified.
    fn scaled(
        &mut self,
        value: NodeId,
        scale: &Scale,
    ) -> NodeId {
        if scale.is_one() {
            return value;
        }
        let factor = self.scale_term(scale);
        let product = self.cx.graph.node(core::MUL, &[value, factor]);
        self.cx.simplify(product)
    }

    fn unit_term(
        &mut self,
        unit: &UnitMap,
    ) -> NodeId {
        if unit.is_empty() {
            return self.cx.graph.int(1);
        }
        let mut factors = Vec::new();
        for (name, &e) in unit {
            let symbol = self.cx.graph.sym(name);
            factors.push(if e.is_one() {
                symbol
            } else {
                let exponent = self.cx.graph.num(Number::fraction(i64::from(*e.numer()), i64::from(*e.denom())).unwrap_or_else(|| Number::from(1)));
                self.cx.graph.node(core::POW, &[symbol, exponent])
            });
        }
        match factors.as_slice() {
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, &factors),
        }
    }

    fn quantity_term(
        &mut self,
        value: NodeId,
        unit: &UnitMap,
    ) -> NodeId {
        let unit = self.unit_term(unit);
        let quantity = self.quantity;
        self.cx.graph.node(quantity, &[value, unit])
    }

    fn rational(
        &self,
        node: NodeId,
    ) -> Option<BigRational> {
        self.cx.graph.number_of(node).and_then(Number::to_rational)
    }

    /// Reads `quantity(v, u)` without evaluating anything else.
    fn lone_quantity(
        &mut self,
        node: NodeId,
    ) -> Option<Val> {
        let term = best(self.cx.graph, node)?;
        if self.cx.graph.op(term) != self.quantity {
            return None;
        }
        let &[value, unit] = self.cx.graph.children(term) else {
            return None;
        };
        let unit = parse_unit(self.cx, unit, 0)?;
        let resolved = resolve_map(&unit)?;
        let Val::Scalar(value) = self.eval(value, 0)? else {
            return None;
        };
        // Affine units are converted through kelvin.
        if let Some(offset) = resolved.offset {
            let scaled = self.scaled(value, &resolved.scale);
            let shift = self.cx.graph.num(Number::rat(offset));
            let sum = self.cx.graph.node(core::ADD, &[scaled, shift]);
            let sum = self.cx.simplify(sum);
            let mut kelvin = UnitMap::new();
            kelvin.insert("K".to_owned(), Exp::one());
            return Some(Val::Quant { value: sum, unit: kelvin });
        }
        Some(Val::Quant { value, unit })
    }

    /// A dimensionless quantity as the number it stands for.
    fn demote(
        &mut self,
        val: Val,
    ) -> Option<Val> {
        match val {
            | Val::Quant { value, unit } => self.collapse(value, unit),
            | scalar @ Val::Scalar(_) => Some(scalar),
        }
    }

    /// A product of unit symbols as a scalar when it is dimensionless.
    fn collapse(
        &mut self,
        value: NodeId,
        unit: UnitMap,
    ) -> Option<Val> {
        let resolved = resolve_map(&unit)?;
        if unit.is_empty() || resolved.dim.is_zero() {
            return Some(Val::Scalar(self.scaled(value, &resolved.scale)));
        }
        Some(Val::Quant { value, unit })
    }

    #[allow(clippy::too_many_lines)]
    fn eval(
        &mut self,
        node: NodeId,
        depth: usize,
    ) -> Option<Val> {
        if depth > 48 {
            return None;
        }
        if let Some(q) = self.lone_quantity(node) {
            return Some(q);
        }
        let term = best(self.cx.graph, node)?;
        let op = self.cx.graph.op(term);
        let children = self.cx.graph.children(term).to_vec();
        if children.is_empty() {
            return Some(Val::Scalar(term));
        }
        if op == self.quantity {
            // A malformed quantity or one with an unknown unit.
            return None;
        }
        match op {
            | core::ADD => {
                let vals: Vec<Val> = children.iter().map(|&c| self.eval(c, depth + 1).and_then(|v| self.demote(v))).collect::<Option<_>>()?;
                if vals.iter().all(|v| matches!(v, Val::Scalar(_))) {
                    let nodes: Vec<NodeId> = vals.iter().map(|v| if let Val::Scalar(n) = v { *n } else { term }).collect();
                    return Some(Val::Scalar(self.cx.graph.node(core::ADD, &nodes)));
                }
                let Some(Val::Quant { unit: first_unit, .. }) = vals.iter().find(|v| matches!(v, Val::Quant { .. })).cloned() else {
                    return None;
                };
                let first = resolve_map(&first_unit)?;
                let mut terms = Vec::new();
                for v in vals {
                    let Val::Quant { value, unit } = v else {
                        return None; // a scalar added to a quantity
                    };
                    let r = resolve_map(&unit)?;
                    if r.dim != first.dim {
                        return None;
                    }
                    terms.push(self.scaled(value, &r.scale.over(&first.scale)));
                }
                let sum = self.cx.graph.node(core::ADD, &terms);
                Some(Val::Quant { value: self.cx.simplify(sum), unit: first_unit })
            },
            | core::MUL => {
                let mut values = Vec::new();
                let mut unit = UnitMap::new();
                let mut any = false;
                for c in children {
                    match self.eval(c, depth + 1)? {
                        | Val::Scalar(n) => values.push(n),
                        | Val::Quant { value, unit: u } => {
                            any = true;
                            values.push(value);
                            for (name, k) in u {
                                *unit.entry(name).or_insert_with(Exp::zero) += k;
                            }
                        },
                    }
                }
                let value = self.cx.graph.node(core::MUL, &values);
                if !any {
                    return Some(Val::Scalar(value));
                }
                self.collapse(value, clean(unit))
            },
            | core::POW => {
                let &[base, exponent] = children.as_slice() else {
                    return None;
                };
                let (b, e) = (self.eval(base, depth + 1)?, self.eval(exponent, depth + 1)?);
                let e = self.demote(e)?;
                match (b, e) {
                    | (Val::Scalar(b), Val::Scalar(e)) => Some(Val::Scalar(self.cx.graph.node(core::POW, &[b, e]))),
                    | (Val::Quant { value, unit }, Val::Scalar(e)) => {
                        let r = self.rational(e)?;
                        let k = Exp::new(r.numer().to_i32()?, r.denom().to_i32()?);
                        let value = self.cx.graph.node(core::POW, &[value, e]);
                        self.collapse(value, clean(unit.into_iter().map(|(n, x)| (n, x * k)).collect()))
                    },
                    | _ => None,
                }
            },
            | _ => {
                // A function of scalars.
                let mut args = Vec::with_capacity(children.len());
                for c in children {
                    let value = self.eval(c, depth + 1)?;
                    let Val::Scalar(n) = self.demote(value)? else {
                        return None;
                    };
                    args.push(n);
                }
                Some(Val::Scalar(self.cx.graph.try_node(op, &args)?))
            },
        }
    }

    fn finish(
        &mut self,
        val: &Val,
    ) -> NodeId {
        match val {
            | Val::Scalar(n) => self.cx.simplify(*n),
            | Val::Quant { value, unit } => {
                let value = self.cx.simplify(*value);
                self.quantity_term(value, unit)
            },
        }
    }

    /// The value in coherent SI units.
    fn to_si(
        &mut self,
        val: &Val,
    ) -> Option<NodeId> {
        let Val::Quant { value, unit } = val else {
            return Some(self.finish(val));
        };
        let r = resolve_map(unit)?;
        let value = self.scaled(*value, &r.scale);
        if r.dim.is_zero() {
            return Some(value);
        }
        let coherent = coherent_unit(&r.dim);
        Some(self.quantity_term(value, &coherent))
    }

    fn convert(
        &mut self,
        q: NodeId,
        target: NodeId,
    ) -> Option<NodeId> {
        let val = self.eval(q, 0)?;
        let Val::Quant { value, unit } = val else {
            return None;
        };
        let target_unit = parse_unit(self.cx, target, 0)?;
        let source = resolve_map(&unit)?;
        let destination = resolve_map(&target_unit)?;
        if source.dim != destination.dim {
            return None;
        }
        let si = self.scaled(value, &source.scale);
        let shifted = match &destination.offset {
            | Some(offset) => {
                let shift = self.cx.graph.num(Number::rat(-offset.clone()));
                let sum = self.cx.graph.node(core::ADD, &[si, shift]);
                self.cx.simplify(sum)
            },
            | None => si,
        };
        let out = self.scaled(shifted, &Scale::one().over(&destination.scale));
        Some(self.quantity_term(out, &target_unit))
    }

    fn dimension(
        &mut self,
        q: NodeId,
    ) -> Option<NodeId> {
        let val = self.eval(q, 0)?;
        let dim = match val {
            | Val::Scalar(_) => Dim::zero(),
            | Val::Quant { unit, .. } => resolve_map(&unit)?.dim,
        };
        let mut factors = Vec::new();
        for (name, &e) in DIMENSIONS.iter().zip(&dim.0) {
            if e.is_zero() {
                continue;
            }
            let symbol = self.cx.graph.sym(name);
            factors.push(if e.is_one() {
                symbol
            } else {
                let exponent = self.cx.graph.num(Number::fraction(i64::from(*e.numer()), i64::from(*e.denom())).unwrap_or_else(|| Number::from(1)));
                self.cx.graph.node(core::POW, &[symbol, exponent])
            });
        }
        let product = match factors.as_slice() {
            | [] => self.cx.graph.int(1),
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, &factors),
        };
        Some(product)
    }
}

/// The coherent SI unit of a dimension: a named derived unit when there
/// is one, else the product of base units.
fn coherent_unit(dim: &Dim) -> UnitMap {
    const NAMED: [(&str, [i32; 7]); 13] = [
        ("N", [1, 1, -2, 0, 0, 0, 0]),
        ("Pa", [-1, 1, -2, 0, 0, 0, 0]),
        ("J", [2, 1, -2, 0, 0, 0, 0]),
        ("W", [2, 1, -3, 0, 0, 0, 0]),
        ("Hz", [0, 0, -1, 0, 0, 0, 0]),
        ("C", [0, 0, 1, 1, 0, 0, 0]),
        ("V", [2, 1, -3, -1, 0, 0, 0]),
        ("F", [-2, -1, 4, 2, 0, 0, 0]),
        ("ohm", [2, 1, -3, -2, 0, 0, 0]),
        ("S", [-2, -1, 3, 2, 0, 0, 0]),
        ("Wb", [2, 1, -2, -1, 0, 0, 0]),
        ("T", [0, 1, -2, -1, 0, 0, 0]),
        ("H", [2, 1, -2, -2, 0, 0, 0]),
    ];
    let mut map = UnitMap::new();
    if let Some((name, _)) = NAMED.iter().find(|(_, d)| Dim::of(*d) == *dim) {
        map.insert((*name).to_owned(), Exp::one());
        return map;
    }
    for (name, &e) in ["m", "kg", "s", "A", "K", "mol", "cd"].iter().zip(&dim.0) {
        if !e.is_zero() {
            map.insert((*name).to_owned(), e);
        }
    }
    map
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Unify,
    Convert,
    Dimension,
    Simplify,
}

struct Request1 {
    op: OpId,
    request: Request,
    quantity: OpId,
}

impl Kernel for Request1 {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let mut units = Units { cx, quantity: self.quantity };
        let result = match (self.request, args.as_slice()) {
            | (Request::Unify, &[e]) => units.eval(e, 0).map(|v| units.finish(&v)),
            | (Request::Simplify, &[e]) => units.eval(e, 0).and_then(|v| units.to_si(&v)),
            | (Request::Convert, &[q, u]) => units.convert(q, u),
            | (Request::Dimension, &[q]) => units.dimension(q),
            | _ => None,
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let quantity = i.op(OpDescriptor::new("quantity", Arity::Fixed(2)))?;
    for (name, arity, request) in [
        ("unify_expression", 1, Request::Unify),
        ("convert", 2, Request::Convert),
        ("dimension", 1, Request::Dimension),
        ("simplify_units", 1, Request::Simplify),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity)).flags(OpFlags::HEAVY).cost(100))?;
        i.kernel(&format!("units/{name}"), Tier::Reduce, Request1 { op, request, quantity });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[units()], src)
    }

    /// The text and whether the request was reduced.
    fn attempt(src: &str) -> (String, bool) {
        reduce_with(&[units()], src, &[])
    }

    #[test]
    fn quantities_are_inert_terms() {
        assert_eq!(run("quantity(5, m)"), "quantity(5, m)");
        assert_eq!(run("quantity(x, kg*m/s^2)"), "quantity(x, kg*m/s^2)");
    }

    #[test]
    fn sums_convert_to_the_first_unit() {
        assert_eq!(run("unify_expression(quantity(5, m) + quantity(3, m))"), "quantity(8, m)");
        assert_eq!(run("unify_expression(quantity(1, km) + quantity(500, m))"), "quantity(3/2, km)");
        assert_eq!(run("unify_expression(quantity(500, m) + quantity(1, km))"), "quantity(1500, m)");
        assert_eq!(run("unify_expression(quantity(2, h) - quantity(30, min))"), "quantity(3/2, h)");
        assert_eq!(run("unify_expression(quantity(1, ft) + quantity(12, in))"), "quantity(2, ft)");
    }

    #[test]
    fn mismatched_dimensions_do_not_reduce() {
        for src in [
            "unify_expression(quantity(5, m) + quantity(3, s))",
            "unify_expression(quantity(5, m) + 3)",
            "unify_expression(quantity(5, m) + x)",
            "unify_expression(quantity(1, kg) - quantity(1, N))",
            "convert(quantity(5, m), s)",
            "unify_expression(sin(quantity(5, m)))",
            "unify_expression(quantity(5, parsec))",
        ] {
            let (text, reduced) = attempt(src);
            assert!(!reduced, "{src} reduced to {text}");
        }
    }

    #[test]
    fn products_quotients_and_powers() {
        assert_eq!(run("unify_expression(quantity(10, m) / quantity(2, s))"), "quantity(5, m/s)");
        assert_eq!(run("unify_expression(quantity(3, m) * quantity(4, m))"), "quantity(12, m^2)");
        assert_eq!(run("unify_expression(quantity(3, m)^2)"), "quantity(9, m^2)");
        assert_eq!(run("unify_expression(quantity(9, m^2)^(1/2))"), "quantity(3, m)");
        assert_eq!(run("unify_expression(2 * quantity(3, kg))"), "quantity(6, kg)");
        assert_eq!(run("unify_expression(quantity(6, kg) / 3)"), "quantity(2, kg)");
        assert_eq!(run("unify_expression(1 / quantity(4, s))"), "quantity(1/4, 1/s)");
        assert_eq!(run("unify_expression(x * quantity(2, m))"), "quantity(2*x, m)");
        // force = mass * acceleration
        assert_eq!(run("unify_expression(quantity(2, kg) * quantity(3, m/s^2))"), "quantity(6, kg*m/s^2)");
    }

    #[test]
    fn dimensionless_results_become_numbers() {
        assert_eq!(run("unify_expression(quantity(1, m) / quantity(1, cm))"), "100");
        assert_eq!(run("unify_expression(quantity(6, m) / quantity(2, m))"), "3");
        assert_eq!(run("unify_expression(quantity(3, s) * quantity(2, Hz))"), "6");
        assert_eq!(run("unify_expression(sin(quantity(3, m) / quantity(3, m)))"), "sin(1)");
        assert_eq!(run("unify_expression(quantity(5, m) + quantity(3, m) - quantity(8, m))"), "quantity(0, m)");
    }

    #[test]
    fn conversion() {
        assert_eq!(run("convert(quantity(5, km), m)"), "quantity(5000, m)");
        assert_eq!(run("convert(quantity(250, cm), m)"), "quantity(5/2, m)");
        assert_eq!(run("convert(quantity(1, h), s)"), "quantity(3600, s)");
        assert_eq!(run("convert(quantity(36, km/h), m/s)"), "quantity(10, m/s)");
        assert_eq!(run("convert(quantity(1, mi), ft)"), "quantity(5280, ft)");
        assert_eq!(run("convert(quantity(1, lb), g)"), "quantity(45359237/100000, g)");
        assert_eq!(run("convert(quantity(1, N), kg*m/s^2)"), "quantity(1, kg*m/s^2)");
        assert_eq!(run("convert(quantity(1, kWh), J)"), "quantity(3600000, J)");
        assert_eq!(run("convert(quantity(1, L), cm^3)"), "quantity(1000, cm^3)");
        assert_eq!(run("convert(quantity(180, deg), rad)"), "quantity(pi, rad)");
        assert_eq!(run("convert(quantity(x, mm), m)"), "quantity(1/1000*x, m)");
        assert_eq!(run("convert(quantity(2, GHz), Hz)"), "quantity(2000000000, Hz)");
        assert_eq!(run("convert(quantity(1, ug), kg)"), "quantity(1/1000000000, kg)");
        assert_eq!(run("convert(quantity(1, kilogram), gram)"), "quantity(1000, gram)");
        // expressions are evaluated first
        assert_eq!(run("convert(quantity(1, km) + quantity(500, m), m)"), "quantity(1500, m)");
    }

    #[test]
    fn temperatures_are_affine() {
        assert_eq!(run("convert(quantity(0, degC), K)"), "quantity(5463/20, K)");
        assert_eq!(run("convert(quantity(100, degC), degF)"), "quantity(212, degF)");
        assert_eq!(run("convert(quantity(32, degF), degC)"), "quantity(0, degC)");
        assert_eq!(run("convert(quantity(20, degC), degC)"), "quantity(20, degC)");
        assert_eq!(run("convert(quantity(300, K), degC)"), "quantity(537/20, degC)");
        // they cannot be combined multiplicatively
        assert!(!attempt("unify_expression(quantity(1, degC*m))").1);
    }

    #[test]
    fn dimensions() {
        assert_eq!(run("dimension(quantity(5, m))"), "Length");
        assert_eq!(run("dimension(quantity(5, m/s))"), "Length/Time");
        assert_eq!(run("dimension(quantity(5, N))"), "Length*Mass/Time^2");
        assert_eq!(run("dimension(quantity(5, J))"), "Length^2*Mass/Time^2");
        assert_eq!(run("dimension(quantity(5, km^2))"), "Length^2");
        assert_eq!(run("dimension(quantity(5, mol/L))"), "Amount/Length^3");
        assert_eq!(run("dimension(quantity(5, rad))"), "1");
        assert_eq!(run("dimension(7)"), "1");
        assert_eq!(run("dimension(quantity(1, kg) * quantity(1, m) / quantity(1, s)^2)"), "Length*Mass/Time^2");
        // equal dimensions: N m and J
        assert_eq!(run("dimension(quantity(1, N*m))"), run("dimension(quantity(1, J))"));
    }

    #[test]
    fn simplification_to_si() {
        assert_eq!(run("simplify_units(quantity(36, km/h))"), "quantity(10, m/s)");
        assert_eq!(run("simplify_units(quantity(2, kg) * quantity(3, m/s^2))"), "quantity(6, N)");
        assert_eq!(run("simplify_units(quantity(1, N) * quantity(2, m))"), "quantity(2, J)");
        assert_eq!(run("simplify_units(quantity(1, J) / quantity(2, s))"), "quantity(1/2, W)");
        assert_eq!(run("simplify_units(quantity(1, V) / quantity(2, A))"), "quantity(1/2, ohm)");
        assert_eq!(run("simplify_units(quantity(1, cm) + quantity(1, m))"), "quantity(101/100, m)");
        assert_eq!(run("simplify_units(quantity(90, deg) + quantity(1, m) / quantity(1, m))"), "1/2*pi + 1");
        assert_eq!(run("simplify_units(quantity(1, kW) * quantity(1, h))"), "quantity(3600000, J)");
        assert_eq!(run("simplify_units(quantity(4, s) * quantity(5, Hz))"), "20");
        assert_eq!(run("simplify_units(quantity(5, m))"), "quantity(5, m)");
        assert_eq!(run("simplify_units(quantity(1, ft)^3)"), "quantity(55306341/1953125000, m^3)");
    }

    #[test]
    fn legacy_units_are_accepted() {
        assert_eq!(run("unify_expression(quantity(5, m) + quantity(300, cm))"), "quantity(8, m)");
        assert_eq!(run("unify_expression(quantity(1, kg) + quantity(500, g))"), "quantity(3/2, kg)");
        assert_eq!(run("unify_expression(quantity(60, s) + quantity(1, min))"), "quantity(120, s)");
        assert_eq!(run("convert(quantity(2, sqm), cm^2)"), "quantity(20000, cm^2)");
        assert_eq!(run("convert(quantity(3, mps), km/h)"), "quantity(54/5, km/h)");
        assert_eq!(run("unify_expression(quantity(10, meter) / quantity(2, second))"), "quantity(5, meter/second)");
        assert_eq!(run("unify_expression(quantity(5, m) * quantity(5, m))"), "quantity(25, m^2)");
    }

    #[test]
    fn every_unit_resolves_and_is_consistent() {
        // exact table entries
        for name in table().keys() {
            assert!(resolve_symbol(name).is_some(), "{name}");
        }
        // prefixes
        for (prefix, power) in PREFIXES {
            let r = resolve_symbol(&format!("{prefix}m")).unwrap_or_else(|| panic!("{prefix}m"));
            assert_eq!(r.scale.rat, BigRational::from_integer(BigInt::from(10)).pow(power), "{prefix}");
            assert_eq!(r.dim, Dim::of([1, 0, 0, 0, 0, 0, 0]));
        }
        // kg is the base unit of mass
        let kg = resolve_symbol("kg").unwrap_or_else(|| panic!("kg"));
        assert!(kg.scale.is_one());
        // derived units are consistent with their definitions
        let check = |name: &str, product: &[(&str, i32)]| {
            let map: UnitMap = product.iter().map(|&(n, e)| (n.to_owned(), Exp::from_integer(e))).collect();
            let a = resolve_symbol(name).unwrap_or_else(|| panic!("{name}"));
            let b = resolve_map(&map).unwrap_or_else(|| panic!("{name}"));
            assert_eq!(a.dim, b.dim, "{name}");
            assert_eq!(a.scale.rat, b.scale.rat, "{name}");
        };
        check("N", &[("kg", 1), ("m", 1), ("s", -2)]);
        check("J", &[("N", 1), ("m", 1)]);
        check("W", &[("J", 1), ("s", -1)]);
        check("Pa", &[("N", 1), ("m", -2)]);
        check("V", &[("W", 1), ("A", -1)]);
        check("ohm", &[("V", 1), ("A", -1)]);
        check("F", &[("C", 1), ("V", -1)]);
        check("Wb", &[("V", 1), ("s", 1)]);
        check("T", &[("Wb", 1), ("m", -2)]);
        check("H", &[("Wb", 1), ("A", -1)]);
        check("S", &[("ohm", -1)]);
        // hp = 550 ft lbf / s and lbf = lb g0 with g0 = 9.80665 m/s^2
        let hp = resolve_symbol("hp").unwrap_or_else(|| panic!("hp"));
        let ft_lbf = resolve_symbol("ft").zip(resolve_symbol("lbf")).map(|(a, b)| a.scale.rat * b.scale.rat);
        assert_eq!(Some(hp.scale.rat), ft_lbf.map(|r| r * BigRational::from_integer(BigInt::from(550))));
        let lbf = resolve_symbol("lbf").unwrap_or_else(|| panic!("lbf"));
        assert_eq!(lbf.scale.rat, dec("0.45359237") * dec("9.80665"));
        // non-unit symbols are rejected
        assert!(resolve_symbol("parsec").is_none());
        assert!(resolve_symbol("xm").is_none());
    }

    #[test]
    fn unit_names_do_not_collide_with_operators() {
        let mut graph = crate::graph::Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard()).unwrap_or_else(|e| panic!("{e}"));
        let mut clashes: Vec<&str> = table().keys().copied().filter(|name| graph.ops().lookup(name).is_some()).collect();
        clashes.extend(ALIASES.iter().map(|(a, _)| *a).filter(|name| graph.ops().lookup(name).is_some()));
        clashes.sort_unstable();
        assert!(clashes.is_empty(), "unit names that are also operators: {clashes:?}");
    }

    #[test]
    fn numeric_values_and_symbolic_mix() {
        assert_eq!(run("unify_expression(quantity(2.5, m) + quantity(50, cm))"), "quantity(3, m)");
        assert_eq!(run("unify_expression(quantity(a, m) + quantity(b, m))"), "quantity(a + b, m)");
        assert_eq!(run("unify_expression(quantity(a, km) + quantity(b, m))"), "quantity(a + 1/1000*b, km)");
        // sums of products
        assert_eq!(
            run("unify_expression(quantity(1, m) * quantity(2, m) + quantity(3, m^2))"),
            "quantity(5, m^2)"
        );
        assert_eq!(
            run("unify_expression(quantity(1, kg) * quantity(10, m/s^2) * quantity(2, m))"),
            "quantity(20, kg*m^2/s^2)"
        );
    }
}
