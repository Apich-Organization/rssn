//! Units: the canonical registry, with exact scales.
//!
//! A unit is a name, a dimension, an exact scale to coherent SI units and,
//! for a few units, an exact additive offset. Everything here is exact: the
//! scale of a foot is the rational `381/1250`, the scale of a degree is
//! `pi/180`, the offset of `degC` is `27315/100`. A float would lose the
//! distinction between `1 in` and `0.0254 m` after a few conversions; a
//! rational does not.
//!
//! # One table
//!
//! This module owns the only unit table in rssn. The `units` rule set
//! (which evaluates `quantity`, `convert`, `dimension`, `unify_expression`)
//! reads it, the dimension layer reads it to turn unit text into a
//! [`DimTerm`](crate::graph::dim::DimTerm), and the
//! [manifest](crate::kernel::manifest) lists it so that a host can generate
//! its own `Units.SI` instead of maintaining a second table.
//!
//! # Extending it
//!
//! [`register`] adds a unit to a process-wide overlay, which
//! [`UnitRegistry::register`] does through a [`Registry`] when a plugin
//! declares one. Registration is additive and checked: a name or alias that
//! already belongs to a different unit is refused, so a plugin cannot
//! silently redefine the metre.
//!
//! # Limitations
//!
//! The table holds ASCII spellings only; the Unicode forms of units (`Å` for
//! `angstrom`, `Ω` for `ohm`, `°C` for `degC`) are the host's aliases, because
//! the parser never sees them (a host generates its `Units` module from the
//! manifest and maps them there). The information units (`bit`, `B`, `KiB`)
//! are dimensionless by design — `spec/05-units.md` §5.3.3 keeps them so that
//! entropy and probability formulas that use bits stay dimensionally
//! transparent; a host that wants a base dimension declares `Information` and
//! registers its own units on it. A compound unit whose scale has a
//! non-integer exponent of a non-coherent unit (`ft^(1/2)`) does not resolve,
//! which is deliberate: only coherent units compose under a root.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::OnceLock;
use std::sync::RwLock;

use serde::Deserialize;
use serde::Serialize;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_rational::Ratio;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::dim::DimArena;
use crate::graph::dim::DimError;
use crate::graph::dim::DimTerm;
use crate::graph::op::core;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp as TypeExp;
use crate::kernel::signature::Rational;

/// A rational exponent in a unit product.
pub(crate) type Exp = Ratio<i32>;

/// A named unit or an immutable definition carried across the value ABI.
/// Ordering is structural: two definitions with the same name but different
/// dimensions/scales must remain distinct factors.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum UnitAtom {
    Named(String),
    Literal(Arc<crate::kernel::value::UnitValue>),
}

impl From<String> for UnitAtom {
    fn from(name: String) -> Self {
        Self::Named(name)
    }
}

impl UnitAtom {
    pub(crate) fn term(
        &self,
        graph: &mut Graph,
    ) -> NodeId {
        match self {
            | Self::Named(name) => graph.sym(name),
            | Self::Literal(unit) => graph.lit(crate::graph::Payload::Unit(unit.clone())),
        }
    }
}

/// A unit as a product of atoms with rational exponents.
pub(crate) type UnitMap = BTreeMap<UnitAtom, Exp>;

/// Names of the seven base dimensions, in the order of [`Dim`].
pub(crate) const DIMENSIONS: [&str; 7] = [
    "Length",
    "Mass",
    "Time",
    "Current",
    "Temperature",
    "Amount",
    "Luminosity",
];

/// Exponents of the base dimensions a unit is measured in.
///
/// The seven SI bases are kept in a fixed array because they are on every
/// path; anything else — a currency, a count of items, whatever a host
/// registers — goes in a sorted side table. That split is what lets
/// `unit EUR : Currency` exist without giving up the constant-time common
/// case.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Dim {
    /// Exponents of the seven SI bases, in `DIMENSIONS` order.
    si: [Exp; 7],
    /// Exponents of every other base, sorted by name.
    extra: Vec<(String, Exp)>,
}

impl Dim {
    pub(crate) fn zero() -> Self {
        Self {
            si: [Exp::zero(); 7],
            extra: Vec::new(),
        }
    }

    pub(crate) fn of(e: [i32; 7]) -> Self {
        Self {
            si: e.map(Exp::from_integer),
            extra: Vec::new(),
        }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.si.iter().all(Zero::is_zero) && self.extra.is_empty()
    }

    pub(crate) fn add_scaled(
        &mut self,
        other: &Self,
        k: Exp,
    ) {
        for (a, b) in self.si.iter_mut().zip(&other.si) {
            *a += *b * k;
        }
        for (name, exponent) in &other.extra {
            self.add_extra(name, *exponent * k);
        }
    }

    /// Adds `exponent` to the exponent of the base dimension `name`.
    ///
    /// A name the seven SI bases do not have lands in the side table, which is
    /// how a custom base dimension reaches a unit at all.
    pub(crate) fn add_base(
        &mut self,
        name: &str,
        exponent: Exp,
    ) {
        if let Some(index) = DIMENSIONS.iter().position(|known| *known == name) {
            if let Some(slot) = self.si.get_mut(index) {
                *slot += exponent;
            }
            return;
        }
        self.add_extra(name, exponent);
    }

    /// Adds `exponent` to the side-table base `name`.
    fn add_extra(
        &mut self,
        name: &str,
        exponent: Exp,
    ) {
        if exponent.is_zero() {
            return;
        }
        match self
            .extra
            .binary_search_by(|(known, _)| known.as_str().cmp(name))
        {
            | Ok(index) => {
                if let Some(slot) = self.extra.get_mut(index) {
                    slot.1 += exponent;
                }
                self.extra.retain(|(_, e)| !e.is_zero());
            },
            | Err(index) => self.extra.insert(index, (name.to_owned(), exponent)),
        }
    }

    /// The seven SI exponents.
    pub(crate) const fn si(&self) -> &[Exp; 7] {
        &self.si
    }

    /// The non-SI base dimensions, sorted by name.
    pub(crate) fn extra(&self) -> &[(String, Exp)] {
        &self.extra
    }
}

/// An exact scale factor `rat * pi^pi`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Scale {
    /// The rational factor.
    pub rat: BigRational,
    /// The power of pi multiplying it.
    pub pi: i32,
}

impl Scale {
    /// The scale `1`.
    #[must_use]
    pub fn one() -> Self {
        Self {
            rat: BigRational::one(),
            pi: 0,
        }
    }

    /// Whether the scale is exactly `1`.
    #[must_use]
    pub fn is_one(&self) -> bool {
        self.rat.is_one() && self.pi == 0
    }

    /// The product of two scales.
    #[must_use]
    pub fn times(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            rat: &self.rat * &other.rat,
            pi: self.pi + other.pi,
        }
    }

    /// `self^k`.
    #[must_use]
    pub fn powi(
        &self,
        k: i32,
    ) -> Self {
        Self {
            rat: self.rat.pow(k),
            pi: self.pi * k,
        }
    }

    /// `self / other`.
    #[must_use]
    pub fn over(
        &self,
        other: &Self,
    ) -> Self {
        self.times(&other.powi(-1))
    }

    /// The scale as a float, when it is representable.
    #[must_use]
    pub fn to_f64(&self) -> Option<f64> {
        let value = self.rat.to_f64()?;
        Some(value * std::f64::consts::PI.powi(self.pi))
    }
}

/// What the table knows about one unit symbol.
#[derive(Clone, Debug)]
pub(crate) struct Resolved {
    pub(crate) dim: Dim,
    pub(crate) scale: Scale,
    /// `SI value = value * scale + offset` for affine units.
    pub(crate) offset: Option<BigRational>,
}

pub(crate) struct Def {
    dim: Dim,
    scale: Scale,
    offset: Option<BigRational>,
    prefixable: bool,
}

/// An exact decimal or fraction (`0.3048`, `1.602e-19`, `5/9`).
pub(crate) fn dec(text: &str) -> BigRational {
    if let Some((n, d)) = text.split_once('/') {
        return dec(n) / dec(d);
    }
    let (mantissa, exponent) = text
        .split_once('e')
        .map_or((text, 0_i32), |(m, e)| (m, e.parse().unwrap_or(0)));
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits: BigInt = format!("{whole}{fraction}")
        .parse()
        .unwrap_or_else(|_| BigInt::zero());
    let places = i32::try_from(fraction.len()).unwrap_or(0) - exponent;
    let ten = BigInt::from(10);
    if places >= 0 {
        BigRational::new(digits, ten.pow(places.unsigned_abs()))
    } else {
        BigRational::from_integer(digits * ten.pow(places.unsigned_abs()))
    }
}

pub(crate) fn table() -> &'static HashMap<&'static str, Def> {
    static TABLE: OnceLock<HashMap<&'static str, Def>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = HashMap::new();
        let mut add = |name: &'static str, dim: [i32; 7], scale: &str, prefixable: bool| {
            t.insert(
                name,
                Def {
                    dim: Dim::of(dim),
                    scale: Scale {
                        rat: dec(scale),
                        pi: 0,
                    },
                    offset: None,
                    prefixable,
                },
            );
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
        // The Aevia design standard (`spec/05-units.md` §5.3.3) declares these
        // and the table did not have them. The dalton is exact by the 2018
        // CODATA relation to the atomic mass constant (1 u = 1.66053906660e-27
        // kg); the knot and mile per hour are defined velocities; the acre is a
        // survey measure; an information bit is dimensionless, a byte is eight
        // bits and a kibibyte is 1024 bytes.
        add("Da", [0, 1, 0, 0, 0, 0, 0], "1.66053906660e-27", true);
        add("kn", [1, 0, -1, 0, 0, 0, 0], "1852/3600", false);
        add("mph", [1, 0, -1, 0, 0, 0, 0], "1609.344/3600", false);
        add("acre", [2, 0, 0, 0, 0, 0, 0], "4046.8564224", false);
        add("bit", [0; 7], "1", true);
        add("B", [0; 7], "8", true);
        add("KiB", [0; 7], "8192", false);
        add("eV", [2, 1, -2, 0, 0, 0, 0], "1.602176634e-19", true);
        add("cal", [2, 1, -2, 0, 0, 0, 0], "4.184", true);
        add("Wh", [2, 1, -2, 0, 0, 0, 0], "3600", true);
        add("lbf", [1, 1, -2, 0, 0, 0, 0], "4.4482216152605", false);
        // hp = 550 ft lbf / s
        add("hp", [2, 1, -3, 0, 0, 0, 0], "745.69987158227022", false);
        // Astronomy (IAU 2012 constants): the astronomical unit is exact, the
        // light-year follows from the defining speed of light, and the parsec
        // from the exact au and the IAU 2015 parallax definition.
        add("au", [1, 0, 0, 0, 0, 0, 0], "149597870700", false);
        add("ly", [1, 0, 0, 0, 0, 0, 0], "9460730472580800", false);
        add("pc", [1, 0, 0, 0, 0, 0, 0], "30856775814913673", true);
        add("Rsun", [1, 0, 0, 0, 0, 0, 0], "695700000", false);
        add("Msun", [0, 1, 0, 0, 0, 0, 0], "1.98892e30", false);
        add("Mearth", [0, 1, 0, 0, 0, 0, 0], "5.9722e24", false);
        add("jup", [0, 1, 0, 0, 0, 0, 0], "1.89813e27", false);
        // CGS: the Gaussian system, still the working unit set of much of
        // astronomy, geophysics and spectroscopy.
        add("dyn", [1, 1, -2, 0, 0, 0, 0], "1/100000", false);
        add("erg", [2, 1, -2, 0, 0, 0, 0], "1/10000000", false);
        add("Gal", [1, 0, -2, 0, 0, 0, 0], "1/100", false);
        add("P", [-1, 1, -1, 0, 0, 0, 0], "1/10", false);
        add("St", [2, 0, -1, 0, 0, 0, 0], "1/10000", false);
        add("G", [0, 1, -2, -1, 0, 0, 0], "1/10000", false);
        add("Mx", [2, 1, -2, -1, 0, 0, 0], "1/100000000", false);
        add("Ba", [-1, 1, -1, 0, 0, 0, 0], "1/10", false);
        // The oersted is `1000/(4 pi) A/m`, so its scale carries a power of pi
        // rather than a rational alone.
        t.insert(
            "Oe",
            Def {
                dim: Dim::of([-1, 0, 0, 1, 0, 0, 0]),
                scale: Scale {
                    rat: dec("250"),
                    pi: -1,
                },
                offset: None,
                prefixable: false,
            },
        );
        let angle = |name: &'static str, t: &mut HashMap<&'static str, Def>, denominator: &str| {
            t.insert(
                name,
                Def {
                    dim: Dim::zero(),
                    scale: Scale {
                        rat: dec(&format!("1/{denominator}")),
                        pi: 1,
                    },
                    offset: None,
                    prefixable: false,
                },
            );
        };
        angle("deg", &mut t, "180");
        angle("arcmin", &mut t, "10800");
        angle("arcsec", &mut t, "648000");
        t.insert(
            "degC",
            Def {
                dim: Dim::of([0, 0, 0, 0, 1, 0, 0]),
                scale: Scale::one(),
                offset: Some(dec("273.15")),
                prefixable: false,
            },
        );
        t.insert(
            "degF",
            Def {
                dim: Dim::of([0, 0, 0, 0, 1, 0, 0]),
                scale: Scale {
                    rat: dec("5/9"),
                    pi: 0,
                },
                offset: Some(dec("45967/180")),
                prefixable: false,
            },
        );
        t
    })
}

/// SI prefixes and their powers of ten.
pub(crate) const PREFIXES: [(&str, i32); 22] = [
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
pub(crate) const ALIASES: [(&str, &[(&str, i32)]); 21] = [
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
    // The standard names the dalton `Da` / `u` (§5.3.3); a byte is the
    // customary name of `B`.
    ("u", &[("Da", 1)]),
    ("byte", &[("B", 1)]),
];

pub(crate) fn resolve_symbol(name: &str) -> Option<Resolved> {
    // A unit the host registered wins over the built-in table: `register`
    // refuses a name that would shadow a different unit, so a hit here is the
    // unit the host meant.
    if let Some(resolved) = resolve_user(name) {
        return Some(resolved);
    }
    if let Some((_, parts)) = ALIASES.iter().find(|(alias, _)| *alias == name) {
        let map: UnitMap = parts
            .iter()
            .map(|&(n, e)| (n.to_owned().into(), Exp::from_integer(e)))
            .collect();
        return resolve_map(&map);
    }
    let table = table();
    if let Some(d) = table.get(name) {
        return Some(Resolved {
            dim: d.dim.clone(),
            scale: d.scale.clone(),
            offset: d.offset.clone(),
        });
    }
    for (prefix, power) in PREFIXES {
        let Some(rest) = name.strip_prefix(prefix) else {
            continue;
        };
        if let Some(d) = table.get(rest).filter(|d| d.prefixable) {
            let factor = Scale {
                rat: BigRational::from_integer(BigInt::from(10)).pow(power),
                pi: 0,
            };
            return Some(Resolved {
                dim: d.dim.clone(),
                scale: d.scale.times(&factor),
                offset: None,
            });
        }
    }
    None
}

/// Dimension, scale and affine offset of a unit product.
pub(crate) fn resolve_map(map: &UnitMap) -> Option<Resolved> {
    let mut dim = Dim::zero();
    let mut scale = Scale::one();
    let mut offset = None;
    for (name, &e) in map {
        let r = match name {
            | UnitAtom::Named(name) => resolve_symbol(name)?,
            | UnitAtom::Literal(unit) => resolve_value(unit)?,
        };
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

/// Resolves a by-value definition without interpreting its descriptive name.
/// Malformed ABI data is rejected before constructing a rational or entering
/// the rewrite engine. Offsets currently have an exact rational representation.
pub(crate) fn resolve_value(unit: &crate::kernel::value::UnitValue) -> Option<Resolved> {
    if unit.scale.denom <= BigInt::zero()
        || unit.scale.numer.is_zero()
        || unit.affine != unit.offset.is_some()
    {
        return None;
    }
    let offset = match &unit.offset {
        | Some(offset) => {
            if offset.denom <= BigInt::zero() || offset.pi != 0 {
                return None;
            }
            Some(BigRational::new(offset.numer.clone(), offset.denom.clone()))
        },
        | None => None,
    };
    Some(Resolved {
        dim: expr_to_dim(&unit.dimension).ok()?,
        scale: Scale {
            rat: BigRational::new(unit.scale.numer.clone(), unit.scale.denom.clone()),
            pi: unit.scale.pi,
        },
        offset,
    })
}

/// Maps a resolved unit to the dimension layer's term.
pub(crate) fn dimension_term(
    dim: &Dim,
    arena: &crate::graph::dim::DimArena,
) -> Option<crate::graph::dim::DimTerm> {
    let mut term = crate::graph::dim::DimTerm::dimensionless();
    for (index, exponent) in dim.si().iter().enumerate() {
        if exponent.is_zero() {
            continue;
        }
        let name = DIMENSIONS.get(index)?;
        let base = arena.base_id(name)?;
        let rational = BigRational::new(
            BigInt::from(*exponent.numer()),
            BigInt::from(*exponent.denom()),
        );
        term = term.mul(&crate::graph::dim::DimTerm::base(base, rational));
    }
    // A base the arena has never seen is not an error here: the graph defines
    // it when the host declares the dimension, and until then the unit simply
    // has no term form.
    for (name, exponent) in dim.extra() {
        let base = arena.base_id(name)?;
        let rational = BigRational::new(
            BigInt::from(*exponent.numer()),
            BigInt::from(*exponent.denom()),
        );
        term = term.mul(&crate::graph::dim::DimTerm::base(base, rational));
    }
    Some(term)
}

pub(crate) fn unit_info(
    resolved: &Resolved,
    arena: &crate::graph::dim::DimArena,
) -> Option<crate::graph::dim::UnitInfo> {
    Some(crate::graph::dim::UnitInfo {
        dimension: dimension_term(&resolved.dim, arena)?,
        scale: resolved.scale.to_f64()?,
        affine: resolved.offset.is_some(),
    })
}

/// Resolves a unit the host registered, by name or alias.
fn resolve_user(name: &str) -> Option<Resolved> {
    // Take a copy under the guard and let it go: nothing below needs the lock,
    // and holding it across the conversion would block a concurrent register.
    let def = {
        let extra = user_units().read().ok()?;
        extra
            .get(name)
            .or_else(|| {
                extra
                    .values()
                    .find(|unit| unit.aliases.iter().any(|a| a == name))
            })?
            .clone()
    };
    let dim = expr_to_dim(&def.dimension).ok()?;
    Some(Resolved {
        dim,
        scale: def.scale,
        offset: def.offset,
    })
}

/// Resolves the unit term of `quantity(value, unit)` for the dimension layer.
pub(crate) fn resolve_unit_node(
    graph: &Graph,
    node: NodeId,
    arena: &crate::graph::dim::DimArena,
) -> Option<crate::graph::dim::UnitInfo> {
    let children = graph.children(node);
    let unit = *children.get(1)?;
    let unit = parse_unit(graph, unit, 0)?;
    let resolved = resolve_map(&unit)?;
    unit_info(&resolved, arena)
}

/// Resolves a unit written as text (`"m/s"`, `"km/h"`).
pub(crate) fn resolve_unit_text(
    _graph: &Graph,
    text: &str,
    arena: &crate::graph::dim::DimArena,
) -> Option<crate::graph::dim::UnitInfo> {
    let unit = parse_unit_text(text)?;
    let resolved = resolve_map(&unit)?;
    unit_info(&resolved, arena)
}

/// Parses a unit product written as text.
pub(crate) fn parse_unit_text(text: &str) -> Option<UnitMap> {
    let text = text.trim();
    if text.is_empty() {
        return None;
    }
    let (numerator, denominator) = crate::graph::dim::split_unit_fraction(text);
    let mut map = UnitMap::new();
    for (part, sign) in [(numerator, 1_i32), (denominator, -1_i32)] {
        let part = part.trim();
        if part.is_empty() {
            if sign > 0 {
                return None;
            }
            continue;
        }
        if part == "1" {
            continue;
        }
        for factor in part.split('*') {
            let factor = factor.trim();
            let (name, exponent) = match factor.split_once('^') {
                | Some((name, exponent)) => (name.trim(), parse_unit_exponent(exponent)?),
                | None => (factor, Exp::new(1, 1)),
            };
            if name.is_empty() {
                return None;
            }
            let value = exponent * Exp::new(sign, 1);
            *map.entry(name.to_string().into())
                .or_insert_with(|| Exp::new(0, 1)) += value;
        }
    }
    if map.is_empty() {
        return None;
    }
    Some(clean(map))
}

pub(crate) fn parse_unit_exponent(text: &str) -> Option<Exp> {
    let text = text.trim();
    let text = text
        .strip_prefix('(')
        .and_then(|rest| rest.strip_suffix(')'))
        .unwrap_or(text);
    // Reject an empty or non-numeric exponent instead of letting `dec`
    // silently read it as zero (`m^` and `m^x` would otherwise be `m^0`).
    if text.is_empty()
        || text.matches('/').count() > 1
        || !text
            .chars()
            .all(|c| c.is_ascii_digit() || matches!(c, '+' | '-' | '.' | '/'))
    {
        return None;
    }
    let (numerator, denominator) = text.split_once('/').map_or((text, "1"), |(a, b)| (a, b));
    let denominator = dec(denominator);
    if denominator.is_zero() {
        return None;
    }
    let value = dec(numerator) / denominator;
    Some(Exp::new(value.numer().to_i32()?, value.denom().to_i32()?))
}

/// The unit product a term denotes.
pub(crate) fn parse_unit(
    graph: &Graph,
    node: NodeId,
    depth: usize,
) -> Option<UnitMap> {
    if depth > 16 {
        return None;
    }
    if let Some(n) = graph.number_of(node) {
        return n.is_one().then(UnitMap::new);
    }
    if let Some(symbol) = graph.symbol_of(node) {
        let mut map = UnitMap::new();
        map.insert(
            graph.interner().symbol_name(symbol).to_owned().into(),
            Exp::one(),
        );
        return Some(map);
    }
    if let Some(crate::graph::Payload::Unit(unit)) = graph.payload(node) {
        return Some([(UnitAtom::Literal(unit.clone()), Exp::one())].into());
    }
    // Several forms may share a class: take the first that reads as a unit.
    graph.enodes(graph.find(node)).find_map(|e| {
        let children = graph.children(e).to_vec();
        match graph.op(e) {
            | core::MUL => {
                let mut map = UnitMap::new();
                for c in children {
                    for (name, k) in parse_unit(graph, c, depth + 1)? {
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
                let map = parse_unit(graph, base, depth + 1)?;
                Some(clean(map.into_iter().map(|(n, k)| (n, k * e)).collect()))
            },
            | _ => None,
        }
    })
}

pub(crate) fn clean(mut map: UnitMap) -> UnitMap {
    map.retain(|_, e| !e.is_zero());
    map
}

// ---------------------------------------------------------------------------
// The public registry
// ---------------------------------------------------------------------------

/// One unit, as data.
///
/// The dimension is a [`DimExpr`] — names, not arena-local ids — so a
/// `UnitDef` can be written to a manifest, sent to another process and read
/// back without the arena it came from.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnitDef {
    /// Symbolic name (`m`, `km`, `degC`).
    pub name: String,
    /// Other names that resolve to this unit.
    pub aliases: Vec<String>,
    /// The dimension it measures.
    pub dimension: DimExpr,
    /// Exact scale to coherent SI units: `SI = value * scale + offset`.
    pub scale: Scale,
    /// Exact additive offset, for affine units.
    pub offset: Option<BigRational>,
    /// Whether SI prefixes may be attached.
    pub prefixable: bool,
    /// Whether the unit is affine (a point scale such as `degC`).
    pub affine: bool,
    /// How it is printed.
    pub display: String,
    /// Rough category (`length`, `time`, `derived`, ...).
    pub category: String,
    /// Version it appeared in.
    pub since: String,
}

impl UnitDef {
    /// A linear unit.
    #[must_use]
    pub fn linear(
        name: &str,
        dimension: DimExpr,
        scale: Scale,
    ) -> Self {
        Self {
            name: name.to_owned(),
            aliases: Vec::new(),
            dimension,
            scale,
            offset: None,
            prefixable: false,
            affine: false,
            display: name.to_owned(),
            category: String::new(),
            since: String::new(),
        }
    }

    /// Declares the unit affine, with `offset`.
    #[must_use]
    pub fn affine(
        mut self,
        offset: BigRational,
    ) -> Self {
        self.offset = Some(offset);
        self.affine = true;
        self
    }

    /// Allows SI prefixes.
    #[must_use]
    pub const fn prefixable(mut self) -> Self {
        self.prefixable = true;
        self
    }

    /// Adds an alias.
    #[must_use]
    pub fn alias(
        mut self,
        alias: &str,
    ) -> Self {
        self.aliases.push(alias.to_owned());
        self
    }

    /// Sets the printed form.
    #[must_use]
    pub fn display(
        mut self,
        display: &str,
    ) -> Self {
        display.clone_into(&mut self.display);
        self
    }

    /// Sets the category.
    #[must_use]
    pub fn category(
        mut self,
        category: &str,
    ) -> Self {
        category.clone_into(&mut self.category);
        self
    }

    /// The unit as the value type the dynamic ABI carries.
    #[must_use]
    pub fn to_value(&self) -> crate::kernel::value::UnitValue {
        crate::kernel::value::UnitValue::from_def(self)
    }
}

/// A unit expression resolved to a dimension and an exact scale.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedUnit {
    /// The dimension of the expression.
    pub dimension: DimExpr,
    /// Its exact scale to coherent SI units.
    pub scale: Scale,
    /// Its additive offset, when it is a single affine unit.
    pub offset: Option<BigRational>,
    /// Whether the expression is a point scale.
    pub affine: bool,
}

/// A snapshot of the registry, for caching and comparison.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnitRegistrySnapshot {
    /// Every unit, sorted by name.
    pub units: Vec<UnitDef>,
    /// BLAKE3 of the deterministic encoding of `units`.
    pub hash: [u8; 32],
}

/// The operations a host needs from a unit registry.
///
/// Implemented by [`StandardRegistry`], which reads the built-in table plus
/// anything [`register`]ed. The trait exists so that a host can supply its
/// own registry — a database of company-specific units, say — without
/// changing the dimension layer.
pub trait UnitRegistry {
    /// Every known unit, sorted by name.
    fn list(&self) -> Vec<UnitDef>;

    /// The unit named `name`, prefix or alias included.
    fn lookup(
        &self,
        name: &str,
    ) -> Option<UnitDef>;

    /// Resolves a unit expression (`"m/s"`, `"km/h"`, `"kg*m^2/s^3"`).
    ///
    /// # Errors
    /// Fails when a factor is not a known unit.
    fn resolve(
        &self,
        expr: &str,
    ) -> Result<ResolvedUnit, DimError>;

    /// Adds a unit.
    ///
    /// # Errors
    /// Fails when the name or an alias already belongs to another unit.
    fn register(
        &mut self,
        unit: UnitDef,
    ) -> Result<(), DimError>;

    /// Prints `value` with the coherent unit of `dim`, when there is one.
    fn format(
        &self,
        dim: &DimTerm,
        value: &crate::kernel::value::Value,
    ) -> String;

    /// A snapshot of the registry.
    fn snapshot(&self) -> UnitRegistrySnapshot;
}

/// The registry rssn ships: the built-in table plus user registrations.
#[derive(Copy, Clone, Debug, Default)]
pub struct StandardRegistry;

impl UnitRegistry for StandardRegistry {
    fn list(&self) -> Vec<UnitDef> {
        all()
    }

    fn lookup(
        &self,
        name: &str,
    ) -> Option<UnitDef> {
        lookup(name)
    }

    fn resolve(
        &self,
        expr: &str,
    ) -> Result<ResolvedUnit, DimError> {
        resolve(expr)
    }

    fn register(
        &mut self,
        unit: UnitDef,
    ) -> Result<(), DimError> {
        register(unit)
    }

    fn format(
        &self,
        dim: &DimTerm,
        value: &crate::kernel::value::Value,
    ) -> String {
        let arena = DimArena::new();
        format_quantity(dim, &arena, value)
    }

    fn snapshot(&self) -> UnitRegistrySnapshot {
        let units = all();
        UnitRegistrySnapshot {
            hash: crate::kernel::codec::digest(&units),
            units,
        }
    }
}

/// The units registered at run time, by name.
fn user_units() -> &'static RwLock<BTreeMap<String, UnitDef>> {
    static USER: OnceLock<RwLock<BTreeMap<String, UnitDef>>> = OnceLock::new();
    USER.get_or_init(|| RwLock::new(BTreeMap::new()))
}

/// Every known unit, sorted by name.
#[must_use]
pub fn all() -> Vec<UnitDef> {
    let mut out: Vec<UnitDef> = table().keys().map(|name| def_of(name)).collect::<Vec<_>>();
    if let Ok(extra) = user_units().read() {
        for unit in extra.values() {
            if let Some(slot) = out.iter_mut().find(|u| u.name == unit.name) {
                *slot = unit.clone();
            } else {
                out.push(unit.clone());
            }
        }
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    out
}

/// The unit named `name`, including prefixed forms and aliases.
#[must_use]
pub fn lookup(name: &str) -> Option<UnitDef> {
    if let Ok(extra) = user_units().read() {
        if let Some(unit) = extra.get(name) {
            return Some(unit.clone());
        }
        for unit in extra.values() {
            if unit.aliases.iter().any(|a| a == name) {
                return Some(unit.clone());
            }
        }
    }
    resolve_symbol(name).map(|resolved| {
        let mut unit = UnitDef::linear(name, dim_to_expr(&resolved.dim), resolved.scale.clone());
        if let Some(offset) = resolved.offset.clone() {
            unit = unit.affine(offset);
        }
        unit.category = category_of(&resolved.dim);
        "0.1".clone_into(&mut unit.since);
        unit
    })
}

/// Resolves a unit expression to its dimension and exact scale.
///
/// # Errors
/// Fails when a factor is not a known unit, or when an affine unit is
/// combined with anything else.
pub fn resolve(expr: &str) -> Result<ResolvedUnit, DimError> {
    let map = parse_unit_text(expr).ok_or_else(|| {
        DimError::UnknownUnit {
            unit: expr.to_owned(),
            context: None,
        }
    })?;
    let resolved = resolve_map(&map).ok_or_else(|| {
        DimError::UnknownUnit {
            unit: expr.to_owned(),
            context: None,
        }
    })?;
    Ok(ResolvedUnit {
        dimension: dim_to_expr(&resolved.dim),
        scale: resolved.scale,
        offset: resolved.offset.clone(),
        affine: resolved.offset.is_some(),
    })
}

/// Adds a unit to the process-wide registry.
///
/// # Errors
/// Fails when the name or one of the aliases already belongs to a different
/// unit, or when the dimension is not a product of base dimensions.
pub fn register(unit: UnitDef) -> Result<(), DimError> {
    if unit.name.trim().is_empty() {
        return Err(DimError::Parse {
            text: unit.name,
            reason: String::from("a unit needs a name"),
        });
    }
    expr_to_dim(&unit.dimension)?;
    let Ok(mut extra) = user_units().write() else {
        return Err(DimError::Parse {
            text: unit.name,
            reason: String::from("the unit registry is poisoned"),
        });
    };
    let clash = |other: &UnitDef| other.name != unit.name && other.dimension != unit.dimension;
    if let Some(existing) = extra.get(&unit.name)
        && clash(existing)
    {
        return Err(DimError::Parse {
            text: unit.name.clone(),
            reason: String::from("a different unit already has this name"),
        });
    }
    for alias in &unit.aliases {
        if extra.get(alias).is_some_and(clash) {
            return Err(DimError::Parse {
                text: alias.clone(),
                reason: String::from("a different unit already has this alias"),
            });
        }
    }
    extra.insert(unit.name.clone(), unit);
    Ok(())
}

/// Prints `value` with the coherent unit of `dim`.
#[must_use]
pub fn format_quantity(
    dim: &DimTerm,
    arena: &DimArena,
    value: &crate::kernel::value::Value,
) -> String {
    match coherent_symbol(dim, arena) {
        | Some(symbol) => format!("{value} {symbol}"),
        | None => format!("{value} {}", dim.describe(arena)),
    }
}

/// The symbol of the coherent unit of `dim`, when the table has one.
fn coherent_symbol(
    dim: &DimTerm,
    arena: &DimArena,
) -> Option<String> {
    let wanted = term_to_dim(dim, arena)?;
    let mut best: Option<&str> = None;
    for (name, def) in table() {
        if def.dim == wanted && def.scale.is_one() && def.offset.is_none() {
            let better = best.is_none_or(|current| name.len() < current.len());
            if better {
                best = Some(name);
            }
        }
    }
    best.map(str::to_owned)
}

/// The unit definition of a table entry.
fn def_of(name: &str) -> UnitDef {
    let Some(def) = table().get(name) else {
        return UnitDef::linear(name, DimExpr::Dimless, Scale::one());
    };
    let aliases: Vec<String> = ALIASES
        .iter()
        .filter(|(_, parts)| {
            parts.len() == 1 && parts.first().is_some_and(|&(n, e)| n == name && e == 1)
        })
        .map(|(alias, _)| (*alias).to_owned())
        .collect();
    UnitDef {
        name: name.to_owned(),
        aliases,
        dimension: dim_to_expr(&def.dim),
        scale: def.scale.clone(),
        offset: def.offset.clone(),
        prefixable: def.prefixable,
        affine: def.offset.is_some(),
        display: display_of(name),
        category: category_of(&def.dim),
        since: "0.1".to_owned(),
    }
}

/// How a unit is printed.
fn display_of(name: &str) -> String {
    match name {
        | "degC" => "\u{b0}C".to_owned(),
        | "degF" => "\u{b0}F".to_owned(),
        | "deg" => "\u{b0}".to_owned(),
        | other => other.to_owned(),
    }
}

/// A rough category for a dimension.
fn category_of(dim: &Dim) -> String {
    let non_zero: Vec<usize> = dim
        .si()
        .iter()
        .enumerate()
        .filter(|(_, e)| !e.is_zero())
        .map(|(i, _)| i)
        .collect();
    match (non_zero.as_slice(), dim.extra()) {
        | ([], []) => "dimensionless".to_owned(),
        | ([i], []) => {
            DIMENSIONS
                .get(*i)
                .map_or_else(|| "derived".to_owned(), |n| n.to_ascii_lowercase())
        },
        | ([], [(name, _)]) => name.to_ascii_lowercase(),
        | _ => "derived".to_owned(),
    }
}

/// The type-level dimension of an internal `Dim`.
pub(crate) fn dim_to_expr(dim: &Dim) -> DimExpr {
    let mut factors = Vec::new();
    for (index, exponent) in dim.si().iter().enumerate() {
        if exponent.is_zero() {
            continue;
        }
        let Some(name) = DIMENSIONS.get(index) else {
            continue;
        };
        let (numer, denom) = (*exponent.numer(), *exponent.denom());
        factors.push(DimExpr::base(name).pow(TypeExp::Rat(Rational::new(
            i64::from(numer),
            i64::from(denom),
        ))));
    }
    for (name, exponent) in dim.extra() {
        let (numer, denom) = (*exponent.numer(), *exponent.denom());
        factors.push(DimExpr::base(name).pow(TypeExp::Rat(Rational::new(
            i64::from(numer),
            i64::from(denom),
        ))));
    }
    match factors.len() {
        | 0 => DimExpr::Dimless,
        | 1 => factors.into_iter().next().unwrap_or(DimExpr::Dimless),
        | _ => DimExpr::Mul(factors),
    }
}

/// The internal `Dim` of a type-level dimension expression.
///
/// # Errors
/// Fails when the expression mentions a dimension variable, `any` or a
/// non-metric type: none of those is a concrete unit dimension.
fn expr_to_dim(dim: &DimExpr) -> Result<Dim, DimError> {
    fn walk(
        dim: &DimExpr,
        factors: &mut BTreeMap<String, BigRational>,
        power: &BigRational,
    ) -> Option<()> {
        match dim {
            | DimExpr::Dimless => {},
            | DimExpr::Base(name) => *factors.entry(name.clone()).or_default() += power,
            | DimExpr::Mul(parts) => {
                for part in parts {
                    walk(part, factors, power)?;
                }
            },
            | DimExpr::Pow(base, TypeExp::Rat(r)) => {
                if r.denom <= 0 {
                    return None;
                }
                let exponent = BigRational::new(r.numer.into(), r.denom.into());
                walk(base, factors, &(power * exponent))?;
            },
            | DimExpr::Inverse(base) => walk(base, factors, &-power)?,
            | _ => return None,
        }
        Some(())
    }
    let bad = || {
        DimError::Parse {
            text: crate::kernel::signature::dim_name(dim),
            reason:
                "unit dimensions must be concrete rational products within the unit exponent range"
                    .to_owned(),
        }
    };
    let mut factors = BTreeMap::new();
    walk(dim, &mut factors, &BigRational::one()).ok_or_else(bad)?;
    let mut out = Dim::zero();
    for (name, exponent) in factors {
        let numer = exponent.numer().to_i32().ok_or_else(bad)?;
        let denom = exponent.denom().to_i32().ok_or_else(bad)?;
        out.add_base(&name, Exp::new(numer, denom));
    }
    Ok(out)
}

/// The internal `Dim` of a dimension term, when it is a plain product of
/// base dimensions with rational exponents.
fn term_to_dim(
    term: &DimTerm,
    arena: &DimArena,
) -> Option<Dim> {
    let mut out = Dim::zero();
    for (name, exponent) in term.factors(arena) {
        let numer = exponent.numer().to_i32()?;
        let denom = exponent.denom().to_i32()?;
        out.add_base(&name, Exp::new(numer, denom));
    }
    term.vars(arena).is_empty().then_some(out)
}

impl crate::kernel::value::UnitValue {
    /// The value form of a unit definition.
    #[must_use]
    pub fn from_def(def: &UnitDef) -> Self {
        let (numer, denom) = (def.scale.rat.numer().clone(), def.scale.rat.denom().clone());
        Self {
            name: def.name.clone(),
            display: def.display.clone(),
            dimension: def.dimension.clone(),
            scale: crate::kernel::value::ScaleValue {
                numer,
                denom,
                pi: def.scale.pi,
            },
            offset: def.offset.as_ref().map(|offset| {
                let (numer, denom) = (offset.numer().clone(), offset.denom().clone());
                crate::kernel::value::ScaleValue { numer, denom, pi: 0 }
            }),
            affine: def.affine,
            aliases: def.aliases.clone(),
            prefixable: def.prefixable,
            category: def.category.clone(),
            since: def.since.clone(),
        }
    }
}
