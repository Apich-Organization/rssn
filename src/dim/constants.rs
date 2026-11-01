//! The physical constants, with their dimensions and exact values.
//!
//! # Why these live in the dimension layer
//!
//! A constant is a number **with a dimension**: `c` is `Length/Time`, `G` is
//! `Length³/(Mass·Time²)`. Storing it as a bare `f64` throws that away, and a
//! host that wants to generate `Units.Constants` then has to maintain a second
//! table of dimensions — exactly the duplication the unit registry exists to
//! prevent.
//!
//! So a constant is a [`Scale`] (an exact rational times a power of π, the same
//! representation a unit uses) plus a [`DimExpr`] plus a flag saying whether
//! the value is **exact by definition** or **measured**.
//!
//! # Exact and measured
//!
//! Since the 2019 SI revision, seven constants have exact values — `c`, `h`,
//! `e`, `k_B`, `N_A`, `Δν_Cs` and `K_cd` — and everything defined from them is
//! exact too (`R = N_A k_B`, `F = N_A e`, `Φ_0 = h/(2e)`, …). A measured
//! constant carries a standard uncertainty and is never exact, whatever the
//! digits look like. The two are distinguished here rather than left for the
//! caller to guess, because "6.67430e-11" reads like a fact and is not one.
//!
//! # Values
//!
//! CODATA 2022, as published by the BIPM. A constant whose value is an exact
//! decimal is stored as that decimal, so `c` really is `299792458` and not
//! `2.99792458e8` rounded into a float.
//!
//! # Not to be confused with the gauss
//!
//! `G` here is the gravitational constant. The unit `G` in the unit table is
//! the gauss. They are in different namespaces on purpose and a host should
//! keep them in different ones too.

use std::sync::OnceLock;

use num_bigint::BigInt;
use num_rational::BigRational;

use crate::dim::units::Scale;
use crate::dim::units::dec;
use crate::kernel::signature::DimExpr;
use crate::kernel::value::ConstantValue;
use crate::kernel::value::ScaleValue;

/// One physical constant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhysicalConstant {
    /// What it is called (`speed of light in vacuum`).
    pub name: &'static str,
    /// The symbol it is written with (`c`).
    pub symbol: &'static str,
    /// Its exact value, in coherent SI units.
    pub value: Scale,
    /// What it measures.
    pub dimension: DimExpr,
    /// The coherent SI unit it is expressed in (`m/s`).
    pub unit: &'static str,
    /// Whether the value is exact by definition.
    pub exact: bool,
    /// The standard uncertainty, when the value is measured.
    pub uncertainty: Option<Scale>,
    /// Where the value comes from.
    pub source: &'static str,
}

impl PhysicalConstant {
    /// The manifest form of this constant.
    #[must_use]
    pub fn to_value(&self) -> ConstantValue {
        ConstantValue {
            name: self.name.to_owned(),
            symbol: self.symbol.to_owned(),
            value: scale_value(&self.value),
            dimension: self.dimension.clone(),
            unit: self.unit.to_owned(),
            exact: self.exact,
            uncertainty: self.uncertainty.as_ref().map(scale_value),
            source: self.source.to_owned(),
        }
    }
}

/// A [`Scale`] in the manifest's wire form.
fn scale_value(scale: &Scale) -> ScaleValue {
    ScaleValue {
        numer: scale.rat.numer().clone(),
        denom: scale.rat.denom().clone(),
        pi: scale.pi,
    }
}

/// An exact value from a decimal string.
fn exact(text: &str) -> Scale {
    Scale {
        rat: dec(text),
        pi: 0,
    }
}

/// An exact value that carries a power of π.
fn exact_pi(
    text: &str,
    pi: i32,
) -> Scale {
    Scale {
        rat: dec(text),
        pi,
    }
}

/// A measured value with its standard uncertainty, both as decimals.
fn measured(
    text: &str,
    uncertainty: &str,
) -> (Scale, Option<Scale>) {
    (exact(text), Some(exact(uncertainty)))
}

/// The dimension `Length`.
fn length() -> DimExpr {
    DimExpr::base("Length")
}

/// A product of named bases.
fn dim(factors: &[(&str, i64)]) -> DimExpr {
    let mut out = Vec::with_capacity(factors.len());
    for (name, exponent) in factors {
        out.push(DimExpr::base(name).pow(crate::kernel::signature::Exp::from(*exponent)));
    }
    match out.len() {
        | 0 => DimExpr::Dimless,
        | 1 => out.into_iter().next().unwrap_or(DimExpr::Dimless),
        | _ => DimExpr::Mul(out),
    }
}

/// Every constant, sorted by symbol.
#[must_use]
pub fn all() -> &'static [PhysicalConstant] {
    static TABLE: OnceLock<Vec<PhysicalConstant>> = OnceLock::new();
    TABLE.get_or_init(build).as_slice()
}

/// The constant written `symbol`, if there is one.
#[must_use]
pub fn lookup(symbol: &str) -> Option<&'static PhysicalConstant> {
    all().iter().find(|constant| constant.symbol == symbol)
}

/// The constants with the given dimension.
#[must_use]
pub fn with_dimension(dimension: &DimExpr) -> Vec<&'static PhysicalConstant> {
    all()
        .iter()
        .filter(|constant| &constant.dimension == dimension)
        .collect()
}

/// Builds the table.
///
/// Sorted by symbol so that the manifest is deterministic without sorting
/// again, and so that a reader can find an entry by eye.
fn build() -> Vec<PhysicalConstant> {
    const CODATA: &str = "CODATA 2022";
    const SI: &str = "SI 2019 definition";
    let mut out = Vec::with_capacity(32);
    let mut add = |name: &'static str,
                   symbol: &'static str,
                   value: Scale,
                   dimension: DimExpr,
                   unit: &'static str,
                   exact: bool,
                   uncertainty: Option<Scale>,
                   source: &'static str| {
        out.push(PhysicalConstant {
            name,
            symbol,
            value,
            dimension,
            unit,
            exact,
            uncertainty,
            source,
        });
    };
    // The seven defining constants of the SI, all exact.
    add(
        "speed of light in vacuum",
        "c",
        exact("299792458"),
        dim(&[("Length", 1), ("Time", -1)]),
        "m/s",
        true,
        None,
        SI,
    );
    add(
        "Planck constant",
        "h",
        exact("6.62607015e-34"),
        dim(&[("Length", 2), ("Mass", 1), ("Time", -1)]),
        "J*s",
        true,
        None,
        SI,
    );
    add(
        "reduced Planck constant",
        "hbar",
        // h / (2 pi), exactly.
        exact_pi("3.313035075e-34", -1),
        dim(&[("Length", 2), ("Mass", 1), ("Time", -1)]),
        "J*s",
        true,
        None,
        SI,
    );
    add(
        "elementary charge",
        "e",
        exact("1.602176634e-19"),
        dim(&[("Time", 1), ("Current", 1)]),
        "C",
        true,
        None,
        SI,
    );
    add(
        "Boltzmann constant",
        "k_B",
        exact("1.380649e-23"),
        dim(&[("Length", 2), ("Mass", 1), ("Time", -2), ("Temperature", -1)]),
        "J/K",
        true,
        None,
        SI,
    );
    add(
        "Avogadro constant",
        "N_A",
        exact("6.02214076e23"),
        dim(&[("Amount", -1)]),
        "1/mol",
        true,
        None,
        SI,
    );
    add(
        "caesium hyperfine frequency",
        "dnu_Cs",
        exact("9192631770"),
        dim(&[("Time", -1)]),
        "Hz",
        true,
        None,
        SI,
    );
    add(
        "luminous efficacy",
        "K_cd",
        exact("683"),
        dim(&[("Length", -2), ("Mass", -1), ("Time", 3), ("Luminosity", 1)]),
        "lm/W",
        true,
        None,
        SI,
    );
    // Exact consequences of the seven.
    add(
        "molar gas constant",
        "R",
        exact("8.31446261815324"),
        dim(&[
            ("Length", 2),
            ("Mass", 1),
            ("Time", -2),
            ("Temperature", -1),
            ("Amount", -1),
        ]),
        "J/(mol*K)",
        true,
        None,
        SI,
    );
    add(
        "Faraday constant",
        "F",
        exact("96485.3321233100184"),
        dim(&[("Time", 1), ("Current", 1), ("Amount", -1)]),
        "C/mol",
        true,
        None,
        SI,
    );
    add(
        "magnetic flux quantum",
        "Phi_0",
        exact("2.067833848e-15"),
        dim(&[("Length", 2), ("Mass", 1), ("Time", -2), ("Current", -1)]),
        "Wb",
        true,
        None,
        SI,
    );
    add(
        "Josephson constant",
        "K_J",
        exact("483597.848416984e9"),
        dim(&[("Length", -2), ("Mass", -1), ("Time", 2), ("Current", 1)]),
        "Hz/V",
        true,
        None,
        SI,
    );
    add(
        "von Klitzing constant",
        "R_K",
        exact("25812.8074593045"),
        dim(&[("Length", 2), ("Mass", 1), ("Time", -3), ("Current", -2)]),
        "ohm",
        true,
        None,
        SI,
    );
    add(
        "conductance quantum",
        "G_0",
        exact("7.748091729e-5"),
        dim(&[("Length", -2), ("Mass", -1), ("Time", 3), ("Current", 2)]),
        "S",
        true,
        None,
        SI,
    );
    add(
        "molar mass constant",
        "M_u",
        exact("1e-3"),
        dim(&[("Mass", 1), ("Amount", -1)]),
        "kg/mol",
        true,
        None,
        SI,
    );
    add(
        "standard acceleration of gravity",
        "g_n",
        exact("9.80665"),
        dim(&[("Length", 1), ("Time", -2)]),
        "m/s^2",
        true,
        None,
        "CGPM 1901",
    );
    add(
        "standard atmosphere",
        "atm",
        exact("101325"),
        dim(&[("Length", -1), ("Mass", 1), ("Time", -2)]),
        "Pa",
        true,
        None,
        "CGPM 1954",
    );
    add(
        "triple point of water",
        "T_tp",
        exact("273.16"),
        dim(&[("Temperature", 1)]),
        "K",
        true,
        None,
        "CGPM 1954",
    );
    // Measured constants, with their standard uncertainties.
    let (g, g_u) = measured("6.67430e-11", "0.00015e-11");
    add(
        "Newtonian constant of gravitation",
        "G",
        g,
        dim(&[("Length", 3), ("Mass", -1), ("Time", -2)]),
        "m^3/(kg*s^2)",
        false,
        g_u,
        CODATA,
    );
    let (me, me_u) = measured("9.1093837139e-31", "0.0000000028e-31");
    add(
        "electron mass",
        "m_e",
        me,
        dim(&[("Mass", 1)]),
        "kg",
        false,
        me_u,
        CODATA,
    );
    let (mp, mp_u) = measured("1.67262192595e-27", "0.00000000052e-27");
    add(
        "proton mass",
        "m_p",
        mp,
        dim(&[("Mass", 1)]),
        "kg",
        false,
        mp_u,
        CODATA,
    );
    let (mn, mn_u) = measured("1.67492750056e-27", "0.00000000085e-27");
    add(
        "neutron mass",
        "m_n",
        mn,
        dim(&[("Mass", 1)]),
        "kg",
        false,
        mn_u,
        CODATA,
    );
    let (mu, mu_u) = measured("1.66053906892e-27", "0.00000000052e-27");
    add(
        "atomic mass constant",
        "m_u",
        mu,
        dim(&[("Mass", 1)]),
        "kg",
        false,
        mu_u,
        CODATA,
    );
    let (alpha, alpha_u) = measured("7.2973525643e-3", "0.0000000011e-3");
    add(
        "fine-structure constant",
        "alpha",
        alpha,
        DimExpr::Dimless,
        "1",
        false,
        alpha_u,
        CODATA,
    );
    let (a0, a0_u) = measured("5.29177210544e-11", "0.00000000082e-11");
    add(
        "Bohr radius",
        "a_0",
        a0,
        length(),
        "m",
        false,
        a0_u,
        CODATA,
    );
    let (ry, ry_u) = measured("10973731.568157", "0.000012");
    add(
        "Rydberg constant",
        "R_inf",
        ry,
        dim(&[("Length", -1)]),
        "1/m",
        false,
        ry_u,
        CODATA,
    );
    let (mub, mub_u) = measured("9.2740100657e-24", "0.0000000029e-24");
    add(
        "Bohr magneton",
        "mu_B",
        mub,
        dim(&[("Length", 2), ("Current", 1)]),
        "J/T",
        false,
        mub_u,
        CODATA,
    );
    let (sigma, sigma_u) = measured("5.670374419e-8", "0.000000000e-8");
    add(
        "Stefan-Boltzmann constant",
        "sigma_SB",
        sigma,
        dim(&[
            ("Mass", 1),
            ("Time", -3),
            ("Temperature", -4),
        ]),
        "W/(m^2*K^4)",
        false,
        sigma_u,
        CODATA,
    );
    let (eps0, eps0_u) = measured("8.8541878188e-12", "0.0000000014e-12");
    add(
        "vacuum electric permittivity",
        "epsilon_0",
        eps0,
        dim(&[
            ("Length", -3),
            ("Mass", -1),
            ("Time", 4),
            ("Current", 2),
        ]),
        "F/m",
        false,
        eps0_u,
        CODATA,
    );
    let (mu0, mu0_u) = measured("1.25663706127e-6", "0.00000000020e-6");
    add(
        "vacuum magnetic permeability",
        "mu_0",
        mu0,
        dim(&[("Length", 1), ("Mass", 1), ("Time", -2), ("Current", -2)]),
        "N/A^2",
        false,
        mu0_u,
        CODATA,
    );
    // Keep the table sorted by symbol: the manifest is hashed, and a stable
    // order is what makes two builds agree.
    out.sort_by(|a, b| a.symbol.cmp(b.symbol));
    out
}

/// A constant as a quantity in a `Value`.
///
/// The value is exact when the constant is, which is the point of storing it
/// as a rational: `c` is an integer, not `2.99792458e8`.
#[must_use]
pub fn to_quantity(
    constant: &PhysicalConstant,
) -> crate::kernel::value::QuantityValue {
    let magnitude = if constant.value.pi == 0 {
        crate::kernel::value::Value::rational(constant.value.rat.clone())
    } else {
        // A power of π is not rational; evaluate it to the nearest float and
        // say so by using a float value rather than an exact one.
        crate::kernel::value::Value::Float(
            rational_to_f64(&constant.value.rat) * std::f64::consts::PI.powi(constant.value.pi),
        )
    };
    crate::kernel::value::QuantityValue {
        value: Box::new(magnitude),
        unit: crate::kernel::value::UnitValue {
            name: constant.symbol.to_owned(),
            display: constant.symbol.to_owned(),
            dimension: constant.dimension.clone(),
            scale: scale_value(&constant.value),
            offset: None,
            affine: false,
            aliases: Vec::new(),
            prefixable: false,
            category: String::from("constant"),
            since: String::new(),
        },
    }
}

/// A rational as the nearest `f64`.
fn rational_to_f64(value: &BigRational) -> f64 {
    use num_traits::ToPrimitive;
    let numer = value.numer().to_f64().unwrap_or(f64::NAN);
    let denom = value.denom().to_f64().unwrap_or(1.0);
    numer / denom
}

/// The value of `constant` as a `f64`, for numeric use.
#[must_use]
pub fn value_of(constant: &PhysicalConstant) -> f64 {
    rational_to_f64(&constant.value.rat) * std::f64::consts::PI.powi(constant.value.pi)
}

/// The uncertainty of `constant` as a `f64`, or zero when it is exact.
#[must_use]
pub fn uncertainty_of(constant: &PhysicalConstant) -> f64 {
    constant
        .uncertainty
        .as_ref()
        .map_or(0.0, |scale| rational_to_f64(&scale.rat))
}

/// An exact integer, for tests and for callers that want the defining value.
#[must_use]
pub fn exact_integer(constant: &PhysicalConstant) -> Option<BigInt> {
    (constant.value.pi == 0 && constant.value.rat.is_integer())
        .then(|| constant.value.rat.to_integer())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_defining_constants_are_exact() {
        for symbol in ["c", "h", "e", "k_B", "N_A", "dnu_Cs", "K_cd"] {
            let constant = lookup(symbol).unwrap_or_else(|| panic!("{symbol}"));
            assert!(constant.exact, "{symbol} must be exact");
            assert!(constant.uncertainty.is_none(), "{symbol} has an uncertainty");
        }
    }

    #[test]
    fn the_speed_of_light_is_the_integer_it_is_defined_to_be() {
        let c = lookup("c").unwrap_or_else(|| panic!("c"));
        assert_eq!(exact_integer(c), Some(BigInt::from(299_792_458)));
        assert_eq!(value_of(c), 299_792_458.0);
    }

    #[test]
    fn a_measured_constant_carries_its_uncertainty() {
        let g = lookup("G").unwrap_or_else(|| panic!("G"));
        assert!(!g.exact);
        let uncertainty = uncertainty_of(g);
        assert!(uncertainty > 0.0, "G must have a non-zero uncertainty");
        // CODATA 2022: 6.67430(15)e-11.
        assert!((value_of(g) / 6.674_30e-11 - 1.0).abs() < 1e-12);
        assert!((uncertainty / 1.5e-15 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn the_derived_constants_are_consistent_with_their_definitions() {
        // R = N_A k_B
        let r = value_of(lookup("R").unwrap_or_else(|| panic!("R")));
        let na = value_of(lookup("N_A").unwrap_or_else(|| panic!("N_A")));
        let kb = value_of(lookup("k_B").unwrap_or_else(|| panic!("k_B")));
        assert!((r / (na * kb) - 1.0).abs() < 1e-14, "{r}");
        // F = N_A e
        let f = value_of(lookup("F").unwrap_or_else(|| panic!("F")));
        let e = value_of(lookup("e").unwrap_or_else(|| panic!("e")));
        assert!((f / (na * e) - 1.0).abs() < 1e-14, "{f}");
        // hbar = h / (2 pi)
        let hbar = value_of(lookup("hbar").unwrap_or_else(|| panic!("hbar")));
        let h = value_of(lookup("h").unwrap_or_else(|| panic!("h")));
        assert!((hbar / (h / (2.0 * std::f64::consts::PI)) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn the_table_is_sorted_and_unique() {
        let symbols: Vec<&str> = all().iter().map(|c| c.symbol).collect();
        let mut sorted = symbols.clone();
        sorted.sort_unstable();
        assert_eq!(symbols, sorted, "the table must be sorted by symbol");
        let mut unique = symbols.clone();
        unique.dedup();
        assert_eq!(symbols.len(), unique.len(), "two constants share a symbol");
    }

    #[test]
    fn every_constant_has_a_dimension_and_a_unit() {
        for constant in all() {
            assert!(!constant.unit.is_empty(), "{} has no unit", constant.symbol);
            assert!(
                !matches!(constant.dimension, DimExpr::Any | DimExpr::NonMetric),
                "{} has a non-physical dimension",
                constant.symbol
            );
        }
    }

    #[test]
    fn a_constant_becomes_a_quantity_with_its_dimension() {
        let c = lookup("c").unwrap_or_else(|| panic!("c"));
        let quantity = to_quantity(c);
        assert_eq!(quantity.unit.dimension, c.dimension);
        assert_eq!(quantity.value.as_f64(), Some(299_792_458.0));
    }
}
