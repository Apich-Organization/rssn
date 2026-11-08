//! Unit systems: natural, Planck, atomic and geometrized units.
//!
//! # What a unit system is
//!
//! A unit system is a set of physical quantities that are **declared equal
//! to one**. Setting `c = 1` does not change any physics; it identifies the
//! dimension of `c` (`Length/Time`) with the dimensionless one, which merges
//! `Length` and `Time` into one class. Setting `c = ħ = k_B = ε₀ = 1` merges
//! every mechanical, thermal and electromagnetic dimension into powers of a
//! single one, energy — the particle physicist's "everything is in eV".
//!
//! Formally the unity quantities `q₁ … q_k` span a sublattice `K` of the
//! dimension group, and the system works in the quotient `D/K`. Each class
//! of the quotient is named by a **representative** built from the system's
//! *generators* (energy for natural units, length for geometrized units,
//! nothing at all for Planck units). The decomposition
//!
//! ```text
//! d = Σ λᵢ · dim(qᵢ) + Σ μⱼ · gⱼ
//! ```
//!
//! is unique when the unity dimensions and the generators form a basis of
//! the seven SI base dimensions, which [`UnitSystem::validate`] checks. Then
//!
//! * the **representative** of `d` is `Π gⱼ^μⱼ` — an ordinary SI dimension;
//! * a quantity `x` of dimension `d` has the system value
//!   `x / Π qᵢ^λᵢ`, which is an ordinary SI quantity of the representative
//!   dimension.
//!
//! The second point is the design decision that keeps natural units sound
//! in a dimension-typed host: a value "in natural units" is not a new kind
//! of number. `1 m` in natural units is `1 m / (ħ c)`, an SI quantity of
//! dimension `1/J`; a host that types it as `Energy^-1` is typing it
//! correctly in SI too. So values can cross between a natural-units module
//! and an SI module without any reinterpretation — the boundary only has to
//! refuse dimensions that are *not* representatives, which ordinary
//! dimension checking already does.
//!
//! # Exactness
//!
//! The exponents `λᵢ` are rationals: a charge in Heaviside–Lorentz natural
//! units is `Q / (ε₀ ħ c)^{1/2}`. The conversion factor is therefore
//! returned as factored powers of the unity quantities ([`Reduction::factors`])
//! rather than as a [`Scale`], which can only hold a rational times a power of
//! π. [`Reduction::exact_scale`] gives the [`Scale`] when every exponent is an
//! integer; a host with a symbolic layer builds the general factor from the
//! factored form, and [`Reduction::factor_f64`] evaluates it numerically.

use std::sync::OnceLock;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::dim::SI_BASES;
use crate::dim::constants;
use crate::dim::units::Scale;
use crate::kernel::signature::DimExpr;
use crate::kernel::signature::Exp;
use crate::kernel::signature::Rational;
use crate::kernel::value::ScaleValue;

/// One quantity a system sets to one: `multiplier · constant`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unity {
    /// The symbol of the physical constant ([`constants::lookup`]).
    pub symbol: &'static str,
    /// An exact numeric multiplier (`4π` for `4π ε₀ = 1`).
    pub multiplier: Scale,
    /// How the system writes it (`4 pi epsilon_0`).
    pub label: &'static str,
}

/// One generator of a system's representatives.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Generator {
    /// What the generator measures (`Energy`).
    pub name: &'static str,
    /// Its SI dimension.
    pub dimension: DimExpr,
    /// The unit values of this generator are displayed in (`eV`).
    pub display: &'static str,
}

/// A unit system.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnitSystem {
    /// The name a host selects it by (`natural`).
    pub name: &'static str,
    /// One-line description.
    pub description: &'static str,
    /// The quantities set to one.
    pub unity: Vec<Unity>,
    /// The generators of the representatives.
    pub generators: Vec<Generator>,
}

/// The decomposition of one SI dimension in a system.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reduction {
    /// The representative: an SI dimension built from the generators.
    pub dimension: DimExpr,
    /// The exponent `μⱼ` of each generator, by generator name; zero
    /// exponents are omitted.
    pub generators: Vec<(String, BigRational)>,
    /// The exponent `λᵢ` of each unity quantity, by constant symbol; zero
    /// exponents are omitted. `x_SI = x_system · Π (mᵢ qᵢ)^λᵢ`.
    pub factors: Vec<(String, BigRational)>,
}

/// Why a dimension or a system cannot be reduced.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SystemError {
    /// No system has this name.
    UnknownSystem(String),
    /// A unity quantity names a constant the table does not have.
    UnknownConstant(String),
    /// The unity dimensions and generators are not a basis.
    NotABasis(String),
    /// The dimension has a variable, a symbolic exponent, or a non-SI base.
    NotReducible(String),
}

impl std::fmt::Display for SystemError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::UnknownSystem(name) => write!(f, "no unit system is called `{name}`"),
            | Self::UnknownConstant(symbol) => write!(f, "no physical constant `{symbol}`"),
            | Self::NotABasis(name) => {
                write!(
                    f,
                    "unit system `{name}`: its unity quantities and generators are not a basis of \
                     the SI dimensions"
                )
            },
            | Self::NotReducible(what) => write!(f, "cannot reduce the dimension {what}"),
        }
    }
}

impl std::error::Error for SystemError {}

/// The wire form of a unit system, carried in the kernel manifest.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnitSystemValue {
    /// The name a host selects it by.
    pub name: String,
    /// One-line description.
    pub description: String,
    /// `(constant symbol, exact multiplier, label)` for each quantity set to
    /// one.
    pub unity: Vec<(String, ScaleValue, String)>,
    /// `(name, SI dimension, display unit)` for each generator.
    pub generators: Vec<(String, DimExpr, String)>,
}

type Vector = [BigRational; 7];

fn zero_vector() -> Vector {
    std::array::from_fn(|_| BigRational::zero())
}

/// The SI exponent vector of a dimension expression.
///
/// # Errors
/// [`SystemError::NotReducible`] for a variable, a symbolic exponent, a
/// non-SI base, or a non-metric dimension.
pub fn si_vector(dimension: &DimExpr) -> Result<Vector, SystemError> {
    let mut out = zero_vector();
    accumulate(dimension, &BigRational::one(), &mut out)?;
    Ok(out)
}

fn accumulate(
    dimension: &DimExpr,
    power: &BigRational,
    out: &mut Vector,
) -> Result<(), SystemError> {
    match dimension {
        | DimExpr::Dimless => Ok(()),
        | DimExpr::Base(name) => {
            let index = SI_BASES
                .iter()
                .position(|base| base == name)
                .ok_or_else(|| SystemError::NotReducible(format!("with base `{name}`")))?;
            out[index] += power;
            Ok(())
        },
        | DimExpr::Mul(factors) => {
            for factor in factors {
                accumulate(factor, power, out)?;
            }
            Ok(())
        },
        | DimExpr::Pow(base, Exp::Rat(exponent)) => {
            let exponent = BigRational::new(exponent.numer.into(), exponent.denom.into());
            accumulate(base, &(power * exponent), out)
        },
        | DimExpr::Inverse(inner) => accumulate(inner, &-power, out),
        | other => Err(SystemError::NotReducible(format!("{other:?}"))),
    }
}

/// The dimension expression of an SI exponent vector.
#[must_use]
pub fn from_si_vector(vector: &Vector) -> DimExpr {
    let mut factors = Vec::new();
    for (name, exponent) in SI_BASES.iter().zip(vector) {
        if exponent.is_zero() {
            continue;
        }
        let base = DimExpr::base(name);
        factors.push(if exponent.is_one() {
            base
        } else {
            base.pow(Exp::Rat(rational(exponent)))
        });
    }
    match factors.len() {
        | 0 => DimExpr::Dimless,
        | 1 => factors.pop().unwrap_or(DimExpr::Dimless),
        | _ => DimExpr::Mul(factors),
    }
}

/// A small rational for the wire form; system exponents are tiny.
fn rational(value: &BigRational) -> Rational {
    Rational::new(
        value.numer().to_i64().unwrap_or(0),
        value.denom().to_i64().unwrap_or(1),
    )
}

/// Solves `A x = b` over the rationals, `A` given by columns.
fn solve(
    columns: &[Vector],
    rhs: &Vector,
) -> Option<Vec<BigRational>> {
    let n = 7;
    if columns.len() != n {
        return None;
    }
    // Row-major augmented matrix.
    let mut m: Vec<Vec<BigRational>> = (0..n)
        .map(|row| {
            let mut line: Vec<BigRational> = columns.iter().map(|c| c[row].clone()).collect();
            line.push(rhs[row].clone());
            line
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n).find(|&row| !m[row][col].is_zero())?;
        m.swap(col, pivot);
        let lead = m[col][col].clone();
        for value in &mut m[col] {
            *value /= &lead;
        }
        for row in 0..n {
            if row != col && !m[row][col].is_zero() {
                let factor = m[row][col].clone();
                for k in col..=n {
                    let delta = &factor * &m[col][k];
                    m[row][k] -= delta;
                }
            }
        }
    }
    Some(
        m.into_iter()
            .map(|mut row| row.pop().unwrap_or_default())
            .collect(),
    )
}

impl UnitSystem {
    /// The SI dimension of each unity quantity, in order.
    fn unity_vectors(&self) -> Result<Vec<Vector>, SystemError> {
        self.unity
            .iter()
            .map(|unity| {
                let constant = constants::lookup(unity.symbol)
                    .ok_or_else(|| SystemError::UnknownConstant(unity.symbol.to_owned()))?;
                si_vector(&constant.dimension)
            })
            .collect()
    }

    fn columns(&self) -> Result<Vec<Vector>, SystemError> {
        let mut columns = self.unity_vectors()?;
        for generator in &self.generators {
            columns.push(si_vector(&generator.dimension)?);
        }
        Ok(columns)
    }

    /// Checks that the unity quantities and the generators form a basis.
    ///
    /// # Errors
    /// [`SystemError::NotABasis`] or [`SystemError::UnknownConstant`].
    pub fn validate(&self) -> Result<(), SystemError> {
        let columns = self.columns()?;
        // A basis solves every unit vector.
        for index in 0..7 {
            let mut target = zero_vector();
            target[index] = BigRational::one();
            if solve(&columns, &target).is_none() {
                return Err(SystemError::NotABasis(self.name.to_owned()));
            }
        }
        Ok(())
    }

    /// Decomposes an SI dimension into its representative and the powers of
    /// the unity quantities.
    ///
    /// # Errors
    /// [`SystemError::NotReducible`] for a dimension that is not a plain SI
    /// product, [`SystemError::NotABasis`] for an inconsistent system.
    pub fn reduce(
        &self,
        dimension: &DimExpr,
    ) -> Result<Reduction, SystemError> {
        let target = si_vector(dimension)?;
        self.reduce_vector(&target)
    }

    /// [`UnitSystem::reduce`] for an exponent vector.
    ///
    /// # Errors
    /// As [`UnitSystem::reduce`].
    pub fn reduce_vector(
        &self,
        target: &Vector,
    ) -> Result<Reduction, SystemError> {
        let columns = self.columns()?;
        let solution =
            solve(&columns, target).ok_or_else(|| SystemError::NotABasis(self.name.to_owned()))?;
        let (lambdas, mus) = solution.split_at(self.unity.len());
        let mut representative = zero_vector();
        let mut generators = Vec::new();
        for (generator, mu) in self.generators.iter().zip(mus) {
            if mu.is_zero() {
                continue;
            }
            let vector = si_vector(&generator.dimension)?;
            for (slot, component) in representative.iter_mut().zip(vector) {
                *slot += component * mu;
            }
            generators.push((generator.name.to_owned(), mu.clone()));
        }
        let factors = self
            .unity
            .iter()
            .zip(lambdas)
            .filter(|(_, lambda)| !lambda.is_zero())
            .map(|(unity, lambda)| (unity.symbol.to_owned(), lambda.clone()))
            .collect();
        Ok(Reduction {
            dimension: from_si_vector(&representative),
            generators,
            factors,
        })
    }

    /// The multiplier of a unity quantity, by constant symbol.
    fn multiplier(
        &self,
        symbol: &str,
    ) -> Option<&Scale> {
        self.unity
            .iter()
            .find(|unity| unity.symbol == symbol)
            .map(|unity| &unity.multiplier)
    }

    /// The manifest form of this system.
    #[must_use]
    pub fn to_value(&self) -> UnitSystemValue {
        UnitSystemValue {
            name: self.name.to_owned(),
            description: self.description.to_owned(),
            unity: self
                .unity
                .iter()
                .map(|unity| {
                    (
                        unity.symbol.to_owned(),
                        ScaleValue {
                            numer: unity.multiplier.rat.numer().clone(),
                            denom: unity.multiplier.rat.denom().clone(),
                            pi: unity.multiplier.pi,
                        },
                        unity.label.to_owned(),
                    )
                })
                .collect(),
            generators: self
                .generators
                .iter()
                .map(|g| (g.name.to_owned(), g.dimension.clone(), g.display.to_owned()))
                .collect(),
        }
    }
}

impl Reduction {
    /// The conversion factor `Π (mᵢ qᵢ)^λᵢ` as an exact [`Scale`], when every
    /// exponent is an integer.
    #[must_use]
    pub fn exact_scale(
        &self,
        system: &UnitSystem,
    ) -> Option<Scale> {
        let mut scale = Scale::one();
        for (symbol, lambda) in &self.factors {
            if !lambda.is_integer() {
                return None;
            }
            let power = lambda.to_integer().to_i32()?;
            let constant = constants::lookup(symbol)?;
            let base = system
                .multiplier(symbol)
                .map_or_else(|| constant.value.clone(), |m| m.times(&constant.value));
            scale = scale.times(&base.powi(power));
        }
        Some(scale)
    }

    /// The conversion factor as a float.
    #[must_use]
    pub fn factor_f64(
        &self,
        system: &UnitSystem,
    ) -> f64 {
        let mut log = 0.0_f64;
        let mut sign_negative = false;
        for (symbol, lambda) in &self.factors {
            let Some(constant) = constants::lookup(symbol) else {
                return f64::NAN;
            };
            let base = system
                .multiplier(symbol)
                .map_or_else(|| constant.value.clone(), |m| m.times(&constant.value));
            let value = base.to_f64().unwrap_or(f64::NAN);
            if value < 0.0 {
                sign_negative = !sign_negative;
            }
            log += value.abs().ln() * lambda.to_f64().unwrap_or(f64::NAN);
        }
        let magnitude = log.exp();
        if sign_negative {
            -magnitude
        } else {
            magnitude
        }
    }

    /// Whether the system leaves this dimension unchanged (no factors).
    #[must_use]
    pub fn is_identity(&self) -> bool {
        self.factors.is_empty()
    }
}

fn four_pi() -> Scale {
    Scale {
        rat: BigRational::from_integer(BigInt::from(4)),
        pi: 1,
    }
}

fn base(name: &str) -> DimExpr {
    DimExpr::base(name)
}

fn energy() -> DimExpr {
    DimExpr::Mul(vec![
        base("Length").pow(Exp::from(2)),
        base("Mass"),
        base("Time").pow(Exp::from(-2)),
    ])
}

fn unity(
    symbol: &'static str,
    label: &'static str,
) -> Unity {
    Unity {
        symbol,
        multiplier: Scale::one(),
        label,
    }
}

fn build() -> Vec<UnitSystem> {
    let amount = Generator {
        name: "Amount",
        dimension: base("Amount"),
        display: "mol",
    };
    let luminosity = Generator {
        name: "Luminosity",
        dimension: base("Luminosity"),
        display: "cd",
    };
    let mut systems = vec![
        UnitSystem {
            name: "natural",
            description: "particle-physics natural units: c = hbar = k_B = epsilon_0 = 1 \
                          (Heaviside-Lorentz); every mechanical, thermal and electromagnetic \
                          dimension is a power of energy",
            unity: vec![
                unity("c", "c"),
                unity("hbar", "hbar"),
                unity("k_B", "k_B"),
                unity("epsilon_0", "epsilon_0"),
            ],
            generators: vec![
                Generator {
                    name: "Energy",
                    dimension: energy(),
                    display: "eV",
                },
                amount.clone(),
                luminosity.clone(),
            ],
        },
        UnitSystem {
            name: "planck",
            description: "Planck units: c = hbar = G = k_B = 4 pi epsilon_0 = 1 (Gaussian \
                          charge); mechanical, thermal and electromagnetic quantities are pure \
                          numbers",
            unity: vec![
                unity("c", "c"),
                unity("hbar", "hbar"),
                unity("G", "G"),
                unity("k_B", "k_B"),
                Unity {
                    symbol: "epsilon_0",
                    multiplier: four_pi(),
                    label: "4 pi epsilon_0",
                },
            ],
            generators: vec![amount.clone(), luminosity.clone()],
        },
        UnitSystem {
            name: "atomic",
            description: "Hartree atomic units: hbar = m_e = e = 4 pi epsilon_0 = 1; lengths \
                          are in bohr and energies in hartree",
            unity: vec![
                unity("hbar", "hbar"),
                unity("m_e", "m_e"),
                unity("e", "e"),
                Unity {
                    symbol: "epsilon_0",
                    multiplier: four_pi(),
                    label: "4 pi epsilon_0",
                },
            ],
            generators: vec![
                Generator {
                    name: "Temperature",
                    dimension: base("Temperature"),
                    display: "K",
                },
                amount.clone(),
                luminosity.clone(),
            ],
        },
        UnitSystem {
            name: "geometrized",
            description: "geometrized units of general relativity: c = G = 1; mass and time are \
                          lengths",
            unity: vec![unity("c", "c"), unity("G", "G")],
            generators: vec![
                Generator {
                    name: "Length",
                    dimension: base("Length"),
                    display: "m",
                },
                Generator {
                    name: "Current",
                    dimension: base("Current"),
                    display: "A",
                },
                Generator {
                    name: "Temperature",
                    dimension: base("Temperature"),
                    display: "K",
                },
                amount,
                luminosity,
            ],
        },
    ];
    systems.sort_by_key(|system| system.name);
    systems
}

/// Every unit system, sorted by name.
#[must_use]
pub fn all() -> &'static [UnitSystem] {
    static SYSTEMS: OnceLock<Vec<UnitSystem>> = OnceLock::new();
    SYSTEMS.get_or_init(build)
}

/// The system called `name`.
#[must_use]
pub fn lookup(name: &str) -> Option<&'static UnitSystem> {
    all().iter().find(|system| system.name == name)
}

/// Decomposes `dimension` in the system called `name`.
///
/// # Errors
/// [`SystemError::UnknownSystem`] or the errors of [`UnitSystem::reduce`].
pub fn reduce(
    name: &str,
    dimension: &DimExpr,
) -> Result<Reduction, SystemError> {
    lookup(name)
        .ok_or_else(|| SystemError::UnknownSystem(name.to_owned()))?
        .reduce(dimension)
}

/// Converts an SI value of `dimension` to the system called `name`.
///
/// # Errors
/// As [`reduce`].
pub fn to_system(
    name: &str,
    dimension: &DimExpr,
    value: f64,
) -> Result<f64, SystemError> {
    let system = lookup(name).ok_or_else(|| SystemError::UnknownSystem(name.to_owned()))?;
    Ok(value / system.reduce(dimension)?.factor_f64(system))
}

/// Converts a value in the system called `name` back to SI `dimension`.
///
/// # Errors
/// As [`reduce`].
pub fn from_system(
    name: &str,
    dimension: &DimExpr,
    value: f64,
) -> Result<f64, SystemError> {
    let system = lookup(name).ok_or_else(|| SystemError::UnknownSystem(name.to_owned()))?;
    Ok(value * system.reduce(dimension)?.factor_f64(system))
}

#[allow(clippy::unwrap_used, clippy::panic)]
#[cfg(test)]
mod tests {
    use super::*;

    fn rat(
        n: i64,
        d: i64,
    ) -> BigRational {
        BigRational::new(n.into(), d.into())
    }

    fn close(
        a: f64,
        b: f64,
        tol: f64,
    ) -> bool {
        ((a - b) / b).abs() < tol
    }

    #[test]
    fn every_system_is_a_basis() {
        for system in all() {
            system.validate().unwrap_or_else(|error| panic!("{error}"));
        }
        assert_eq!(all().len(), 4);
    }

    #[test]
    fn natural_units_reduce_everything_to_energy() {
        let natural = lookup("natural").unwrap();
        let length = natural.reduce(&base("Length")).unwrap();
        assert_eq!(length.generators, vec![("Energy".to_owned(), rat(-1, 1))]);
        assert_eq!(
            length.factors,
            vec![("c".to_owned(), rat(1, 1)), ("hbar".to_owned(), rat(1, 1))]
        );
        let mass = natural.reduce(&base("Mass")).unwrap();
        assert_eq!(mass.generators, vec![("Energy".to_owned(), rat(1, 1))]);
        assert_eq!(mass.factors, vec![("c".to_owned(), rat(-2, 1))]);
        let time = natural.reduce(&base("Time")).unwrap();
        assert_eq!(time.factors, vec![("hbar".to_owned(), rat(1, 1))]);
        // Velocity is dimensionless: v/c.
        let velocity = natural
            .reduce(&DimExpr::Mul(vec![
                base("Length"),
                base("Time").pow(Exp::from(-1)),
            ]))
            .unwrap();
        assert!(velocity.generators.is_empty());
        assert_eq!(velocity.dimension, DimExpr::Dimless);
        // Charge is dimensionless in Heaviside-Lorentz units, with a
        // half-integer factor (epsilon_0 hbar c)^(1/2).
        let charge = natural
            .reduce(&DimExpr::Mul(vec![base("Current"), base("Time")]))
            .unwrap();
        assert_eq!(charge.dimension, DimExpr::Dimless);
        assert!(charge.factors.iter().all(|(_, l)| *l == rat(1, 2)));
        assert!(charge.exact_scale(natural).is_none());
    }

    #[test]
    fn the_electron_mass_is_half_an_mev() {
        let natural = lookup("natural").unwrap();
        let m_e = constants::value_of(constants::lookup("m_e").unwrap());
        let joules = to_system("natural", &base("Mass"), m_e).unwrap();
        let ev = constants::value_of(constants::lookup("e").unwrap());
        assert!(
            close(joules / ev, 0.510_998_950_69e6, 1e-9),
            "{}",
            joules / ev
        );
        // The exact scale of a mass is c^-2.
        let scale = natural
            .reduce(&base("Mass"))
            .unwrap()
            .exact_scale(natural)
            .unwrap();
        assert_eq!(
            scale.rat,
            BigRational::from_integer(299_792_458_i64.into()).pow(-2)
        );
    }

    #[test]
    fn a_fermi_is_about_one_over_197_mev() {
        // hbar c = 197.3269804 MeV fm.
        let in_joule_inverse = to_system("natural", &base("Length"), 1e-15).unwrap();
        let ev = constants::value_of(constants::lookup("e").unwrap());
        let per_mev = in_joule_inverse * ev * 1e6;
        assert!(close(per_mev, 1.0 / 197.326_980_4, 1e-9), "{per_mev}");
        let back = from_system("natural", &base("Length"), in_joule_inverse).unwrap();
        assert!(close(back, 1e-15, 1e-14));
    }

    #[test]
    fn the_fine_structure_constant_is_e_squared_over_four_pi() {
        // In Heaviside-Lorentz natural units alpha = e^2 / (4 pi).
        let charge = DimExpr::Mul(vec![base("Current"), base("Time")]);
        let e = constants::value_of(constants::lookup("e").unwrap());
        let e_natural = to_system("natural", &charge, e).unwrap();
        let alpha = e_natural * e_natural / (4.0 * std::f64::consts::PI);
        assert!(close(alpha, 7.297_352_5643e-3, 1e-8), "{alpha}");
    }

    #[test]
    fn planck_and_atomic_units() {
        // The Planck length is sqrt(hbar G / c^3) ~ 1.616255e-35 m.
        let one = from_system("planck", &base("Length"), 1.0).unwrap();
        assert!(close(one, 1.616_255e-35, 1e-5), "{one}");
        assert_eq!(
            lookup("planck")
                .unwrap()
                .reduce(&base("Mass"))
                .unwrap()
                .dimension,
            DimExpr::Dimless
        );
        // The bohr radius is 5.29177210544e-11 m and the hartree 4.3597447222e-18 J.
        let bohr = from_system("atomic", &base("Length"), 1.0).unwrap();
        assert!(close(bohr, 5.291_772_105_44e-11, 1e-9), "{bohr}");
        let hartree = from_system("atomic", &energy(), 1.0).unwrap();
        assert!(close(hartree, 4.359_744_722_2e-18, 1e-9), "{hartree}");
    }

    #[test]
    fn geometrized_units_make_mass_a_length() {
        let geometrized = lookup("geometrized").unwrap();
        let mass = geometrized.reduce(&base("Mass")).unwrap();
        assert_eq!(mass.dimension, base("Length"));
        // The Sun's mass is about 1.477 km.
        let sun = to_system("geometrized", &base("Mass"), 1.988_47e30).unwrap();
        assert!(close(sun, 1476.6, 1e-3), "{sun}");
    }

    #[test]
    fn a_dimension_variable_is_not_reducible() {
        assert!(matches!(
            reduce("natural", &DimExpr::Var("D".into())),
            Err(SystemError::NotReducible(_))
        ));
        assert!(matches!(
            reduce("nope", &DimExpr::Dimless),
            Err(SystemError::UnknownSystem(_))
        ));
    }

    #[test]
    fn manifest_form_round_trips_through_json() {
        let value = lookup("atomic").unwrap().to_value();
        let text = serde_json::to_string(&value).unwrap();
        let back: UnitSystemValue = serde_json::from_str(&text).unwrap();
        assert_eq!(back, value);
        assert_eq!(back.unity.len(), 4);
    }
}
