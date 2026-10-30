//! Exponent expressions: rational functions over atoms.
//!
//! An exponent of a base dimension is not always a number. `x^(a+b)` has
//! exponent `a + b`, `x^(1/n)` has exponent `1/n`, and `product(f(k), k=1..n)`
//! has exponent `n`. The exponent ring is therefore the field of rational
//! functions over `Q` in *atoms*: symbols and opaque term forms that cannot be
//! reduced to arithmetic. Equality is decided by cross-multiplication, so no
//! polynomial gcd is needed; the polynomials are tiny (a handful of terms).

use std::collections::BTreeMap;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

/// Identifier of an atom of the exponent field (a symbol or an opaque term).
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct AtomId(pub(crate) u32);

/// A monomial: atoms with non-zero integer exponents, sorted by id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Monomial {
    powers: Vec<(AtomId, i32)>,
}

impl Monomial {
    const fn one() -> Self {
        Self { powers: Vec::new() }
    }

    fn atom(atom: AtomId) -> Self {
        Self {
            powers: vec![(atom, 1)],
        }
    }

    const fn is_one(&self) -> bool {
        self.powers.is_empty()
    }

    fn mul(
        &self,
        other: &Self,
    ) -> Self {
        let mut powers = self.powers.clone();
        for &(atom, exponent) in &other.powers {
            match powers.binary_search_by_key(&atom, |&(a, _)| a) {
                | Ok(index) => {
                    if let Some(slot) = powers.get_mut(index) {
                        slot.1 = slot.1.saturating_add(exponent);
                    }
                },
                | Err(index) => powers.insert(index, (atom, exponent)),
            }
        }
        powers.retain(|&(_, exponent)| exponent != 0);
        Self { powers }
    }

    /// The componentwise minimum of two monomials: their common factor.
    fn gcd(
        &self,
        other: &Self,
    ) -> Self {
        let mut powers = Vec::new();
        for &(atom, exponent) in &self.powers {
            if let Ok(index) = other.powers.binary_search_by_key(&atom, |&(a, _)| a)
                && let Some(&(_, other_exponent)) = other.powers.get(index)
            {
                let common = exponent.min(other_exponent);
                if common != 0 {
                    powers.push((atom, common));
                }
            }
        }
        Self { powers }
    }

    fn display(
        &self,
        atoms: &[String],
    ) -> String {
        if self.powers.is_empty() {
            return String::from("1");
        }
        let mut factors = Vec::with_capacity(self.powers.len());
        for &(atom, exponent) in &self.powers {
            let name = atoms.get(atom.0 as usize).map_or("?", String::as_str);
            if exponent == 1 {
                factors.push(name.to_string());
            } else {
                factors.push(format!("{name}^{exponent}"));
            }
        }
        factors.join("*")
    }
}

/// A polynomial over `Q` in the atoms.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Poly {
    terms: BTreeMap<Monomial, BigRational>,
}

impl Poly {
    fn zero() -> Self {
        Self::default()
    }

    fn one() -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(Monomial::one(), BigRational::one());
        Self { terms }
    }

    fn constant(value: BigRational) -> Self {
        if value.is_zero() {
            return Self::zero();
        }
        let mut terms = BTreeMap::new();
        terms.insert(Monomial::one(), value);
        Self { terms }
    }

    fn atom(atom: AtomId) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(Monomial::atom(atom), BigRational::one());
        Self { terms }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn is_one(&self) -> bool {
        self.terms.len() == 1
            && self
                .terms
                .iter()
                .next()
                .is_some_and(|(mono, coeff)| mono.is_one() && coeff.is_one())
    }

    fn add(
        &self,
        other: &Self,
    ) -> Self {
        let mut terms = self.terms.clone();
        for (monomial, coefficient) in &other.terms {
            match terms.get_mut(monomial) {
                | Some(slot) => {
                    *slot += coefficient;
                    if slot.is_zero() {
                        terms.remove(monomial);
                    }
                },
                | None => {
                    terms.insert(monomial.clone(), coefficient.clone());
                },
            }
        }
        Self { terms }
    }

    fn neg(&self) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(mono, coeff)| (mono.clone(), -coeff))
                .collect(),
        }
    }

    fn mul(
        &self,
        other: &Self,
    ) -> Self {
        let mut terms: BTreeMap<Monomial, BigRational> = BTreeMap::new();
        for (a_mono, a_coeff) in &self.terms {
            for (b_mono, b_coeff) in &other.terms {
                let monomial = a_mono.mul(b_mono);
                let coefficient = a_coeff * b_coeff;
                match terms.get_mut(&monomial) {
                    | Some(slot) => {
                        *slot += &coefficient;
                        if slot.is_zero() {
                            terms.remove(&monomial);
                        }
                    },
                    | None => {
                        terms.insert(monomial, coefficient);
                    },
                }
            }
        }
        Self { terms }
    }

    fn scalar_mul(
        &self,
        factor: &BigRational,
    ) -> Self {
        if factor.is_zero() {
            return Self::zero();
        }
        Self {
            terms: self
                .terms
                .iter()
                .map(|(mono, coeff)| (mono.clone(), coeff * factor))
                .collect(),
        }
    }

    fn leading_coefficient(&self) -> Option<&BigRational> {
        self.terms.values().next_back()
    }

    /// The monomial dividing every term, if the polynomial is non-zero.
    fn common_monomial(&self) -> Option<Monomial> {
        let mut iter = self.terms.keys();
        let first = iter.next()?.clone();
        Some(iter.fold(first, |acc, monomial| acc.gcd(monomial)))
    }

    fn divide_monomial(
        &self,
        divisor: &Monomial,
    ) -> Self {
        if divisor.is_one() {
            return self.clone();
        }
        let mut terms = BTreeMap::new();
        for (monomial, coefficient) in &self.terms {
            let mut powers = monomial.powers.clone();
            for &(atom, exponent) in &divisor.powers {
                if let Ok(index) = powers.binary_search_by_key(&atom, |&(a, _)| a)
                    && let Some(slot) = powers.get_mut(index)
                {
                    slot.1 -= exponent;
                }
            }
            powers.retain(|&(_, e)| e != 0);
            terms.insert(Monomial { powers }, coefficient.clone());
        }
        Self { terms }
    }

    fn display(
        &self,
        atoms: &[String],
    ) -> String {
        if self.terms.is_empty() {
            return String::from("0");
        }
        let mut parts = Vec::with_capacity(self.terms.len());
        for (monomial, coefficient) in self.terms.iter().rev() {
            let number = if coefficient.is_integer() {
                coefficient.to_integer().to_string()
            } else {
                format!(
                    "{}/{}",
                    coefficient.numer(),
                    coefficient.denom()
                )
            };
            if monomial.is_one() {
                parts.push(number);
            } else if coefficient.is_one() {
                parts.push(monomial.display(atoms));
            } else if coefficient == &BigRational::from_integer(BigInt::from(-1)) {
                parts.push(format!("-{}", monomial.display(atoms)));
            } else {
                parts.push(format!("{number}*{}", monomial.display(atoms)));
            }
        }
        // Join with explicit signs where possible.
        let mut out = String::new();
        for (index, part) in parts.iter().enumerate() {
            if index == 0 {
                out.push_str(part);
            } else if let Some(rest) = part.strip_prefix('-') {
                out.push_str(" - ");
                out.push_str(rest);
            } else {
                out.push_str(" + ");
                out.push_str(part);
            }
        }
        out
    }
}

/// An exponent: a rational function `num / den` over the atoms.
#[derive(Clone, Debug)]
pub(crate) struct Expo {
    num: Poly,
    den: Poly,
}

impl Expo {
    pub(crate) fn zero() -> Self {
        Self {
            num: Poly::zero(),
            den: Poly::one(),
        }
    }

    pub(crate) fn one() -> Self {
        Self {
            num: Poly::one(),
            den: Poly::one(),
        }
    }

    pub(crate) fn from_rat(value: BigRational) -> Self {
        Self {
            num: Poly::constant(value),
            den: Poly::one(),
        }
    }

    pub(crate) fn atom(atom: AtomId) -> Self {
        Self {
            num: Poly::atom(atom),
            den: Poly::one(),
        }
    }

    fn normalize(
        mut num: Poly,
        mut den: Poly,
    ) -> Self {
        if num.is_zero() {
            return Self::zero();
        }
        // A zero denominator cannot occur from the constructors; guard anyway.
        if den.is_zero() {
            return Self::zero();
        }
        // Scalar normalization: make the leading denominator coefficient one.
        if let Some(leading) = den.leading_coefficient().cloned()
            && !leading.is_zero()
            && !leading.is_one()
        {
            let inverse = leading.recip();
            num = num.scalar_mul(&inverse);
            den = den.scalar_mul(&inverse);
        }
        // Cancel a common monomial factor.
        if let (Some(num_common), Some(den_common)) = (num.common_monomial(), den.common_monomial())
        {
            let common = num_common.gcd(&den_common);
            if !common.is_one() {
                num = num.divide_monomial(&common);
                den = den.divide_monomial(&common);
            }
        }
        Self { num, den }
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.num.is_zero()
    }

    #[cfg(test)]
    pub(crate) fn is_one(&self) -> bool {
        self.den.is_one() && self.num.is_one()
    }

    pub(crate) fn is_constant(&self) -> Option<BigRational> {
        if !self.den.is_one() {
            return None;
        }
        if self.num.terms.is_empty() {
            return Some(BigRational::zero());
        }
        if self.num.terms.len() != 1 {
            return None;
        }
        let (monomial, coefficient) = self.num.terms.iter().next()?;
        if monomial.is_one() {
            Some(coefficient.clone())
        } else {
            None
        }
    }

    pub(crate) fn as_integer(&self) -> Option<BigInt> {
        let constant = self.is_constant()?;
        if constant.is_integer() {
            Some(constant.to_integer())
        } else {
            None
        }
    }

    pub(crate) fn add(
        &self,
        other: &Self,
    ) -> Self {
        let num = self.num.mul(&other.den).add(&other.num.mul(&self.den));
        let den = self.den.mul(&other.den);
        Self::normalize(num, den)
    }

    pub(crate) fn neg(&self) -> Self {
        Self {
            num: self.num.neg(),
            den: self.den.clone(),
        }
    }

    pub(crate) fn sub(
        &self,
        other: &Self,
    ) -> Self {
        self.add(&other.neg())
    }

    pub(crate) fn mul(
        &self,
        other: &Self,
    ) -> Self {
        Self::normalize(
            self.num.mul(&other.num),
            self.den.mul(&other.den),
        )
    }

    #[cfg(test)]
    pub(crate) fn inverse(&self) -> Option<Self> {
        if self.num.is_zero() {
            return None;
        }
        Some(Self::normalize(self.den.clone(), self.num.clone()))
    }

    pub(crate) fn pow_i32(
        &self,
        exponent: i32,
    ) -> Option<Self> {
        if exponent == 0 {
            return Some(Self::one());
        }
        let magnitude = exponent.unsigned_abs();
        // A dimension with a huge literal exponent is not worth expanding;
        // treat it as an opaque atom by refusing the exact power.
        if magnitude > 256 {
            return None;
        }
        let mut num = Poly::one();
        let mut den = Poly::one();
        for _ in 0..magnitude {
            num = num.mul(&self.num);
            den = den.mul(&self.den);
        }
        if exponent > 0 {
            Some(Self::normalize(num, den))
        } else {
            if num.is_zero() {
                return None;
            }
            Some(Self::normalize(den, num))
        }
    }

    /// Equality in the field of rational functions, by cross-multiplication.
    pub(crate) fn equals(
        &self,
        other: &Self,
    ) -> bool {
        self.num.mul(&other.den) == other.num.mul(&self.den)
    }

    pub(crate) fn display(
        &self,
        atoms: &[String],
    ) -> String {
        if self.den.is_one() {
            return self.num.display(atoms);
        }
        if self.num.is_one() {
            return format!("1/({})", self.den.display(atoms));
        }
        format!(
            "({})/({})",
            self.num.display(atoms),
            self.den.display(atoms)
        )
    }
}

impl PartialEq for Expo {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.equals(other)
    }
}

impl Eq for Expo {}