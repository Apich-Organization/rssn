//! Partial fraction decomposition over the rationals.
//!
//! `n / d = q + sum_i sum_{j <= m_i} a_ij / f_i^j`, where `d = c prod f_i^m_i`
//! is the factorisation of the denominator into irreducible factors over `Q`
//! and every numerator `a_ij` has degree below `deg f_i`. For each factor
//! power `F_i = f_i^m_i`, the polynomial Chinese remainder theorem gives the
//! local numerator `r * (d/F_i)^(-1) mod F_i`, where `r` is the remainder of
//! `n` divided by `d`. Repeated division by `f_i` gives the individual
//! numerators for that block. Large decompositions can compute these blocks
//! independently on the active CPU pool.
//!
//! The full irreducible factorisation handles repeated factors and higher
//! degrees uniformly. `1/(x^2 (x + 1))` yields `-1/x + 1/x^2 + 1/(x + 1)`.
//! `1/(x^3 + 1)` yields a constant over `x + 1` plus a linear numerator over
//! `x^2 - x + 1`.
//!
//! `split_real_quadratics` goes one step further and rewrites such a
//! decomposition over `R`, adjoining `sqrt(discriminant)` for each real
//! irreducible quadratic factor.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::kernels::par;

use super::univariate;
use super::univariate::QPoly;

/// One term `numerator / factor^power` of a decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Piece {
    /// An irreducible factor of the denominator, primitive over `Z`.
    pub factor: QPoly,
    /// Its power in this term.
    pub power: u32,
    /// A numerator of degree below `deg factor` (never zero).
    pub numerator: QPoly,
}

/// A complete decomposition of a rational function.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Apart {
    /// The polynomial part.
    pub quotient: QPoly,
    /// The proper fractions.
    pub pieces: Vec<Piece>,
}

/// A term `(real + radical*sqrt(radicand)) / (x - pole)^power` of a
/// decomposition taken over a real quadratic extension of `Q`.
///
/// The pole is `offset + scale*sqrt(radicand)`; for a term with a rational
/// numerator and a rational pole `radicand`, `radical` and `scale` are zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RadicalPiece {
    /// Rational part of the numerator.
    pub real: BigRational,
    /// Coefficient of `sqrt(radicand)` in the numerator.
    pub radical: BigRational,
    /// Positive rational whose square root is adjoined; zero for a term
    /// that is entirely rational.
    pub radicand: BigRational,
    /// Rational part of the pole.
    pub offset: BigRational,
    /// Coefficient of `sqrt(radicand)` in the pole.
    pub scale: BigRational,
    /// Power of the linear denominator.
    pub power: u32,
}

/// Decomposes `numer / denom` over the rationals.
///
/// Input coefficients are in
/// ascending degree. The result gives the polynomial quotient and nonzero
/// fraction numerators, in factor order and increasing denominator power.
/// Returns `None` for a zero denominator or an incomplete factorisation.
#[must_use]
pub fn apart(
    numer: &[BigRational],
    denom: &[BigRational],
) -> Option<Apart> {
    let trim = |p: &[BigRational]| {
        let mut p = p.to_vec();
        while p.last().is_some_and(Zero::is_zero) {
            p.pop();
        }
        p
    };
    let (numer, denom) = (trim(numer), trim(denom));
    if denom.is_empty() {
        return None;
    }
    let (quotient, remainder) = univariate::divrem(&numer, &denom)?;
    if remainder.is_empty() || denom.len() < 2 {
        return Some(Apart {
            quotient,
            pieces: Vec::new(),
        });
    }
    let (_, factors) = univariate::factor(&denom);
    let factors: Vec<_> = factors
        .into_iter()
        .map(|(factor, power)| {
            (
                factor.into_iter().map(BigRational::from_integer).collect(),
                power,
            )
        })
        .collect();
    let pieces = decompose_proper(&remainder, &denom, &factors)?;
    Some(Apart { quotient, pieces })
}

fn decompose_proper(
    remainder: &[BigRational],
    denominator: &[BigRational],
    factors: &[(QPoly, u32)],
) -> Option<Vec<Piece>> {
    const FACTOR_GRAIN: usize = 4;
    let degree = denominator.len().checked_sub(1)?;
    let factored_degree = factors
        .iter()
        .try_fold(0_usize, |total, (factor, multiplicity)| {
            let factor_degree = factor.len().checked_sub(1)?;
            if factor_degree == 0 || *multiplicity == 0 {
                return None;
            }
            total.checked_add(factor_degree.checked_mul(usize::try_from(*multiplicity).ok()?)?)
        })?;
    if factored_degree != degree {
        return None;
    }

    if degree.saturating_mul(factors.len()) >= 192 {
        let blocks = par::try_map_indexed(factors.len(), FACTOR_GRAIN, |index| {
            let (factor, multiplicity) = factors.get(index).ok_or(())?;
            factor_block(remainder, denominator, factor, *multiplicity).ok_or(())
        })
        .ok()?;
        return Some(blocks.into_iter().flatten().collect());
    }
    let mut pieces = Vec::new();
    for (factor, multiplicity) in factors {
        pieces.extend(factor_block(remainder, denominator, factor, *multiplicity)?);
    }
    Some(pieces)
}

fn factor_block(
    remainder: &[BigRational],
    denominator: &[BigRational],
    factor: &[BigRational],
    multiplicity: u32,
) -> Option<Vec<Piece>> {
    let power = polynomial_power(factor, multiplicity);
    let (cofactor, defect) = univariate::divrem(denominator, &power)?;
    if !defect.is_empty() {
        return None;
    }
    let (_, cofactor) = univariate::divrem(&cofactor, &power)?;
    let inverse = inverse_modulo(&cofactor, &power)?;
    let (_, reduced_remainder) = univariate::divrem(remainder, &power)?;
    let (_, mut local) =
        univariate::divrem(&univariate::mul(&reduced_remainder, &inverse), &power)?;

    let mut pieces = Vec::new();
    for exponent in (1..=multiplicity).rev() {
        if local.is_empty() {
            break;
        }
        let (quotient, numerator) = univariate::divrem(&local, factor)?;
        local = quotient;
        if !numerator.is_empty() {
            pieces.push(Piece {
                factor: factor.to_vec(),
                power: exponent,
                numerator,
            });
        }
    }
    pieces.reverse();
    Some(pieces)
}

fn polynomial_power(
    factor: &[BigRational],
    exponent: u32,
) -> QPoly {
    match exponent {
        | 0 => vec![BigRational::one()],
        | 1 => factor.to_vec(),
        | _ => {
            let half = polynomial_power(factor, exponent / 2);
            let squared = univariate::mul(&half, &half);
            if exponent.is_multiple_of(2) {
                squared
            } else {
                univariate::mul(&squared, factor)
            }
        },
    }
}

/// Inverts a reduced residue modulo `modulus` with the Euclidean algorithm.
fn inverse_modulo(
    residue: &[BigRational],
    modulus: &[BigRational],
) -> Option<QPoly> {
    let (mut previous, mut current) = (modulus.to_vec(), residue.to_vec());
    let (mut previous_weight, mut current_weight) = (Vec::new(), vec![BigRational::one()]);
    while current.len() > 1 {
        let (quotient, remainder) = univariate::divrem(&previous, &current)?;
        let next_weight = univariate::sub(
            &previous_weight,
            &univariate::mul(&quotient, &current_weight),
        );
        (previous, current) = (current, remainder);
        (previous_weight, current_weight) = (current_weight, next_weight);
    }
    let [unit] = current.as_slice() else {
        return None;
    };
    Some(
        current_weight
            .into_iter()
            .map(|value| value / unit)
            .collect(),
    )
}

/// Rewrites the proper fractions of `parts` as terms with linear
/// denominators, adjoining `sqrt(discriminant)` for every real irreducible
/// quadratic factor.
///
/// `1/(x^2 - 2)` decomposes as `1/(x^2 - 2)` over `Q` and here becomes
/// `(sqrt(2)/4)/(x - sqrt(2)) - (sqrt(2)/4)/(x + sqrt(2))`, carried exactly as
/// `real + radical*sqrt(radicand)` with `radicand = 8`.
///
/// Returns `None`, rather than a wrong answer, when a factor has degree three
/// or more, when a quadratic factor has no real roots (it is irreducible over
/// `R` as well), or when a real quadratic factor occurs with power two or more
/// — a higher-order pole over the extension needs the full expansion there,
/// which this does not attempt. The polynomial part `parts.quotient` is
/// unchanged and is not repeated in the result.
#[must_use]
pub fn split_real_quadratics(parts: &Apart) -> Option<Vec<RadicalPiece>> {
    let mut out = Vec::new();
    for piece in &parts.pieces {
        match piece.factor.len() {
            | 2 => {
                let lead = piece.factor.get(1)?.clone();
                let constant = piece.factor.first()?.clone();
                out.push(RadicalPiece {
                    real: piece
                        .numerator
                        .first()
                        .cloned()
                        .unwrap_or_else(BigRational::zero),
                    radical: BigRational::zero(),
                    radicand: BigRational::zero(),
                    offset: -constant / lead,
                    scale: BigRational::zero(),
                    power: piece.power,
                });
            },
            | 3 => {
                if piece.power != 1 {
                    return None;
                }
                let a = piece.factor.get(2)?.clone();
                let b = piece.factor.get(1)?.clone();
                let c = piece.factor.first()?.clone();
                let two_a = BigRational::from_integer(BigInt::from(2)) * &a;
                let disc = &b * &b - BigRational::from_integer(BigInt::from(4)) * &a * &c;
                if disc <= BigRational::zero() {
                    return None;
                }
                let alpha = piece
                    .numerator
                    .get(1)
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                let beta = piece
                    .numerator
                    .first()
                    .cloned()
                    .unwrap_or_else(BigRational::zero);
                // Residues at the two roots are conjugates: `p +- q sqrt(disc)`
                // with `p = alpha / (2a)` and `q = (2a beta - alpha b)/(2a disc)`.
                let real = &alpha / &two_a;
                let radical = (&two_a * &beta - &alpha * &b) / (&two_a * &disc);
                let offset = -b / &two_a;
                let scale = BigRational::one() / &two_a;
                out.push(RadicalPiece {
                    real: real.clone(),
                    radical: radical.clone(),
                    radicand: disc.clone(),
                    offset: offset.clone(),
                    scale: scale.clone(),
                    power: 1,
                });
                out.push(RadicalPiece {
                    real,
                    radical: -radical,
                    radicand: disc,
                    offset,
                    scale: -scale,
                    power: 1,
                });
            },
            | _ => return None,
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;

    fn q(values: &[i64]) -> QPoly {
        values
            .iter()
            .map(|&v| BigRational::from_integer(BigInt::from(v)))
            .collect()
    }

    fn frac(
        numer: i64,
        denom: i64,
    ) -> BigRational {
        BigRational::new(BigInt::from(numer), BigInt::from(denom))
    }

    /// Recombines a decomposition over the original denominator.
    fn recombine(
        parts: &Apart,
        denom: &[BigRational],
    ) -> QPoly {
        let mut total = univariate::mul(&parts.quotient, denom);
        for piece in &parts.pieces {
            let mut power: QPoly = vec![BigRational::one()];
            for _ in 0..piece.power {
                power = univariate::mul(&power, &piece.factor);
            }
            let (cofactor, rest) = univariate::divrem(denom, &power).unwrap_or_default();
            assert_eq!(rest, []);
            total = univariate::add(&total, &univariate::mul(&cofactor, &piece.numerator));
        }
        total
    }

    /// `recombine - numer` after cross-multiplication; empty exactly when the
    /// decomposition is the original rational function.
    fn recombination_defect(
        parts: &Apart,
        numer: &[BigRational],
        denom: &[BigRational],
    ) -> QPoly {
        univariate::sub(&recombine(parts, denom), numer)
    }

    /// Taylor coefficients of `p(r + t)` up to `t^(n - 1)`, by repeated
    /// division by `x - r`.
    fn taylor(
        p: &[BigRational],
        r: &BigRational,
        n: usize,
    ) -> Vec<BigRational> {
        let divisor = vec![-r.clone(), BigRational::one()];
        let mut rest: QPoly = p.to_vec();
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            match univariate::divrem(&rest, &divisor) {
                | Some((quotient, remainder)) => {
                    out.push(remainder.first().cloned().unwrap_or_else(BigRational::zero));
                    rest = quotient;
                },
                | None => out.push(BigRational::zero()),
            }
        }
        out
    }

    /// The first `n` coefficients of `a(t) / b(t)` at `t = 0`, `b(0) != 0`.
    fn series(
        a: &[BigRational],
        b: &[BigRational],
        n: usize,
    ) -> Vec<BigRational> {
        let lead = b.first().cloned().unwrap_or_else(BigRational::one);
        let mut out: Vec<BigRational> = Vec::with_capacity(n);
        for k in 0..n {
            let mut acc = a.get(k).cloned().unwrap_or_else(BigRational::zero);
            for i in 1..=k {
                let bi = b.get(i).cloned().unwrap_or_else(BigRational::zero);
                let si = out.get(k - i).cloned().unwrap_or_else(BigRational::zero);
                acc -= bi * si;
            }
            out.push(acc / &lead);
        }
        out
    }

    /// Multiplicity of `root` in `denom`, and `denom / (x - root)^m`.
    fn local(
        denom: &[BigRational],
        root: &BigRational,
    ) -> (u32, QPoly) {
        let divisor = vec![-root.clone(), BigRational::one()];
        let mut rest: QPoly = denom.to_vec();
        let mut power = 0;
        while rest.len() >= 2 {
            match univariate::divrem(&rest, &divisor) {
                | Some((quotient, remainder)) if remainder.is_empty() => {
                    rest = quotient;
                    power += 1;
                },
                | _ => break,
            }
        }
        (power, rest)
    }

    type Ext = (BigRational, BigRational);

    fn ext_mul(
        a: &Ext,
        b: &Ext,
        d: &BigRational,
    ) -> Ext {
        (&a.0 * &b.0 + &a.1 * &b.1 * d, &a.0 * &b.1 + &a.1 * &b.0)
    }

    fn ext_add(
        a: &Ext,
        b: &Ext,
    ) -> Ext {
        (&a.0 + &b.0, &a.1 + &b.1)
    }

    /// Asserts `c0/(x - r0) + c1/(x - r1) == numer / factor` exactly, for
    /// `factor = a (x - r0)(x - r1)`, by cross-multiplying in `Q(sqrt d)`:
    /// the sum of the conjugate pair is rational, and its cleared numerator
    /// must be the original numerator with the leading coefficient removed.
    fn quadratic_pair(
        split: &[RadicalPiece],
        numer: &[BigRational],
        factor: &[BigRational],
    ) {
        assert_eq!(split.len(), 2, "a quadratic factor splits into two terms");
        let d = split[0].radicand.clone();
        let a = factor.last().cloned().unwrap_or_else(BigRational::one);
        let coef: Vec<Ext> = split
            .iter()
            .map(|p| (p.real.clone(), p.radical.clone()))
            .collect();
        let pole: Vec<Ext> = split
            .iter()
            .map(|p| (p.offset.clone(), p.scale.clone()))
            .collect();
        // N = c0 (x - r1) + c1 (x - r0), D = (x - r0)(x - r1).
        let n_x = ext_add(&coef[0], &coef[1]);
        let n_0 = {
            let sum = ext_add(
                &ext_mul(&coef[0], &pole[1], &d),
                &ext_mul(&coef[1], &pole[0], &d),
            );
            (-sum.0, -sum.1)
        };
        let d_x = (-(&pole[0].0 + &pole[1].0), -(&pole[0].1 + &pole[1].1));
        let d_0 = ext_mul(&pole[0], &pole[1], &d);
        let zero = BigRational::zero();
        assert_eq!(n_x.1, zero, "numerator is not rational");
        assert_eq!(n_0.1, zero, "numerator is not rational");
        assert_eq!(d_x.1, zero, "denominator is not rational");
        assert_eq!(d_0.1, zero, "denominator is not rational");
        assert_eq!(
            &a * &n_x.0,
            numer.get(1).cloned().unwrap_or_else(BigRational::zero)
        );
        assert_eq!(
            &a * &n_0.0,
            numer.first().cloned().unwrap_or_else(BigRational::zero)
        );
        assert_eq!(
            d_x.0,
            &factor.get(1).cloned().unwrap_or_else(BigRational::zero) / &a
        );
        assert_eq!(
            d_0.0,
            &factor.first().cloned().unwrap_or_else(BigRational::zero) / &a
        );
    }

    #[test]
    fn decompositions_recombine() {
        let cases: [(&[i64], &[i64]); 5] = [
            (&[1], &[-1, 0, 1]),            // 1/(x^2 - 1)
            (&[0, 0, 0, 1], &[-1, 0, 1]),   // x^3/(x^2 - 1)
            (&[1, 2], &[0, 0, 1, 1]),       // (2x + 1)/(x^2 (x + 1))
            (&[3, 0, 1], &[1, 0, 2, 0, 1]), // (x^2 + 3)/(x^2 + 1)^2
            (&[1], &[-1, 0, 0, 1]),         // 1/(x^3 - 1)
        ];
        for (n, d) in cases {
            let (n, d) = (q(n), q(d));
            let parts = apart(&n, &d).unwrap_or_else(|| panic!("no decomposition"));
            assert_eq!(recombine(&parts, &d), n, "{parts:?}");
            for piece in &parts.pieces {
                assert!(piece.numerator.len() < piece.factor.len());
            }
        }
    }

    /// The cross-multiplied numerator is identically zero for every case,
    /// including repeated factors, irreducible quadratics, an irreducible
    /// cubic, a cancelling pair and an improper fraction.
    #[test]
    fn recombination_is_exactly_zero() {
        let cases: [(&[i64], &[i64]); 9] = [
            (&[1], &[-1, 0, 1]),          // 1/(x^2 - 1)
            (&[0, 0, 0, 1], &[-1, 0, 1]), // x^3/(x^2 - 1)
            (&[1, 2], &[0, 0, 1, 1]),     // (2x + 1)/(x^2 (x + 1))
            (&[1], &[0, 0, 1, 1]),        // 1/(x^2 (x + 1))
            (&[1], &[1, 2, 1]),           // 1/(x + 1)^2
            (&[1], &[-1, 0, 0, 1]),       // 1/(x^3 - 1)
            (&[1], &[1, 0, 0, 1]),        // 1/(x^3 + 1)
            (&[2], &[4, 0, 0, 0, 1]),     // 2/(x^4 + 4)
            (&[1, 1], &[2, 3, 1]),        // (x + 1)/((x + 1)(x + 2))
        ];
        for (n, d) in cases {
            let (n, d) = (q(n), q(d));
            let parts = apart(&n, &d).unwrap_or_else(|| panic!("no decomposition for {n:?}/{d:?}"));
            let defect = recombination_defect(&parts, &n, &d);
            assert!(defect.is_empty(), "{n:?}/{d:?} leaves {defect:?}");
            for piece in &parts.pieces {
                assert_ne!(piece.numerator, []);
                assert!(piece.numerator.len() < piece.factor.len());
                assert!(piece.power >= 1);
            }
        }
    }

    /// A repeated irreducible quadratic keeps both powers, each with a
    /// constant numerator: `(x^2 + 3)/(x^2 + 1)^2 = 1/(x^2 + 1) + 2/(x^2 + 1)^2`.
    #[test]
    fn repeated_irreducible_quadratic() {
        let parts = apart(&q(&[3, 0, 1]), &q(&[1, 0, 2, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[1, 0, 1]),
                    power: 1,
                    numerator: q(&[1]),
                },
                Piece {
                    factor: q(&[1, 0, 1]),
                    power: 2,
                    numerator: q(&[2]),
                },
            ]
        );
    }

    #[test]
    fn rational_unit_and_non_monic_repeated_factors_keep_order_and_omit_zero_pieces() {
        let linear = q(&[1, 2]);
        let quadratic = q(&[1, 0, 1]);
        let denominator = univariate::mul(
            &[frac(-3, 7)],
            &univariate::mul(
                &polynomial_power(&linear, 3),
                &polynomial_power(&quadratic, 2),
            ),
        );
        let expected = Apart {
            quotient: q(&[-2, 3]),
            pieces: vec![
                Piece {
                    factor: linear.clone(),
                    power: 1,
                    numerator: vec![frac(1, 3)],
                },
                Piece {
                    factor: linear,
                    power: 3,
                    numerator: vec![frac(5, 2)],
                },
                Piece {
                    factor: quadratic.clone(),
                    power: 1,
                    numerator: q(&[-3, 1]),
                },
                Piece {
                    factor: quadratic,
                    power: 2,
                    numerator: q(&[4, -2]),
                },
            ],
        };
        let mut numerator = recombine(&expected, &denominator);
        numerator.push(BigRational::zero());
        let mut padded_denominator = denominator;
        padded_denominator.push(BigRational::zero());
        assert_eq!(apart(&numerator, &padded_denominator), Some(expected));
    }

    #[test]
    fn independent_factor_blocks_preserve_factor_and_power_order_across_pools() {
        let factors: Vec<_> = (1..=12)
            .map(|k| (q(&[k, 1]), u32::try_from(k % 2 + 1).unwrap()))
            .collect();
        let denominator = factors
            .iter()
            .fold(q(&[1]), |product, (factor, multiplicity)| {
                univariate::mul(&product, &polynomial_power(factor, *multiplicity))
            });
        let expected = Apart {
            quotient: Vec::new(),
            pieces: factors
                .iter()
                .enumerate()
                .flat_map(|(index, (factor, multiplicity))| {
                    (1..=*multiplicity).map(move |power| {
                        Piece {
                            factor: factor.clone(),
                            power,
                            numerator: vec![frac(
                                i64::try_from(index).unwrap() + 1,
                                i64::from(power),
                            )],
                        }
                    })
                })
                .collect(),
        };
        let remainder = recombine(&expected, &denominator);
        for workers in [1, 2, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            assert_eq!(
                pool.install(|| decompose_proper(&remainder, &denominator, &factors)),
                Some(expected.pieces.clone())
            );
        }
    }

    #[test]
    fn invalid_factor_products_have_no_decomposition() {
        let denominator = q(&[1, 2, 1]);
        let remainder = q(&[1]);
        for factors in [
            vec![(q(&[1, 1]), 1)],
            vec![(q(&[1, 1]), 1), (q(&[1, 1]), 1)],
            vec![(q(&[0, 1]), 1), (q(&[1, 1]), 1)],
            vec![(q(&[1, 1]), 2), (q(&[2, 1]), 0)],
        ] {
            assert_eq!(decompose_proper(&remainder, &denominator, &factors), None);
        }
    }

    /// `1/(x^2 (x + 1)) = -1/x + 1/x^2 + 1/(x + 1)`, asserted as an exact
    /// term list.
    #[test]
    fn repeated_factor_terms_are_exact() {
        let parts = apart(&q(&[1]), &q(&[0, 0, 1, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(parts.quotient, []);
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[0, 1]),
                    power: 1,
                    numerator: q(&[-1]),
                },
                Piece {
                    factor: q(&[0, 1]),
                    power: 2,
                    numerator: q(&[1]),
                },
                Piece {
                    factor: q(&[1, 1]),
                    power: 1,
                    numerator: q(&[1]),
                },
            ]
        );
    }

    /// `1/(x^3 + 1) = 1/3/(x + 1) + (2/3 - x/3)/(x^2 - x + 1)`: a linear term
    /// plus a proper term over a genuinely irreducible quadratic.
    #[test]
    fn irreducible_quadratic_terms_are_exact() {
        let parts = apart(&q(&[1]), &q(&[1, 0, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(parts.quotient, []);
        assert_eq!(
            parts.pieces,
            vec![
                Piece {
                    factor: q(&[1, 1]),
                    power: 1,
                    numerator: vec![frac(1, 3)],
                },
                Piece {
                    factor: q(&[1, -1, 1]),
                    power: 1,
                    numerator: vec![frac(2, 3), frac(-1, 3)],
                },
            ]
        );
        let quadratic = &parts.pieces[1].factor;
        assert!(
            univariate::rational_roots(quadratic).is_empty(),
            "x^2 - x + 1 must have no rational root"
        );
        assert_eq!(quadratic.len(), 3);
        assert_eq!(parts.pieces[1].numerator.len(), 2);
    }

    /// At every rational pole, the coefficient of `(x - r)^-k` equals the
    /// standard local limit: the Taylor coefficient of order `m - k` of
    /// `numer(r + t) / (denom / (x - r)^m)` evaluated at `t = 0`.
    #[test]
    fn pole_coefficients_are_the_local_limits() {
        let cases: [(&[i64], &[i64]); 4] = [
            (&[1], &[-1, 0, 1]),      // 1/(x^2 - 1)
            (&[1], &[0, 0, 1, 1]),    // 1/(x^2 (x + 1))
            (&[1, 2], &[0, 0, 1, 1]), // (2x + 1)/(x^2 (x + 1))
            (&[1], &[-1, 0, 0, 1]),   // 1/(x^3 - 1)
        ];
        for (n, d) in cases {
            let (numer, denom) = (q(n), q(d));
            let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
            let roots = univariate::rational_roots(&denom);
            assert!(!roots.is_empty(), "{d:?} has rational poles");
            for root in roots {
                let (power, cofactor) = local(&denom, &root);
                assert!(power >= 1);
                let order = power as usize;
                let expanded = series(
                    &taylor(&numer, &root, order),
                    &taylor(&cofactor, &root, order),
                    order,
                );
                for k in 1..=power {
                    let expected = expanded[order - k as usize].clone();
                    let got = parts
                        .pieces
                        .iter()
                        .find(|p| {
                            p.power == k
                                && p.factor.len() == 2
                                && univariate::eval(&p.factor, &root).is_zero()
                        })
                        .and_then(|p| p.numerator.first().cloned());
                    assert_eq!(got, Some(expected), "pole {root} order {k} of {n:?}/{d:?}");
                }
            }
        }
    }

    /// `1/(x^2 - 2) = (sqrt 2/4)/(x - sqrt 2) - (sqrt 2/4)/(x + sqrt 2)`,
    /// exactly, with the recombination cross-multiplied in `Q(sqrt 8)`.
    #[test]
    fn real_quadratic_splits_over_the_extension() {
        let numer = q(&[1]);
        let denom = q(&[-2, 0, 1]);
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            parts.pieces,
            vec![Piece {
                factor: denom.clone(),
                power: 1,
                numerator: q(&[1]),
            }]
        );
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].radicand, frac(8, 1));
        assert_eq!(split[0].real, BigRational::zero());
        assert_eq!(split[0].radical, frac(1, 8));
        assert_eq!(split[0].offset, BigRational::zero());
        assert_eq!(split[0].scale, frac(1, 2));
        assert_eq!(split[1].radical, frac(-1, 8));
        assert_eq!(split[1].scale, frac(-1, 2));
        quadratic_pair(&split, &numer, &denom);
    }

    /// The same for a quadratic with a non-zero linear term.
    #[test]
    fn shifted_real_quadratic_splits() {
        let numer = q(&[1]);
        let denom = q(&[-1, 1, 1]); // x^2 + x - 1
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].radicand, frac(5, 1));
        assert_eq!(split[0].radical, frac(1, 5));
        assert_eq!(split[0].offset, frac(-1, 2));
        assert_eq!(split[0].scale, frac(1, 2));
        assert_eq!(split[1].radical, frac(-1, 5));
        assert_eq!(split[1].scale, frac(-1, 2));
        quadratic_pair(&split, &numer, &denom);
        assert_eq!(univariate::rational_roots(&denom), []);
    }

    /// A rational factor keeps its exact residue when a real quadratic sits
    /// beside it: `(x + 1)/((x + 2)(x^2 - 3))` has residue `-1` at `x = -2`.
    #[test]
    fn linear_and_quadratic_factors_mix() {
        let numer = q(&[1, 1]);
        let denom = univariate::mul(&q(&[2, 1]), &q(&[-3, 0, 1]));
        let parts = apart(&numer, &denom).unwrap_or_else(|| panic!("none"));
        let split = split_real_quadratics(&parts).unwrap_or_else(|| panic!("no split"));
        assert_eq!(split.len(), 3);
        assert_eq!(split[0].offset, frac(-2, 1));
        assert_eq!(split[0].scale, BigRational::zero());
        assert_eq!(split[0].radicand, BigRational::zero());
        assert_eq!(split[0].real, frac(-1, 1));
        assert_eq!(split[0].power, 1);
        assert_eq!(split[1].radicand, frac(12, 1));
        assert_eq!(split[1].scale, frac(1, 2));
        assert_eq!(split[2].scale, frac(-1, 2));
        quadratic_pair(&split[1..], &q(&[-1, 1]), &q(&[-3, 0, 1]));
    }

    /// Refusal is explicit: a cubic that is irreducible over `Q`, a repeated
    /// real quadratic, and a zero denominator are all reported as `None`
    /// rather than answered with a wrong structure.
    #[test]
    fn refusal_is_honest() {
        assert_eq!(apart(&q(&[1]), &[]), None);
        assert_eq!(apart(&q(&[1]), &q(&[0])), None);
        // x^3 - 2 has no rational root, so it stays one cubic factor.
        let cubic = apart(&q(&[1]), &q(&[-2, 0, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(cubic.pieces.len(), 1);
        assert_eq!(cubic.pieces[0].factor, q(&[-2, 0, 0, 1]));
        assert_eq!(split_real_quadratics(&cubic), None);
        // (x^2 - 2)^2 is a double pole over Q(sqrt 2), which is not attempted.
        let repeated = apart(&q(&[1]), &q(&[4, 0, -4, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(
            repeated.pieces,
            vec![Piece {
                factor: q(&[-2, 0, 1]),
                power: 2,
                numerator: q(&[1]),
            }]
        );
        assert_eq!(split_real_quadratics(&repeated), None);
        // A quadratic with complex roots is already irreducible over R.
        let complex = apart(&q(&[1]), &q(&[1, 0, 1])).unwrap_or_else(|| panic!("none"));
        assert_eq!(split_real_quadratics(&complex), None);
    }
}
