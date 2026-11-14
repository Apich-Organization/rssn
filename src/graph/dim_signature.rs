//! # Dimension signatures for plugin operators
//!
//! A built-in operator states how its arguments' dimensions determine its
//! value's with one of the fixed [`DimRule`](crate::graph::DimRule)s
//! (`Same(0)`, `Product`, …). An operator a plugin defines from a typed
//! source function knows something more precise: its declared parameter and
//! result types, such as
//!
//! ```text
//! area     : (Length^?d)           -> Length^(2 ?d)
//! kinetic  : (Mass, Length/Time)   -> Mass Length^2 / Time^2
//! softplus : (1)                   -> 1
//! ```
//!
//! A [`DimSignature`] is that statement: every parameter and the result is a
//! *monomial* `C · Π v_k^(e_k)` of a constant dimension `C` and signature
//! variables `v_k` with rational exponents. [`DimSignature::apply`] solves
//! the variables from the argument dimensions (a parameter that holds one
//! unsolved variable determines it: `v = (arg / C)^(1/e)`), checks every
//! other parameter against its solved form, and returns the result monomial.
//! An argument whose dimension is unknown (`Any`, or still a unification
//! variable of the inference) neither determines nor contradicts anything;
//! a result that depends on an unsolved variable is `Any`.
//!
//! The signature is attached to the operator as an attribute and read by the
//! generic [`signature_rule`], which is what
//! [`Registry::register_dim_signature`](crate::kernel::plugin::Registry::register_dim_signature)
//! installs as the operator's [`DimRule::Custom`](crate::graph::DimRule::Custom)
//! rule. One rule function serves every such operator, so no per-operator
//! closure has to live in a `fn` pointer.

use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;

use crate::dim::DimExpr;
use crate::dim::Factor;
use crate::graph::DimError;
use crate::graph::DimTerm;
use crate::graph::dim::DimArena;
use crate::graph::dim::DimCtx;
use crate::kernel::signature::Exp;

/// `C · Π v_k^(e_k)`.
#[derive(Clone, Debug)]
struct Monomial {
    constant: DimTerm,
    vars: Vec<(usize, BigRational)>,
}

/// The dimensions of an operator's parameters and result; see the
/// [module documentation](self).
#[derive(Clone, Debug)]
pub struct DimSignature {
    params: Vec<Monomial>,
    result: Monomial,
    /// Names of the signature variables, for messages.
    names: Vec<String>,
}

fn rational(exponent: &Exp) -> Result<BigRational, DimError> {
    match exponent {
        | Exp::Rat(r) => {
            if r.denom == 0 {
                return Err(DimError::Mismatch {
                    left: format!("{}/{}", r.numer, r.denom),
                    right: "a rational exponent".to_owned(),
                    context: "dimension signature".to_owned(),
                });
            }
            Ok(BigRational::new(r.numer.into(), r.denom.into()))
        },
        | Exp::Symbol(name) => {
            Err(DimError::SymbolicConstraint {
                constraint: format!("the exponent `{name}` is symbolic"),
                context: "a dimension signature needs rational exponents".to_owned(),
            })
        },
    }
}

impl DimSignature {
    /// Builds a signature; base names are interned into `arena`.
    ///
    /// # Errors
    /// A parameter or result that is not a product of factors with rational
    /// exponents ([`DimExpr::Any`] and [`DimExpr::NonMetric`] are not).
    pub fn new(
        params: &[DimExpr],
        result: &DimExpr,
        arena: &mut DimArena,
    ) -> Result<Self, DimError> {
        let mut names: Vec<String> = Vec::new();
        let mut monomial = |expr: &DimExpr| -> Result<Monomial, DimError> {
            let mut constant = DimTerm::dimensionless();
            let mut vars: Vec<(usize, BigRational)> = Vec::new();
            for (factor, exponent) in expr.factors()? {
                let power = rational(&exponent)?;
                match factor {
                    | Factor::Base(name) => {
                        let base = arena.define_base(&name);
                        constant = constant.mul(&DimTerm::base(base, power));
                    },
                    | Factor::Var(name) => {
                        let index = match names.iter().position(|known| *known == name) {
                            | Some(index) => index,
                            | None => {
                                names.push(name);
                                names.len().saturating_sub(1)
                            },
                        };
                        match vars.iter_mut().find(|(v, _)| *v == index) {
                            | Some((_, e)) => *e += power,
                            | None => vars.push((index, power)),
                        }
                    },
                }
            }
            vars.retain(|(_, e)| !e.is_zero());
            Ok(Monomial { constant, vars })
        };
        let params = params
            .iter()
            .map(&mut monomial)
            .collect::<Result<Vec<_>, _>>()?;
        let result = monomial(result)?;
        Ok(Self {
            params,
            result,
            names,
        })
    }

    /// Number of parameters.
    #[must_use]
    pub const fn arity(&self) -> usize {
        self.params.len()
    }

    /// The value's dimension for arguments of dimensions `args`.
    ///
    /// # Errors
    /// [`DimError::Mismatch`] when an argument contradicts the signature
    /// (`expected` and `found` are rendered with `arena` when given).
    pub fn apply(
        &self,
        args: &[DimTerm],
        arena: Option<&DimArena>,
    ) -> Result<DimTerm, DimError> {
        let render = |term: &DimTerm| {
            arena.map_or_else(|| format!("{term:?}"), |arena| term.to_string(arena))
        };
        let mut solved: Vec<Option<DimTerm>> = vec![None; self.names.len()];
        let known = |term: &DimTerm| !term.is_any() && !term.has_unknown();
        // Each pass solves every parameter with exactly one open variable;
        // a parameter with several waits for the others.
        let mut open: Vec<usize> = (0..self.params.len().min(args.len())).collect();
        loop {
            let mut progress = false;
            let mut still = Vec::new();
            for index in open {
                let (Some(param), Some(arg)) = (self.params.get(index), args.get(index)) else {
                    continue;
                };
                if !known(arg) {
                    continue;
                }
                let unsolved: Vec<&(usize, BigRational)> = param
                    .vars
                    .iter()
                    .filter(|(v, _)| solved.get(*v).is_some_and(Option::is_none))
                    .collect();
                match unsolved.as_slice() {
                    | [] => {
                        let expected = evaluate(param, &solved);
                        if let Some(expected) = expected
                            && expected != *arg
                        {
                            return Err(DimError::Mismatch {
                                left: render(&expected),
                                right: render(arg),
                                context: format!("argument {}", index.saturating_add(1)),
                            });
                        }
                    },
                    | [(var, exponent)] => {
                        // arg = C · Π(solved) · v^e  ⇒  v = (arg / (C · Π(solved)))^(1/e)
                        let mut rest = param.constant.clone();
                        for (other, e) in &param.vars {
                            if other == var {
                                continue;
                            }
                            if let Some(Some(value)) = solved.get(*other) {
                                rest = rest.mul(&value.pow(e));
                            }
                        }
                        let value = arg.mul(&rest.inverse()).pow(&exponent.recip());
                        if let Some(slot) = solved.get_mut(*var) {
                            *slot = Some(value);
                        }
                        progress = true;
                    },
                    | _ => still.push(index),
                }
            }
            open = still;
            if !progress || open.is_empty() {
                break;
            }
        }
        Ok(evaluate(&self.result, &solved).unwrap_or_else(DimTerm::any))
    }
}

/// `C · Π v^e` with the solved variables, or `None` when one is open.
fn evaluate(
    monomial: &Monomial,
    solved: &[Option<DimTerm>],
) -> Option<DimTerm> {
    let mut value = monomial.constant.clone();
    for (var, exponent) in &monomial.vars {
        let solution = solved.get(*var)?.as_ref()?;
        if exponent.is_one() {
            value = value.mul(solution);
        } else {
            value = value.mul(&solution.pow(exponent));
        }
    }
    Some(value)
}

/// The [`DimRule::Custom`](crate::graph::DimRule::Custom) rule of every
/// operator with a [`DimSignature`] attribute.
///
/// An operator without the attribute (a misconfiguration) is `Any`, the
/// same answer as for a polymorphic operator, rather than an error.
///
/// # Errors
/// What [`DimSignature::apply`] reports.
pub fn signature_rule(
    ctx: &DimCtx<'_>,
    args: &[DimTerm],
) -> Result<DimTerm, DimError> {
    let op = ctx.graph.op(ctx.node);
    match ctx.graph.ops().attr::<DimSignature>(op) {
        | Some(signature) => signature.apply(args, Some(ctx.arena)),
        | None => Ok(DimTerm::any()),
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::kernel::signature::Rational;

    fn length() -> DimExpr {
        DimExpr::base("Length")
    }

    fn time() -> DimExpr {
        DimExpr::base("Time")
    }

    #[test]
    fn constant_signatures_check_and_return() {
        let mut arena = DimArena::default();
        let sig = DimSignature::new(
            &[length(), time()],
            &DimExpr::Mul(vec![length(), DimExpr::Inverse(Box::new(time()))]),
            &mut arena,
        )
        .unwrap();
        let l = length().to_term(&mut arena).unwrap();
        let t = time().to_term(&mut arena).unwrap();
        let v = sig.apply(&[l.clone(), t.clone()], Some(&arena)).unwrap();
        assert_eq!(v, l.mul(&t.inverse()));
        assert!(sig.apply(&[t.clone(), t], Some(&arena)).is_err());
    }

    #[test]
    fn a_variable_is_solved_from_its_argument() {
        let mut arena = DimArena::default();
        let d = DimExpr::Var("d".to_owned());
        // (d) -> d^2
        let sig =
            DimSignature::new(std::slice::from_ref(&d), &d.clone().pow(2_i64), &mut arena).unwrap();
        let l = length().to_term(&mut arena).unwrap();
        let area = sig.apply(std::slice::from_ref(&l), Some(&arena)).unwrap();
        assert_eq!(area, l.pow(&BigRational::from_integer(2.into())));
        // (d^(1/2), d) -> d: the first argument determines d, the second must agree.
        let half = DimExpr::Pow(
            Box::new(d.clone()),
            Exp::Rat(Rational { numer: 1, denom: 2 }),
        );
        let sig = DimSignature::new(&[half, d.clone()], &d, &mut arena).unwrap();
        let area_term = l.pow(&BigRational::from_integer(2.into()));
        let ok = sig
            .apply(&[l.clone(), area_term.clone()], Some(&arena))
            .unwrap();
        assert_eq!(ok, area_term);
        assert!(sig.apply(&[l.clone(), l], Some(&arena)).is_err());
    }

    #[test]
    fn unknown_arguments_leave_the_result_open() {
        let mut arena = DimArena::default();
        let d = DimExpr::Var("d".to_owned());
        let sig = DimSignature::new(std::slice::from_ref(&d), &d, &mut arena).unwrap();
        assert!(sig.apply(&[DimTerm::any()], Some(&arena)).unwrap().is_any());
        let fixed = DimSignature::new(&[d], &time(), &mut arena).unwrap();
        let t = time().to_term(&mut arena).unwrap();
        assert_eq!(fixed.apply(&[DimTerm::any()], Some(&arena)).unwrap(), t);
    }

    #[test]
    fn symbolic_exponents_are_refused() {
        let mut arena = DimArena::default();
        let n = DimExpr::Pow(Box::new(length()), Exp::Symbol("n".to_owned()));
        assert!(DimSignature::new(&[n], &length(), &mut arena).is_err());
    }
}
