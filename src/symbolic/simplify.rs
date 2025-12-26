//! # Symbolic Expression Simplification (Legacy)
//!
//! **⚠️ DEPRECATION NOTICE**: This module is deprecated in favor of
//! [`simplify_dag`](crate::symbolic::simplify_dag), which provides better performance
//! and more comprehensive simplification through DAG-based algorithms. This module is
//! maintained for backward compatibility and will continue to receive bug fixes.
//!
//! ## Overview
//!
//! This module provides AST-based symbolic expression simplification. It includes:
//!
//! - **`simplify`**: Core simplification function applying deterministic algebraic rules
//! - **`heuristic_simplify`**: Pattern-based simplification using rewrite rules
//! - **Utility functions**: Term collection, rational expression simplification
//!
//! ## Simplification Strategy
//!
//! The simplification process works in multiple phases:
//!
//! 1. **Recursive simplification**: Bottom-up traversal simplifying children first
//! 2. **Rule application**: Apply algebraic identities and arithmetic evaluations
//! 3. **Term collection**: Collect like terms in sums and products
//! 4. **Rational simplification**: Cancel common factors in fractions
//!
//! ## Supported Operations
//!
//! The simplifier handles:
//! - **Arithmetic**: Addition, subtraction, multiplication, division, power
//! - **Trigonometric**: sin, cos, tan with basic identities
//! - **Exponential/Logarithmic**: exp, log with basic rules
//! - **Algebraic**: Polynomial expansion, factoring, rational expressions
//! - **N-ary operations**: `AddList`, `MulList` (new)
//! - **Dynamic operations**: `UnaryList`, `BinaryList`, `NaryList` (new)
//!
//! ## Examples
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::simplify::simplify;
//!
//! // Basic arithmetic simplification
//! let expr = Expr::new_add(
//!     Expr::new_constant(2.0),
//!     Expr::new_constant(3.0),
//! );
//!
//! let result = simplify(expr);
//! // result is Expr::Constant(5.0)
//! ```
//!
//! ## Migration Guide
//!
//! For new code, prefer [`simplify_dag::simplify`](crate::symbolic::simplify_dag::simplify):
//!
//! ```rust
//! 
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::simplify_dag::simplify;
//!
//! let expr = Expr::new_add(
//!     Expr::new_variable("x"),
//!     Expr::new_constant(0.0),
//! );
//!
//! let result = simplify(&expr); // DAG-based simplification
//! ```
//!
//! ## Performance Notes
//!
//! - This AST-based simplifier may duplicate work on shared subexpressions
//! - For complex expressions, consider using `simplify_dag` for better performance
//! - Caching is used internally to avoid redundant simplifications
//!
//! ## See Also
//!
//! - [`simplify_dag`](crate::symbolic::simplify_dag) - Modern DAG-based simplification (recommended)
//! - [`core`](crate::symbolic::core) - Core expression types and operations
//! - [`calculus`](crate::symbolic::calculus) - Symbolic calculus operations

#![allow(deprecated)]

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::sync::Arc;

use num_bigint::BigInt;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::symbolic::calculus::substitute;
use crate::symbolic::core::DagNode;
use crate::symbolic::core::DagOp;
use crate::symbolic::core::Expr;

pub(crate) fn simplify_dag_node(
    node: &Arc<DagNode>,
    cache: &mut HashMap<u64, Expr>,
) -> Expr {

    if let Some(simplified) =
        cache.get(&node.hash)
    {

        return simplified.clone();
    }

    let simplified_children = node
        .children
        .iter()
        .map(|child| {

            simplify_dag_node(
                child, cache,
            )
        })
        .collect::<Vec<Expr>>();

    let new_expr =
        build_expr_from_op_and_children(
            &node.op,
            simplified_children,
        );

    let simplified_expr =
        apply_rules(new_expr); // apply_rules from simplify.rs

    cache.insert(
        node.hash,
        simplified_expr.clone(),
    );

    simplified_expr
}

pub(crate) fn build_expr_from_op_and_children(
    op: &DagOp,
    children: Vec<Expr>,
) -> Expr {

    macro_rules! arc {
        ($idx:expr) => {

            Arc::new(
                children[$idx].clone(),
            )
        };
    }

    match op {
        | DagOp::Constant(c) => Expr::Constant(c.into_inner()),
        | DagOp::BigInt(i) => Expr::BigInt(i.clone()),
        | DagOp::Rational(r) => Expr::Rational(r.clone()),
        | DagOp::Boolean(b) => Expr::Boolean(*b),
        | DagOp::Variable(s) => Expr::Variable(s.clone()),
        | DagOp::Pattern(s) => Expr::Pattern(s.clone()),
        | DagOp::Domain(s) => Expr::Domain(s.clone()),
        | DagOp::Pi => Expr::Pi,
        | DagOp::E => Expr::E,
        | DagOp::Infinity => Expr::Infinity,
        | DagOp::NegativeInfinity => Expr::NegativeInfinity,
        | DagOp::InfiniteSolutions => Expr::InfiniteSolutions,
        | DagOp::NoSolution => Expr::NoSolution,
        | DagOp::Derivative(s) => Expr::Derivative(arc!(0), s.clone()),
        | DagOp::DerivativeN(s) => {
            Expr::DerivativeN(
                arc!(0),
                s.clone(),
                arc!(1),
            )
        },
        | DagOp::Limit(s) => {
            Expr::Limit(
                arc!(0),
                s.clone(),
                arc!(1),
            )
        },
        | DagOp::Solve(s) => Expr::Solve(arc!(0), s.clone()),
        | DagOp::ConvergenceAnalysis(s) => Expr::ConvergenceAnalysis(arc!(0), s.clone()),
        | DagOp::ForAll(s) => Expr::ForAll(s.clone(), arc!(0)),
        | DagOp::Exists(s) => Expr::Exists(s.clone(), arc!(0)),
        | DagOp::Substitute(s) => {
            Expr::Substitute(
                arc!(0),
                s.clone(),
                arc!(1),
            )
        },
        | DagOp::Ode {
            func,
            var,
        } => {
            Expr::Ode {
                equation : arc!(0),
                func : func.clone(),
                var : var.clone(),
            }
        },
        | DagOp::Pde {
            func,
            vars,
        } => {
            Expr::Pde {
                equation : arc!(0),
                func : func.clone(),
                vars : vars.clone(),
            }
        },
        | DagOp::Predicate {
            name,
        } => {
            Expr::Predicate {
                name : name.clone(),
                args : children,
            }
        },
        | DagOp::Path(pt) => {
            Expr::Path(
                pt.clone(),
                arc!(0),
                arc!(1),
            )
        },
        | DagOp::Interval(incl_lower, incl_upper) => {
            Expr::Interval(
                arc!(0),
                arc!(1),
                *incl_lower,
                *incl_upper,
            )
        },
        | DagOp::RootOf {
            index,
        } => {
            Expr::RootOf {
                poly : arc!(0),
                index : *index,
            }
        },
        | DagOp::SparsePolynomial(p) => Expr::SparsePolynomial(p.clone()),
        | DagOp::QuantityWithValue(u) => Expr::QuantityWithValue(arc!(0), u.clone()),
        | DagOp::Add => Expr::Add(arc!(0), arc!(1)),
        | DagOp::Sub => Expr::Sub(arc!(0), arc!(1)),
        | DagOp::Mul => Expr::Mul(arc!(0), arc!(1)),
        | DagOp::Div => Expr::Div(arc!(0), arc!(1)),
        | DagOp::Neg => Expr::Neg(arc!(0)),
        | DagOp::Power => Expr::Power(arc!(0), arc!(1)),
        | DagOp::Sin => Expr::Sin(arc!(0)),
        | DagOp::Cos => Expr::Cos(arc!(0)),
        | DagOp::Tan => Expr::Tan(arc!(0)),
        | DagOp::Exp => Expr::Exp(arc!(0)),
        | DagOp::Log => Expr::Log(arc!(0)),
        | DagOp::Abs => Expr::Abs(arc!(0)),
        | DagOp::Sqrt => Expr::Sqrt(arc!(0)),
        | DagOp::Eq => Expr::Eq(arc!(0), arc!(1)),
        | DagOp::Lt => Expr::Lt(arc!(0), arc!(1)),
        | DagOp::Gt => Expr::Gt(arc!(0), arc!(1)),
        | DagOp::Le => Expr::Le(arc!(0), arc!(1)),
        | DagOp::Ge => Expr::Ge(arc!(0), arc!(1)),
        | DagOp::Matrix {
            rows: _,
            cols,
        } => {

            let reconstructed_matrix : Vec<Vec<Expr>> = children
                .chunks(*cols)
                .map(<[Expr]>::to_vec)
                .collect();

            Expr::Matrix(reconstructed_matrix)
        },
        | DagOp::Vector => Expr::Vector(children),
        | DagOp::Complex => Expr::Complex(arc!(0), arc!(1)),
        | DagOp::Transpose => Expr::Transpose(arc!(0)),
        | DagOp::MatrixMul => Expr::MatrixMul(arc!(0), arc!(1)),
        | DagOp::MatrixVecMul => Expr::MatrixVecMul(arc!(0), arc!(1)),
        | DagOp::Inverse => Expr::Inverse(arc!(0)),
        | DagOp::Integral => {
            Expr::Integral {
                integrand : arc!(0),
                var : arc!(1),
                lower_bound : arc!(2),
                upper_bound : arc!(3),
            }
        },
        | DagOp::VolumeIntegral => {
            Expr::VolumeIntegral {
                scalar_field : arc!(0),
                volume : arc!(1),
            }
        },
        | DagOp::SurfaceIntegral => {
            Expr::SurfaceIntegral {
                vector_field : arc!(0),
                surface : arc!(1),
            }
        },
        | DagOp::Sum => {
            Expr::Sum {
                body : arc!(0),
                var : arc!(1),
                from : arc!(2),
                to : arc!(3),
            }
        },
        | DagOp::Series(s) => {
            Expr::Series(
                arc!(0),
                s.clone(),
                arc!(1),
                arc!(2),
            )
        },
        | DagOp::Summation(s) => {
            Expr::Summation(
                arc!(0),
                s.clone(),
                arc!(1),
                arc!(2),
            )
        },
        | DagOp::Product(s) => {
            Expr::Product(
                arc!(0),
                s.clone(),
                arc!(1),
                arc!(2),
            )
        },
        | DagOp::AsymptoticExpansion(s) => {
            Expr::AsymptoticExpansion(
                arc!(0),
                s.clone(),
                arc!(1),
                arc!(2),
            )
        },
        | DagOp::Sec => Expr::Sec(arc!(0)),
        | DagOp::Csc => Expr::Csc(arc!(0)),
        | DagOp::Cot => Expr::Cot(arc!(0)),
        | DagOp::ArcSin => Expr::ArcSin(arc!(0)),
        | DagOp::ArcCos => Expr::ArcCos(arc!(0)),
        | DagOp::ArcTan => Expr::ArcTan(arc!(0)),
        | DagOp::ArcSec => Expr::ArcSec(arc!(0)),
        | DagOp::ArcCsc => Expr::ArcCsc(arc!(0)),
        | DagOp::ArcCot => Expr::ArcCot(arc!(0)),
        | DagOp::Sinh => Expr::Sinh(arc!(0)),
        | DagOp::Cosh => Expr::Cosh(arc!(0)),
        | DagOp::Tanh => Expr::Tanh(arc!(0)),
        | DagOp::Sech => Expr::Sech(arc!(0)),
        | DagOp::Csch => Expr::Csch(arc!(0)),
        | DagOp::Coth => Expr::Coth(arc!(0)),
        | DagOp::ArcSinh => Expr::ArcSinh(arc!(0)),
        | DagOp::ArcCosh => Expr::ArcCosh(arc!(0)),
        | DagOp::ArcTanh => Expr::ArcTanh(arc!(0)),
        | DagOp::ArcSech => Expr::ArcSech(arc!(0)),
        | DagOp::ArcCsch => Expr::ArcCsch(arc!(0)),
        | DagOp::ArcCoth => Expr::ArcCoth(arc!(0)),
        | DagOp::LogBase => Expr::LogBase(arc!(0), arc!(1)),
        | DagOp::Atan2 => Expr::Atan2(arc!(0), arc!(1)),
        | DagOp::Binomial => Expr::Binomial(arc!(0), arc!(1)),
        | DagOp::Factorial => Expr::Factorial(arc!(0)),
        | DagOp::Permutation => Expr::Permutation(arc!(0), arc!(1)),
        | DagOp::Combination => Expr::Combination(arc!(0), arc!(1)),
        | DagOp::FallingFactorial => Expr::FallingFactorial(arc!(0), arc!(1)),
        | DagOp::RisingFactorial => Expr::RisingFactorial(arc!(0), arc!(1)),
        | DagOp::Boundary => Expr::Boundary(arc!(0)),
        | DagOp::Gamma => Expr::Gamma(arc!(0)),
        | DagOp::Beta => Expr::Beta(arc!(0), arc!(1)),
        | DagOp::Erf => Expr::Erf(arc!(0)),
        | DagOp::Erfc => Expr::Erfc(arc!(0)),
        | DagOp::Erfi => Expr::Erfi(arc!(0)),
        | DagOp::Zeta => Expr::Zeta(arc!(0)),
        | DagOp::BesselJ => Expr::BesselJ(arc!(0), arc!(1)),
        | DagOp::BesselY => Expr::BesselY(arc!(0), arc!(1)),
        | DagOp::LegendreP => Expr::LegendreP(arc!(0), arc!(1)),
        | DagOp::LaguerreL => Expr::LaguerreL(arc!(0), arc!(1)),
        | DagOp::HermiteH => Expr::HermiteH(arc!(0), arc!(1)),
        | DagOp::Digamma => Expr::Digamma(arc!(0)),
        | DagOp::KroneckerDelta => Expr::KroneckerDelta(arc!(0), arc!(1)),
        | DagOp::And => Expr::And(children),
        | DagOp::Or => Expr::Or(children),
        | DagOp::Not => Expr::Not(arc!(0)),
        | DagOp::Xor => Expr::Xor(arc!(0), arc!(1)),
        | DagOp::Implies => Expr::Implies(arc!(0), arc!(1)),
        | DagOp::Equivalent => Expr::Equivalent(arc!(0), arc!(1)),
        | DagOp::Union => Expr::Union(children),
        | DagOp::Polynomial => Expr::Polynomial(children),
        | DagOp::Floor => Expr::Floor(arc!(0)),
        | DagOp::IsPrime => Expr::IsPrime(arc!(0)),
        | DagOp::Gcd => Expr::Gcd(arc!(0), arc!(1)),
        | DagOp::Mod => Expr::Mod(arc!(0), arc!(1)),
        | DagOp::System => Expr::System(children),
        | DagOp::Solutions => Expr::Solutions(children),
        | DagOp::ParametricSolution => {
            Expr::ParametricSolution {
                x : arc!(0),
                y : arc!(1),
            }
        },
        | DagOp::GeneralSolution => Expr::GeneralSolution(arc!(0)),
        | DagOp::ParticularSolution => Expr::ParticularSolution(arc!(0)),
        | DagOp::Fredholm => {
            Expr::Fredholm(
                arc!(0),
                arc!(1),
                arc!(2),
                arc!(3),
            )
        },
        | DagOp::Volterra => {
            Expr::Volterra(
                arc!(0),
                arc!(1),
                arc!(2),
                arc!(3),
            )
        },
        | DagOp::Apply => Expr::Apply(arc!(0), arc!(1)),
        | DagOp::Tuple => Expr::Tuple(children),
        | DagOp::Distribution => {
            Expr::Distribution(
                children[0]
                    .clone_box_dist()
                    .expect("Dag Distribution"),
            )
        },
        | DagOp::Max => Expr::Max(arc!(0), arc!(1)),
        | DagOp::Quantity => {
            Expr::Quantity(
                children[0]
                    .clone_box_quant()
                    .expect("Dag Quatity"),
            )
        },
        | DagOp::UnaryList(s) => Expr::UnaryList(s.clone(), arc!(0)),
        | DagOp::BinaryList(s) => {
            Expr::BinaryList(
                s.clone(),
                arc!(0),
                arc!(1),
            )
        },
        | DagOp::NaryList(s) => Expr::NaryList(s.clone(), children),
        | _ => {
            Expr::CustomString(format!(
                "Unimplemented: {op:?}"
            ))
        },
    }
}

#[deprecated(
    since = "0.1.10",
    note = "Please use `simplify_dag` \
            instead."
)]
/// The main simplification function.
/// It recursively simplifies an expression tree by applying deterministic algebraic rules.
///
/// This function performs a deep simplification, traversing the expression tree
/// and applying various algebraic identities and arithmetic evaluations.
/// It also includes a step for simplifying rational expressions by canceling common factors.
///
/// # Arguments
/// * `expr` - The expression to simplify.
///
/// # Returns
/// A new, simplified `Expr`.
#[must_use]

pub fn simplify(expr: Expr) -> Expr {

    if let Expr::Dag(node) = expr {

        let mut cache = HashMap::new();

        simplify_dag_node(
            &node,
            &mut cache,
        )
    } else {

        let mut cache = HashMap::new();

        simplify_with_cache(
            &expr,
            &mut cache,
        )
    }
}

#[inline]
#[must_use]

pub fn is_zero(expr: &Expr) -> bool {

    match expr {
        | Expr::Dag(node) => {
            is_zero(
                &node
                    .to_expr()
                    .expect(
                        "Dag is Zero",
                    ),
            )
        },
        | Expr::Constant(val)
            if *val == 0.0 =>
        {
            true
        },
        | Expr::BigInt(val)
            if val.is_zero() =>
        {
            true
        },
        | Expr::Rational(val)
            if val.is_zero() =>
        {
            true
        },
        | _ => false,
    }
}

#[inline]
#[must_use]

pub fn is_one(expr: &Expr) -> bool {

    match expr {
        | Expr::Dag(node) => {
            is_one(
                &node
                    .to_expr()
                    .expect(
                        "Dag is One",
                    ),
            )
        },
        | Expr::Constant(val)
            if *val == 1.0 =>
        {
            true
        },
        | Expr::BigInt(val)
            if val.is_one() =>
        {
            true
        },
        | Expr::Rational(val)
            if val.is_one() =>
        {
            true
        },
        | _ => false,
    }
}

#[inline]
#[must_use]

pub fn as_f64(
    expr: &Expr
) -> Option<f64> {

    match expr {
        | Expr::Dag(node) => {
            as_f64(
                &node
                    .to_expr()
                    .expect(
                        "Dat is f64",
                    ),
            )
        },
        | Expr::Constant(val) => {
            Some(*val)
        },
        | Expr::BigInt(val) => {
            val.to_f64()
        },
        | Expr::Rational(val) => {
            val.to_f64()
        },
        | _ => None,
    }
}

/// The main simplification function with caching.
/// It recursively simplifies an expression tree by applying deterministic algebraic rules.

pub(crate) fn simplify_with_cache(
    expr: &Expr,
    cache: &mut HashMap<Expr, Expr>,
) -> Expr {

    if let Some(cached_result) =
        cache.get(expr)
    {

        return cached_result.clone();
    }

    let result = {

        let simplified_children_expr = match expr {
            | Expr::Add(a, b) => {
                Expr::new_add(
                    simplify_with_cache(a, cache),
                    simplify_with_cache(b, cache),
                )
            },
            | Expr::Sub(a, b) => {
                Expr::new_sub(
                    simplify_with_cache(a, cache),
                    simplify_with_cache(b, cache),
                )
            },
            | Expr::Mul(a, b) => {
                Expr::new_mul(
                    simplify_with_cache(a, cache),
                    simplify_with_cache(b, cache),
                )
            },
            | Expr::Div(a, b) => {
                Expr::new_div(
                    simplify_with_cache(a, cache),
                    simplify_with_cache(b, cache),
                )
            },
            | Expr::Power(b, e) => {
                Expr::new_pow(
                    simplify_with_cache(b, cache),
                    simplify_with_cache(e, cache),
                )
            },
            | Expr::Sin(arg) => {
                Expr::new_sin(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Cos(arg) => {
                Expr::new_cos(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Tan(arg) => {
                Expr::new_tan(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Exp(arg) => {
                Expr::new_exp(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Log(arg) => {
                Expr::new_log(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Neg(arg) => {
                Expr::new_neg(simplify_with_cache(
                    arg, cache,
                ))
            },
            | Expr::Sum {
                body,
                var,
                from,
                to,
            } => {
                Expr::Sum {
                    body : Arc::new(simplify_with_cache(
                        body, cache,
                    )),
                    var : Arc::new(simplify_with_cache(
                        var, cache,
                    )),
                    from : Arc::new(simplify_with_cache(
                        from, cache,
                    )),
                    to : Arc::new(simplify_with_cache(
                        to, cache,
                    )),
                }
            },
            // N-ary list variants
            | Expr::AddList(terms) => {

                let simplified_terms : Vec<Expr> = terms
                    .iter()
                    .map(|t| simplify_with_cache(t, cache))
                    .collect();

                Expr::AddList(simplified_terms)
            },
            | Expr::MulList(factors) => {

                let simplified_factors : Vec<Expr> = factors
                    .iter()
                    .map(|f| simplify_with_cache(f, cache))
                    .collect();

                Expr::MulList(simplified_factors)
            },
            // Generic list variants
            | Expr::UnaryList(name, arg) => {
                Expr::UnaryList(
                    name.clone(),
                    Arc::new(simplify_with_cache(
                        arg, cache,
                    )),
                )
            },
            | Expr::BinaryList(name, a, b) => {
                Expr::BinaryList(
                    name.clone(),
                    Arc::new(simplify_with_cache(
                        a, cache,
                    )),
                    Arc::new(simplify_with_cache(
                        b, cache,
                    )),
                )
            },
            | Expr::NaryList(name, args) => {

                let simplified_args : Vec<Expr> = args
                    .iter()
                    .map(|arg| simplify_with_cache(arg, cache))
                    .collect();

                Expr::NaryList(
                    name.clone(),
                    simplified_args,
                )
            },
            | _ => expr.clone(),
        };

        let simplified_expr = apply_rules(simplified_children_expr);

        simplify_rational_expression(
            &simplified_expr,
        )
    };

    cache.insert(
        expr.clone(),
        result.clone(),
    );

    result
}

/// Applies a set of deterministic simplification rules to an expression.
#[allow(clippy::unnecessary_to_owned)]

pub(crate) fn apply_rules(
    expr: Expr
) -> Expr {

    match expr {
        | Expr::Add(a, b) => {
            match simplify_add(
                (*a).clone(),
                (*b).clone(),
            ) {
                | Ok(value) => value,
                | Err(value) => value,
            }
        },
        | Expr::Sub(a, b) => {

            if let Some(value) =
                simplify_sub(&a, &b)
            {

                return value;
            }

            Expr::new_sub(a, b)
        },
        | Expr::Mul(a, b) => {

            if let Some(value) =
                simplify_mul(&a, &b)
            {

                return value;
            }

            Expr::new_mul(a, b)
        },
        | Expr::Div(a, b) => {

            if let Some(value) =
                simplify_div(&a, &b)
            {

                return value;
            }

            Expr::new_div(a, b)
        },
        | Expr::Power(b, e) => {

            if let Some(value) =
                simplify_power(&b, &e)
            {

                return value;
            }

            Expr::new_pow(b, e)
        },
        | Expr::Sqrt(arg) => {
            simplify_sqrt(
                (*arg).clone(),
            )
        },
        | Expr::Neg(mut arg) => {

            if matches!(
                *arg,
                Expr::Neg(_)
            ) && crate::is_exclusive(
                &arg,
            ) {

                let temp_arg = arg;

                match Arc::try_unwrap(
                    temp_arg,
                ) {
                    | Ok(
                        Expr::Neg(
                            inner_arc,
                        ),
                    ) => {

                        return Arc::try_unwrap(inner_arc).unwrap_or_else(|a| (*a).clone());
                    },
                    | Ok(other) => {

                        arg = Arc::new(
                            other,
                        );
                    },
                    | Err(
                        reclaimed_arg,
                    ) => {

                        arg = reclaimed_arg;
                    },
                }
            }

            if let Expr::Neg(
                ref inner_arg,
            ) = *arg
            {

                return inner_arg
                    .as_ref()
                    .clone();
            }

            if let Some(v) =
                as_f64(&arg)
            {

                return Expr::Constant(
                    -v,
                );
            }

            Expr::new_neg(arg)
        },
        | Expr::Log(arg) => {

            if let Some(value) =
                simplify_log(&arg)
            {

                return value;
            }

            Expr::new_log(arg)
        },
        | Expr::Exp(arg) => {

            if let Expr::Log(
                ref inner,
            ) = *arg
            {

                return inner
                    .as_ref()
                    .clone();
            }

            if is_zero(&arg) {

                return Expr::BigInt(
                    BigInt::one(),
                );
            }

            Expr::new_exp(arg)
        },
        | Expr::Sin(arg) => {

            if *arg == Expr::Pi {

                return Expr::BigInt(
                    BigInt::zero(),
                );
            }

            if let Expr::Neg(
                ref inner_arg,
            ) = *arg
            {

                return simplify(
                    Expr::new_neg(
                        Expr::new_sin(
                            inner_arg
                                .clone(
                                ),
                        ),
                    ),
                );
            }

            Expr::new_sin(arg)
        },
        | Expr::Cos(arg) => {

            if *arg == Expr::Pi {

                return Expr::new_neg(
                    Expr::BigInt(
                        BigInt::one(),
                    ),
                );
            }

            if let Expr::Neg(
                ref inner_arg,
            ) = *arg
            {

                return simplify(
                    Expr::new_cos(
                        inner_arg
                            .clone(),
                    ),
                );
            }

            Expr::new_cos(arg)
        },
        | Expr::Tan(arg) => {

            if *arg == Expr::Pi {

                return Expr::BigInt(
                    BigInt::zero(),
                );
            }

            if let Expr::Neg(
                ref inner_arg,
            ) = *arg
            {

                return simplify(
                    Expr::new_neg(
                        Expr::new_tan(
                            inner_arg
                                .clone(
                                ),
                        ),
                    ),
                );
            }

            Expr::new_tan(arg)
        },
        | Expr::Sum {
            body,
            var,
            from,
            to,
        } => {

            if let (
                Some(start),
                Some(end),
            ) = (
                as_f64(&from),
                as_f64(&to),
            ) {

                let mut total =
                    Expr::Constant(0.0);

                for i in (start.round()
                    as i64)
                    ..= (end.round()
                        as i64)
                {

                    let i_expr =
                        Expr::Constant(
                            i as f64,
                        );

                    if let Expr::Variable(ref v) = *var {

                        let term = substitute(&body, v, &i_expr);

                        total = simplify(Expr::new_add(
                            total, term,
                        ));
                    } else {

                        return Expr::Sum {
                            body,
                            var,
                            from,
                            to,
                        };
                    }
                }

                total
            } else {

                Expr::Sum {
                    body,
                    var,
                    from,
                    to,
                }
            }
        },
        // Handle N-ary list variants
        | Expr::AddList(terms) => {

            // Simplify each term
            let simplified_terms: Vec<
                Expr,
            > = terms
                .iter()
                .map(|t| {

                    simplify(t.clone())
                })
                .collect();

            // Flatten nested AddLists
            let mut flattened =
                Vec::new();

            for term in simplified_terms
            {

                if let Expr::AddList(
                    sub_terms,
                ) = term
                {

                    flattened.extend(
                        sub_terms,
                    );
                } else {

                    flattened
                        .push(term);
                }
            }

            // Filter out zeros and combine constants
            let mut constant_sum = 0.0;

            let mut non_constants =
                Vec::new();

            for term in flattened {

                if is_zero(&term) {

                    continue;
                }

                if let Some(val) =
                    as_f64(&term)
                {

                    constant_sum += val;
                } else {

                    non_constants
                        .push(term);
                }
            }

            // Build result
            if !non_constants.is_empty()
            {

                if constant_sum != 0.0 {

                    non_constants.insert(
                        0,
                        Expr::Constant(constant_sum),
                    );
                }

                if non_constants.len()
                    == 1
                {

                    non_constants[0]
                        .clone()
                } else {

                    Expr::AddList(
                        non_constants,
                    )
                }
            } else if constant_sum
                != 0.0
            {

                Expr::Constant(
                    constant_sum,
                )
            } else {

                Expr::BigInt(
                    BigInt::zero(),
                )
            }
        },
        | Expr::MulList(factors) => {

            // Simplify each factor
            let simplified_factors : Vec<Expr> = factors
                .iter()
                .map(|f| simplify(f.clone()))
                .collect();

            // Flatten nested MulLists
            let mut flattened =
                Vec::new();

            for factor in
                simplified_factors
            {

                if let Expr::MulList(
                    sub_factors,
                ) = factor
                {

                    flattened.extend(
                        sub_factors,
                    );
                } else {

                    flattened
                        .push(factor);
                }
            }

            // Filter out ones, check for zeros, and combine constants
            let mut constant_product =
                1.0;

            let mut non_constants =
                Vec::new();

            for factor in flattened {

                if is_zero(&factor) {

                    return Expr::BigInt(BigInt::zero());
                }

                if is_one(&factor) {

                    continue;
                }

                if let Some(val) =
                    as_f64(&factor)
                {

                    constant_product *=
                        val;
                } else {

                    non_constants
                        .push(factor);
                }
            }

            // Build result
            if !non_constants.is_empty()
            {

                if constant_product
                    != 1.0
                {

                    non_constants.insert(
                        0,
                        Expr::Constant(constant_product),
                    );
                }

                if non_constants.len()
                    == 1
                {

                    non_constants[0]
                        .clone()
                } else {

                    Expr::MulList(
                        non_constants,
                    )
                }
            } else if constant_product
                != 1.0
            {

                Expr::Constant(
                    constant_product,
                )
            } else {

                Expr::BigInt(
                    BigInt::one(),
                )
            }
        },
        // Generic list variants - simplify children
        | Expr::UnaryList(
            name,
            arg,
        ) => {
            Expr::UnaryList(
                name,
                Arc::new(simplify(
                    arg.as_ref()
                        .clone(),
                )),
            )
        },
        | Expr::BinaryList(
            name,
            a,
            b,
        ) => {
            Expr::BinaryList(
                name,
                Arc::new(simplify(
                    a.as_ref().clone(),
                )),
                Arc::new(simplify(
                    b.as_ref().clone(),
                )),
            )
        },
        | Expr::NaryList(
            name,
            args,
        ) => {

            let simplified_args: Vec<
                Expr,
            > = args
                .iter()
                .map(|arg| {

                    simplify(
                        arg.clone(),
                    )
                })
                .collect();

            Expr::NaryList(
                name,
                simplified_args,
            )
        },
        | _ => expr,
    }
}

#[inline]

pub(crate) fn simplify_log(
    arg: &Expr
) -> Option<Expr> {

    if let Expr::Complex(re, im) = &arg
    {

        let magnitude_sq =
            Expr::new_add(
                Expr::new_pow(
                    re.clone(),
                    Expr::Constant(2.0),
                ),
                Expr::new_pow(
                    im.clone(),
                    Expr::Constant(2.0),
                ),
            );

        let magnitude = Expr::new_sqrt(
            magnitude_sq,
        );

        let real_part =
            Expr::new_log(magnitude);

        let imag_part = Expr::new_atan2(
            im.clone(),
            re.clone(),
        );

        return Some(simplify(
            Expr::new_complex(
                real_part,
                imag_part,
            ),
        ));
    }

    if matches!(arg, Expr::E) {

        return Some(Expr::BigInt(
            BigInt::one(),
        ));
    }

    if let Expr::Exp(inner) = arg {

        return Some(
            inner
                .as_ref()
                .clone(),
        );
    }

    if is_one(arg) {

        return Some(Expr::BigInt(
            BigInt::zero(),
        ));
    }

    if let Expr::Power(base, exp) = arg
    {

        return Some(simplify(
            Expr::new_mul(
                exp.clone(),
                Expr::new_log(
                    base.clone(),
                ),
            ),
        ));
    }

    None
}

#[inline]

pub(crate) fn simplify_sqrt(
    arg: Expr
) -> Expr {

    let simplified_arg = simplify(arg);

    let denested = crate::symbolic::radicals::denest_sqrt(&Expr::new_sqrt(
        simplified_arg.clone(),
    ));

    if let Expr::Sqrt(_) = denested {

        if let Expr::Power(
            ref b,
            ref e,
        ) = simplified_arg
        {

            if let Some(val) = as_f64(e)
            {

                return simplify(
                    Expr::new_pow(
                        b.clone(),
                        Expr::Constant(
                            val / 2.0,
                        ),
                    ),
                );
            }
        }

        Expr::new_sqrt(simplified_arg)
    } else {

        denested
    }
}

#[inline]

pub(crate) fn simplify_power(
    b: &Expr,
    e: &Expr,
) -> Option<Expr> {

    if let (Some(vb), Some(ve)) =
        (as_f64(b), as_f64(e))
    {

        return Some(Expr::Constant(
            vb.powf(ve),
        ));
    }

    if is_zero(e) {

        return Some(Expr::BigInt(
            BigInt::one(),
        ));
    }

    if is_one(e) {

        return Some(b.clone());
    }

    if is_zero(b) {

        return Some(Expr::BigInt(
            BigInt::zero(),
        ));
    }

    if is_one(b) {

        return Some(Expr::BigInt(
            BigInt::one(),
        ));
    }

    if let Expr::Power(
        inner_b,
        inner_e,
    ) = b
    {

        return Some(simplify(
            Expr::new_pow(
                inner_b.clone(),
                Expr::new_mul(
                    inner_e.clone(),
                    e.clone(),
                ),
            ),
        ));
    }

    if let Expr::Exp(base_inner) = b {

        return Some(simplify(
            Expr::new_exp(
                Expr::new_mul(
                    base_inner.clone(),
                    e.clone(),
                ),
            ),
        ));
    }

    None
}

#[inline]

pub(crate) fn simplify_div(
    a: &Expr,
    b: &Expr,
) -> Option<Expr> {

    if let (Some(va), Some(vb)) =
        (as_f64(a), as_f64(b))
    {

        if vb != 0.0 {

            return Some(
                Expr::Constant(va / vb),
            );
        }
    }

    if is_zero(a) {

        return Some(Expr::BigInt(
            BigInt::zero(),
        ));
    }

    if is_one(b) {

        return Some(a.clone());
    }

    if *a == *b {

        return Some(Expr::BigInt(
            BigInt::one(),
        ));
    }

    None
}

#[inline]

pub(crate) fn simplify_mul(
    a: &Expr,
    b: &Expr,
) -> Option<Expr> {

    if let (Some(va), Some(vb)) =
        (as_f64(a), as_f64(b))
    {

        return Some(Expr::Constant(
            va * vb,
        ));
    }

    if is_zero(a) || is_zero(b) {

        return Some(Expr::BigInt(
            BigInt::zero(),
        ));
    }

    if is_one(a) {

        return Some(b.clone());
    }

    if is_one(b) {

        return Some(a.clone());
    }

    if let (
        Expr::Exp(a_inner),
        Expr::Exp(b_inner),
    ) = (&a, &b)
    {

        return Some(simplify(
            Expr::new_exp(
                Expr::new_add(
                    a_inner.clone(),
                    b_inner.clone(),
                ),
            ),
        ));
    }

    if let (
        Expr::Power(base1, exp1),
        Expr::Power(base2, exp2),
    ) = (&a, &b)
    {

        if base1 == base2 {

            return Some(simplify(
                Expr::new_pow(
                    base1.clone(),
                    Expr::new_add(
                        exp1.clone(),
                        exp2.clone(),
                    ),
                ),
            ));
        }
    }

    if let Expr::Add(b_inner, c_inner) =
        b
    {

        return Some(simplify(
            Expr::new_add(
                Expr::new_mul(
                    a.clone(),
                    b_inner.clone(),
                ),
                Expr::new_mul(
                    a.clone(),
                    c_inner.clone(),
                ),
            ),
        ));
    }

    None
}

#[inline]
#[allow(unused_allocation)]

pub(crate) fn simplify_sub(
    a: &Expr,
    b: &Expr,
) -> Option<Expr> {

    if let (Some(va), Some(vb)) =
        (as_f64(a), as_f64(b))
    {

        return Some(Expr::Constant(
            va - vb,
        ));
    }

    if is_zero(b) {

        return Some(a.clone());
    }

    if *a == *b {

        return Some(Expr::BigInt(
            BigInt::zero(),
        ));
    }

    if is_one(a) {

        if let Expr::Power(base, exp) =
            b
        {

            let two = Expr::BigInt(
                BigInt::from(2),
            );

            let two_f =
                Expr::Constant(2.0);

            if *exp == Arc::new(two)
                || *exp
                    == Arc::new(two_f)
            {

                if let Expr::Cos(arg) =
                    &**base
                {

                    return Some(simplify(
                        Expr::new_pow(
                            Expr::new_sin(arg.clone()),
                            Expr::Constant(2.0),
                        ),
                    ));
                }

                if let Expr::Sin(arg) =
                    &**base
                {

                    return Some(simplify(
                        Expr::new_pow(
                            Expr::new_cos(arg.clone()),
                            Expr::Constant(2.0),
                        ),
                    ));
                }
            }
        }
    }

    None
}

#[inline]

pub(crate) fn simplify_add(
    a: Expr,
    b: Expr,
) -> Result<Expr, Expr> {

    if let (
        Expr::BigInt(ia),
        Expr::BigInt(ib),
    ) = (&a, &b)
    {

        return Err(Expr::BigInt(
            ia + ib,
        ));
    }

    if let (
        Expr::Rational(ra),
        Expr::Rational(rb),
    ) = (&a, &b)
    {

        return Err(Expr::Rational(
            ra + rb,
        ));
    }

    if let (Some(va), Some(vb)) = (
        as_f64(&a),
        as_f64(&b),
    ) {

        return Err(Expr::Constant(
            va + vb,
        ));
    }

    let original_expr =
        Expr::new_add(a, b);

    let (constant_term, terms) =
        collect_and_order_terms(
            &original_expr,
        );

    let mut term_iter = terms
        .into_iter()
        .filter(|(_, coeff)| {

            !is_zero(coeff)
        });

    let mut result_expr =
        match term_iter.next() {
            | Some((base, coeff)) => {

                let first_term =
                    if is_one(&coeff) {

                        base
                    } else {

                        simplify(Expr::new_mul(
                    coeff, base,
                ))
                    };

                if is_zero(
                    &constant_term,
                ) {

                    first_term
                } else {

                    simplify(Expr::new_add(
                    constant_term,
                    first_term,
                ))
                }
            },
            | None => constant_term,
        };

    for (base, coeff) in term_iter {

        let term = if is_one(&coeff) {

            base
        } else {

            simplify(Expr::new_mul(
                coeff, base,
            ))
        };

        result_expr =
            simplify(Expr::new_add(
                result_expr,
                term,
            ));
    }

    Ok(result_expr)
}

pub struct RewriteRule {
    name: &'static str,
    pattern: Expr,
    replacement: Expr,
}

#[must_use]

pub fn get_name(
    rule: &RewriteRule
) -> String {

    println!("{}", rule.name);

    rule.name
        .to_string()
}

pub(crate) fn get_default_rules(
) -> Vec<RewriteRule> {

    vec![
        RewriteRule {
            name : "factor_common_term",
            pattern : Expr::Add(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern(
                        "a".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "b".to_string(),
                    )),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern(
                        "a".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "c".to_string(),
                    )),
                )),
            ),
            replacement : Expr::Mul(
                Arc::new(Expr::Pattern(
                    "a".to_string(),
                )),
                Arc::new(Expr::Add(
                    Arc::new(Expr::Pattern(
                        "b".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "c".to_string(),
                    )),
                )),
            ),
        },
        RewriteRule {
            name : "distribute_mul_add",
            pattern : Expr::Mul(
                Arc::new(Expr::Pattern(
                    "a".to_string(),
                )),
                Arc::new(Expr::Add(
                    Arc::new(Expr::Pattern(
                        "b".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "c".to_string(),
                    )),
                )),
            ),
            replacement : Expr::Add(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern(
                        "a".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "b".to_string(),
                    )),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Pattern(
                        "a".to_string(),
                    )),
                    Arc::new(Expr::Pattern(
                        "c".to_string(),
                    )),
                )),
            ),
        },
        RewriteRule {
            name : "tan_to_sin_cos",
            pattern : Expr::Tan(Arc::new(
                Expr::Pattern("x".to_string()),
            )),
            replacement : Expr::Div(
                Arc::new(Expr::Sin(Arc::new(
                    Expr::Pattern("x".to_string()),
                ))),
                Arc::new(Expr::Cos(Arc::new(
                    Expr::Pattern("x".to_string()),
                ))),
            ),
        },
        RewriteRule {
            name : "sin_cos_to_tan",
            pattern : Expr::Div(
                Arc::new(Expr::Sin(Arc::new(
                    Expr::Pattern("x".to_string()),
                ))),
                Arc::new(Expr::Cos(Arc::new(
                    Expr::Pattern("x".to_string()),
                ))),
            ),
            replacement : Expr::Tan(Arc::new(
                Expr::Pattern("x".to_string()),
            )),
        },
        RewriteRule {
            name : "double_angle_sin",
            pattern : Expr::Mul(
                Arc::new(Expr::BigInt(
                    BigInt::from(2),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::Sin(Arc::new(
                        Expr::Pattern("x".to_string()),
                    ))),
                    Arc::new(Expr::Cos(Arc::new(
                        Expr::Pattern("x".to_string()),
                    ))),
                )),
            ),
            replacement : Expr::Sin(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(
                    BigInt::from(2),
                )),
                Arc::new(Expr::Pattern(
                    "x".to_string(),
                )),
            ))),
        },
        RewriteRule {
            name : "double_angle_cos_1",
            pattern : Expr::Sub(
                Arc::new(Expr::Power(
                    Arc::new(Expr::Cos(Arc::new(
                        Expr::Pattern("x".to_string()),
                    ))),
                    Arc::new(Expr::BigInt(
                        BigInt::from(2),
                    )),
                )),
                Arc::new(Expr::Power(
                    Arc::new(Expr::Sin(Arc::new(
                        Expr::Pattern("x".to_string()),
                    ))),
                    Arc::new(Expr::BigInt(
                        BigInt::from(2),
                    )),
                )),
            ),
            replacement : Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(
                    BigInt::from(2),
                )),
                Arc::new(Expr::Pattern(
                    "x".to_string(),
                )),
            ))),
        },
        RewriteRule {
            name : "double_angle_cos_2",
            pattern : Expr::Sub(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::BigInt(
                        BigInt::from(2),
                    )),
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Cos(Arc::new(
                            Expr::Pattern("x".to_string()),
                        ))),
                        Arc::new(Expr::BigInt(
                            BigInt::from(2),
                        )),
                    )),
                )),
                Arc::new(Expr::BigInt(
                    BigInt::from(1),
                )),
            ),
            replacement : Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(
                    BigInt::from(2),
                )),
                Arc::new(Expr::Pattern(
                    "x".to_string(),
                )),
            ))),
        },
        RewriteRule {
            name : "double_angle_cos_3",
            pattern : Expr::Sub(
                Arc::new(Expr::BigInt(
                    BigInt::from(1),
                )),
                Arc::new(Expr::Mul(
                    Arc::new(Expr::BigInt(
                        BigInt::from(2),
                    )),
                    Arc::new(Expr::Power(
                        Arc::new(Expr::Sin(Arc::new(
                            Expr::Pattern("x".to_string()),
                        ))),
                        Arc::new(Expr::BigInt(
                            BigInt::from(2),
                        )),
                    )),
                )),
            ),
            replacement : Expr::Cos(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(
                    BigInt::from(2),
                )),
                Arc::new(Expr::Pattern(
                    "x".to_string(),
                )),
            ))),
        },
    ]
}

#[must_use]

pub fn substitute_patterns(
    template: &Expr,
    assignments: &HashMap<String, Expr>,
) -> Expr {

    match template {
        | Expr::Pattern(name) => {
            assignments
                .get(name)
                .cloned()
                .unwrap_or_else(|| {

                    template.clone()
                })
        },
        | Expr::Add(a, b) => {
            Expr::new_add(
                substitute_patterns(
                    a,
                    assignments,
                ),
                substitute_patterns(
                    b,
                    assignments,
                ),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                substitute_patterns(
                    a,
                    assignments,
                ),
                substitute_patterns(
                    b,
                    assignments,
                ),
            )
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                substitute_patterns(
                    a,
                    assignments,
                ),
                substitute_patterns(
                    b,
                    assignments,
                ),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                substitute_patterns(
                    a,
                    assignments,
                ),
                substitute_patterns(
                    b,
                    assignments,
                ),
            )
        },
        | Expr::Power(b, e) => {
            Expr::new_pow(
                substitute_patterns(
                    b,
                    assignments,
                ),
                substitute_patterns(
                    e,
                    assignments,
                ),
            )
        },
        | Expr::Sin(arg) => {
            Expr::new_sin(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | Expr::Cos(arg) => {
            Expr::new_cos(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | Expr::Tan(arg) => {
            Expr::new_tan(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | Expr::Exp(arg) => {
            Expr::new_exp(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | Expr::Log(arg) => {
            Expr::new_log(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | Expr::Neg(arg) => {
            Expr::new_neg(
                substitute_patterns(
                    arg,
                    assignments,
                ),
            )
        },
        | _ => template.clone(),
    }
}

pub(crate) fn apply_rules_recursively(
    expr: &Expr,
    rules: &[RewriteRule],
) -> (Expr, bool) {

    let mut current_expr = expr.clone();

    let mut changed = false;

    let simplified_children =
        match &current_expr {
            | Expr::Add(a, b) => {

                let (na, ca) = apply_rules_recursively(a, rules);

                let (nb, cb) = apply_rules_recursively(b, rules);

                if ca || cb {

                    Some(Expr::new_add(
                        na, nb,
                    ))
                } else {

                    None
                }
            },
            | Expr::Sub(a, b) => {

                let (na, ca) = apply_rules_recursively(a, rules);

                let (nb, cb) = apply_rules_recursively(b, rules);

                if ca || cb {

                    Some(Expr::new_sub(
                        na, nb,
                    ))
                } else {

                    None
                }
            },
            | Expr::Mul(a, b) => {

                let (na, ca) = apply_rules_recursively(a, rules);

                let (nb, cb) = apply_rules_recursively(b, rules);

                if ca || cb {

                    Some(Expr::new_mul(
                        na, nb,
                    ))
                } else {

                    None
                }
            },
            | Expr::Div(a, b) => {

                let (na, ca) = apply_rules_recursively(a, rules);

                let (nb, cb) = apply_rules_recursively(b, rules);

                if ca || cb {

                    Some(Expr::new_div(
                        na, nb,
                    ))
                } else {

                    None
                }
            },
            | Expr::Power(b, e) => {

                let (nb, cb) = apply_rules_recursively(b, rules);

                let (ne, ce) = apply_rules_recursively(e, rules);

                if cb || ce {

                    Some(Expr::new_pow(
                        nb, ne,
                    ))
                } else {

                    None
                }
            },
            | Expr::Sin(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_sin(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | Expr::Cos(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_cos(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | Expr::Tan(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_tan(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | Expr::Exp(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_exp(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | Expr::Log(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_log(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | Expr::Neg(arg) => {

                let (narg, carg) = apply_rules_recursively(arg, rules);

                if carg {

                    Some(Expr::new_neg(
                        narg,
                    ))
                } else {

                    None
                }
            },
            | _ => None,
        };

    if let Some(new_expr) =
        simplified_children
    {

        current_expr = new_expr;

        changed = true;
    }

    for rule in rules {

        if let Some(assignments) =
            pattern_match(
                &current_expr,
                &rule.pattern,
            )
        {

            let new_expr =
                substitute_patterns(
                    &rule.replacement,
                    &assignments,
                );

            let simplified_new_expr =
                simplify(new_expr);

            if complexity(
                &simplified_new_expr,
            ) < complexity(
                &current_expr,
            ) {

                current_expr =
                    simplified_new_expr;

                changed = true;
            }
        }
    }

    (
        current_expr,
        changed,
    )
}

/// Applies a set of heuristic transformations to find a simpler form of an expression.
///
/// This function uses pattern matching and rewrite rules to transform the expression.
/// It iteratively applies rules until a fixed point is reached or a maximum number
/// of iterations is exceeded. After each pass of rule application, it performs a
/// deterministic simplification using `simplify`.
///
/// # Arguments
/// * `expr` - The expression to heuristically simplify.
///
/// # Returns
/// A new, heuristically simplified `Expr`.
#[must_use]

pub fn heuristic_simplify(
    expr: Expr
) -> Expr {

    let mut current_expr = expr;

    let rules = get_default_rules();

    const MAX_ITERATIONS: usize = 10;

    for _ in 0 .. MAX_ITERATIONS {

        let (next_expr, changed) =
            apply_rules_recursively(
                &current_expr,
                &rules,
            );

        current_expr =
            simplify(next_expr);

        if !changed {

            break;
        }
    }

    current_expr
}

pub(crate) fn complexity(
    expr: &Expr
) -> usize {

    match expr {
        | Expr::BigInt(_) => 1,
        | Expr::Rational(_) => 2,
        | Expr::Constant(_) => 3,
        | Expr::Variable(_)
        | Expr::Pattern(_) => 5,
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b) => {
            complexity(a)
                + complexity(b)
                + 1
        },
        | Expr::Power(a, b) => {
            complexity(a)
                + complexity(b)
                + 2
        },
        | Expr::Sin(a)
        | Expr::Cos(a)
        | Expr::Tan(a)
        | Expr::Exp(a)
        | Expr::Log(a)
        | Expr::Neg(a) => {
            complexity(a) + 3
        },
        | _ => 100,
    }
}

/// Attempts to match an expression against a pattern.
///
/// If a match is found, it returns a `HashMap` containing the assignments
/// for the pattern variables. Pattern variables are represented by `Expr::Pattern(name)`.
///
/// # Arguments
/// * `expr` - The expression to match.
/// * `pattern` - The pattern to match against.
///
/// # Returns
/// `Some(HashMap<String, Expr>)` with variable assignments if a match is found,
/// `None` otherwise.
#[inline]
#[must_use]

pub fn pattern_match(
    expr: &Expr,
    pattern: &Expr,
) -> Option<HashMap<String, Expr>> {

    let mut assignments =
        HashMap::new();

    if pattern_match_recursive(
        expr,
        pattern,
        &mut assignments,
    ) {

        Some(assignments)
    } else {

        None
    }
}

pub(crate) fn pattern_match_recursive(
    expr: &Expr,
    pattern: &Expr,
    assignments: &mut HashMap<
        String,
        Expr,
    >,
) -> bool {

    match (expr, pattern) {
        | (_, Expr::Pattern(name)) => {

            if let Some(existing) =
                assignments.get(name)
            {

                return existing
                    == expr;
            }

            assignments.insert(
                name.clone(),
                expr.clone(),
            );

            true
        },
        | (
            Expr::Add(e1, e2),
            Expr::Add(p1, p2),
        )
        | (
            Expr::Mul(e1, e2),
            Expr::Mul(p1, p2),
        ) => {

            let original_assignments =
                assignments.clone();

            if pattern_match_recursive(e1, p1, assignments)
                && pattern_match_recursive(e2, p2, assignments)
            {

                return true;
            }

            *assignments =
                original_assignments;

            pattern_match_recursive(
                e1,
                p2,
                assignments,
            ) && pattern_match_recursive(
                e2,
                p1,
                assignments,
            )
        },
        | (
            Expr::Sub(e1, e2),
            Expr::Sub(p1, p2),
        )
        | (
            Expr::Div(e1, e2),
            Expr::Div(p1, p2),
        )
        | (
            Expr::Power(e1, e2),
            Expr::Power(p1, p2),
        ) => {
            pattern_match_recursive(
                e1,
                p1,
                assignments,
            ) && pattern_match_recursive(
                e2,
                p2,
                assignments,
            )
        },
        | (
            Expr::Sin(e),
            Expr::Sin(p),
        )
        | (
            Expr::Cos(e),
            Expr::Cos(p),
        )
        | (
            Expr::Tan(e),
            Expr::Tan(p),
        )
        | (
            Expr::Exp(e),
            Expr::Exp(p),
        )
        | (
            Expr::Log(e),
            Expr::Log(p),
        )
        | (
            Expr::Neg(e),
            Expr::Neg(p),
        ) => {
            pattern_match_recursive(
                e,
                p,
                assignments,
            )
        },
        | _ => expr == pattern,
    }
}

#[must_use]

pub fn collect_and_order_terms(
    expr: &Expr
) -> (
    Expr,
    Vec<(Expr, Expr)>,
) {

    /// Collects terms from an expression and orders them by complexity.
    ///
    /// This function is useful for canonicalizing expressions, especially sums and differences.
    /// It extracts a constant term and a vector of `(base, coefficient)` pairs for other terms.
    /// Terms are ordered heuristically by their complexity.
    ///
    /// # Arguments
    /// * `expr` - The expression to collect terms from.
    ///
    /// # Returns
    /// A tuple `(constant_term, terms)` where `constant_term` is an `Expr` and `terms` is a
    /// `Vec<(Expr, Expr)>` of `(base, coefficient)` pairs.
    let mut terms = BTreeMap::new();

    collect_terms_recursive(
        expr,
        &Expr::BigInt(BigInt::one()),
        &mut terms,
    );

    let mut sorted_terms: Vec<(
        Expr,
        Expr,
    )> = terms
        .into_iter()
        .collect();

    sorted_terms.sort_by(
        |(b1, _), (b2, _)| {

            complexity(b2)
                .cmp(&complexity(b1))
        },
    );

    let constant_term =
        if let Some(pos) = sorted_terms
            .iter()
            .position(|(b, _)| {

                is_one(b)
            })
        {

            let (_, c) = sorted_terms
                .remove(pos);

            c
        } else {

            Expr::BigInt(BigInt::zero())
        };

    (
        constant_term,
        sorted_terms,
    )
}

pub(crate) fn fold_constants(
    expr: Expr
) -> Expr {

    let expr = match expr {
        | Expr::Add(a, b) => {
            Expr::new_add(
                fold_constants(
                    a.as_ref().clone(),
                ),
                fold_constants(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Sub(a, b) => {
            Expr::new_sub(
                fold_constants(
                    a.as_ref().clone(),
                ),
                fold_constants(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Mul(a, b) => {
            Expr::new_mul(
                fold_constants(
                    a.as_ref().clone(),
                ),
                fold_constants(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Div(a, b) => {
            Expr::new_div(
                fold_constants(
                    a.as_ref().clone(),
                ),
                fold_constants(
                    b.as_ref().clone(),
                ),
            )
        },
        | Expr::Power(base, exp) => {
            Expr::new_pow(
                fold_constants(
                    (*base).clone(),
                ),
                fold_constants(
                    (*exp).clone(),
                ),
            )
        },
        | Expr::Neg(arg) => {
            Expr::new_neg(
                fold_constants(
                    (*arg).clone(),
                ),
            )
        },
        | _ => expr,
    };

    match expr {
        | Expr::Add(a, b) => {

            if let (
                Some(va),
                Some(vb),
            ) = (
                as_f64(&a),
                as_f64(&b),
            ) {

                Expr::Constant(va + vb)
            } else {

                Expr::new_add(a, b)
            }
        },
        | Expr::Sub(a, b) => {

            if let (
                Some(va),
                Some(vb),
            ) = (
                as_f64(&a),
                as_f64(&b),
            ) {

                Expr::Constant(va - vb)
            } else {

                Expr::new_sub(a, b)
            }
        },
        | Expr::Mul(a, b) => {

            if let (
                Some(va),
                Some(vb),
            ) = (
                as_f64(&a),
                as_f64(&b),
            ) {

                Expr::Constant(va * vb)
            } else {

                Expr::new_mul(a, b)
            }
        },
        | Expr::Div(a, b) => {

            if let (
                Some(va),
                Some(vb),
            ) = (
                as_f64(&a),
                as_f64(&b),
            ) {

                if vb == 0.0 {

                    Expr::new_div(a, b)
                } else {

                    Expr::Constant(
                        va / vb,
                    )
                }
            } else {

                Expr::new_div(a, b)
            }
        },
        | Expr::Power(b, e) => {

            if let (
                Some(vb),
                Some(ve),
            ) = (
                as_f64(&b),
                as_f64(&e),
            ) {

                Expr::Constant(
                    vb.powf(ve),
                )
            } else {

                Expr::new_pow(b, e)
            }
        },
        | Expr::Neg(arg) => {
            if let Some(v) =
                as_f64(&arg)
            {

                Expr::Constant(-v)
            } else {

                Expr::new_neg(arg)
            }
        },
        | _ => expr,
    }
}

#[must_use]

pub const fn is_numeric(
    expr: &Expr
) -> bool {

    matches!(
        expr,
        Expr::Constant(_)
            | Expr::BigInt(_)
            | Expr::Rational(_)
    )
}

pub(crate) fn collect_terms_recursive(
    expr: &Expr,
    coeff: &Expr,
    terms: &mut BTreeMap<Expr, Expr>,
) {

    let mut stack = vec![(
        expr.clone(),
        coeff.clone(),
    )];

    while let Some((
        current_expr,
        current_coeff,
    )) = stack.pop()
    {

        match &current_expr {
            | Expr::Add(a, b) => {

                stack.push((
                    a.as_ref().clone(),
                    current_coeff
                        .clone(),
                ));

                stack.push((
                    b.as_ref().clone(),
                    current_coeff,
                ));
            },
            | Expr::AddList(
                terms_list,
            ) => {

                // Flatten AddList by pushing all terms with the same coefficient
                for term in terms_list {

                    stack.push((
                        term.clone(),
                        current_coeff
                            .clone(),
                    ));
                }
            },
            | Expr::Sub(a, b) => {

                stack.push((
                    a.as_ref().clone(),
                    current_coeff
                        .clone(),
                ));

                stack.push((
                    b.as_ref().clone(),
                    fold_constants(Expr::new_neg(
                        current_coeff,
                    )),
                ));
            },
            | Expr::Mul(a, b) => {
                if is_numeric(a) {

                    stack.push((
                        b.as_ref().clone(),
                        fold_constants(Expr::new_mul(
                            current_coeff,
                            a.as_ref().clone(),
                        )),
                    ));
                } else if is_numeric(b)
                {

                    stack.push((
                        a.as_ref().clone(),
                        fold_constants(Expr::new_mul(
                            current_coeff,
                            b.as_ref().clone(),
                        )),
                    ));
                } else {

                    let base =
                        current_expr;

                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                    *entry = fold_constants(Expr::new_add(
                        entry.clone(),
                        current_coeff,
                    ));
                }
            },
            | Expr::MulList(
                factors,
            ) => {

                // Try to extract numeric coefficient from MulList
                let mut numeric_part =
                    Expr::BigInt(
                        BigInt::one(),
                    );

                let mut
                non_numeric_parts =
                    Vec::new();

                for factor in factors {

                    if is_numeric(
                        factor,
                    ) {

                        numeric_part = fold_constants(Expr::new_mul(
                            numeric_part,
                            factor.clone(),
                        ));
                    } else {

                        non_numeric_parts.push(factor.clone());
                    }
                }

                if non_numeric_parts
                    .is_empty()
                {

                    // All factors are numeric
                    let base =
                        Expr::BigInt(
                            BigInt::one(
                            ),
                        );

                    let new_coeff = fold_constants(Expr::new_mul(
                        current_coeff,
                        numeric_part,
                    ));

                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                    *entry = fold_constants(Expr::new_add(
                        entry.clone(),
                        new_coeff,
                    ));
                } else {

                    let base = if non_numeric_parts.len() == 1 {

                        non_numeric_parts[0].clone()
                    } else {

                        Expr::MulList(non_numeric_parts)
                    };

                    let new_coeff = fold_constants(Expr::new_mul(
                        current_coeff,
                        numeric_part,
                    ));

                    let entry = terms
                        .entry(base)
                        .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                    *entry = fold_constants(Expr::new_add(
                        entry.clone(),
                        new_coeff,
                    ));
                }
            },
            | _ => {

                let base = current_expr;

                let entry = terms
                    .entry(base)
                    .or_insert_with(|| Expr::BigInt(BigInt::zero()));

                *entry = fold_constants(
                    Expr::new_add(
                        entry.clone(),
                        current_coeff,
                    ),
                );
            },
        }
    }
}

#[inline]

pub(crate) fn as_rational(
    expr: &Expr
) -> (Expr, Expr) {

    if let Expr::Div(num, den) = expr {

        (
            num.as_ref().clone(),
            den.as_ref().clone(),
        )
    } else {

        (
            expr.clone(),
            Expr::Constant(1.0),
        )
    }
}

pub(crate) fn simplify_rational_expression(
    expr: &Expr
) -> Expr {

    if let Expr::Add(a, b)
    | Expr::Sub(a, b)
    | Expr::Mul(a, b)
    | Expr::Div(a, b) = expr
    {

        let (num1, den1) =
            as_rational(a);

        let (num2, den2) =
            as_rational(b);

        let (
            new_num_expr,
            new_den_expr,
        ) = match expr {
            | Expr::Add(_, _) => (
                apply_rules(
                    Expr::new_add(
                        Expr::new_mul(
                            num1,
                            den2.clone(
                            ),
                        ),
                        Expr::new_mul(
                            num2,
                            den1.clone(
                            ),
                        ),
                    ),
                ),
                apply_rules(
                    Expr::new_mul(
                        den1, den2,
                    ),
                ),
            ),
            | Expr::Sub(_, _) => (
                apply_rules(
                    Expr::new_sub(
                        Expr::new_mul(
                            num1,
                            den2.clone(
                            ),
                        ),
                        Expr::new_mul(
                            num2,
                            den1.clone(
                            ),
                        ),
                    ),
                ),
                apply_rules(
                    Expr::new_mul(
                        den1, den2,
                    ),
                ),
            ),
            | Expr::Mul(_, _) => {
                (
                    apply_rules(
                        Expr::new_mul(
                            num1, num2,
                        ),
                    ),
                    apply_rules(
                        Expr::new_mul(
                            den1, den2,
                        ),
                    ),
                )
            },
            | Expr::Div(_, _) => {
                (
                    apply_rules(
                        Expr::new_mul(
                            num1, den2,
                        ),
                    ),
                    apply_rules(
                        Expr::new_mul(
                            den1, num2,
                        ),
                    ),
                )
            },
            | _ => unreachable!(),
        };

        if is_one(&new_den_expr) {

            return new_num_expr;
        }

        if is_zero(&new_num_expr) {

            return Expr::Constant(0.0);
        }

        let var = "x";

        let p_num = crate::symbolic::polynomial::expr_to_sparse_poly(
            &new_num_expr,
            &[var],
        );

        let p_den = crate::symbolic::polynomial::expr_to_sparse_poly(
            &new_den_expr,
            &[var],
        );

        let common_divisor = crate::symbolic::polynomial::gcd(
            p_num.clone(),
            p_den.clone(),
            var,
        );

        if common_divisor.degree(var)
            > 0
        {

            let final_num_poly = p_num
                .long_division(
                    common_divisor
                        .clone(),
                    var,
                )
                .0;

            let final_den_poly = p_den
                .long_division(
                    common_divisor,
                    var,
                )
                .0;

            let final_num = crate::symbolic::polynomial::sparse_poly_to_expr(&final_num_poly);

            let final_den = crate::symbolic::polynomial::sparse_poly_to_expr(&final_den_poly);

            if is_one(&final_den) {

                return final_num;
            }

            return Expr::new_div(
                final_num,
                final_den,
            );
        }

        return Expr::new_div(
            new_num_expr,
            new_den_expr,
        );
    }

    expr.clone()
}
