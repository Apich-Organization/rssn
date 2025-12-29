//! # Unit Unification and Dimensional Analysis
//!
//! This module provides functionality for handling physical units and performing
//! dimensional analysis within symbolic expressions. It supports parsing unit strings,
//! unifying expressions involving quantities with units, and performing arithmetic
//! operations while respecting dimensional consistency.
//!
//! # Supported Units
//! Currently, the module supports basic SI units and some common derivatives:
//! - **Length**: meter (m), centimeter (cm)
//! - **Mass**: kilogram (kg), gram (g)
//! - **Time**: second (s), minute (min)
//! - **Area**: square meter (m², sqm)
//! - **Velocity**: meter per second (m/s, mps)
//!
//! # Examples
//! ```
//! 
//! use std::sync::Arc;
//!
//! use rssn::symbolic::core::Expr;
//! use rssn::symbolic::unit_unification::unify_expression;
//!
//! // Define quantities: 5m and 3m
//! let q1 = Expr::QuantityWithValue(
//!     Arc::new(Expr::new_constant(
//!         5.0,
//!     )),
//!     "m".to_string(),
//! );
//!
//! let q2 = Expr::QuantityWithValue(
//!     Arc::new(Expr::new_constant(
//!         3.0,
//!     )),
//!     "m".to_string(),
//! );
//!
//! // Add them: 5m + 3m
//! let sum_expr = Expr::new_add(q1, q2);
//!
//! // Unify to get the result: 8m
//! let result = unify_expression(&sum_expr).unwrap();
//! ```

use std::fmt::Debug;
use std::hash::Hash;
use std::hash::Hasher;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;
use std::sync::Arc;

use ordered_float;
use uom::si::area;
use uom::si::f64::Area;
use uom::si::f64::Length;
use uom::si::f64::Mass;
use uom::si::f64::Time;
use uom::si::f64::Velocity;
use uom::si::length;
use uom::si::mass;
use uom::si::time;
use uom::si::velocity;

use crate::symbolic::core::Expr;

/// Represents a physical quantity with a specific dimension.
///
/// This enum wraps the `uom` crate's types to provide a unified interface for
/// different physical dimensions.
#[derive(
    Debug,
    Clone,
    PartialEq,
    serde::Serialize,
    serde::Deserialize,
)]
#[repr(C)]

pub enum SupportedQuantity {
    /// Represents a length quantity.
    Length(Length),
    /// Represents a mass quantity.
    Mass(Mass),
    /// Represents a time quantity.
    Time(Time),
    /// Represents an area quantity.
    Area(Area),
    /// Represents a velocity quantity.
    Velocity(Velocity),
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Add for SupportedQuantity {
    type Output = Result<Self, String>;

    fn add(
        self,
        rhs: Self,
    ) -> Self::Output {

        match (self, rhs) {
            | (
                Self::Length(l1),
                Self::Length(l2),
            ) => {
                Ok(Self::Length(
                    l1 + l2,
                ))
            },
            | (
                Self::Mass(m1),
                Self::Mass(m2),
            ) => {
                Ok(Self::Mass(m1 + m2))
            },
            | (
                Self::Time(t1),
                Self::Time(t2),
            ) => {
                Ok(Self::Time(t1 + t2))
            },
            | (
                Self::Area(a1),
                Self::Area(a2),
            ) => {
                Ok(Self::Area(a1 + a2))
            },
            | (
                Self::Velocity(v1),
                Self::Velocity(v2),
            ) => {
                Ok(Self::Velocity(
                    v1 + v2,
                ))
            },
            | _ => Err(
                "Incompatible types \
                 for addition"
                    .to_string(),
            ),
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Sub for SupportedQuantity {
    type Output = Result<Self, String>;

    fn sub(
        self,
        rhs: Self,
    ) -> Self::Output {

        match (self, rhs) {
            | (
                Self::Length(l1),
                Self::Length(l2),
            ) => {
                Ok(Self::Length(
                    l1 - l2,
                ))
            },
            | (
                Self::Mass(m1),
                Self::Mass(m2),
            ) => {
                Ok(Self::Mass(m1 - m2))
            },
            | (
                Self::Time(t1),
                Self::Time(t2),
            ) => {
                Ok(Self::Time(t1 - t2))
            },
            | (
                Self::Area(a1),
                Self::Area(a2),
            ) => {
                Ok(Self::Area(a1 - a2))
            },
            | (
                Self::Velocity(v1),
                Self::Velocity(v2),
            ) => {
                Ok(Self::Velocity(
                    v1 - v2,
                ))
            },
            | _ => Err(
                "Incompatible types \
                 for subtraction"
                    .to_string(),
            ),
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Neg for SupportedQuantity {
    type Output = Self;

    fn neg(self) -> Self::Output {

        match self {
            | Self::Length(l) => {
                Self::Length(-l)
            },
            | Self::Mass(m) => {
                Self::Mass(-m)
            },
            | Self::Time(t) => {
                Self::Time(-t)
            },
            | Self::Area(a) => {
                Self::Area(-a)
            },
            | Self::Velocity(v) => {
                Self::Velocity(-v)
            },
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Mul for SupportedQuantity {
    type Output = Result<Self, String>;

    fn mul(
        self,
        rhs: Self,
    ) -> Self::Output {

        match (self, rhs) {
            | (
                Self::Length(l1),
                Self::Length(l2),
            ) => {
                Ok(Self::Area(l1 * l2))
            },
            | _ => {
                Err("Unsupported or \
                     complex quantity \
                     combination for \
                     multiplication"
                    .to_string())
            },
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Mul<f64> for SupportedQuantity {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self::Output {

        match self {
            | Self::Length(l) => {
                Self::Length(l * rhs)
            },
            | Self::Mass(m) => {
                Self::Mass(m * rhs)
            },
            | Self::Time(t) => {
                Self::Time(t * rhs)
            },
            | Self::Area(a) => {
                Self::Area(a * rhs)
            },
            | Self::Velocity(v) => {
                Self::Velocity(v * rhs)
            },
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Div for SupportedQuantity {
    type Output = Result<Self, String>;

    fn div(
        self,
        rhs: Self,
    ) -> Self::Output {

        match (self, rhs) {
            | (
                Self::Length(l),
                Self::Time(t),
            ) => {
                Ok(Self::Velocity(
                    l / t,
                ))
            },
            | (
                Self::Length(_l),
                Self::Length(_l2),
            ) => Err(
                "Division resulting \
                 in dimensionless \
                 scalar is not yet \
                 supported"
                    .to_string(),
            ),
            | _ => {
                Err("Unsupported or \
                     complex quantity \
                     combination for \
                     division"
                    .to_string())
            },
        }
    }
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Div<f64> for SupportedQuantity {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self::Output {

        match self {
            | Self::Length(l) => {
                Self::Length(l / rhs)
            },
            | Self::Mass(m) => {
                Self::Mass(m / rhs)
            },
            | Self::Time(t) => {
                Self::Time(t / rhs)
            },
            | Self::Area(a) => {
                Self::Area(a / rhs)
            },
            | Self::Velocity(v) => {
                Self::Velocity(v / rhs)
            },
        }
    }
}

/// A wrapper struct for `SupportedQuantity` to be used within `Expr`.
#[derive(
    Clone,
    Debug,
    PartialEq,
    serde::Serialize,
    serde::Deserialize,
)]

pub struct UnitQuantity(
    pub SupportedQuantity,
);

#[allow(
    clippy::arithmetic_side_effects
)]

impl Eq for UnitQuantity {
}

#[allow(
    clippy::arithmetic_side_effects
)]

impl Hash for UnitQuantity {
    fn hash<H: Hasher>(
        &self,
        state: &mut H,
    ) {

        match &self.0 {
            | SupportedQuantity::Length(l) => {

                (
                    "Length",
                    ordered_float::OrderedFloat(l.value),
                )
                    .hash(state);
            },
            | SupportedQuantity::Mass(m) => {

                (
                    "Mass",
                    ordered_float::OrderedFloat(m.value),
                )
                    .hash(state);
            },
            | SupportedQuantity::Time(t) => {

                (
                    "Time",
                    ordered_float::OrderedFloat(t.value),
                )
                    .hash(state);
            },
            | SupportedQuantity::Area(a) => {

                (
                    "Area",
                    ordered_float::OrderedFloat(a.value),
                )
                    .hash(state);
            },
            | SupportedQuantity::Velocity(v) => {

                (
                    "Velocity",
                    ordered_float::OrderedFloat(v.value),
                )
                    .hash(state);
            },
        }
    }
}

/// Parses a value and a unit string into a `SupportedQuantity` enum.
#[inline]

pub(crate) fn parse_quantity(
    value: f64,
    unit: &str,
) -> Result<SupportedQuantity, String> {

    let unit_lower =
        unit.to_lowercase();

    match unit_lower.as_str() {
        | "m" | "meter" => {
            Ok(
                SupportedQuantity::Length(Length::new::<
                    length::meter,
                >(
                    value
                )),
            )
        },
        | "cm" | "centimeter" => {
            Ok(
                SupportedQuantity::Length(Length::new::<
                    length::centimeter,
                >(
                    value
                )),
            )
        },
        | "kg" | "kilogram" => {
            Ok(
                SupportedQuantity::Mass(Mass::new::<
                    mass::kilogram,
                >(
                    value
                )),
            )
        },
        | "g" | "gram" => {
            Ok(
                SupportedQuantity::Mass(Mass::new::<
                    mass::gram,
                >(
                    value
                )),
            )
        },
        | "s" | "second" => {
            Ok(
                SupportedQuantity::Time(Time::new::<
                    time::second,
                >(
                    value
                )),
            )
        },
        | "min" | "minute" => {
            Ok(
                SupportedQuantity::Time(Time::new::<
                    time::minute,
                >(
                    value
                )),
            )
        },
        | "m2" | "sqm" => {
            Ok(
                SupportedQuantity::Area(Area::new::<
                    area::square_meter,
                >(
                    value
                )),
            )
        },
        | "m/s" | "mps" => {
            Ok(
                SupportedQuantity::Velocity(Velocity::new::<
                    velocity::meter_per_second,
                >(
                    value
                )),
            )
        },
        | _ => {
            Err(format!(
                "Unknown or unsupported unit: {unit}"
            ))
        },
    }
}

/// Converts a numeric `Expr` variant into an `f64`.
///
/// This function handles `Expr::Constant`, `Expr::BigInt`, `Expr::Rational`, and `Expr::Dag`
/// (by converting it to an AST expression first).
///
/// # Arguments
/// * `expr` - The expression to convert.
///
/// # Returns
/// A `Result` containing the `f64` value or an error string.
#[inline]

pub(crate) fn expr_to_f64(
    expr: &Expr
) -> Result<f64, String> {

    let expr_ast =
        if let Expr::Dag(node) = expr {

            node.to_expr()
                .map_err(|e| {

                    format!(
                    "DAG conversion \
                     error: {e}"
                )
                })?
        } else {

            expr.clone()
        };

    expr_ast
        .to_f64()
        .ok_or_else(|| {

            format!(
                "Expression cannot be \
                 converted to a \
                 numeric value: \
                 {expr:?}"
            )
        })
}

/// Unifies an expression containing quantities with units.
///
/// This function recursively traverses the expression tree, resolving `QuantityWithValue`
/// nodes into `Quantity` nodes with concrete `UnitQuantity` values. It then performs
/// arithmetic operations on these quantities, checking for dimensional consistency.
///
/// # Arguments
/// * `expr` - The expression to unify.
///
/// # Returns
/// A `Result` containing the unified `Expr` or an error string if unification fails
/// (e.g., due to incompatible units).
///
/// # Examples
/// ```
/// 
/// use std::sync::Arc;
///
/// use rssn::symbolic::core::Expr;
/// use rssn::symbolic::unit_unification::unify_expression;
///
/// // 10m / 2s
/// let dist = Expr::QuantityWithValue(
///     Arc::new(Expr::new_constant(
///         10.0,
///     )),
///     "m".to_string(),
/// );
///
/// let time = Expr::QuantityWithValue(
///     Arc::new(Expr::new_constant(
///         2.0,
///     )),
///     "s".to_string(),
/// );
///
/// let velocity_expr = Expr::new_div(dist, time);
///
/// let result = unify_expression(&velocity_expr).unwrap();
/// // Result is a Quantity representing 5 m/s
/// ```

/// # Panics
///
/// Panics if a `Dag` node cannot be converted to an `Expr`, which indicates an
/// internal inconsistency in the expression representation. This should ideally
/// not happen in a well-formed expression DAG.
///
/// # Errors
///
/// This function will return an error if:
/// - A numerical value cannot be extracted from a quantity expression.
/// - Units cannot be parsed from the provided strings.
/// - Addition, subtraction, multiplication, or division of quantities fails
///   due to incompatible dimensions.
pub fn unify_expression(
    expr: &Expr
) -> Result<Expr, String> {

    match expr {
        | Expr::Dag(node) => {
            unify_expression(
                &node
                    .to_expr()
                    .expect(
                        "Dag Unify",
                    ),
            )
        },
        | Expr::QuantityWithValue(
            val_expr,
            unit_str,
        ) => {

            let value =
                expr_to_f64(val_expr)?;

            let quantity =
                parse_quantity(
                    value,
                    unit_str,
                )?;

            Ok(Expr::Quantity(
                Arc::new(UnitQuantity(
                    quantity,
                )),
            ))
        },
        | Expr::Add(a, b) => {

            let unified_a =
                unify_expression(a)?;

            let unified_b =
                unify_expression(b)?;

            match (unified_a, unified_b)
            {
                | (
                    Expr::Quantity(qa),
                    Expr::Quantity(qb),
                ) => {

                    let result = qa
                        .0
                        .clone()
                        + qb.0.clone();

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(
                            result?,
                        )),
                    ))
                },
                | (a, b) => {
                    Ok(Expr::new_add(
                        a, b,
                    ))
                },
            }
        },
        | Expr::Sub(a, b) => {

            let unified_a =
                unify_expression(a)?;

            let unified_b =
                unify_expression(b)?;

            match (unified_a, unified_b)
            {
                | (
                    Expr::Quantity(qa),
                    Expr::Quantity(qb),
                ) => {

                    let result = qa
                        .0
                        .clone()
                        - qb.0.clone();

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(
                            result?,
                        )),
                    ))
                },
                | (a, b) => {
                    Ok(Expr::new_sub(
                        a, b,
                    ))
                },
            }
        },
        | Expr::Mul(a, b) => {

            let unified_a =
                unify_expression(a)?;

            let unified_b =
                unify_expression(b)?;

            match (unified_a, unified_b)
            {
                | (
                    Expr::Quantity(qa),
                    Expr::Quantity(qb),
                ) => {

                    let result = (qa
                        .0
                        .clone()
                        * qb.0
                            .clone())?;

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(result)),
                    ))
                },
                | (
                    Expr::Quantity(qa),
                    Expr::Constant(
                        scalar_f64,
                    ),
                )
                | (
                    Expr::Constant(
                        scalar_f64,
                    ),
                    Expr::Quantity(qa),
                ) => {

                    let result = qa
                        .0
                        .clone()
                        * scalar_f64;

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(result)),
                    ))
                },
                | (a, b) => {
                    Ok(Expr::new_mul(
                        a, b,
                    ))
                },
            }
        },
        | Expr::Div(a, b) => {

            let unified_a =
                unify_expression(a)?;

            let unified_b =
                unify_expression(b)?;

            match (unified_a, unified_b)
            {
                | (
                    Expr::Quantity(qa),
                    Expr::Quantity(qb),
                ) => {

                    let result = (qa
                        .0
                        .clone()
                        / qb.0
                            .clone())?;

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(result)),
                    ))
                },
                #[allow(clippy::arithmetic_side_effects)]
                | (
                    Expr::Quantity(qa),
                    Expr::Constant(
                        scalar_f64,
                    ),
                ) => {

                    if !scalar_f64
                        .is_normal()
                    {

                        return Err(format!(
                            "Error: Division scalar must be a non-zero, finite number. Received: \
                             {scalar_f64}"
                        ));
                    }

                    let result = qa
                        .0
                        .clone()
                        / scalar_f64;

                    Ok(Expr::Quantity(
                        Arc::new(UnitQuantity(result)),
                    ))
                },
                | (
                    Expr::Constant(_),
                    Expr::Quantity(_),
                ) => Err(
                    "Error: Scalar \
                     divided by \
                     Quantity (S / Q) \
                     is not yet \
                     supported as it \
                     requires new \
                     reciprocal \
                     dimension types \
                     (like Frequency \
                     or Reciprocal \
                     Length)."
                        .to_string(),
                ),
                | (a, b) => {
                    Ok(Expr::new_div(
                        a, b,
                    ))
                },
            }
        },
        | Expr::Neg(a) => {

            let unified_a =
                unify_expression(a)?;

            if let Expr::Quantity(qa) =
                unified_a
            {

                Ok(Expr::Quantity(
                    Arc::new(
                        UnitQuantity(
                            qa.0.clone(
                            )
                            .neg(),
                        ),
                    ),
                ))
            } else {

                Ok(Expr::new_neg(
                    unified_a,
                ))
            }
        },
        | _ => Ok(expr.clone()),
    }
}
