use crate::symbolic::core::Expr;
use std::fmt::{self, Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::ops::{Add, Div, Mul, Neg, Sub};
// NOTE: Added 'Area' and 'Velocity' to f64 imports for Mul/Div results.
use ordered_float;
use uom::si::f64::{Area, Length, Mass, Time, Velocity};
use uom::si::{area, length, mass, time, velocity}; // NOTE: Required for safe hashing of f64 values.

// 1. The enum to hold different, but specific, quantity types.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum SupportedQuantity {
    Length(Length),
    Mass(Mass),
    Time(Time),
    Area(Area), // New: Needed for Length * Length
    Velocity(Velocity), // New: Needed for Length / Time
                // Add other quantities like Force, Energy, etc. here as needed
}

// 2. Implement arithmetic operations for the enum.

impl Add for SupportedQuantity {
    type Output = Result<Self, String>;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Length(l1), Self::Length(l2)) => Ok(Self::Length(l1 + l2)),
            (Self::Mass(m1), Self::Mass(m2)) => Ok(Self::Mass(m1 + m2)),
            (Self::Time(t1), Self::Time(t2)) => Ok(Self::Time(t1 + t2)),
            // Note: Area and Velocity are only compatible with themselves
            (Self::Area(a1), Self::Area(a2)) => Ok(Self::Area(a1 + a2)),
            (Self::Velocity(v1), Self::Velocity(v2)) => Ok(Self::Velocity(v1 + v2)),
            _ => Err("Incompatible types for addition".to_string()),
        }
    }
}

impl Sub for SupportedQuantity {
    type Output = Result<Self, String>;
    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Length(l1), Self::Length(l2)) => Ok(Self::Length(l1 - l2)),
            (Self::Mass(m1), Self::Mass(m2)) => Ok(Self::Mass(m1 - m2)),
            (Self::Time(t1), Self::Time(t2)) => Ok(Self::Time(t1 - t2)),
            (Self::Area(a1), Self::Area(a2)) => Ok(Self::Area(a1 - a2)),
            (Self::Velocity(v1), Self::Velocity(v2)) => Ok(Self::Velocity(v1 - v2)),
            _ => Err("Incompatible types for subtraction".to_string()),
        }
    }
}

// Implementation for Negation
impl Neg for SupportedQuantity {
    type Output = Self; // Negation does not change the type, so it doesn't need to return Result
    fn neg(self) -> Self::Output {
        match self {
            Self::Length(l) => Self::Length(-l),
            Self::Mass(m) => Self::Mass(-m),
            Self::Time(t) => Self::Time(-t),
            Self::Area(a) => Self::Area(-a),
            Self::Velocity(v) => Self::Velocity(-v),
        }
    }
}

// Implementation for Multiplication (Key Dimensional Changes)
impl Mul for SupportedQuantity {
    type Output = Result<Self, String>;
    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            // Length * Length = Area
            (Self::Length(l1), Self::Length(l2)) => Ok(Self::Area(l1 * l2)),
            // NOTE: Only simple, dimensionally clean products are supported for now.
            _ => Err("Unsupported or complex quantity combination for multiplication".to_string()),
        }
    }
}

// NEW: Implementation for SupportedQuantity * f64 (Scalar Multiplication)
impl Mul<f64> for SupportedQuantity {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        match self {
            Self::Length(l) => Self::Length(l * rhs),
            Self::Mass(m) => Self::Mass(m * rhs),
            Self::Time(t) => Self::Time(t * rhs),
            Self::Area(a) => Self::Area(a * rhs),
            Self::Velocity(v) => Self::Velocity(v * rhs),
        }
    }
}

// Implementation for Division (Key Dimensional Changes)
impl Div for SupportedQuantity {
    type Output = Result<Self, String>;
    fn div(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            // Length / Time = Velocity
            (Self::Length(l), Self::Time(t)) => Ok(Self::Velocity(l / t)),

            // Length / Length = dimensionless (Requires Scalar/Dimensionless variant)
            (Self::Length(_l), Self::Length(_l2)) => {
                Err("Division resulting in dimensionless scalar is not yet supported".to_string())
            }

            _ => Err("Unsupported or complex quantity combination for division".to_string()),
        }
    }
}

// NEW: Implementation for SupportedQuantity / f64 (Scalar Division)
impl Div<f64> for SupportedQuantity {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        match self {
            Self::Length(l) => Self::Length(l / rhs),
            Self::Mass(m) => Self::Mass(m / rhs),
            Self::Time(t) => Self::Time(t / rhs),
            Self::Area(a) => Self::Area(a / rhs),
            Self::Velocity(v) => Self::Velocity(v / rhs),
        }
    }
}

// 3. The wrapper to be stored in the Expr enum.
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct UnitQuantity(pub SupportedQuantity);

impl Eq for UnitQuantity {}

impl Hash for UnitQuantity {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // NOTE: Uses OrderedFloat for reliable hashing of f64 values.
        match &self.0 {
            SupportedQuantity::Length(l) => {
                ("Length", ordered_float::OrderedFloat(l.value)).hash(state)
            }
            SupportedQuantity::Mass(m) => {
                ("Mass", ordered_float::OrderedFloat(m.value)).hash(state)
            }
            SupportedQuantity::Time(t) => {
                ("Time", ordered_float::OrderedFloat(t.value)).hash(state)
            }
            SupportedQuantity::Area(a) => {
                ("Area", ordered_float::OrderedFloat(a.value)).hash(state)
            }
            SupportedQuantity::Velocity(v) => {
                ("Velocity", ordered_float::OrderedFloat(v.value)).hash(state)
            }
        }
    }
}

/// Parses a value and a unit string into a SupportedQuantity enum.
fn parse_quantity(value: f64, unit: &str) -> Result<SupportedQuantity, String> {
    let unit_lower = unit.to_lowercase();
    match unit_lower.as_str() {
        "m" | "meter" => Ok(SupportedQuantity::Length(Length::new::<length::meter>(
            value,
        ))),
        "cm" | "centimeter" => Ok(SupportedQuantity::Length(
            Length::new::<length::centimeter>(value),
        )),
        "kg" | "kilogram" => Ok(SupportedQuantity::Mass(Mass::new::<mass::kilogram>(value))),
        "g" | "gram" => Ok(SupportedQuantity::Mass(Mass::new::<mass::gram>(value))),
        "s" | "second" => Ok(SupportedQuantity::Time(Time::new::<time::second>(value))),
        "min" | "minute" => Ok(SupportedQuantity::Time(Time::new::<time::minute>(value))),
        // New unit parsing for added types (optional, as they can be results of operations)
        "m2" | "sqm" => Ok(SupportedQuantity::Area(Area::new::<area::square_meter>(
            value,
        ))),
        "m/s" | "mps" => Ok(SupportedQuantity::Velocity(Velocity::new::<
            velocity::meter_per_second,
        >(value))),
        _ => Err(format!("Unknown or unsupported unit: {}", unit)),
    }
}

/// Converts a numeric `Expr` variant into an `f64`.
/// NOTE: Assumes the `Expr` struct (not shown) has a method `to_f64()`.
fn expr_to_f64(expr: &Expr) -> Result<f64, String> {
    expr.to_f64().ok_or_else(|| {
        format!(
            "Expression cannot be converted to a numeric value: {:?}",
            expr
        )
    })
}

/// The main unification function.
pub fn unify_expression(expr: &Expr) -> Result<Expr, String> {
    match expr {
        Expr::QuantityWithValue(val_expr, unit_str) => {
            let value = expr_to_f64(val_expr)?;
            let quantity = parse_quantity(value, unit_str)?;
            Ok(Expr::Quantity(Box::new(UnitQuantity(quantity))))
        }
        Expr::Add(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;

            // FIX: Use 'match' to correctly handle the move.
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = (*qa).0 + (*qb).0;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result?))))
                }
                (a, b) => Ok(Expr::Add(Box::new(a), Box::new(b))),
            }
        }
        Expr::Sub(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;

            // FIX: Use 'match' to correctly handle the move.
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = (*qa).0 - (*qb).0;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result?))))
                }
                (a, b) => Ok(Expr::Sub(Box::new(a), Box::new(b))),
            }
        }
        Expr::Mul(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;

            // NEW: Multiplication logic
            match (unified_a, unified_b) {
                // Case 1: Quantity * Quantity (Q * Q)
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = ((*qa).0 * (*qb).0)?;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result))))
                }

                // Case 2: Quantity * Scalar (Q * S or S * Q)
                (Expr::Quantity(qa), Expr::Constant(scalar_f64))
                | (Expr::Constant(scalar_f64), Expr::Quantity(qa)) => {
                    // FIX: Removed '*' - scalar_f64 is already f64
                    let result = qa.0 * scalar_f64;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result))))
                }

                // Case 3: Re-wrap other combinations
                (a, b) => Ok(Expr::Mul(Box::new(a), Box::new(b))),
            }
        }
        Expr::Div(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;

            // NEW: Division logic
            match (unified_a, unified_b) {
                // Case 1: Quantity / Quantity (Q / Q)
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = ((*qa).0 / (*qb).0)?;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result))))
                }

                // Case 2: Quantity / Scalar (Q / S)
                (Expr::Quantity(qa), Expr::Constant(scalar_f64)) => {
                    // FIX: Removed '*' - scalar_f64 is already f64
                    let result = qa.0 / scalar_f64;
                    Ok(Expr::Quantity(Box::new(UnitQuantity(result))))
                }

                // Case 3: Scalar / Quantity (S / Q) - Resulting in reciprocal units
                (Expr::Constant(_), Expr::Quantity(_)) => {
                    Err("Error: Scalar divided by Quantity (S / Q) is not yet supported as it requires new reciprocal dimension types (like Frequency or Reciprocal Length).".to_string())
                }

                // Case 4: Re-wrap other combinations
                (a, b) => Ok(Expr::Div(Box::new(a), Box::new(b))),
            }
        }
        Expr::Neg(a) => {
            let unified_a = unify_expression(a)?;

            // NEW: Negation logic
            if let Expr::Quantity(qa) = unified_a {
                // Negation on quantity is straightforward
                Ok(Expr::Quantity(Box::new(UnitQuantity(qa.0.neg()))))
            } else {
                Ok(Expr::Neg(Box::new(unified_a)))
            }
        }
        // Pass through other expressions
        _ => Ok(expr.clone()),
    }
}
