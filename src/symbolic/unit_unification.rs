use crate::symbolic::core::Expr;
use ordered_float;
use std::fmt::Debug;
use std::hash::{Hash, Hasher};
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::sync::Arc;
use uom::si::f64::{Area, Length, Mass, Time, Velocity};
use uom::si::{area, length, mass, time, velocity};
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum SupportedQuantity {
    Length(Length),
    Mass(Mass),
    Time(Time),
    Area(Area),
    Velocity(Velocity),
}
#[allow(clippy::arithmetic_side_effects)]
impl Add for SupportedQuantity {
    type Output = Result<Self, String>;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Length(l1), Self::Length(l2)) => Ok(Self::Length(l1 + l2)),
            (Self::Mass(m1), Self::Mass(m2)) => Ok(Self::Mass(m1 + m2)),
            (Self::Time(t1), Self::Time(t2)) => Ok(Self::Time(t1 + t2)),
            (Self::Area(a1), Self::Area(a2)) => Ok(Self::Area(a1 + a2)),
            (Self::Velocity(v1), Self::Velocity(v2)) => Ok(Self::Velocity(v1 + v2)),
            _ => Err("Incompatible types for addition".to_string()),
        }
    }
}
#[allow(clippy::arithmetic_side_effects)]
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
#[allow(clippy::arithmetic_side_effects)]
impl Neg for SupportedQuantity {
    type Output = Self;
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
#[allow(clippy::arithmetic_side_effects)]
impl Mul for SupportedQuantity {
    type Output = Result<Self, String>;
    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Length(l1), Self::Length(l2)) => Ok(Self::Area(l1 * l2)),
            _ => Err("Unsupported or complex quantity combination for multiplication".to_string()),
        }
    }
}
#[allow(clippy::arithmetic_side_effects)]
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
#[allow(clippy::arithmetic_side_effects)]
impl Div for SupportedQuantity {
    type Output = Result<Self, String>;
    fn div(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Length(l), Self::Time(t)) => Ok(Self::Velocity(l / t)),
            (Self::Length(_l), Self::Length(_l2)) => {
                Err("Division resulting in dimensionless scalar is not yet supported".to_string())
            }
            _ => Err("Unsupported or complex quantity combination for division".to_string()),
        }
    }
}
#[allow(clippy::arithmetic_side_effects)]
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
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct UnitQuantity(pub SupportedQuantity);
#[allow(clippy::arithmetic_side_effects)]
impl Eq for UnitQuantity {}
#[allow(clippy::arithmetic_side_effects)]
impl Hash for UnitQuantity {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match &self.0 {
            SupportedQuantity::Length(l) => {
                ("Length", ordered_float::OrderedFloat(l.value)).hash(state);
            }
            SupportedQuantity::Mass(m) => {
                ("Mass", ordered_float::OrderedFloat(m.value)).hash(state);
            }
            SupportedQuantity::Time(t) => {
                ("Time", ordered_float::OrderedFloat(t.value)).hash(state);
            }
            SupportedQuantity::Area(a) => {
                ("Area", ordered_float::OrderedFloat(a.value)).hash(state);
            }
            SupportedQuantity::Velocity(v) => {
                ("Velocity", ordered_float::OrderedFloat(v.value)).hash(state);
            }
        }
    }
}
/// Parses a value and a unit string into a SupportedQuantity enum.
#[inline]
pub(crate) fn parse_quantity(value: f64, unit: &str) -> Result<SupportedQuantity, String> {
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
#[inline]
pub(crate) fn expr_to_f64(expr: &Expr) -> Result<f64, String> {
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
            Ok(Expr::Quantity(Arc::new(UnitQuantity(quantity))))
        }
        Expr::Add(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = qa.0.clone() + qb.0.clone();
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result?))))
                }
                (a, b) => Ok(Expr::new_add(a, b)),
            }
        }
        Expr::Sub(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = qa.0.clone() - qb.0.clone();
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result?))))
                }
                (a, b) => Ok(Expr::new_sub(a, b)),
            }
        }
        Expr::Mul(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = (qa.0.clone() * qb.0.clone())?;
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result))))
                }
                (Expr::Quantity(qa), Expr::Constant(scalar_f64))
                | (Expr::Constant(scalar_f64), Expr::Quantity(qa)) => {
                    let result = qa.0.clone() * scalar_f64;
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result))))
                }
                (a, b) => Ok(Expr::new_mul(a, b)),
            }
        }
        Expr::Div(a, b) => {
            let unified_a = unify_expression(a)?;
            let unified_b = unify_expression(b)?;
            match (unified_a, unified_b) {
                (Expr::Quantity(qa), Expr::Quantity(qb)) => {
                    let result = (qa.0.clone() / qb.0.clone())?;
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result))))
                }
                #[allow(clippy::arithmetic_side_effects)]
                (Expr::Quantity(qa), Expr::Constant(scalar_f64)) => {
                    if !scalar_f64.is_normal() {
                        return Err(
                            format!(
                                "Error: Division scalar must be a non-zero, finite number. Received: {}",
                                scalar_f64
                            ),
                        );
                    }
                    let result = qa.0.clone() / scalar_f64;
                    Ok(Expr::Quantity(Arc::new(UnitQuantity(result))))
                }
                (Expr::Constant(_), Expr::Quantity(_)) => {
                    Err(
                        "Error: Scalar divided by Quantity (S / Q) is not yet supported as it requires new reciprocal dimension types (like Frequency or Reciprocal Length)."
                            .to_string(),
                    )
                }
                (a, b) => Ok(Expr::new_div(a, b)),
            }
        }
        Expr::Neg(a) => {
            let unified_a = unify_expression(a)?;
            if let Expr::Quantity(qa) = unified_a {
                Ok(Expr::Quantity(Arc::new(UnitQuantity(qa.0.clone().neg()))))
            } else {
                Ok(Expr::new_neg(unified_a))
            }
        }
        _ => Ok(expr.clone()),
    }
}
