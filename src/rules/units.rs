//! Units of measurement and dimensional analysis.
//!
//! A physical quantity is the term `quantity(value, unit)`. The value is
//! any scalar term (exact, float or symbolic); the unit is a product of
//! unit symbols and rational powers of them: `m`, `km/h`, `kg*m/s^2`,
//! `m^2`, `1` (dimensionless). Unit symbols are case sensitive, as in the
//! SI: `m` is the metre, `s` the second, `N` the newton.
//!
//! | operator | value |
//! |---|---|
//! | `quantity(v, u)` | the quantity itself (an inert term) |
//! | `unify_expression(e)` | evaluates an expression of quantities: sums need equal dimensions and are expressed in the first operand's unit, products and powers combine units, a dimensionless result becomes a plain number; a **dimension mismatch leaves the request unreduced** |
//! | `convert(q, u)` | the quantity `q` expressed in the unit `u` (same dimension required) |
//! | `dimension(q)` | the dimension as a product of `Length`, `Mass`, `Time`, `Current`, `Temperature`, `Amount`, `Luminosity` (`1` if dimensionless) |
//! | `simplify_units(e)` | like `unify_expression`, with the result in coherent SI units: a named derived unit (`N`, `Pa`, `J`, `W`, `Hz`, `C`, `V`, `F`, `ohm`, `S`, `Wb`, `T`, `H`) when the dimension has one, otherwise a product of `m`, `kg`, `s`, `A`, `K`, `mol`, `cd` |
//!
//! # Units
//!
//! * SI base units `m`, `g` (so `kg` is `k` + `g`), `s`, `A`, `K`, `mol`,
//!   `cd`; derived `rad`, `sr`, `Hz`, `N`, `Pa`, `J`, `W`, `C`, `V`, `F`,
//!   `ohm`, `S`, `Wb`, `T`, `H`, `lm`, `lx`, `Bq`, `Gy`, `Sv`, `kat`.
//! * Every prefix from `Y` (10^24) to `y` (10^-24), including `da`, `u` or
//!   `µ` for micro, on the SI units and on `L`, `l`, `t`, `bar`, `eV`, `cal`,
//!   `Wh`, `rad`.
//! * Other units: `min`, `h`, `d`, `week`, `yr`, `deg`, `arcmin`, `arcsec`,
//!   `in`, `ft`, `yd`, `mi`, `nmi`, `angstrom`, `lb`, `oz`, `L`, `gal`,
//!   `atm`, `mmHg`, `psi`, `lbf`, `hp`, and the affine temperature units
//!   `degC` and `degF`.
//! * Long names: `meter`, `metre`, `centimeter`, `kilometer`, `kilogram`,
//!   `gram`, `second`, `minute`, `hour`, `day`, `liter`, `litre`, `kelvin`,
//!   `celsius`, `fahrenheit`, and the legacy abbreviations `sqm` (`m^2`) and
//!   `mps` (`m/s`).
//!
//! All scale factors are exact rationals (degrees carry an exact factor of
//! `pi`), so exact values convert exactly: `convert(quantity(5, km), m)` is
//! `quantity(5000, m)`. `degC` and `degF` are affine: they can only be used
//! alone, a lone `quantity(t, degC)` is converted through kelvin when it
//! takes part in arithmetic, and `convert` goes between them exactly.
//! Values with no known unit symbol (or an unknown one inside a unit) leave
//! every request unreduced.

use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::Arity;
use crate::graph::Cx;
use crate::graph::Kernel;
use crate::graph::NodeId;
use crate::graph::Number;
use crate::graph::OpDescriptor;
use crate::graph::OpFlags;
use crate::graph::OpId;
use crate::graph::Outcome;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::op::core;
use crate::graph::rule::Installer;
use crate::rules::arith::arith;
use crate::rules::elementary::elementary;
use crate::rules::poly::best;

use crate::dim::units::DIMENSIONS;
use crate::dim::units::Dim;
use crate::dim::units::Exp;
use crate::dim::units::Scale;
use crate::dim::units::UnitMap;
use crate::dim::units::clean;
use crate::dim::units::parse_unit;
use crate::dim::units::resolve_map;
use crate::dim::units::resolve_unit_node;
use crate::dim::units::resolve_unit_text;

/// The units rule set.
#[must_use]
pub fn units() -> RuleSet {
    RuleSet::new("units", install)
        .needs(arith())
        .needs(elementary())
}

/// A scalar or a quantity with a unit.
#[derive(Clone, Debug)]
enum Val {
    Scalar(NodeId),
    Quant { value: NodeId, unit: UnitMap },
}

struct Units<'c, 'a> {
    cx: &'c mut Cx<'a>,
    quantity: OpId,
}

impl Units<'_, '_> {
    fn scale_term(
        &mut self,
        scale: &Scale,
    ) -> NodeId {
        let rat = self.cx.graph.num(Number::rat(scale.rat.clone()));
        if scale.pi == 0 {
            return rat;
        }
        let pi = self
            .cx
            .graph
            .ops()
            .lookup("pi")
            .map(|p| self.cx.graph.node(p, &[]));
        match pi {
            | Some(pi) => {
                let e = self.cx.graph.int(i64::from(scale.pi));
                let power = self.cx.graph.node(core::POW, &[pi, e]);
                self.cx.graph.node(core::MUL, &[rat, power])
            },
            | None => rat,
        }
    }

    /// `value * scale`, simplified.
    fn scaled(
        &mut self,
        value: NodeId,
        scale: &Scale,
    ) -> NodeId {
        if scale.is_one() {
            return value;
        }
        let factor = self.scale_term(scale);
        let product = self.cx.graph.node(core::MUL, &[value, factor]);
        self.cx.simplify(product)
    }

    fn unit_term(
        &mut self,
        unit: &UnitMap,
    ) -> NodeId {
        if unit.is_empty() {
            return self.cx.graph.int(1);
        }
        let mut factors = Vec::new();
        for (name, &e) in unit {
            let symbol = name.term(self.cx.graph);
            factors.push(if e.is_one() {
                symbol
            } else {
                let exponent = self.cx.graph.num(
                    Number::fraction(i64::from(*e.numer()), i64::from(*e.denom()))
                        .unwrap_or_else(|| Number::from(1)),
                );
                self.cx.graph.node(core::POW, &[symbol, exponent])
            });
        }
        match factors.as_slice() {
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, &factors),
        }
    }

    fn quantity_term(
        &mut self,
        value: NodeId,
        unit: &UnitMap,
    ) -> NodeId {
        let unit = self.unit_term(unit);
        let quantity = self.quantity;
        self.cx.graph.node(quantity, &[value, unit])
    }

    fn rational(
        &self,
        node: NodeId,
    ) -> Option<BigRational> {
        self.cx.graph.number_of(node).and_then(Number::to_rational)
    }

    /// Reads `quantity(v, u)` without evaluating anything else.
    fn lone_quantity(
        &mut self,
        node: NodeId,
    ) -> Option<Val> {
        let term = best(self.cx.graph, node)?;
        if self.cx.graph.op(term) != self.quantity {
            return None;
        }
        let &[value, unit] = self.cx.graph.children(term) else {
            return None;
        };
        let unit = parse_unit(self.cx.graph, unit, 0)?;
        let resolved = resolve_map(&unit)?;
        let Val::Scalar(value) = self.eval(value, 0)? else {
            return None;
        };
        // Affine units are converted through kelvin.
        if let Some(offset) = resolved.offset {
            let scaled = self.scaled(value, &resolved.scale);
            let shift = self.cx.graph.num(Number::rat(offset));
            let sum = self.cx.graph.node(core::ADD, &[scaled, shift]);
            let sum = self.cx.simplify(sum);
            let kelvin = coherent_unit(&resolved.dim);
            return Some(Val::Quant {
                value: sum,
                unit: kelvin,
            });
        }
        Some(Val::Quant { value, unit })
    }

    /// Whether `node` is a lone affine quantity (an absolute temperature
    /// point such as `quantity(20, degC)`).
    ///
    /// Affine points have no multiplicative structure and two of them cannot
    /// be added: `unify_expression` refuses such requests instead of
    /// silently reading them as kelvin. A point plus a difference written in
    /// `K` is still allowed (the difference is the vector unit in SI).
    fn is_affine_point(
        &mut self,
        node: NodeId,
    ) -> bool {
        let Some(term) = best(self.cx.graph, node) else {
            return false;
        };
        if self.cx.graph.op(term) != self.quantity {
            return false;
        }
        let &[_, unit] = self.cx.graph.children(term) else {
            return false;
        };
        let Some(unit) = parse_unit(self.cx.graph, unit, 0) else {
            return false;
        };
        resolve_map(&unit).is_some_and(|resolved| resolved.offset.is_some())
    }

    /// A dimensionless quantity as the number it stands for.
    fn demote(
        &mut self,
        val: Val,
    ) -> Option<Val> {
        match val {
            | Val::Quant { value, unit } => self.collapse(value, unit),
            | scalar @ Val::Scalar(_) => Some(scalar),
        }
    }

    /// A product of unit symbols as a scalar when it is dimensionless.
    fn collapse(
        &mut self,
        value: NodeId,
        unit: UnitMap,
    ) -> Option<Val> {
        let resolved = resolve_map(&unit)?;
        if unit.is_empty() || resolved.dim.is_zero() {
            return Some(Val::Scalar(self.scaled(value, &resolved.scale)));
        }
        Some(Val::Quant { value, unit })
    }

    #[allow(clippy::too_many_lines)]
    fn eval(
        &mut self,
        node: NodeId,
        depth: usize,
    ) -> Option<Val> {
        if depth > 48 {
            return None;
        }
        if let Some(q) = self.lone_quantity(node) {
            return Some(q);
        }
        let term = best(self.cx.graph, node)?;
        let op = self.cx.graph.op(term);
        let children = self.cx.graph.children(term).to_vec();
        if children.is_empty() {
            return Some(Val::Scalar(term));
        }
        if op == self.quantity {
            // A malformed quantity or one with an unknown unit.
            return None;
        }
        match op {
            | core::ADD => {
                // point + point (two absolute temperatures) is not an
                // absolute temperature: refuse it rather than convert both
                // through kelvin and add.
                if children
                    .iter()
                    .filter(|&&c| self.is_affine_point(c))
                    .count()
                    >= 2
                {
                    return None;
                }
                let vals: Vec<Val> = children
                    .iter()
                    .map(|&c| self.eval(c, depth + 1).and_then(|v| self.demote(v)))
                    .collect::<Option<_>>()?;
                if vals.iter().all(|v| matches!(v, Val::Scalar(_))) {
                    let nodes: Vec<NodeId> = vals
                        .iter()
                        .map(|v| {
                            if let Val::Scalar(n) = v {
                                *n
                            } else {
                                term
                            }
                        })
                        .collect();
                    return Some(Val::Scalar(self.cx.graph.node(core::ADD, &nodes)));
                }
                let Some(Val::Quant { unit: first_unit, .. }) = vals
                    .iter()
                    .find(|v| matches!(v, Val::Quant { .. }))
                    .cloned()
                else {
                    return None;
                };
                let first = resolve_map(&first_unit)?;
                let mut terms = Vec::new();
                for v in vals {
                    let Val::Quant { value, unit } = v else {
                        return None; // a scalar added to a quantity
                    };
                    let r = resolve_map(&unit)?;
                    if r.dim != first.dim {
                        return None;
                    }
                    terms.push(self.scaled(value, &r.scale.over(&first.scale)));
                }
                let sum = self.cx.graph.node(core::ADD, &terms);
                Some(Val::Quant {
                    value: self.cx.simplify(sum),
                    unit: first_unit,
                })
            },
            | core::MUL => {
                // A point is affine: `degC * 2` is not a temperature.
                if children.len() > 1 && children.iter().any(|&c| self.is_affine_point(c)) {
                    return None;
                }
                let mut values = Vec::new();
                let mut unit = UnitMap::new();
                let mut any = false;
                for c in children {
                    match self.eval(c, depth + 1)? {
                        | Val::Scalar(n) => values.push(n),
                        | Val::Quant { value, unit: u } => {
                            any = true;
                            values.push(value);
                            for (name, k) in u {
                                *unit.entry(name).or_insert_with(Exp::zero) += k;
                            }
                        },
                    }
                }
                let value = self.cx.graph.node(core::MUL, &values);
                if !any {
                    return Some(Val::Scalar(value));
                }
                self.collapse(value, clean(unit))
            },
            | core::POW => {
                let &[base, exponent] = children.as_slice() else {
                    return None;
                };
                // A power of an affine point is not a temperature.
                if self.is_affine_point(base) {
                    return None;
                }
                let (b, e) = (self.eval(base, depth + 1)?, self.eval(exponent, depth + 1)?);
                let e = self.demote(e)?;
                match (b, e) {
                    | (Val::Scalar(b), Val::Scalar(e)) => {
                        Some(Val::Scalar(self.cx.graph.node(core::POW, &[b, e])))
                    },
                    | (Val::Quant { value, unit }, Val::Scalar(e)) => {
                        let r = self.rational(e)?;
                        let k = Exp::new(r.numer().to_i32()?, r.denom().to_i32()?);
                        let value = self.cx.graph.node(core::POW, &[value, e]);
                        self.collapse(
                            value,
                            clean(unit.into_iter().map(|(n, x)| (n, x * k)).collect()),
                        )
                    },
                    | _ => None,
                }
            },
            | _ => {
                // A function of scalars.
                let mut args = Vec::with_capacity(children.len());
                for c in children {
                    let value = self.eval(c, depth + 1)?;
                    let Val::Scalar(n) = self.demote(value)? else {
                        return None;
                    };
                    args.push(n);
                }
                Some(Val::Scalar(self.cx.graph.try_node(op, &args)?))
            },
        }
    }

    fn finish(
        &mut self,
        val: &Val,
    ) -> NodeId {
        match val {
            | Val::Scalar(n) => self.cx.simplify(*n),
            | Val::Quant { value, unit } => {
                let value = self.cx.simplify(*value);
                self.quantity_term(value, unit)
            },
        }
    }

    /// The value in coherent SI units.
    fn coherent(
        &mut self,
        val: &Val,
    ) -> Option<NodeId> {
        let Val::Quant { value, unit } = val else {
            return Some(self.finish(val));
        };
        let r = resolve_map(unit)?;
        let value = self.scaled(*value, &r.scale);
        if r.dim.is_zero() {
            return Some(value);
        }
        let coherent = coherent_unit(&r.dim);
        Some(self.quantity_term(value, &coherent))
    }

    fn convert(
        &mut self,
        q: NodeId,
        target: NodeId,
    ) -> Option<NodeId> {
        let val = self.eval(q, 0)?;
        let Val::Quant { value, unit } = val else {
            return None;
        };
        let target_unit = parse_unit(self.cx.graph, target, 0)?;
        let source = resolve_map(&unit)?;
        let destination = resolve_map(&target_unit)?;
        if source.dim != destination.dim {
            return None;
        }
        let si = self.scaled(value, &source.scale);
        let shifted = match &destination.offset {
            | Some(offset) => {
                let shift = self.cx.graph.num(Number::rat(-offset.clone()));
                let sum = self.cx.graph.node(core::ADD, &[si, shift]);
                self.cx.simplify(sum)
            },
            | None => si,
        };
        let out = self.scaled(shifted, &Scale::one().over(&destination.scale));
        Some(self.quantity_term(out, &target_unit))
    }

    fn dimension(
        &mut self,
        q: NodeId,
    ) -> Option<NodeId> {
        let val = self.eval(q, 0)?;
        let dim = match val {
            | Val::Scalar(_) => Dim::zero(),
            | Val::Quant { unit, .. } => resolve_map(&unit)?.dim,
        };
        let mut factors = Vec::new();
        for (name, &e) in DIMENSIONS
            .iter()
            .copied()
            .zip(dim.si())
            .chain(dim.extra().iter().map(|(name, e)| (name.as_str(), e)))
        {
            if e.is_zero() {
                continue;
            }
            let symbol = self.cx.graph.sym(name);
            factors.push(if e.is_one() {
                symbol
            } else {
                let exponent = self.cx.graph.num(
                    Number::fraction(i64::from(*e.numer()), i64::from(*e.denom()))
                        .unwrap_or_else(|| Number::from(1)),
                );
                self.cx.graph.node(core::POW, &[symbol, exponent])
            });
        }
        let product = match factors.as_slice() {
            | [] => self.cx.graph.int(1),
            | [only] => *only,
            | _ => self.cx.graph.node(core::MUL, &factors),
        };
        Some(product)
    }
}

/// The coherent SI unit of a dimension: a named derived unit when there
/// is one, else the product of base units.
fn coherent_unit(dim: &Dim) -> UnitMap {
    const NAMED: [(&str, [i32; 7]); 13] = [
        ("N", [1, 1, -2, 0, 0, 0, 0]),
        ("Pa", [-1, 1, -2, 0, 0, 0, 0]),
        ("J", [2, 1, -2, 0, 0, 0, 0]),
        ("W", [2, 1, -3, 0, 0, 0, 0]),
        ("Hz", [0, 0, -1, 0, 0, 0, 0]),
        ("C", [0, 0, 1, 1, 0, 0, 0]),
        ("V", [2, 1, -3, -1, 0, 0, 0]),
        ("F", [-2, -1, 4, 2, 0, 0, 0]),
        ("ohm", [2, 1, -3, -2, 0, 0, 0]),
        ("S", [-2, -1, 3, 2, 0, 0, 0]),
        ("Wb", [2, 1, -2, -1, 0, 0, 0]),
        ("T", [0, 1, -2, -1, 0, 0, 0]),
        ("H", [2, 1, -2, -2, 0, 0, 0]),
    ];
    let mut map = UnitMap::new();
    if !dim.extra().is_empty() {
        let unit = crate::dim::UnitDef::linear(
            "coherent",
            crate::dim::units::dim_to_expr(dim),
            Scale::one(),
        )
        .to_value();
        map.insert(
            crate::dim::units::UnitAtom::Literal(std::sync::Arc::new(unit)),
            Exp::one(),
        );
        return map;
    }
    if let Some((name, _)) = NAMED.iter().find(|(_, d)| Dim::of(*d) == *dim) {
        map.insert((*name).to_owned().into(), Exp::one());
        return map;
    }
    for (name, &e) in ["m", "kg", "s", "A", "K", "mol", "cd"].iter().zip(dim.si()) {
        if !e.is_zero() {
            map.insert((*name).to_owned().into(), e);
        }
    }
    map
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Request {
    Unify,
    Convert,
    Dimension,
    Simplify,
}

struct Request1 {
    op: OpId,
    request: Request,
    quantity: OpId,
}

impl Kernel for Request1 {
    fn ops(&self) -> Vec<OpId> {
        vec![self.op]
    }

    fn reduce(
        &self,
        cx: &mut Cx<'_>,
        node: NodeId,
    ) -> Outcome {
        let args = cx.graph.children(node).to_vec();
        let mut units = Units {
            cx,
            quantity: self.quantity,
        };
        let result = match (self.request, args.as_slice()) {
            | (Request::Unify, &[e]) => units.eval(e, 0).map(|v| units.finish(&v)),
            | (Request::Simplify, &[e]) => units.eval(e, 0).and_then(|v| units.coherent(&v)),
            | (Request::Convert, &[q, u]) => units.convert(q, u),
            | (Request::Dimension, &[q]) => units.dimension(q),
            | _ => None,
        };
        result.map_or(Outcome::Pass, Outcome::Equal)
    }
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let quantity = i.op(OpDescriptor::new("quantity", Arity::Fixed(2)))?;
    i.graph()
        .ops_mut()
        .set_attr(quantity, crate::graph::dim::UnitResolver(resolve_unit_node));
    i.graph().ops_mut().set_attr(
        quantity,
        crate::graph::dim::UnitTextResolver(resolve_unit_text),
    );
    for (name, arity, request) in [
        ("unify_expression", 1, Request::Unify),
        ("convert", 2, Request::Convert),
        ("dimension", 1, Request::Dimension),
        ("simplify_units", 1, Request::Simplify),
    ] {
        let op = i.op(OpDescriptor::new(name, Arity::Fixed(arity))
            .flags(OpFlags::HEAVY)
            .cost(100))?;
        i.kernel(
            &format!("units/{name}"),
            Tier::Reduce,
            Request1 {
                op,
                request,
                quantity,
            },
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use crate::dim::units::ALIASES;
    use crate::dim::units::PREFIXES;
    use crate::dim::units::dec;
    use crate::dim::units::resolve_symbol;
    use crate::dim::units::table;

    use super::*;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> String {
        simplify(&[units()], src)
    }

    /// The text and whether the request was reduced.
    fn attempt(src: &str) -> (String, bool) {
        reduce_with(&[units()], src, &[])
    }

    #[test]
    fn quantities_are_inert_terms() {
        assert_eq!(run("quantity(5, m)"), "quantity(5, m)");
        assert_eq!(run("quantity(x, kg*m/s^2)"), "quantity(x, kg*m/s^2)");
    }

    #[test]
    fn sums_convert_to_the_first_unit() {
        assert_eq!(
            run("unify_expression(quantity(5, m) + quantity(3, m))"),
            "quantity(8, m)"
        );
        assert_eq!(
            run("unify_expression(quantity(1, km) + quantity(500, m))"),
            "quantity(3/2, km)"
        );
        assert_eq!(
            run("unify_expression(quantity(500, m) + quantity(1, km))"),
            "quantity(1500, m)"
        );
        assert_eq!(
            run("unify_expression(quantity(2, h) - quantity(30, min))"),
            "quantity(3/2, h)"
        );
        assert_eq!(
            run("unify_expression(quantity(1, ft) + quantity(12, in))"),
            "quantity(2, ft)"
        );
    }

    #[test]
    fn mismatched_dimensions_do_not_reduce() {
        for src in [
            "unify_expression(quantity(5, m) + quantity(3, s))",
            "unify_expression(quantity(5, m) + 3)",
            "unify_expression(quantity(5, m) + x)",
            "unify_expression(quantity(1, kg) - quantity(1, N))",
            "convert(quantity(5, m), s)",
            "unify_expression(sin(quantity(5, m)))",
            "unify_expression(quantity(5, parsec))",
        ] {
            let (text, reduced) = attempt(src);
            assert!(!reduced, "{src} reduced to {text}");
        }
    }

    #[test]
    fn products_quotients_and_powers() {
        assert_eq!(
            run("unify_expression(quantity(10, m) / quantity(2, s))"),
            "quantity(5, m/s)"
        );
        assert_eq!(
            run("unify_expression(quantity(3, m) * quantity(4, m))"),
            "quantity(12, m^2)"
        );
        assert_eq!(
            run("unify_expression(quantity(3, m)^2)"),
            "quantity(9, m^2)"
        );
        assert_eq!(
            run("unify_expression(quantity(9, m^2)^(1/2))"),
            "quantity(3, m)"
        );
        assert_eq!(
            run("unify_expression(2 * quantity(3, kg))"),
            "quantity(6, kg)"
        );
        assert_eq!(
            run("unify_expression(quantity(6, kg) / 3)"),
            "quantity(2, kg)"
        );
        assert_eq!(
            run("unify_expression(1 / quantity(4, s))"),
            "quantity(1/4, 1/s)"
        );
        assert_eq!(
            run("unify_expression(x * quantity(2, m))"),
            "quantity(2*x, m)"
        );
        // force = mass * acceleration
        assert_eq!(
            run("unify_expression(quantity(2, kg) * quantity(3, m/s^2))"),
            "quantity(6, kg*m/s^2)"
        );
    }

    #[test]
    fn dimensionless_results_become_numbers() {
        assert_eq!(
            run("unify_expression(quantity(1, m) / quantity(1, cm))"),
            "100"
        );
        assert_eq!(
            run("unify_expression(quantity(6, m) / quantity(2, m))"),
            "3"
        );
        assert_eq!(
            run("unify_expression(quantity(3, s) * quantity(2, Hz))"),
            "6"
        );
        assert_eq!(
            run("unify_expression(sin(quantity(3, m) / quantity(3, m)))"),
            "sin(1)"
        );
        assert_eq!(
            run("unify_expression(quantity(5, m) + quantity(3, m) - quantity(8, m))"),
            "quantity(0, m)"
        );
    }

    #[test]
    fn conversion() {
        assert_eq!(run("convert(quantity(5, km), m)"), "quantity(5000, m)");
        assert_eq!(run("convert(quantity(250, cm), m)"), "quantity(5/2, m)");
        assert_eq!(run("convert(quantity(1, h), s)"), "quantity(3600, s)");
        assert_eq!(run("convert(quantity(36, km/h), m/s)"), "quantity(10, m/s)");
        assert_eq!(run("convert(quantity(1, mi), ft)"), "quantity(5280, ft)");
        assert_eq!(
            run("convert(quantity(1, lb), g)"),
            "quantity(45359237/100000, g)"
        );
        assert_eq!(
            run("convert(quantity(1, N), kg*m/s^2)"),
            "quantity(1, kg*m/s^2)"
        );
        assert_eq!(run("convert(quantity(1, kWh), J)"), "quantity(3600000, J)");
        assert_eq!(run("convert(quantity(1, L), cm^3)"), "quantity(1000, cm^3)");
        assert_eq!(run("convert(quantity(180, deg), rad)"), "quantity(pi, rad)");
        assert_eq!(run("convert(quantity(x, mm), m)"), "quantity(1/1000*x, m)");
        assert_eq!(
            run("convert(quantity(2, GHz), Hz)"),
            "quantity(2000000000, Hz)"
        );
        assert_eq!(
            run("convert(quantity(1, ug), kg)"),
            "quantity(1/1000000000, kg)"
        );
        assert_eq!(
            run("convert(quantity(1, kilogram), gram)"),
            "quantity(1000, gram)"
        );
        // expressions are evaluated first
        assert_eq!(
            run("convert(quantity(1, km) + quantity(500, m), m)"),
            "quantity(1500, m)"
        );
    }

    #[test]
    fn temperatures_are_affine() {
        assert_eq!(run("convert(quantity(0, degC), K)"), "quantity(5463/20, K)");
        assert_eq!(
            run("convert(quantity(100, degC), degF)"),
            "quantity(212, degF)"
        );
        assert_eq!(
            run("convert(quantity(32, degF), degC)"),
            "quantity(0, degC)"
        );
        assert_eq!(
            run("convert(quantity(20, degC), degC)"),
            "quantity(20, degC)"
        );
        assert_eq!(
            run("convert(quantity(300, K), degC)"),
            "quantity(537/20, degC)"
        );
        // they cannot be combined multiplicatively
        assert!(!attempt("unify_expression(quantity(1, degC*m))").1);
        // an affine point is not a vector: it cannot be scaled, powered or
        // added to another point. These used to be quietly converted to K.
        for source in [
            "unify_expression(quantity(20, degC) * 2)",
            "unify_expression(2 * quantity(20, degC))",
            "unify_expression(quantity(20, degC) / 2)",
            "unify_expression(quantity(20, degC)^2)",
            "unify_expression(quantity(20, degC) + quantity(30, degC))",
        ] {
            let (text, reduced) = attempt(source);
            assert!(!reduced, "{source} reduced to {text}");
        }
        // a point plus a difference expressed in kelvin is still allowed
        assert_eq!(
            run("unify_expression(quantity(20, degC) + quantity(2, K))"),
            "quantity(5903/20, K)"
        );
    }

    #[test]
    fn dimensions() {
        assert_eq!(run("dimension(quantity(5, m))"), "Length");
        assert_eq!(run("dimension(quantity(5, m/s))"), "Length/Time");
        assert_eq!(run("dimension(quantity(5, N))"), "Length*Mass/Time^2");
        assert_eq!(run("dimension(quantity(5, J))"), "Length^2*Mass/Time^2");
        assert_eq!(run("dimension(quantity(5, km^2))"), "Length^2");
        assert_eq!(run("dimension(quantity(5, mol/L))"), "Amount/Length^3");
        assert_eq!(run("dimension(quantity(5, rad))"), "1");
        assert_eq!(run("dimension(7)"), "1");
        assert_eq!(
            run("dimension(quantity(1, kg) * quantity(1, m) / quantity(1, s)^2)"),
            "Length*Mass/Time^2"
        );
        // equal dimensions: N m and J
        assert_eq!(
            run("dimension(quantity(1, N*m))"),
            run("dimension(quantity(1, J))")
        );
    }

    #[test]
    fn simplification_to_si() {
        assert_eq!(
            run("simplify_units(quantity(36, km/h))"),
            "quantity(10, m/s)"
        );
        assert_eq!(
            run("simplify_units(quantity(2, kg) * quantity(3, m/s^2))"),
            "quantity(6, N)"
        );
        assert_eq!(
            run("simplify_units(quantity(1, N) * quantity(2, m))"),
            "quantity(2, J)"
        );
        assert_eq!(
            run("simplify_units(quantity(1, J) / quantity(2, s))"),
            "quantity(1/2, W)"
        );
        assert_eq!(
            run("simplify_units(quantity(1, V) / quantity(2, A))"),
            "quantity(1/2, ohm)"
        );
        assert_eq!(
            run("simplify_units(quantity(1, cm) + quantity(1, m))"),
            "quantity(101/100, m)"
        );
        assert_eq!(
            run("simplify_units(quantity(90, deg) + quantity(1, m) / quantity(1, m))"),
            "1/2*pi + 1"
        );
        assert_eq!(
            run("simplify_units(quantity(1, kW) * quantity(1, h))"),
            "quantity(3600000, J)"
        );
        assert_eq!(
            run("simplify_units(quantity(4, s) * quantity(5, Hz))"),
            "20"
        );
        assert_eq!(run("simplify_units(quantity(5, m))"), "quantity(5, m)");
        assert_eq!(
            run("simplify_units(quantity(1, ft)^3)"),
            "quantity(55306341/1953125000, m^3)"
        );
    }

    #[test]
    fn legacy_units_are_accepted() {
        assert_eq!(
            run("unify_expression(quantity(5, m) + quantity(300, cm))"),
            "quantity(8, m)"
        );
        assert_eq!(
            run("unify_expression(quantity(1, kg) + quantity(500, g))"),
            "quantity(3/2, kg)"
        );
        assert_eq!(
            run("unify_expression(quantity(60, s) + quantity(1, min))"),
            "quantity(120, s)"
        );
        assert_eq!(
            run("convert(quantity(2, sqm), cm^2)"),
            "quantity(20000, cm^2)"
        );
        assert_eq!(
            run("convert(quantity(3, mps), km/h)"),
            "quantity(54/5, km/h)"
        );
        assert_eq!(
            run("unify_expression(quantity(10, meter) / quantity(2, second))"),
            "quantity(5, meter/second)"
        );
        assert_eq!(
            run("unify_expression(quantity(5, m) * quantity(5, m))"),
            "quantity(25, m^2)"
        );
    }

    #[test]
    fn every_unit_resolves_and_is_consistent() {
        // exact table entries
        for name in table().keys() {
            assert!(resolve_symbol(name).is_some(), "{name}");
        }
        // prefixes
        for (prefix, power) in PREFIXES {
            let r = resolve_symbol(&format!("{prefix}m")).unwrap_or_else(|| panic!("{prefix}m"));
            assert_eq!(
                r.scale.rat,
                BigRational::from_integer(BigInt::from(10)).pow(power),
                "{prefix}"
            );
            assert_eq!(r.dim, Dim::of([1, 0, 0, 0, 0, 0, 0]));
        }
        // kg is the base unit of mass
        let kg = resolve_symbol("kg").unwrap_or_else(|| panic!("kg"));
        assert!(kg.scale.is_one());
        // derived units are consistent with their definitions
        let check = |name: &str, product: &[(&str, i32)]| {
            let map: UnitMap = product
                .iter()
                .map(|&(n, e)| (n.to_owned().into(), Exp::from_integer(e)))
                .collect();
            let a = resolve_symbol(name).unwrap_or_else(|| panic!("{name}"));
            let b = resolve_map(&map).unwrap_or_else(|| panic!("{name}"));
            assert_eq!(a.dim, b.dim, "{name}");
            assert_eq!(a.scale.rat, b.scale.rat, "{name}");
        };
        check("N", &[("kg", 1), ("m", 1), ("s", -2)]);
        check("J", &[("N", 1), ("m", 1)]);
        check("W", &[("J", 1), ("s", -1)]);
        check("Pa", &[("N", 1), ("m", -2)]);
        check("V", &[("W", 1), ("A", -1)]);
        check("ohm", &[("V", 1), ("A", -1)]);
        check("F", &[("C", 1), ("V", -1)]);
        check("Wb", &[("V", 1), ("s", 1)]);
        check("T", &[("Wb", 1), ("m", -2)]);
        check("H", &[("Wb", 1), ("A", -1)]);
        check("S", &[("ohm", -1)]);
        // hp = 550 ft lbf / s and lbf = lb g0 with g0 = 9.80665 m/s^2
        let hp = resolve_symbol("hp").unwrap_or_else(|| panic!("hp"));
        let ft_lbf = resolve_symbol("ft")
            .zip(resolve_symbol("lbf"))
            .map(|(a, b)| a.scale.rat * b.scale.rat);
        assert_eq!(
            Some(hp.scale.rat),
            ft_lbf.map(|r| r * BigRational::from_integer(BigInt::from(550)))
        );
        let lbf = resolve_symbol("lbf").unwrap_or_else(|| panic!("lbf"));
        assert_eq!(lbf.scale.rat, dec("0.45359237") * dec("9.80665"));
        // non-unit symbols are rejected
        assert!(resolve_symbol("parsec").is_none());
        assert!(resolve_symbol("xm").is_none());
    }

    #[test]
    fn unit_names_do_not_collide_with_operators() {
        let mut graph = crate::graph::Graph::new();
        crate::graph::Engine::install(&mut graph, &crate::rules::standard())
            .unwrap_or_else(|e| panic!("{e}"));
        let mut clashes: Vec<&str> = table()
            .keys()
            .copied()
            .filter(|name| graph.ops().lookup(name).is_some())
            .collect();
        clashes.extend(
            ALIASES
                .iter()
                .map(|(a, _)| *a)
                .filter(|name| graph.ops().lookup(name).is_some()),
        );
        clashes.sort_unstable();
        assert!(
            clashes.is_empty(),
            "unit names that are also operators: {clashes:?}"
        );
    }

    #[test]
    fn numeric_values_and_symbolic_mix() {
        assert_eq!(
            run("unify_expression(quantity(2.5, m) + quantity(50, cm))"),
            "quantity(3, m)"
        );
        assert_eq!(
            run("unify_expression(quantity(a, m) + quantity(b, m))"),
            "quantity(a + b, m)"
        );
        assert_eq!(
            run("unify_expression(quantity(a, km) + quantity(b, m))"),
            "quantity(a + 1/1000*b, km)"
        );
        // sums of products
        assert_eq!(
            run("unify_expression(quantity(1, m) * quantity(2, m) + quantity(3, m^2))"),
            "quantity(5, m^2)"
        );
        assert_eq!(
            run("unify_expression(quantity(1, kg) * quantity(10, m/s^2) * quantity(2, m))"),
            "quantity(20, kg*m^2/s^2)"
        );
    }
}
