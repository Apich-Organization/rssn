use num_bigint::BigInt;
use num_traits::{One, Zero};
use ordered_float::OrderedFloat;

use super::Rule;
use crate::symbolic::egraph::egraph::EGraph;
use crate::symbolic::egraph::enode::ENode;
use crate::symbolic::egraph::id::Id;

/// Differentiation rules (Tier 0: Top priority de-cocooning).
/// Transforms `Derivative(f, x)` into its algebraic derivative components.
#[derive(Debug)]
pub struct DifferentiationRule;

impl Rule for DifferentiationRule {
    fn name(&self) -> &str {
        "calculus::differentiation"
    }

    fn tier(&self) -> u8 {
        0
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        // Step 1: Scan read-only and collect all (class_id, body_id, var) targets
        let mut targets: Vec<(Id, Id, String)> = Vec::new();

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            // If this class already contains a reduced (non-Derivative) node, it has already been de-cocooned
            if class.nodes.iter().any(|n| !matches!(n, ENode::Derivative(..))) {
                continue;
            }
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                if let ENode::Derivative(body, var) = node {
                    targets.push((class_id, *body, var.clone()));
                    break;
                }
            }
        }

        if targets.is_empty() {
            return 0;
        }

        let zero_id = egraph.add_node(ENode::BigInt(BigInt::zero()));
        let one_id = egraph.add_node(ENode::BigInt(BigInt::one()));
        let two_id = egraph.add_node(ENode::BigInt(BigInt::from(2)));

        let mut applied = 0;

        // Step 2: For each target, examine body nodes (by snapshotting) and generate derivation identities
        for (class_id, body_id, var) in targets {
            let canonical_body = egraph.find_immut(body_id);
            let (body_nodes, is_pure_constant) = match egraph.get_class(canonical_body) {
                Some(c) => (c.nodes.clone(), c.free_vars.is_empty()),
                None => continue,
            };
            let is_numeric = body_nodes.iter().any(|n| {
                matches!(n, ENode::Constant(_) | ENode::BigInt(_) | ENode::Rational(_) | ENode::Pi | ENode::E | ENode::Infinity | ENode::NegativeInfinity)
            });
            if is_numeric || is_pure_constant {
                if egraph.union(class_id, zero_id) {
                    applied += 1;
                }
                continue;
            }

            // Differentiate the best structural representation of the body
            let chosen_node = body_nodes.iter().cloned().find(|n| {
                matches!(n, ENode::Variable(v) if v == &var)
            }).or_else(|| body_nodes.iter().cloned().find(|n| {
                matches!(n, ENode::Add(..) | ENode::Sub(..) | ENode::Div(..))
            })).or_else(|| body_nodes.iter().cloned().find(|n| {
                matches!(n, ENode::Mul(..) | ENode::Power(..) | ENode::Sin(_) | ENode::Cos(_) | ENode::Tan(_) | ENode::Exp(_) | ENode::Log(_) | ENode::Neg(_))
            })).or_else(|| body_nodes.first().cloned());

            if let Some(body_node) = chosen_node {
                let derived_id = match body_node {
                    // Base case 2: Variable x -> 1; independent variables -> 0; unknown functions preserved
                    ENode::Variable(name) => {
                        if name == var {
                            one_id
                        } else if matches!(name.as_str(), "y" | "z" | "w" | "u" | "v" | "psi" | "phi" | "q" | "f" | "g" | "h")
                            || (var == "t" && matches!(name.as_str(), "x" | "y" | "z" | "theta")) {
                            continue;
                        } else {
                            zero_id
                        }
                    }

                    // Constants
                    ENode::Constant(_)
                    | ENode::BigInt(_)
                    | ENode::Rational(_)
                    | ENode::Pi
                    | ENode::E
                    | ENode::Infinity
                    | ENode::NegativeInfinity => zero_id,

                    // Linearity: d/dx(u + v) = d/dx(u) + d/dx(v)
                    ENode::Add(u, v) => {
                        let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                        let dv = egraph.add_node(ENode::Derivative(v, var.clone()));
                        egraph.add_node(ENode::Add(du, dv))
                    }

                    // Linearity: d/dx(u - v) = d/dx(u) - d/dx(v)
                    ENode::Sub(u, v) => {
                        let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                        let dv = egraph.add_node(ENode::Derivative(v, var.clone()));
                        egraph.add_node(ENode::Sub(du, dv))
                    }

                    // Linearity: d/dx(-u) = -(d/dx(u))
                    ENode::Neg(u) => {
                        let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                        egraph.add_node(ENode::Neg(du))
                    }

                    // Product Rule: d/dx(u * v) = u'*v + u*v'
                    ENode::Mul(u, v) => {
                        let u_has_var = egraph.get_class(u).is_some_and(|c| c.contains_var(&var));
                        let v_has_var = egraph.get_class(v).is_some_and(|c| c.contains_var(&var));

                        if !u_has_var && !v_has_var {
                            zero_id
                        } else if !u_has_var && v_has_var {
                            let dv = egraph.add_node(ENode::Derivative(v, var.clone()));
                            egraph.add_node(ENode::Mul(u, dv))
                        } else if u_has_var && !v_has_var {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            egraph.add_node(ENode::Mul(du, v))
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let dv = egraph.add_node(ENode::Derivative(v, var.clone()));
                            let t1 = egraph.add_node(ENode::Mul(du, v));
                            let t2 = egraph.add_node(ENode::Mul(u, dv));
                            egraph.add_node(ENode::Add(t1, t2))
                        }
                    }

                    // Quotient Rule: d/dx(u / v) = (u'*v - u*v') / v^2
                    ENode::Div(u, v) => {
                        let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                        let dv = egraph.add_node(ENode::Derivative(v, var.clone()));
                        let num1 = egraph.add_node(ENode::Mul(du, v));
                        let num2 = egraph.add_node(ENode::Mul(u, dv));
                        let num = egraph.add_node(ENode::Sub(num1, num2));
                        let den = egraph.add_node(ENode::Power(v, two_id));
                        egraph.add_node(ENode::Div(num, den))
                    }

                    // Power Rule: d/dx(u^n)
                    ENode::Power(base, exp) => {
                        let exp_has_var = egraph.get_class(exp).is_some_and(|c| c.contains_var(&var));
                        let base_has_var = egraph.get_class(base).is_some_and(|c| c.contains_var(&var));

                        let is_exp_zero = egraph.get_class(exp).is_some_and(|c| c.is_zero());
                        let is_exp_one = egraph.get_class(exp).is_some_and(|c| c.is_one());

                        if is_exp_zero || (!exp_has_var && !base_has_var) {
                            zero_id
                        } else if is_exp_one {
                            egraph.add_node(ENode::Derivative(base, var.clone()))
                        } else if !exp_has_var && base_has_var {
                            let dbase = egraph.add_node(ENode::Derivative(base, var.clone()));
                            let n_minus_1 = egraph.add_node(ENode::Sub(exp, one_id));
                            let new_pow = egraph.add_node(ENode::Power(base, n_minus_1));
                            let prod1 = egraph.add_node(ENode::Mul(exp, new_pow));
                            egraph.add_node(ENode::Mul(prod1, dbase))
                        } else if exp_has_var && !base_has_var {
                            let dexp = egraph.add_node(ENode::Derivative(exp, var.clone()));
                            let ln_base = egraph.add_node(ENode::Log(base));
                            let pow = egraph.add_node(ENode::Power(base, exp));
                            let t1 = egraph.add_node(ENode::Mul(pow, ln_base));
                            egraph.add_node(ENode::Mul(t1, dexp))
                        } else {
                            continue;
                        }
                    }

                    // Chain Rule: d/dx(sin(u)) = cos(u) * u'
                    ENode::Sin(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let cos_u = egraph.add_node(ENode::Cos(u));
                            egraph.add_node(ENode::Mul(cos_u, du))
                        }
                    }

                    // Chain Rule: d/dx(cos(u)) = -sin(u) * u'
                    ENode::Cos(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let sin_u = egraph.add_node(ENode::Sin(u));
                            let neg_sin = egraph.add_node(ENode::Neg(sin_u));
                            egraph.add_node(ENode::Mul(neg_sin, du))
                        }
                    }

                    // Chain Rule: d/dx(exp(u)) = exp(u) * u'
                    ENode::Exp(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let exp_u = egraph.add_node(ENode::Exp(u));
                            egraph.add_node(ENode::Mul(exp_u, du))
                        }
                    }

                    // Chain Rule: d/dx(ln(u)) = u' / u
                    ENode::Log(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            egraph.add_node(ENode::Div(du, u))
                        }
                    }

                    // Chain Rule: d/dx(tan(u)) = sec(u)^2 * u'
                    ENode::Tan(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let sec_u = egraph.add_node(ENode::Sec(u));
                            let sec2 = egraph.add_node(ENode::Power(sec_u, two_id));
                            egraph.add_node(ENode::Mul(sec2, du))
                        }
                    }

                    // Chain Rule: d/dx(sqrt(u)) = u' / (2 * sqrt(u))
                    ENode::Sqrt(u) => {
                        if !egraph.get_class(u).is_some_and(|c| c.contains_var(&var)) {
                            zero_id
                        } else {
                            let du = egraph.add_node(ENode::Derivative(u, var.clone()));
                            let sqrt_u = egraph.add_node(ENode::Sqrt(u));
                            let den = egraph.add_node(ENode::Mul(two_id, sqrt_u));
                            egraph.add_node(ENode::Div(du, den))
                        }
                    }

                    // Nested derivative: d/dx(d/dy(body))
                    // If same variable: d/dx(d/dx(body)) = d²/dx²(body)
                    // If different variable: d/dx(d/dy(body)) = Derivative(Derivative(body, y), x)
                    ENode::Derivative(inner_body, inner_var) => {
                        if inner_var == var {
                            // Same variable: produce second derivative
                            let n2 = egraph.add_node(ENode::Constant(OrderedFloat(2.0)));
                            egraph.add_node(ENode::DerivativeN(inner_body, var.clone(), n2))
                        } else {
                            // Mixed partial: wrap as nested derivative
                            let inner_deriv = egraph.add_node(ENode::Derivative(inner_body, inner_var));
                            egraph.add_node(ENode::Derivative(inner_deriv, var.clone()))
                        }
                    }

                    // Nested n-th derivative: d/dx(d^n/dx^n(body)) = d^(n+1)/dx^(n+1)(body)
                    ENode::DerivativeN(inner_body, inner_var, n_id) if inner_var == var => {
                        let n_plus_1 = egraph.add_node(ENode::Add(n_id, one_id));
                        egraph.add_node(ENode::DerivativeN(inner_body, var.clone(), n_plus_1))
                    }

                    _ => continue,
                };

                if egraph.union(class_id, derived_id) {
                    applied += 1;
                }
            }
        }

        applied
    }
}
