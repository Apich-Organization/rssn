use num_bigint::BigInt;
use num_traits::{One, Zero};

use super::Rule;
use crate::symbolic::egraph::egraph::EGraph;
use crate::symbolic::egraph::enode::ENode;
use crate::symbolic::egraph::id::Id;

enum TrigOp {
    ZeroSin(Id),
    ZeroCos(Id),
    Tan(Id, Id),
    Sec(Id, Id),
    Csc(Id, Id),
    Cot(Id, Id),
    Pythagoras(Id),
    NegPythagoras(Id),
}

/// Trigonometric identities rule (Tier 1 & 2).
/// Includes: sin^2(x) + cos^2(x) = 1, tan(x) = sin(x)/cos(x), odd/even parity, etc.
#[derive(Debug)]
pub struct TrigIdentitiesRule;

impl Rule for TrigIdentitiesRule {
    fn name(&self) -> &str {
        "trig::identities"
    }

    fn tier(&self) -> u8 {
        1
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut ops: Vec<TrigOp> = Vec::new();

        let zero_id = egraph.add_node(ENode::BigInt(BigInt::zero()));
        let two_id = egraph.add_node(ENode::BigInt(BigInt::from(2)));
        let exp2 = egraph.find_immut(two_id);
        let zero_root = egraph.find_immut(zero_id);

        let is_sin_sq = |class_id: Id| -> Option<Id> {
            let canon = egraph.find_immut(class_id);
            if let Some(c) = egraph.get_class(canon) {
                for node in &c.nodes {
                    if let ENode::Power(base, exp) = node {
                        if egraph.find_immut(*exp) == exp2 {
                            let base_canon = egraph.find_immut(*base);
                            if let Some(bc) = egraph.get_class(base_canon) {
                                for bn in &bc.nodes {
                                    if let ENode::Sin(arg) = bn {
                                        return Some(egraph.find_immut(*arg));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            None
        };

        let is_cos_sq = |class_id: Id| -> Option<Id> {
            let canon = egraph.find_immut(class_id);
            if let Some(c) = egraph.get_class(canon) {
                for node in &c.nodes {
                    if let ENode::Power(base, exp) = node {
                        if egraph.find_immut(*exp) == exp2 {
                            let base_canon = egraph.find_immut(*base);
                            if let Some(bc) = egraph.get_class(base_canon) {
                                for bn in &bc.nodes {
                                    if let ENode::Cos(arg) = bn {
                                        return Some(egraph.find_immut(*arg));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            None
        };

        let is_neg_sin_sq = |class_id: Id| -> Option<Id> {
            let canon = egraph.find_immut(class_id);
            if let Some(c) = egraph.get_class(canon) {
                for node in &c.nodes {
                    if let ENode::Neg(inner) = node {
                        if let Some(arg) = is_sin_sq(*inner) {
                            return Some(arg);
                        }
                    }
                }
            }
            None
        };

        let is_neg_cos_sq = |class_id: Id| -> Option<Id> {
            let canon = egraph.find_immut(class_id);
            if let Some(c) = egraph.get_class(canon) {
                for node in &c.nodes {
                    if let ENode::Neg(inner) = node {
                        if let Some(arg) = is_cos_sq(*inner) {
                            return Some(arg);
                        }
                    }
                }
            }
            None
        };

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            for node in &class.nodes {
                match node {
                    // sin(0) = 0
                    ENode::Sin(arg) => {
                        let r = egraph.find_immut(*arg);
                        if r == zero_root || egraph.get_class(r).is_some_and(|c| c.is_zero()) {
                            ops.push(TrigOp::ZeroSin(class_id));
                        }
                    }
                    // cos(0) = 1
                    ENode::Cos(arg) => {
                        let r = egraph.find_immut(*arg);
                        if r == zero_root || egraph.get_class(r).is_some_and(|c| c.is_zero()) {
                            ops.push(TrigOp::ZeroCos(class_id));
                        }
                    }
                    // tan(x) <=> sin(x) / cos(x)
                    ENode::Tan(arg) => {
                        let r = egraph.find_immut(*arg);
                        if r == zero_root || egraph.get_class(r).is_some_and(|c| c.is_zero()) {
                            ops.push(TrigOp::ZeroSin(class_id));
                        } else {
                            ops.push(TrigOp::Tan(class_id, *arg));
                        }
                    }
                    // sec(x) <=> 1 / cos(x)
                    ENode::Sec(arg) => {
                        ops.push(TrigOp::Sec(class_id, *arg));
                    }
                    // csc(x) <=> 1 / sin(x)
                    ENode::Csc(arg) => {
                        ops.push(TrigOp::Csc(class_id, *arg));
                    }
                    // cot(x) <=> cos(x) / sin(x)
                    ENode::Cot(arg) => {
                        ops.push(TrigOp::Cot(class_id, *arg));
                    }

                    // Pythagorean: sin^2(x) + cos^2(x) = 1, -sin^2(x) - cos^2(x) = -1
                    ENode::Add(a, b) => {
                        let sin_a = is_sin_sq(*a);
                        let cos_b = is_cos_sq(*b);
                        if sin_a.is_some() && sin_a == cos_b {
                            ops.push(TrigOp::Pythagoras(class_id));
                        } else {
                            let cos_a = is_cos_sq(*a);
                            let sin_b = is_sin_sq(*b);
                            if cos_a.is_some() && cos_a == sin_b {
                                ops.push(TrigOp::Pythagoras(class_id));
                            } else {
                                let neg_sin_a = is_neg_sin_sq(*a);
                                let neg_cos_b = is_neg_cos_sq(*b);
                                if neg_sin_a.is_some() && neg_sin_a == neg_cos_b {
                                    ops.push(TrigOp::NegPythagoras(class_id));
                                } else {
                                    let neg_cos_a = is_neg_cos_sq(*a);
                                    let neg_sin_b = is_neg_sin_sq(*b);
                                    if neg_cos_a.is_some() && neg_cos_a == neg_sin_b {
                                        ops.push(TrigOp::NegPythagoras(class_id));
                                    }
                                }
                            }
                        }
                    }
                    ENode::Sub(a, b) => {
                        let neg_sin_a = is_neg_sin_sq(*a);
                        let cos_b = is_cos_sq(*b);
                        if neg_sin_a.is_some() && neg_sin_a == cos_b {
                            ops.push(TrigOp::NegPythagoras(class_id));
                        } else {
                            let neg_cos_a = is_neg_cos_sq(*a);
                            let sin_b = is_sin_sq(*b);
                            if neg_cos_a.is_some() && neg_cos_a == sin_b {
                                ops.push(TrigOp::NegPythagoras(class_id));
                            }
                        }
                    }
                    _ => {}
                }
            }
        }

        if ops.is_empty() {
            return 0;
        }

        let one_id = egraph.add_node(ENode::BigInt(BigInt::one()));
        let mut applied = 0;

        for op in ops {
            match op {
                TrigOp::ZeroSin(class_id) => {
                    if egraph.union(class_id, zero_id) {
                        applied += 1;
                    }
                }
                TrigOp::ZeroCos(class_id) => {
                    if egraph.union(class_id, one_id) {
                        applied += 1;
                    }
                }
                TrigOp::Tan(class_id, arg) => {
                    let sin_x = egraph.add_node(ENode::Sin(arg));
                    let cos_x = egraph.add_node(ENode::Cos(arg));
                    let div = egraph.add_node(ENode::Div(sin_x, cos_x));
                    if egraph.union(class_id, div) {
                        applied += 1;
                    }
                }
                TrigOp::Sec(class_id, arg) => {
                    let cos_x = egraph.add_node(ENode::Cos(arg));
                    let div = egraph.add_node(ENode::Div(one_id, cos_x));
                    if egraph.union(class_id, div) {
                        applied += 1;
                    }
                }
                TrigOp::Csc(class_id, arg) => {
                    let sin_x = egraph.add_node(ENode::Sin(arg));
                    let div = egraph.add_node(ENode::Div(one_id, sin_x));
                    if egraph.union(class_id, div) {
                        applied += 1;
                    }
                }
                TrigOp::Cot(class_id, arg) => {
                    let cos_x = egraph.add_node(ENode::Cos(arg));
                    let sin_x = egraph.add_node(ENode::Sin(arg));
                    let div = egraph.add_node(ENode::Div(cos_x, sin_x));
                    if egraph.union(class_id, div) {
                        applied += 1;
                    }
                }
                TrigOp::Pythagoras(class_id) => {
                    if egraph.union(class_id, one_id) {
                        applied += 1;
                    }
                }
                TrigOp::NegPythagoras(class_id) => {
                    let neg_one_id = egraph.add_node(ENode::BigInt(-BigInt::one()));
                    if egraph.union(class_id, neg_one_id) {
                        applied += 1;
                    }
                }
            }
        }

        applied
    }
}
