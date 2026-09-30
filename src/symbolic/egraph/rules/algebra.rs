use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero, ToPrimitive};
use ordered_float::OrderedFloat;

use super::Rule;
use crate::symbolic::egraph::eclass::EClass;
use crate::symbolic::egraph::egraph::EGraph;
use crate::symbolic::egraph::enode::ENode;
use crate::symbolic::egraph::id::Id;

#[derive(Clone, Debug, PartialEq)]
enum NumVal {
    Int(BigInt),
    Rat(BigRational),
    Float(f64),
}

impl NumVal {
    fn to_enode(&self) -> ENode {
        match self {
            NumVal::Int(i) => ENode::BigInt(i.clone()),
            NumVal::Rat(r) => {
                if r.is_integer() {
                    ENode::BigInt(r.to_integer())
                } else {
                    ENode::Rational(r.clone())
                }
            }
            NumVal::Float(f) => ENode::Constant(OrderedFloat(*f)),
        }
    }

    fn add(&self, other: &Self) -> Self {
        match (self, other) {
            (NumVal::Int(a), NumVal::Int(b)) => NumVal::Int(a + b),
            (NumVal::Int(a), NumVal::Rat(b)) => NumVal::Rat(BigRational::from(a.clone()) + b),
            (NumVal::Rat(a), NumVal::Int(b)) => NumVal::Rat(a + BigRational::from(b.clone())),
            (NumVal::Rat(a), NumVal::Rat(b)) => NumVal::Rat(a + b),
            (NumVal::Float(a), NumVal::Float(b)) => NumVal::Float(a + b),
            (NumVal::Float(a), NumVal::Int(b)) => NumVal::Float(a + b.to_f64().unwrap_or(0.0)),
            (NumVal::Int(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) + b),
            (NumVal::Float(a), NumVal::Rat(b)) => NumVal::Float(a + b.to_f64().unwrap_or(0.0)),
            (NumVal::Rat(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) + b),
        }
    }

    fn sub(&self, other: &Self) -> Self {
        match (self, other) {
            (NumVal::Int(a), NumVal::Int(b)) => NumVal::Int(a - b),
            (NumVal::Int(a), NumVal::Rat(b)) => NumVal::Rat(BigRational::from(a.clone()) - b),
            (NumVal::Rat(a), NumVal::Int(b)) => NumVal::Rat(a - BigRational::from(b.clone())),
            (NumVal::Rat(a), NumVal::Rat(b)) => NumVal::Rat(a - b),
            (NumVal::Float(a), NumVal::Float(b)) => NumVal::Float(a - b),
            (NumVal::Float(a), NumVal::Int(b)) => NumVal::Float(a - b.to_f64().unwrap_or(0.0)),
            (NumVal::Int(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) - b),
            (NumVal::Float(a), NumVal::Rat(b)) => NumVal::Float(a - b.to_f64().unwrap_or(0.0)),
            (NumVal::Rat(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) - b),
        }
    }

    fn mul(&self, other: &Self) -> Self {
        match (self, other) {
            (NumVal::Int(a), NumVal::Int(b)) => NumVal::Int(a * b),
            (NumVal::Int(a), NumVal::Rat(b)) => NumVal::Rat(BigRational::from(a.clone()) * b),
            (NumVal::Rat(a), NumVal::Int(b)) => NumVal::Rat(a * BigRational::from(b.clone())),
            (NumVal::Rat(a), NumVal::Rat(b)) => NumVal::Rat(a * b),
            (NumVal::Float(a), NumVal::Float(b)) => NumVal::Float(a * b),
            (NumVal::Float(a), NumVal::Int(b)) => NumVal::Float(a * b.to_f64().unwrap_or(0.0)),
            (NumVal::Int(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) * b),
            (NumVal::Float(a), NumVal::Rat(b)) => NumVal::Float(a * b.to_f64().unwrap_or(0.0)),
            (NumVal::Rat(a), NumVal::Float(b)) => NumVal::Float(a.to_f64().unwrap_or(0.0) * b),
        }
    }

    fn div(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (_, NumVal::Int(b)) if b.is_zero() => None,
            (_, NumVal::Rat(b)) if b.is_zero() => None,
            (_, NumVal::Float(b)) if *b == 0.0 => None,
            (NumVal::Int(a), NumVal::Int(b)) => {
                if a % b == BigInt::zero() {
                    Some(NumVal::Int(a / b))
                } else {
                    Some(NumVal::Rat(BigRational::new(a.clone(), b.clone())))
                }
            }
            (NumVal::Int(a), NumVal::Rat(b)) => Some(NumVal::Rat(BigRational::from(a.clone()) / b)),
            (NumVal::Rat(a), NumVal::Int(b)) => Some(NumVal::Rat(a / BigRational::from(b.clone()))),
            (NumVal::Rat(a), NumVal::Rat(b)) => Some(NumVal::Rat(a / b)),
            (NumVal::Float(a), NumVal::Float(b)) => Some(NumVal::Float(a / b)),
            (NumVal::Float(a), NumVal::Int(b)) => Some(NumVal::Float(a / b.to_f64().unwrap_or(1.0))),
            (NumVal::Int(a), NumVal::Float(b)) => Some(NumVal::Float(a.to_f64().unwrap_or(0.0) / b)),
            (NumVal::Float(a), NumVal::Rat(b)) => Some(NumVal::Float(a / b.to_f64().unwrap_or(1.0))),
            (NumVal::Rat(a), NumVal::Float(b)) => Some(NumVal::Float(a.to_f64().unwrap_or(0.0) / b)),
        }
    }

    fn pow(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (NumVal::Int(a), NumVal::Int(b)) => {
                if let Some(exp_u32) = b.to_u32() {
                    Some(NumVal::Int(a.pow(exp_u32)))
                } else if let Some(exp_i32) = b.to_i32() {
                    if exp_i32 < 0 && !a.is_zero() {
                        let inv = BigRational::new(BigInt::one(), a.pow((-exp_i32) as u32));
                        Some(NumVal::Rat(inv))
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            (NumVal::Float(a), NumVal::Float(b)) => Some(NumVal::Float(a.powf(*b))),
            (NumVal::Float(a), NumVal::Int(b)) => Some(NumVal::Float(a.powf(b.to_f64().unwrap_or(0.0)))),
            (NumVal::Int(a), NumVal::Float(b)) => Some(NumVal::Float(a.to_f64().unwrap_or(0.0).powf(*b))),
            (NumVal::Rat(a), NumVal::Int(b)) => {
                if let Some(exp_i32) = b.to_i32() {
                    if exp_i32 >= 0 {
                        Some(NumVal::Rat(a.pow(exp_i32)))
                    } else {
                        Some(NumVal::Rat(a.recip().pow(-exp_i32)))
                    }
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    fn neg(&self) -> Self {
        match self {
            NumVal::Int(a) => NumVal::Int(-a.clone()),
            NumVal::Rat(a) => NumVal::Rat(-a.clone()),
            NumVal::Float(a) => NumVal::Float(-a),
        }
    }

    fn sqrt(&self) -> Option<Self> {
        match self {
            NumVal::Int(a) => {
                if a >= &BigInt::zero() {
                    let s = a.sqrt();
                    if &s * &s == *a {
                        Some(NumVal::Int(s))
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            NumVal::Rat(a) => {
                if a >= &BigRational::zero() {
                    let num = a.numer();
                    let den = a.denom();
                    let s_num = num.sqrt();
                    let s_den = den.sqrt();
                    if &s_num * &s_num == *num && &s_den * &s_den == *den {
                        Some(NumVal::Rat(BigRational::new(s_num, s_den)))
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            NumVal::Float(a) => {
                if *a >= 0.0 {
                    let s = a.sqrt();
                    let r = s.round();
                    if (r - s).abs() < 1e-12 && ((r * r) - *a).abs() < 1e-12 {
                        Some(NumVal::Int(BigInt::from(r as i64)))
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
        }
    }

    fn abs(&self) -> Self {
        match self {
            NumVal::Int(a) => {
                if a < &BigInt::zero() {
                    NumVal::Int(-a.clone())
                } else {
                    NumVal::Int(a.clone())
                }
            }
            NumVal::Rat(a) => {
                if a < &BigRational::zero() {
                    NumVal::Rat(-a.clone())
                } else {
                    NumVal::Rat(a.clone())
                }
            }
            NumVal::Float(a) => NumVal::Float(a.abs()),
        }
    }
}

fn get_num(class: &EClass) -> Option<NumVal> {
    for node in &class.nodes {
        match node {
            ENode::BigInt(b) => return Some(NumVal::Int(b.clone())),
            ENode::Rational(r) => return Some(NumVal::Rat(r.clone())),
            ENode::Constant(c) => return Some(NumVal::Float(c.0)),
            _ => {}
        }
    }
    None
}

/// Constant folding rule (Tier 1: strictly reduces complexity).
#[derive(Debug)]
pub struct ConstantFoldingRule;

impl Rule for ConstantFoldingRule {
    fn name(&self) -> &str {
        "algebra::constant_folding"
    }

    fn tier(&self) -> u8 {
        1
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut rewrites = Vec::new();

        let zero_id = egraph.add_node(ENode::BigInt(BigInt::zero()));
        let one_id = egraph.add_node(ENode::BigInt(BigInt::one()));

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            // Canonicalize 0 and 1 representations
            if class.is_zero() && class_id != zero_id {
                rewrites.push((class_id, ENode::BigInt(BigInt::zero())));
            }
            if class.is_one() && class_id != one_id {
                rewrites.push((class_id, ENode::BigInt(BigInt::one())));
            }

            for node in &class.nodes {
                match node {
                    ENode::Add(a, b) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        let nb = egraph.get_class(*b).and_then(get_num);
                        if let (Some(v1), Some(v2)) = (na, nb) {
                            rewrites.push((class_id, v1.add(&v2).to_enode()));
                        }
                    }
                    ENode::Sub(a, b) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        let nb = egraph.get_class(*b).and_then(get_num);
                        if let (Some(v1), Some(v2)) = (na, nb) {
                            rewrites.push((class_id, v1.sub(&v2).to_enode()));
                        }
                    }
                    ENode::Mul(a, b) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        let nb = egraph.get_class(*b).and_then(get_num);
                        if let (Some(v1), Some(v2)) = (na, nb) {
                            rewrites.push((class_id, v1.mul(&v2).to_enode()));
                        }
                    }
                    ENode::Div(a, b) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        let nb = egraph.get_class(*b).and_then(get_num);
                        if let (Some(v1), Some(v2)) = (na, nb) {
                            if let Some(res) = v1.div(&v2) {
                                rewrites.push((class_id, res.to_enode()));
                            }
                        }
                    }
                    ENode::Power(a, b) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        let nb = egraph.get_class(*b).and_then(get_num);
                        if let (Some(v1), Some(v2)) = (na, nb) {
                            if let Some(res) = v1.pow(&v2) {
                                rewrites.push((class_id, res.to_enode()));
                            }
                        }
                    }
                    ENode::Neg(a) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        if let Some(v) = na {
                            rewrites.push((class_id, v.neg().to_enode()));
                        }
                    }
                    ENode::Sqrt(a) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        if let Some(v) = na {
                            if let Some(res) = v.sqrt() {
                                rewrites.push((class_id, res.to_enode()));
                            }
                        }
                    }
                    ENode::Abs(a) => {
                        let na = egraph.get_class(*a).and_then(get_num);
                        if let Some(v) = na {
                            rewrites.push((class_id, v.abs().to_enode()));
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut applied = 0;
        for (class_id, new_node) in rewrites {
            let new_id = egraph.add_node(new_node);
            if egraph.union(class_id, new_id) {
                applied += 1;
            }
        }
        applied
    }
}

/// Identity reduction rule (Tier 1: strictly reduces complexity).
/// Covers: x+0->x, x*1->x, x*0->0, x^1->x, x^0->1, x-x->0, x/x->1, etc.
#[derive(Debug)]
pub struct IdentityReductionRule;

impl Rule for IdentityReductionRule {
    fn name(&self) -> &str {
        "algebra::identity_reduction"
    }

    fn tier(&self) -> u8 {
        1
    }

    fn apply(&self, egraph: &mut EGraph) -> usize {
        let mut unions = Vec::new();
        let mut new_nodes = Vec::new();
        let mut neg_muls = Vec::new();
        let mut assoc_muls = Vec::new();
        let mut commutative_muls: Vec<(Id, Id, Id)> = Vec::new(); // (class_id, a, b) → Mul(b, a)
        // (class_id, inner_enode, outer_id) → Mul(outer_id, add_node(inner_enode))
        let mut deferred_nested_muls: Vec<(Id, ENode, Id)> = Vec::new();
        // (class_id, base_id, exp_enode) → Power(base_id, add_node(exp_enode))
        let mut deferred_powers: Vec<(Id, Id, ENode)> = Vec::new();

        let zero_id = egraph.add_node(ENode::BigInt(BigInt::zero()));
        let one_id = egraph.add_node(ENode::BigInt(BigInt::one()));
        let two_id = egraph.add_node(ENode::BigInt(BigInt::from(2)));
        let neg_one_id = egraph.add_node(ENode::BigInt(BigInt::from(-1)));
        let r_zero = egraph.find_immut(zero_id);
        let r_one = egraph.find_immut(one_id);
        let r_neg_one = egraph.find_immut(neg_one_id);

        for class_idx in 0..egraph.classes.len() {
            let class = match &egraph.classes[class_idx] {
                Some(c) => c,
                None => continue,
            };
            let class_id = Id::from_usize(class_idx);

            // Canonicalize 0, 1, -1
            if class.is_zero() && class_id != zero_id {
                unions.push((class_id, zero_id));
            }
            if class.is_one() && class_id != one_id {
                unions.push((class_id, one_id));
            }
            if class.is_neg_one() && class_id != neg_one_id {
                unions.push((class_id, neg_one_id));
            }

            for node in &class.nodes {
                match node {
                    ENode::Add(a, b) => {
                        let ra = egraph.find_immut(*a);
                        let rb = egraph.find_immut(*b);
                        let is_a_zero = ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero());
                        let is_b_zero = rb == r_zero || egraph.get_class(rb).is_some_and(|c| c.is_zero());

                        if is_a_zero {
                            unions.push((class_id, rb));
                        } else if is_b_zero {
                            unions.push((class_id, ra));
                        } else if ra == rb {
                            new_nodes.push((class_id, ENode::Mul(two_id, *a)));
                        } else {
                            let a_is_neg_b = egraph.get_class(ra).is_some_and(|c| {
                                c.nodes.iter().any(|n| match n {
                                    ENode::Neg(inner) => egraph.find_immut(*inner) == rb,
                                    _ => false,
                                })
                            });
                            let b_is_neg_a = egraph.get_class(rb).is_some_and(|c| {
                                c.nodes.iter().any(|n| match n {
                                    ENode::Neg(inner) => egraph.find_immut(*inner) == ra,
                                    _ => false,
                                })
                            });
                            if a_is_neg_b || b_is_neg_a {
                                unions.push((class_id, r_zero));
                            }

                            // (u - b) + b => u
                            if let Some(ca) = egraph.get_class(ra) {
                                for n in &ca.nodes {
                                    if let ENode::Sub(u, v) = n {
                                        if egraph.find_immut(*v) == rb {
                                            unions.push((class_id, *u));
                                        }
                                    }
                                }
                            }
                            // a + (u - a) => u
                            if let Some(cb) = egraph.get_class(rb) {
                                for n in &cb.nodes {
                                    if let ENode::Sub(u, v) = n {
                                        if egraph.find_immut(*v) == ra {
                                            unions.push((class_id, *u));
                                        }
                                    }
                                }
                            }
                            // a + (-b) => a - b
                            // (-a) + b => b - a
                            if let Some(ca) = egraph.get_class(ra) {
                                for na in &ca.nodes {
                                    if let ENode::Neg(inner) = na {
                                        new_nodes.push((class_id, ENode::Sub(*b, *inner)));
                                    }
                                }
                            }
                            if let Some(cb) = egraph.get_class(rb) {
                                for nb in &cb.nodes {
                                    if let ENode::Neg(inner) = nb {
                                        new_nodes.push((class_id, ENode::Sub(*a, *inner)));
                                    }
                                }
                            }
                        }
                    }
                    ENode::Sub(a, b) => {
                        let ra = egraph.find_immut(*a);
                        let rb = egraph.find_immut(*b);
                        let is_a_zero = ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero());
                        let is_b_zero = rb == r_zero || egraph.get_class(rb).is_some_and(|c| c.is_zero());

                        if ra == rb {
                            unions.push((class_id, r_zero));
                        } else if is_b_zero {
                            unions.push((class_id, ra));
                        } else if is_a_zero {
                            new_nodes.push((class_id, ENode::Neg(*b)));
                        } else {
                            // a - (-b) => a + b
                            if let Some(cb) = egraph.get_class(rb) {
                                for nb in &cb.nodes {
                                    if let ENode::Neg(inner) = nb {
                                        new_nodes.push((class_id, ENode::Add(*a, *inner)));
                                        // (u - v) - (-w) where v == w => u
                                        if let Some(ca) = egraph.get_class(ra) {
                                            for na in &ca.nodes {
                                                if let ENode::Sub(u, v) = na {
                                                    if egraph.find_immut(*v) == egraph.find_immut(*inner) {
                                                        unions.push((class_id, *u));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            // (u - v) - b: if u == b => -v
                            if let Some(ca) = egraph.get_class(ra) {
                                for na in &ca.nodes {
                                    if let ENode::Sub(u, v) = na {
                                        if egraph.find_immut(*u) == rb {
                                            new_nodes.push((class_id, ENode::Neg(*v)));
                                        }
                                    }
                                }
                            }
                            // a - (x - y): if a == x => y
                            if let Some(cb) = egraph.get_class(rb) {
                                for nb in &cb.nodes {
                                    if let ENode::Sub(x, y) = nb {
                                        if ra == egraph.find_immut(*x) {
                                            unions.push((class_id, *y));
                                        }
                                    }
                                }
                            }

                            // (u + v) - b
                            if let Some(ca) = egraph.get_class(ra) {
                                for na in &ca.nodes {
                                    if let ENode::Add(u, v) = na {
                                        let ru = egraph.find_immut(*u);
                                        let rv = egraph.find_immut(*v);
                                        if ru == rb {
                                            unions.push((class_id, *v));
                                        } else if rv == rb {
                                            unions.push((class_id, *u));
                                        }
                                    }
                                }
                            }
                            // Term cancellation: (u + v) - (x + y)
                            if let (Some(ca), Some(cb)) = (egraph.get_class(ra), egraph.get_class(rb)) {
                                for na in &ca.nodes {
                                    if let ENode::Add(u, v) = na {
                                        let ru = egraph.find_immut(*u);
                                        let rv = egraph.find_immut(*v);
                                        for nb in &cb.nodes {
                                            if let ENode::Add(x, y) = nb {
                                                let rx = egraph.find_immut(*x);
                                                let ry = egraph.find_immut(*y);
                                                if ru == rx {
                                                    new_nodes.push((class_id, ENode::Sub(*v, *y)));
                                                } else if ru == ry {
                                                    new_nodes.push((class_id, ENode::Sub(*v, *x)));
                                                } else if rv == rx {
                                                    new_nodes.push((class_id, ENode::Sub(*u, *y)));
                                                } else if rv == ry {
                                                    new_nodes.push((class_id, ENode::Sub(*u, *x)));
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Distributive factoring: Mul(a, x) - Mul(b, x) → Mul(Sub(a, b), x)
                            // Also: Mul(x, a) - Mul(x, b) → Mul(x, Sub(a, b))
                            if let (Some(ca), Some(cb)) = (egraph.get_class(ra), egraph.get_class(rb)) {
                                for na in &ca.nodes {
                                    if let ENode::Mul(a1, a2) = na {
                                        let ra1 = egraph.find_immut(*a1);
                                        let ra2 = egraph.find_immut(*a2);
                                        for nb in &cb.nodes {
                                            if let ENode::Mul(b1, b2) = nb {
                                                let rb1 = egraph.find_immut(*b1);
                                                let rb2 = egraph.find_immut(*b2);
                                                // Mul(a1, a2) - Mul(b1, b2) where a2==b2 → Mul(a1-b1, a2)
                                                if ra2 == rb2 {
                                                    deferred_nested_muls.push((class_id, ENode::Sub(*a1, *b1), *a2));
                                                }
                                                // Mul(a1, a2) - Mul(b1, b2) where a1==b1 → Mul(a1, a2-b2)
                                                if ra1 == rb1 {
                                                    deferred_nested_muls.push((class_id, ENode::Sub(*a2, *b2), *a1));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ENode::Mul(a, b) => {
                        let ra = egraph.find_immut(*a);
                        let rb = egraph.find_immut(*b);
                        let is_a_zero = ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero());
                        let is_b_zero = rb == r_zero || egraph.get_class(rb).is_some_and(|c| c.is_zero());
                        let is_a_one = ra == r_one || egraph.get_class(ra).is_some_and(|c| c.is_one());
                        let is_b_one = rb == r_one || egraph.get_class(rb).is_some_and(|c| c.is_one());
                        let is_a_neg_one = ra == r_neg_one || egraph.get_class(ra).is_some_and(|c| c.is_neg_one());
                        let is_b_neg_one = rb == r_neg_one || egraph.get_class(rb).is_some_and(|c| c.is_neg_one());

                        if is_a_zero || is_b_zero {
                            unions.push((class_id, r_zero));
                        } else if is_a_one {
                            unions.push((class_id, rb));
                        } else if is_b_one {
                            unions.push((class_id, ra));
                        } else if is_a_neg_one {
                            new_nodes.push((class_id, ENode::Neg(*b)));
                        } else if is_b_neg_one {
                            new_nodes.push((class_id, ENode::Neg(*a)));
                        }

                        // Commutativity: a * b ≡ b * a (canonical orientation: numbers to left, else ra > rb)
                        let a_is_num = egraph.get_class(ra).and_then(get_num).is_some();
                        let b_is_num = egraph.get_class(rb).and_then(get_num).is_some();
                        let should_commute = if !a_is_num && b_is_num {
                            true
                        } else if a_is_num && !b_is_num {
                            false
                        } else {
                            ra > rb
                        };

                        if should_commute && !is_a_zero && !is_b_zero && !is_a_one && !is_b_one {
                            let already_has = egraph.get_class(class_id).is_some_and(|c| {
                                c.nodes.iter().any(|n| match n {
                                    ENode::Mul(x, y) => egraph.find_immut(*x) == rb && egraph.find_immut(*y) == ra,
                                    _ => false,
                                })
                            });
                            if !already_has {
                                commutative_muls.push((class_id, rb, ra));
                            }
                        }

                        // a * a => a^2
                        if ra == rb && !is_a_zero && !is_a_one && !is_a_neg_one {
                            new_nodes.push((class_id, ENode::Power(*a, two_id)));
                        }

                        // Power combinations: a * (a^exp) => a^(exp+1)
                        if let Some(cb) = egraph.get_class(rb) {
                            for nb in &cb.nodes {
                                if let ENode::Power(base, exp) = nb {
                                    if ra == egraph.find_immut(*base) {
                                        deferred_powers.push((class_id, *base, ENode::Add(*exp, one_id)));
                                    }
                                }
                            }
                        }
                        if let Some(ca) = egraph.get_class(ra) {
                            for na in &ca.nodes {
                                if let ENode::Power(base, exp) = na {
                                    if rb == egraph.find_immut(*base) {
                                        deferred_powers.push((class_id, *base, ENode::Add(*exp, one_id)));
                                    }
                                }
                            }
                        }

                        // a * (u / a) => u, (w * a) * (u / a) => w * u, a * a^(-1) => 1
                        if let Some(cb) = egraph.get_class(rb) {
                            for nb in &cb.nodes {
                                if let ENode::Div(u, v) = nb {
                                    let rv = egraph.find_immut(*v);
                                    if ra == rv {
                                        unions.push((class_id, *u));
                                    } else if let Some(ca) = egraph.get_class(ra) {
                                        for na in &ca.nodes {
                                            if let ENode::Mul(w, z) = na {
                                                if egraph.find_immut(*z) == rv {
                                                    new_nodes.push((class_id, ENode::Mul(*w, *u)));
                                                } else if egraph.find_immut(*w) == rv {
                                                    new_nodes.push((class_id, ENode::Mul(*z, *u)));
                                                }
                                            }
                                        }
                                    }
                                } else if let ENode::Power(base, exp) = nb {
                                    if egraph.find_immut(*exp) == r_neg_one && ra == egraph.find_immut(*base) {
                                        unions.push((class_id, r_one));
                                    }
                                }
                            }
                        }
                        if let Some(ca) = egraph.get_class(ra) {
                            for na in &ca.nodes {
                                if let ENode::Div(u, v) = na {
                                    let rv = egraph.find_immut(*v);
                                    if rb == rv {
                                        unions.push((class_id, *u));
                                    } else if let Some(cb) = egraph.get_class(rb) {
                                        for nb in &cb.nodes {
                                            if let ENode::Mul(w, z) = nb {
                                                if egraph.find_immut(*z) == rv {
                                                    new_nodes.push((class_id, ENode::Mul(*w, *u)));
                                                } else if egraph.find_immut(*w) == rv {
                                                    new_nodes.push((class_id, ENode::Mul(*z, *u)));
                                                }
                                            }
                                        }
                                    }
                                } else if let ENode::Power(base, exp) = na {
                                    if egraph.find_immut(*exp) == r_neg_one && rb == egraph.find_immut(*base) {
                                        unions.push((class_id, r_one));
                                    }
                                }
                            }
                        }

                        // Numeric associativity: (c1 * x) * c2 => (c1 * c2) * x
                        if !is_a_zero && !is_b_zero && !is_a_one && !is_b_one {
                            if let Some(nb_val) = egraph.get_class(rb).and_then(get_num) {
                                if let Some(ca) = egraph.get_class(ra) {
                                    for na in &ca.nodes {
                                        if let ENode::Mul(u, v) = na {
                                            if let Some(nu_val) = egraph.get_class(*u).and_then(get_num) {
                                                let n_prod = nb_val.mul(&nu_val).to_enode();
                                                assoc_muls.push((class_id, n_prod, *v));
                                            } else if let Some(nv_val) = egraph.get_class(*v).and_then(get_num) {
                                                let n_prod = nb_val.mul(&nv_val).to_enode();
                                                assoc_muls.push((class_id, n_prod, *u));
                                            }
                                        }
                                    }
                                }
                            }
                            if let Some(na_val) = egraph.get_class(ra).and_then(get_num) {
                                if let Some(cb) = egraph.get_class(rb) {
                                    for nb in &cb.nodes {
                                        if let ENode::Mul(u, v) = nb {
                                            if let Some(nu_val) = egraph.get_class(*u).and_then(get_num) {
                                                let n_prod = na_val.mul(&nu_val).to_enode();
                                                assoc_muls.push((class_id, n_prod, *v));
                                            } else if let Some(nv_val) = egraph.get_class(*v).and_then(get_num) {
                                                let n_prod = na_val.mul(&nv_val).to_enode();
                                                assoc_muls.push((class_id, n_prod, *u));
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if let Some(ca) = egraph.get_class(ra) {
                            for na in &ca.nodes {
                                if let ENode::Neg(inner) = na {
                                    neg_muls.push((class_id, *inner, *b));
                                }
                            }
                        }
                        if let Some(cb) = egraph.get_class(rb) {
                            for nb in &cb.nodes {
                                if let ENode::Neg(inner) = nb {
                                    neg_muls.push((class_id, *a, *inner));
                                }
                            }
                        }
                    }
                    ENode::Div(a, b) => {
                        let ra = egraph.find_immut(*a);
                        let rb = egraph.find_immut(*b);
                        let is_a_zero = ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero());
                        let is_b_zero = rb == r_zero || egraph.get_class(rb).is_some_and(|c| c.is_zero());
                        let is_b_one = rb == r_one || egraph.get_class(rb).is_some_and(|c| c.is_one());

                        if ra == rb && !is_a_zero {
                            unions.push((class_id, r_one));
                        } else if is_b_one {
                            unions.push((class_id, ra));
                        } else if is_a_zero && !is_b_zero {
                            unions.push((class_id, r_zero));
                        } else if !is_b_zero {
                            // (u * b) / b => u or (b * v) / b => v
                            if let Some(ca) = egraph.get_class(ra) {
                                for na in &ca.nodes {
                                    if let ENode::Mul(u, v) = na {
                                        if egraph.find_immut(*u) == rb {
                                            unions.push((class_id, *v));
                                        } else if egraph.find_immut(*v) == rb {
                                            unions.push((class_id, *u));
                                        }

                                        let ru = egraph.find_immut(*u);
                                        let rv = egraph.find_immut(*v);
                                        if let Some(cu) = egraph.get_class(ru) {
                                            for nu in &cu.nodes {
                                                if let ENode::Mul(x, y) = nu {
                                                    if egraph.find_immut(*x) == rb {
                                                        new_nodes.push((class_id, ENode::Mul(*y, *v)));
                                                    } else if egraph.find_immut(*y) == rb {
                                                        new_nodes.push((class_id, ENode::Mul(*x, *v)));
                                                    }
                                                }
                                            }
                                        }
                                        if let Some(cv) = egraph.get_class(rv) {
                                            for nv in &cv.nodes {
                                                if let ENode::Mul(x, y) = nv {
                                                    if egraph.find_immut(*x) == rb {
                                                        new_nodes.push((class_id, ENode::Mul(*u, *y)));
                                                    } else if egraph.find_immut(*y) == rb {
                                                        new_nodes.push((class_id, ENode::Mul(*u, *x)));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    ENode::Power(a, b) => {
                        let ra = egraph.find_immut(*a);
                        let rb = egraph.find_immut(*b);
                        let is_a_one = ra == r_one || egraph.get_class(ra).is_some_and(|c| c.is_one());
                        let is_b_zero = rb == r_zero || egraph.get_class(rb).is_some_and(|c| c.is_zero());
                        let is_b_one = rb == r_one || egraph.get_class(rb).is_some_and(|c| c.is_one());

                        if is_b_zero {
                            unions.push((class_id, r_one));
                        } else if is_b_one {
                            unions.push((class_id, ra));
                        } else if is_a_one {
                            unions.push((class_id, r_one));
                        }

                        // i^2 -> -1
                        let is_b_two = egraph.get_class(rb).is_some_and(|c| {
                            c.nodes.iter().any(|n| match n {
                                ENode::BigInt(i) => *i == BigInt::from(2),
                                ENode::Constant(val) => val.0 == 2.0,
                                _ => false,
                            })
                        });
                        let is_a_i = egraph.get_class(ra).is_some_and(|c| {
                            c.nodes.iter().any(|n| match n {
                                ENode::Variable(v) => v == "i",
                                _ => false,
                            })
                        });
                        if is_a_i && is_b_two {
                            unions.push((class_id, neg_one_id));
                        }

                        // sqrt(x)^2 -> x
                        if is_b_two {
                            if let Some(c) = egraph.get_class(ra) {
                                for n in &c.nodes {
                                    if let ENode::Sqrt(x) = n {
                                        unions.push((class_id, *x));
                                    }
                                }
                            }
                        }
                    }
                    ENode::Neg(a) => {
                        let ra = egraph.find_immut(*a);
                        let is_a_zero = ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero());
                        if is_a_zero {
                            unions.push((class_id, r_zero));
                        } else if let Some(target_class) = egraph.get_class(ra) {
                            for target_node in &target_class.nodes {
                                if let ENode::Neg(inner) = target_node {
                                    // -(-x) => x
                                    unions.push((class_id, *inner));
                                } else if let ENode::Mul(u, v) = target_node {
                                    // -(x * -1) => x
                                    if egraph.find_immut(*v) == r_neg_one || egraph.get_class(*v).is_some_and(|c| c.is_neg_one()) {
                                        unions.push((class_id, *u));
                                    } else if egraph.find_immut(*u) == r_neg_one || egraph.get_class(*u).is_some_and(|c| c.is_neg_one()) {
                                        unions.push((class_id, *v));
                                    }
                                }
                            }
                        }
                    }
                    ENode::Exp(a) => {
                        let ra = egraph.find_immut(*a);
                        if ra == r_zero || egraph.get_class(ra).is_some_and(|c| c.is_zero()) {
                            unions.push((class_id, r_one));
                        }
                    }
                    ENode::Log(a) => {
                        let ra = egraph.find_immut(*a);
                        if ra == r_one || egraph.get_class(ra).is_some_and(|c| c.is_one()) {
                            unions.push((class_id, r_zero));
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut applied = 0;
        for (class_id, u, v) in neg_muls {
            let mul = egraph.add_node(ENode::Mul(u, v));
            let neg = egraph.add_node(ENode::Neg(mul));
            if egraph.union(class_id, neg) {
                applied += 1;
            }
        }
        for (class_id, node) in new_nodes {
            let new_id = egraph.add_node(node);
            if egraph.union(class_id, new_id) {
                applied += 1;
            }
        }
        for (class_id, prod_node, var_id) in assoc_muls {
            let n_id = egraph.add_node(prod_node);
            let mul_node = egraph.add_node(ENode::Mul(n_id, var_id));
            if egraph.union(class_id, mul_node) {
                applied += 1;
            }
        }
        for (class_id, b, a) in commutative_muls {
            let comm = egraph.add_node(ENode::Mul(b, a));
            if egraph.union(class_id, comm) {
                applied += 1;
            }
        }
        for (class_id, inner_node, outer_id) in deferred_nested_muls {
            let inner_id = egraph.add_node(inner_node);
            let outer = egraph.add_node(ENode::Mul(outer_id, inner_id));
            if egraph.union(class_id, outer) {
                applied += 1;
            }
        }
        for (class_id, base_id, exp_node) in deferred_powers {
            let exp_id = egraph.add_node(exp_node);
            let power = egraph.add_node(ENode::Power(base_id, exp_id));
            if egraph.union(class_id, power) {
                applied += 1;
            }
        }
        for (c1, c2) in unions {
            if egraph.union(c1, c2) {
                applied += 1;
            }
        }
        applied
    }
}
