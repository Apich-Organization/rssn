//! Dimension inference over the term DAG.

use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::RwLock;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::ToPrimitive;
use num_traits::Zero;

use crate::graph::ClassId;
use crate::graph::DimVarId;
use crate::graph::Graph;
use crate::graph::NodeId;
use crate::graph::SymbolId;
use crate::graph::dim::DimArena;
use crate::graph::dim::DimCtx;
use crate::graph::dim::DimError;
use crate::graph::dim::DimRepr;
use crate::graph::dim::DimRule;
use crate::graph::dim::DimTerm;
use crate::graph::dim::Expo;
use crate::graph::dim::UnitResolver;
use crate::graph::dim::exact_exponent;
use crate::graph::op::core;

/// Why two dimensions could not be unified.
#[derive(Clone, Debug)]
pub enum UnifyError {
    /// The dimensions differ.
    Mismatch(DimTerm, DimTerm),
    /// A symbolic exponent would have to satisfy an equation.
    Constraint(DimTerm),
}

/// The substitutions produced while unifying dimension variables.
#[derive(Clone, Debug, Default)]
pub struct Unifier {
    subst: Vec<(DimVarId, DimTerm)>,
}

impl Unifier {
    /// Applies the current substitutions to `term`.
    pub(crate) fn resolve(
        &self,
        term: &DimTerm,
    ) -> DimTerm {
        self.resolve_depth(term, 0)
    }

    fn resolve_depth(
        &self,
        term: &DimTerm,
        depth: usize,
    ) -> DimTerm {
        if depth > 32 {
            return DimTerm::any();
        }
        match &term.repr {
            | DimRepr::Scalar { bases, vars } => {
                let mut out_bases = bases.clone();
                let mut out_vars: Vec<(DimVarId, Expo)> = Vec::new();
                for (id, exponent) in vars {
                    let bound = self
                        .subst
                        .iter()
                        .find(|(variable, _)| variable == id)
                        .map(|(_, value)| value);
                    if let Some(bound) = bound {
                        let value = self.resolve_depth(bound, depth + 1).pow_expo(exponent);
                        merge_term(&mut out_bases, &mut out_vars, &value);
                    } else {
                        out_vars.push((*id, exponent.clone()));
                    }
                }
                DimTerm::scalar(out_bases, out_vars)
            },
            | _ => term.clone(),
        }
    }

    fn bind(
        &mut self,
        id: DimVarId,
        term: &DimTerm,
    ) -> Result<(), UnifyError> {
        let term = self.resolve(term);
        if term_has_var(&term, id, 0) {
            return Err(UnifyError::Constraint(term));
        }
        if let Some(index) = self.subst.iter().position(|(variable, _)| *variable == id) {
            let existing = self.subst.get(index).map(|(_, value)| value.clone());
            if let Some(existing) = existing {
                let merged = self.unify(&existing, &term)?;
                if let Some(slot) = self.subst.get_mut(index) {
                    slot.1 = merged;
                }
                return Ok(());
            }
        }
        self.subst.push((id, term));
        Ok(())
    }

    /// Unifies two dimensions, solving dimension variables when their
    /// exponents are rational constants.
    pub(crate) fn unify(
        &mut self,
        left: &DimTerm,
        right: &DimTerm,
    ) -> Result<DimTerm, UnifyError> {
        let mut left = self.resolve(left);
        let mut right = self.resolve(right);
        for _ in 0..64 {
            match (&left.repr, &right.repr) {
                | (DimRepr::Any, _) => return Ok(right),
                | (_, DimRepr::Any) | (DimRepr::NonMetric, DimRepr::NonMetric) => {
                    return Ok(left);
                },
                | (DimRepr::NonMetric, _) | (_, DimRepr::NonMetric) => {
                    return Err(UnifyError::Mismatch(left, right));
                },
                | (DimRepr::Scalar { .. }, DimRepr::Scalar { .. }) => {},
            }
            let quotient = left.mul(&right.inverse());
            let quotient = self.resolve(&quotient);
            match &quotient.repr {
                | DimRepr::Scalar { bases, vars } => {
                    if let Some((id, exponent)) = solvable(vars) {
                        let rest = without_var(bases, vars, id);
                        let inverse = BigRational::one() / exponent;
                        let solution = rest.pow(&-inverse);
                        self.bind(id, &solution)?;
                        left = self.resolve(&left);
                        right = self.resolve(&right);
                        continue;
                    }
                    if !vars.is_empty() {
                        return Err(UnifyError::Constraint(quotient));
                    }
                    if bases.is_empty() {
                        return Ok(self.resolve(&left));
                    }
                    return Err(if bases.iter().any(|(_, e)| e.is_constant().is_none()) {
                        UnifyError::Constraint(quotient)
                    } else {
                        UnifyError::Mismatch(left, right)
                    });
                },
                | _ => return Ok(left),
            }
        }
        Err(UnifyError::Constraint(left))
    }
}

fn solvable(vars: &[(DimVarId, Expo)]) -> Option<(DimVarId, BigRational)> {
    for (id, exponent) in vars {
        if let Some(rational) = exponent.is_constant()
            && !rational.is_zero()
        {
            return Some((*id, rational));
        }
    }
    None
}

fn without_var(
    bases: &[(crate::graph::BaseId, Expo)],
    vars: &[(DimVarId, Expo)],
    id: DimVarId,
) -> DimTerm {
    let vars = vars
        .iter()
        .filter(|(other, _)| *other != id)
        .cloned()
        .collect();
    DimTerm::scalar(bases.to_vec(), vars)
}

fn term_has_var(
    term: &DimTerm,
    id: DimVarId,
    depth: usize,
) -> bool {
    if depth > 32 {
        return true;
    }
    match &term.repr {
        | DimRepr::Scalar { vars, .. } => vars.iter().any(|(other, _)| *other == id),
        | _ => false,
    }
}

fn merge_term(
    bases: &mut Vec<(crate::graph::BaseId, Expo)>,
    vars: &mut Vec<(DimVarId, Expo)>,
    term: &DimTerm,
) {
    if let DimRepr::Scalar {
        bases: other_bases,
        vars: other_vars,
    } = &term.repr
    {
        for (id, exponent) in other_bases {
            merge_base(bases, *id, exponent);
        }
        for (id, exponent) in other_vars {
            merge_var(vars, *id, exponent);
        }
    }
}

fn merge_base(
    factors: &mut Vec<(crate::graph::BaseId, Expo)>,
    id: crate::graph::BaseId,
    exponent: &Expo,
) {
    if let Some(index) = factors.iter().position(|(other, _)| *other == id)
        && let Some(slot) = factors.get_mut(index)
    {
        slot.1 = slot.1.add(exponent);
    } else {
        factors.push((id, exponent.clone()));
    }
}

fn merge_var(
    factors: &mut Vec<(DimVarId, Expo)>,
    id: DimVarId,
    exponent: &Expo,
) {
    if let Some(index) = factors.iter().position(|(other, _)| *other == id)
        && let Some(slot) = factors.get_mut(index)
    {
        slot.1 = slot.1.add(exponent);
    } else {
        factors.push((id, exponent.clone()));
    }
}

/// The inference engine. It borrows the graph's arena through a lock so the
/// graph itself stays `Sync`.
pub struct DimInfer<'a> {
    graph: &'a Graph,
    arena: &'a RwLock<DimArena>,
    assumptions: &'a HashMap<SymbolId, DimTerm>,
    unifier: Unifier,
    memo: HashMap<ClassId, DimTerm>,
    visiting: HashSet<ClassId>,
    errors: Vec<DimError>,
    unannotated: HashSet<String>,
}

impl<'a> DimInfer<'a> {
    pub(crate) fn new(graph: &'a Graph) -> Self {
        Self {
            graph,
            arena: graph.dims.as_ref(),
            assumptions: &graph.dim_assumptions,
            unifier: Unifier::default(),
            memo: HashMap::new(),
            visiting: HashSet::new(),
            errors: Vec::new(),
            unannotated: HashSet::new(),
        }
    }

    pub(crate) fn run(
        &mut self,
        node: NodeId,
    ) -> DimTerm {
        let class = self.graph.find(node);
        let term = self.class_dim(class);
        self.unifier.resolve(&term)
    }

    pub(crate) fn into_errors(self) -> Vec<DimError> {
        self.errors
    }

    fn class_dim(
        &mut self,
        class: ClassId,
    ) -> DimTerm {
        if let Some(cached) = self.memo.get(&class) {
            return cached.clone();
        }
        if !self.visiting.insert(class) {
            return DimTerm::any();
        }
        let members: Vec<NodeId> = self.graph.members(class).collect();
        let mut result: Option<DimTerm> = None;
        for member in members {
            let dimension = self.derive(member);
            match result.take() {
                | None => result = Some(dimension),
                | Some(current) => {
                    let left = self.unifier.resolve(&current);
                    match self.unifier.unify(&left, &dimension) {
                        | Ok(unified) => result = Some(unified),
                        | Err(error) => {
                            self.record_unify(error, member);
                            result = Some(current);
                        },
                    }
                },
            }
        }
        self.visiting.remove(&class);
        let dimension = result.unwrap_or_else(DimTerm::any);
        let resolved = self.unifier.resolve(&dimension);
        self.memo.insert(class, resolved.clone());
        resolved
    }

    fn derive(
        &mut self,
        node: NodeId,
    ) -> DimTerm {
        if let Some(number) = self.graph.as_number(node) {
            return if number.is_zero() {
                DimTerm::any()
            } else {
                DimTerm::dimensionless()
            };
        }
        if let Some(symbol) = self.graph.as_symbol(node) {
            return self.symbol_dimension(symbol);
        }
        let op = self.graph.op(node);
        let children: Vec<NodeId> = self.graph.children(node).to_vec();
        let args: Vec<DimTerm> = children
            .iter()
            .map(|&child| {
                let class = self.graph.find(child);
                self.class_dim(class)
            })
            .collect();
        let name = self.graph.ops().get(op).name.to_string();
        let rule = self
            .graph
            .ops()
            .attr::<DimRule>(op)
            .copied()
            .or_else(|| crate::graph::dim::default_rule(&name))
            .unwrap_or(DimRule::Unspecified);
        self.apply_rule(rule, &name, node, &children, &args)
    }

    fn symbol_dimension(
        &self,
        symbol: SymbolId,
    ) -> DimTerm {
        if let Some(known) = self.assumptions.get(&symbol) {
            return self.unifier.resolve(known);
        }
        let label = self.graph.interner().symbol_name(symbol).to_string();
        let id = self.intern_var(&label);
        DimTerm::var(id)
    }

    #[allow(clippy::too_many_lines)]
    fn apply_rule(
        &mut self,
        rule: DimRule,
        name: &str,
        node: NodeId,
        children: &[NodeId],
        args: &[DimTerm],
    ) -> DimTerm {
        match rule {
            | DimRule::Unspecified => {
                if self.unannotated.insert(name.to_string()) {
                    self.errors
                        .push(DimError::Unannotated { op: name.to_string() });
                }
                let id = self.intern_var(&format!("{name}(...)"));
                DimTerm::var(id)
            },
            | DimRule::NonMetric => DimTerm::non_metric(),
            | DimRule::Any => DimTerm::any(),
            | DimRule::Same(index) => args.get(index).cloned().unwrap_or_else(DimTerm::any),
            | DimRule::Sum => self.sum(args, node),
            | DimRule::Product => {
                args.iter()
                    .fold(DimTerm::dimensionless(), |acc, arg| acc.mul(arg))
            },
            | DimRule::Quotient => {
                match (args.first(), args.get(1)) {
                    | (Some(left), Some(right)) => left.mul(&right.inverse()),
                    | _ => DimTerm::any(),
                }
            },
            | DimRule::Power => self.power(node, children, args),
            | DimRule::Root(degree) => {
                let exponent = BigRational::new(BigInt::one(), BigInt::from(degree));
                args.first()
                    .map_or_else(DimTerm::any, |base| base.pow(&exponent))
            },
            | DimRule::Dimensionless => {
                for (index, arg) in args.iter().enumerate() {
                    if let Err(error) = self.unifier.unify(arg, &DimTerm::dimensionless()) {
                        self.record_unify(error, node);
                        self.errors.push(DimError::NonDimensionless {
                            op: name.to_string(),
                            argument: children
                                .get(index)
                                .map_or_else(String::new, |&child| self.graph.display(child)),
                            found: self.text(arg),
                        });
                    }
                }
                DimTerm::dimensionless()
            },
            | DimRule::Inverse => args.first().map_or_else(DimTerm::any, DimTerm::inverse),
            | DimRule::Square => {
                args.first().map_or_else(DimTerm::any, |base| {
                    base.pow(&BigRational::from_integer(BigInt::from(2)))
                })
            },
            | DimRule::List => self.list(args, node),
            | DimRule::Compare => {
                self.compare(args, node);
                DimTerm::non_metric()
            },
            | DimRule::Angle => {
                self.compare(args, node);
                DimTerm::dimensionless()
            },
            | DimRule::Diff { f, x } => {
                let function = args.get(f).cloned().unwrap_or_else(DimTerm::any);
                let variable = args.get(x).cloned().unwrap_or_else(DimTerm::any);
                self.check_not_non_metric(&function, node, name);
                function.mul(&variable.inverse())
            },
            | DimRule::DiffN { f, x, n } => {
                let function = args.get(f).cloned().unwrap_or_else(DimTerm::any);
                let variable = args.get(x).cloned().unwrap_or_else(DimTerm::any);
                let degree = children
                    .get(n)
                    .and_then(|&child| self.expo_of(child))
                    .unwrap_or_else(Expo::one);
                self.check_not_non_metric(&function, node, name);
                function.mul(&variable.pow_expo(&degree).inverse())
            },
            | DimRule::Integral { f, x } => self.integral(name, node, children, args, f, x),
            | DimRule::SumOver { f, index } => {
                self.bind_bounds(name, node, children, args, index);
                args.get(f).cloned().unwrap_or_else(DimTerm::any)
            },
            | DimRule::ProductOver { f, index } => {
                self.bind_bounds(name, node, children, args, index);
                let factor = args.get(f).cloned().unwrap_or_else(DimTerm::any);
                let count = self.count_exponent(children, index);
                factor.pow_expo(&count)
            },
            | DimRule::Quantity => self.quantity(node),
            | DimRule::Custom(custom) => {
                let result = {
                    let arena = self.read_arena();
                    custom(
                        &DimCtx {
                            graph: self.graph,
                            arena: &arena,
                            node,
                        },
                        args,
                    )
                };
                match result {
                    | Ok(dimension) => dimension,
                    | Err(error) => {
                        self.errors.push(error);
                        DimTerm::any()
                    },
                }
            },
            | DimRule::UnknownFunction => {
                let id = self.intern_var(&self.graph.display(node));
                DimTerm::var(id)
            },
        }
    }

    fn sum(
        &mut self,
        args: &[DimTerm],
        node: NodeId,
    ) -> DimTerm {
        let mut result: Option<DimTerm> = None;
        for arg in args {
            match result.take() {
                | None => result = Some(arg.clone()),
                | Some(current) => {
                    let left = self.unifier.resolve(&current);
                    match self.unifier.unify(&left, arg) {
                        | Ok(unified) => result = Some(unified),
                        | Err(error) => {
                            self.record_unify(error, node);
                            result = Some(current);
                        },
                    }
                },
            }
        }
        result.unwrap_or_else(DimTerm::any)
    }

    fn compare(
        &mut self,
        args: &[DimTerm],
        node: NodeId,
    ) {
        let mut result: Option<DimTerm> = None;
        for arg in args {
            match result.take() {
                | None => result = Some(arg.clone()),
                | Some(current) => {
                    let left = self.unifier.resolve(&current);
                    match self.unifier.unify(&left, arg) {
                        | Ok(unified) => result = Some(unified),
                        | Err(error) => {
                            self.record_unify(error, node);
                            result = Some(current);
                        },
                    }
                },
            }
        }
    }

    fn list(
        &mut self,
        args: &[DimTerm],
        node: NodeId,
    ) -> DimTerm {
        if args.is_empty() {
            return DimTerm::non_metric();
        }
        if args.iter().all(DimTerm::is_non_metric) {
            return DimTerm::non_metric();
        }
        let mut result: Option<DimTerm> = None;
        for arg in args {
            if arg.is_non_metric() {
                continue;
            }
            match result.take() {
                | None => result = Some(arg.clone()),
                | Some(current) => {
                    let left = self.unifier.resolve(&current);
                    match self.unifier.unify(&left, arg) {
                        | Ok(unified) => result = Some(unified),
                        | Err(_) => {
                            self.errors.push(DimError::MixedList {
                                context: self.graph.display(node),
                            });
                            return DimTerm::any();
                        },
                    }
                },
            }
        }
        result.unwrap_or_else(DimTerm::any)
    }

    fn power(
        &mut self,
        node: NodeId,
        children: &[NodeId],
        args: &[DimTerm],
    ) -> DimTerm {
        let Some(base) = args.first().cloned() else {
            return DimTerm::any();
        };
        if base.is_non_metric() {
            self.errors.push(DimError::NonDimensionless {
                op: String::from("pow"),
                argument: self.graph.display(node),
                found: self.text(&base),
            });
            return DimTerm::any();
        }
        if let Some(exponent) = args.get(1)
            && let Err(error) = self.unifier.unify(exponent, &DimTerm::dimensionless())
        {
            self.record_unify(error, node);
            self.errors.push(DimError::NonDimensionless {
                op: String::from("pow"),
                argument: self.graph.display(node),
                found: self.text(exponent),
            });
        }
        let exponent = children
            .get(1)
            .and_then(|&child| self.expo_of(child))
            .unwrap_or_else(Expo::one);
        base.pow_expo(&exponent)
    }

    fn integral(
        &mut self,
        name: &str,
        node: NodeId,
        children: &[NodeId],
        args: &[DimTerm],
        f: usize,
        x: usize,
    ) -> DimTerm {
        let function = args.get(f).cloned().unwrap_or_else(DimTerm::any);
        let variable = args.get(x).cloned().unwrap_or_else(DimTerm::any);
        if name == "defint" {
            for limit in [2_usize, 3] {
                if let Some(limit_dim) = args.get(limit)
                    && let Err(error) = self.unifier.unify(limit_dim, &variable)
                {
                    self.record_unify(error, node);
                }
            }
        }
        let _ = children;
        function.mul(&variable)
    }

    fn bind_bounds(
        &mut self,
        _name: &str,
        node: NodeId,
        _children: &[NodeId],
        args: &[DimTerm],
        index: usize,
    ) {
        let Some(index_dim) = args.get(index).cloned() else {
            return;
        };
        for bound in [2_usize, 3] {
            if let Some(bound_dim) = args.get(bound)
                && let Err(error) = self.unifier.unify(bound_dim, &index_dim)
            {
                self.record_unify(error, node);
            }
        }
    }

    fn count_exponent(
        &mut self,
        children: &[NodeId],
        _index: usize,
    ) -> Expo {
        let upper = children.get(3).and_then(|&child| self.expo_of(child));
        let lower = children.get(2).and_then(|&child| self.expo_of(child));
        let (Some(upper), Some(lower)) = (upper, lower) else {
            return Expo::one();
        };
        upper.sub(&lower).add(&Expo::one())
    }

    fn quantity(
        &self,
        node: NodeId,
    ) -> DimTerm {
        let Some(op) = self.graph.ops().lookup("quantity") else {
            return DimTerm::any();
        };
        let Some(resolver) = self.graph.ops().attr::<UnitResolver>(op) else {
            return DimTerm::any();
        };
        let arena = self.read_arena();
        (resolver.0)(self.graph, node, &arena).map_or_else(DimTerm::any, |info| info.dimension)
    }

    fn check_not_non_metric(
        &mut self,
        dimension: &DimTerm,
        node: NodeId,
        name: &str,
    ) {
        if dimension.is_non_metric() {
            self.errors.push(DimError::NonDimensionless {
                op: name.to_string(),
                argument: self.graph.display(node),
                found: self.text(dimension),
            });
        }
    }

    fn expo_of(
        &mut self,
        node: NodeId,
    ) -> Option<Expo> {
        if let Some(number) = self.graph.as_number(node) {
            return match exact_exponent(number) {
                | Some(value) => Some(Expo::from_rat(value)),
                | None => Some(self.atom_of(node)),
            };
        }
        if let Some(symbol) = self.graph.as_symbol(node) {
            let label = self.graph.interner().symbol_name(symbol).to_string();
            let atom = self.intern_atom(&label);
            return Some(Expo::atom(atom));
        }
        let op = self.graph.op(node);
        if op == core::ADD {
            let children: Vec<NodeId> = self.graph.children(node).to_vec();
            let mut total = Expo::zero();
            for child in children {
                total = total.add(&self.expo_of(child)?);
            }
            return Some(total);
        }
        if op == core::MUL {
            let children: Vec<NodeId> = self.graph.children(node).to_vec();
            let mut product = Expo::one();
            for child in children {
                product = product.mul(&self.expo_of(child)?);
            }
            return Some(product);
        }
        if op == core::POW {
            let children: Vec<NodeId> = self.graph.children(node).to_vec();
            let base = children.first().and_then(|&child| self.expo_of(child))?;
            let exponent = children.get(1).and_then(|&child| self.expo_of(child))?;
            if let Some(integer) = exponent.as_integer()
                && let Some(k) = integer.to_i32()
                && let Some(value) = base.pow_i32(k)
            {
                return Some(value);
            }
            return Some(self.atom_of(node));
        }
        Some(self.atom_of(node))
    }

    fn atom_of(
        &self,
        node: NodeId,
    ) -> Expo {
        let label = self.graph.display(node);
        let atom = self.intern_atom(&label);
        Expo::atom(atom)
    }

    fn intern_atom(
        &self,
        label: &str,
    ) -> crate::graph::dim::AtomId {
        let mut arena = match self.arena.write() {
            | Ok(guard) => guard,
            | Err(poisoned) => poisoned.into_inner(),
        };
        arena.intern_atom(label)
    }

    fn intern_var(
        &self,
        label: &str,
    ) -> DimVarId {
        let mut arena = match self.arena.write() {
            | Ok(guard) => guard,
            | Err(poisoned) => poisoned.into_inner(),
        };
        arena.intern_var(label)
    }

    fn read_arena(&self) -> std::sync::RwLockReadGuard<'_, DimArena> {
        match self.arena.read() {
            | Ok(guard) => guard,
            | Err(poisoned) => poisoned.into_inner(),
        }
    }

    fn text(
        &self,
        dimension: &DimTerm,
    ) -> String {
        let arena = self.read_arena();
        dimension.describe(&arena)
    }

    fn record_unify(
        &mut self,
        error: UnifyError,
        node: NodeId,
    ) {
        if self.errors.len() >= 64 {
            return;
        }
        let context = self.graph.display(node);
        match error {
            | UnifyError::Mismatch(left, right) => {
                let left = self.text(&left);
                let right = self.text(&right);
                self.errors
                    .push(DimError::Mismatch { left, right, context });
            },
            | UnifyError::Constraint(quotient) => {
                let constraint = format!("{} = 1", self.text(&quotient));
                self.errors
                    .push(DimError::SymbolicConstraint { constraint, context });
            },
        }
    }
}
