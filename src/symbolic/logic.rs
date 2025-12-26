//! # Symbolic Logic Module
//!
//! This module provides functions for symbolic manipulation of logical expressions.
//! It includes capabilities for simplifying logical formulas, converting them to
//! normal forms (CNF, DNF), and a basic SAT solver for quantifier-free predicate logic.

use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Arc;

use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;

/// Checks if a variable occurs freely in an expression.

pub(crate) fn free_vars(
    expr: &Expr,
    free: &mut BTreeSet<String>,
    bound: &mut BTreeSet<String>,
) {

    match expr {
        | Expr::Dag(node) => {

            free_vars(
                &node
                    .to_expr()
                    .expect(
                        "Free Vars",
                    ),
                free,
                bound,
            );
        },
        | Expr::Variable(s) => {

            if !bound.contains(s) {

                free.insert(s.clone());
            }
        },
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b)
        | Expr::Eq(a, b)
        | Expr::Lt(a, b)
        | Expr::Gt(a, b)
        | Expr::Le(a, b)
        | Expr::Ge(a, b)
        | Expr::Xor(a, b)
        | Expr::Implies(a, b)
        | Expr::Equivalent(a, b) => {

            free_vars(a, free, bound);

            free_vars(b, free, bound);
        },
        | Expr::Neg(a)
        | Expr::Not(a) => {

            free_vars(a, free, bound);
        },
        | Expr::And(v)
        | Expr::Or(v) => {

            for sub_expr in v {

                free_vars(
                    sub_expr,
                    free,
                    bound,
                );
            }
        },
        | Expr::ForAll(var, body)
        | Expr::Exists(var, body) => {

            bound.insert(var.clone());

            free_vars(
                body, free, bound,
            );

            bound.remove(var);
        },
        | Expr::Predicate {
            args,
            ..
        } => {
            for arg in args {

                free_vars(
                    arg, free, bound,
                );
            }
        },
        | _ => {},
    }
}

/// Helper to check if an expression contains a specific free variable.

pub(crate) fn has_free_var(
    expr: &Expr,
    var: &str,
) -> bool {

    let mut free = BTreeSet::new();

    let mut bound = BTreeSet::new();

    free_vars(
        expr,
        &mut free,
        &mut bound,
    );

    free.contains(var)
}

/// Simplifies a logical expression by applying a set of transformation rules.
///
/// This function recursively traverses the expression tree and applies rules such as:
/// - **Double Negation**: `Not(Not(P))` -> `P`
/// - **De Morgan's Laws**: `Not(ForAll(x, P(x)))` -> `Exists(x, Not(P(x)))`
/// - **Constant Folding**: `A And False` -> `False`, `A Or True` -> `True`
/// - **Identity and Idempotence**: `A And True` -> `A`, `A Or A` -> `A`
/// - **Contradiction/Tautology**: `A And Not(A)` -> `False`, `A Or Not(A)` -> `True`
/// - **Quantifier Reduction**: Removes redundant quantifiers where the variable is not free in the body.
/// - **Quantifier Pushing**: Moves quantifiers inwards to narrow their scope, e.g.,
///   `ForAll(x, P(x) And Q(y))` -> `(ForAll(x, P(x))) And Q(y)`.
///
/// # Arguments
/// * `expr` - The logical expression to simplify.
///
/// # Returns
/// A new, simplified logical expression.
#[must_use]

pub fn simplify_logic(
    expr: &Expr
) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            simplify_logic(
                &node
                    .to_expr()
                    .expect("Simplify Logic"),
            )
        },
        | Expr::Not(inner) => {
            match simplify_logic(inner) {
                | Expr::Boolean(b) => Expr::Boolean(!b),
                | Expr::Not(sub) => (*sub).clone(),
                | Expr::ForAll(var, body) => {
                    Expr::Exists(
                        var,
                        Arc::new(simplify_logic(
                            &Expr::new_not(body),
                        )),
                    )
                },
                | Expr::Exists(var, body) => {
                    Expr::ForAll(
                        var,
                        Arc::new(simplify_logic(
                            &Expr::new_not(body),
                        )),
                    )
                },
                | simplified_inner => Expr::new_not(simplified_inner),
            }
        },
        | Expr::And(v) => {

            let mut new_terms = Vec::new();

            for term in v {

                let simplified = simplify_logic(term);

                match simplified {
                    | Expr::Boolean(false) => return Expr::Boolean(false),
                    | Expr::Boolean(true) => continue,
                    | Expr::And(mut sub_terms) => new_terms.append(&mut sub_terms),
                    | _ => new_terms.push(simplified),
                }
            }

            let mut unique_terms = BTreeSet::new();

            for term in new_terms {

                unique_terms.insert(term);
            }

            for term in &unique_terms {

                if unique_terms.contains(&Expr::new_not(
                    term.clone(),
                )) {

                    return Expr::Boolean(false);
                }
            }

            if unique_terms.is_empty() {

                Expr::Boolean(true)
            } else if unique_terms.len() == 1 {

                unique_terms
                    .into_iter()
                    .next()
                    .expect("Unique Term Parsing Failed")
            } else {

                Expr::And(
                    unique_terms
                        .into_iter()
                        .collect(),
                )
            }
        },
        | Expr::Or(v) => {

            let mut new_terms = Vec::new();

            for term in v {

                let simplified = simplify_logic(term);

                match simplified {
                    | Expr::Boolean(true) => return Expr::Boolean(true),
                    | Expr::Boolean(false) => continue,
                    | Expr::Or(mut sub_terms) => new_terms.append(&mut sub_terms),
                    | _ => new_terms.push(simplified),
                }
            }

            let mut unique_terms = BTreeSet::new();

            for term in new_terms {

                unique_terms.insert(term);
            }

            for term in &unique_terms {

                if unique_terms.contains(&Expr::new_not(
                    term.clone(),
                )) {

                    return Expr::Boolean(true);
                }
            }

            if unique_terms.is_empty() {

                Expr::Boolean(false)
            } else if unique_terms.len() == 1 {

                unique_terms
                    .into_iter()
                    .next()
                    .expect("Unique Term Parsing Failed")
            } else {

                Expr::Or(
                    unique_terms
                        .into_iter()
                        .collect(),
                )
            }
        },
        | Expr::Implies(a, b) => {
            simplify_logic(&Expr::Or(vec![
                Expr::Not(Arc::new(
                    a.as_ref().clone(),
                )),
                b.as_ref().clone(),
            ]))
        },
        | Expr::Equivalent(a, b) => {
            simplify_logic(&Expr::And(vec![
                Expr::Implies(a.clone(), b.clone()),
                Expr::Implies(b.clone(), a.clone()),
            ]))
        },
        | Expr::Xor(a, b) => {
            simplify_logic(&Expr::And(vec![
                Expr::Or(vec![
                    a.as_ref().clone(),
                    b.as_ref().clone(),
                ]),
                Expr::Not(Arc::new(Expr::And(
                    vec![
                        a.as_ref().clone(),
                        b.as_ref().clone(),
                    ],
                ))),
            ]))
        },
        | Expr::ForAll(var, body) => {

            let simplified_body = simplify_logic(body);

            if !has_free_var(
                &simplified_body,
                var,
            ) {

                return simplified_body;
            }

            if let Expr::And(terms) = &simplified_body {

                let mut with_var = vec![];

                let mut without_var = vec![];

                for term in terms {

                    if has_free_var(term, var) {

                        with_var.push(term.clone());
                    } else {

                        without_var.push(term.clone());
                    }
                }

                if !without_var.is_empty() {

                    let forall_part = if with_var.is_empty() {

                        Expr::Boolean(true)
                    } else {

                        Expr::ForAll(
                            var.clone(),
                            Arc::new(Expr::And(with_var)),
                        )
                    };

                    without_var.push(simplify_logic(
                        &forall_part,
                    ));

                    return simplify_logic(&Expr::And(
                        without_var,
                    ));
                }
            }

            Expr::ForAll(
                var.clone(),
                Arc::new(simplified_body),
            )
        },
        | Expr::Exists(var, body) => {

            let simplified_body = simplify_logic(body);

            if !has_free_var(
                &simplified_body,
                var,
            ) {

                return simplified_body;
            }

            if let Expr::Or(terms) = &simplified_body {

                let mut with_var = vec![];

                let mut without_var = vec![];

                for term in terms {

                    if has_free_var(term, var) {

                        with_var.push(term.clone());
                    } else {

                        without_var.push(term.clone());
                    }
                }

                if !without_var.is_empty() {

                    let exists_part = if with_var.is_empty() {

                        Expr::Boolean(false)
                    } else {

                        Expr::Exists(
                            var.clone(),
                            Arc::new(Expr::Or(with_var)),
                        )
                    };

                    without_var.push(simplify_logic(
                        &exists_part,
                    ));

                    return simplify_logic(&Expr::Or(
                        without_var,
                    ));
                }
            }

            Expr::Exists(
                var.clone(),
                Arc::new(simplified_body),
            )
        },
        | Expr::Predicate { name, args } => {
            Expr::Predicate {
                name: name.clone(),
                args: args
                    .iter()
                    .map(|expr: &Expr| simplify(&expr.clone()))
                    .collect(),
            }
        },
        | _ => expr.clone(),
    }
}

pub(crate) fn to_basic_logic_ops(
    expr: &Expr
) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            to_basic_logic_ops(
                &node
                    .to_expr()
                    .expect("To Basic Logic Ops"),
            )
        },
        | Expr::Implies(a, b) => {
            Expr::Or(vec![
                Expr::Not(Arc::new(
                    to_basic_logic_ops(a),
                )),
                to_basic_logic_ops(b),
            ])
        },
        | Expr::Equivalent(a, b) => {
            Expr::And(vec![
                Expr::Or(vec![
                    Expr::Not(Arc::new(
                        to_basic_logic_ops(a),
                    )),
                    to_basic_logic_ops(b),
                ]),
                Expr::Or(vec![
                    Expr::Not(Arc::new(
                        to_basic_logic_ops(b),
                    )),
                    to_basic_logic_ops(a),
                ]),
            ])
        },
        | Expr::Xor(a, b) => {
            Expr::And(vec![
                Expr::Or(vec![
                    to_basic_logic_ops(a),
                    to_basic_logic_ops(b),
                ]),
                Expr::Not(Arc::new(Expr::And(
                    vec![
                        to_basic_logic_ops(a),
                        to_basic_logic_ops(b),
                    ],
                ))),
            ])
        },
        | Expr::And(v) => {
            Expr::And(
                v.iter()
                    .map(to_basic_logic_ops)
                    .collect(),
            )
        },
        | Expr::Or(v) => {
            Expr::Or(
                v.iter()
                    .map(to_basic_logic_ops)
                    .collect(),
            )
        },
        | Expr::Not(a) => {
            Expr::new_not(to_basic_logic_ops(
                a,
            ))
        },
        | _ => expr.clone(),
    }
}

pub(crate) fn move_not_inwards(
    expr: &Expr
) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            move_not_inwards(
                &node
                    .to_expr()
                    .expect("Move not Inwards"),
            )
        },
        | Expr::Not(a) => {
            match &**a {
                | Expr::And(v) => {
                    Expr::Or(
                        v.iter()
                            .map(|e| {

                                move_not_inwards(&Expr::new_not(
                                    e.clone(),
                                ))
                            })
                            .collect(),
                    )
                },
                | Expr::Or(v) => {
                    Expr::And(
                        v.iter()
                            .map(|e| {

                                move_not_inwards(&Expr::new_not(
                                    e.clone(),
                                ))
                            })
                            .collect(),
                    )
                },
                | Expr::Not(b) => move_not_inwards(b),
                | Expr::ForAll(var, body) => {
                    Expr::Exists(
                        var.clone(),
                        Arc::new(move_not_inwards(
                            &Expr::new_not(body.clone()),
                        )),
                    )
                },
                | Expr::Exists(var, body) => {
                    Expr::ForAll(
                        var.clone(),
                        Arc::new(move_not_inwards(
                            &Expr::new_not(body.clone()),
                        )),
                    )
                },
                | _ => expr.clone(),
            }
        },
        | Expr::And(v) => {
            Expr::And(
                v.iter()
                    .map(move_not_inwards)
                    .collect(),
            )
        },
        | Expr::Or(v) => {
            Expr::Or(
                v.iter()
                    .map(move_not_inwards)
                    .collect(),
            )
        },
        | _ => expr.clone(),
    }
}

pub(crate) fn distribute_or_over_and(
    expr: &Expr
) -> Expr {

    match expr {
        | Expr::Dag(node) => {
            distribute_or_over_and(
                &node
                    .to_expr()
                    .expect("Distribute or Over"),
            )
        },
        | Expr::Or(v) => {

            let v_dist: Vec<Expr> = v
                .iter()
                .map(distribute_or_over_and)
                .collect();

            if let Some(pos) = v_dist
                .iter()
                .position(|e| matches!(e, Expr::And(_)))
            {

                let and_clause = v_dist[pos].clone();

                let other_terms: Vec<Expr> = v_dist
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != pos)
                    .map(|(_, e)| e.clone())
                    .collect();

                if let Expr::And(and_terms) = and_clause {

                    let new_clauses = and_terms
                        .iter()
                        .map(|term| {

                            let mut new_or_list = other_terms.clone();

                            new_or_list.push(term.clone());

                            distribute_or_over_and(&Expr::Or(
                                new_or_list,
                            ))
                        })
                        .collect();

                    return Expr::And(new_clauses);
                }
            }

            Expr::Or(v_dist)
        },
        | Expr::And(v) => {
            Expr::And(
                v.iter()
                    .map(distribute_or_over_and)
                    .collect(),
            )
        },
        | _ => expr.clone(),
    }
}

/// Converts a logical expression into Conjunctive Normal Form (CNF).
///
/// CNF is a standardized representation of a logical formula which is a conjunction
/// of one or more clauses, where each clause is a disjunction of literals.
/// The conversion process involves three main steps:
/// 1.  Eliminating complex logical operators like `Implies`, `Equivalent`, and `Xor`.
/// 2.  Moving all `Not` operators inwards using De Morgan's laws.
/// 3.  Distributing `Or` over `And` to achieve the final CNF structure.
///
/// # Arguments
/// * `expr` - The logical expression to convert.
///
/// # Returns
/// An equivalent expression in Conjunctive Normal Form.
#[must_use]

pub fn to_cnf(expr: &Expr) -> Expr {

    let simplified =
        simplify_logic(expr);

    let basic_ops =
        to_basic_logic_ops(&simplified);

    let not_inwards =
        move_not_inwards(&basic_ops);

    let distributed =
        distribute_or_over_and(
            &not_inwards,
        );

    simplify_logic(&distributed)
}

/// Converts a logical expression into Disjunctive Normal Form (DNF).
///
/// DNF is a standardized representation of a logical formula which is a disjunction
/// of one or more clauses, where each clause is a conjunction of literals.
/// This implementation cleverly achieves the conversion by using the `to_cnf` function:
/// 1.  The input expression `expr` is negated: `Not(expr)`.
/// 2.  The negated expression is converted to CNF: `cnf(Not(expr))`.
/// 3.  The resulting CNF is negated again, and De Morgan's laws are applied implicitly
///     by `simplify_logic`, resulting in the DNF of the original expression.
///
/// # Arguments
/// * `expr` - The logical expression to convert.
///
/// # Returns
/// An equivalent expression in Disjunctive Normal Form.
#[must_use]

pub fn to_dnf(expr: &Expr) -> Expr {

    let not_expr = simplify_logic(
        &Expr::new_not(expr.clone()),
    );

    let cnf_of_not = to_cnf(&not_expr);

    simplify_logic(&Expr::new_not(
        cnf_of_not,
    ))
}

/// Determines if a quantifier-free logical formula is satisfiable using the DPLL algorithm.
///
/// This function first checks if the expression contains any quantifiers (`ForAll`, `Exists`).
/// If it does, the problem is generally undecidable, and the function returns `None`.
///
/// For quantifier-free formulas, it proceeds by:
/// 1.  Converting the expression to Conjunctive Normal Form (CNF).
/// 2.  Applying the recursive DPLL (Davis-Putnam-Logemann-Loveland) algorithm to the CNF clauses.
///
/// The DPLL algorithm attempts to find a satisfying assignment for the propositional variables
/// (in this case, predicate instances like `P(x)`) by using unit propagation, pure literal
/// elimination (implicitly), and recursive branching on variable assignments.
///
/// # Arguments
/// * `expr` - A logical expression, which should be quantifier-free for a definitive result.
///
/// # Returns
/// * `Some(true)` if the formula is satisfiable.
/// * `Some(false)` if the formula is unsatisfiable.
/// * `None` if the formula contains quantifiers, as this solver does not handle them.
#[must_use]

pub fn is_satisfiable(
    expr: &Expr
) -> Option<bool> {

    if contains_quantifier(expr) {

        return None;
    }

    let cnf = to_cnf(expr);

    if let Expr::Boolean(b) = cnf {

        return Some(b);
    }

    let mut clauses =
        extract_clauses(&cnf);

    let mut assignments =
        HashMap::new();

    Some(dpll(
        &mut clauses,
        &mut assignments,
    ))
}

/// A literal is an atomic proposition (e.g., P(x)) or its negation.
#[derive(
    Debug,
    Clone,
    PartialEq,
    Eq,
    Hash,
    PartialOrd,
    Ord,
)]

pub enum Literal {
    Positive(Expr),
    Negative(Expr),
}

pub(crate) const fn get_atom(
    literal: &Literal
) -> &Expr {

    match literal {
        | Literal::Positive(atom) => {
            atom
        },
        | Literal::Negative(atom) => {
            atom
        },
    }
}

pub(crate) fn extract_clauses(
    cnf_expr: &Expr
) -> Vec<HashSet<Literal>> {

    let mut clauses = Vec::new();

    if let Expr::And(conjuncts) =
        cnf_expr
    {

        for clause_expr in conjuncts {

            clauses.push(extract_literals_from_clause(clause_expr));
        }
    } else {

        clauses.push(extract_literals_from_clause(cnf_expr));
    }

    clauses
}

pub(crate) fn extract_literals_from_clause(
    clause_expr: &Expr
) -> HashSet<Literal> {

    let mut literals = HashSet::new();

    if let Expr::Or(disjuncts) =
        clause_expr
    {

        for literal_expr in disjuncts {

            if let Expr::Not(atom) =
                literal_expr
            {

                literals.insert(
                    Literal::Negative(
                        atom.as_ref()
                            .clone(),
                    ),
                );
            } else {

                literals.insert(
                    Literal::Positive(
                        literal_expr
                            .clone(),
                    ),
                );
            }
        }
    } else if let Expr::Not(atom) =
        clause_expr
    {

        literals.insert(
            Literal::Negative(
                atom.as_ref()
                    .clone(),
            ),
        );
    } else {

        literals.insert(
            Literal::Positive(
                clause_expr.clone(),
            ),
        );
    }

    literals
}

pub(crate) fn dpll(
    clauses: &mut Vec<HashSet<Literal>>,
    assignments: &mut HashMap<
        Expr,
        bool,
    >,
) -> bool {

    while let Some(unit_literal) =
        find_unit_clause(clauses)
    {

        let (atom, value) =
            match unit_literal {
                | Literal::Positive(
                    a,
                ) => (a, true),
                | Literal::Negative(
                    a,
                ) => (a, false),
            };

        assignments.insert(
            atom.clone(),
            value,
        );

        simplify_clauses(
            clauses,
            &atom,
            value,
        );

        if clauses.is_empty() {

            return true;
        }

        if clauses
            .iter()
            .any(HashSet::is_empty)
        {

            return false;
        }
    }

    if clauses.is_empty() {

        return true;
    }

    let atom_to_branch =
        match get_unassigned_atom(
            clauses,
            assignments,
        ) {
            | Some(v) => v,
            | _none => return true,
        };

    let mut clauses_true =
        clauses.clone();

    let mut assignments_true =
        assignments.clone();

    assignments_true.insert(
        atom_to_branch.clone(),
        true,
    );

    simplify_clauses(
        &mut clauses_true,
        &atom_to_branch,
        true,
    );

    if dpll(
        &mut clauses_true,
        &mut assignments_true,
    ) {

        return true;
    }

    let mut clauses_false =
        clauses.clone();

    let mut assignments_false =
        assignments.clone();

    assignments_false.insert(
        atom_to_branch.clone(),
        false,
    );

    simplify_clauses(
        &mut clauses_false,
        &atom_to_branch,
        false,
    );

    if dpll(
        &mut clauses_false,
        &mut assignments_false,
    ) {

        return true;
    }

    false
}

pub(crate) fn find_unit_clause(
    clauses: &[HashSet<Literal>]
) -> Option<Literal> {

    clauses
        .iter()
        .find(|c| c.len() == 1)
        .and_then(|c| {

            c.iter()
                .next()
                .cloned()
        })
}

pub(crate) fn simplify_clauses(
    clauses: &mut Vec<HashSet<Literal>>,
    atom: &Expr,
    value: bool,
) {

    clauses.retain(|clause| {

        !clause
            .iter()
            .any(|lit| {

                match lit {
                    | Literal::Positive(a) => a == atom && value,
                    | Literal::Negative(a) => a == atom && !value,
                }
            })
    });

    let opposite_literal = if value {

        Literal::Negative(atom.clone())
    } else {

        Literal::Positive(atom.clone())
    };

    for clause in clauses {

        clause
            .remove(&opposite_literal);
    }
}

pub(crate) fn get_unassigned_atom(
    clauses: &[HashSet<Literal>],
    assignments: &HashMap<Expr, bool>,
) -> Option<Expr> {

    for clause in clauses {

        for literal in clause {

            let atom =
                get_atom(literal);

            if !assignments
                .contains_key(atom)
            {

                return Some(
                    atom.clone(),
                );
            }
        }
    }

    None
}

pub(crate) fn contains_quantifier(
    expr: &Expr
) -> bool {

    match expr {
        | Expr::Dag(node) => {
            contains_quantifier(
                &node
                    .to_expr()
                    .expect("Contains Quantifier"),
            )
        },
        | Expr::ForAll(_, _) | Expr::Exists(_, _) => true,
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b)
        | Expr::Eq(a, b)
        | Expr::Lt(a, b)
        | Expr::Gt(a, b)
        | Expr::Le(a, b)
        | Expr::Ge(a, b)
        | Expr::Xor(a, b)
        | Expr::Implies(a, b)
        | Expr::Equivalent(a, b) => contains_quantifier(a) || contains_quantifier(b),
        | Expr::Neg(a) | Expr::Not(a) => contains_quantifier(a),
        | Expr::And(v) | Expr::Or(v) => {
            v.iter()
                .any(contains_quantifier)
        },
        | Expr::Predicate { args, .. } => {
            args.iter()
                .any(contains_quantifier)
        },
        | _ => false,
    }
}
