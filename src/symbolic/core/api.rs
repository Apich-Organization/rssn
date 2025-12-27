#![allow(deprecated)]

use std::cmp::Ordering;
use std::collections::hash_map::Entry;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::convert::AsRef;
use std::fmt::Debug;
use std::fmt::Write;
use std::fmt::{
    self,
};
use std::hash::Hash;
use std::hash::Hasher;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::RwLock;

use lazy_static::lazy_static;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use ordered_float::OrderedFloat;

use super::ast_impl::*;
use super::dag_mgr::*;
use super::expr::*;
use super::expr_impl::*;
use super::to_expr::*;
use crate::symbolic::unit_unification::UnitQuantity;

impl AsRef<Self> for Expr {
    fn as_ref(&self) -> &Self {

        self
    }
}

// --- Helper Macros ---
macro_rules! unary_constructor {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.

        pub fn $name<A>(a : A) -> Expr
        where
            A : AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    vec![dag_a],
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

macro_rules! binary_constructor {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.

        pub fn $name<A, B>(
            a : A,
            b : B,
        ) -> Expr
        where
            A : AsRef<Expr>,
            B : AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let dag_b = DAG_MANAGER
                .get_or_create(b.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    vec![dag_a, dag_b],
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

macro_rules! n_ary_constructor {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.

        pub fn $name<I, T>(elements : I) -> Expr
        where
            I : IntoIterator<Item = T>,
            T : AsRef<Expr>,
        {

            let children_nodes = elements
                .into_iter()
                .map(|child| {

                    DAG_MANAGER
                        .get_or_create(child.as_ref())
                        .expect("DAG manager get_or_create failed")
                })
                .collect::<Vec<_>>();

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    children_nodes,
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the \
            'UnaryList' variant \
            instead."
)]

macro_rules! unary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'UnaryList' variant instead."
        )]

        pub fn $name<A>(a : A) -> Expr
        where
            A : AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    vec![dag_a],
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the \
            'BinaryList' variant \
            instead."
)]

macro_rules! binary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'BinaryList' variant instead."
        )]

        pub fn $name<A, B>(
            a : A,
            b : B,
        ) -> Expr
        where
            A : AsRef<Expr>,
            B : AsRef<Expr>,
        {

            let dag_a = DAG_MANAGER
                .get_or_create(a.as_ref())
                .expect("DAG manager get_or_create failed");

            let dag_b = DAG_MANAGER
                .get_or_create(b.as_ref())
                .expect("DAG manager get_or_create failed");

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    vec![dag_a, dag_b],
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

#[deprecated(
    since = "0.1.18",
    note = "Please use the 'NaryList' \
            variant instead."
)]

macro_rules! n_ary_constructor_deprecated {
    ($name:ident, $op:ident) => {
        /// Creates a new
        #[doc = stringify!($op)]
        /// expression, managed by the DAG.
        #[deprecated(
            since = "0.1.18",
            note = "Please use the 'NaryList' variant instead."
        )]

        pub fn $name<I, T>(elements : I) -> Expr
        where
            I : IntoIterator<Item = T>,
            T : AsRef<Expr>,
        {

            let children_nodes = elements
                .into_iter()
                .map(|child| {

                    DAG_MANAGER
                        .get_or_create(child.as_ref())
                        .expect("DAG manager get_or_create failed")
                })
                .collect::<Vec<_>>();

            let node = DAG_MANAGER
                .get_or_create_normalized(
                    DagOp::$op,
                    children_nodes,
                )
                .expect("DAG manager get_or_create_normalized failed");

            Expr::Dag(node)
        }
    };
}

// --- Smart Constructors ---
impl Expr {
    // --- Unary Operator Constructors ---
    unary_constructor!(new_sin, Sin);

    unary_constructor!(new_cos, Cos);

    unary_constructor!(new_tan, Tan);

    unary_constructor!(new_exp, Exp);

    unary_constructor!(new_log, Log);

    unary_constructor!(new_neg, Neg);

    unary_constructor!(new_abs, Abs);

    unary_constructor!(new_sqrt, Sqrt);

    unary_constructor!(
        new_transpose,
        Transpose
    );

    unary_constructor!(
        new_inverse,
        Inverse
    );

    unary_constructor!(new_sec, Sec);

    unary_constructor!(new_csc, Csc);

    unary_constructor!(new_cot, Cot);

    unary_constructor!(
        new_arcsin,
        ArcSin
    );

    unary_constructor!(
        new_arccos,
        ArcCos
    );

    unary_constructor!(
        new_arctan,
        ArcTan
    );

    unary_constructor!(
        new_arcsec,
        ArcSec
    );

    unary_constructor!(
        new_arccsc,
        ArcCsc
    );

    unary_constructor!(
        new_arccot,
        ArcCot
    );

    unary_constructor!(new_sinh, Sinh);

    unary_constructor!(new_cosh, Cosh);

    unary_constructor!(new_tanh, Tanh);

    unary_constructor!(new_sech, Sech);

    unary_constructor!(new_csch, Csch);

    unary_constructor!(new_coth, Coth);

    unary_constructor!(
        new_arcsinh,
        ArcSinh
    );

    unary_constructor!(
        new_arccosh,
        ArcCosh
    );

    unary_constructor!(
        new_arctanh,
        ArcTanh
    );

    unary_constructor!(
        new_arcsech,
        ArcSech
    );

    unary_constructor!(
        new_arccsch,
        ArcCsch
    );

    unary_constructor!(
        new_arccoth,
        ArcCoth
    );

    unary_constructor!(new_not, Not);

    unary_constructor!(
        new_floor,
        Floor
    );

    unary_constructor!(
        new_gamma,
        Gamma
    );

    unary_constructor!(new_erf, Erf);

    unary_constructor!(new_erfc, Erfc);

    unary_constructor!(new_erfi, Erfi);

    unary_constructor!(new_zeta, Zeta);

    unary_constructor!(
        new_digamma,
        Digamma
    );

    // --- Binary Operator Constructors ---
    binary_constructor!(new_add, Add);

    binary_constructor!(new_sub, Sub);

    binary_constructor!(new_mul, Mul);

    binary_constructor!(new_div, Div);

    binary_constructor!(new_pow, Power);

    binary_constructor!(
        new_complex,
        Complex
    );

    binary_constructor!(
        new_matrix_mul,
        MatrixMul
    );

    binary_constructor!(
        new_matrix_vec_mul,
        MatrixVecMul
    );

    binary_constructor!(
        new_log_base,
        LogBase
    );

    binary_constructor!(
        new_atan2,
        Atan2
    );

    binary_constructor!(new_xor, Xor);

    binary_constructor!(
        new_implies,
        Implies
    );

    binary_constructor!(
        new_equivalent,
        Equivalent
    );

    binary_constructor!(new_beta, Beta);

    binary_constructor!(
        new_bessel_j,
        BesselJ
    );

    binary_constructor!(
        new_bessel_y,
        BesselY
    );

    binary_constructor!(
        new_legendre_p,
        LegendreP
    );

    binary_constructor!(
        new_laguerre_l,
        LaguerreL
    );

    binary_constructor!(
        new_hermite_h,
        HermiteH
    );

    binary_constructor!(
        new_kronecker_delta,
        KroneckerDelta
    );

    binary_constructor!(
        new_apply,
        Apply
    );

    // --- N-ary Constructors ---
    n_ary_constructor!(
        new_vector,
        Vector
    );

    n_ary_constructor!(new_and, And);

    n_ary_constructor!(new_or, Or);

    n_ary_constructor!(
        new_union,
        Union
    );

    n_ary_constructor!(
        new_polynomial,
        Polynomial
    );

    n_ary_constructor!(
        new_tuple,
        Tuple
    );

    unary_constructor_deprecated!(
        new_custom_arc_one,
        CustomArcOne
    );

    binary_constructor_deprecated!(
        new_custom_arc_two,
        CustomArcTwo
    );

    n_ary_constructor_deprecated!(
        new_custom_vec_one,
        CustomVecOne
    );

    n_ary_constructor_deprecated!(
        new_custom_vec_two,
        CustomVecTwo
    );

    n_ary_constructor_deprecated!(
        new_custom_vec_three,
        CustomVecThree
    );

    n_ary_constructor_deprecated!(
        new_custom_vec_four,
        CustomVecFour
    );

    n_ary_constructor_deprecated!(
        new_custom_vec_five,
        CustomVecFive
    );

    // --- Leaf Node Constructors ---
    /// Creates a new Constant expression, managed by the DAG.
    #[must_use]

    pub fn new_constant(
        c: f64
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Constant(
                    OrderedFloat(c),
                ),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Variable expression, managed by the DAG.
    #[must_use]

    pub fn new_variable(
        name: &str
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Variable(
                    name.to_string(),
                ),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new `BigInt` expression, managed by the DAG.
    #[must_use]

    pub fn new_bigint(
        i: BigInt
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::BigInt(i),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Rational expression, managed by the DAG.
    #[must_use]

    pub fn new_rational(
        r: BigRational
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Rational(r),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Pi expression, managed by the DAG.
    #[must_use]

    pub fn new_pi() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Pi,
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new E expression, managed by the DAG.
    #[must_use]

    pub fn new_e() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::E,
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new Infinity expression, managed by the DAG.
    #[must_use]

    pub fn new_infinity() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Infinity,
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    /// Creates a new `NegativeInfinity` expression, managed by the DAG.
    #[must_use]

    pub fn new_negative_infinity(
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::NegativeInfinity,
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    // --- Special Constructors ---
    /// Creates a new Matrix expression, managed by the DAG.

    pub fn new_matrix<I, J, T>(
        elements: I
    ) -> Self
    where
        I: IntoIterator<Item = J>,
        J: IntoIterator<Item = T>,
        T: AsRef<Self>,
    {

        let mut flat_children_nodes =
            Vec::new();

        let mut rows = 0;

        let mut cols = 0;

        for row_iter in elements {

            rows += 1;

            let mut current_cols = 0;

            for element in row_iter {

                let node = DAG_MANAGER
                    .get_or_create(element.as_ref())
                    .expect("Value is valid");

                flat_children_nodes
                    .push(node);

                current_cols += 1;
            }

            if cols == 0 {

                cols = current_cols;
            } else if current_cols
                != cols
            {

                panic!(
                    "Matrix rows must \
                     have consistent \
                     length"
                );
            }
        }

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Matrix {
                    rows,
                    cols,
                },
                flat_children_nodes,
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_predicate<I, T>(
        name: &str,
        args: I,
    ) -> Self
    where
        I: IntoIterator<Item = T>,
        T: AsRef<Self>,
    {

        let children_nodes = args
            .into_iter()
            .map(|child| {

                DAG_MANAGER
                    .get_or_create(child.as_ref())
                    .expect("Value is valid")
            })
            .collect::<Vec<_>>();

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Predicate {
                    name: name
                        .to_string(),
                },
                children_nodes,
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_forall<A>(
        var: &str,
        expr: A,
    ) -> Self
    where
        A: AsRef<Self>,
    {

        let child_node = DAG_MANAGER
            .get_or_create(
                expr.as_ref(),
            )
            .expect("Value is valid");

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::ForAll(
                    var.to_string(),
                ),
                vec![child_node],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_exists<A>(
        var: &str,
        expr: A,
    ) -> Self
    where
        A: AsRef<Self>,
    {

        let child_node = DAG_MANAGER
            .get_or_create(
                expr.as_ref(),
            )
            .expect("Value is valid");

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Exists(
                    var.to_string(),
                ),
                vec![child_node],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    pub fn new_interval<A, B>(
        lower: A,
        upper: B,
        incl_lower: bool,
        incl_upper: bool,
    ) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
    {

        let dag_lower = DAG_MANAGER
            .get_or_create(
                lower.as_ref(),
            )
            .expect("Value is valid");

        let dag_upper = DAG_MANAGER
            .get_or_create(
                upper.as_ref(),
            )
            .expect("Value is valid");

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::Interval(
                    incl_lower,
                    incl_upper,
                ),
                vec![
                    dag_lower,
                    dag_upper,
                ],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[must_use]

    pub fn new_sparse_polynomial(
        p: SparsePolynomial
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::SparsePolynomial(
                    p,
                ),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    // --- Custom Constructors ---
    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    #[must_use]

    pub fn new_custom_zero() -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::CustomZero,
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'UnaryList' variant \
                instead."
    )]
    #[must_use]

    pub fn new_custom_string(
        s: &str
    ) -> Self {

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::CustomString(
                    s.to_string(),
                ),
                vec![],
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]

    pub fn new_custom_arc_three<
        A,
        B,
        C,
    >(
        a: A,
        b: B,
        c: C,
    ) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
    {

        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");

        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");

        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");

        let children =
            vec![dag_a, dag_b, dag_c];

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::CustomArcThree,
                children,
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]

    pub fn new_custom_arc_four<
        A,
        B,
        C,
        D,
    >(
        a: A,
        b: B,
        c: C,
        d: D,
    ) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
        D: AsRef<Self>,
    {

        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");

        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");

        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");

        let dag_d = DAG_MANAGER
            .get_or_create(d.as_ref())
            .expect("Value is valid");

        let children = vec![
            dag_a, dag_b, dag_c, dag_d,
        ];

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::CustomArcFour,
                children,
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    #[deprecated(
        since = "0.1.18",
        note = "Please use the \
                'NaryList' variant \
                instead."
    )]

    pub fn new_custom_arc_five<
        A,
        B,
        C,
        D,
        E,
    >(
        a: A,
        b: B,
        c: C,
        d: D,
        e: E,
    ) -> Self
    where
        A: AsRef<Self>,
        B: AsRef<Self>,
        C: AsRef<Self>,
        D: AsRef<Self>,
        E: AsRef<Self>,
    {

        let dag_a = DAG_MANAGER
            .get_or_create(a.as_ref())
            .expect("Value is valid");

        let dag_b = DAG_MANAGER
            .get_or_create(b.as_ref())
            .expect("Value is valid");

        let dag_c = DAG_MANAGER
            .get_or_create(c.as_ref())
            .expect("Value is valid");

        let dag_d = DAG_MANAGER
            .get_or_create(d.as_ref())
            .expect("Value is valid");

        let dag_e = DAG_MANAGER
            .get_or_create(e.as_ref())
            .expect("Value is valid");

        let children = vec![
            dag_a, dag_b, dag_c, dag_d,
            dag_e,
        ];

        let node = DAG_MANAGER
            .get_or_create_normalized(
                DagOp::CustomArcFive,
                children,
            )
            .expect("Value is valid");

        Self::Dag(node)
    }

    // --- AST to DAG Migration Utilities ---

    /// Checks if this expression is in DAG form.
    ///
    /// # Returns
    /// * `true` if the expression is `Expr::Dag`, `false` otherwise
    ///
    /// # Examples
    /// ```
    /// 
    /// use rssn::symbolic::core::Expr;
    ///
    /// let dag_expr =
    ///     Expr::new_variable("x");
    ///
    /// assert!(dag_expr.is_dag());
    ///
    /// let ast_expr = Expr::Constant(1.0);
    ///
    /// assert!(!ast_expr.is_dag());
    /// ```
    #[inline]
    #[must_use]

    pub const fn is_dag(&self) -> bool {

        matches!(self, Self::Dag(_))
    }

    /// Converts this expression to DAG form if not already.
    ///
    /// This is a key function for the AST→DAG migration. It ensures that any
    /// expression, whether in old AST form or new DAG form, is converted to
    /// a canonical DAG representation.
    ///
    /// # Returns
    /// * `Ok(Expr::Dag)` - The expression in DAG form
    /// * `Err(String)` - If conversion fails
    ///
    /// # Examples
    /// ```
    /// 
    /// use std::sync::Arc;
    ///
    /// use rssn::symbolic::core::Expr;
    ///
    /// // Old AST form
    /// let ast = Expr::Add(
    ///     Arc::new(Expr::Variable(
    ///         "x".to_string(),
    ///     )),
    ///     Arc::new(Expr::Constant(1.0)),
    /// );
    ///
    /// // Convert to DAG
    /// let dag = ast
    ///     .to_dag()
    ///     .unwrap();
    ///
    /// assert!(dag.is_dag());
    /// ```

    pub fn to_dag(
        &self
    ) -> Result<Self, String> {

        match self {
            // Already in DAG form, just clone
            | Self::Dag(_) => {
                Ok(self.clone())
            },

            // Convert AST to DAG
            | _ => {

                let dag_node =
                    DAG_MANAGER
                        .get_or_create(
                            self,
                        )?;

                Ok(Self::Dag(dag_node))
            },
        }
    }

    /// Converts this expression to DAG form in-place.
    ///
    /// This is a convenience method that converts the expression to DAG form
    /// and replaces the current value.
    ///
    /// # Examples
    /// ```
    /// 
    /// use rssn::symbolic::core::Expr;
    ///
    /// let mut expr = Expr::Constant(1.0);
    ///
    /// assert!(!expr.is_dag());
    ///
    /// expr.to_dag_form();
    ///
    /// assert!(expr.is_dag());
    /// ```

    pub fn to_dag_form(&mut self) {

        if let Ok(dag) = self.to_dag() {

            *self = dag;
        }
    }

    /// Converts a DAG expression back to AST form (for serialization).
    ///
    /// This is used internally for backward-compatible serialization.
    /// Note: This may fail for expressions that cannot be represented in AST form.
    ///
    /// # Returns
    /// * `Ok(Expr)` - The expression in AST form
    /// * `Err(String)` - If conversion fails

    pub fn to_ast(
        &self
    ) -> Result<Self, String> {

        match self {
            | Self::Dag(node) => {
                node.to_expr()
            },
            | _ => Ok(self.clone()),
        }
    }
}

// --- Dynamic Operation Registry ---
//
// This registry allows for runtime registration of custom operations without
// modifying the core `Expr` enum. It's designed to support:
// 1. Plugin systems that need to add new operations
// 2. Domain-specific extensions without core changes
// 3. Future-proofing the expression system
//
// Operations registered here can be used with UnaryList, BinaryList, and NaryList
// variants, and their properties (associativity, commutativity) are used during
// normalization and simplification.

/// Properties of a dynamically registered operation.
///
/// This struct defines the characteristics of a custom operation that can be
/// registered at runtime. These properties are used by the simplification engine
/// to apply appropriate transformations.
#[derive(Debug, Clone, Default)]

pub struct DynamicOpProperties {
    /// The name of the operation (must be unique).
    pub name: String,
    /// A human-readable description of what the operation does.
    pub description: String,
    /// Whether the operation is associative: f(f(a,b),c) = f(a,f(b,c))
    pub is_associative: bool,
    /// Whether the operation is commutative: f(a,b) = f(b,a)
    pub is_commutative: bool,
}

lazy_static! {
    /// Global registry for dynamically registered operations.
    ///
    /// This registry maps operation names to their properties. It's thread-safe
    /// and can be accessed from multiple threads simultaneously.
    pub static ref DYNAMIC_OP_REGISTRY: RwLock<HashMap<String, DynamicOpProperties>> = RwLock::new(HashMap::new());
}

/// Registers a dynamic operation with the global registry.
///
/// This function allows you to add custom operations at runtime without modifying
/// the core `Expr` enum. Once registered, operations can be used with `UnaryList`,
/// `BinaryList`, or `NaryList` variants.
///
/// # Arguments
/// * `name` - The unique name of the operation
/// * `props` - The properties of the operation (associativity, commutativity, etc.)
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::register_dynamic_op;
/// use rssn::symbolic::core::DynamicOpProperties;
///
/// register_dynamic_op(
///     "my_custom_op",
///     DynamicOpProperties {
///         name : "my_custom_op".to_string(),
///         description : "A custom commutative operation"
///             .to_string(),
///         is_associative : true,
///         is_commutative : true,
///     },
/// );
/// ```

pub fn register_dynamic_op(
    name: &str,
    props: DynamicOpProperties,
) {

    let mut registry =
        DYNAMIC_OP_REGISTRY
            .write()
            .unwrap();

    registry.insert(
        name.to_string(),
        props,
    );
}

/// Retrieves the properties of a dynamically registered operation.
///
/// This function looks up an operation in the global registry and returns its
/// properties if found.
///
/// # Arguments
/// * `name` - The name of the operation to look up
///
/// # Returns
/// * `Some(DynamicOpProperties)` if the operation is registered
/// * `None` if the operation is not found
///
/// # Examples
/// ```
/// 
/// use rssn::symbolic::core::get_dynamic_op_properties;
/// use rssn::symbolic::core::register_dynamic_op;
/// use rssn::symbolic::core::DynamicOpProperties;
///
/// register_dynamic_op(
///     "test_op",
///     DynamicOpProperties {
///         name : "test_op".to_string(),
///         description : "Test operation".to_string(),
///         is_associative : false,
///         is_commutative : true,
///     },
/// );
///
/// let props = get_dynamic_op_properties("test_op");
///
/// assert!(props.is_some());
///
/// assert!(
///     props
///         .unwrap()
///         .is_commutative
/// );
/// ```
#[must_use]

pub fn get_dynamic_op_properties(
    name: &str
) -> Option<DynamicOpProperties> {

    let registry = DYNAMIC_OP_REGISTRY
        .read()
        .unwrap();

    registry
        .get(name)
        .cloned()
}
