use std::marker::PhantomData;

use super::cost::Extractor;
use super::egraph::EGraph;
use super::heuristics::HeuristicBudget;
use super::rules::Rule;
use super::rules::algebra::{ConstantFoldingRule, IdentityReductionRule};
use super::rules::calculus::DifferentiationRule;
use super::rules::oracles::{
    CombinatoricsOracleRule, ElementaryOracleRule, EulerLagrangeOracleRule, FactorOracleRule,
    FunctionalAnalysisOracleRule, GrobnerOracleRule, IndefiniteSumOracleRule,
    IntegralEquationsOracleRule, IntegralOracleRule, LieAlgebraOracleRule, LimitOracleRule,
    LogicOracleRule, MatrixOracleRule, NumberTheoryOracleRule, OdeOracleRule,
    OptimizationOracleRule, PdeOracleRule, PolynomialOracleRule, ProductOracleRule,
    QuantumOracleRule, RadicalsOracleRule, ResidueOracleRule, SeriesOracleRule, SolveOracleRule,
    SpecialFunctionsOracleRule, StatsOracleRule, SumOracleRule, TransformOracleRule,
    VectorCalculusOracleRule,
};
use super::rules::trig::TrigIdentitiesRule;
use crate::symbolic::core::Expr;

// --- Typestate Markers ---

/// Initial typestate: No rules or configuration configured.
pub struct Unconfigured;

/// Typestate indicating that rules have been registered.
pub struct HasRules;

/// Final typestate: Pipeline is compiled, validated, and ready for execution.
pub struct Ready;

/// A chained typestate builder for configuring and compiling E-Graph saturation pipelines.
pub struct PipelineBuilder<State> {
    budget: HeuristicBudget,
    rules: Vec<Box<dyn Rule>>,
    _state: PhantomData<State>,
}

impl Default for PipelineBuilder<Unconfigured> {
    fn default() -> Self {
        Self::new()
    }
}

impl PipelineBuilder<Unconfigured> {
    /// Starts constructing a new E-Graph pipeline with the typestate machine.
    #[must_use]
    pub fn new() -> Self {
        Self {
            budget: HeuristicBudget::default(),
            rules: Vec::new(),
            _state: PhantomData,
        }
    }
}

impl<S> PipelineBuilder<S> {
    /// Enables algebraic simplification and constant folding rules.
    pub fn with_algebra(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(ConstantFoldingRule));
        self.rules.push(Box::new(IdentityReductionRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Enables calculus rules (differentiation chain rules, product rules, etc.).
    pub fn with_calculus(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(DifferentiationRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Enables trigonometric and hyperbolic identities.
    pub fn with_trigonometry(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(TrigIdentitiesRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Ordinary Differential Equation solver oracle (bridges `ode.rs`).
    pub fn with_ode_solver(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(OdeOracleRule::new()));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the algebraic equation solver oracle (bridges `solve.rs`).
    pub fn with_equation_solver(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(SolveOracleRule::new()));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Gröbner basis relations simplifier oracle (bridges `cas_foundations.rs`).
    pub fn with_grobner_simplifier(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(GrobnerOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Symbolic Integration solver oracle (bridges `calculus::integrate_internal`).
    pub fn with_integral_solver(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(IntegralOracleRule::default()));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Limit solver oracle (bridges `calculus::limit_internal`).
    pub fn with_limit_solver(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(LimitOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the PDE solver oracle (bridges `pde::solve_pde`).
    pub fn with_pde_solver(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(PdeOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Polynomial Factorization oracle (bridges `cas_foundations::factorize`).
    pub fn with_factorizer(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(FactorOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Series expansion oracle (bridges `series::taylor_series`).
    pub fn with_series_expansion(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(SeriesOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Summation oracle (bridges `series::summation`).
    pub fn with_summation(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(SumOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Product oracle (bridges `series::product`).
    pub fn with_product(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(ProductOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Indefinite Summation oracle (bridges closed-form antidifference solver).
    pub fn with_indefinite_sum(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(IndefiniteSumOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Radicals Denesting oracle (bridges `radicals::denest_sqrt`).
    pub fn with_radicals_denesting(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(RadicalsOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Matrix Inversion oracle (bridges `matrix::inverse_matrix`).
    pub fn with_matrix_inversion(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(MatrixOracleRule::default()));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Transforms oracle (bridges Laplace, Fourier, Z-transforms).
    pub fn with_transforms(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(TransformOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Vector Calculus oracle (bridges Gradient, Divergence, Curl, Laplacian).
    pub fn with_vector_calculus(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(VectorCalculusOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Euler-Lagrange variational oracle.
    pub fn with_euler_lagrange(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(EulerLagrangeOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Quantum operator commutator oracle.
    pub fn with_quantum(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(QuantumOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Complex Analysis residue oracle.
    pub fn with_residue(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(ResidueOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Combinatorics and Number Theory oracle.
    pub fn with_combinatorics(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(CombinatoricsOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Propositional Logic oracle.
    pub fn with_logic(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(LogicOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Number Theory oracle (bridges Mod, Gcd, Lcm, IsPrime).
    pub fn with_number_theory(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(NumberTheoryOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Elementary Algebra oracle (bridges expand, pow_expand, etc.).
    pub fn with_elementary(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(ElementaryOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Special Functions oracle (bridges Gamma, Zeta, Erf, Bessel).
    pub fn with_special_functions(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(SpecialFunctionsOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Polynomial Algebra oracle (bridges polynomial GCD and division).
    pub fn with_polynomial_ops(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(PolynomialOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Optimization oracle (bridges Hessian matrix).
    pub fn with_optimization(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(OptimizationOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Functional Analysis oracle (bridges Hilbert and Banach spaces).
    pub fn with_functional_analysis(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(FunctionalAnalysisOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Lie Algebra oracle (bridges Lie bracket commutator).
    pub fn with_lie_algebra(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(LieAlgebraOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Integral Equations oracle (bridges airfoil equation solver).
    pub fn with_integral_equations(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(IntegralEquationsOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects the Statistics and Information Theory oracle (bridges mean, variance, entropy, etc.).
    pub fn with_stats(mut self) -> PipelineBuilder<HasRules> {
        self.rules.push(Box::new(StatsOracleRule));
        PipelineBuilder {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }

    /// Injects all standard advanced algorithm oracles.
    pub fn with_all_oracles(self) -> PipelineBuilder<HasRules> {
        self.with_ode_solver()
            .with_equation_solver()
            .with_grobner_simplifier()
            .with_integral_solver()
            .with_limit_solver()
            .with_pde_solver()
            .with_factorizer()
            .with_series_expansion()
            .with_summation()
            .with_product()
            .with_indefinite_sum()
            .with_radicals_denesting()
            .with_matrix_inversion()
            .with_transforms()
            .with_vector_calculus()
            .with_euler_lagrange()
            .with_quantum()
            .with_residue()
            .with_combinatorics()
            .with_logic()
            .with_number_theory()
            .with_elementary()
            .with_special_functions()
            .with_polynomial_ops()
            .with_optimization()
            .with_functional_analysis()
            .with_lie_algebra()
            .with_integral_equations()
            .with_stats()
    }

    /// Sets the heuristic saturation budget and limits.
    #[must_use]
    pub fn with_heuristic_budget(mut self, budget: HeuristicBudget) -> Self {
        self.budget = budget;
        self
    }
}

impl PipelineBuilder<HasRules> {
    /// Compiles the configured rules and budget into an executable `EGraphPipeline<Ready>`.
    /// The typestate transition guarantees at compile-time that rules are registered.
    #[must_use]
    pub fn build(self) -> EGraphPipeline<Ready> {
        EGraphPipeline {
            budget: self.budget,
            rules: self.rules,
            _state: PhantomData,
        }
    }
}

/// The executable E-Graph pipeline in the `Ready` typestate.
pub struct EGraphPipeline<State> {
    budget: HeuristicBudget,
    rules: Vec<Box<dyn Rule>>,
    _state: PhantomData<State>,
}

impl EGraphPipeline<Ready> {
    /// Creates a builder to configure a custom pipeline.
    #[must_use]
    pub fn builder() -> PipelineBuilder<Unconfigured> {
        PipelineBuilder::new()
    }

    /// Constructs a standard, fully-equipped pipeline with all core rules and oracles enabled.
    #[must_use]
    pub fn standard() -> Self {
        PipelineBuilder::new()
            .with_algebra()
            .with_calculus()
            .with_trigonometry()
            .with_all_oracles()
            .build()
    }

    /// Simplifies and evaluates an `Expr` using heuristic equality saturation,
    /// returning the minimal-cost canonical `Expr` with DAG structural sharing preserved.
    pub fn simplify(&self, expr: &Expr) -> Expr {
        let mut egraph = EGraph::new();
        let root = egraph.add_expr(expr);
        egraph.rebuild();

        // Separate rules by priority tiers
        let tier0_rules: Vec<&dyn Rule> = self.rules.iter().filter(|r| r.tier() == 0).map(AsRef::as_ref).collect();
        let tier1_rules: Vec<&dyn Rule> = self.rules.iter().filter(|r| r.tier() == 1).map(AsRef::as_ref).collect();
        let other_rules: Vec<&dyn Rule> = self.rules.iter().filter(|r| r.tier() >= 2).map(AsRef::as_ref).collect();

        for _iter in 0..self.budget.max_iterations {
            let mut iteration_changes = 0;

            // Phase 1: High-Weight De-cocooning (Derivatives, Integrals, ODEs, Solves)
            if self.budget.prioritize_high_weight_ops {
                for rule in &tier0_rules {
                    iteration_changes += rule.apply(&mut egraph);
                }
                egraph.rebuild();
            }

            // Phase 2: Constant folding & Identity reduction
            for rule in &tier1_rules {
                iteration_changes += rule.apply(&mut egraph);
            }
            egraph.rebuild();

            // Phase 3: Structural transformations (if within budget)
            if egraph.total_nodes() < self.budget.max_nodes && egraph.total_classes() < self.budget.max_classes {
                for rule in &other_rules {
                    iteration_changes += rule.apply(&mut egraph);
                }
                egraph.rebuild();
            }

            // Stop if reached a fixpoint (no changes occurred in this iteration)
            if iteration_changes == 0 {
                break;
            }

            // Stop if exceeded memory/node bounds
            if egraph.total_nodes() >= self.budget.max_nodes || egraph.total_classes() >= self.budget.max_classes {
                break;
            }
        }

        let extractor = Extractor::new(&egraph);
        extractor.extract(&egraph, root)
    }
}
