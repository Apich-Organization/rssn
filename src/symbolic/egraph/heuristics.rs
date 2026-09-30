/// Budget and heuristic scheduling configuration for the E-Graph engine.
#[derive(Clone, Debug)]
pub struct HeuristicBudget {
    /// Maximum search depth for chained rule applications.
    pub max_depth: usize,
    /// Maximum number of saturation iterations before forcing extraction.
    pub max_iterations: usize,
    /// Upper limit on the total number of nodes in the E-Graph to prevent memory bloat.
    pub max_nodes: usize,
    /// Upper limit on the total number of classes in the E-Graph.
    pub max_classes: usize,
    /// Whether to prioritize de-cocooning high-weight operators (Derivative, Integral, Ode, etc.)
    /// before permitting explosive general rewrites (associativity, distributivity).
    pub prioritize_high_weight_ops: bool,
}

impl Default for HeuristicBudget {
    fn default() -> Self {
        Self {
            max_depth: 10,
            max_iterations: 25,
            max_nodes: 50_000,
            max_classes: 20_000,
            prioritize_high_weight_ops: true,
        }
    }
}

impl HeuristicBudget {
    /// Creates a new default heuristic budget.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the maximum iteration depth.
    #[must_use]
    pub fn max_depth(mut self, depth: usize) -> Self {
        self.max_depth = depth;
        self
    }

    /// Sets the maximum number of iterations.
    #[must_use]
    pub fn max_iterations(mut self, iters: usize) -> Self {
        self.max_iterations = iters;
        self
    }

    /// Sets the node budget.
    #[must_use]
    pub fn max_nodes(mut self, nodes: usize) -> Self {
        self.max_nodes = nodes;
        self
    }

    /// Sets the class budget.
    #[must_use]
    pub fn max_classes(mut self, classes: usize) -> Self {
        self.max_classes = classes;
        self
    }

    /// Configures high-weight de-cocooning priority.
    #[must_use]
    pub fn prioritize_high_weight_ops(mut self, enable: bool) -> Self {
        self.prioritize_high_weight_ops = enable;
        self
    }
}
