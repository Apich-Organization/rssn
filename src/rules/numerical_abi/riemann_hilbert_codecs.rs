// A fixed row-major tensor is the portable complex matrix representation.
use crate::kernels::riemann_hilbert::Matrix2 as RhMatrix;
use crate::kernels::riemann_hilbert::MatrixRhSolution as RhSolution;

impl Arg for RhMatrix {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(2),
            shape: Some(vec![NatExpr::Lit(2), NatExpr::Lit(2)]),
            elem: Box::new(TypeExpr::Complex),
        }
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(tensor) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if tensor.shape != [2, 2] {
            return Err(OpError::Domain(
                "Riemann–Hilbert matrices must have shape [2, 2]".into(),
            ));
        }
        let cells = tensor.data.as_slice().try_into().map_err(|_| {
            OpError::Domain("Riemann–Hilbert matrix must contain four entries".into())
        })?;
        Ok(Self {
            m: decode_fixed_array(cells, index)?,
        })
    }
}

impl Output for RhMatrix {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        Ok(Value::Tensor(TensorValue::of(
            vec![2, 2],
            self.m
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            TypeExpr::Complex,
        )))
    }
}

const RH_MAX_NODES: usize = 256;

impl Arg for RhSolution {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("nodes".into(), <Vec<f64> as Arg>::ty()),
            ("weights".into(), <Vec<f64> as Arg>::ty()),
            (
                "density".into(),
                TypeExpr::Tensor {
                    rank: NatExpr::Lit(3),
                    shape: None,
                    elem: Box::new(TypeExpr::Complex),
                },
            ),
        ])
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Record(record) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        let field = |name| {
            record
                .get(name)
                .ok_or_else(|| OpError::Domain(format!("missing Riemann–Hilbert field {name}")))
        };
        let nodes = field("nodes")?;
        let weights = field("weights")?;
        let Value::Tensor(density) = field("density")? else {
            return Err(OpError::Domain("jump density must be a tensor".into()));
        };
        let [n, 2, 2] = density.shape.as_slice() else {
            return Err(OpError::Domain(
                "jump density must have shape [nodes, 2, 2]".into(),
            ));
        };
        if *n > RH_MAX_NODES {
            return Err(OpError::BudgetExceeded);
        }
        if *n == 0 || density.data.len() != 4 * n {
            return Err(OpError::Domain(
                "invalid Riemann–Hilbert density size".into(),
            ));
        }
        for vector in [nodes, weights] {
            if !matches!(vector, Value::Tensor(tensor) if tensor.shape == [*n] && tensor.data.len() == *n)
            {
                return Err(OpError::Domain(
                    "nodes, weights and density lengths must agree".into(),
                ));
            }
        }
        let nodes = Vec::<f64>::decode(nodes, index)?;
        let weights = Vec::<f64>::decode(weights, index)?;
        if nodes.iter().any(|&x| !(-1.0 < x && x < 1.0))
            || nodes.windows(2).any(|pair| pair[0] >= pair[1])
            || weights.iter().any(|&w| !w.is_finite() || w <= 0.0)
        {
            return Err(OpError::Domain(
                "expected increasing nodes inside (-1, 1) and positive finite quadrature weights"
                    .into(),
            ));
        }
        let (matrices, _) = density.data.as_chunks::<4>();
        let density = matrices
            .iter()
            .map(|entries| {
                Ok(RhMatrix {
                    m: decode_fixed_array(entries, index)?,
                })
            })
            .collect::<Result<_, OpError>>()?;
        Ok(Self {
            nodes,
            weights,
            density,
        })
    }
}

impl Output for RhSolution {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let density = Value::Tensor(TensorValue::of(
            vec![self.density.len(), 2, 2],
            self.density
                .into_iter()
                .flat_map(|matrix| matrix.m)
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            TypeExpr::Complex,
        ));
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("nodes".into(), self.nodes.encode()?),
                ("weights".into(), self.weights.encode()?),
                ("density".into(), density),
            ],
        }))
    }
}

fn rh_solution_value(
    solution: &RhSolution,
    point: num_complex::Complex64,
) -> Result<RhMatrix, OpError> {
    if point.im == 0.0 && (-1.0..=1.0).contains(&point.re) {
        return Err(OpError::Domain(
            "Riemann–Hilbert evaluation requires a point off the closed arc [-1, 1]".into(),
        ));
    }
    crate::cancellation::checkpoint()?;
    let result = solution.eval(point);
    crate::cancellation::checkpoint()?;
    Ok(result)
}

binding!(rh_matrix_new, "num_riemann_hilbert_matrix_new", "numerical_matrix", "Constructs a complex 2x2 matrix from four row-major entries of the same dimension.", "crate::kernels::riemann_hilbert::Matrix2::new", DimensionPolicy::Preserve,
    (entries: [num_complex::Complex64;4]) -> super::RhMatrix,
    |[a,b,c,d]: [num_complex::Complex64;4]| super::RhMatrix::new(a,b,c,d));
binding!(rh_matrix_identity, "num_riemann_hilbert_matrix_identity", "numerical_matrix", "Returns the dimensionless 2x2 identity matrix.", "crate::kernels::riemann_hilbert::Matrix2::identity", DimensionPolicy::Dimensionless,
    () -> super::RhMatrix, super::RhMatrix::identity);
binding!(rh_matrix_zero, "num_riemann_hilbert_matrix_zero", "numerical_matrix", "Returns the dimensionless 2x2 zero matrix.", "crate::kernels::riemann_hilbert::Matrix2::zero", DimensionPolicy::Dimensionless,
    () -> super::RhMatrix, super::RhMatrix::zero);
binding!(rh_matrix_add, "num_riemann_hilbert_matrix_add", "numerical_matrix", "Adds two 2x2 matrices with the same entry dimension.", "crate::kernels::riemann_hilbert::Matrix2::add", DimensionPolicy::EqualPreserve,
    (first: super::RhMatrix, second: super::RhMatrix) -> super::RhMatrix, super::RhMatrix::add);
binding!(rh_matrix_sub, "num_riemann_hilbert_matrix_sub", "numerical_matrix", "Subtracts two 2x2 matrices with the same entry dimension.", "crate::kernels::riemann_hilbert::Matrix2::sub", DimensionPolicy::EqualPreserve,
    (first: super::RhMatrix, second: super::RhMatrix) -> super::RhMatrix, super::RhMatrix::sub);
binding!(rh_matrix_mul, "num_riemann_hilbert_matrix_mul", "numerical_matrix", "Multiplies two 2x2 matrices. Result entries have the product dimension.", "crate::kernels::riemann_hilbert::Matrix2::mul", DimensionPolicy::Product,
    (first: super::RhMatrix, second: super::RhMatrix) -> super::RhMatrix, super::RhMatrix::mul);
binding!(rh_matrix_scale, "num_riemann_hilbert_matrix_scale", "numerical_matrix", "Scales a 2x2 matrix by a complex scalar. Result entries have the product dimension.", "crate::kernels::riemann_hilbert::Matrix2::scale", DimensionPolicy::Product,
    (matrix: super::RhMatrix, scalar: num_complex::Complex64) -> super::RhMatrix, super::RhMatrix::scale);
binding!(rh_matrix_inv, "num_riemann_hilbert_matrix_inv", "numerical_matrix", "Returns the inverse in reciprocal units. None denotes a numerically singular matrix or an inverse outside binary64 range.", "crate::kernels::riemann_hilbert::Matrix2::inv", DimensionPolicy::Inverse,
    (matrix: super::RhMatrix) -> Option<super::RhMatrix>, super::RhMatrix::inv);
binding!(rh_matrix_norm, "num_riemann_hilbert_matrix_norm", "numerical_matrix", "Returns a scale-stable Frobenius norm retaining the matrix entry dimension.", "crate::kernels::riemann_hilbert::Matrix2::norm", DimensionPolicy::Preserve,
    (matrix: super::RhMatrix) -> f64, super::RhMatrix::norm);
binding!(rh_solution_eval, "num_riemann_hilbert_solution_eval", "numerical_matrix", "Evaluates a portable Riemann–Hilbert collocation record at a dimensionless complex point off [-1, 1]. Nodes and positive weights have shape [n], and jump density has shape [n, 2, 2].", "crate::kernels::riemann_hilbert::MatrixRhSolution::eval", DimensionPolicy::Dimensionless,
    (solution: super::RhSolution, point: num_complex::Complex64) -> Result<super::RhMatrix, OpError>,
    |solution: super::RhSolution, point: num_complex::Complex64| super::rh_solution_value(&solution, point));
