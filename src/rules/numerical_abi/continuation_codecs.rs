// Typed, dimensionless adapters for complex Taylor data and sequence sums.
use crate::kernels::analytic_continuation as continuation;

struct ComplexPade {
    numerator: Vec<Complex64>,
    denominator: Vec<Complex64>,
}
record!(ComplexPade {
    numerator: Vec<Complex64>,
    denominator: Vec<Complex64>
});

impl Arg for ComplexPade {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        use super::native_object_api::field;

        let Value::Record(record) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        if record.fields.len() != 2 {
            return Err(OpError::Domain(
                "expected a Padé numerator and denominator".into(),
            ));
        }
        Ok(Self {
            numerator: Vec::<Complex64>::decode(field(value, "numerator")?, index)?,
            denominator: Vec::<Complex64>::decode(field(value, "denominator")?, index)?,
        })
    }
}

binding!(
    continuation_pade,
    "num_continuation_pade",
    "numerical_special",
    "Forms a complex [m/n] Padé approximant from Taylor coefficients in ascending order. Supply at least m+n+1 dimensionless coefficients. The denominator has constant coefficient one and can drop in degree for a rank-deficient system. The denominator degree may not exceed 64.",
    "crate::kernels::analytic_continuation::pade_approximant_complex",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<super::Complex64>, m: usize, n: usize) -> Result<super::ComplexPade, OpError>,
    |coefficients: Vec<super::Complex64>, m: usize, n: usize| {
        if n > 64 {
            return Err(OpError::Domain("Padé denominator degree exceeds 64".into()));
        }
        let (numerator, denominator) = super::continuation::pade_approximant_complex(&coefficients, m, n)
            .ok_or_else(|| OpError::Domain("insufficient Taylor coefficients for the requested Padé order".into()))?;
        Ok(super::ComplexPade { numerator, denominator })
    }
);

binding!(
    continuation_eval_pade,
    "num_continuation_eval_pade",
    "numerical_special",
    "Evaluates the dimensionless Padé numerator divided by its denominator at a complex point. Pass the record returned by num_continuation_pade. A zero denominator or unrepresentable result returns a domain error.",
    "crate::kernels::analytic_continuation::eval_pade_complex",
    DimensionPolicy::Dimensionless,
    (pade: super::ComplexPade, point: super::Complex64) -> super::Complex64,
    |pade: super::ComplexPade, point: super::Complex64| {
        super::continuation::eval_pade_complex(&pade.numerator, &pade.denominator, point)
    }
);

binding!(
    continuation_wynn,
    "num_continuation_wynn",
    "numerical_special",
    "Accelerates convergence of up to 4096 dimensionless complex partial sums with Wynn's epsilon algorithm. With no finite accelerated estimate, returns the last sum. An empty sequence returns zero.",
    "crate::kernels::analytic_continuation::wynn_epsilon_complex",
    DimensionPolicy::Dimensionless,
    (partial_sums: Vec<super::Complex64>) -> Result<super::Complex64, OpError>,
    |partial_sums: Vec<super::Complex64>| {
        if partial_sums.len() > 4096 {
            return Err(OpError::Domain("Wynn's epsilon algorithm exceeds the sequence work budget".into()));
        }
        Ok(super::continuation::wynn_epsilon_complex(&partial_sums))
    }
);

binding!(
    continuation_recenter,
    "num_continuation_recenter",
    "numerical_special",
    "Re-centers ascending complex Taylor coefficients along the supplied chain of dimensionless centers. The first center corresponds to the input coefficients and the last to the output. Empty inputs and work above eight million coefficient pairs are rejected.",
    "crate::kernels::analytic_continuation::weierstrass_chain_continue",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<super::Complex64>, centers: Vec<super::Complex64>) -> Result<Vec<super::Complex64>, OpError>,
    |coefficients: Vec<super::Complex64>, centers: Vec<super::Complex64>| {
        let n = coefficients.len();
        let steps = centers.len().saturating_sub(1);
        let pairs = n.saturating_mul(n.saturating_add(1)) / 2;
        if pairs.saturating_mul(steps) > 8_000_000 {
            return Err(OpError::Domain("Taylor recentering exceeds the coefficient work budget".into()));
        }
        super::continuation::weierstrass_chain_continue(&coefficients, &centers)
            .ok_or_else(|| OpError::Domain("Taylor recentering needs coefficients and at least one center".into()))
    }
);

struct HolonomicMatrix(Vec<Vec<Complex64>>);

impl Arg for HolonomicMatrix {
    fn ty() -> TypeExpr {
        <Self as Output>::ty()
    }

    fn decode(
        value: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(matrix) = value else {
            return Err(OpError::Type {
                expected: Box::new(<Self as Arg>::ty()),
                got: Box::new(value.type_expr()),
                arg: index,
            });
        };
        let [rows, width] = matrix.shape.as_slice() else {
            return Err(OpError::Domain(
                "holonomic coefficients must be a matrix".into(),
            ));
        };
        if *rows == 0 || *width == 0 || rows.checked_mul(*width) != Some(matrix.data.len()) {
            return Err(OpError::Domain(
                "holonomic coefficients need finite rectangular rows".into(),
            ));
        }
        if *rows > 17 || *width > 33 {
            return Err(OpError::BudgetExceeded);
        }
        let coefficients = matrix
            .data
            .chunks_exact(*width)
            .map(|row| {
                row.iter()
                    .map(|coefficient| Complex64::decode(coefficient, index))
                    .collect()
            })
            .collect::<Result<Vec<Vec<_>>, OpError>>()?;
        if coefficients
            .iter()
            .flatten()
            .any(|coefficient| !coefficient.is_finite())
        {
            return Err(OpError::Domain(
                "holonomic coefficients must be finite".into(),
            ));
        }
        Ok(Self(coefficients))
    }
}

impl Output for HolonomicMatrix {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: crate::kernel::signature::NatExpr::Lit(2),
            shape: None,
            elem: Box::new(TypeExpr::Complex),
        }
    }

    fn encode(self) -> Result<Value, OpError> {
        let rows = self.0.len();
        let columns = self.0.first().map_or(0, Vec::len);
        let data = self
            .0
            .into_iter()
            .flatten()
            .map(Output::encode)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Value::Tensor(crate::kernel::value::TensorValue::of(
            vec![rows, columns],
            data,
            TypeExpr::Complex,
        )))
    }
}

fn validate_holonomic_path(
    path: &[Complex64],
    order: usize,
) -> Result<(), OpError> {
    if path.len() > 1024 || order > 128 {
        return Err(OpError::BudgetExceeded);
    }
    if path.is_empty() || order == 0 || path.iter().any(|point| !point.is_finite()) {
        return Err(OpError::Domain(
            "holonomic continuation requires finite path points and a positive Taylor order".into(),
        ));
    }
    Ok(crate::cancellation::checkpoint()?)
}

binding!(
    continuation_holonomic_order,
    "num_continuation_holonomic_order",
    "numerical_special",
    "Returns the order of an ODE described by rows p_0 through p_r with polynomial coefficients in ascending powers, padded with zeros.",
    "crate::kernels::analytic_continuation::HolonomicOde::order",
    DimensionPolicy::Dimensionless,
    (coefficients: super::HolonomicMatrix) -> usize,
    |coefficients: super::HolonomicMatrix| super::continuation::HolonomicOde { coeffs: coefficients.0 }.order()
);

binding!(
    continuation_holonomic_follow,
    "num_continuation_holonomic_follow",
    "numerical_special",
    "Continues the derivative state [y, y', ...] of a dimensionless polynomial-coefficient ODE along a complex polygonal path. Returns the endpoint state, or an error for a singular or unresolved step.",
    "crate::kernels::analytic_continuation::continue_holonomic_along_path",
    DimensionPolicy::Dimensionless,
    (coefficients: super::HolonomicMatrix, path: Vec<super::Complex64>, initial: Vec<super::Complex64>, taylor_order: usize) -> Result<Vec<super::Complex64>, OpError>,
    |coefficients: super::HolonomicMatrix, path: Vec<super::Complex64>, initial: Vec<super::Complex64>, taylor_order: usize| {
        super::validate_holonomic_path(&path,taylor_order)?;
        let ode = super::continuation::HolonomicOde { coeffs: coefficients.0 };
        let state = super::continuation::continue_holonomic_along_path(&ode, &path, &initial, taylor_order);
        crate::cancellation::checkpoint()?;
        state.ok_or_else(|| OpError::NotReduced("the holonomic path reaches a singular point or lacks a converged Taylor step".into()))
    }
);

binding!(
    continuation_holonomic_monodromy,
    "num_continuation_holonomic_monodromy",
    "numerical_special",
    "Continues the derivative-state identity basis around the closed complex loop; an omitted closing point is appended. Matrix columns correspond to the initial basis vectors.",
    "crate::kernels::analytic_continuation::holonomic_monodromy_matrix",
    DimensionPolicy::Dimensionless,
    (coefficients: super::HolonomicMatrix, path: Vec<super::Complex64>, taylor_order: usize) -> Result<super::HolonomicMatrix, OpError>,
    |coefficients: super::HolonomicMatrix, path: Vec<super::Complex64>, taylor_order: usize| {
        super::validate_holonomic_path(&path,taylor_order)?;
        let ode = super::continuation::HolonomicOde { coeffs: coefficients.0 };
        let monodromy = super::continuation::holonomic_monodromy_matrix(&ode, &path, taylor_order);
        crate::cancellation::checkpoint()?;
        monodromy.map(super::HolonomicMatrix).ok_or_else(|| OpError::NotReduced("monodromy reached a singular point or exceeded the Taylor refinement limit".into()))
    }
);

binding!(
    continuation_borel_pade,
    "num_continuation_borel_pade",
    "numerical_special",
    "Borel–Padé resummation of an ascending dimensionless asymptotic series through order m+n. Poles on the positive Borel axis select the upper lateral integral.",
    "crate::kernels::analytic_continuation::borel_pade_resum",
    DimensionPolicy::Dimensionless,
    (coefficients: Vec<super::Complex64>, point: super::Complex64, m: usize, n: usize) -> Result<super::Complex64, OpError>,
    |coefficients: Vec<super::Complex64>, point: super::Complex64, m: usize, n: usize| {
        let requested = m.checked_add(n).ok_or(OpError::BudgetExceeded)?;
        if n > 64 || requested > 128 {
            return Err(OpError::BudgetExceeded);
        }
        if coefficients.len() <= requested || !point.is_finite() || coefficients.iter().any(|coefficient| !coefficient.is_finite()) {
            return Err(OpError::Domain("Borel–Padé needs a finite point and m+n+1 finite coefficients".into()));
        }
        crate::cancellation::checkpoint()?;
        let value = super::continuation::borel_pade_resum(&coefficients,point,m,n);
        crate::cancellation::checkpoint()?;
        value.ok_or_else(|| OpError::NotReduced("the rational Borel integral did not converge".into()))
    }
);
