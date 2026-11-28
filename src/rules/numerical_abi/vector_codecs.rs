// Vector adapters which reuse decoded storage or check output size before calling the kernel.
binding!(
    vector_axpy,
    "num_vector_axpy",
    "numerical_vector",
    "`num_vector_axpy(alpha, x, y)` computes alpha*x + y with the native SIMD kernel. The scaled increment and the retained input y must have the same dimension. Returns a new vector.",
    "crate::kernels::vector::axpy",
    DimensionPolicy::Axpy,
    (alpha: f64, x: Vec<f64>, y: Vec<f64>) -> Result<Vec<f64>, OpError>,
    |alpha: f64, x: Vec<f64>, mut y: Vec<f64>| {
        crate::kernels::vector::axpy(alpha, &x, &mut y).map_err(OpError::Domain)?;
        Ok(y)
    }
);
binding!(
    vector_axpby,
    "num_vector_axpby",
    "numerical_vector",
    "`num_vector_axpby(alpha, x, beta, y)` computes alpha*x + beta*y with the native SIMD kernel. Both scaled terms must have the same dimension. Returns a new vector in those units.",
    "crate::kernels::vector::axpby",
    DimensionPolicy::Axpby,
    (alpha: f64, x: Vec<f64>, beta: f64, y: Vec<f64>) -> Result<Vec<f64>, OpError>,
    |alpha: f64, x: Vec<f64>, beta: f64, mut y: Vec<f64>| {
        crate::kernels::vector::axpby(alpha, &x, beta, &mut y).map_err(OpError::Domain)?;
        Ok(y)
    }
);
binding!(
    vector_outer_product,
    "num_vector_outer_product",
    "numerical_vector",
    "`num_vector_outer_product(left, right)` returns a row-major matrix with shape [left length, right length] and product units. The result may contain at most 8388608 entries.",
    "crate::kernels::vector::outer_product",
    DimensionPolicy::Product,
    (left: Vec<f64>, right: Vec<f64>) -> Result<dense::Mat, OpError>,
    |left: Vec<f64>, right: Vec<f64>| {
        if left.len().saturating_mul(right.len()) > 8_388_608 {
            return Err(OpError::Domain("outer product exceeds the result allocation limit".into()));
        }
        Ok(crate::kernels::vector::outer_product(&left, &right))
    }
);
