type NativeDenseFactors<T> = (Vec<T>, Vec<T>, Vec<T>);
// Explicit normalized enum/value codecs for native linear APIs.
impl Arg for crate::kernels::dense::MatrixNorm {
    fn ty() -> TypeExpr {
        TypeExpr::String
    }
    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        match <String as Arg>::decode(v, i)?.as_str() {
            | "One" => Ok(Self::One),
            | "Inf" => Ok(Self::Inf),
            | "Frobenius" => Ok(Self::Frobenius),
            | "Spectral" => Ok(Self::Spectral),
            | _ => Err(OpError::Domain(
                "matrix norm must be One, Inf, Frobenius or Spectral".into(),
            )),
        }
    }
}
record!(crate::kernels::dense::Eigen{values:Vec<(f64,f64)>,vectors:Vec<Vec<num_complex::Complex64>>,residuals:Vec<f64>,defective:bool});
record!(crate::kernels::dense::GeneralizedEigen{values:Vec<(f64,f64)>,vectors:Option<dense::Mat>,b_positive_definite:bool});
impl Output for Vec<dense::Mat> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: None,
            elem: Box::new(<dense::Mat as Output>::ty()),
        }
    }
    fn encode(self) -> Result<Value, OpError> {
        if self.len() > 64
            || self
                .iter()
                .try_fold(0usize, |n, m| n.checked_add(m.data.len()))
                .is_none_or(|n| n > 262144)
        {
            return Err(OpError::Domain(
                "matrix family allocation budget exceeded".into(),
            ));
        }
        Ok(Value::Tensor(TensorValue {
            shape: vec![self.len()],
            data: self
                .into_iter()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: <dense::Mat as Output>::ty(),
        }))
    }
}
binding!(linear_codec_dense_norm,"num_linear_dense_norm","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::norm",DimensionPolicy::Preserve,(a:dense::Mat,kind:crate::kernels::dense::MatrixNorm)->Result<f64,OpError>,|a:dense::Mat,kind:crate::kernels::dense::MatrixNorm|{if a.rows>256||a.cols>256{return Err(OpError::Domain("dense spectral work budget exceeded".into()));}crate::kernels::dense::norm(&a,kind).map_err(dense_error)});
binding!(linear_codec_dense_cond_norm,"num_linear_dense_cond_norm","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::cond_norm",DimensionPolicy::UnitlessResult,(a:dense::Mat,kind:crate::kernels::dense::MatrixNorm)->Result<f64,OpError>,|a:dense::Mat,kind:crate::kernels::dense::MatrixNorm|{if a.rows>256||a.cols>256{return Err(OpError::Domain("dense spectral work budget exceeded".into()));}crate::kernels::dense::cond_norm(&a,kind).map_err(dense_error)});
binding!(linear_codec_dense_eigen,"num_linear_dense_eigen","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::eigen",DimensionPolicy::Dimensionless,(a:dense::Mat)->Result<crate::kernels::dense::Eigen,OpError>,|a:dense::Mat|{if a.rows>256||a.cols>256{return Err(OpError::Domain("dense spectral work budget exceeded".into()));}crate::kernels::dense::eigen(&a).map_err(dense_error)});
binding!(linear_codec_dense_eigenvalues,"num_linear_dense_eigenvalues","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::eigenvalues",DimensionPolicy::Dimensionless,(a:dense::Mat)->Result<Vec<(f64,f64)>,OpError>,|a:dense::Mat|{if a.rows>256||a.cols>256{return Err(OpError::Domain("dense spectral work budget exceeded".into()));}crate::kernels::dense::eigenvalues(&a).map_err(dense_error)});
binding!(linear_codec_dense_eigen_generalized_symmetric,"num_linear_dense_eigen_generalized_symmetric","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::eigen_generalized_symmetric",DimensionPolicy::Dimensionless,(a:dense::Mat,b:dense::Mat)->Result<crate::kernels::dense::GeneralizedEigen,OpError>,|a:dense::Mat,b:dense::Mat|{if a.rows>256||a.cols>256{return Err(OpError::Domain("dense spectral work budget exceeded".into()));}crate::kernels::dense::eigen_generalized_symmetric(&a,&b).map_err(dense_error)});
binding!(linear_codec_phi_matrices,"num_linear_phi_matrices","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::matfun::phi_matrices",DimensionPolicy::Dimensionless,(a:dense::Mat,p:usize)->Result<Vec<dense::Mat>,OpError>,|a:dense::Mat,p:usize|{if p>16||a.rows>64||a.cols>64||a.rows.checked_mul(p+1).is_none_or(|n|n>256){return Err(OpError::Domain("augmented exponential work budget exceeded".into()));}crate::kernels::matfun::phi_matrices(&a,p).map_err(dense_error)});
binding!(linear_codec_mat_zeros,"num_linear_mat_zeros","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::zeros",DimensionPolicy::Dimensionless,(rows:usize,cols:usize)->Result<dense::Mat,OpError>,|rows:usize,cols:usize|{if rows.checked_mul(cols).is_none_or(|n|n>262144){return Err(OpError::Domain("matrix allocation budget exceeded".into()));}Ok(dense::Mat::zeros(rows,cols))});
binding!(linear_codec_mat_identity,"num_linear_mat_identity","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::identity",DimensionPolicy::Dimensionless,(n:usize)->Result<dense::Mat,OpError>,|n:usize|{if n>512{return Err(OpError::Domain("matrix allocation budget exceeded".into()));}Ok(dense::Mat::identity(n))});
binding!(linear_codec_mat_from_rows,"num_linear_mat_from_rows","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::from_rows",DimensionPolicy::Dimensionless,(rows:Vec<Vec<f64>>)->Result<dense::Mat,OpError>,|rows:Vec<Vec<f64>>|{dense::Mat::from_rows(&rows).map_err(dense_error)});
binding!(linear_codec_mat_at,"num_linear_mat_at","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::at",DimensionPolicy::Preserve,(a:dense::Mat,i:usize,j:usize)->Result<f64,OpError>,|a:dense::Mat,i:usize,j:usize|{if i>=a.rows||j>=a.cols{return Err(OpError::Domain("matrix index out of range".into()));}Ok(a.at(i,j))});
binding!(linear_codec_mat_set,"num_linear_mat_set","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::set",DimensionPolicy::Dimensionless,(a:dense::Mat,i:usize,j:usize,value:f64)->Result<dense::Mat,OpError>,|a:dense::Mat,i:usize,j:usize,value:f64|{let mut a=a;if i>=a.rows||j>=a.cols{return Err(OpError::Domain("matrix index out of range".into()));}a.set(i,j,value);Ok(a)});
binding!(linear_codec_mat_is_finite,"num_linear_mat_is_finite","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::is_finite",DimensionPolicy::UnitlessResult,(a:dense::Mat)->bool,|a:dense::Mat|{a.is_finite()});
binding!(linear_codec_mat_norm_fro,"num_linear_mat_norm_fro","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::norm_fro",DimensionPolicy::Preserve,(a:dense::Mat)->f64,|a:dense::Mat|{a.norm_fro()});
binding!(linear_codec_mat_matvec,"num_linear_mat_matvec","numerical_matrix","Checked native linear calculation in normalized reference coordinates. Physical tensors must be explicitly nondimensionalized; normalized outputs use the same declared reference basis.","crate::kernels::dense::Mat::matvec",DimensionPolicy::Product,(a:dense::Mat,x:Vec<f64>)->Result<Vec<f64>,OpError>,|a:dense::Mat,x:Vec<f64>|{if x.len()!=a.cols{return Err(OpError::Domain("matrix-vector shape mismatch".into()));}Ok(a.matvec(&x))});
pub(super) static LINEAR_CODEC_EXPORTS: &[&CodecExport] = &[
    &linear_codec_optimized_det::EXPORT,
    &linear_codec_optimized_inv::EXPORT,
    &linear_codec_optimized_lstsq::EXPORT,
    &linear_codec_optimized_lu::EXPORT,
    &linear_codec_optimized_matmul::EXPORT,
    &linear_codec_optimized_qr::EXPORT,
    &linear_codec_optimized_solve::EXPORT,
    &linear_codec_reference_det::EXPORT,
    &linear_codec_reference_inv::EXPORT,
    &linear_codec_reference_lu::EXPORT,
    &linear_codec_reference_matmul::EXPORT,
    &linear_codec_reference_qr::EXPORT,
    &linear_codec_reference_solve::EXPORT,
    &linear_codec_reference_tril::EXPORT,
    &linear_codec_reference_triu::EXPORT,
    &linear_codec_exact_sqrt::EXPORT,
    &linear_codec_exact_lu::EXPORT,
    &linear_codec_exact_qr::EXPORT,
    &linear_codec_dense_norm::EXPORT,
    &linear_codec_dense_cond_norm::EXPORT,
    &linear_codec_dense_eigen::EXPORT,
    &linear_codec_dense_eigenvalues::EXPORT,
    &linear_codec_dense_eigen_generalized_symmetric::EXPORT,
    &linear_codec_phi_matrices::EXPORT,
    &linear_codec_mat_zeros::EXPORT,
    &linear_codec_mat_identity::EXPORT,
    &linear_codec_mat_from_rows::EXPORT,
    &linear_codec_mat_at::EXPORT,
    &linear_codec_mat_set::EXPORT,
    &linear_codec_mat_is_finite::EXPORT,
    &linear_codec_mat_norm_fro::EXPORT,
    &linear_codec_mat_matvec::EXPORT,
];
binding!(linear_codec_optimized_det,"num_linear_optimized_det","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::det",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<f64,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::optimized::det(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_optimized_inv,"num_linear_optimized_inv","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::inv",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::optimized::inv(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_optimized_lstsq,"num_linear_optimized_lstsq","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::lstsq",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}
if b.len()!=rows{return Err(OpError::Domain("RHS length differs from matrix rows".into()));}crate::kernels::dense::optimized::lstsq(rows,cols,&data,&b).map_err(dense_error)});
binding!(linear_codec_optimized_lu,"num_linear_optimized_lu","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::lu",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<crate::rules::numerical_abi::codecs::NativeDenseFactors<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::optimized::lu(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_optimized_matmul,"num_linear_optimized_matmul","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::matmul",DimensionPolicy::Dimensionless,(rows:usize,inner:usize,cols:usize,a:Vec<f64>,b:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,inner:usize,cols:usize,a:Vec<f64>,b:Vec<f64>|{if rows>256||inner>256||cols>256||rows.checked_mul(inner)!=Some(a.len())||inner.checked_mul(cols)!=Some(b.len()){return Err(OpError::Domain("product shape/work budget exceeded".into()));}crate::kernels::dense::optimized::matmul(rows,inner,cols,&a,&b).map_err(dense_error)});
binding!(linear_codec_optimized_qr,"num_linear_optimized_qr","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::qr",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<(Vec<f64>,Vec<f64>),OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::optimized::qr(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_optimized_solve,"num_linear_optimized_solve","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::optimized::solve",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}
if b.len()!=rows{return Err(OpError::Domain("RHS length differs from matrix rows".into()));}crate::kernels::dense::optimized::solve(rows,cols,&data,&b).map_err(dense_error)});
binding!(linear_codec_reference_det,"num_linear_reference_det","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::det",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<f64,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::reference::det(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_reference_inv,"num_linear_reference_inv","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::inv",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::reference::inv(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_reference_lu,"num_linear_reference_lu","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::lu",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<crate::rules::numerical_abi::codecs::NativeDenseFactors<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::reference::lu(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_reference_matmul,"num_linear_reference_matmul","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::matmul",DimensionPolicy::Dimensionless,(rows:usize,inner:usize,cols:usize,a:Vec<f64>,b:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,inner:usize,cols:usize,a:Vec<f64>,b:Vec<f64>|{if rows>256||inner>256||cols>256||rows.checked_mul(inner)!=Some(a.len())||inner.checked_mul(cols)!=Some(b.len()){return Err(OpError::Domain("product shape/work budget exceeded".into()));}Ok(crate::kernels::dense::reference::matmul(rows,inner,cols,&a,&b))});
binding!(linear_codec_reference_qr,"num_linear_reference_qr","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::qr",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<(Vec<f64>,Vec<f64>),OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}crate::kernels::dense::reference::qr(rows,cols,&data).map_err(dense_error)});
binding!(linear_codec_reference_solve,"num_linear_reference_solve","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::solve",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>,b:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}
if b.len()!=rows{return Err(OpError::Domain("RHS length differs from matrix rows".into()));}crate::kernels::dense::reference::solve(rows,cols,&data,&b).map_err(dense_error)});
binding!(linear_codec_reference_tril,"num_linear_reference_tril","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::tril",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}Ok(crate::kernels::dense::reference::tril(rows,cols,&data))});
binding!(linear_codec_reference_triu,"num_linear_reference_triu","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::reference::triu",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<f64>)->Result<Vec<f64>,OpError>,|rows:usize,cols:usize,data:Vec<f64>|{if rows>256||cols>256||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("row-major matrix shape/allocation budget exceeded".into()));}Ok(crate::kernels::dense::reference::triu(rows,cols,&data))});
binding!(linear_codec_exact_sqrt,"num_linear_exact_sqrt","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::exact::exact_sqrt",DimensionPolicy::Dimensionless,(value:num_rational::BigRational)->Option<num_rational::BigRational>,|value:num_rational::BigRational|{crate::kernels::dense::exact::exact_sqrt(&value)});
binding!(linear_codec_exact_lu,"num_linear_exact_lu","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::exact::lu_exact",DimensionPolicy::Dimensionless,(n:usize,data:Vec<num_rational::BigRational>)->Result<crate::rules::numerical_abi::codecs::NativeDenseFactors<num_rational::BigRational>,OpError>,|n:usize,data:Vec<num_rational::BigRational>|{if n>32||n.checked_mul(n)!=Some(data.len()){return Err(OpError::Domain("exact LU shape/work budget exceeded".into()));}crate::kernels::dense::exact::lu_exact(n,&data).map_err(|e|OpError::Domain(e.to_string()))});
binding!(linear_codec_exact_qr,"num_linear_exact_qr","numerical_matrix","Native row-major linear calculation in explicitly normalized reference coordinates. Shapes, coefficient/extent bounds and factor work limits are checked before native allocation.","crate::kernels::dense::exact::qr_exact",DimensionPolicy::Dimensionless,(rows:usize,cols:usize,data:Vec<num_rational::BigRational>)->Result<(Vec<num_rational::BigRational>,Vec<num_rational::BigRational>),OpError>,|rows:usize,cols:usize,data:Vec<num_rational::BigRational>|{if rows>32||cols>32||rows.checked_mul(cols)!=Some(data.len()){return Err(OpError::Domain("exact QR shape/work budget exceeded".into()));}crate::kernels::dense::exact::qr_exact(rows,cols,&data).map_err(|e|OpError::Domain(e.to_string()))});
