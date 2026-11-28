// Native normalized enums and value objects. Enum spellings are exact Rust variant names.
macro_rules! analysis_enum_arg {($ty:path; $($variant:ident),+ $(,)?)=>{
impl Arg for $ty {
 fn ty()->TypeExpr{TypeExpr::String}
 fn decode(value:&Value,index:usize)->Result<Self,OpError>{let s=<String as Arg>::decode(value,index)?;match s.as_str(){$(stringify!($variant)=>Ok(Self::$variant),)+_=>Err(OpError::Domain(format!("unknown {} variant: {}",stringify!($ty),s)))}}
}
};}
analysis_enum_arg!(crate::kernels::fft::Norm;Backward,Forward,Unitary);
analysis_enum_arg!(crate::kernels::convergence::LevinVariant;U,T,V);
analysis_enum_arg!(crate::kernels::coordinates::CoordinateSystem;Cartesian,Polar,Cylindrical,Spherical);
analysis_enum_arg!(crate::kernels::wavelets::Boundary;Periodic,Zero,Symmetric,Reflect);
analysis_enum_arg!(crate::kernels::wavelets::Wavelet;Haar,Db2,Db3,Db4,Db5,Db6,Db7,Db8,Db9,Db10,Sym4,Sym5,Sym6,Bior1_1,Bior1_3,Bior1_5,Bior2_2,Bior2_4,Bior3_3,Bior3_5,Bior4_4);
integer_codec!(isize);
impl Arg for crate::kernels::complex_scalar::Complex {
    fn ty() -> TypeExpr {
        TypeExpr::Complex
    }

    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        let z = <num_complex::Complex64 as Arg>::decode(v, i)?;
        Ok(Self::new(z.re, z.im))
    }
}
impl Output for crate::kernels::complex_scalar::Complex {
    fn ty() -> TypeExpr {
        TypeExpr::Complex
    }

    fn encode(self) -> Result<Value, OpError> {
        num_complex::Complex64::new(self.re, self.im).encode()
    }
}
impl Arg for crate::kernels::complex::Mobius {
    fn ty() -> TypeExpr {
        TypeExpr::Record(
            ["a", "b", "c", "d"]
                .into_iter()
                .map(|n| (n.into(), TypeExpr::Complex))
                .collect(),
        )
    }

    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        let f = |n| {
            <num_complex::Complex64 as Arg>::decode(model_field(v, n, &["a", "b", "c", "d"])?, i)
        };
        Ok(Self::new(f("a")?, f("b")?, f("c")?, f("d")?))
    }
}
record!(crate::kernels::complex::Mobius {
    a: num_complex::Complex64,
    b: num_complex::Complex64,
    c: num_complex::Complex64,
    d: num_complex::Complex64
});
record!(crate::kernels::signal::PowerSpectrum {
    frequencies: Vec<f64>,
    density: Vec<f64>
});
record!(crate::kernels::signal::CrossSpectrum {
    frequencies: Vec<f64>,
    density: Vec<num_complex::Complex64>,
    coherence: Vec<f64>
});
impl Arg for Vec<Vec<num_complex::Complex64>> {
    fn ty() -> TypeExpr {
        TypeExpr::Tensor {
            rank: NatExpr::Lit(2),
            shape: None,
            elem: Box::new(TypeExpr::Complex),
        }
    }

    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        let Value::Tensor(t) = v else {
            return Err(OpError::Domain("expected complex matrix".into()));
        };
        if t.shape.len() != 2
            || t.shape.iter().try_fold(1usize, |a, &b| a.checked_mul(b)) != Some(t.data.len())
            || t.data.len() > 8_388_608
        {
            return Err(OpError::Domain(
                "invalid complex matrix shape or allocation budget".into(),
            ));
        }
        let rows = t.shape[0];
        let cols = t.shape[1];
        if cols == 0 {
            return Ok(vec![vec![]; rows]);
        }
        t.data
            .chunks(cols)
            .map(|row| {
                row.iter()
                    .map(|x| <num_complex::Complex64 as Arg>::decode(x, i))
                    .collect()
            })
            .collect()
    }
}
impl Output for Vec<Vec<num_complex::Complex64>> {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let n = self.len();
        let m = self.first().map_or(0, Vec::len);
        if n.checked_mul(m).is_none_or(|n| n > 8_388_608) || self.iter().any(|r| r.len() != m) {
            return Err(OpError::Domain("ragged or excessive complex matrix".into()));
        }
        Ok(Value::Tensor(TensorValue {
            shape: vec![n, m],
            data: self
                .into_iter()
                .flatten()
                .map(Output::encode)
                .collect::<Result<_, _>>()?,
            elem: TypeExpr::Complex,
        }))
    }
}
// Ragged coefficient bands are encoded as flat samples plus checked offsets.
struct AnalysisWaveletBands(Vec<Vec<f64>>);
impl Arg for AnalysisWaveletBands {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("data".into(), <Vec<f64> as Arg>::ty()),
            ("offsets".into(), <Vec<usize> as Arg>::ty()),
        ])
    }

    fn decode(
        v: &Value,
        i: usize,
    ) -> Result<Self, OpError> {
        let names = ["data", "offsets"];
        let data = <Vec<f64> as Arg>::decode(model_field(v, "data", &names)?, i)?;
        let offsets = <Vec<usize> as Arg>::decode(model_field(v, "offsets", &names)?, i)?;
        if offsets.len() > 66
            || offsets.first() != Some(&0)
            || offsets.last() != Some(&data.len())
            || offsets.windows(2).any(|w| w[0] > w[1])
            || data.len() > 8_388_608
        {
            return Err(OpError::Domain(
                "invalid wavelet band offsets or execution budget".into(),
            ));
        }
        Ok(Self(
            offsets
                .windows(2)
                .map(|w| data[w[0]..w[1]].to_vec())
                .collect(),
        ))
    }
}
impl Output for AnalysisWaveletBands {
    fn ty() -> TypeExpr {
        <Self as Arg>::ty()
    }

    fn encode(self) -> Result<Value, OpError> {
        let mut offsets = vec![0usize];
        let mut data = Vec::new();
        for row in self.0 {
            if data
                .len()
                .checked_add(row.len())
                .is_none_or(|n| n > 8_388_608)
            {
                return Err(OpError::Domain(
                    "wavelet bands exceed allocation budget".into(),
                ));
            }
            data.extend(row);
            offsets.push(data.len());
        }
        Ok(Value::Record(RecordValue {
            fields: vec![
                ("data".into(), data.encode()?),
                ("offsets".into(), offsets.encode()?),
            ],
        }))
    }
}
binding!(analysis_complex_scalar_arg,"num_analysis_complex_scalar_arg","numerical_vector","Native complex_scalar arg; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::arg",DimensionPolicy::UnitlessResult,(value:crate::kernels::complex_scalar::Complex)->f64,|value:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::arg(value));
binding!(analysis_complex_scalar_magnitude,"num_analysis_complex_scalar_magnitude","numerical_vector","Native complex_scalar magnitude; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::magnitude",DimensionPolicy::Preserve,(value:crate::kernels::complex_scalar::Complex)->f64,|value:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::magnitude(value));
binding!(analysis_complex_scalar_norm_squared,"num_analysis_complex_scalar_norm_squared","numerical_vector","Native complex_scalar norm_squared; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::norm_squared",DimensionPolicy::Square,(value:crate::kernels::complex_scalar::Complex)->f64,|value:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::norm_squared(value));
binding!(analysis_complex_scalar_cos,"num_analysis_complex_scalar_cos","numerical_vector","Native complex_scalar cos; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::cos",DimensionPolicy::Dimensionless,(value:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|value:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::cos(value));
binding!(analysis_complex_scalar_sin,"num_analysis_complex_scalar_sin","numerical_vector","Native complex_scalar sin; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::sin",DimensionPolicy::Dimensionless,(value:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|value:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::sin(value));
binding!(analysis_complex_scalar_pow,"num_analysis_complex_scalar_pow","numerical_vector","Native complex_scalar pow; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::pow",DimensionPolicy::Dimensionless,(value:crate::kernels::complex_scalar::Complex, exponent:i64)->crate::kernels::complex_scalar::Complex,|value:crate::kernels::complex_scalar::Complex, exponent:i64|crate::kernels::complex_scalar::pow(value,exponent));
binding!(analysis_complex_scalar_powc,"num_analysis_complex_scalar_powc","numerical_vector","Native complex_scalar powc; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::powc",DimensionPolicy::Dimensionless,(base:crate::kernels::complex_scalar::Complex, exponent:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|base:crate::kernels::complex_scalar::Complex, exponent:crate::kernels::complex_scalar::Complex|crate::kernels::complex_scalar::powc(base,exponent));
binding!(analysis_complex_scalar_complex,"num_analysis_complex_scalar_complex","numerical_vector","Native complex_scalar complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::complex",DimensionPolicy::EqualPreserve,(re:f64, im:f64)->crate::kernels::complex_scalar::Complex,|re:f64, im:f64|crate::kernels::complex_scalar::complex(re,im));
binding!(analysis_complex_scalar_method_new,"num_analysis_complex_scalar_method_new","numerical_vector","Native complex_scalar method_new; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::new",DimensionPolicy::EqualPreserve,(re:f64, im:f64)->crate::kernels::complex_scalar::Complex,|re:f64, im:f64|crate::kernels::complex_scalar::Complex::new(re,im));
binding!(analysis_complex_scalar_method_abs,"num_analysis_complex_scalar_method_abs","numerical_vector","Native complex_scalar method_abs; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::abs",DimensionPolicy::Preserve,(z:crate::kernels::complex_scalar::Complex)->f64,|z:crate::kernels::complex_scalar::Complex|z.abs());
binding!(analysis_complex_scalar_method_arg,"num_analysis_complex_scalar_method_arg","numerical_vector","Native complex_scalar method_arg; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::arg",DimensionPolicy::UnitlessResult,(z:crate::kernels::complex_scalar::Complex)->f64,|z:crate::kernels::complex_scalar::Complex|z.arg());
binding!(analysis_complex_scalar_method_conj,"num_analysis_complex_scalar_method_conj","numerical_vector","Native complex_scalar method_conj; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::conj",DimensionPolicy::Preserve,(z:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex|z.conj());
binding!(analysis_complex_scalar_method_div,"num_analysis_complex_scalar_method_div","numerical_vector","Native complex_scalar method_div; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::div",DimensionPolicy::Quotient,(z:crate::kernels::complex_scalar::Complex, other:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex, other:crate::kernels::complex_scalar::Complex|z.div(other));
binding!(
    analysis_complex_scalar_method_recip,
    "num_analysis_complex_scalar_method_recip",
    "numerical_vector",
    "Computes the complex reciprocal 1/z with a scaled evaluation when intermediate products leave the ordinary floating-point range. Returns the inverse physical dimension of the input. Zero or a nonfinite result produces a domain error.",
    "crate::kernels::complex_scalar::Complex::recip",
    DimensionPolicy::Inverse,
    (z: crate::kernels::complex_scalar::Complex) -> crate::kernels::complex_scalar::Complex,
    crate::kernels::complex_scalar::Complex::recip
);
binding!(analysis_complex_scalar_method_mul,"num_analysis_complex_scalar_method_mul","numerical_vector","Native complex_scalar method_mul; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::mul",DimensionPolicy::Product,(z:crate::kernels::complex_scalar::Complex, other:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex, other:crate::kernels::complex_scalar::Complex|z.mul(other));
binding!(analysis_complex_scalar_method_exp,"num_analysis_complex_scalar_method_exp","numerical_vector","Native complex_scalar method_exp; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::exp",DimensionPolicy::Dimensionless,(z:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex|z.exp());
binding!(analysis_complex_scalar_method_ln,"num_analysis_complex_scalar_method_ln","numerical_vector","Native complex_scalar method_ln; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex_scalar::Complex::ln",DimensionPolicy::Dimensionless,(z:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex|z.ln());
binding!(analysis_complex_scalar_method_sqrt,"num_analysis_complex_scalar_method_sqrt","numerical_vector","Returns the principal complex square root, with the imaginary sign selecting the side of the negative real axis. Rescaling protects finite large inputs and the quotient preserves representable small components. The result carries the square root of the input dimension.","crate::kernels::complex_scalar::Complex::sqrt",DimensionPolicy::Root,(z:crate::kernels::complex_scalar::Complex)->crate::kernels::complex_scalar::Complex,|z:crate::kernels::complex_scalar::Complex|z.sqrt());
binding!(analysis_complex_mobius_new,"num_analysis_complex_mobius_new","numerical_vector","Native complex mobius_new; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex::Mobius::new",DimensionPolicy::Dimensionless,(a:num_complex::Complex64, b:num_complex::Complex64, c:num_complex::Complex64, d:num_complex::Complex64)->crate::kernels::complex::Mobius,|a:num_complex::Complex64, b:num_complex::Complex64, c:num_complex::Complex64, d:num_complex::Complex64|crate::kernels::complex::Mobius::new(a,b,c,d));
binding!(analysis_complex_mobius_identity,"num_analysis_complex_mobius_identity","numerical_vector","Native complex mobius_identity; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex::Mobius::identity",DimensionPolicy::Dimensionless,()->crate::kernels::complex::Mobius,crate::kernels::complex::Mobius::identity);
binding!(analysis_complex_mobius_apply,"num_analysis_complex_mobius_apply","numerical_vector","Native complex mobius_apply; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex::Mobius::apply",DimensionPolicy::Dimensionless,(map:crate::kernels::complex::Mobius, z:num_complex::Complex64)->num_complex::Complex64,|map:crate::kernels::complex::Mobius, z:num_complex::Complex64|map.apply(z));
binding!(analysis_complex_mobius_compose,"num_analysis_complex_mobius_compose","numerical_vector","Native complex mobius_compose; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex::Mobius::compose",DimensionPolicy::Dimensionless,(map:crate::kernels::complex::Mobius, other:crate::kernels::complex::Mobius)->crate::kernels::complex::Mobius,|map:crate::kernels::complex::Mobius, other:crate::kernels::complex::Mobius|map.compose(&other));
binding!(analysis_complex_mobius_inverse,"num_analysis_complex_mobius_inverse","numerical_vector","Native complex mobius_inverse; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::complex::Mobius::inverse",DimensionPolicy::Dimensionless,(map:crate::kernels::complex::Mobius)->crate::kernels::complex::Mobius,|map:crate::kernels::complex::Mobius|map.inverse());
binding!(analysis_fft_fft_norm,"num_analysis_fft_fft_norm","numerical_vector","Native fft fft_norm; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::fft::fft_norm",DimensionPolicy::Preserve,(input:Vec<num_complex::Complex64>, norm:crate::kernels::fft::Norm)->Vec<num_complex::Complex64>,|input:Vec<num_complex::Complex64>, norm:crate::kernels::fft::Norm|crate::kernels::fft::fft_norm(&input,norm));
binding!(analysis_fft_ifft_norm,"num_analysis_fft_ifft_norm","numerical_vector","Native fft ifft_norm; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::fft::ifft_norm",DimensionPolicy::Preserve,(input:Vec<num_complex::Complex64>, norm:crate::kernels::fft::Norm)->Vec<num_complex::Complex64>,|input:Vec<num_complex::Complex64>, norm:crate::kernels::fft::Norm|crate::kernels::fft::ifft_norm(&input,norm));
binding!(analysis_fft_fft_axes,"num_analysis_fft_fft_axes","numerical_vector","Native fft fft_axes; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::fft::fft_axes",DimensionPolicy::Preserve,(input:Vec<num_complex::Complex64>, shape:Vec<usize>, axes:Vec<usize>, inverse:bool, norm:crate::kernels::fft::Norm)->Result<Vec<num_complex::Complex64>, OpError>,|input:Vec<num_complex::Complex64>, shape:Vec<usize>, axes:Vec<usize>, inverse:bool, norm:crate::kernels::fft::Norm|crate::kernels::fft::fft_axes(&input,&shape,&axes,inverse,norm).map_err(OpError::Domain));
binding!(analysis_fftconv_fft,"num_analysis_fftconv_fft","numerical_vector","Native fftconv fft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::fftconv::fft",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>, inverse:bool)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>, inverse:bool|crate::kernels::fftconv::fft(&x.into_iter().map(|z|(z.re,z.im)).collect::<Vec<_>>(),inverse).into_iter().map(|(re,im)|num_complex::Complex64::new(re,im)).collect());
binding!(analysis_signal_fft,"num_analysis_signal_fft","numerical_vector","Native signal fft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::signal::fft",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::signal::fft(&mut x)});
binding!(analysis_signal_ifft,"num_analysis_signal_ifft","numerical_vector","Native signal ifft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::signal::ifft",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::signal::ifft(&mut x)});
binding!(analysis_signal_stft,"num_analysis_signal_stft","numerical_vector","Native signal stft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::signal::stft",DimensionPolicy::Preserve,(x:Vec<f64>, window:Vec<f64>, hop:usize)->Result<Vec<Vec<num_complex::Complex64>>,OpError>,|x:Vec<f64>, window:Vec<f64>, hop:usize|{if hop>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(crate::kernels::signal::stft(&x,&window,hop))});
binding!(analysis_transforms_fft,"num_analysis_transforms_fft","numerical_vector","Native transforms fft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::fft",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::transforms::fft(&mut x);x});
binding!(analysis_transforms_ifft,"num_analysis_transforms_ifft","numerical_vector","Native transforms ifft; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::ifft",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::transforms::ifft(&mut x);x});
binding!(analysis_transforms_fft_slice,"num_analysis_transforms_fft_slice","numerical_vector","Native transforms fft_slice; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::fft_slice",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::transforms::fft_slice(&mut x);x});
binding!(analysis_transforms_ifft_slice,"num_analysis_transforms_ifft_slice","numerical_vector","Native transforms ifft_slice; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::ifft_slice",DimensionPolicy::Preserve,(x:Vec<num_complex::Complex64>)->Vec<num_complex::Complex64>,|x:Vec<num_complex::Complex64>|{let mut x=x;crate::kernels::transforms::ifft_slice(&mut x);x});
binding!(analysis_transforms_fft_many_par,"num_analysis_transforms_fft_many_par","numerical_vector","Native transforms fft_many_par; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::fft_many_par",DimensionPolicy::Preserve,(rows:Vec<Vec<num_complex::Complex64>>, grain:usize)->Result<Vec<Vec<num_complex::Complex64>>,OpError>,|rows:Vec<Vec<num_complex::Complex64>>, grain:usize|{if grain>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(crate::kernels::transforms::fft_many_par(&rows,grain))});
binding!(analysis_transforms_ifft_many_par,"num_analysis_transforms_ifft_many_par","numerical_vector","Native transforms ifft_many_par; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::transforms::ifft_many_par",DimensionPolicy::Preserve,(rows:Vec<Vec<num_complex::Complex64>>, grain:usize)->Result<Vec<Vec<num_complex::Complex64>>,OpError>,|rows:Vec<Vec<num_complex::Complex64>>, grain:usize|{if grain>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(crate::kernels::transforms::ifft_many_par(&rows,grain))});
binding!(analysis_hypergeometric_gamma_c,"num_analysis_hypergeometric_gamma_c","numerical_vector","Native hypergeometric gamma_c; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::gamma_c",DimensionPolicy::Dimensionless,(z:num_complex::Complex64)->num_complex::Complex64,|z:num_complex::Complex64|crate::kernels::hypergeometric::gamma_c(z));
binding!(analysis_hypergeometric_ln_gamma_c,"num_analysis_hypergeometric_ln_gamma_c","numerical_vector","Native hypergeometric ln_gamma_c; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::ln_gamma_c",DimensionPolicy::Dimensionless,(z:num_complex::Complex64)->num_complex::Complex64,|z:num_complex::Complex64|crate::kernels::hypergeometric::ln_gamma_c(z));
binding!(analysis_hypergeometric_hyp1f1_cplx,"num_analysis_hypergeometric_hyp1f1_cplx","numerical_vector","Native hypergeometric hyp1f1_cplx; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::hyp1f1_cplx",DimensionPolicy::Dimensionless,(a:num_complex::Complex64, b:num_complex::Complex64, z:num_complex::Complex64)->num_complex::Complex64,|a:num_complex::Complex64, b:num_complex::Complex64, z:num_complex::Complex64|crate::kernels::hypergeometric::hyp1f1_cplx(a,b,z));
binding!(analysis_hypergeometric_hyp2f1_cplx,"num_analysis_hypergeometric_hyp2f1_cplx","numerical_vector","Native hypergeometric hyp2f1_cplx; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::hyp2f1_cplx",DimensionPolicy::Dimensionless,(a:num_complex::Complex64, b:num_complex::Complex64, c:num_complex::Complex64, z:num_complex::Complex64)->num_complex::Complex64,|a:num_complex::Complex64, b:num_complex::Complex64, c:num_complex::Complex64, z:num_complex::Complex64|crate::kernels::hypergeometric::hyp2f1_cplx(a,b,c,z));
binding!(analysis_hypergeometric_hyp_u_cplx,"num_analysis_hypergeometric_hyp_u_cplx","numerical_vector","Native hypergeometric hyp_u_cplx; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::hyp_u_cplx",DimensionPolicy::Dimensionless,(a:num_complex::Complex64, b:num_complex::Complex64, z:num_complex::Complex64)->num_complex::Complex64,|a:num_complex::Complex64, b:num_complex::Complex64, z:num_complex::Complex64|crate::kernels::hypergeometric::hyp_u_cplx(a,b,z));
binding!(analysis_hypergeometric_hyp_pfq,"num_analysis_hypergeometric_hyp_pfq","numerical_vector","Native hypergeometric hyp_pfq; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::hyp_pfq",DimensionPolicy::Dimensionless,(a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64)->num_complex::Complex64,|a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64|crate::kernels::hypergeometric::hyp_pfq(&a,&b,z));
binding!(analysis_hypergeometric_hyp_pfq_checked,"num_analysis_hypergeometric_hyp_pfq_checked","numerical_vector","Native hypergeometric hyp_pfq_checked; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::hypergeometric::hyp_pfq_checked",DimensionPolicy::Dimensionless,(a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64)->Result<num_complex::Complex64, OpError>,|a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64|crate::kernels::hypergeometric::hyp_pfq_checked(&a,&b,z).map_err(|e|OpError::Domain(e.to_string())));
binding!(analysis_special_complex_bessel_i_complex,"num_analysis_special_complex_bessel_i_complex","numerical_vector","Native special_complex bessel_i_complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::special_complex::bessel_i_complex",DimensionPolicy::Dimensionless,(nu:f64, z:num_complex::Complex64)->num_complex::Complex64,|nu:f64, z:num_complex::Complex64|crate::kernels::special_complex::bessel_i_complex(nu,z));
binding!(analysis_special_complex_bessel_j_complex,"num_analysis_special_complex_bessel_j_complex","numerical_vector","Native special_complex bessel_j_complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::special_complex::bessel_j_complex",DimensionPolicy::Dimensionless,(nu:f64, z:num_complex::Complex64)->num_complex::Complex64,|nu:f64, z:num_complex::Complex64|crate::kernels::special_complex::bessel_j_complex(nu,z));
binding!(analysis_special_complex_digamma_complex,"num_analysis_special_complex_digamma_complex","numerical_vector","Native special_complex digamma_complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::special_complex::digamma_complex",DimensionPolicy::Dimensionless,(z:num_complex::Complex64)->num_complex::Complex64,|z:num_complex::Complex64|crate::kernels::special_complex::digamma_complex(z));
binding!(analysis_special_complex_gamma_complex,"num_analysis_special_complex_gamma_complex","numerical_vector","Native special_complex gamma_complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::special_complex::gamma_complex",DimensionPolicy::Dimensionless,(z:num_complex::Complex64)->num_complex::Complex64,|z:num_complex::Complex64|crate::kernels::special_complex::gamma_complex(z));
binding!(analysis_special_complex_ln_gamma_complex,"num_analysis_special_complex_ln_gamma_complex","numerical_vector","Native special_complex ln_gamma_complex; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::special_complex::ln_gamma_complex",DimensionPolicy::Dimensionless,(z:num_complex::Complex64)->num_complex::Complex64,|z:num_complex::Complex64|crate::kernels::special_complex::ln_gamma_complex(z));
binding!(analysis_special_complex_lambert_w_complex,"num_analysis_special_complex_lambert_w_complex","numerical_vector","Return the dimensionless complex solution to w exp(w) = z on integer branch k. When the imaginary component is signed zero on the negative real axis, its sign selects the side of the branch cut. The principal branch preserves small real and imaginary components near the origin.","crate::kernels::special_complex::lambert_w_complex",DimensionPolicy::Dimensionless,(z:num_complex::Complex64, k:i32)->num_complex::Complex64,|z:num_complex::Complex64, k:i32|crate::kernels::special_complex::lambert_w_complex(z,k));
binding!(analysis_meijer_meijer_g,"num_analysis_meijer_meijer_g","numerical_vector","Native meijer meijer_g; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::meijer::meijer_g",DimensionPolicy::Dimensionless,(m:usize, n:usize, a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64)->Result<Option<num_complex::Complex64>,OpError>,|m:usize, n:usize, a:Vec<num_complex::Complex64>, b:Vec<num_complex::Complex64>, z:num_complex::Complex64|{if m>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if n>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(crate::kernels::meijer::meijer_g(m,n,&a,&b,z))});
binding!(analysis_meijer_meijer_g_real,"num_analysis_meijer_meijer_g_real","numerical_vector","Native meijer meijer_g_real; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::meijer::meijer_g_real",DimensionPolicy::Dimensionless,(m:usize, n:usize, a:Vec<f64>, b:Vec<f64>, z:f64)->Result<Option<f64>,OpError>,|m:usize, n:usize, a:Vec<f64>, b:Vec<f64>, z:f64|{if m>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if n>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(crate::kernels::meijer::meijer_g_real(m,n,&a,&b,z))});
binding!(analysis_convergence_levin_condition_number,"num_analysis_convergence_levin_condition_number","numerical_vector","Native convergence levin_condition_number; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::convergence::levin_condition_number",DimensionPolicy::UnitlessResult,(terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64)->f64,|terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64|crate::kernels::convergence::levin_condition_number(&terms,variant,beta));
binding!(
    analysis_convergence_levin_condition_numbers,
    "num_analysis_convergence_levin_condition_numbers",
    "numerical_vector",
    "Returns one dimensionless Levin condition number per available order, from order zero. Infinite values indicate an undefined or unbounded condition number. Measured terms may share any dimension, beta is dimensionless, and the input budget is 4096 terms.",
    "crate::kernels::convergence::levin_condition_numbers",
    DimensionPolicy::UnitlessResult,
    (terms: Vec<f64>, variant: crate::kernels::convergence::LevinVariant, beta: f64) -> Result<Vec<f64>, OpError>,
    |terms: Vec<f64>, variant: crate::kernels::convergence::LevinVariant, beta: f64| {
        if terms.len() > 4096 {
            return Err(OpError::Domain("Levin condition profile exceeds its term budget".into()));
        }
        Ok(crate::kernels::convergence::levin_condition_numbers(&terms, variant, beta))
    }
);
binding!(analysis_convergence_levin_transform,"num_analysis_convergence_levin_transform","numerical_vector","Native convergence levin_transform; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::convergence::levin_transform",DimensionPolicy::Preserve,(terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64)->Vec<f64>,|terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64|crate::kernels::convergence::levin_transform(&terms,variant,beta));
binding!(analysis_convergence_levin_transform_terms,"num_analysis_convergence_levin_transform_terms","numerical_vector","Native convergence levin_transform_terms; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::convergence::levin_transform_terms",DimensionPolicy::Preserve,(terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64)->Vec<f64>,|terms:Vec<f64>, variant:crate::kernels::convergence::LevinVariant, beta:f64|crate::kernels::convergence::levin_transform_terms(&terms,variant,beta));
binding!(analysis_convergence_levin_sum_lazy,"num_analysis_convergence_levin_sum_lazy","numerical_vector","Native convergence levin_sum_lazy; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::convergence::levin_sum_lazy",DimensionPolicy::WaveletSplit,(terms:Vec<f64>, max_terms:usize, variant:crate::kernels::convergence::LevinVariant, beta:f64)->Result<(f64,f64), OpError>,|terms:Vec<f64>, max_terms:usize, variant:crate::kernels::convergence::LevinVariant, beta:f64|{if max_terms>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::convergence::levin_sum_lazy(terms.into_iter(),max_terms,variant,beta).map_err(OpError::Domain)});
binding!(
    analysis_coordinates_cartesian_to_spherical,
    "num_analysis_coordinates_cartesian_to_spherical",
    "numerical_vector",
    "Converts the fixed shape [x, y, z] to [radius, polar angle from +z, azimuth from +x]. Express all lengths in one reference and use radians for both angles. The polar angle is zero at the origin.",
    "crate::kernels::coordinates::cartesian_to_spherical",
    DimensionPolicy::Dimensionless,
    (point: [f64; 3]) -> [f64; 3],
    crate::kernels::coordinates::cartesian_to_spherical
);
binding!(
    analysis_coordinates_spherical_to_cartesian,
    "num_analysis_coordinates_spherical_to_cartesian",
    "numerical_vector",
    "Converts [radius, polar angle from +z, azimuth from +x] to [x, y, z]. The incoming angles are radians and the result uses the same length reference as the radius. A negative radius reflects the point through the origin.",
    "crate::kernels::coordinates::spherical_to_cartesian",
    DimensionPolicy::Dimensionless,
    (point: [f64; 3]) -> [f64; 3],
    crate::kernels::coordinates::spherical_to_cartesian
);
binding!(analysis_coordinates_from_cartesian,"num_analysis_coordinates_from_cartesian","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::from_cartesian",DimensionPolicy::Dimensionless,(point:Vec<f64>, to:crate::kernels::coordinates::CoordinateSystem)->Result<Vec<f64>, OpError>,|point:Vec<f64>, to:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::from_cartesian(&point,to).map_err(OpError::Domain));
binding!(analysis_coordinates_from_cartesian_pure,"num_analysis_coordinates_from_cartesian_pure","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::from_cartesian_pure",DimensionPolicy::Dimensionless,(point:Vec<f64>, to:crate::kernels::coordinates::CoordinateSystem)->Vec<f64>,|point:Vec<f64>, to:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::from_cartesian_pure(&point,to));
binding!(analysis_coordinates_to_cartesian,"num_analysis_coordinates_to_cartesian","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::to_cartesian",DimensionPolicy::Dimensionless,(point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem)->Result<Vec<f64>, OpError>,|point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::to_cartesian(&point,from).map_err(OpError::Domain));
binding!(analysis_coordinates_to_cartesian_pure,"num_analysis_coordinates_to_cartesian_pure","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::to_cartesian_pure",DimensionPolicy::Dimensionless,(point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem)->Vec<f64>,|point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::to_cartesian_pure(&point,from));
binding!(analysis_coordinates_transform_point,"num_analysis_coordinates_transform_point","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::transform_point",DimensionPolicy::Dimensionless,(point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem)->Result<Vec<f64>, OpError>,|point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::transform_point(&point,from,to).map_err(OpError::Domain));
binding!(analysis_coordinates_transform_point_pure,"num_analysis_coordinates_transform_point_pure","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::transform_point_pure",DimensionPolicy::Dimensionless,(point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem)->Vec<f64>,|point:Vec<f64>, from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem|crate::kernels::coordinates::transform_point_pure(&point,from,to));
binding!(analysis_coordinates_jacobian_determinant,"num_analysis_coordinates_jacobian_determinant","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::jacobian_determinant",DimensionPolicy::Dimensionless,(from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem, at_point:Vec<f64>)->Result<f64, OpError>,|from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem, at_point:Vec<f64>|crate::kernels::coordinates::jacobian_determinant(from,to,&at_point).map_err(OpError::Domain));
binding!(analysis_coordinates_numerical_jacobian,"num_analysis_coordinates_numerical_jacobian","numerical_vector","Normalized coordinate conversion: radial/Cartesian coordinates are expressed in one chosen positive reference length, angles in radians. Result coordinates and Jacobians use the same normalization; mixed physical units are not erased silently.","crate::kernels::coordinates::numerical_jacobian",DimensionPolicy::Dimensionless,(from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem, at_point:Vec<f64>)->Result<Vec<Vec<f64>>, OpError>,|from:crate::kernels::coordinates::CoordinateSystem, to:crate::kernels::coordinates::CoordinateSystem, at_point:Vec<f64>|crate::kernels::coordinates::numerical_jacobian(from,to,&at_point).map_err(OpError::Domain));
binding!(analysis_wavelets_method_analysis_high,"num_analysis_wavelets_method_analysis_high","numerical_vector","Native wavelets method_analysis_high; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::analysis_high",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.analysis_high());
binding!(analysis_wavelets_method_analysis_low,"num_analysis_wavelets_method_analysis_low","numerical_vector","Native wavelets method_analysis_low; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::analysis_low",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.analysis_low());
binding!(analysis_wavelets_method_high_pass,"num_analysis_wavelets_method_high_pass","numerical_vector","Native wavelets method_high_pass; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::high_pass",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.high_pass());
binding!(analysis_wavelets_method_low_pass,"num_analysis_wavelets_method_low_pass","numerical_vector","Native wavelets method_low_pass; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::low_pass",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.low_pass());
binding!(analysis_wavelets_method_synthesis_high,"num_analysis_wavelets_method_synthesis_high","numerical_vector","Native wavelets method_synthesis_high; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::synthesis_high",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.synthesis_high());
binding!(analysis_wavelets_method_synthesis_low,"num_analysis_wavelets_method_synthesis_low","numerical_vector","Native wavelets method_synthesis_low; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::synthesis_low",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->Vec<f64>,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.synthesis_low());
binding!(analysis_wavelets_method_coeff_len,"num_analysis_wavelets_method_coeff_len","numerical_vector","Native wavelets method_coeff_len; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::coeff_len",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet, n:usize, boundary:crate::kernels::wavelets::Boundary)->Result<usize,OpError>,|wavelet:crate::kernels::wavelets::Wavelet, n:usize, boundary:crate::kernels::wavelets::Boundary|{if n>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}Ok(wavelet.coeff_len(n,boundary))});
binding!(analysis_wavelets_method_taps,"num_analysis_wavelets_method_taps","numerical_vector","Native wavelets method_taps; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::taps",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->usize,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.taps());
binding!(analysis_wavelets_method_is_orthogonal,"num_analysis_wavelets_method_is_orthogonal","numerical_vector","Native wavelets method_is_orthogonal; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::is_orthogonal",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->bool,|wavelet:crate::kernels::wavelets::Wavelet|wavelet.is_orthogonal());
binding!(analysis_wavelets_method_vanishing_moments,"num_analysis_wavelets_method_vanishing_moments","numerical_vector","Native wavelets method_vanishing_moments; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::Wavelet::vanishing_moments",DimensionPolicy::Dimensionless,(wavelet:crate::kernels::wavelets::Wavelet)->(usize,usize),|wavelet:crate::kernels::wavelets::Wavelet|wavelet.vanishing_moments());
binding!(analysis_wavelets_dwt,"num_analysis_wavelets_dwt","numerical_vector","Native wavelets dwt; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::dwt",DimensionPolicy::WaveletSplit,(x:Vec<f64>, w:crate::kernels::wavelets::Wavelet)->Result<(Vec<f64>,Vec<f64>), OpError>,|x:Vec<f64>, w:crate::kernels::wavelets::Wavelet|crate::kernels::wavelets::dwt(&x,w).map_err(|e|OpError::Domain(format!("{e:?}"))));
binding!(analysis_wavelets_dwt_with,"num_analysis_wavelets_dwt_with","numerical_vector","Native wavelets dwt_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::dwt_with",DimensionPolicy::WaveletSplit,(x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary)->Result<(Vec<f64>,Vec<f64>), OpError>,|x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary|crate::kernels::wavelets::dwt_with(&x,w,b).map_err(|e|OpError::Domain(format!("{e:?}"))));
binding!(analysis_wavelets_idwt,"num_analysis_wavelets_idwt","numerical_vector","Native wavelets idwt; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::idwt",DimensionPolicy::EqualPreserve,(a:Vec<f64>, d:Vec<f64>, w:crate::kernels::wavelets::Wavelet)->Result<Vec<f64>, OpError>,|a:Vec<f64>, d:Vec<f64>, w:crate::kernels::wavelets::Wavelet|crate::kernels::wavelets::idwt(&a,&d,w).map_err(|e|OpError::Domain(format!("{e:?}"))));
binding!(analysis_wavelets_idwt_with,"num_analysis_wavelets_idwt_with","numerical_vector","Native wavelets idwt_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::idwt_with",DimensionPolicy::EqualPreserve,(a:Vec<f64>, d:Vec<f64>, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary, n:usize)->Result<Vec<f64>, OpError>,|a:Vec<f64>, d:Vec<f64>, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary, n:usize|{if n>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::idwt_with(&a,&d,w,b,n).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_dwt2,"num_analysis_wavelets_dwt2","numerical_vector","Native wavelets dwt2; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::dwt2",DimensionPolicy::Preserve,(x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary)->Result<Vec<f64>, OpError>,|x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::dwt2(&x,rows,cols,w,b).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_idwt2,"num_analysis_wavelets_idwt2","numerical_vector","Native wavelets idwt2; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::idwt2",DimensionPolicy::Preserve,(c:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary)->Result<Vec<f64>, OpError>,|c:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::idwt2(&c,rows,cols,w,b).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_extend,"num_analysis_wavelets_extend","numerical_vector","Native wavelets extend; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::extend",DimensionPolicy::Preserve,(x:Vec<f64>, i:isize, b:crate::kernels::wavelets::Boundary)->f64,|x:Vec<f64>, i:isize, b:crate::kernels::wavelets::Boundary|crate::kernels::wavelets::extend(&x,i,b));
binding!(analysis_wavelets_wavedec,"num_analysis_wavelets_wavedec","numerical_vector","Native wavelets wavedec; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::wavedec",DimensionPolicy::AnalysisWaveletBandsBuild,(x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, levels:usize)->Result<super::AnalysisWaveletBands, OpError>,|x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, levels:usize|{if levels>64 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::wavedec(&x,w,levels).map(super::AnalysisWaveletBands).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_wavedec_with,"num_analysis_wavelets_wavedec_with","numerical_vector","Native wavelets wavedec_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::wavedec_with",DimensionPolicy::AnalysisWaveletBandsBuild,(x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, levels:usize, b:crate::kernels::wavelets::Boundary)->Result<super::AnalysisWaveletBands, OpError>,|x:Vec<f64>, w:crate::kernels::wavelets::Wavelet, levels:usize, b:crate::kernels::wavelets::Boundary|{if levels>64 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::wavedec_with(&x,w,levels,b).map(super::AnalysisWaveletBands).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_wavedec2,"num_analysis_wavelets_wavedec2","numerical_vector","Native wavelets wavedec2; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::wavedec2",DimensionPolicy::AnalysisWaveletBandsBuild,(x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, levels:usize)->Result<super::AnalysisWaveletBands, OpError>,|x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, levels:usize|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if levels>64 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::wavedec2(&x,rows,cols,w,levels).map(super::AnalysisWaveletBands).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_wavedec2_with,"num_analysis_wavelets_wavedec2_with","numerical_vector","Native wavelets wavedec2_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::wavedec2_with",DimensionPolicy::AnalysisWaveletBandsBuild,(x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, levels:usize, b:crate::kernels::wavelets::Boundary)->Result<super::AnalysisWaveletBands, OpError>,|x:Vec<f64>, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, levels:usize, b:crate::kernels::wavelets::Boundary|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if levels>64 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::wavedec2_with(&x,rows,cols,w,levels,b).map(super::AnalysisWaveletBands).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_waverec,"num_analysis_wavelets_waverec","numerical_vector","Native wavelets waverec; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::waverec",DimensionPolicy::AnalysisWaveletBandsRecover,(coeffs:super::AnalysisWaveletBands, w:crate::kernels::wavelets::Wavelet)->Result<Vec<f64>, OpError>,|coeffs:super::AnalysisWaveletBands, w:crate::kernels::wavelets::Wavelet|crate::kernels::wavelets::waverec(&coeffs.0,w).map_err(|e|OpError::Domain(format!("{e:?}"))));
binding!(analysis_wavelets_waverec_with,"num_analysis_wavelets_waverec_with","numerical_vector","Native wavelets waverec_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::waverec_with",DimensionPolicy::AnalysisWaveletBandsRecover,(coeffs:super::AnalysisWaveletBands, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary, n:usize)->Result<Vec<f64>, OpError>,|coeffs:super::AnalysisWaveletBands, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary, n:usize|{if n>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::waverec_with(&coeffs.0,w,b,n).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_waverec2,"num_analysis_wavelets_waverec2","numerical_vector","Native wavelets waverec2; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::waverec2",DimensionPolicy::AnalysisWaveletBandsRecover,(coeffs:super::AnalysisWaveletBands, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet)->Result<Vec<f64>, OpError>,|coeffs:super::AnalysisWaveletBands, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::waverec2(&coeffs.0,rows,cols,w).map_err(|e|OpError::Domain(format!("{e:?}")))});
binding!(analysis_wavelets_waverec2_with,"num_analysis_wavelets_waverec2_with","numerical_vector","Native wavelets waverec2_with; normalized dimensionless arguments unless the declared contract preserves measured signal units. Enum arguments use exact Rust variant names.","crate::kernels::wavelets::waverec2_with",DimensionPolicy::AnalysisWaveletBandsRecover,(coeffs:super::AnalysisWaveletBands, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary)->Result<Vec<f64>, OpError>,|coeffs:super::AnalysisWaveletBands, rows:usize, cols:usize, w:crate::kernels::wavelets::Wavelet, b:crate::kernels::wavelets::Boundary|{if rows>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}
if cols>1000000 {return Err(OpError::Domain("analysis count exceeds execution budget".into()));}crate::kernels::wavelets::waverec2_with(&coeffs.0,rows,cols,w,b).map_err(|e|OpError::Domain(format!("{e:?}")))});
include!("polynomial.rs");
include!("optimization.rs");
include!("homotopy.rs");

analysis_enum_arg!(crate::kernels::rootfind::LeastSquaresStatus;ResidualConverged,GradientConverged,IterationLimit,EvaluationLimit,Stagnated);
impl Arg for crate::kernels::rootfind::LeastSquaresResult {
    fn ty() -> TypeExpr {
        TypeExpr::Record(vec![
            ("x".into(), <Vec<f64> as Arg>::ty()),
            ("residual_norm".into(), TypeExpr::Float),
            ("iterations".into(), TypeExpr::Int),
            ("evaluations".into(), TypeExpr::Int),
            ("jacobian_evaluations".into(), TypeExpr::Int),
            ("accepted_steps".into(), TypeExpr::Int),
            ("rejected_steps".into(), TypeExpr::Int),
            ("status".into(), TypeExpr::String),
        ])
    }

    fn decode(
        v: &Value,
        index: usize,
    ) -> Result<Self, OpError> {
        let field = |name| crate::rules::numerical_abi::native_object_api::field(v, name);
        let n = |name| <usize as Arg>::decode(field(name)?, index);
        Ok(Self {
            x: <Vec<f64> as Arg>::decode(field("x")?, index)?,
            residual_norm: <f64 as Arg>::decode(field("residual_norm")?, index)?,
            iterations: n("iterations")?,
            evaluations: n("evaluations")?,
            jacobian_evaluations: n("jacobian_evaluations")?,
            accepted_steps: n("accepted_steps")?,
            rejected_steps: n("rejected_steps")?,
            status: <crate::kernels::rootfind::LeastSquaresStatus as Arg>::decode(
                field("status")?,
                index,
            )?,
        })
    }
}
binding!(analysis_lsq_status_converged,"num_analysis_lsq_status_converged","numerical_vector","Actual native least-squares status predicate; budget/stagnation statuses never imply convergence.","crate::kernels::rootfind::LeastSquaresStatus::converged",DimensionPolicy::Dimensionless,(status:crate::kernels::rootfind::LeastSquaresStatus)->bool,|status:crate::kernels::rootfind::LeastSquaresStatus|status.converged());
binding!(analysis_lsq_result_converged,"num_analysis_lsq_result_converged","numerical_vector","Actual native normalized least-squares report convergence predicate. All original report fields/status are supplied explicitly; this does not synthesize convergence from an iteration count.","crate::kernels::rootfind::LeastSquaresResult::converged",DimensionPolicy::Dimensionless,(report:crate::kernels::rootfind::LeastSquaresResult)->bool,|report:crate::kernels::rootfind::LeastSquaresResult|report.converged());
pub(super) static ANALYSIS_EXPORTS: &[&CodecExport] = &[
    &analysis_lsq_status_converged::EXPORT,
    &analysis_lsq_result_converged::EXPORT,
    &analysis_complex_scalar_arg::EXPORT,
    &analysis_complex_scalar_magnitude::EXPORT,
    &analysis_complex_scalar_norm_squared::EXPORT,
    &analysis_complex_scalar_cos::EXPORT,
    &analysis_complex_scalar_sin::EXPORT,
    &analysis_complex_scalar_pow::EXPORT,
    &analysis_complex_scalar_powc::EXPORT,
    &analysis_complex_scalar_complex::EXPORT,
    &analysis_complex_scalar_method_new::EXPORT,
    &analysis_complex_scalar_method_abs::EXPORT,
    &analysis_complex_scalar_method_arg::EXPORT,
    &analysis_complex_scalar_method_conj::EXPORT,
    &analysis_complex_scalar_method_div::EXPORT,
    &analysis_complex_scalar_method_recip::EXPORT,
    &analysis_complex_scalar_method_mul::EXPORT,
    &analysis_complex_scalar_method_exp::EXPORT,
    &analysis_complex_scalar_method_ln::EXPORT,
    &analysis_complex_scalar_method_sqrt::EXPORT,
    &analysis_complex_mobius_new::EXPORT,
    &analysis_complex_mobius_identity::EXPORT,
    &analysis_complex_mobius_apply::EXPORT,
    &analysis_complex_mobius_compose::EXPORT,
    &analysis_complex_mobius_inverse::EXPORT,
    &analysis_fft_fft_norm::EXPORT,
    &analysis_fft_ifft_norm::EXPORT,
    &analysis_fft_fft_axes::EXPORT,
    &analysis_fftconv_fft::EXPORT,
    &analysis_signal_fft::EXPORT,
    &analysis_signal_ifft::EXPORT,
    &analysis_signal_stft::EXPORT,
    &analysis_transforms_fft::EXPORT,
    &analysis_transforms_ifft::EXPORT,
    &analysis_transforms_fft_slice::EXPORT,
    &analysis_transforms_ifft_slice::EXPORT,
    &analysis_transforms_fft_many_par::EXPORT,
    &analysis_transforms_ifft_many_par::EXPORT,
    &analysis_hypergeometric_gamma_c::EXPORT,
    &analysis_hypergeometric_ln_gamma_c::EXPORT,
    &analysis_hypergeometric_hyp1f1_cplx::EXPORT,
    &analysis_hypergeometric_hyp2f1_cplx::EXPORT,
    &analysis_hypergeometric_hyp_u_cplx::EXPORT,
    &analysis_hypergeometric_hyp_pfq::EXPORT,
    &analysis_hypergeometric_hyp_pfq_checked::EXPORT,
    &analysis_special_complex_bessel_i_complex::EXPORT,
    &analysis_special_complex_bessel_j_complex::EXPORT,
    &analysis_special_complex_digamma_complex::EXPORT,
    &analysis_special_complex_gamma_complex::EXPORT,
    &analysis_special_complex_ln_gamma_complex::EXPORT,
    &analysis_special_complex_lambert_w_complex::EXPORT,
    &analysis_meijer_meijer_g::EXPORT,
    &analysis_meijer_meijer_g_real::EXPORT,
    &analysis_convergence_levin_condition_number::EXPORT,
    &analysis_convergence_levin_condition_numbers::EXPORT,
    &analysis_convergence_levin_transform::EXPORT,
    &analysis_convergence_levin_transform_terms::EXPORT,
    &analysis_convergence_levin_sum_lazy::EXPORT,
    &analysis_coordinates_cartesian_to_spherical::EXPORT,
    &analysis_coordinates_spherical_to_cartesian::EXPORT,
    &analysis_coordinates_from_cartesian::EXPORT,
    &analysis_coordinates_from_cartesian_pure::EXPORT,
    &analysis_coordinates_to_cartesian::EXPORT,
    &analysis_coordinates_to_cartesian_pure::EXPORT,
    &analysis_coordinates_transform_point::EXPORT,
    &analysis_coordinates_transform_point_pure::EXPORT,
    &analysis_coordinates_jacobian_determinant::EXPORT,
    &analysis_coordinates_numerical_jacobian::EXPORT,
    &analysis_wavelets_method_analysis_high::EXPORT,
    &analysis_wavelets_method_analysis_low::EXPORT,
    &analysis_wavelets_method_high_pass::EXPORT,
    &analysis_wavelets_method_low_pass::EXPORT,
    &analysis_wavelets_method_synthesis_high::EXPORT,
    &analysis_wavelets_method_synthesis_low::EXPORT,
    &analysis_wavelets_method_coeff_len::EXPORT,
    &analysis_wavelets_method_taps::EXPORT,
    &analysis_wavelets_method_is_orthogonal::EXPORT,
    &analysis_wavelets_method_vanishing_moments::EXPORT,
    &analysis_wavelets_dwt::EXPORT,
    &analysis_wavelets_dwt_with::EXPORT,
    &analysis_wavelets_idwt::EXPORT,
    &analysis_wavelets_idwt_with::EXPORT,
    &analysis_wavelets_dwt2::EXPORT,
    &analysis_wavelets_idwt2::EXPORT,
    &analysis_wavelets_extend::EXPORT,
    &analysis_wavelets_wavedec::EXPORT,
    &analysis_wavelets_wavedec_with::EXPORT,
    &analysis_wavelets_wavedec2::EXPORT,
    &analysis_wavelets_wavedec2_with::EXPORT,
    &analysis_wavelets_waverec::EXPORT,
    &analysis_wavelets_waverec_with::EXPORT,
    &analysis_wavelets_waverec2::EXPORT,
    &analysis_wavelets_waverec2_with::EXPORT,
    &analysis_rootfind_polynomial_roots_companion::EXPORT,
    &analysis_rootfind_polynomial_roots_jenkins_traub::EXPORT,
    &analysis_rootfind_homotopy_univariate::EXPORT,
    &analysis_rootfind_polynomial_roots_jenkins_traub_complex::EXPORT,
    &analysis_poly_eval_complex::EXPORT,
    &analysis_qp::EXPORT,
    &analysis_simplex::EXPORT,
    &analysis_simplex_with::EXPORT,
    &analysis_kkt::EXPORT,
    &analysis_homotopy_system::EXPORT,
    &analysis_homotopy_gamma::EXPORT,
];
