use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
pub fn factors(expr: &DimExpr) -> std::collections::BTreeMap<String, num_rational::BigRational> {
    fn collect(
        expr: &DimExpr,
        power: num_rational::BigRational,
        out: &mut std::collections::BTreeMap<String, num_rational::BigRational>,
    ) {
        match expr {
            | DimExpr::Dimless => {},
            | DimExpr::Base(name) => *out.entry(name.clone()).or_default() += power,
            | DimExpr::Mul(parts) => {
                for part in parts {
                    collect(part, power.clone(), out);
                }
            },
            | DimExpr::Inverse(inner) => collect(inner, -power, out),
            | DimExpr::Pow(inner, Exp::Rat(r)) => {
                collect(
                    inner,
                    power * num_rational::BigRational::new(r.numer.into(), r.denom.into()),
                    out,
                )
            },
            | other => panic!("nonconcrete result dimension {other:?}"),
        }
    }
    let mut result = std::collections::BTreeMap::new();
    collect(
        expr,
        num_rational::BigRational::from_integer(1.into()),
        &mut result,
    );
    result.retain(|_, power| power != &num_rational::BigRational::from_integer(0.into()));
    result
}
