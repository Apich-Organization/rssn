//! The graph and the manifest share the cumulative probability dimension rule.

use rssn::graph::DimensionMode;
use rssn::kernel::signature::DimRuleInfo;
use rssn::prelude::ComputeError;
use rssn::prelude::Config;
use rssn::prelude::Session;

#[test]
fn cdf_manifest_preserves_the_dimension_check_in_its_portable_rule() {
    let session = Session::new();
    let manifest = session.manifest();
    let cdf = manifest.op("cdf").unwrap();
    let rule = cdf.dim_rule.as_ref().unwrap();
    assert_eq!(*rule, DimRuleInfo::CompatibleDimensionless);
    assert_eq!(DimRuleInfo::from(rule.to_rule()), *rule);
    let wire = serde_json::to_string(rule).unwrap();
    assert_eq!(wire, "\"compatible_dimensionless\"");
    assert_eq!(serde_json::from_str::<DimRuleInfo>(&wire).unwrap(), *rule);
}

#[test]
fn cumulative_probability_uses_the_point_dimension_but_returns_one() {
    let session = Session::new();
    let config = Config::new().dimensions(DimensionMode::Strict);
    let distribution = "normal(quantity(1, m), quantity(2, m))";
    for point in ["quantity(3, m)", "x"] {
        let term = session
            .parse(&format!("cdf({distribution}, {point})"))
            .unwrap();
        let result = session.compute(term, &config).unwrap();
        assert_eq!(result.dimension.unwrap().to_string(), "1");
    }
    for operator in ["cdf", "pdf"] {
        let term = session
            .parse(&format!("{operator}({distribution}, quantity(3, s))"))
            .unwrap();
        assert!(matches!(
            session.compute(term, &config),
            Err(ComputeError::Dimension(_))
        ));
    }
}
