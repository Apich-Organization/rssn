//! Explicit pinned-call diagnostics and compatibility.

#[test]
fn live_key_diagnostics_do_not_dispatch_to_a_sibling_overload() {
    use rssn::api::Session;
    use rssn::dim::Scale;
    use rssn::dim::UnitDef;
    use rssn::kernel::call::CallCtx;
    use rssn::kernel::call::Kernel;
    use rssn::kernel::call::OpError;
    use rssn::kernel::signature::ArityInfo;
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::signature::TypeExpr;
    use rssn::kernel::value::Value;

    let session = Session::new();
    let (key, _) = session
        .overloads("quantity")
        .into_iter()
        .find(|(_, sig)| sig.params.get(1) == Some(&TypeExpr::Symbol))
        .expect("quantity's Symbol overload");
    let args = [
        Value::Float(3.0),
        Value::Unit(UnitDef::linear("m", DimExpr::base("Length"), Scale::one()).to_value()),
    ];
    let ctx = CallCtx::new();
    assert!(
        session.call_value("quantity", &args, &ctx).is_ok(),
        "the other overload still accepts a unit"
    );
    assert_eq!(
        session.call_key(&key, &args, &ctx),
        Err(OpError::UnknownOp(key.clone()))
    );
    assert_eq!(
        Kernel::call(&session, &key, &args, &ctx),
        Err(OpError::UnknownOp(key.clone()))
    );
    assert!(matches!(
        session.call_key_diagnosed(&key, &args, &ctx),
        Err(OpError::Type { arg: 1, expected, got })
            if *expected == TypeExpr::Symbol && matches!(*got, TypeExpr::Unit(_))
    ));
    assert_eq!(
        session.call_key_diagnosed(&key, &args[..1], &ctx),
        Err(OpError::Arity {
            expected: ArityInfo::fixed(2),
            got: 1,
        })
    );
    let mut stale = key;
    stale.signature_hash[0] ^= 1;
    assert_eq!(
        session.call_key_diagnosed(&stale, &args, &ctx),
        Err(OpError::UnknownOp(stale))
    );
}
