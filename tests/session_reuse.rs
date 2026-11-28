//! Standard sessions share immutable setup while their units can diverge.

use std::sync::Arc;

use rssn::api::Session;
use rssn::dim::Scale;
use rssn::dim::UnitDef;
use rssn::dim::units;
use rssn::kernel::signature::DimExpr;

#[test]
fn standard_templates_refresh_for_global_units_and_keep_forks_local() {
    let session = Session::new();
    let initial = session.manifest();
    let threaded = std::thread::spawn(|| Session::new().manifest())
        .join()
        .unwrap();
    assert!(Arc::ptr_eq(&initial, &threaded));

    let local_name = "rssn_test_fork_local_unit";
    let mut fork = session.fork();
    assert_ne!(session.id(), fork.id());
    assert!(Arc::ptr_eq(&initial, &fork.manifest()));
    fork.register_unit(UnitDef::linear(local_name, DimExpr::Dimless, Scale::one()))
        .unwrap();
    assert!(
        fork.manifest()
            .units
            .iter()
            .any(|unit| unit.name == local_name)
    );
    assert!(!initial.units.iter().any(|unit| unit.name == local_name));
    assert!(Arc::ptr_eq(&initial, &Session::new().manifest()));

    let global_name = "rssn_test_process_wide_unit";
    units::register(UnitDef::linear(global_name, DimExpr::Dimless, Scale::one())).unwrap();
    let after_registration = Session::new().manifest();
    assert!(
        after_registration
            .units
            .iter()
            .any(|unit| unit.name == global_name)
    );
    assert!(
        !after_registration
            .units
            .iter()
            .any(|unit| unit.name == local_name)
    );
    assert!(
        !session
            .manifest()
            .units
            .iter()
            .any(|unit| unit.name == global_name)
    );
    assert!(Arc::ptr_eq(&after_registration, &Session::new().manifest()));
}
