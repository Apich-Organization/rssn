//! Checked native handles are session-local capabilities with bounded lifetimes.
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::native_objects::NativeObjects;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::value::OpaqueValue;

const KIND: &str = "rssn/test-native/v1";
const MEMORY_LIMIT: usize = 64 * 1024 * 1024;
const OBJECT_LIMIT: usize = 4096;

fn read(
    objects: &NativeObjects,
    handle: &Value,
) -> Result<i64, OpError> {
    objects.with::<i64, _>(handle, KIND, |value| Ok(*value))
}

#[test]
fn handles_reject_wrong_session_kind_type_and_malformed_encoding() {
    let objects = NativeObjects::new();
    let foreign = NativeObjects::new();
    let handle = objects.insert(KIND, vec![], 17i64, 0).unwrap();
    // Both stores deliberately have object number one: scope authentication
    // must distinguish them before any native value can be borrowed.
    let foreign_handle = foreign.insert(KIND, vec![], 99i64, 0).unwrap();
    assert_eq!(read(&objects, &handle).unwrap(), 17);
    assert!(read(&objects, &foreign_handle).is_err());
    assert!(read(&foreign, &handle).is_err());
    assert!(objects.dimensions(&foreign_handle).is_err());
    assert!(objects.release(&foreign_handle).is_err());
    assert!(
        objects
            .with::<i64, _>(&handle, "rssn/wrong-kind/v1", |_| Ok(()))
            .is_err()
    );
    assert!(
        objects
            .with::<String, _>(&handle, KIND, |_| Ok(()))
            .is_err()
    );
    assert!(
        objects
            .with_mut::<String, _>(&handle, KIND, |_| Ok(()))
            .is_err()
    );
    assert!(read(&objects, &Value::int(1)).is_err());

    let Value::Opaque(original) = &handle else {
        panic!("opaque native handle")
    };
    let forged_kind = Value::Opaque(OpaqueValue {
        kind: "rssn/forged-kind/v1".into(),
        bytes: original.bytes.clone(),
    });
    assert!(read(&objects, &forged_kind).is_err());
    assert!(objects.dimensions(&forged_kind).is_err());
    for bytes in [vec![], vec![0], vec![0; 15], vec![0; 16], vec![0; 17]] {
        assert!(
            read(
                &objects,
                &Value::Opaque(OpaqueValue {
                    kind: KIND.into(),
                    bytes
                })
            )
            .is_err()
        );
    }
    let mut unknown_id = original.bytes.clone();
    unknown_id[8..].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(
        read(
            &objects,
            &Value::Opaque(OpaqueValue {
                kind: KIND.into(),
                bytes: unknown_id
            })
        )
        .is_err()
    );
    // Failed authentication and failed downcasts never poison the live value.
    assert_eq!(read(&objects, &handle).unwrap(), 17);
    assert_eq!(objects.usage().unwrap().0, 1);
}

#[test]
fn release_invalidates_all_copies_without_reusing_old_identifiers() {
    let objects = NativeObjects::new();
    let dimensions = vec![DimExpr::Base("Length".into())];
    let handle = objects
        .insert(KIND, dimensions.clone(), 5i64, 1024)
        .unwrap();
    let alias = handle.clone();
    assert_eq!(objects.dimensions(&alias).unwrap(), dimensions);
    assert_eq!(
        objects
            .clone()
            .with::<i64, _>(&alias, KIND, |n| Ok(*n))
            .unwrap(),
        5
    );
    objects.release(&handle).unwrap();
    assert!(read(&objects, &alias).is_err());
    assert!(objects.dimensions(&alias).is_err());
    assert!(objects.release(&alias).is_err());
    assert_eq!(objects.usage().unwrap(), (0, 0));
    let next = objects.insert(KIND, vec![], 9i64, 0).unwrap();
    assert_ne!(next, handle);
    assert!(read(&objects, &alias).is_err());
    assert_eq!(read(&objects, &next).unwrap(), 9);
}

#[test]
fn conflicting_reentrant_borrows_return_errors_and_restore_access() {
    let objects = NativeObjects::new();
    let handle = objects.insert(KIND, vec![], 4i64, 0).unwrap();
    objects
        .with_mut::<i64, _>(&handle, KIND, |value| {
            assert!(read(&objects, &handle).is_err());
            assert!(
                objects
                    .with_mut::<i64, _>(&handle, KIND, |_| Ok(()))
                    .is_err()
            );
            assert!(objects.release(&handle).is_err());
            *value += 1;
            Ok(())
        })
        .unwrap();
    objects
        .with::<i64, _>(&handle, KIND, |value| {
            assert_eq!(read(&objects, &handle).unwrap(), *value);
            assert!(
                objects
                    .with_mut::<i64, _>(&handle, KIND, |_| Ok(()))
                    .is_err()
            );
            Ok(())
        })
        .unwrap();
    let failed: Result<(), OpError> = objects.with_mut::<i64, _>(&handle, KIND, |_| {
        Err(OpError::Host("intentional callback failure".into()))
    });
    assert!(failed.is_err());
    assert_eq!(read(&objects, &handle).unwrap(), 5);
    objects.release(&handle).unwrap();
}

#[test]
fn object_and_memory_budgets_fail_atomically_and_release_restores_capacity() {
    let objects = NativeObjects::new();
    assert!(objects.insert(KIND, vec![], (), MEMORY_LIMIT + 1).is_err());
    assert!(objects.insert(KIND, vec![], (), usize::MAX).is_err());
    assert_eq!(objects.usage().unwrap(), (0, 0));
    let charged = objects.insert(KIND, vec![], (), MEMORY_LIMIT / 2).unwrap();
    let before = objects.usage().unwrap();
    assert!(before.1 >= MEMORY_LIMIT / 2 && before.1 <= MEMORY_LIMIT);
    assert!(
        objects
            .insert(KIND, vec![], (), MEMORY_LIMIT / 2 + 1)
            .is_err()
    );
    assert_eq!(objects.usage().unwrap(), before);
    objects.release(&charged).unwrap();
    assert_eq!(objects.usage().unwrap(), (0, 0));

    let mut handles = Vec::new();
    for _ in 0..OBJECT_LIMIT {
        handles.push(objects.insert(KIND, vec![], (), 0).unwrap());
    }
    assert_eq!(objects.usage().unwrap().0, OBJECT_LIMIT);
    let full = objects.usage().unwrap();
    assert!(objects.insert(KIND, vec![], (), 0).is_err());
    assert_eq!(objects.usage().unwrap(), full);
    objects.release(&handles.pop().unwrap()).unwrap();
    let replacement = objects.insert(KIND, vec![], (), 0).unwrap();
    assert_eq!(objects.usage().unwrap().0, OBJECT_LIMIT);
    objects.release(&replacement).unwrap();
    for handle in handles {
        objects.release(&handle).unwrap();
    }
    assert_eq!(objects.usage().unwrap(), (0, 0));
}

#[test]
fn dependent_handles_guard_parent_release_and_invalid_parent_inserts_are_atomic() {
    let objects = NativeObjects::new();
    let foreign = NativeObjects::new();
    let parent = objects.insert(KIND, vec![], 10i64, 1024).unwrap();
    let child = objects
        .insert_dependent(
            "rssn/test-dependent/v1",
            vec![],
            (),
            128,
            std::slice::from_ref(&parent),
        )
        .unwrap();
    let alias = parent.clone();
    assert!(objects.release(&alias).is_err());
    assert_eq!(read(&objects, &parent).unwrap(), 10);
    let foreign_parent = foreign.insert(KIND, vec![], 0i64, 0).unwrap();
    let before = objects.usage().unwrap();
    assert!(
        objects
            .insert_dependent(KIND, vec![], (), 0, &[foreign_parent])
            .is_err()
    );
    assert_eq!(objects.usage().unwrap(), before);
    objects.release(&child).unwrap();
    assert!(
        objects
            .insert_dependent(KIND, vec![], (), 0, &[child])
            .is_err()
    );
    objects.release(&parent).unwrap();
    assert!(read(&objects, &alias).is_err());
    assert_eq!(objects.usage().unwrap(), (0, 0));
}

#[test]
fn nested_dynamic_scopes_invalidate_escaped_handles_and_preserve_older_objects() {
    let objects = NativeObjects::new();
    let base = objects.insert(KIND, vec![], 1i64, 1000).unwrap();
    let baseline = objects.usage().unwrap();
    let escaped;
    {
        let outer = objects.scoped().unwrap();
        escaped = objects
            .insert_dependent(KIND, vec![], 2i64, 2000, std::slice::from_ref(&base))
            .unwrap();
        let inner_escaped;
        {
            let inner = objects.scoped().unwrap();
            inner_escaped = objects.insert(KIND, vec![], 3i64, 3000).unwrap();
            assert_eq!(read(&objects, &inner_escaped).unwrap(), 3);
            drop(inner);
        }
        assert!(read(&objects, &inner_escaped).is_err());
        assert_eq!(read(&objects, &escaped).unwrap(), 2);
        assert!(objects.release(&base).is_err());
        drop(outer);
    }
    assert!(read(&objects, &escaped).is_err());
    assert!(objects.dimensions(&escaped).is_err());
    assert_eq!(read(&objects, &base).unwrap(), 1);
    assert_eq!(objects.usage().unwrap(), baseline);
    objects.release(&base).unwrap();
    assert_eq!(objects.usage().unwrap(), (0, 0));
}

#[test]
fn native_metadata_rejects_excessive_or_nonconcrete_domains_without_retaining_objects() {
    use rssn::kernel::signature::Exp;
    use rssn::kernel::signature::Rational;
    let objects = NativeObjects::new();
    assert!(objects.insert(&"k".repeat(257), vec![], (), 0).is_err());
    assert!(
        objects
            .insert(KIND, vec![DimExpr::Dimless; 65], (), 0)
            .is_err()
    );
    assert!(
        objects
            .insert(KIND, vec![DimExpr::Base("x".repeat(129))], (), 0)
            .is_err()
    );
    assert!(
        objects
            .insert(KIND, vec![DimExpr::Var("D".into())], (), 0)
            .is_err()
    );
    let mut deep = DimExpr::Base("Length".into());
    for _ in 0..33 {
        deep = DimExpr::Inverse(Box::new(deep));
    }
    assert!(objects.insert(KIND, vec![deep], (), 0).is_err());
    assert!(
        objects
            .insert(KIND, vec![DimExpr::Mul(vec![DimExpr::Dimless; 257])], (), 0)
            .is_err()
    );
    let invalid_power = DimExpr::Pow(
        Box::new(DimExpr::Base("Length".into())),
        Exp::Rat(Rational { numer: 1, denom: 0 }),
    );
    assert!(objects.insert(KIND, vec![invalid_power], (), 0).is_err());
    assert_eq!(objects.usage().unwrap(), (0, 0));
    let parent = objects.insert(KIND, vec![], (), 0).unwrap();
    let before = objects.usage().unwrap();
    assert!(
        objects
            .insert_dependent(KIND, vec![], (), 0, &vec![parent.clone(); 65])
            .is_err()
    );
    assert_eq!(objects.usage().unwrap(), before);
    objects.release(&parent).unwrap();
    assert_eq!(objects.usage().unwrap(), (0, 0));
}

#[test]
fn committed_scopes_keep_new_objects_while_failed_scopes_rollback() {
    let objects = NativeObjects::new();
    let existing = objects.insert(KIND, vec![], 17i64, 0).unwrap();
    let scope = objects.scoped().unwrap();
    let retained = objects.insert(KIND, vec![], 23i64, 0).unwrap();
    scope.commit();
    assert_eq!(read(&objects, &retained).unwrap(), 23);
    let failed = objects.scoped().unwrap();
    let rolled_back = objects.insert(KIND, vec![], 31i64, 0).unwrap();
    drop(failed);
    assert!(read(&objects, &rolled_back).is_err());
    assert_eq!(read(&objects, &existing).unwrap(), 17);
    assert_eq!(read(&objects, &retained).unwrap(), 23);
}
