use std::sync::Arc;

use rssn::is_exclusive;

#[test]

fn test_is_exclusive() {

    let arc = Arc::new(5);

    assert!(is_exclusive(&arc));

    let arc2 = arc.clone();

    assert!(!is_exclusive(&arc));

    assert!(!is_exclusive(&arc2));

    drop(arc2);

    assert!(is_exclusive(&arc));
}
