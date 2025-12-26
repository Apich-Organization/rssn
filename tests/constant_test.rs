use rssn::constant;

#[test]

fn test_get_build_date() {

    let date =
        constant::get_build_date();

    assert!(!date.is_empty());

    println!(
        "Build Date: {}",
        date
    );
}

#[test]

fn test_get_commit_sha() {

    let sha =
        constant::get_commit_sha();

    assert!(!sha.is_empty());

    println!(
        "Commit SHA: {}",
        sha
    );
}

#[test]

fn test_get_rustc_version() {

    let version =
        constant::get_rustc_version();

    assert!(!version.is_empty());

    println!(
        "Rustc Version: {}",
        version
    );
}

#[test]

fn test_get_cargo_target_triple() {

    let triple = constant::get_cargo_target_triple();

    assert!(!triple.is_empty());

    println!(
        "Target Triple: {}",
        triple
    );
}

#[test]

fn test_get_system_info() {

    let info =
        constant::get_system_info();

    assert!(!info.is_empty());

    println!(
        "System Info: {}",
        info
    );
}
