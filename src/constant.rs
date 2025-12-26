pub const BUILD_DATE : &str =
    env!("VERGEN_BUILD_DATE");

pub const COMMIT_SHA : &str =
    env!("VERGEN_GIT_SHA");

pub const RUSTC_VERSION : &str =
    env!("VERGEN_RUSTC_SEMVER");

pub const CARGO_TARGET_TRIPLE : &str =
    env!("VERGEN_CARGO_TARGET_TRIPLE");

pub const SYSTEM_INFO : &str =
    env!("VERGEN_SYSINFO_OS_VERSION");

// --- Getter functions ---

/// Returns the build date (e.g., "2023-10-24").
///
/// # Examples
///
/// ```
/// 
/// use rssn::constant::get_build_date;
///
/// let date = get_build_date();
///
/// assert!(!date.is_empty());
/// ```
#[must_use]

pub const fn get_build_date(
) -> &'static str {

    BUILD_DATE
}

/// Returns the Git short SHA for the current commit.
///
/// # Examples
///
/// ```
/// 
/// use rssn::constant::get_commit_sha;
///
/// let sha = get_commit_sha();
///
/// assert!(!sha.is_empty());
/// ```
#[must_use]

pub const fn get_commit_sha(
) -> &'static str {

    COMMIT_SHA
}

/// Returns the rustc semantic version (e.g., "1.70.0-nightly").
///
/// # Examples
///
/// ```
/// 
/// use rssn::constant::get_rustc_version;
///
/// let version = get_rustc_version();
///
/// assert!(!version.is_empty());
/// ```
#[must_use]

pub const fn get_rustc_version(
) -> &'static str {

    RUSTC_VERSION
}

/// Returns the Cargo target triple (e.g., "x86_64-unknown-linux-gnu").
///
/// # Examples
///
/// ```
/// 
/// use rssn::constant::get_cargo_target_triple;
///
/// let triple = get_cargo_target_triple();
///
/// assert!(!triple.is_empty());
/// ```
#[must_use]

pub const fn get_cargo_target_triple(
) -> &'static str {

    CARGO_TARGET_TRIPLE
}

/// Returns the system information string (e.g., "Linux Arch Linux").
///
/// # Examples
///
/// ```
/// 
/// use rssn::constant::get_system_info;
///
/// let sys_info = get_system_info();
///
/// assert!(!sys_info.is_empty());
/// ```
#[must_use]

pub const fn get_system_info(
) -> &'static str {

    SYSTEM_INFO
}
