/// The date the library was built.

pub const BUILD_DATE: &str =
    env!("VERGEN_BUILD_DATE");

/// The Git commit SHA the library was built from.

pub const COMMIT_SHA: &str =
    env!("VERGEN_GIT_SHA");

/// The version of the Rust compiler used to build the library.

pub const RUSTC_VERSION: &str =
    env!("VERGEN_RUSTC_SEMVER");

/// The target triple for which the library was built.

pub const CARGO_TARGET_TRIPLE: &str =
    env!("VERGEN_CARGO_TARGET_TRIPLE");

/// Operating system and version information of the build environment.

pub const SYSTEM_INFO: &str =
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
#[inline(always)]

pub const fn get_build_date()
-> &'static str {

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
#[inline(always)]

pub const fn get_commit_sha()
-> &'static str {

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
#[inline(always)]

pub const fn get_rustc_version()
-> &'static str {

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
#[inline(always)]

pub const fn get_cargo_target_triple()
-> &'static str {

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
#[inline(always)]

pub const fn get_system_info()
-> &'static str {

    SYSTEM_INFO
}

// --- Math & Physics Constants ---

macro_rules! nist_const {
    ($const_name:ident, $fn_name:ident, $value:expr_2021, $unit:expr_2021, $uncert:expr_2021, $desc:expr_2021) => {
        #[doc = concat!($desc, "\n\n**Value:** ", stringify!($value), " ", $unit, "\n**Uncertainty:** ", $uncert)]
        pub const $const_name: f64 = $value;

        #[doc = concat!("Returns the ", $desc)]
        #[must_use]
        #[inline(always)]
        pub const fn $fn_name() -> f64 {
            $const_name
        }
    };
}

// --- Fundamental Constants ---

nist_const!(
    SPEED_OF_LIGHT,
    get_speed_of_light,
    299_792_458.0,
    "m s⁻¹",
    "exact",
    "speed of light in vacuum"
); // [cite: 16, 18, 19]

nist_const!(
    PLANCK_CONSTANT,
    get_planck_constant,
    6.626_070_15e-34,
    "J Hz⁻¹",
    "exact",
    "Planck constant"
); // [cite: 28, 30, 31]

nist_const!(
    ELEMENTARY_CHARGE,
    get_elementary_charge,
    1.602_176_634e-19,
    "C",
    "exact",
    "elementary charge"
); // [cite: 57, 59, 60]

nist_const!(
    BOLTZMANN_CONSTANT,
    get_boltzmann_constant,
    1.380_649e-23,
    "J K⁻¹",
    "exact",
    "Boltzmann constant"
); // [cite: 227, 229, 230]

nist_const!(
    AVOGADRO_CONSTANT,
    get_avogadro_constant,
    6.022_140_76e23,
    "mol⁻¹",
    "exact",
    "Avogadro constant"
); // [cite: 235, 237, 238]

// --- Measured Constants (with uncertainty) ---

nist_const!(
    GRAVITATIONAL_CONSTANT,
    get_gravitational_constant,
    6.674_30e-11,
    "m³ kg⁻¹ s⁻²",
    "0.00015 x 10⁻¹¹",
    "Newtonian constant of gravitation"
); // [cite: 24, 26, 27]

nist_const!(
    ELECTRON_MASS,
    get_electron_mass,
    9.109_383_713_9e-31,
    "kg",
    "2.8e-40",
    "electron mass"
); // [cite: 147, 149, 150]

nist_const!(
    FINE_STRUCTURE_CONSTANT,
    get_fine_structure_constant,
    7.297_352_564_3e-3,
    "dimensionless",
    "0.000_000_001_1e-3",
    "fine-structure constant"
); // [cite: 121, 123]

nist_const!(
    RYDBERG_CONSTANT,
    get_rydberg_constant,
    10_973_731.568_157,
    "m⁻¹",
    "0.000_000_012",
    "Rydberg constant"
); // [cite: 134, 136, 137]

nist_const!(
    VACUUM_ELECTRIC_PERMITTIVITY,
    get_vacuum_electric_permittivity,
    8.854_187_818_8e-12,
    "F m⁻¹",
    "0.000_000_001_4e-12",
    "vacuum electric permittivity"
); // [cite: 74, 76]

// --- Electromagnetic Constants ---

nist_const!(
    VACUUM_MAGNETIC_PERMEABILITY,
    get_vacuum_magnetic_permeability,
    1.256_637_061_27e-6,
    "N A⁻²",
    "0.000_000_000_20e-6",
    "vacuum magnetic permeability"
); // [cite: 65, 66, 67]

nist_const!(
    JOSEPHSON_CONSTANT,
    get_josephson_constant,
    483_597.848_4e9,
    "Hz V⁻¹",
    "exact (definition based)",
    "Josephson constant (2e/h)"
); // [cite: 82, 84, 85]

nist_const!(
    VON_KLITZING_CONSTANT,
    get_von_klitzing_constant,
    25_812.807_45,
    "Ω",
    "exact (definition based)",
    "von Klitzing constant (h/e²)"
); // [cite: 90, 92, 93]

nist_const!(
    MAGNETIC_FLUX_QUANTUM,
    get_magnetic_flux_quantum,
    2.067_833_848e-15,
    "Wb",
    "exact (definition based)",
    "magnetic flux quantum (h/2e)"
); // [cite: 96, 98, 100]

// --- Atomic & Particle Masses ---

nist_const!(
    PROTON_MASS_KG,
    get_proton_mass_kg,
    1.672_621_925_95e-27,
    "kg",
    "0.000_000_000_52e-27",
    "proton mass"
); // [cite: 35, 36]

nist_const!(
    NEUTRON_MASS_U,
    get_neutron_mass_u,
    1.008_664_916_06,
    "u",
    "0.000_000_000_40",
    "neutron mass in atomic mass units"
); // [cite: 174, 182, 183]

nist_const!(
    ATOMIC_MASS_CONSTANT,
    get_atomic_mass_constant,
    1.660_539_068_92e-27,
    "kg",
    "0.000_000_000_52e-27",
    "atomic mass constant (m_u)"
); // [cite: 243, 244, 245]

nist_const!(
    PROTON_ELECTRON_MASS_RATIO,
    get_proton_electron_mass_ratio,
    1_836.152_673_426,
    "dimensionless",
    "0.000_000_032",
    "proton-electron mass ratio"
); // [cite: 54, 56]

// --- Magnetic Moments & Factors ---

nist_const!(
    BOHR_MAGNETON,
    get_bohr_magneton,
    9.274_010_065_7e-24,
    "J T⁻¹",
    "0.000_000_002_9e-24",
    "Bohr magneton"
); // [cite: 101, 103, 107]

nist_const!(
    NUCLEAR_MAGNETON,
    get_nuclear_magneton,
    5.050_783_739_3e-27,
    "J T⁻¹",
    "0.000_000_001_6e-27",
    "nuclear magneton"
); // [cite: 114, 116, 117]

nist_const!(
    ELECTRON_G_FACTOR,
    get_electron_g_factor,
    -2.002_319_304_360_92,
    "dimensionless",
    "0.000_000_000_000_36",
    "electron g-factor"
); // [cite: 259, 265, 266]

// --- Physico-Chemical Constants ---

nist_const!(
    MOLAR_GAS_CONSTANT,
    get_molar_gas_constant,
    8.314_462_618,
    "J mol⁻¹ K⁻¹",
    "exact (defined by R = Na * k)",
    "molar gas constant"
); // [cite: 261, 262, 263, 264]

nist_const!(
    FARADAY_CONSTANT,
    get_faraday_constant,
    96_485.332_12,
    "C mol⁻¹",
    "exact (defined by F = Na * e)",
    "Faraday constant"
); // [cite: 255, 256, 257, 258]

nist_const!(
    STEFAN_BOLTZMANN_CONSTANT,
    get_stefan_boltzmann_constant,
    5.670_374_419e-8,
    "W m⁻² K⁻⁴",
    "exact (calculated from k, h, c)",
    "Stefan-Boltzmann constant"
); // [cite: 312, 313, 317]

// --- Additional Fundamental & Quantum Constants ---

nist_const!(
    REDUCED_PLANCK_CONSTANT,
    get_reduced_planck_constant,
    1.054_571_817e-34,
    "J s",
    "uncertainty in source",
    "reduced Planck constant (h-bar)"
); //

nist_const!(
    INVERSE_FINE_STRUCTURE_CONSTANT,
    get_inverse_fine_structure_constant,
    137.035_999_177,
    "dimensionless",
    "0.000_000_021",
    "inverse fine-structure constant \
     (1/α)"
); //

nist_const!(
    BOHR_RADIUS,
    get_bohr_radius,
    5.291_772_105_44e-11,
    "m",
    "0.000_000_000_82e-11",
    "Bohr radius (a₀)"
); //

nist_const!(
    HARTREE_ENERGY,
    get_hartree_energy,
    4.359_744_722_206_0e-18,
    "J",
    "0.000_000_000_004_8e-18",
    "Hartree energy (Eh)"
); //

// --- Particle Masses (Atomic Units) ---

nist_const!(
    ELECTRON_MASS_U,
    get_electron_mass_u,
    5.485_799_090_441e-4,
    "u",
    "0.000_000_000_097e-4",
    "electron mass in atomic mass \
     units"
); //

nist_const!(
    PROTON_MASS_U,
    get_proton_mass_u,
    1.007_276_466_578_9,
    "u",
    "0.000_000_000_008_3",
    "proton mass in atomic mass units"
); //

nist_const!(
    DEUTERON_MASS_U,
    get_deuteron_mass_u,
    2.013_553_212_544,
    "u",
    "0.000_000_000_015",
    "deuteron mass in atomic mass \
     units"
); //

nist_const!(
    ALPHA_PARTICLE_MASS_U,
    get_alpha_particle_mass_u,
    4.001_506_179_129,
    "u",
    "0.000_000_000_062",
    "alpha particle mass in atomic \
     mass units"
); //

// --- Electromagnetic Interaction ---

nist_const!(
    CLASSICAL_ELECTRON_RADIUS,
    get_classical_electron_radius,
    2.817_940_320_5e-15,
    "m",
    "0.000_000_001_3e-15",
    "classical electron radius"
); //

nist_const!(
    THOMSON_CROSS_SECTION,
    get_thomson_cross_section,
    6.652_458_705_1e-29,
    "m²",
    "0.000_000_006_2e-29",
    "Thomson cross section"
); //

// --- Radiation Constants ---

nist_const!(
    WIEN_DISPLACEMENT_CONSTANT,
    get_wien_displacement_constant,
    2.897_771_955e-3,
    "m K",
    "exact (calculated)",
    "Wien displacement law constant \
     (b)"
); //

nist_const!(
    FIRST_RADIATION_CONSTANT,
    get_first_radiation_constant,
    3.741_771_852e-16,
    "W m²",
    "exact (calculated)",
    "first radiation constant (c₁)"
); //

nist_const!(
    SECOND_RADIATION_CONSTANT,
    get_second_radiation_constant,
    1.438_776_877e-2,
    "m K",
    "exact (calculated)",
    "second radiation constant (c₂)"
); //

// --- Muon Data ---
nist_const!(
    MUON_G_FACTOR,
    get_muon_g_factor,
    -2.002_331_841_23,
    "dimensionless",
    "0.000_000_000_82",
    "muon g-factor"
); // [cite: 20-23]

nist_const!(
    MUON_MASS_U,
    get_muon_mass_u,
    0.113_428_925_7,
    "u",
    "0.000_000_002_5",
    "muon mass in atomic mass units"
); // [cite: 272-273]

nist_const!(
    MUON_ELECTRON_MASS_RATIO,
    get_muon_electron_mass_ratio,
    206.768_282_7,
    "dimensionless",
    "0.000_004_6",
    "muon-electron mass ratio"
); // [cite: 278]

// --- Magnetic Moments & Shielding ---
nist_const!(
    PROTON_MAGNETIC_MOMENT,
    get_proton_magnetic_moment,
    1.410_606_795_45e-26,
    "J T⁻¹",
    "0.000_000_000_60e-26",
    "proton magnetic moment"
); // [cite: 61-63]

nist_const!(
    NEUTRON_MAGNETIC_MOMENT,
    get_neutron_magnetic_moment,
    -9.662_365_3e-27,
    "J T⁻¹",
    "0.000_002_3e-27",
    "neutron magnetic moment"
); // [cite: 188-189]

nist_const!(
    PROTON_MAGNETIC_SHIELDING_CORRECTION, get_proton_magnetic_shielding_correction,
    2.567_15e-5, "dimensionless", "0.000_41e-5",
    "proton magnetic shielding correction (H2O sphere, 25°C)"
); // [cite: 77-78, 81]

nist_const!(
    SHIELDED_PROTON_GYROMAGNETIC_RATIO, get_shielded_proton_gyromagnetic_ratio,
    2.675_153_194e8, "s⁻¹ T⁻¹", "0.000_000_011e8",
    "shielded proton gyromagnetic ratio (H2O, sphere, 25°C)"
); // [cite: 104-106, 108]

// --- Mass Ratios & Specific Quotients ---
nist_const!(
    NEUTRON_PROTON_MASS_RATIO,
    get_neutron_proton_mass_ratio,
    1.001_378_419_46,
    "dimensionless",
    "0.000_000_000_40",
    "neutron-proton mass ratio"
); // [cite: 186-187]

nist_const!(
    ELECTRON_MUON_MASS_RATIO,
    get_electron_muon_mass_ratio,
    4.836_331_70e-3,
    "dimensionless",
    "0.000_000_11e-3",
    "electron-muon mass ratio"
); // [cite: 158-160]

nist_const!(
    DEUTERON_PROTON_MASS_RATIO,
    get_deuteron_proton_mass_ratio,
    1.999_007_501_269_9,
    "dimensionless",
    "0.000_000_000_008_4",
    "deuteron-proton mass ratio"
); // [cite: 200-201]

nist_const!(
    ELECTRON_CHARGE_TO_MASS_QUOTIENT, get_electron_charge_to_mass_quotient,
    -1.758_820_008_38e11, "C kg⁻¹", "0.000_000_000_55e11",
    "electron charge to mass quotient"
); // [cite: 163-164, 170-171]

// --- Physicochemical Data ---
nist_const!(
    MOLAR_VOLUME_IDEAL_GAS,
    get_molar_volume_ideal_gas,
    22.413_969_54e-3,
    "m³ mol⁻¹",
    "exact (at 273.15 K, 101.325 kPa)",
    "molar volume of ideal gas"
); // [cite: 308-311]

nist_const!(
    MUON_MAGNETIC_MOMENT,
    get_muon_magnetic_moment,
    -4.490_448_30e-26,
    "J T⁻¹",
    "0.000_000_10e-26",
    "muon magnetic moment"
);
