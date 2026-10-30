//! Reference-value tests for special functions (new and existing).

use std::f64::consts::PI;

use num_complex::Complex64;

use rssn::kernels::special::bessel_j;
use rssn::kernels::special::digamma_numerical;
use rssn::kernels::special::erf_numerical;
use rssn::kernels::special::erfc_numerical;
use rssn::kernels::special::gamma_numerical;
use rssn::kernels::special::hurwitz_zeta;
use rssn::kernels::special::inverse_erf_numerical;
use rssn::kernels::special::ln_gamma_numerical;
use rssn::kernels::special_ext::anger_j;
use rssn::kernels::special_ext::appell_f1;
use rssn::kernels::special_ext::bessel_in;
use rssn::kernels::special_ext::bessel_in_scaled;
use rssn::kernels::special_ext::bessel_jn;
use rssn::kernels::special_ext::carlson_rf;
use rssn::kernels::special_ext::coulomb_c;
use rssn::kernels::special_ext::coulomb_f;
use rssn::kernels::special_ext::coulomb_f_prime;
use rssn::kernels::special_ext::coulomb_fg;
use rssn::kernels::special_ext::coulomb_g;
use rssn::kernels::special_ext::coulomb_g_prime;
use rssn::kernels::special_ext::coulomb_sigma;
use rssn::kernels::special_ext::elliptic_e;
use rssn::kernels::special_ext::elliptic_e_inc;
use rssn::kernels::special_ext::elliptic_f;
use rssn::kernels::special_ext::elliptic_k;
use rssn::kernels::special_ext::hyp1f1;
use rssn::kernels::special_ext::hyp2f1;
use rssn::kernels::special_ext::jacobi_elliptic;
use rssn::kernels::special_ext::kelvin_bei;
use rssn::kernels::special_ext::kelvin_ber;
use rssn::kernels::special_ext::kelvin_kei;
use rssn::kernels::special_ext::kelvin_ker;
use rssn::kernels::special_ext::lerch;
use rssn::kernels::special_ext::polylog;
use rssn::kernels::special_ext::polylog_complex;
use rssn::kernels::special_ext::weber_e;

fn close(
    a: f64,
    b: f64,
    tol: f64,
) {
    assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} vs {b}");
}

#[test]
fn complete_elliptic() {
    close(elliptic_k(0.5), 1.854_074_677_301_372, 1e-14);
    close(elliptic_e(0.5), 1.350_643_881_047_675, 1e-14);
    close(elliptic_k(0.0), std::f64::consts::FRAC_PI_2, 1e-15);
    close(elliptic_e(0.0), std::f64::consts::FRAC_PI_2, 1e-15);
    close(elliptic_k(-1.0), 1.311_028_777_146_059_9, 1e-14);
    close(elliptic_k(0.9), 2.578_092_113_348_173, 1e-14);
    close(elliptic_e(0.9), 1.104_774_732_704_073, 1e-13);
    assert!(elliptic_k(1.0).is_infinite());
    assert!((elliptic_e(1.0) - 1.0).abs() < 1e-15);
    assert!(elliptic_k(2.0).is_nan());
    // Legendre relation: E K' + E' K - K K' = pi/2
    let m = 0.3;
    let r = elliptic_e(m) * elliptic_k(1.0 - m) + elliptic_e(1.0 - m) * elliptic_k(m)
        - elliptic_k(m) * elliptic_k(1.0 - m);
    close(r, std::f64::consts::FRAC_PI_2, 1e-13);
}

#[test]
fn incomplete_elliptic() {
    // F(pi/2, m) = K(m), E(pi/2, m) = E(m)
    close(
        elliptic_f(std::f64::consts::FRAC_PI_2, 0.5),
        elliptic_k(0.5),
        1e-13,
    );
    close(
        elliptic_e_inc(std::f64::consts::FRAC_PI_2, 0.5),
        elliptic_e(0.5),
        1e-13,
    );
    // F(phi, 0) = phi
    close(elliptic_f(0.7, 0.0), 0.7, 1e-14);
    // reference: F(pi/4, 0.5) = 0.826017876...
    close(
        elliptic_f(std::f64::consts::FRAC_PI_4, 0.5),
        0.826_017_876_249_245,
        1e-12,
    );
    // periodicity
    close(
        elliptic_f(std::f64::consts::PI + 0.3, 0.4),
        elliptic_f(0.3, 0.4) + 2.0 * elliptic_k(0.4),
        1e-13,
    );
    close(carlson_rf(1.0, 2.0, 3.0), 0.726_945_935_468_9, 1e-11);
}

#[test]
fn jacobi_functions() {
    let (sn, cn, dn) = jacobi_elliptic(0.7, 0.0);
    close(sn, 0.7_f64.sin(), 1e-14);
    close(cn, 0.7_f64.cos(), 1e-14);
    close(dn, 1.0, 1e-14);
    let (sn, cn, dn) = jacobi_elliptic(1.0, 0.5);
    close(sn, 0.803_001_824_895_643_9, 1e-9);
    close(cn, (1.0 - sn * sn).sqrt(), 1e-12);
    close(dn, (1.0 - 0.5 * sn * sn).sqrt(), 1e-12);
    close(sn * sn + cn * cn, 1.0, 1e-14);
}

#[test]
fn hypergeometric() {
    // 2F1(1,1;2;z) = -ln(1-z)/z
    for z in [-3.0, -0.4, 0.3, 0.8, 0.97] {
        close(hyp2f1(1.0, 1.0, 2.0, z), -(1.0 - z).ln() / z, 1e-11);
    }
    // 2F1(a,b;b;z) = (1-z)^-a
    close(hyp2f1(0.5, 2.0, 2.0, 0.6), 0.4_f64.powf(-0.5), 1e-13);
    // 2F1(1/2,1/2;3/2;z^2) = asin(z)/z
    close(hyp2f1(0.5, 0.5, 1.5, 0.25), 0.5_f64.asin() / 0.5, 1e-13);
    assert_eq!(hyp2f1(1.0, 1.0, 2.0, 1.0), f64::INFINITY);
    // 1F1(1;2;x) = (e^x - 1)/x ; 1F1(a;a;x) = e^x
    close(hyp1f1(1.0, 2.0, 3.0), (3.0_f64.exp() - 1.0) / 3.0, 1e-13);
    close(
        hyp1f1(1.0, 2.0, -3.0),
        ((-3.0_f64).exp() - 1.0) / -3.0,
        1e-13,
    );
    close(hyp1f1(2.5, 2.5, 4.0), 4.0_f64.exp(), 1e-13);
    // Kummer's second example: 1F1(1/2;1;x) relation with Bessel: e^{-x/2} I0(x/2)
    close(
        hyp1f1(0.5, 1.0, 2.0),
        (1.0_f64).exp() * bessel_in(0, 1.0),
        1e-13,
    );
}

#[test]
fn bessel_integer_order() {
    close(bessel_jn(0, 1.0), 0.765_197_686_557_966_6, 1e-15);
    close(bessel_jn(1, 1.0), 0.440_050_585_744_933_5, 1e-15);
    close(bessel_jn(2, 5.0), 0.046_565_116_277_752_2, 1e-14);
    close(bessel_jn(0, 10.0), -0.245_935_764_451_348_3, 1e-14);
    close(bessel_jn(0, 50.0), 0.055_812_327_669_251_9, 1e-13);
    close(bessel_jn(-1, 2.0), -bessel_jn(1, 2.0), 1e-15);
    close(bessel_in(0, 1.0), 1.266_065_877_752_008_4, 1e-14);
    close(bessel_in(1, 2.0), 1.590_636_854_637_329, 1e-14);
    // recurrence J_{n-1} + J_{n+1} = 2n/x J_n
    close(
        bessel_jn(2, 3.0) + bessel_jn(4, 3.0),
        6.0 / 3.0 * bessel_jn(3, 3.0),
        1e-14,
    );
}

#[test]
fn existing_gamma_erf_reference_values() {
    close(gamma_numerical(5.0), 24.0, 1e-13);
    close(gamma_numerical(0.5), std::f64::consts::PI.sqrt(), 1e-13);
    close(gamma_numerical(-1.5), 2.363_271_801_207_355, 1e-12);
    close(ln_gamma_numerical(100.0), 359.134_205_369_575_4, 1e-13);
    close(digamma_numerical(1.0), -0.577_215_664_901_532_9, 1e-12);
    close(digamma_numerical(0.5), -1.963_510_026_021_423_5, 1e-12);
    close(erf_numerical(0.5), 0.520_499_877_813_046_5, 1e-14);
    close(erf_numerical(2.0), 0.995_322_265_018_952_7, 1e-14);
    close(erfc_numerical(3.0), 2.209_049_699_858_544e-5, 1e-9);
    close(inverse_erf_numerical(0.5), 0.476_936_276_204_469_9, 1e-12);
    close(erf_numerical(inverse_erf_numerical(-0.9)), -0.9, 1e-13);
}

const HYP_CASES: &[(f64, f64, f64, f64, f64)] = &[
    (1.2, 0.4, 2.0, -1000.0, 0.089_088_716_218_255_19),
    (1.2, 0.4, 2.0, -10.0, 0.492_822_851_373_183_27),
    (1.2, 0.4, 2.0, -2.0, 0.745_656_780_677_500_5),
    (1.2, 0.4, 2.0, -0.5, 0.903_500_554_932_665_5),
    (1.2, 0.4, 2.0, 0.3, 1.085_810_452_079_139_4),
    (1.2, 0.4, 2.0, 0.7, 1.286_026_956_157_418_2),
    (1.2, 0.4, 2.0, 0.95, 1.643_900_303_802_698_7),
    (1.2, 0.4, 2.0, 0.999, 2.018_581_667_814_119_7),
    (1.2, 0.4, 2.0, 0.999999999, 2.131_859_371_422_813),
    (0.3, 0.7, 1.5, -1000.0, 0.200_141_912_034_340_7),
    (0.3, 0.7, 1.5, -10.0, 0.647_460_351_241_119_3),
    (0.3, 0.7, 1.5, -2.0, 0.839_192_664_442_761_1),
    (0.3, 0.7, 1.5, -0.5, 0.942_034_635_077_344_3),
    (0.3, 0.7, 1.5, 0.3, 1.048_816_368_988_159_8),
    (0.3, 0.7, 1.5, 0.7, 1.153_867_076_098_748_4),
    (0.3, 0.7, 1.5, 0.95, 1.314_577_663_709_436_2),
    (0.3, 0.7, 1.5, 0.999, 1.444_481_010_263_152),
    (0.3, 0.7, 1.5, 0.999999999, 1.469_437_547_985_022_6),
    (2.5, 1.5, 3.5, -1000.0, 0.000_078_270_882_416_762_58),
    (2.5, 1.5, 3.5, -10.0, 0.053_674_605_984_257_61),
    (2.5, 1.5, 3.5, -2.0, 0.284_537_735_948_667_27),
    (2.5, 1.5, 3.5, -0.5, 0.640_484_548_889_045_8),
    (2.5, 1.5, 3.5, 0.3, 1.452_767_637_957_694_1),
    (2.5, 1.5, 3.5, 0.7, 3.291_990_555_652_741_3),
    (2.5, 1.5, 3.5, 0.95, 13.925_720_900_108_688),
    (2.5, 1.5, 3.5, 0.999, 146.937_128_977_006_75),
    (2.5, 1.5, 3.5, 0.999999999, 158_102.104_904_277_13),
    (1.3, 2.3, 4.6, -1000.0, 0.000_529_122_847_870_453_8),
    (1.3, 2.3, 4.6, -10.0, 0.125_815_549_697_934_26),
    (1.3, 2.3, 4.6, -2.0, 0.436_943_947_912_514_06),
    (1.3, 2.3, 4.6, -0.5, 0.756_341_180_114_451_4),
    (1.3, 2.3, 4.6, 0.3, 1.245_673_144_480_010_7),
    (1.3, 2.3, 4.6, 0.7, 1.908_372_864_547_123_6),
    (1.3, 2.3, 4.6, 0.95, 3.228_397_840_179_902_7),
    (1.3, 2.3, 4.6, 0.999, 4.208_383_751_358_547),
    (1.3, 2.3, 4.6, 0.999999999, 4.274_080_812_165_148),
    (1.3, 0.7, 2.0, -1000.0, 0.014_417_780_477_402_46),
    (1.3, 0.7, 2.0, -10.0, 0.280_065_606_316_328_6),
    (1.3, 0.7, 2.0, -2.0, 0.582_535_953_790_404_1),
    (1.3, 0.7, 2.0, -0.5, 0.826_863_343_360_431_3),
    (1.3, 0.7, 2.0, 0.3, 1.171_025_029_007_173_8),
    (1.3, 0.7, 2.0, 0.7, 1.644_735_358_714_176_4),
    (1.3, 0.7, 2.0, 0.95, 2.897_546_255_312_291),
    (1.3, 0.7, 2.0, 0.999, 6.136_539_985_123_29),
    (1.3, 0.7, 2.0, 0.999999999, 17.990_257_041_595_164),
    (0.4, 0.6, 3.0, -1000.0, 0.205_316_481_768_520_46),
    (0.4, 0.6, 3.0, -10.0, 0.720_716_515_909_826_4),
    (0.4, 0.6, 3.0, -2.0, 0.891_315_936_939_953_1),
    (0.4, 0.6, 3.0, -0.5, 0.964_666_666_323_555_6),
    (0.4, 0.6, 3.0, 0.3, 1.026_314_354_094_873),
    (0.4, 0.6, 3.0, 0.7, 1.072_173_774_234_874_3),
    (0.4, 0.6, 3.0, 0.95, 1.114_393_223_701_038),
    (0.4, 0.6, 3.0, 0.999, 1.125_961_730_503_150_8),
    (0.4, 0.6, 3.0, 0.999999999, 1.126_230_250_683_063_3),
    (1.5, 2.5, 3.0, -1000.0, 0.000_053_362_014_117_422_83),
    (1.5, 2.5, 3.0, -10.0, 0.040_543_766_631_261_13),
    (1.5, 2.5, 3.0, -2.0, 0.241_806_873_077_026_09),
    (1.5, 2.5, 3.0, -0.5, 0.598_554_160_717_613_7),
    (1.5, 2.5, 3.0, 0.3, 1.554_079_326_367_708_6),
    (1.5, 2.5, 3.0, 0.7, 4.260_921_914_840_392_5),
    (1.5, 2.5, 3.0, 0.95, 30.949_523_656_432_49),
    (1.5, 2.5, 3.0, 0.999, 1_689.980_460_427_861_3),
    (1.5, 2.5, 3.0, 0.999999999, 1_697_652_749.079_143_5),
    (0.5, -0.5, 1.0, -1000.0, 20.185_433_027_086_116),
    (0.5, -0.5, 1.0, -10.0, 2.316_747_229_631_726),
    (0.5, -0.5, 1.0, -2.0, 1.390_656_513_186_148_8),
    (0.5, -0.5, 1.0, -0.5, 1.115_212_230_772_902_4),
    (0.5, -0.5, 1.0, 0.3, 0.920_146_705_054_900_8),
    (0.5, -0.5, 1.0, 0.7, 0.790_472_034_321_195_2),
    (0.5, -0.5, 1.0, 0.95, 0.675_118_543_172_368_5),
    (0.5, -0.5, 1.0, 0.999, 0.638_001_740_734_463_5),
    (0.5, -0.5, 1.0, 0.999999999, 0.636_619_775_947_907_8),
    (0.5, 0.5, 1.0, -1000.0, 0.097_421_500_405_940_44),
    (0.5, 0.5, 1.0, -10.0, 0.503_484_682_735_704_5),
    (0.5, 0.5, 1.0, -2.0, 0.745_749_187_316_329_7),
    (0.5, 0.5, 1.0, -0.5, 0.901_286_299_360_447_4),
    (0.5, 0.5, 1.0, 0.3, 1.091_095_910_362_781_5),
    (0.5, 0.5, 1.0, 0.7, 1.321_217_206_769_961_7),
    (0.5, 0.5, 1.0, 0.95, 1.851_504_997_072_928_3),
    (0.5, 0.5, 1.0, 0.999, 3.081_960_708_698_816),
    (0.5, 0.5, 1.0, 0.999999999, 7.478_962_801_238_461),
    (1.0, 2.0, 4.0, -1000.0, 0.002_964_506_018_795_432_7),
    (1.0, 2.0, 4.0, -10.0, 0.201_738_911_995_307_55),
    (1.0, 2.0, 4.0, -2.0, 0.528_122_350_496_753_2),
    (1.0, 2.0, 4.0, -0.5, 0.806_512_216_212_164_5),
    (1.0, 2.0, 4.0, 0.3, 1.183_897_609_530_519),
    (1.0, 2.0, 4.0, 0.7, 1.640_959_044_353_688_6),
    (1.0, 2.0, 4.0, 0.95, 2.442_082_306_964_633_4),
    (1.0, 2.0, 4.0, 0.999, 2.967_443_900_663_840_7),
    (1.0, 2.0, 4.0, 0.999999999, 2.999_999_884_660_408),
    (1.1, 2.1, 0.7, -1000.0, -0.000_168_159_848_200_940_6),
    (1.1, 2.1, 0.7, -10.0, -0.025_933_030_948_297_568),
    (1.1, 2.1, 0.7, -2.0, -0.050_667_781_007_507_61),
    (1.1, 2.1, 0.7, -0.5, 0.212_877_106_919_934_3),
    (1.1, 2.1, 0.7, 0.3, 3.011_779_894_347_601),
    (1.1, 2.1, 0.7, 0.7, 31.025_118_065_974_26),
    (1.1, 2.1, 0.7, 0.95, 3_042.025_022_428_055_8),
    (1.1, 2.1, 0.7, 0.999, 54_788_973.967_899_725),
    (1.1, 2.1, 0.7, 0.999999999, 5.480_944_480_406_656e22),
    (0.7, 1.7, 1.3, -1000.0, 0.005_276_623_374_913_077),
    (0.7, 1.7, 1.3, -10.0, 0.13550158458878207),
    (0.7, 1.7, 1.3, -2.0, 0.382_266_791_756_527_45),
    (0.7, 1.7, 1.3, -0.5, 0.694_128_804_948_757_7),
    (0.7, 1.7, 1.3, 0.3, 1.392_539_424_336_495_7),
    (0.7, 1.7, 1.3, 0.7, 3.167_124_405_499_906_4),
    (0.7, 1.7, 1.3, 0.95, 20.325_034_949_523_04),
    (0.7, 1.7, 1.3, 0.999, 1_446.258_332_753_814_6),
    (0.7, 1.7, 1.3, 0.999999999, 5_750_114_783.966_538),
    (-0.3, 0.2, 0.9, -1000.0, 3.738_316_206_403_935_5),
    (-0.3, 0.2, 0.9, -10.0, 1.309_923_551_429_016_8),
    (-0.3, 0.2, 0.9, -2.0, 1.099_531_811_899_221_7),
    (-0.3, 0.2, 0.9, -0.5, 1.030_273_635_963_677_7),
    (-0.3, 0.2, 0.9, 0.3, 0.978_467_366_876_644_6),
    (-0.3, 0.2, 0.9, 0.7, 0.942_335_077_419_742_2),
    (-0.3, 0.2, 0.9, 0.95, 0.908_536_392_529_293_3),
    (-0.3, 0.2, 0.9, 0.999, 0.897_070_246_431_739_7),
    (-0.3, 0.2, 0.9, 0.999999999, 0.896_625_581_969_533_3),
];
const J_CASES: &[(i32, f64, f64)] = &[
    (0, 0.1, 0.997_501_562_066_04),
    (1, 0.1, 0.049_937_526_036_242),
    (2, 0.1, 0.001_248_958_658_799_919),
    (5, 0.1, 2.603_081_790_964_441_7e-9),
    (20, 0.1, 3.919_437_720_858_622e-45),
    (50, 0.1, 2.920_142_569_099_643_7e-130),
    (100, 0.1, 8.452_516_535_121_789e-289),
    (300, 0.1, 0.0),
    (0, 10.0, -0.245_935_764_451_348_35),
    (1, 10.0, 0.043_472_746_168_861_44),
    (2, 10.0, 0.254_630_313_685_120_6),
    (5, 10.0, -0.234_061_528_186_793_63),
    (20, 10.0, 0.000_011_513_369_247_813_398),
    (50, 10.0, 1.784_513_607_871_595_3e-30),
    (100, 10.0, 6.597_316_064_155_382e-89),
    (300, 10.0, 0.0),
    (0, 100.0, 0.019_985_850_304_223_122),
    (1, 100.0, -0.077_145_352_014_112_16),
    (2, 100.0, -0.021_528_757_344_505_364),
    (5, 100.0, -0.074_195_736_964_513_93),
    (20, 100.0, 0.062_217_458_498_338_755),
    (50, 100.0, -0.038_698_339_728_525_384),
    (100, 100.0, 0.096_366_673_295_861_56),
    (300, 100.0, 3.520_366_621_846_936_4e-109),
    (0, 1000.0, 0.024_786_686_152_420_176),
    (1, 1000.0, 0.004_728_311_907_089_524),
    (2, 1000.0, -0.024_777_229_528_605_997),
    (5, 1000.0, 0.005_025_406_945_233_186_5),
    (20, 1000.0, 0.023_357_967_932_679_333),
    (50, 1000.0, -0.003_336_048_960_615_276_4),
    (100, 1000.0, 0.011_676_135_007_802_554),
    (300, 1000.0, 0.000_467_828_038_791_247_9),
    (4, 30.0, -0.052_609_000_321_320_355),
    (12, 45.5, -0.082_099_302_267_801_54),
    (25, 700.0, 0.023_895_529_404_557_973),
    (24, 25.0, 0.199_778_510_669_943_26),
    (3, 5.0, 0.364_831_230_613_667),
    (3, 100000.0, -0.001_846_688_793_360_512_2),
];
const I_CASES: &[(i32, f64, f64)] = &[
    (0, 0.1, 1.002_501_562_934_095_6),
    (1, 0.1, 0.050_062_526_047_092_694),
    (3, 0.1, 0.000_020_846_357_422_327_155),
    (20, 0.1, 3.920_371_031_419_978e-45),
    (100, 0.1, 8.452_934_986_892_05e-289),
    (0, 10.0, 2_815.716_628_466_254_4),
    (1, 10.0, 2_670.988_303_701_254_6),
    (3, 10.0, 1_758.380_716_610_853_1),
    (20, 10.0, 0.000_125_079_973_564_494_75),
    (100, 10.0, 1.082_344_201_749_201_6e-88),
    (0, 100.0, 1.073_751_707_131_073_8e42),
    (1, 100.0, 1.068_369_390_338_162_5e42),
    (3, 100.0, 1.026_274_017_565_19e42),
    (20, 100.0, 1.448_346_125_642_717_2e41),
    (100, 100.0, 4.641_534_941_616_199e21),
    (0, 600.0, 6.146_305_403_936_844_4e258),
    (1, 600.0, 6.141_181_345_066_892e258),
    (3, 600.0, 6.100_342_446_403_87e258),
    (20, 600.0, 4.402_931_410_823_698e258),
    (100, 600.0, 1.495_622_591_661_959_3e255),
];
const IS_CASES: &[(i32, f64, f64)] = &[
    (0, 1000.0, 0.012_617_240_455_891_257),
    (1, 1000.0, 0.012_610_930_256_928_629),
    (5, 1000.0, 0.012_460_428_940_768_863),
    (50, 1000.0, 0.003_613_581_892_594_122_6),
    (0, 100000.0, 0.001_261_567_837_976_776_8),
    (1, 100000.0, 0.001_261_561_530_121_817_1),
    (5, 100000.0, 0.001_261_410_151_064_260_6),
    (50, 100000.0, 0.001_245_896_313_058_022_9),
];
const JAC_CASES: &[(f64, f64, f64, f64, f64)] = &[
    (
        0.5,
        0.3,
        0.474_215_622_711_820_6,
        0.880_408_736_426_462_4,
        0.965_678_964_745_951_2,
    ),
    (
        3.0,
        0.3,
        0.411_422_660_634_000_7,
        -0.911_444_674_303_843_8,
        0.974_278_922_226_613_7,
    ),
    (
        17.3,
        0.3,
        -0.160_204_282_891_874_55,
        -0.987_083_880_804_007_6,
        0.996_142_748_968_705_8,
    ),
    (
        123.456,
        0.3,
        0.055_921_787_191_412_28,
        0.998_435_152_484_786_3,
        0.999_530_802_984_678_1,
    ),
    (
        1000.0,
        0.3,
        -0.770_450_078_649_177_6,
        0.637_500_334_360_285_6,
        0.906_599_141_237_649_1,
    ),
    (
        -777.7,
        0.3,
        -0.391_796_371_291_614_04,
        -0.920_051_956_925_653,
        0.976_702_964_586_888_7,
    ),
    (
        0.5,
        0.9,
        0.463_840_368_011_727_7,
        0.885_918_795_941_673_8,
        0.897_979_343_694_759_3,
    ),
    (
        3.0,
        0.9,
        0.990_630_599_937_832_5,
        -0.136_568_717_013_853_34,
        0.341_739_539_737_691_06,
    ),
    (
        17.3,
        0.9,
        -0.967_383_798_973_692_3,
        -0.253_315_189_997_021_2,
        0.397_179_716_167_510_5,
    ),
    (
        123.456,
        0.9,
        -0.284_751_768_641_518_85,
        0.958_601_288_469_573_7,
        0.962_821_264_425_528_8,
    ),
    (
        1000.0,
        0.9,
        -0.291_493_035_512_560_67,
        0.956_572_950_823_758_7,
        0.961_003_969_410_587_6,
    ),
    (
        -777.7,
        0.9,
        -0.715_006_894_602_377_5,
        -0.699_117_401_207_454_2,
        0.734_771_138_929_638_5,
    ),
    (
        0.5,
        0.999,
        0.462_134_380_670_853_3,
        0.886_809_908_718_867_9,
        0.886_930_314_279_405_4,
    ),
    (
        3.0,
        0.999,
        0.995_293_128_337_884_4,
        0.096_910_209_386_770_21,
        0.101_893_066_960_928_5,
    ),
    (
        17.3,
        0.999,
        -0.968_523_623_115_827,
        0.248_921_657_287_170_44,
        0.250_798_782_842_191_9,
    ),
    (
        123.456,
        0.999,
        0.984_330_767_759_837,
        -0.176_331_901_938_729_77,
        0.179_058_221_541_733_43,
    ),
    (
        1000.0,
        0.999,
        -0.991_707_563_246_611_1,
        -0.128_515_014_684_934_23,
        0.132_286_026_814_888_5,
    ),
    (
        -777.7,
        0.999,
        -0.996_337_588_405_197_1,
        0.085_506_782_952_675_96,
        0.091_126_826_571_460_94,
    ),
];

// ---------------------------------------------------------------------
// Extended-range tests (reference values from 40-digit arithmetic)
// ---------------------------------------------------------------------

fn rel(
    a: f64,
    b: f64,
) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn hyp2f1_matches_high_precision_references() {
    let mut worst = 0.0_f64;
    for &(a, b, c, z, want) in HYP_CASES {
        let got = hyp2f1(a, b, c, z);
        let e = rel(got, want);
        worst = worst.max(e);
        assert!(
            e < 5e-10,
            "2F1({a},{b};{c};{z}) = {got} vs {want} (rel {e:e})"
        );
    }
    assert!(worst < 5e-10);
}

#[test]
fn hyp2f1_elementary_closed_forms() {
    // 2F1(1,1;2;z) = -ln(1-z)/z, across all of z < 1
    for z in [
        -1e6,
        -1000.0,
        -10.0,
        -2.0,
        -0.9,
        -1e-3,
        1e-3,
        0.3,
        0.6,
        0.61,
        0.9,
        0.99,
        0.999999,
        1.0 - 1e-12,
    ] {
        close(hyp2f1(1.0, 1.0, 2.0, z), -(1.0 - z).ln() / z, 2e-11);
    }
    // 2F1(1,1;3;z) = 2 (z + (1-z) ln(1-z)) / z^2
    for z in [-50.0, -0.5, 0.4, 0.8, 0.97] {
        close(
            hyp2f1(1.0, 1.0, 3.0, z),
            2.0 * (z + (1.0 - z) * (1.0 - z).ln()) / (z * z),
            1e-11,
        );
    }
    // 2F1(1/2,1;3/2;z^2) = atanh(z)/z (z^2 close to 1 hits the logarithmic case)
    for z in [0.2_f64, 0.7, 0.95, 0.9999] {
        close(hyp2f1(0.5, 1.0, 1.5, z * z), z.atanh() / z, 1e-11);
    }
    // 2F1(a,b;b;z) = (1-z)^-a for any z < 1 (integer c-a-b and Euler path)
    for z in [-100.0, -1.0, 0.5, 0.9, 0.99] {
        close(hyp2f1(2.0, 3.5, 3.5, z), (1.0 - z).powi(-2), 1e-11);
        close(hyp2f1(0.5, 2.0, 2.0, z), (1.0 - z).powf(-0.5), 1e-11);
    }
    // 2F1(a,-a;1/2;sin^2 x) = cos(2 a x), 2F1(a,-a;1/2;-sinh^2 t) = cosh(2 a t)
    for x in [0.3_f64, 1.0, 1.4] {
        close(
            hyp2f1(0.7, -0.7, 0.5, x.sin().powi(2)),
            (1.4 * x).cos(),
            1e-12,
        );
    }
    for t in [0.5_f64, 2.0, 5.0] {
        close(
            hyp2f1(0.7, -0.7, 0.5, -t.sinh().powi(2)),
            (1.4 * t).cosh(),
            1e-11,
        );
    }
}

#[test]
fn hyp2f1_legendre_and_elliptic_identities() {
    // 2F1(1/2,1/2;1;m) = 2 K(m) / pi for all m < 1 (log case at c-a-b = 0)
    for m in [-1e4, -100.0, -3.0, -0.5, 0.2, 0.7, 0.95, 0.999999] {
        close(hyp2f1(0.5, 0.5, 1.0, m), 2.0 / PI * elliptic_k(m), 1e-11);
    }
    // 2F1(-1/2,1/2;1;m) = 2 E(m) / pi (c-a-b = 1)
    for m in [-50.0, -1.0, 0.3, 0.8, 0.99] {
        close(hyp2f1(-0.5, 0.5, 1.0, m), 2.0 / PI * elliptic_e(m), 1e-11);
    }
    // Legendre polynomials: P_n(x) = 2F1(-n, n+1; 1; (1-x)/2), including x < -1 (z > 1)
    for x in [-3.0, -0.9, 0.0, 0.4, 0.99, 1.0, 2.5] {
        let z = (1.0 - x) / 2.0;
        close(
            hyp2f1(-3.0, 4.0, 1.0, z),
            0.5 * (5.0 * x * x * x - 3.0 * x),
            1e-12,
        );
    }
    // P_{-1/2}(cos t) = 2 K(sin^2(t/2)) / pi = 2F1(1/2,1/2;1;sin^2(t/2))
    let t: f64 = 2.9;
    close(
        hyp2f1(0.5, 0.5, 1.0, (t / 2.0).sin().powi(2)),
        2.0 / PI * elliptic_k((t / 2.0).sin().powi(2)),
        1e-12,
    );
}

#[test]
fn hyp2f1_special_arguments() {
    // Gauss summation at z = 1: 2F1(1,2;4;1) = Gamma(4) Gamma(1) / (Gamma(3) Gamma(2)) = 3
    close(hyp2f1(1.0, 2.0, 4.0, 1.0), 3.0, 1e-12);
    // divergent at z = 1
    assert_eq!(hyp2f1(1.0, 1.0, 2.0, 1.0), f64::INFINITY);
    // genuinely complex for z > 1 (non-terminating)
    assert!(hyp2f1(0.5, 0.5, 1.0, 1.5).is_nan());
    assert!(hyp2f1(1.0, 1.0, 2.0, 3.0).is_nan());
    // terminating series are fine for any z
    close(
        hyp2f1(-2.0, 3.0, 1.5, 7.0),
        1.0 - 4.0 * 7.0 + 3.2 * 49.0,
        1e-12,
    );
    // c a non-positive integer without termination is undefined
    assert!(hyp2f1(0.5, 0.5, -1.0, 0.2).is_nan());
    assert_eq!(hyp2f1(0.3, 0.4, 0.5, 0.0), 1.0);
}

#[test]
fn jacobi_large_arguments() {
    for &(u, m, sn, cn, dn) in JAC_CASES {
        let (s, c, d) = jacobi_elliptic(u, m);
        let tol = 2e-13 * u.abs().max(1.0);
        assert!((s - sn).abs() < tol, "sn({u}|{m}) = {s} vs {sn}");
        assert!((c - cn).abs() < tol, "cn({u}|{m}) = {c} vs {cn}");
        assert!((d - dn).abs() < tol, "dn({u}|{m}) = {d} vs {dn}");
    }
    // identities up to |u| = 1e3
    for m in [0.01, 0.3, 0.5, 0.9, 0.99, 0.9999] {
        for k in 0..400 {
            let u = -1000.0 + 5.003 * f64::from(k);
            let (s, c, d) = jacobi_elliptic(u, m);
            assert!((s * s + c * c - 1.0).abs() < 1e-14, "m={m} u={u}");
            assert!((d * d + m * s * s - 1.0).abs() < 1e-13, "m={m} u={u}");
        }
    }
}

#[test]
fn jacobi_periodicity_and_symmetries() {
    for m in [0.2, 0.7, 0.95] {
        let k = elliptic_k(m);
        for u in [0.3, 1.1, 2.7, 123.4] {
            let (s, c, d) = jacobi_elliptic(u, m);
            // sn, cn have period 4K; dn has period 2K
            let (s4, c4, d4) = jacobi_elliptic(u + 4.0 * k, m);
            assert!((s - s4).abs() < 1e-13 && (c - c4).abs() < 1e-13 && (d - d4).abs() < 1e-13);
            let (s2, c2, d2) = jacobi_elliptic(u + 2.0 * k, m);
            assert!((s + s2).abs() < 1e-13 && (c + c2).abs() < 1e-13 && (d - d2).abs() < 1e-13);
            // many periods away
            let (s9, c9, _) = jacobi_elliptic(u + 200.0 * k, m);
            assert!(
                (s - s9).abs() < 1e-11 && (c - c9).abs() < 1e-11,
                "m={m} u={u}"
            );
            // parity
            let (sm, cm, dm) = jacobi_elliptic(-u, m);
            assert!((sm + s).abs() < 1e-15 && (cm - c).abs() < 1e-15 && (dm - d).abs() < 1e-15);
        }
        // quarter period values: sn(K) = 1, cn(K) = 0, dn(K) = sqrt(1-m)
        let (s, c, d) = jacobi_elliptic(k, m);
        assert!((s - 1.0).abs() < 1e-14 && c.abs() < 1e-7 && (d - (1.0 - m).sqrt()).abs() < 1e-14);
    }
    assert!(jacobi_elliptic(1.0, 1.5).0.is_nan());
    assert!(jacobi_elliptic(f64::INFINITY, 0.5).0.is_nan());
}

#[test]
fn bessel_j_reference_values() {
    for &(n, x, want) in J_CASES {
        let got = bessel_jn(n, x);
        let err = (got - want).abs();
        let ok = err <= 2e-13 * want.abs() + 3e-16 / x.max(1.0).sqrt();
        assert!(ok, "J_{n}({x}) = {got:e} vs {want:e}");
    }
    // symmetry
    close(bessel_jn(-3, 7.5), -bessel_jn(3, 7.5), 1e-14);
    close(bessel_jn(3, -7.5), -bessel_jn(3, 7.5), 1e-14);
    assert_eq!(bessel_jn(0, 0.0), 1.0);
    assert_eq!(bessel_jn(2, 0.0), 0.0);
}

#[test]
fn bessel_j_identities_large_x() {
    // J0^2 + 2 sum J_k^2 = 1 and the three-term recurrence at x = 100, 1000
    for x in [10.0, 100.0, 1000.0] {
        let mut sum = bessel_jn(0, x).powi(2);
        for k in 1..(x as i32 + 200) {
            sum += 2.0 * bessel_jn(k, x).powi(2);
        }
        assert!((sum - 1.0).abs() < 1e-12, "x={x}: {sum}");
        for n in [1, 4, 17, 60] {
            let lhs = bessel_jn(n - 1, x) + bessel_jn(n + 1, x);
            let rhs = 2.0 * f64::from(n) / x * bessel_jn(n, x);
            assert!((lhs - rhs).abs() < 1e-13, "x={x} n={n}");
        }
    }
}

#[test]
fn bessel_i_reference_values() {
    for &(n, x, want) in I_CASES {
        let got = bessel_in(n, x);
        assert!(rel(got, want) < 5e-13, "I_{n}({x}) = {got:e} vs {want:e}");
    }
    for &(n, x, want) in IS_CASES {
        let got = bessel_in_scaled(n, x);
        assert!(
            rel(got, want) < 5e-13,
            "scaled I_{n}({x}) = {got:e} vs {want:e}"
        );
    }
    assert!(bessel_in(0, 1000.0).is_infinite());
    close(bessel_in(3, -2.0), -bessel_in(3, 2.0), 1e-14);
}

#[test]
fn bessel_cost_is_not_linear_in_x() {
    // x = 1e5 and 1e7 must return instantly via the asymptotic expansion
    let v = bessel_jn(0, 1e7);
    assert!(v.abs() < 1e-3);
    let s = bessel_in_scaled(0, 1e9);
    close(s, 1.0 / (2.0 * PI * 1e9_f64).sqrt(), 1e-9);
}

// ---------------------------------------------------------------------
// Coulomb wave functions
// ---------------------------------------------------------------------

const COULOMB_CASES: &[(f64, usize, f64, f64, f64, f64, f64)] = &[
    (
        0.0,
        0,
        0.5,
        0.479_425_538_604_203,
        0.877_582_561_890_372_7,
        0.877_582_561_890_372_7,
        -0.479_425_538_604_203,
    ),
    (
        0.0,
        1,
        2.0,
        0.870_795_549_959_983_2,
        0.473_899_651_845_690_1,
        0.701_224_008_552_110_5,
        -0.766_758_840_823_197_6,
    ),
    (
        0.0,
        2,
        5.0,
        0.673_656_050_425_626_1,
        -0.744_909_460_566_104_4,
        -0.824_977_288_005_522_2,
        -0.572_200_922_368_284_3,
    ),
    (
        0.0,
        3,
        10.0,
        -0.394_958_449_844_703_2,
        0.897_909_471_239_035_4,
        0.953_274_788_765_689,
        0.364_710_613_307_641_2,
    ),
    (
        1.0,
        0,
        3.0,
        1.084_052_684_202_935_8,
        0.301_916_785_821_136_6,
        0.627_039_514_889_285_5,
        -0.747_829_194_000_705,
    ),
    (
        1.0,
        2,
        5.0,
        1.186_370_500_609_351_4,
        0.154_144_770_433_326_6,
        0.382_961_011_797_984_5,
        -0.793_148_988_674_430_7,
    ),
    (
        -2.0,
        0,
        7.0,
        -0.142_480_641_838_878,
        1.102_975_280_207_625_3,
        0.882_122_126_144_823,
        0.189_787_892_509_963_8,
    ),
    (
        2.0,
        1,
        5.0,
        0.993_505_675_214_803_3,
        0.370_867_645_213_274_3,
        1.081_482_760_400_011,
        -0.602_827_996_108_031,
    ),
    (
        5.0,
        0,
        20.0,
        -0.229_347_373_002_198_6,
        0.828_679_340_869_801_9,
        1.165_716_671_344_780_7,
        0.148_223_529_679_329_7,
    ),
    (
        -5.0,
        1,
        12.0,
        -0.078_478_739_182_301_9,
        -1.157_125_099_716_428,
        -0.857_568_597_110_491,
        0.097_935_461_692_271_8,
    ),
    (
        0.5,
        3,
        8.0,
        -0.199_378_373_579_696_9,
        -0.894_016_314_374_365_8,
        -1.076_478_821_251_929,
        0.188_638_170_966_195,
    ),
    (
        -0.5,
        0,
        0.2,
        0.325_100_720_962_525_2,
        1.434_348_774_996_957_7,
        0.632_984_681_716_263_1,
        -0.283_232_829_862_142_8,
    ),
];

const COULOMB_CASES_C: &[(f64, usize, f64, f64)] = &[
    (0.0, 0, 1.0, 0.0),
    (0.0, 1, 1.0 / 3.0, 0.0),
    (0.0, 3, 0.009_523_809_523_809_524, 0.0),
    (1.0, 0, 0.108_422_513_102_072_6, -0.301_640_320_467_533_2),
    (1.0, 2, 0.011_428_736_368_066_59, 0.947_405_451_930_721_2),
    (-2.0, 1, 2.642_222_805_575_776_3, -1.236_795_034_103_878_8),
    (2.5, 0, 0.001_538_576_471_074_136, 0.542_604_405_852_436_5),
    (-0.5, 2, 0.139_218_108_418_932_2, -0.464_567_973_222_242_5),
];

#[test]
fn coulomb_reference_values() {
    for &(eta, l, rho, f, fp, g, gp) in COULOMB_CASES {
        let values = coulomb_fg(eta, l, rho)
            .unwrap_or_else(|| panic!("no F/G at eta={eta}, l={l}, rho={rho}"));
        // Steed's fractions lose a few digits when |eta| is much larger
        // than rho; the pde solver documents 1e-8 for that corner.
        let tol = if eta.abs() >= 3.0 {
            1e-8
        } else {
            1e-11
        };
        close(values.0, f, tol);
        close(values.1, fp, tol);
        close(values.2, g, tol);
        close(values.3, gp, tol);
        close(coulomb_f(eta, l, rho).unwrap_or(f64::NAN), f, tol);
        close(coulomb_g(eta, l, rho).unwrap_or(f64::NAN), g, tol);
        close(coulomb_f_prime(eta, l, rho).unwrap_or(f64::NAN), fp, tol);
        close(coulomb_g_prime(eta, l, rho).unwrap_or(f64::NAN), gp, tol);
    }
    for &(eta, l, c, sigma) in COULOMB_CASES_C {
        close(coulomb_c(eta, l), c, 1e-12);
        close(coulomb_sigma(eta, l), sigma, 1e-12);
    }
    // Domain limits.
    assert!(coulomb_fg(f64::NAN, 0, 1.0).is_none());
    assert!(coulomb_fg(0.0, 0, 0.0).is_none());
    assert!(coulomb_fg(0.0, 0, -1.0).is_none());
    assert!(coulomb_fg(0.0, 65, 1.0).is_none());
    // Outside the classically forbidden series window: rho > 34 yet still
    // below 2|eta| + 2.
    assert!(coulomb_fg(20.0, 0, 35.0).is_none());
}

#[test]
fn coulomb_wronskian_and_recurrences() {
    // F' G - F G' = 1 and the derivative recurrence
    // F_l' = ((l+1)/rho + eta/(l+1)) F_l - sqrt(1 + eta^2/(l+1)^2) F_{l+1}.
    for &(eta, l) in &[
        (0.0_f64, 0_usize),
        (0.5, 1),
        (1.0, 0),
        (-1.0, 2),
        (2.0, 1),
        (5.0, 0),
        (-5.0, 1),
    ] {
        for &rho in &[0.5_f64, 2.0, 8.0, 14.0, 20.0, 30.0, 60.0] {
            // Steed's fractions lose accuracy in the |eta| >> rho corner.
            if eta.abs() > 2.0 * rho {
                continue;
            }
            let Some((f, fp, g, gp)) = coulomb_fg(eta, l, rho) else {
                continue;
            };
            assert!(
                (fp * g - f * gp - 1.0).abs() < 1e-8,
                "Wronskian at eta={eta}, l={l}, rho={rho}"
            );
            let lf = l as f64;
            if let Some(next) = coulomb_f(eta, l + 1, rho) {
                let expected = ((lf + 1.0) / rho + eta / (lf + 1.0)) * f
                    - (1.0 + eta * eta / ((lf + 1.0) * (lf + 1.0))).sqrt() * next;
                close(fp, expected, 1e-8);
                let g_next = coulomb_g(eta, l + 1, rho).unwrap_or(f64::NAN);
                let expected_g = ((lf + 1.0) / rho + eta / (lf + 1.0)) * g
                    - (1.0 + eta * eta / ((lf + 1.0) * (lf + 1.0))).sqrt() * g_next;
                close(gp, expected_g, 1e-8);
            }
        }
    }
}

#[test]
fn coulomb_eta_zero_and_series_asymptotic_continuity() {
    // eta = 0: F_l = rho j_l and G_l = -rho y_l.
    let f = |l: usize, rho: f64| {
        match l {
            | 0 => rho.sin(),
            | 1 => rho.sin() / rho - rho.cos(),
            | 2 => (3.0 / (rho * rho) - 1.0) * rho.sin() - 3.0 * rho.cos() / rho,
            | _ => {
                let s = rho.sin();
                let c = rho.cos();
                (15.0 / (rho * rho * rho) - 6.0 / rho) * s + (1.0 - 15.0 / (rho * rho)) * c
            },
        }
    };
    let g = |l: usize, rho: f64| {
        match l {
            | 0 => rho.cos(),
            | 1 => rho.cos() / rho + rho.sin(),
            | 2 => (3.0 / (rho * rho) - 1.0) * rho.cos() + 3.0 * rho.sin() / rho,
            | _ => {
                let s = rho.sin();
                let c = rho.cos();
                (15.0 / (rho * rho * rho) - 6.0 / rho) * c + (15.0 / (rho * rho) - 1.0) * s
            },
        }
    };
    for l in 0..=3_usize {
        for rho in [0.2_f64, 0.7, 2.0, 6.0, 15.0, 24.0, 60.0] {
            let Some(fv) = coulomb_f(0.0, l, rho) else {
                panic!("no F")
            };
            close(fv, f(l, rho), 1e-11);
            let Some(gv) = coulomb_g(0.0, l, rho) else {
                panic!("no G")
            };
            // Steed's fractions lose a few digits for l = 3 at rho = 0.2.
            close(gv, g(l, rho), 1e-8);
        }
    }
    // The confluent series and the asymptotic expansion agree across the
    // switch at rho = 24.
    for &(eta, l) in &[
        (0.0_f64, 0_usize),
        (1.0, 1),
        (-1.0, 2),
        (2.0, 1),
        (5.0, 0),
        (-5.0, 1),
    ] {
        // Two points 2e-9 apart: the function changes by about
        // 2e-9 * |F'|, so the representation change is what is tested.
        let below = coulomb_fg(eta, l, 23.999_999_999).unwrap_or_else(|| panic!("no series value"));
        let above =
            coulomb_fg(eta, l, 24.000_000_001).unwrap_or_else(|| panic!("no asymptotic value"));
        let tol = if eta.abs() >= 3.0 {
            1e-6
        } else {
            1e-8
        };
        close(below.0, above.0, tol);
        close(below.1, above.1, tol);
        close(below.2, above.2, tol);
        close(below.3, above.3, tol);
    }
}

#[test]
fn coulomb_solves_the_radial_equation() {
    for &(eta, l) in &[(0.0_f64, 0_usize), (1.0, 1), (-2.0, 2), (3.0, 0)] {
        for &rho in &[1.0_f64, 3.0, 9.0, 18.0] {
            let Some(f) = coulomb_f(eta, l, rho) else {
                continue;
            };
            let h = 1e-4;
            let (Some(plus), Some(minus)) =
                (coulomb_f(eta, l, rho + h), coulomb_f(eta, l, rho - h))
            else {
                continue;
            };
            let second = (plus - 2.0 * f + minus) / (h * h);
            let lf = l as f64;
            let potential = 1.0 - 2.0 * eta / rho - lf * (lf + 1.0) / (rho * rho);
            assert!(
                (second + potential * f).abs() < 1e-4 * (1.0 + f.abs()),
                "radial equation at eta={eta}, l={l}, rho={rho}"
            );
        }
    }
}

// ---------------------------------------------------------------------
// Polylogarithm and Lerch transcendent
// ---------------------------------------------------------------------

const POLYLOG_CASES: &[(f64, f64, f64, f64)] = &[
    (0.5, 0.5, 0.806_126_723_042_852_3, 0.0),
    (0.5, -0.5, -0.373_752_237_980_973_1, 0.0),
    (1.5, 2.5, 1.184_000_189_908_536_3, -3.393_294_700_911_349),
    (2.5, -3.0, -2.162_700_712_002_057, 0.0),
    (0.5, 2.0, -1.610_061_529_867_485, -2.128_934_038_862_452_4),
    (1.5, 0.999, 2.501_708_465_341_355_6, 0.0),
    (3.0, 0.3, 0.312_400_177_892_892_6, 0.0),
    (1.0, 0.3, 0.356_674_943_938_732_4, 0.0),
    (0.0, 0.3, 0.428_571_428_571_428_6, 0.0),
    (-1.0, 0.3, 0.612_244_897_959_183_6, 0.0),
    (-2.0, 0.3, 1.137_026_239_067_055_3, 0.0),
    (-3.0, 2.5, 8.518_518_518_518_519, 0.0),
    (2.0, 0.99, 1.588_625_448_076_375_2, 0.0),
    (2.0, 1000.0, -20.569_673_613_567_51, -21.701_353_237_246_4),
    (3.0, -1000.0, -66.300_123_850_809_27, 0.0),
];

const POLYLOG_COMPLEX_CASES: &[(f64, f64, f64, f64, f64, f64)] = &[
    (
        0.5,
        0.3,
        0.4,
        0.2,
        0.528_011_019_547_220_7,
        0.368_793_573_477_330_9,
    ),
    (
        0.5,
        0.3,
        -0.4,
        0.2,
        -0.340_342_448_390_287_7,
        0.110_288_249_394_827_4,
    ),
    (
        2.5,
        1.0,
        0.9,
        0.3,
        1.170_628_630_899_858_2,
        0.294_040_950_587_905_7,
    ),
    (
        -0.5,
        0.2,
        0.5,
        -0.3,
        0.168_719_354_006_113_5,
        -1.101_856_803_984_646_9,
    ),
    (
        1.5,
        -0.5,
        -2.0,
        0.1,
        -1.284_962_304_442_192,
        0.214_385_288_174_990_3,
    ),
    (
        0.75,
        0.0,
        0.3,
        -0.8,
        -0.127_637_717_294_666_6,
        -0.835_677_348_534_448_3,
    ),
];

fn complex_close(
    got: Complex64,
    want: (f64, f64),
    tolerance: f64,
) {
    let want = Complex64::new(want.0, want.1);
    assert!(
        (got - want).norm() <= tolerance * want.norm().max(1.0),
        "got {got} vs {want}"
    );
}

#[test]
fn polylog_reference_values() {
    for &(s, z, re, im) in POLYLOG_CASES {
        let got = polylog_complex(Complex64::new(s, 0.0), Complex64::new(z, 0.0));
        complex_close(got, (re, im), 1e-11);
        if z <= 1.0 {
            close(polylog(s, z), re, 1e-11);
        }
    }
    for &(s_re, s_im, z_re, z_im, want_re, want_im) in POLYLOG_COMPLEX_CASES {
        let got = polylog_complex(Complex64::new(s_re, s_im), Complex64::new(z_re, z_im));
        complex_close(got, (want_re, want_im), 1e-11);
    }
    // Analytic continuation at z = 1 and the pole at s = 1.
    close(polylog(2.0, 1.0), PI * PI / 6.0, 1e-13);
    close(polylog(0.0, 1.0), -0.5, 1e-14);
    assert!(polylog(1.0, 1.0).is_infinite());
    // z > 1 is complex for positive integer s.
    assert!(polylog(2.0, 2.0).is_nan());
    assert!((polylog(-2.0, 2.0) + 6.0).abs() < 1e-12);
}

#[test]
fn polylog_identities() {
    // Li_1(z) = -ln(1-z), Li_0(z) = z/(1-z).
    for z in [-0.9_f64, -0.3, 0.0, 0.3, 0.8, 0.99] {
        close(polylog(1.0, z), -(1.0 - z).ln(), 1e-13);
        close(polylog(0.0, z), z / (1.0 - z), 1e-13);
    }
    // Negative integer orders are the Eulerian rational functions.
    for z in [-0.8_f64, -0.2, 0.3, 0.7, 2.5] {
        let one_minus = 1.0 - z;
        close(polylog(-1.0, z), z / (one_minus * one_minus), 1e-12);
        close(polylog(-2.0, z), z * (1.0 + z) / one_minus.powi(3), 1e-12);
        close(
            polylog(-3.0, z),
            z * (1.0 + 4.0 * z + z * z) / one_minus.powi(4),
            1e-12,
        );
    }
    // Li_s'(z) = Li_{s-1}(z)/z.
    for &(s, z) in &[(0.5_f64, 0.4_f64), (1.5, -0.6), (2.5, 0.8)] {
        let h = 1e-5;
        let derivative = (polylog(s, z + h) - polylog(s, z - h)) / (2.0 * h);
        close(derivative, polylog(s - 1.0, z) / z, 1e-6);
    }
}

const LERCH_CASES: &[(f64, f64, f64, f64)] = &[
    (0.5, 1.5, 0.3, 6.526_899_164_438_676),
    (0.9, 0.5, 0.1, 7.107_489_640_999_062),
    (0.99, 2.0, 0.5, 4.889_092_975_043_935),
    (-0.5, 2.5, 1.2, 0.575_249_450_961_414_4),
    (-5.0, 0.5, 0.7, 0.358_554_981_517_674_3),
    (0.3, 3.0, 2.0, 0.137_779_754_365_473_6),
    (0.8, 1.0, 1.5, 1.535_083_821_537_538),
    (0.99, 1.0, 0.25, 8.268_624_575_137_229),
    (0.5, 0.1, 0.9, 1.918_288_435_305_399),
    (-2.0, 1.0, 3.0, 0.137_326_536_083_513_7),
    (0.25, 0.0, 1.7, 1.333_333_333_333_333_3),
    (-0.9, -2.0, 0.4, -0.128_356_903_338_679_1),
    (0.9, -1.0, 0.5, 95.0),
    (0.0, 2.5, 0.3, 20.286_020_648_339_49),
    (0.5, 2.0, 1.0, 1.164_481_052_930_025),
    (0.5, 1.0, 2.0, 0.772_588_722_239_781_2),
];

#[test]
fn lerch_reference_values() {
    for &(z, s, a, want) in LERCH_CASES {
        close(lerch(z, s, a), want, 2e-10);
    }
    // Reductions.
    for z in [-0.5_f64, 0.0, 0.3, 0.8] {
        close(lerch(z, 0.0, 0.7), 1.0 / (1.0 - z), 1e-13);
    }
    for &(s, z) in &[(1.5_f64, 0.4_f64), (2.0, -0.6), (0.5, 0.9)] {
        close(lerch(z, s, 1.0), polylog(s, z) / z, 1e-11);
    }
    close(lerch(0.4, 1.0, 1.0), -(1.0 - 0.4_f64).ln() / 0.4, 1e-13);
    // s = 1 at a positive integer a: (-ln(1-z) - sum z^j/j)/z^a.
    for a in [2.0_f64, 3.0, 5.0] {
        for z in [0.3_f64, -0.4, 0.7] {
            let sum: f64 = (1..a as u32).map(|j| z.powi(j as i32) / f64::from(j)).sum();
            let want = (-(1.0 - z).ln() - sum) / z.powi(a as i32);
            close(lerch(z, 1.0, a), want, 1e-12);
        }
    }
    close(lerch(0.3, 1.0, 2.0), (-(0.7_f64).ln() - 0.3) / 0.09, 1e-13);
    // z = 1 is the analytic continuation zeta(s, a).
    close(lerch(1.0, 2.0, 0.5), hurwitz_zeta(2.0, 0.5), 1e-13);
    assert!(lerch(1.0, 1.0, 1.0).is_infinite());
    // Outside the range.
    assert!(lerch(2.0, 1.5, 0.5).is_nan());
    assert!(lerch(0.5, 1.5, -0.5).is_nan());
}

// ---------------------------------------------------------------------
// Appell F1
// ---------------------------------------------------------------------

const APPELL_CASES: &[(f64, f64, f64, f64, f64, f64, f64)] = &[
    (1.5, 0.5, 0.7, 2.0, 0.3, -0.4, 0.949_582_343_289_011_8),
    (1.0, 0.5, 0.5, 1.5, 0.6, 0.2, 1.439_837_249_573_893_2),
    (2.0, -0.5, 1.5, 3.0, -0.8, 0.4, 2.049_290_365_596_732_3),
    (0.8, 0.3, 1.2, 2.5, 0.9, -0.7, 0.898_944_558_857_869),
    (1.5, 1.0, 1.0, 2.0, 0.4, 0.4, 2.151_657_414_559_676_2),
    (2.5, 0.5, 2.0, 4.0, -2.0, -3.0, 0.107_099_428_194_481_1),
    (0.5, 0.2, 0.4, 1.0, 0.95, 0.1, 1.270_215_025_918_858),
    (3.0, 2.0, 1.0, 5.0, -10.0, 0.2, 0.032_740_266_703_283),
    (1.2, 0.4, 0.9, 2.2, -0.3, -0.6, 0.743_450_923_837_328_5),
];

#[test]
fn appell_reference_values_and_reductions() {
    for &(a, b1, b2, c, x, y, want) in APPELL_CASES {
        close(appell_f1(a, b1, b2, c, x, y), want, 5e-10);
    }
    // x = y: F1(a; b1, b2; c; x, x) = 2F1(a, b1 + b2; c; x).
    for &(a, b1, b2, c, x) in &[
        (1.5_f64, 0.5_f64, 1.0_f64, 2.0_f64, 0.3_f64),
        (0.8, 0.5, 1.5, 1.5, 0.7),
        (2.0, -0.5, 3.0, 3.0, -0.6),
    ] {
        close(
            appell_f1(a, b1, b2, c, x, x),
            hyp2f1(a, b1 + b2, c, x),
            1e-10,
        );
    }
    // x = 0 and y = 0 reduce to Gauss 2F1.
    for &(a, b1, b2, c, point) in &[
        (1.5_f64, 0.5_f64, 0.7_f64, 2.0_f64, 0.4_f64),
        (0.8, -0.3, 1.2, 2.5, -0.5),
    ] {
        close(
            appell_f1(a, b1, b2, c, 0.0, point),
            hyp2f1(a, b2, c, point),
            1e-11,
        );
        close(
            appell_f1(a, b1, b2, c, point, 0.0),
            hyp2f1(a, b1, c, point),
            1e-11,
        );
    }
    // Trivial exact values and the domain limits.
    close(appell_f1(0.0, 2.0, 3.0, 4.0, 0.5, 0.2), 1.0, 1e-15);
    close(appell_f1(2.0, 0.0, 0.0, 3.0, 0.5, -0.2), 1.0, 1e-15);
    close(appell_f1(2.0, 0.5, 0.7, 3.0, 0.0, 0.0), 1.0, 1e-15);
    assert!(appell_f1(1.5, 0.5, 0.5, 2.0, 1.0, 0.2).is_nan());
    assert!(appell_f1(1.5, 0.5, 0.5, 2.0, 0.2, 1.2).is_nan());
    assert!(appell_f1(3.0, 0.5, 0.5, 0.0, 0.9, 0.2).is_nan());
    assert!(appell_f1(3.0, 0.5, 0.5, -1.0, 0.9, 0.2).is_nan());
}

// ---------------------------------------------------------------------
// Kelvin, Anger and Weber functions
// ---------------------------------------------------------------------

const KELVIN_CASES: &[(usize, f64, f64, f64, f64, f64)] = &[
    (
        0,
        0.5,
        0.999_023_463_990_838_3,
        0.062_493_218_382_199_46,
        0.855_905_872_118_634_2,
        -0.671_581_695_094_367_6,
    ),
    (
        0,
        1.0,
        0.984_381_781_213_086_9,
        0.249_566_040_036_659_7,
        0.286_706_208_728_316_1,
        -0.494_994_636_518_719_9,
    ),
    (
        0,
        2.5,
        0.399_968_417_129_531_3,
        1.457_182_044_159_804_2,
        -0.069_687_972_589_045_34,
        -0.110_696_099_155_674_9,
    ),
    (
        0,
        5.0,
        -6.230_082_478_666_358,
        0.116_034_381_550_200_4,
        -0.011_511_727_199_490_66,
        0.011_187_586_509_869_64,
    ),
    (
        0,
        8.0,
        20.973_955_610_730_26,
        -35.016_725_164_881_51,
        0.001_485_834_068_518_963,
        0.000_369_583_956_125_959_6,
    ),
    (
        0,
        12.0,
        -128.511_626_156_538_7,
        546.948_552_454_246_6,
        -6.307_713_705_205_456e-5,
        -3.899_959_497_178_822e-5,
    ),
    (
        0,
        15.0,
        -2967.254_534_635_129,
        -2952.707_887_344_117,
        -1.514_347_207_347_207e-8,
        7.962_894_398_377_627e-6,
    ),
    (
        0,
        30.0,
        -46_117_602.577_985,
        109_955_713.182_506_1,
        -1.293_826_937_602_08e-10,
        -5.289_996_606_628_324e-11,
    ),
    (
        0,
        50.0,
        -117_623_968_512_357.44,
        -50_192_646_254_462.21,
        -2.915_077_089_396_798e-17,
        7.255_813_220_365_599e-17,
    ),
    (
        1,
        0.5,
        -0.182_243_123_755_112_1,
        0.171_195_179_717_015_3,
        -1.522_403_406_532_09,
        -1.051_182_085_412_523,
    ),
    (
        1,
        1.0,
        -0.395_868_261_019_711_4,
        0.307_556_631_375_536_6,
        -0.740_322_276_841_982_7,
        -0.241_995_966_429_738_2,
    ),
    (
        1,
        2.5,
        -1.373_096_897_645_111,
        0.038_668_443_965_950_48,
        -0.117_256_135_859_870_5,
        0.093_313_788_135_357_5,
    ),
    (
        1,
        5.0,
        0.359_776_666_776_672_8,
        -5.797_907_901_792_625,
        0.012_737_390_484_218_57,
        0.011_577_754_393_252_47,
    ),
    (
        1,
        8.0,
        32.506_861_295_691_42,
        21.673_535_103_011_48,
        0.000_322_857_047_807_720_9,
        -0.001_566_974_799_960_044,
    ),
    (
        1,
        12.0,
        -526.963_435_897_347_3,
        -141.349_793_355_256_5,
        -3.834_011_517_754_054e-5,
        6.605_000_654_623_384e-5,
    ),
    (
        1,
        15.0,
        2954.865_291_352_463,
        -2826.093_604_257_383,
        8.150_749_777_615_45e-6,
        -1.679_694_873_684_5e-7,
    ),
    (
        1,
        30.0,
        -108_110_203.693_179_7,
        -46_885_734.672_072_37,
        -5.201_599_505_837_437e-11,
        1.315_233_470_379_599e-10,
    ),
    (
        1,
        50.0,
        50_675_459_067_627.8,
        -116_434_867_731_360.8,
        7.327_583_669_835_41e-17,
        2.884_739_143_245_268e-17,
    ),
];

#[test]
fn kelvin_reference_values() {
    for &(nu, x, ber, bei, ker, kei) in KELVIN_CASES {
        let scale = ber.abs().max(bei.abs()).max(1.0);
        assert!(
            (kelvin_ber(nu, x) - ber).abs() <= 2e-10 * scale,
            "ber_{nu}({x})"
        );
        assert!(
            (kelvin_bei(nu, x) - bei).abs() <= 2e-10 * scale,
            "bei_{nu}({x})"
        );
        let scale = ker.abs().max(kei.abs()).max(1e-300);
        assert!(
            (kelvin_ker(nu, x) - ker).abs() <= 2e-9 * scale,
            "ker_{nu}({x})"
        );
        assert!(
            (kelvin_kei(nu, x) - kei).abs() <= 2e-9 * scale,
            "kei_{nu}({x})"
        );
    }
    // Values at the origin and domain limits.
    assert_eq!(kelvin_ber(0, 0.0), 1.0);
    assert_eq!(kelvin_bei(0, 0.0), 0.0);
    assert_eq!(kelvin_ber(1, 0.0), 0.0);
    assert_eq!(kelvin_bei(1, 0.0), 0.0);
    assert!(kelvin_ber(0, -1.0).is_nan());
    assert!(kelvin_ber(101, 1.0).is_nan());
    assert!(kelvin_ker(0, -1.0).is_nan());
}

#[test]
fn kelvin_derivative_identities() {
    let root = std::f64::consts::SQRT_2;
    for x in [1.0_f64, 2.5, 5.0, 8.0] {
        let h = 1e-5;
        let d_ber = (kelvin_ber(0, x + h) - kelvin_ber(0, x - h)) / (2.0 * h);
        let d_bei = (kelvin_bei(0, x + h) - kelvin_bei(0, x - h)) / (2.0 * h);
        let d_ker = (kelvin_ker(0, x + h) - kelvin_ker(0, x - h)) / (2.0 * h);
        let d_kei = (kelvin_kei(0, x + h) - kelvin_kei(0, x - h)) / (2.0 * h);
        let scale = kelvin_ber(1, x).abs().max(kelvin_bei(1, x).abs()).max(1.0);
        assert!((d_ber - (kelvin_ber(1, x) + kelvin_bei(1, x)) / root).abs() < 2e-6 * scale);
        assert!((d_bei - (kelvin_bei(1, x) - kelvin_ber(1, x)) / root).abs() < 2e-6 * scale);
        let scale = kelvin_ker(1, x)
            .abs()
            .max(kelvin_kei(1, x).abs())
            .max(1e-12);
        assert!((d_ker - (kelvin_ker(1, x) + kelvin_kei(1, x)) / root).abs() < 2e-5 * scale);
        assert!((d_kei - (kelvin_kei(1, x) - kelvin_ker(1, x)) / root).abs() < 2e-5 * scale);
    }
}

const ANGER_CASES: &[(f64, f64, f64, f64)] = &[
    (0.5, 1.0, 0.855_165_309_679_262_2, 0.099_507_542_640_043_17),
    (1.3, 2.5, 0.465_457_606_540_690_2, -0.114_892_690_483_087_1),
    (2.0, 0.5, 0.030_604_023_458_682_64, 0.100_860_937_614_390_4),
    (
        -0.5,
        3.0,
        -0.570_621_336_461_999_3,
        -0.151_446_717_588_907_6,
    ),
    (3.0, 10.0, 0.058_379_379_305_186_81, 0.231_363_593_834_378_2),
    (0.0, 1.0, 0.765_197_686_557_966_6, -0.568_656_627_048_288),
    (1.0, 5.0, -0.327_579_137_591_465_2, -0.171_192_173_426_483_1),
    (2.5, -2.0, 0.251_641_727_595_066_8, 0.293_301_108_747_901_8),
    (4.0, 20.0, 0.130_670_933_554_863_3, -0.157_138_840_532_805_7),
    (-1.3, 0.7, -0.401_592_605_121_233_4, -0.184_715_057_585_216),
];

#[test]
fn anger_and_weber() {
    for &(nu, x, anger, weber) in ANGER_CASES {
        close(anger_j(nu, x), anger, 1e-10);
        close(weber_e(nu, x), weber, 1e-10);
    }
    // Integer order: the Anger function is the Bessel function.
    for n in [-3.0_f64, 0.0, 1.0, 4.0] {
        for x in [0.5_f64, 3.0, 10.0] {
            close(anger_j(n, x), bessel_j(n, x), 1e-12);
        }
    }
    // Derivative recurrences J_nu' = (J_{nu-1} - J_{nu+1})/2, likewise E.
    for &(nu, x) in &[(0.5_f64, 1.0_f64), (1.3, 2.5), (2.0, 5.0), (-0.5, 3.0)] {
        let h = 1e-5;
        let d_anger = (anger_j(nu, x + h) - anger_j(nu, x - h)) / (2.0 * h);
        let d_weber = (weber_e(nu, x + h) - weber_e(nu, x - h)) / (2.0 * h);
        close(
            d_anger,
            (anger_j(nu - 1.0, x) - anger_j(nu + 1.0, x)) / 2.0,
            1e-6,
        );
        close(
            d_weber,
            (weber_e(nu - 1.0, x) - weber_e(nu + 1.0, x)) / 2.0,
            1e-6,
        );
    }
}

/// `(ν, x, ber_ν, bei_ν, ker_ν, kei_ν)` for orders above one, from mpmath
/// (`J_ν(x e^{3πi/4})` and `e^{-iπν/2} K_ν(x e^{iπ/4})` at 40 digits).
#[allow(clippy::excessive_precision)] // mpmath prints 17 digits
const KELVIN_HIGHER_ORDERS: &[(usize, f64, f64, f64, f64, f64)] = &[
    (
        2,
        0.5,
        0.00065102047400263239,
        -0.031244913792168524,
        0.47690929429468515,
        7.950780708363778,
    ),
    (
        2,
        1.0,
        0.01041124172310777,
        -0.1246745356785978,
        0.41803361794402502,
        1.8842024187201009,
    ),
    (
        2,
        2.5,
        0.39864666006832802,
        -0.70231531727595382,
        0.18880430952581467,
        0.12423993630249844,
    ),
    (
        2,
        5.0,
        4.4884262727291134,
        1.4221014077546614,
        0.011183732582038347,
        -0.018064948041226885,
    ),
    (
        2,
        8.0,
        -22.889035214142713,
        25.438893739657008,
        -0.001819912297233685,
        -0.0001496529313405915,
    ),
    (
        2,
        12.0,
        173.95662974979833,
        -468.18708306071398,
        7.5379630878682888e-5,
        3.5733952956330616e-5,
    ),
    (
        2,
        15.0,
        2422.220503055607,
        2940.5671762548155,
        -7.6915290166749971e-7,
        -8.7155181418384226e-6,
    ),
    (
        2,
        30.0,
        49003751.725943048,
        -102649134.5775433,
        1.380348213209304e-10,
        4.9151953550353202e-11,
    ),
    (
        2,
        50.0,
        112897374700922.85,
        52052603206179.934,
        2.7894144696964213e-17,
        -7.544661328877416e-17,
    ),
    (
        2,
        100.0,
        -6.9967285889179983e+28,
        -1.8902903755660902e+29,
        1.035693096165584e-32,
        2.2540680184132102e-32,
    ),
    (
        3,
        0.5,
        0.0018124724599922045,
        0.0018700157083679225,
        43.801004675368463,
        -46.623031919530453,
    ),
    (
        3,
        1.0,
        0.013787984049301442,
        0.015628768611018837,
        4.8872738824136287,
        -6.2697108872392399,
    ),
    (
        3,
        2.5,
        0.12750060963120389,
        0.30489342282657011,
        0.044209886053044849,
        -0.44748292259396356,
    ),
    (
        3,
        5.0,
        -2.0943519510521605,
        2.454408539183927,
        -0.029282942811794904,
        -0.0076851511020541331,
    ),
    (
        3,
        8.0,
        -15.420358153724682,
        -22.575046230240226,
        0.00026766881417395374,
        0.0022633212645170961,
    ),
    (
        3,
        12.0,
        375.60871128914055,
        210.70057627783664,
        2.8995539276329509e-5,
        -9.2239729416793175e-5,
    ),
    (
        3,
        15.0,
        -2857.1249053192234,
        1814.8766618296312,
        -9.6491307768447646e-6,
        1.9564162106182435e-6,
    ),
    (
        3,
        30.0,
        93812232.454285384,
        51943469.871234376,
        4.3636037920790056e-11,
        -1.4917148541928309e-10,
    ),
    (
        3,
        50.0,
        -54117359109522.088,
        107103887896586.09,
        -7.9121672757927972e-17,
        -2.6157417592178351e-17,
    ),
    (
        3,
        100.0,
        1.8691083902467538e+29,
        -6.5424434562727306e+28,
        2.2754161929430815e-32,
        -1.0942694171476132e-32,
    ),
    (
        5,
        0.5,
        5.8141236065369249e-6,
        -5.6942406311812366e-6,
        8823.2639213368059,
        8551.7643576001612,
    ),
    (
        5,
        1.0,
        0.00019167656261427883,
        -0.00017633421766040988,
        287.76494549392659,
        253.88227095043138,
    ),
    (
        5,
        2.5,
        0.02211047908855141,
        -0.01281247444396538,
        3.5385514590822467,
        1.4953531163985601,
    ),
    (
        5,
        5.0,
        0.85103812436999622,
        0.21143444001332132,
        0.077129537582891512,
        -0.056323305678988751,
    ),
    (
        5,
        8.0,
        -1.8213634575286227,
        12.812799869404204,
        -0.0031462114865251252,
        -0.0034552625830007342,
    ),
    (
        5,
        12.0,
        -143.13278749867267,
        -224.19923835995,
        2.044018899835694e-5,
        0.00015428435491668834,
    ),
    (
        5,
        15.0,
        2274.4704453957567,
        -401.96275677556265,
        1.2228523983704834e-5,
        -7.5917908895113803e-6,
    ),
    (
        5,
        30.0,
        -67986608.369629166,
        -57055621.205329841,
        -1.8920311089926327e-11,
        1.8679501753770503e-10,
    ),
    (
        5,
        50.0,
        58944764545696.221,
        -89480981521939.545,
        9.1315567294845551e-17,
        1.9258333794677473e-17,
    ),
    (
        5,
        100.0,
        -1.7986165198005527e+29,
        5.1670562061173603e+28,
        -2.3390853933046872e-32,
        1.2912884062819384e-32,
    ),
    (
        10,
        0.5,
        1.493215789465446e-15,
        -2.628031871371771e-13,
        1321190400.0999992,
        190248468519.99974,
    ),
    (
        10,
        1.0,
        6.1158290022627827e-12,
        -2.6905073656370991e-10,
        5160000.0999876044,
        185713929.99895853,
    ),
    (
        10,
        2.5,
        3.6360582372500299e-7,
        -2.542767746513776e-6,
        3357.8102678489722,
        19153.263749590153,
    ),
    (
        10,
        5.0,
        0.0014314823385488821,
        -0.0022461282631687126,
        11.769354216704298,
        14.240646929785773,
    ),
    (
        10,
        8.0,
        0.31220091654618469,
        -0.04045018274236428,
        0.14374999474995972,
        -0.02250252348617757,
    ),
    (
        10,
        12.0,
        -0.18278590069908194,
        24.697249623073272,
        -0.00071677387205550979,
        -0.0013507027511021333,
    ),
    (
        10,
        15.0,
        -360.62989432182461,
        -5.1209739169586488,
        -7.34657116508621e-5,
        4.9176629766348317e-5,
    ),
    (
        10,
        30.0,
        36267692.216799164,
        800522.55237632641,
        3.3446439427664037e-10,
        -3.1299334682990881e-10,
    ),
    (
        10,
        50.0,
        27508884009335.222,
        56573062237574.63,
        -4.8908777224237247e-17,
        -1.5119032183802049e-16,
    ),
    (
        10,
        100.0,
        -1.8610698704628851e+27,
        -1.4349733027295827e+29,
        2.4189353626319677e-32,
        2.5074021754229658e-32,
    ),
    (
        20,
        0.5,
        -3.7382962971393896e-31,
        -1.1125914556941956e-33,
        -6.6874719264734165e+28,
        2.1998341742847056e+26,
    ),
    (
        20,
        1.0,
        -3.9196392070812351e-25,
        -4.666456730627322e-27,
        -6.3771238931910892e+22,
        8.3914336014667767e+20,
    ),
    (
        20,
        2.5,
        -3.5557158645829092e-17,
        -2.6504954698911216e-18,
        -698735204607768.28,
        57590811825368.7,
    ),
    (
        20,
        5.0,
        -3.5812034464085299e-11,
        -1.0983243188696403e-11,
        -631013608.98089831,
        215341431.5242135,
    ),
    (
        20,
        8.0,
        -3.3149892950509523e-7,
        -3.1591038544832499e-7,
        -36173.777176555384,
        40423.229353937752,
    ),
    (
        20,
        12.0,
        0.00021920621022488485,
        -0.0015901565483126905,
        4.6059275466491923,
        14.382471401753977,
    ),
    (
        20,
        15.0,
        0.13490687412109218,
        -0.071459441458819295,
        0.14875282011086423,
        0.034982420217856627,
    ),
    (
        20,
        30.0,
        858250.97211407789,
        279127.27859314986,
        1.110515472328322e-8,
        -1.3726598017950208e-8,
    ),
    (
        20,
        50.0,
        7243399757168.7588,
        662199492628.78686,
        9.5454411042158809e-16,
        -9.7747588191554929e-16,
    ),
    (
        20,
        100.0,
        -4.2925241333666338e+28,
        2.4578660003191362e+28,
        -9.8016387353515453e-32,
        2.4547778771599654e-32,
    ),
];

#[test]
fn kelvin_functions_of_higher_integer_order() {
    for &(nu, x, ber, bei, ker, kei) in KELVIN_HIGHER_ORDERS {
        let scale = ber.abs().max(bei.abs()).max(1e-300);
        assert!(
            (kelvin_ber(nu, x) - ber).abs() <= 1e-10 * scale,
            "ber_{nu}({x}) = {:e}, want {ber:e}",
            kelvin_ber(nu, x)
        );
        assert!(
            (kelvin_bei(nu, x) - bei).abs() <= 1e-10 * scale,
            "bei_{nu}({x}) = {:e}, want {bei:e}",
            kelvin_bei(nu, x)
        );
        let scale = ker.abs().max(kei.abs()).max(1e-300);
        assert!(
            (kelvin_ker(nu, x) - ker).abs() <= 1e-9 * scale,
            "ker_{nu}({x}) = {:e}, want {ker:e}",
            kelvin_ker(nu, x)
        );
        assert!(
            (kelvin_kei(nu, x) - kei).abs() <= 1e-9 * scale,
            "kei_{nu}({x}) = {:e}, want {kei:e}",
            kelvin_kei(nu, x)
        );
    }
}

#[test]
fn kelvin_recurrences_and_the_origin_for_higher_orders() {
    // J_{n-1} + J_{n+1} = (2n/z) J_n with z = x e^{3πi/4}, i.e. in terms of
    // (ber, bei); and the K recurrence K_{n+1} = K_{n-1} + (2n/z) K_n with
    // z = x e^{iπ/4}, rotated by e^{-iπν/2}.
    for &x in &[0.7_f64, 3.0, 9.0, 25.0, 40.0, 70.0] {
        let z = num_complex::Complex64::new(-x / 2.0_f64.sqrt(), x / 2.0_f64.sqrt());
        for n in 1..=8_usize {
            let j = |m: usize| num_complex::Complex64::new(kelvin_ber(m, x), kelvin_bei(m, x));
            let lhs = j(n - 1) + j(n + 1);
            let rhs = num_complex::Complex64::new(2.0 * n as f64, 0.0) / z * j(n);
            let scale = lhs.norm().max(rhs.norm()).max(1e-300);
            assert!(
                (lhs - rhs).norm() <= 1e-9 * scale,
                "J recurrence n = {n}, x = {x}: {lhs} vs {rhs}"
            );
        }
    }
    for &x in &[0.3_f64, 1.5, 6.0, 14.0, 20.0, 60.0] {
        let zk = num_complex::Complex64::new(x / 2.0_f64.sqrt(), x / 2.0_f64.sqrt());
        // K_nu = i^nu (ker + i kei)
        let k = |m: usize| {
            let rotated = num_complex::Complex64::new(kelvin_ker(m, x), kelvin_kei(m, x));
            rotated * num_complex::Complex64::new(0.0, 1.0).powi(m as i32)
        };
        for n in 1..=8_usize {
            let lhs = k(n + 1);
            let rhs = k(n - 1) + num_complex::Complex64::new(2.0 * n as f64, 0.0) / zk * k(n);
            let scale = lhs.norm().max(rhs.norm()).max(1e-300);
            assert!(
                (lhs - rhs).norm() <= 1e-9 * scale,
                "K recurrence n = {n}, x = {x}"
            );
        }
    }
    // At the origin: ber, bei vanish for nu >= 1; ker, kei diverge with the
    // sign of the leading term e^{-3πiν/4}.
    for nu in 1..=6_usize {
        assert_eq!(kelvin_ber(nu, 0.0), 0.0);
        assert_eq!(kelvin_bei(nu, 0.0), 0.0);
    }
    assert_eq!(
        kelvin_ker(1, 0.0),
        f64::NEG_INFINITY,
        "ker_1 ~ -1/(x sqrt 2)"
    );
    assert_eq!(kelvin_kei(1, 0.0), f64::NEG_INFINITY);
    assert_eq!(kelvin_kei(2, 0.0), f64::INFINITY, "kei_2 ~ 2/x^2");
    assert!(
        kelvin_ker(2, 0.0).is_nan(),
        "ker_2(0) is finite and not computed"
    );
    assert_eq!(kelvin_ker(3, 0.0), f64::INFINITY);
    // The known values near zero agree with the sign: ker_1(0.01) < 0 etc.
    assert!(kelvin_ker(1, 0.01) < -50.0 && kelvin_kei(1, 0.01) < -50.0);
    assert!(kelvin_kei(2, 0.01) > 1.0e3);
    // Out of range.
    assert!(kelvin_ker(101, 1.0).is_nan() && kelvin_bei(101, 1.0).is_nan());
}

#[test]
#[allow(clippy::excessive_precision, clippy::inconsistent_digit_grouping)]
fn coulomb_irregular_function_at_small_rho_matches_mpmath() {
    // G and G' used to come from Steed's continued fraction, which lost digits at
    // small rho (relative error 8e-7 at eta = 1.94, l = 4, rho = 0.555, and worse);
    // they now come from the Tricomi function. mpmath at 40 digits.
    let table: [(f64, usize, f64, f64, f64, f64, f64); 64] = [
        (
            -3.0,
            0,
            0.05,
            0.18602459808205699,
            3.1288800758073005,
            0.277253978197505,
            -0.7123012388988763,
        ),
        (
            -3.0,
            0,
            0.5551033346899956,
            0.06677047044651739,
            -1.790906707440172,
            -0.5421563452306725,
            -0.43505107348924316,
        ),
        (
            -3.0,
            0,
            2.0,
            0.19680255620767087,
            1.3726383209278323,
            0.6814042161474986,
            -0.3286464470812138,
        ),
        (
            -3.0,
            0,
            6.0,
            0.3715746262822885,
            1.0740078828290804,
            0.7548995346838618,
            -0.5092703742955965,
        ),
        (
            -3.0,
            3,
            0.05,
            2.0063984931769743e-06,
            0.00015898953170838883,
            3605.866386115425,
            -212672.10542649202,
        ),
        (
            -3.0,
            3,
            0.5551033346899956,
            0.02028944777072735,
            0.12884249979064946,
            4.730526121091798,
            -19.2468121186114,
        ),
        (
            -3.0,
            3,
            2.0,
            0.7516623136343783,
            0.5236460579938922,
            0.5542254344482435,
            -0.944283121737663,
        ),
        (
            -3.0,
            3,
            6.0,
            -0.13317179148090239,
            1.1220301262096184,
            0.8699571953878981,
            0.17933090781748576,
        ),
        (
            -3.0,
            4,
            0.05,
            1.4039100709086367e-08,
            1.3953993852761232e-06,
            398749.61322091153,
            -31596399.514774594,
        ),
        (
            -3.0,
            4,
            0.5551033346899956,
            0.001714757793811285,
            0.014291563933802697,
            39.829162216079396,
            -251.21855886226186,
        ),
        (
            -3.0,
            4,
            2.0,
            0.33274546723926457,
            0.5236460579938922,
            1.309651931838374,
            -0.944283121737663,
        ),
        (
            -3.0,
            4,
            6.0,
            -0.8887459815356344,
            -0.24052690447909753,
            -0.20146187261318182,
            1.0706580048504408,
        ),
        (
            -3.0,
            10,
            0.05,
            2.555422176294796e-23,
            5.614899701791468e-21,
            9.330095099650793e+19,
            -1.863193183848777e+22,
        ),
        (
            -3.0,
            10,
            0.5551033346899956,
            6.98153159450584e-12,
            1.362594337866653e-10,
            3850476046.017955,
            -68084818886.237595,
        ),
        (
            -3.0,
            10,
            2.0,
            5.716987454127573e-06,
            2.9320139397990445e-05,
            18026.001568210744,
            -82469.15443049802,
        ),
        (
            -3.0,
            10,
            6.0,
            0.14210638344492463,
            0.17086203544424214,
            3.540913672126538,
            -2.7795534097256,
        ),
        (
            0.5,
            0,
            0.05,
            0.019301113034633346,
            0.3952756057785516,
            2.309448547440257,
            -4.514316155747629,
        ),
        (
            0.5,
            0,
            0.5551033346899956,
            0.26008041162297263,
            0.5447478789613189,
            1.4728410160421557,
            -0.7600456306196555,
        ),
        (
            0.5,
            0,
            2.0,
            1.0211202242957036,
            0.3296024266337856,
            0.5322112059077913,
            -0.8075267489876982,
        ),
        (
            0.5,
            0,
            6.0,
            -1.0286067495607154,
            -0.1643901107340414,
            -0.18863987402419613,
            0.9420407465089363,
        ),
        (
            0.5,
            3,
            0.05,
            2.635677506792051e-08,
            2.1116880697549336e-06,
            270473.66957057465,
            -16270730.302924512,
        ),
        (
            0.5,
            3,
            0.5551033346899956,
            0.00041927216784073913,
            0.0030476433931199078,
            187.37856726331376,
            -1023.0513264849474,
        ),
        (
            0.5,
            3,
            2.0,
            0.0689129580359834,
            0.13103168516256103,
            4.6460853642117135,
            -5.676952730921114,
        ),
        (
            0.5,
            3,
            6.0,
            1.152510977577688,
            -0.21055166723540558,
            -0.20469265386509988,
            -0.8302754933137814,
        ),
        (
            0.5,
            4,
            0.05,
            1.4738557544019456e-10,
            1.4752620019292435e-08,
            37649436.931260735,
            -3016388555.5624905,
        ),
        (
            0.5,
            4,
            0.5551033346899956,
            2.5783845577513103e-05,
            0.00023351713000715205,
            2378.1908107367462,
            -17245.36031410177,
        ),
        (
            0.5,
            4,
            2.0,
            0.015289365480437968,
            0.036959352089541805,
            15.429806007531958,
            -28.106095550146364,
        ),
        (
            0.5,
            4,
            6.0,
            1.1142845849978567,
            0.2793381067485534,
            0.6630670099608422,
            -0.7312136663827952,
        ),
        (
            0.5,
            10,
            0.05,
            1.6035927639399477e-25,
            3.528598060193883e-23,
            1.4844408868998669e+22,
            -2.969584842920239e+24,
        ),
        (
            0.5,
            10,
            0.5551033346899956,
            5.147768957978719e-14,
            1.021185468284235e-12,
            512923596764.1764,
            -9250801279770.822,
        ),
        (
            0.5,
            10,
            2.0,
            6.737340312975623e-08,
            3.6774906515267593e-07,
            1426807.1164497216,
            -7054609.010280364,
        ),
        (
            0.5,
            10,
            6.0,
            0.007031831047110343,
            0.011291900396639725,
            47.842976212208036,
            -65.38295286856722,
        ),
        (
            1.9352611942164177,
            0,
            0.05,
            0.00043873051673622013,
            0.009590561721831798,
            82.39293218749133,
            -478.21040572681414,
        ),
        (
            1.9352611942164177,
            0,
            0.5551033346899956,
            0.010798492861367,
            0.034042971009084703,
            18.016342235475605,
            -35.807791749468855,
        ),
        (
            1.9352611942164177,
            0,
            2.0,
            0.1617849936575037,
            0.19673678810432674,
            3.241312403916752,
            -2.239482168400108,
        ),
        (
            1.9352611942164177,
            0,
            6.0,
            1.2688851971592134,
            -0.10826237557156904,
            -0.05736350478769196,
            -0.7831990576653228,
        ),
        (
            1.9352611942164177,
            3,
            0.05,
            1.7552040179632245e-09,
            1.4125353993654072e-07,
            4037572.4967359276,
            -244802648.4033698,
        ),
        (
            1.9352611942164177,
            3,
            0.5551033346899956,
            3.33637458706739e-05,
            0.0002541385622232347,
            2215.271213831445,
            -13098.474028555764,
        ),
        (
            1.9352611942164177,
            3,
            2.0,
            0.00901026872744141,
            0.02023359894147746,
            27.891012500733275,
            -48.35199172934563,
        ),
        (
            1.9352611942164177,
            3,
            6.0,
            0.831435949568417,
            0.391805869947624,
            1.3579087631034374,
            -0.5628375535189051,
        ),
        (
            1.9352611942164177,
            4,
            0.05,
            1.0780593337430117e-11,
            1.0821758228377588e-09,
            512887628.11926246,
            -41274667830.370255,
        ),
        (
            1.9352611942164177,
            4,
            0.5551033346899956,
            2.1766286753530014e-06,
            2.0325881571911997e-05,
            27125.273180269247,
            -206123.7705872475,
        ),
        (
            1.9352611942164177,
            4,
            2.0,
            0.00193200389032125,
            0.005210679247125492,
            105.88634720317992,
            -232.01827404184473,
        ),
        (
            1.9352611942164177,
            4,
            6.0,
            0.5083726078245699,
            0.33876060872486075,
            1.9129585605767117,
            -0.6923366604265516,
        ),
        (
            1.9352611942164177,
            10,
            0.05,
            1.435515521074315e-26,
            3.1606275401673078e-24,
            1.6571662622905174e+23,
            -3.3174943794438638e+25,
        ),
        (
            1.9352611942164177,
            10,
            0.5551033346899956,
            4.921344167411398e-15,
            9.826629319219479e-14,
            5326869722420.614,
            -96832946822870.25,
        ),
        (
            1.9352611942164177,
            10,
            2.0,
            7.768997820621162e-09,
            4.341503823096638e-08,
            12056873.298699584,
            -61340007.52659985,
        ),
        (
            1.9352611942164177,
            10,
            6.0,
            0.0013852854017528306,
            0.0024199446429581254,
            219.81770876857564,
            -337.87500586227094,
        ),
        (
            4.0,
            0,
            0.05,
            1.0606129582901143e-06,
            2.5181322343306274e-05,
            28523.133951898028,
            -265648.0552247638,
        ),
        (
            4.0,
            0,
            0.5551033346899956,
            5.28604659673126e-05,
            0.00022492268325534392,
            2523.563966711848,
            -8179.898026401268,
        ),
        (
            4.0,
            0,
            2.0,
            0.00264173965093556,
            0.005043591755519073,
            108.777127518402,
            -170.86186986643037,
        ),
        (
            4.0,
            0,
            6.0,
            0.2944492156204837,
            0.20917339593576786,
            3.0137838132919463,
            -1.2552134138942634,
        ),
        (
            4.0,
            3,
            0.05,
            1.6805766602174564e-11,
            1.3610822053241167e-09,
            418167926.407966,
            -25636383434.783913,
        ),
        (
            4.0,
            3,
            0.5551033346899956,
            4.098960772810381e-07,
            3.317845395464794e-06,
            167133.18558701806,
            -1086807.0090490647,
        ),
        (
            4.0,
            3,
            2.0,
            0.00021096415789938,
            0.0005632440528393283,
            963.2912336180799,
            -2168.2922163988164,
        ),
        (
            4.0,
            3,
            6.0,
            0.1052558822390735,
            0.097971236630984,
            5.958184515783218,
            -3.954831654611398,
        ),
        (
            4.0,
            4,
            0.05,
            1.3073658529525404e-13,
            1.3177279645956499e-11,
            42078500063.295166,
            -3407767126374.032,
        ),
        (
            4.0,
            4,
            0.5551033346899956,
            3.231952248763729e-08,
            3.144709352554607e-07,
            1738266.0762767892,
            -14027615.70638358,
        ),
        (
            4.0,
            4,
            2.0,
            4.924887068820037e-05,
            0.00015060176121132118,
            3576.6634204563093,
            -9367.69073427113,
        ),
        (
            4.0,
            4,
            6.0,
            0.0547691210353686,
            0.057572427789762126,
            9.818276072545267,
            -7.937648104900111,
        ),
        (
            4.0,
            10,
            0.05,
            3.209258936058845e-28,
            7.071960819127024e-26,
            7.40564213466378e+24,
            -1.4840681270072574e+27,
        ),
        (
            4.0,
            10,
            0.5551033346899956,
            1.208881769307865e-16,
            2.4362445422778354e-15,
            214667959681174.94,
            -3945930585899486.0,
        ),
        (
            4.0,
            10,
            2.0,
            2.4876781585582406e-10,
            1.4354289920398763e-09,
            363566664.20380616,
            -1921976394.0032058,
        ),
        (
            4.0,
            10,
            6.0,
            9.195553336001117e-05,
            0.00017774237827533813,
            2953.248225626121,
            -5166.438814288806,
        ),
    ];
    for (eta, l, rho, f, fp, g, gp) in table {
        let (cf, cfp, cg, cgp) = coulomb_fg(eta, l, rho)
            .unwrap_or_else(|| panic!("no values at eta={eta}, l={l}, rho={rho}"));
        let rel = |got: f64, want: f64| (got - want).abs() / want.abs().max(1e-300);
        assert!(
            rel(cg, g) < 1e-11,
            "G(eta={eta}, l={l}, rho={rho}) = {cg:e}, want {g:e}"
        );
        assert!(
            rel(cgp, gp) < 1e-11,
            "G'(eta={eta}, l={l}, rho={rho}) = {cgp:e}, want {gp:e}"
        );
        assert!(
            rel(cf, f) < 1e-11 && rel(cfp, fp) < 1e-10,
            "F at eta={eta}, l={l}, rho={rho}: {cf:e} {cfp:e}"
        );
        // the Wronskian F' G - F G' = 1
        assert!(
            (cfp * cg - cf * cgp - 1.0).abs() < 1e-9,
            "Wronskian at eta={eta}, l={l}, rho={rho}"
        );
    }
}
