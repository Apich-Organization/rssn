//! Combinatorics: factorials, binomial coefficients, counting sequences and
//! the identities among them.
//!
//! Operators are reduced by exact kernels on literal integers. Results grow
//! super-exponentially, so each kernel declines inputs whose answer would
//! run to tens of thousands of bits; the term then stays symbolic.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use statrs::function::gamma::ln_gamma;

use crate::graph::Arity;
use crate::graph::OpDescriptor;
use crate::graph::RuleError;
use crate::graph::RuleSet;
use crate::graph::Tier;
use crate::graph::rule::Installer;

use super::arith::arith;
use super::number_theory::Value;
use super::number_theory::exact;
use super::number_theory::exact_lists;

/// Largest supported argument of `factorial` (a 54,000-bit result).
const MAX_FACTORIAL: u64 = 5000;
/// Bound on the bits a multiplicative kernel may accumulate.
const MAX_BITS: u64 = 60_000;

/// The combinatorics rule set.
#[must_use]
pub fn combinatorics() -> RuleSet {
    RuleSet::new("combinatorics", install).needs(arith())
}

/// The value of `v` if it is a natural number not above `cap`.
fn nat(
    v: &BigInt,
    cap: u64,
) -> Option<u64> {
    v.to_u64().filter(|&x| x <= cap)
}

fn factorial(n: u64) -> BigInt {
    (2..=n).fold(BigInt::one(), |acc, k| acc * k)
}

fn double_factorial(n: u64) -> BigInt {
    (1..=n)
        .rev()
        .step_by(2)
        .fold(BigInt::one(), |acc, k| acc * k)
}

/// `binomial(n, k)` for `0 <= k <= n`, by the multiplicative formula.
fn choose(
    n: &BigInt,
    k: u64,
) -> Option<BigInt> {
    let k = n.to_u64().map_or(k, |n| k.min(n.saturating_sub(k)));
    if k.saturating_mul(n.bits()) > MAX_BITS {
        return None;
    }
    let base = n - k;
    let mut acc = BigInt::one();
    for i in 1..=k {
        // Each partial product is itself a binomial coefficient, so the
        // division is exact.
        acc = acc * (&base + i) / i;
    }
    Some(acc)
}

/// Generalised `binomial(n, k)` for every integer `n`; zero for `k < 0`.
fn binomial(
    n: &BigInt,
    k: &BigInt,
) -> Option<BigInt> {
    if k.is_negative() {
        return Some(BigInt::zero());
    }
    if n.is_negative() {
        // (-m choose k) = (-1)^k (m + k - 1 choose k)
        let k = nat(k, MAX_FACTORIAL)?;
        let c = choose(&(BigInt::from(k) - n - 1), k)?;
        return Some(if k % 2 == 0 { c } else { -c });
    }
    if k > n {
        return Some(BigInt::zero());
    }
    choose(n, k.to_u64()?)
}

/// `x (x - step) ... ` over `k` factors.
fn product_run(
    x: &BigInt,
    k: &BigInt,
    step: i8,
) -> Option<BigInt> {
    let k = nat(k, MAX_FACTORIAL)?;
    if k.saturating_mul(
        x.bits()
            .max(1)
            .saturating_add(u64::from(64 - k.leading_zeros())),
    ) > MAX_BITS
    {
        return None;
    }
    let mut acc = BigInt::one();
    let mut term = x.clone();
    for _ in 0..k {
        acc *= &term;
        term += i32::from(step);
    }
    Some(acc)
}

/// `(F(n), F(n + 1))` by fast doubling.
fn fibonacci_pair(n: u64) -> (BigInt, BigInt) {
    if n == 0 {
        return (BigInt::zero(), BigInt::one());
    }
    let (a, b) = fibonacci_pair(n / 2);
    let c = &a * (&b * 2 - &a);
    let d = &a * &a + &b * &b;
    if n.is_multiple_of(2) {
        (c, d)
    } else {
        (d.clone(), c + d)
    }
}

fn bell(n: u64) -> BigInt {
    // Bell's triangle: each row starts with the end of the previous one.
    let mut row = vec![BigInt::one()];
    for _ in 0..n {
        let mut next = vec![row.last().cloned().unwrap_or_default()];
        for above in &row {
            let value = next.last().cloned().unwrap_or_default() + above;
            next.push(value);
        }
        row = next;
    }
    row.into_iter().next().unwrap_or_default()
}

/// Row `n` of the triangle `t(n, j) = t(n-1, j-1) + weight(n, j) * t(n-1, j)`
/// started from `t(0, 0) = 1`.
fn stirling_row(
    n: u64,
    weight: fn(u64, u64) -> u64,
) -> Vec<BigInt> {
    let mut row = vec![BigInt::one()];
    for m in 1..=n {
        let mut next = vec![BigInt::zero(); row.len() + 1];
        for (j, value) in row.iter().enumerate() {
            let w = weight(m, u64::try_from(j).unwrap_or(0));
            // Contribution to t(m, j + 1) and to t(m, j).
            if let Some(slot) = next.get_mut(j + 1) {
                *slot += value;
            }
            if let Some(slot) = next.get_mut(j) {
                *slot += value * w;
            }
        }
        row = next;
    }
    row
}

fn partitions(n: usize) -> BigInt {
    // Euler's pentagonal number recurrence.
    let mut p = vec![BigInt::one()];
    for m in 1..=n {
        let mut total = BigInt::zero();
        let mut j = 1_usize;
        loop {
            let (g1, g2) = (j * (3 * j - 1) / 2, j * (3 * j + 1) / 2);
            if g1 > m {
                break;
            }
            let mut part = p[m - g1].clone();
            if g2 <= m {
                part += &p[m - g2];
            }
            if j % 2 == 1 {
                total += part;
            } else {
                total -= part;
            }
            j += 1;
        }
        p.push(total);
    }
    p.pop().unwrap_or_default()
}

fn derangements(n: u64) -> BigInt {
    let (mut prev, mut cur) = (BigInt::one(), BigInt::zero());
    if n == 0 {
        return prev;
    }
    for k in 2..=n {
        let next = (&prev + &cur) * (k - 1);
        prev = cur;
        cur = next;
    }
    cur
}

fn harmonic(n: u64) -> BigRational {
    (1..=n).fold(BigRational::zero(), |acc, k| {
        acc + BigRational::new(BigInt::one(), BigInt::from(k))
    })
}

fn multinomial(args: &[Vec<BigInt>]) -> Option<Value> {
    let [parts] = args else {
        return None;
    };
    let mut total = 0_u64;
    let mut acc = BigInt::one();
    for part in parts {
        let k = nat(part, MAX_FACTORIAL)?;
        total = total.checked_add(k).filter(|&t| t <= MAX_FACTORIAL)?;
        acc *= choose(&BigInt::from(total), k)?;
    }
    Some(Value::Int(acc))
}

/// `Gamma(n + 1)` on reals: the value `n!` extends to.
fn factorial_eval(args: &[f64]) -> f64 {
    match args.first() {
        | Some(&x) if x > -1.0 => ln_gamma(x + 1.0).exp(),
        | _ => f64::NAN,
    }
}

/// `Gamma(n + 1) / (Gamma(k + 1) Gamma(n - k + 1))`, zero where the last
/// factor has a pole, and the usual extension to negative integers `n`.
fn binomial_eval(args: &[f64]) -> f64 {
    let (Some(&n), Some(&k)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    if n.is_nan() || k.is_nan() {
        return f64::NAN;
    }
    if k < 0.0 {
        return 0.0;
    }
    if n < 0.0 {
        if n.fract() != 0.0 || k.fract() != 0.0 {
            return f64::NAN;
        }
        let sign = if k.rem_euclid(2.0) == 0.0 {
            1.0
        } else {
            -1.0
        };
        return sign * binomial_eval(&[k - n - 1.0, k]);
    }
    let m = n - k + 1.0;
    if m <= 0.0 && m.fract() == 0.0 {
        return 0.0;
    }
    if m < 0.0 {
        return f64::NAN;
    }
    (ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(m)).exp()
}

/// `Gamma(n + 1) / Gamma(n - k + 1)`, zero where the denominator has a pole.
fn permutations_eval(args: &[f64]) -> f64 {
    let (Some(&n), Some(&k)) = (args.first(), args.get(1)) else {
        return f64::NAN;
    };
    let m = n - k + 1.0;
    if n < 0.0 || m < 0.0 && m.fract() != 0.0 {
        return f64::NAN;
    }
    if m <= 0.0 {
        return 0.0;
    }
    (ln_gamma(n + 1.0) - ln_gamma(m)).exp()
}

fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    let unary = |name: &str| OpDescriptor::new(name, Arity::Fixed(1));
    let binary = |name: &str| OpDescriptor::new(name, Arity::Fixed(2));
    let set = "combinatorics";

    exact(
        i,
        set,
        unary("factorial").eval(factorial_eval),
        |a| match a {
            | [n] => Some(Value::Int(factorial(nat(n, MAX_FACTORIAL)?))),
            | _ => None,
        },
    )?;
    exact(i, set, unary("double_factorial"), |a| match a {
        | [n] => Some(Value::Int(double_factorial(nat(n, 2 * MAX_FACTORIAL)?))),
        | _ => None,
    })?;
    exact(
        i,
        set,
        binary("binomial").eval(binomial_eval),
        |a| match a {
            | [n, k] => binomial(n, k).map(Value::Int),
            | _ => None,
        },
    )?;
    exact_lists(i, set, unary("multinomial"), multinomial)?;
    exact(
        i,
        set,
        binary("permutations").eval(permutations_eval),
        |a| match a {
            | [n, k] if !n.is_negative() && !k.is_negative() => {
                if k > n {
                    Some(Value::Int(BigInt::zero()))
                } else {
                    product_run(n, k, -1).map(Value::Int)
                }
            },
            | _ => None,
        },
    )?;
    exact(i, set, binary("rising"), |a| match a {
        | [x, k] => product_run(x, k, 1).map(Value::Int),
        | _ => None,
    })?;
    exact(i, set, binary("falling"), |a| match a {
        | [x, k] => product_run(x, k, -1).map(Value::Int),
        | _ => None,
    })?;
    exact(i, set, unary("catalan"), |a| match a {
        | [n] => {
            let n = nat(n, MAX_FACTORIAL / 2)?;
            let c = choose(&BigInt::from(2 * n), n)?;
            Some(Value::Int(c / (n + 1)))
        },
        | _ => None,
    })?;
    exact(i, set, unary("fibonacci"), |a| match a {
        | [n] => Some(Value::Int(fibonacci_pair(nat(n, 100_000)?).0)),
        | _ => None,
    })?;
    exact(i, set, unary("lucas"), |a| match a {
        | [n] => {
            let (f, next) = fibonacci_pair(nat(n, 100_000)?);
            Some(Value::Int(next * 2 - f))
        },
        | _ => None,
    })?;
    exact(i, set, unary("bell"), |a| match a {
        | [n] => Some(Value::Int(bell(nat(n, 1000)?))),
        | _ => None,
    })?;
    exact(i, set, binary("stirling1"), |a| match a {
        | [n, k] => {
            let (n, k) = (nat(n, 1000)?, k.to_u64()?);
            let row = stirling_row(n, |m, _| m - 1);
            Some(Value::Int(
                row.get(usize::try_from(k).ok()?)
                    .cloned()
                    .unwrap_or_default(),
            ))
        },
        | _ => None,
    })?;
    exact(i, set, binary("stirling2"), |a| match a {
        | [n, k] => {
            let (n, k) = (nat(n, 1000)?, k.to_u64()?);
            let row = stirling_row(n, |_, j| j);
            Some(Value::Int(
                row.get(usize::try_from(k).ok()?)
                    .cloned()
                    .unwrap_or_default(),
            ))
        },
        | _ => None,
    })?;
    exact(i, set, unary("partitions"), |a| match a {
        | [n] => Some(Value::Int(partitions(usize::try_from(nat(n, 5000)?).ok()?))),
        | _ => None,
    })?;
    exact(i, set, unary("derangements"), |a| match a {
        | [n] => Some(Value::Int(derangements(nat(n, MAX_FACTORIAL)?))),
        | _ => None,
    })?;
    exact(i, set, unary("harmonic"), |a| match a {
        | [n] => Some(Value::Rat(harmonic(nat(n, 1000)?))),
        | _ => None,
    })?;

    // All rules hold for the continuous extensions used as numeric
    // semantics, hence for every real argument satisfying the guards.
    i.rewrites(
        Tier::Normalize,
        &[
            "combinatorics/binomial-0: binomial(?n, 0) => 1",
            "combinatorics/binomial-1: binomial(?n, 1) => ?n if nonnegative(?n)",
            "combinatorics/binomial-self: binomial(?n, ?n) => 1 if nonnegative(?n)",
            "combinatorics/permutations-0: permutations(?n, 0) => 1 if nonnegative(?n)",
            "combinatorics/permutations-1: permutations(?n, 1) => ?n if nonnegative(?n)",
            "combinatorics/permutations-self: permutations(?n, ?n) => factorial(?n) \
             if nonnegative(?n)",
            "combinatorics/binomial-permutations: binomial(?n, ?k) * factorial(?k) \
             => permutations(?n, ?k) if nonnegative(?n), nonnegative(?k)",
            "combinatorics/factorial-ratio: factorial(?n + 1) / factorial(?n) => ?n + 1 \
             if nonnegative(?n)",
        ],
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rules::testing::eval;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn s(src: &str) -> String {
        simplify(&[combinatorics()], src)
    }

    fn value(src: &str) -> BigInt {
        s(src)
            .parse()
            .unwrap_or_else(|_| panic!("`{src}` did not reduce to an integer"))
    }

    #[test]
    fn factorials() {
        assert_eq!(s("factorial(0)"), "1");
        assert_eq!(s("factorial(1)"), "1");
        assert_eq!(s("factorial(5)"), "120");
        assert_eq!(s("factorial(20)"), "2432902008176640000");
        assert_eq!(s("factorial(30)"), "265252859812191058636308480000000");
        assert_eq!(s("factorial(-1)"), "factorial(-1)");
        assert_eq!(s("factorial(1000000)"), "factorial(1000000)");
        assert_eq!(s("double_factorial(0)"), "1");
        assert_eq!(s("double_factorial(7)"), "105");
        assert_eq!(s("double_factorial(8)"), "384");
        assert_eq!(s("double_factorial(10)"), "3840");
        assert_eq!(s("double_factorial(-3)"), "double_factorial(-3)");
    }

    #[test]
    fn binomials() {
        assert_eq!(s("binomial(5, 2)"), "10");
        assert_eq!(s("binomial(5, 0)"), "1");
        assert_eq!(s("binomial(5, 5)"), "1");
        assert_eq!(s("binomial(0, 0)"), "1");
        assert_eq!(s("binomial(5, 7)"), "0");
        assert_eq!(s("binomial(5, -1)"), "0");
        assert_eq!(s("binomial(-3, 2)"), "6");
        assert_eq!(s("binomial(-2, 3)"), "-4");
        assert_eq!(s("binomial(100, 50)"), "100891344545564193334812497256");
        assert_eq!(s("binomial(100, 3)"), "161700");
        assert_eq!(
            s("binomial(10^20, 2)"),
            "4999999999999999999950000000000000000000"
        );
        assert_eq!(
            s("binomial(10^30, 10^9)"),
            "binomial(1000000000000000000000000000000, 1000000000)"
        );
        assert_eq!(s("multinomial(list(2, 3, 4))"), "1260");
        assert_eq!(s("multinomial(list())"), "1");
        assert_eq!(s("multinomial(list(5))"), "1");
        assert_eq!(s("multinomial(list(1, 1, 1, 1))"), "24");
        assert_eq!(s("multinomial(list(2, -1))"), "multinomial(list(2, -1))");
    }

    #[test]
    fn falling_and_rising() {
        assert_eq!(s("permutations(5, 2)"), "20");
        assert_eq!(s("permutations(5, 0)"), "1");
        assert_eq!(s("permutations(5, 5)"), "120");
        assert_eq!(s("permutations(3, 5)"), "0");
        assert_eq!(s("permutations(0, 0)"), "1");
        assert_eq!(s("rising(3, 4)"), "360");
        assert_eq!(s("rising(-2, 4)"), "0");
        assert_eq!(s("rising(7, 0)"), "1");
        assert_eq!(s("falling(5, 3)"), "60");
        assert_eq!(s("falling(-1, 3)"), "-6");
        assert_eq!(s("falling(3, 5)"), "0");
        assert_eq!(s("falling(5, -1)"), "falling(5, -1)");
    }

    #[test]
    fn counting_sequences() {
        let expect = |op: &str, values: &[&str]| {
            for (n, v) in values.iter().enumerate() {
                assert_eq!(s(&format!("{op}({n})")), *v, "{op}({n})");
            }
        };
        expect("catalan", &["1", "1", "2", "5", "14", "42", "132", "429"]);
        expect("fibonacci", &["0", "1", "1", "2", "3", "5", "8", "13"]);
        expect("lucas", &["2", "1", "3", "4", "7", "11", "18", "29"]);
        expect(
            "bell",
            &[
                "1", "1", "2", "5", "15", "52", "203", "877", "4140", "21147", "115975",
            ],
        );
        expect(
            "partitions",
            &["1", "1", "2", "3", "5", "7", "11", "15", "22", "30", "42"],
        );
        expect(
            "derangements",
            &["1", "0", "1", "2", "9", "44", "265", "1854"],
        );
        assert_eq!(s("catalan(30)"), "3814986502092304");
        assert_eq!(s("fibonacci(100)"), "354224848179261915075");
        assert_eq!(
            s("fibonacci(200)"),
            "280571172992510140037611932413038677189525"
        );
        assert_eq!(s("lucas(100)"), "792070839848372253127");
        assert_eq!(s("bell(20)"), "51724158235372");
        assert_eq!(s("partitions(50)"), "204226");
        assert_eq!(s("partitions(100)"), "190569292");
        assert_eq!(s("derangements(10)"), "1334961");
        assert_eq!(s("catalan(-1)"), "catalan(-1)");
        assert_eq!(s("fibonacci(-1)"), "fibonacci(-1)");
        assert_eq!(s("bell(-1)"), "bell(-1)");
        assert_eq!(s("partitions(-1)"), "partitions(-1)");
    }

    #[test]
    fn stirling_numbers() {
        assert_eq!(s("stirling1(0, 0)"), "1");
        assert_eq!(s("stirling1(3, 0)"), "0");
        assert_eq!(s("stirling1(4, 2)"), "11");
        assert_eq!(s("stirling1(5, 2)"), "50");
        assert_eq!(s("stirling1(5, 5)"), "1");
        assert_eq!(s("stirling1(10, 3)"), "1172700");
        assert_eq!(s("stirling1(3, 5)"), "0");
        assert_eq!(s("stirling2(0, 0)"), "1");
        assert_eq!(s("stirling2(4, 2)"), "7");
        assert_eq!(s("stirling2(5, 2)"), "15");
        assert_eq!(s("stirling2(6, 3)"), "90");
        assert_eq!(s("stirling2(10, 3)"), "9330");
        assert_eq!(s("stirling2(3, 5)"), "0");
        assert_eq!(s("stirling2(5, -1)"), "stirling2(5, -1)");
    }

    #[test]
    fn harmonic_numbers() {
        assert_eq!(s("harmonic(0)"), "0");
        assert_eq!(s("harmonic(1)"), "1");
        assert_eq!(s("harmonic(2)"), "3/2");
        assert_eq!(s("harmonic(5)"), "137/60");
        assert_eq!(s("harmonic(10)"), "7381/2520");
    }

    #[test]
    fn symbolic_arguments_stay() {
        let sets = [combinatorics()];
        for (src, expected) in [
            ("factorial(n)", "factorial(n)"),
            ("binomial(n, 2)", "binomial(n, 2)"),
            ("fibonacci(n)", "fibonacci(n)"),
            ("stirling2(3, k)", "stirling2(3, k)"),
            ("multinomial(list(a, 2))", "multinomial(list(a, 2))"),
            ("harmonic(n)", "harmonic(n)"),
        ] {
            let (text, reduced) = reduce_with(&sets, src, &[]);
            assert!(reduced);
            assert_eq!(text, expected);
        }
    }

    #[test]
    fn rewrites_fire() {
        let sets = [combinatorics()];
        let run =
            |src: &str, assume: &[(&str, crate::graph::Facts)]| reduce_with(&sets, src, assume).0;
        let nonneg = [("n", crate::graph::Facts::NONNEGATIVE)];
        assert_eq!(run("binomial(n, 0)", &[]), "1");
        assert_eq!(run("binomial(n, 1)", &nonneg), "n");
        assert_eq!(run("binomial(n, 1)", &[]), "binomial(n, 1)");
        assert_eq!(run("binomial(n, n)", &nonneg), "1");
        assert_eq!(run("binomial(n, n)", &[]), "binomial(n, n)");
        assert_eq!(run("permutations(n, 0)", &nonneg), "1");
        assert_eq!(run("permutations(n, 1)", &nonneg), "n");
        assert_eq!(run("permutations(n, n)", &nonneg), "factorial(n)");
        assert_eq!(run("factorial(n + 1) / factorial(n)", &nonneg), "n + 1");
        assert_eq!(
            run("factorial(n + 1) / factorial(n)", &[]),
            "factorial(n + 1)/factorial(n)"
        );
        let both = [
            ("n", crate::graph::Facts::NONNEGATIVE),
            ("k", crate::graph::Facts::NONNEGATIVE),
        ];
        assert_eq!(
            run("binomial(n, k) * factorial(k)", &both),
            "permutations(n, k)"
        );
    }

    #[test]
    fn float_semantics() {
        let sets = [combinatorics()];
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs().max(1.0);
        assert!(close(
            eval(&sets, "factorial(x)", &[("x", 10.0)]),
            3_628_800.0
        ));
        assert!(close(
            eval(&sets, "factorial(x)", &[("x", 0.5)]),
            0.886_226_925_452_758
        ));
        assert!(eval(&sets, "factorial(x)", &[("x", -2.0)]).is_nan());
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 10.0), ("y", 3.0)]),
            120.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 3.0), ("y", 5.0)]),
            0.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", -3.0), ("y", 2.0)]),
            6.0
        ));
        assert!(close(
            eval(&sets, "binomial(x, y)", &[("x", 4.0), ("y", -1.0)]),
            0.0
        ));
        assert!(close(
            eval(&sets, "permutations(x, y)", &[("x", 5.0), ("y", 2.0)]),
            20.0
        ));
    }

    #[test]
    fn fibonacci_recurrence() {
        for n in 0..60 {
            let sum = value(&format!("fibonacci({n})")) + value(&format!("fibonacci({})", n + 1));
            assert_eq!(sum, value(&format!("fibonacci({})", n + 2)), "n = {n}");
        }
        // d'Ocagne / Cassini: F(n-1) F(n+1) - F(n)^2 = (-1)^n.
        for n in 1..40_i32 {
            let lhs = value(&format!("fibonacci({})", n - 1))
                * value(&format!("fibonacci({})", n + 1))
                - value(&format!("fibonacci({n})")).pow(2);
            assert_eq!(
                lhs,
                BigInt::from(if n % 2 == 0 { 1 } else { -1 }),
                "n = {n}"
            );
        }
        // Lucas numbers satisfy the same recurrence and L(n) = F(n-1) + F(n+1).
        for n in 1..40 {
            assert_eq!(
                value(&format!("lucas({n})")),
                value(&format!("fibonacci({})", n - 1)) + value(&format!("fibonacci({})", n + 1)),
                "n = {n}"
            );
        }
    }

    #[test]
    fn binomial_recurrence() {
        // Pascal's rule.
        for n in 1..14 {
            for k in 0..=n {
                let lhs = value(&format!("binomial({n}, {k})"));
                let rhs = value(&format!("binomial({}, {})", n - 1, k - 1))
                    + value(&format!("binomial({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "C({n}, {k})");
                assert_eq!(lhs, value(&format!("binomial({n}, {})", n - k)), "symmetry");
            }
        }
        // C(n, k) = n! / (k! (n-k)!) and row sums 2^n.
        for n in 0..12 {
            let mut total = BigInt::zero();
            for k in 0..=n {
                let c = value(&format!("binomial({n}, {k})"));
                assert_eq!(
                    c * value(&format!("factorial({k})")) * value(&format!("factorial({})", n - k)),
                    value(&format!("factorial({n})"))
                );
                total += value(&format!("binomial({n}, {k})"));
            }
            assert_eq!(total, BigInt::from(2).pow(n), "row {n}");
        }
    }

    #[test]
    fn stirling_and_bell_recurrences() {
        for n in 1..10 {
            for k in 1..=n {
                let lhs = value(&format!("stirling2({n}, {k})"));
                let rhs = value(&format!("stirling2({}, {})", n - 1, k - 1))
                    + BigInt::from(k) * value(&format!("stirling2({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "S({n}, {k})");
                let lhs = value(&format!("stirling1({n}, {k})"));
                let rhs = value(&format!("stirling1({}, {})", n - 1, k - 1))
                    + BigInt::from(n - 1) * value(&format!("stirling1({}, {k})", n - 1));
                assert_eq!(lhs, rhs, "c({n}, {k})");
            }
            // Bell numbers are row sums of the Stirling numbers of the second kind.
            let total: BigInt = (0..=n)
                .map(|k| value(&format!("stirling2({n}, {k})")))
                .sum();
            assert_eq!(total, value(&format!("bell({n})")));
        }
    }

    #[test]
    fn partition_recurrence() {
        // n p(n) = sum_{k=1..n} sigma(k) p(n-k), with sigma from number theory.
        let sets = [combinatorics(), crate::rules::number_theory()];
        let v = |src: &str| -> BigInt {
            simplify(&sets, src)
                .parse()
                .unwrap_or_else(|_| panic!("`{src}` did not reduce to an integer"))
        };
        for n in 1..40 {
            let rhs: BigInt = (1..=n)
                .map(|k| v(&format!("divisor_sum({k})")) * v(&format!("partitions({})", n - k)))
                .sum();
            assert_eq!(
                BigInt::from(n) * v(&format!("partitions({n})")),
                rhs,
                "n = {n}"
            );
        }
        // Derangements: D(n) = n D(n-1) + (-1)^n.
        for n in 1..30_i32 {
            let sign = if n % 2 == 0 { 1 } else { -1 };
            assert_eq!(
                value(&format!("derangements({n})")),
                BigInt::from(n) * value(&format!("derangements({})", n - 1)) + sign
            );
        }
    }
}
