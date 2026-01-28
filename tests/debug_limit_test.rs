use num_bigint::BigInt;
use rssn::symbolic::calculus::*;
use rssn::symbolic::core::Expr;

#[test]

fn test_lim_sinc_only() {

    let x = Expr::new_variable("x");

    let zero = Expr::new_constant(0.0);

    let sin_x =
        Expr::new_sin(x.clone());

    let sin_x_over_x =
        Expr::new_div(sin_x, x.clone());

    println!(
        "Expr: {}",
        sin_x_over_x
    );

    let lim_sinc = limit(
        &sin_x_over_x,
        "x",
        &zero,
    );

    println!(
        "Limit result variant: {:?}",
        lim_sinc
    );

    println!(
        "Limit result display: {}",
        lim_sinc
    );

    match &lim_sinc {
        | Expr::new_bigint(i) => {

            println!(
                "DEBUG: lim_sinc is \
                 BigInt({})",
                i
            )
        },
        | Expr::new_constant(c) => {

            println!(
                "DEBUG: lim_sinc is \
                 Constant({})",
                c
            )
        },
        | Expr::Dag(_) => {

            println!(
                "DEBUG: lim_sinc is \
                 Dag"
            )
        },
        | _ => {

            println!(
                "DEBUG: lim_sinc is \
                 {:?}",
                lim_sinc
            )
        },
    }

    assert_eq!(
        lim_sinc,
        Expr::new_constant(1.0)
    );
}
