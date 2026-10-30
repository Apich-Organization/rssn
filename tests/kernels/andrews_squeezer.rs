//! Andrews' squeezer mechanism: an index-3 DAE benchmark for [`radau_with`].
//!
//! # Reference
//!
//! The mechanism and the data below are transcribed from the canonical test-set
//! source `andrews.f` (revision 1.2, 2006-10-02) of the Test Set for IVP
//! Solvers:
//!
//! * <http://archimede.dm.uniba.it/~testset/src/problems/andrews.f>
//! * problem description: <http://archimede.uniba.it/~testset/report/andrews.pdf>
//! * IFToMM multibody benchmark:
//!   <https://iftomm-multibody.org/media/problems/andrews'%20squeezer%20mechanism/andrews'%20squeezer%20mechanism.pdf>
//!
//! The same model is described in E. Hairer and G. Wanner, *Solving Ordinary
//! Differential Equations II* (2nd rev. ed., Springer 1996), pp. 530-536. The
//! test-set interval is `[0, 0.03]` (`t(0) = 0`, `t(1) = 3d-2` in `andrews.f`);
//! its `solut` array tabulates the full 27-component state only at the
//! terminal time `t = 0.03`, so that is the reference output time used here.
//! The reference state was computed with PSIDE on a Cray C90 at
//! `rtol = atol = 1e-14` (16 significant digits).
//!
//! # Formulation
//!
//! Seven rigid bodies, six holonomic loop-closure constraints and one spring
//! (stiffness `c0`, natural length `l0`); no gravity; a constant driving torque
//! `mom` acts on the crank. With the generalized coordinates
//! `q = (β, Θ, γ, Φ, δ, Ω, ε)` the equations are
//!
//! ```text
//! q' = v
//! M(q) v' = f(q, v) - G(q)^T λ,    G = ∂g/∂q
//! 0 = g(q)
//! ```
//!
//! and are passed to [`radau_with`] in the first-order Hessenberg form
//! `K y' = φ(y)` of `andrews.f` with the 27-vector
//! `y = (q, v, a, λ)`:
//!
//! ```text
//! y'[0..7]   = v
//! y'[7..14]  = a
//! 0          = -f(q, v) + M(q) a + G(q)^T λ
//! 0          = g(q)
//! ```
//!
//! The mass matrix is `K = diag(I_14, 0_13)`. Exactly as in `andrews.f`
//! (`ind(i) = 1` for `i = 1..7`, `2` for `8..14`, `3` for `15..27`), the
//! index blocks are `n2 = 7` and `n3 = 13`, so the last 13 components carry
//! the `h^2` error weighting of the index-3 Hessenberg form.

use rssn::kernels::dense::Mat;
use rssn::kernels::dense::{
    self,
};
use rssn::kernels::jacobian::fd_jacobian_dense;
use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::ode_adaptive::OdeSolution;
use rssn::kernels::ode_stiff::DaeIndex;
use rssn::kernels::ode_stiff::RadauOptions;
use rssn::kernels::ode_stiff::radau_with;

// ---------------------------------------------------------------- mechanism

/// Body masses [kg] (`m1..m7` in `andrews.f`).
const M1: f64 = 0.04325;
const M2: f64 = 0.00365;
const M3: f64 = 0.02373;
const M4: f64 = 0.00706;
const M5: f64 = 0.07050;
const M6: f64 = 0.00706;
const M7: f64 = 0.05498;

/// Moments of inertia [kg m²] about the centres of mass.
const I1: f64 = 2.194e-6;
const I2: f64 = 4.410e-7;
const I3: f64 = 5.255e-6;
const I4: f64 = 5.667e-7;
const I5: f64 = 1.169e-5;
const I6: f64 = 5.667e-7;
const I7: f64 = 1.912e-5;

/// Fixed pivot coordinates [m].
const XA: f64 = -0.06934;
const YA: f64 = -0.00227;
const XB: f64 = -0.03635;
const YB: f64 = 0.03273;
const XC: f64 = 0.014;
const YC: f64 = 0.072;

/// Spring stiffness [N/m] and natural length [m].
const C0: f64 = 4530.0;
const L0: f64 = 0.07785;

/// Link-length data [m] (`d`, `da`, `e`, `ea`, `rr`, `ra`, `ss`, `sa`, `sb`,
/// `sc`, `sd`, `ta`, `tb`, `u`, `ua`, `ub`, `zf`, `zt`, `fa`).
const D: f64 = 0.028;
const DA: f64 = 0.0115;
const E: f64 = 0.02;
const EA: f64 = 0.01421;
const RR: f64 = 0.007;
const RA: f64 = 0.00092;
const SS: f64 = 0.035;
const SA: f64 = 0.01874;
const SB: f64 = 0.01043;
const SC: f64 = 0.018;
const SD: f64 = 0.02;
const TA: f64 = 0.02308;
const TB: f64 = 0.00916;
const U: f64 = 0.04;
const UA: f64 = 0.01228;
const UB: f64 = 0.00449;
const ZF: f64 = 0.02;
const ZT: f64 = 0.04;
const FA: f64 = 0.01421;

/// Constant driving torque [N m] on the crank `O-F`.
const MOM: f64 = 0.033;

/// Consistent initial positions `q0` from `andrews.f` (`init`).
const Q0: [f64; 7] = [
    -0.0617138900142764496358948458001,
    0.0,
    0.455279819163070380255912382449,
    0.222668390165885884674473185609,
    0.487364979543842550225598953530,
    -0.222668390165885884674473185609,
    1.23054744454982119249735015568,
];

/// Consistent initial accelerations `a0 = q''(0)` (`init`); all initial
/// velocities vanish.
const A0: [f64; 7] = [
    14222.4439199541138705911625887,
    -10666.8329399655854029433719415,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
];

/// Consistent initial multipliers `λ0` (`init`; remaining four are zero).
const LAMBDA0: [f64; 6] = [
    98.5668703962410896057654982170,
    -6.12268834425566265503114393122,
    0.0,
    0.0,
    0.0,
    0.0,
];

/// Reference state at `t = 0.03`, the `solut` array of `andrews.f` (PSIDE,
/// Cray C90, `rtol = atol = 1e-14`). Components 1-7 are `q`, 8-14 `v`,
/// 15-21 `a`, 22-27 `λ`.
#[rustfmt::skip]
const REFERENCE: [f64; 27] = [
    0.1581077119629904e+2, -0.1575637105984298e+2, 0.4082224013073101e-1,
    -0.5347301163226948e+0, 0.5244099658805304e+0, 0.5347301163226948e+0,
    0.1048080741042263e+1,
    0.1139920302151208e+4, -0.1424379294994111e+4, 0.1103291221937134e+2,
    0.1929337464421385e+2, 0.5735699284790808e+0, -0.1929337464421385e+2,
    0.3231791658026955e+0,
    -0.2463176316945196e+5, 0.5185037701610329e+5, 0.3241025686413781e+6,
    0.5667493645176213e+6, 0.1674362929479361e+5, -0.5667493645176222e+6,
    0.9826520791458422e+4,
    0.1991753333731910e+3, -0.2975531228015052e+2, 0.2306654119098399e+2,
    0.3145271365475927e+2, 0.2264249232082739e+2, 0.1161740700019673e+2,
];

/// Initial state of the 27-variable Hessenberg form.
fn initial_state() -> [f64; 27] {
    let mut y = [0.0; 27];
    y[..7].copy_from_slice(&Q0);
    y[14..21].copy_from_slice(&A0);
    y[21..27].copy_from_slice(&LAMBDA0);
    y
}

/// `M(q)`, the 7 x 7 configuration-dependent mass matrix of `andrews.f`
/// (`feval`, lines `m(1,1)` ... `m(7,7)`).
fn mass_matrix(q: &[f64]) -> Mat {
    let (_, b2, _, b4, _, b6, _) = (q[0], q[1], q[2], q[3], q[4], q[5], q[6]);
    let siph = b4.sin();
    let siom = b6.sin();
    let coth = b2.cos();
    let mut m = Mat::zeros(7, 7);
    m.set(
        0,
        0,
        M1 * RA * RA + M2 * (RR * RR - 2.0 * DA * RR * coth + DA * DA) + I1 + I2,
    );
    m.set(1, 0, M2 * (DA * DA - DA * RR * coth) + I2);
    m.set(1, 1, M2 * DA * DA + I2);
    m.set(2, 2, M3 * (SA * SA + SB * SB) + I3);
    m.set(3, 3, M4 * (E - EA) * (E - EA) + I4);
    m.set(4, 3, M4 * ((E - EA) * (E - EA) + ZT * (E - EA) * siph) + I4);
    m.set(
        4,
        4,
        M4 * (ZT * ZT + 2.0 * ZT * (E - EA) * siph + (E - EA) * (E - EA))
            + M5 * (TA * TA + TB * TB)
            + I4
            + I5,
    );
    m.set(5, 5, M6 * (ZF - FA) * (ZF - FA) + I6);
    m.set(
        6,
        5,
        M6 * ((ZF - FA) * (ZF - FA) - U * (ZF - FA) * siom) + I6,
    );
    m.set(
        6,
        6,
        M6 * ((ZF - FA) * (ZF - FA) - 2.0 * U * (ZF - FA) * siom + U * U)
            + M7 * (UA * UA + UB * UB)
            + I6
            + I7,
    );
    for j in 1..7 {
        for i in 0..j {
            m.set(i, j, m.at(j, i));
        }
    }
    m
}

/// The six loop-closure constraints `g(q)` (`feval`, lines `g(1)` ... `g(6)`).
fn constraints(q: &[f64]) -> [f64; 6] {
    let (b1, b2, b3, b4, b5, b6, b7) = (q[0], q[1], q[2], q[3], q[4], q[5], q[6]);
    let sibe = b1.sin();
    let cobe = b1.cos();
    let sibeth = (b1 + b2).sin();
    let cobeth = (b1 + b2).cos();
    let siga = b3.sin();
    let coga = b3.cos();
    let siphde = (b4 + b5).sin();
    let cophde = (b4 + b5).cos();
    let side = b5.sin();
    let code = b5.cos();
    let siomep = (b6 + b7).sin();
    let coomep = (b6 + b7).cos();
    let siep = b7.sin();
    let coep = b7.cos();
    [
        RR * cobe - D * cobeth - SS * siga - XB,
        RR * sibe - D * sibeth + SS * coga - YB,
        RR * cobe - D * cobeth - E * siphde - ZT * code - XA,
        RR * sibe - D * sibeth + E * cophde - ZT * side - YA,
        RR * cobe - D * cobeth - ZF * coomep - U * siep - XA,
        RR * sibe - D * sibeth - ZF * siomep + U * coep - YA,
    ]
}

/// `G(q) = ∂g/∂q`, the 6 x 7 constraint Jacobian (`feval`, `gp`).
fn constraint_jacobian(q: &[f64]) -> Mat {
    let (b1, b2, b3, b4, b5, b6, b7) = (q[0], q[1], q[2], q[3], q[4], q[5], q[6]);
    let sibe = b1.sin();
    let cobe = b1.cos();
    let sibeth = (b1 + b2).sin();
    let cobeth = (b1 + b2).cos();
    let siga = b3.sin();
    let coga = b3.cos();
    let siphde = (b4 + b5).sin();
    let cophde = (b4 + b5).cos();
    let side = b5.sin();
    let code = b5.cos();
    let siomep = (b6 + b7).sin();
    let coomep = (b6 + b7).cos();
    let siep = b7.sin();
    let coep = b7.cos();
    let mut g = Mat::zeros(6, 7);
    g.set(0, 0, -RR * sibe + D * sibeth);
    g.set(0, 1, D * sibeth);
    g.set(0, 2, -SS * coga);
    g.set(1, 0, RR * cobe - D * cobeth);
    g.set(1, 1, -D * cobeth);
    g.set(1, 2, -SS * siga);
    g.set(2, 0, -RR * sibe + D * sibeth);
    g.set(2, 1, D * sibeth);
    g.set(2, 3, -E * cophde);
    g.set(2, 4, -E * cophde + ZT * side);
    g.set(3, 0, RR * cobe - D * cobeth);
    g.set(3, 1, -D * cobeth);
    g.set(3, 3, -E * siphde);
    g.set(3, 4, -E * siphde - ZT * code);
    g.set(4, 0, -RR * sibe + D * sibeth);
    g.set(4, 1, D * sibeth);
    g.set(4, 5, ZF * siomep);
    g.set(4, 6, ZF * siomep - U * coep);
    g.set(5, 0, RR * cobe - D * cobeth);
    g.set(5, 1, -D * cobeth);
    g.set(5, 5, -ZF * coomep);
    g.set(5, 6, -ZF * coomep - U * siep);
    g
}

/// Generalized forces `f(q, v)`, including the spring and the driving torque
/// (`feval`, `f(1)` ... `f(7)`).
fn generalized_forces(
    q: &[f64],
    v: &[f64],
) -> [f64; 7] {
    let (_, b2, b3, b4, _, b6, _) = (q[0], q[1], q[2], q[3], q[4], q[5], q[6]);
    let sith = b2.sin();
    let coph = b4.cos();
    let coom = b6.cos();
    let coga = b3.cos();
    let siga = b3.sin();
    let (bep, thp, php, dep, omp, epp) = (v[0], v[1], v[3], v[4], v[5], v[6]);

    let xd = SD * coga + SC * siga + XB;
    let yd = SD * siga - SC * coga + YB;
    let lang = ((xd - XC) * (xd - XC) + (yd - YC) * (yd - YC)).sqrt();
    let force = -C0 * (lang - L0) / lang;
    let fx = force * (xd - XC);
    let fy = force * (yd - YC);

    [
        MOM - M2 * DA * RR * thp * (thp + 2.0 * bep) * sith,
        M2 * DA * RR * bep * bep * sith,
        fx * (SC * coga - SD * siga) + fy * (SD * coga + SC * siga),
        M4 * ZT * (E - EA) * dep * dep * coph,
        -M4 * ZT * (E - EA) * php * (php + 2.0 * dep) * coph,
        -M6 * U * (ZF - FA) * epp * epp * coom,
        M6 * U * (ZF - FA) * omp * (omp + 2.0 * epp) * coom,
    ]
}

/// Right-hand side `φ(y)` of `K y' = φ(y)`, exactly as in `andrews.f`
/// (`feval`, lines `ff(i)`).
fn rhs(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
) {
    let q = &y[0..7];
    let v = &y[7..14];
    let a = &y[14..21];
    let lam = &y[21..27];
    let m = mass_matrix(q);
    let g = constraint_jacobian(q);
    let f = generalized_forces(q, v);

    for i in 0..7 {
        dy[i] = v[i];
        dy[7 + i] = a[i];
    }
    for i in 0..7 {
        let mut s = -f[i];
        for j in 0..7 {
            s += m.at(i, j) * a[j];
        }
        for j in 0..6 {
            s += g.at(j, i) * lam[j];
        }
        dy[14 + i] = s;
    }
    let gg = constraints(q);
    for j in 0..6 {
        dy[21 + j] = gg[j];
    }
}

/// The approximate Jacobian used by `andrews.f` (`jeval`): derivatives of
/// `f`, `M` and `G` are neglected, only the Hessenberg block structure is
/// retained (Hairer & Wanner, p. 540).
fn jacobian(
    _t: f64,
    y: &[f64],
) -> Mat {
    let q = &y[0..7];
    let m = mass_matrix(q);
    let g = constraint_jacobian(q);
    let mut j = Mat::zeros(27, 27);
    for i in 0..14 {
        j.set(i, i + 7, 1.0);
    }
    for i in 0..7 {
        for k in 0..7 {
            j.set(14 + k, 14 + i, m.at(k, i));
        }
    }
    for i in 0..6 {
        for k in 0..7 {
            j.set(14 + k, 21 + i, g.at(i, k));
        }
    }
    for k in 0..6 {
        for i in 0..7 {
            j.set(21 + k, i, g.at(k, i));
        }
    }
    j
}

/// Mass matrix `K = diag(I_14, 0_13)` of the first-order form.
fn mass_matrix_27() -> Mat {
    let mut k = Mat::identity(27);
    for i in 14..27 {
        k.set(i, i, 0.0);
    }
    k
}

/// Integrate the benchmark from `0` to `t1` with the given tolerances.
fn integrate(
    t1: f64,
    rtol: f64,
    atol: f64,
) -> OdeSolution {
    let y0 = initial_state();
    let mass = mass_matrix_27();
    let jac = |t: f64, y: &[f64]| jacobian(t, y);
    let opts = OdeOptions {
        rtol,
        atol,
        h0: Some(1e-6),
        // The mechanism has a sub-millisecond time scale (IFToMM benchmark
        // statement); capping the step also keeps the simplified Newton
        // iteration well behaved at loose tolerances.
        hmax: 1e-3,
        max_steps: 200_000,
        ..OdeOptions::default()
    };
    let ro = RadauOptions {
        jac: Some(&jac),
        mass: Some(&mass),
        index: DaeIndex { n2: 7, n3: 13 },
        min_stages: 3,
        max_stages: 3,
        ..RadauOptions::default()
    };
    radau_with(rhs, 0.0, t1, &y0, &opts, &ro).expect("Andrews' squeezer integration failed")
}

/// Jacobian of [`rhs`] at a stage value. The linear Hessenberg blocks (the
/// `q' = v`, `v' = a` identities and the `M`, `G^T`, `G` blocks of the
/// algebraic rows) are inserted exactly; the remaining derivatives of the
/// generalized forces and of `M(q) a` are finite-difference approximations.
/// The exact linear blocks avoid differencing the tiny mass-matrix entries
/// (`M ~ 1e-6`), which is what limits a pure finite-difference Jacobian here.
fn stage_jacobian(
    t: f64,
    y: &[f64],
) -> Mat {
    let n = 27;
    let mut f0 = vec![0.0; n];
    rhs(t, y, &mut f0);
    let mut j = fd_jacobian_dense(&rhs, t, y, &f0);
    let q = &y[0..7];
    let m = mass_matrix(q);
    let g = constraint_jacobian(q);
    for i in 0..14 {
        for c in 0..n {
            j.set(i, c, 0.0);
        }
        j.set(i, i + 7, 1.0);
    }
    for i in 0..7 {
        for k in 0..7 {
            j.set(14 + i, 14 + k, m.at(i, k));
        }
        for k in 0..6 {
            j.set(14 + i, 21 + k, g.at(k, i));
        }
    }
    for k in 0..6 {
        for c in 0..7 {
            j.set(21 + k, c, g.at(k, c));
        }
        for c in 7..n {
            j.set(21 + k, c, 0.0);
        }
    }
    j
}

/// Constant-step Radau IIA(5) for the same 27-variable DAE, used only by the
/// convergence test. The three-stage tableau is the one of Hairer & Wanner,
/// *Solving Ordinary Differential Equations II*, Section IV.5 (the method
/// implemented by `radau5`): nodes `c = ((4 - √6)/10, (4 + √6)/10, 1)` and
///
/// ```text
/// A = [ (88 - 7√6)/360   (296 - 169√6)/1800   (-2 + 3√6)/225 ]
///     [ (296 + 169√6)/1800 (88 + 7√6)/360      (-2 - 3√6)/225 ]
///     [ (16 - √6)/36       (16 + √6)/36        1/9          ]
/// ```
///
/// The stage equations are `K Z_i = h Σ_j A_ij φ(y + Z_j)`; because the
/// method is stiffly accurate (`c_3 = 1`) the new state is `y + Z_3`. The
/// Newton iteration uses [`stage_jacobian`] and stops when the stage update
/// reaches the linear-solve floor (the stage matrix is ill-conditioned).
fn fixed_step_radau5(
    t1: f64,
    n_steps: usize,
) -> Vec<f64> {
    let s6 = 6.0_f64.sqrt();
    let c = [(4.0 - s6) / 10.0, (4.0 + s6) / 10.0, 1.0];
    let a = [
        [
            (88.0 - 7.0 * s6) / 360.0,
            (296.0 - 169.0 * s6) / 1800.0,
            (-2.0 + 3.0 * s6) / 225.0,
        ],
        [
            (296.0 + 169.0 * s6) / 1800.0,
            (88.0 + 7.0 * s6) / 360.0,
            (-2.0 - 3.0 * s6) / 225.0,
        ],
        [(16.0 - s6) / 36.0, (16.0 + s6) / 36.0, 1.0 / 9.0],
    ];
    let n = 27;
    let h = t1 / n_steps as f64;
    let mut y = initial_state().to_vec();
    for step in 0..n_steps {
        let t = step as f64 * h;
        let mut z = vec![0.0; 3 * n];
        let mut converged = false;
        for _ in 0..60 {
            let mut stage_y = vec![vec![0.0; n]; 3];
            let mut stage_f = vec![vec![0.0; n]; 3];
            for i in 0..3 {
                stage_y[i].copy_from_slice(&y);
                for k in 0..n {
                    stage_y[i][k] += z[i * n + k];
                }
                rhs(t + c[i] * h, &stage_y[i], &mut stage_f[i]);
            }
            let jfs: Vec<Mat> = (0..3)
                .map(|j| stage_jacobian(t + c[j] * h, &stage_y[j]))
                .collect();
            let mut r = vec![0.0; 3 * n];
            let mut mat = Mat::zeros(3 * n, 3 * n);
            for i in 0..3 {
                for k in 0..n {
                    let mut s = if k < 14 { z[i * n + k] } else { 0.0 };
                    for j in 0..3 {
                        s -= h * a[i][j] * stage_f[j][k];
                    }
                    // Newton solves J dz = -R, so store -R.
                    r[i * n + k] = -s;
                }
                for j in 0..3 {
                    let jf = &jfs[j];
                    for row in 0..n {
                        for col in 0..n {
                            let mut val = -h * a[i][j] * jf.at(row, col);
                            if i == j && row == col && row < 14 {
                                val += 1.0;
                            }
                            mat.set(i * n + row, j * n + col, val);
                        }
                    }
                }
            }
            let lu = dense::lu_factor(&mat).expect("fixed-step stage matrix is singular");
            let dz = lu.solve(&r);
            let dn = dz.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            for (zi, d) in z.iter_mut().zip(&dz) {
                *zi += d;
            }
            // The stage matrix is ill-conditioned (M ~ 1e-6, h small), so the
            // linear solve bottoms out in the 1e-7 ... 1e-6 range; this floor
            // is below the truncation error of the step sizes used for the
            // order measurement.
            if dn <= 1e-6 {
                converged = true;
                break;
            }
        }
        assert!(
            converged,
            "fixed-step Newton did not converge at step {step}"
        );
        for k in 0..n {
            y[k] += z[2 * n + k];
        }
    }
    y
}

/// Total mechanical energy `T + V` and the spring potential.
fn mechanical_energy(
    q: &[f64],
    v: &[f64],
) -> (f64, f64) {
    let m = mass_matrix(q);
    let mut t = 0.0;
    for i in 0..7 {
        for j in 0..7 {
            t += 0.5 * v[i] * m.at(i, j) * v[j];
        }
    }
    let xd = SD * q[2].cos() + SC * q[2].sin() + XB;
    let yd = SD * q[2].sin() - SC * q[2].cos() + YB;
    let lang = ((xd - XC) * (xd - XC) + (yd - YC) * (yd - YC)).sqrt();
    (t, 0.5 * C0 * (lang - L0) * (lang - L0))
}

// ------------------------------------------------------------------- tests

/// At `t = 0.03` the computed state is compared with the 27 values of
/// `andrews.f`'s `solut`. The reference was computed at `rtol = atol = 1e-14`
/// and is tabulated to 16 digits, so it is much more accurate than the
/// `1e-10` integration used here; the asserted bounds are therefore dominated
/// by our solver error, not by the reference. The generalized coordinates
/// (components 1-7) are required to agree to `5e-7` absolute, a few units in
/// the last digit of the reference's `t = 0.03` values; velocities to `1e-4`
/// (relative `~1e-7` on values of order `1e3`). The algebraic components
/// (`a` and `λ`) are deliberately weighted with `h^2` by the index-3 error
/// control, so they are only required to `5e-3` relative; the differential
/// variables are the physically meaningful reference quantities.
#[test]
fn andrews_squeezer_reference() {
    let sol = integrate(0.03, 1e-10, 1e-12);
    let y = sol.y.last().unwrap();
    let mut worst_q = 0.0_f64;
    let mut worst_q_at = 0;
    for i in 0..7 {
        let d = (y[i] - REFERENCE[i]).abs();
        if d > worst_q {
            worst_q = d;
            worst_q_at = i;
        }
    }
    let mut worst_v = 0.0_f64;
    for i in 7..14 {
        worst_v = worst_v.max((y[i] - REFERENCE[i]).abs());
    }
    let mut worst_a = 0.0_f64;
    for i in 14..21 {
        let scale = REFERENCE[i].abs().max(1.0);
        worst_a = worst_a.max((y[i] - REFERENCE[i]).abs() / scale);
    }
    let mut worst_l = 0.0_f64;
    for i in 21..27 {
        let scale = REFERENCE[i].abs().max(1.0);
        worst_l = worst_l.max((y[i] - REFERENCE[i]).abs() / scale);
    }
    println!(
        "reference: steps={} nfev={} rej={} |dq|max={worst_q:.3e} (q{}) |dv|max={worst_v:.3e} \
         |da|rel={worst_a:.3e} |dλ|rel={worst_l:.3e}",
        sol.t.len(),
        sol.nfev,
        sol.rejected,
        worst_q_at + 1
    );
    assert!(
        worst_q < 5e-7,
        "max |Δq| = {worst_q:.3e} at component {}",
        worst_q_at + 1
    );
    assert!(worst_v < 1e-4, "max |Δv| = {worst_v:.3e}");
    assert!(worst_a < 5e-3, "max relative |Δa| = {worst_a:.3e}");
    assert!(worst_l < 5e-3, "max relative |Δλ| = {worst_l:.3e}");
    // The solver must actually reach t = 0.03.
    assert!((sol.t.last().unwrap() - 0.03).abs() < 1e-15);
    assert!(
        sol.orders.iter().all(|&o| o == 5),
        "orders: {:?}",
        sol.orders
    );
}

/// The loop-closure constraints are part of the index-3 system, so they must
/// be satisfied at the final time without any projection. The bound is the
/// integration tolerance scale (`rtol = 1e-9` on coordinates of order 10),
/// i.e. the algebraic equations are solved to solver accuracy.
#[test]
fn andrews_squeezer_constraints() {
    let sol = integrate(0.03, 1e-9, 1e-11);
    let y = sol.y.last().unwrap();
    let g = constraints(&y[..7]);
    let gmax = g.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
    println!("constraints at t=0.03: {g:?} max={gmax:.3e}");
    assert!(gmax < 1e-7, "max |g(q(0.03))| = {gmax:.3e}");
}

/// Energy balance from the IFToMM problem statement. With no gravity the
/// mechanical energy is `T + V`, where `V = c0 (L - l0)^2 / 2` is the spring
/// potential and `T = v^T M(q) v / 2`. The driving torque adds the work
/// `W = mom (q1 - q1(0))`, so `U = T + V - W` must stay constant; the IFToMM
/// benchmark asks for a drift below `1e-5` J over the motion. We integrate
/// with `rtol = 1e-9` and assert that benchmark threshold.
#[test]
fn andrews_squeezer_energy_balance() {
    let sol = integrate(0.03, 1e-9, 1e-11);
    let y0 = initial_state();
    let (t0, v0) = mechanical_energy(&y0[..7], &y0[7..14]);
    let u0 = t0 + v0;
    let mut worst = 0.0_f64;
    for y in &sol.y {
        let (t, v) = mechanical_energy(&y[..7], &y[7..14]);
        let u = t + v - MOM * (y[0] - Q0[0]);
        worst = worst.max((u - u0).abs());
    }
    println!(
        "energy: U0={u0:.15} max drift={worst:.3e} J ({} points)",
        sol.t.len()
    );
    assert!(worst < 1e-5, "max energy-balance drift = {worst:.3e} J");
}

/// Fifth-order convergence of the generalized coordinates. The reference at
/// `t = 0.002` is an adaptive `radau_with` run at `rtol = 1e-10`; the
/// constant-step Radau IIA(5) solutions use `1, 2, 4` steps, so `h` halves
/// each time. Halving `h` must reduce the error in the index-1 coordinates
/// `q` by about `2^5 = 32`; the observed orders below are
/// `log2(e(h) / e(h/2))`. Radau IIA is order 5 and the positions are seen to
/// converge at `~5.0`. The index-2 velocities (components 8-14) show the
/// classical order reduction of the index-3 Hessenberg form (here `~2.7`,
/// consistent with the `O(h^3)` error of the algebraic accelerations being
/// integrated into `v`); this is a property of the DAE block structure, not
/// of the reference comparison, which is made on all 27 state components in
/// [`andrews_squeezer_reference`]. The convergence assertion is therefore on
/// the generalized coordinates, whose fifth-order behaviour is the point of
/// the index-3 error weighting.
#[test]
fn andrews_squeezer_convergence() {
    let t1 = 0.002;
    let reference = integrate(t1, 1e-10, 1e-12);
    let yref = reference.y.last().unwrap();
    let mut eq = Vec::new();
    let mut ev = Vec::new();
    for &n_steps in &[1usize, 2, 4] {
        let y = fixed_step_radau5(t1, n_steps);
        let dq = (0..7).fold(0.0_f64, |m, i| m.max((y[i] - yref[i]).abs()));
        let dv = (7..14).fold(0.0_f64, |m, i| m.max((y[i] - yref[i]).abs()));
        println!(
            "steps={n_steps} h={:.3e} |Δq|={dq:.3e} |Δv|={dv:.3e}",
            t1 / n_steps as f64
        );
        eq.push(dq);
        ev.push(dv);
    }
    let mut order_q = Vec::new();
    let mut order_v = Vec::new();
    for k in 0..eq.len() - 1 {
        order_q.push((eq[k] / eq[k + 1]).log2());
        order_v.push((ev[k] / ev[k + 1]).log2());
        println!("order q={:.2} v={:.2}", order_q[k], order_v[k]);
    }
    let p_q = *order_q.last().unwrap();
    let p_v = *order_v.last().unwrap();
    assert!(p_q > 4.5 && p_q < 6.5, "observed position order {p_q:.2}");
    assert!(
        *eq.last().unwrap() < 1e-7,
        "position error {:e}",
        eq.last().unwrap()
    );
    // The velocity block is index 2; its reduced order is expected and is
    // reported rather than required to be 5.
    assert!(p_v > 2.0 && p_v < 3.5, "observed velocity order {p_v:.2}");
    assert!(
        *ev.last().unwrap() < 1e-2,
        "velocity error {:e}",
        ev.last().unwrap()
    );
}
