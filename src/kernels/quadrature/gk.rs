//! Adaptive Gauss–Kronrod panels and their error-driven partition.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::collections::VecDeque;

use super::GkError;
use super::GkRule;
use super::Integral;

const XGK15: [f64; 8] = [
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.0,
];
const WGK15: [f64; 8] = [
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714,
];
const WG7: [f64; 4] = [
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327,
];

const XGK21: [f64; 11] = [
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.0,
];
const WGK21: [f64; 11] = [
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
];
const WG10: [f64; 5] = [
    0.066671344308688137593568809893332,
    0.149451349150580593145776339657697,
    0.219086362515982043995534934228163,
    0.269266719309996355091226921569469,
    0.295524224714752870173892994651338,
];

/// The panel arithmetic. `eval` receives the centre first and then each pair
/// `c - dx_j`, `c + dx_j`, and stores each result at its node's index.
pub(super) fn panel_with(
    a: f64,
    b: f64,
    rule: GkRule,
    eval: impl FnOnce(&[f64], &mut [f64]),
) -> (f64, f64, usize) {
    let (xgk, wgk, wg): (&[f64], &[f64], &[f64]) = match rule {
        | GkRule::G7K15 => (&XGK15, &WGK15, &WG7),
        | GkRule::G10K21 => (&XGK21, &WGK21, &WG10),
    };
    let pairs = xgk.len() - 1;
    let count = 2 * pairs + 1;
    let centre = f64::midpoint(a, b);
    let width = b - a;
    let half = if width.is_finite() {
        0.5 * width
    } else {
        0.5 * b - 0.5 * a
    };
    let mut nodes = [0.0; 21];
    nodes[0] = centre;
    let (node_pairs, _) = nodes[1..count].as_chunks_mut::<2>();
    for (pair, &node) in node_pairs.iter_mut().zip(xgk) {
        let dx = half * node;
        pair[0] = centre - dx;
        pair[1] = centre + dx;
    }
    let mut values = [0.0; 21];
    eval(&nodes[..count], &mut values[..count]);
    let fc = values[0];
    let mut resk = wgk[pairs] * fc;
    let mut resg = match rule {
        | GkRule::G7K15 => WG7[3] * fc,
        | GkRule::G10K21 => 0.0,
    };
    let mut resabs = wgk[pairs] * fc.abs();
    let (values, _) = values[1..count].as_chunks::<2>();
    for (j, &[f1, f2]) in values.iter().enumerate() {
        resk += wgk[j] * (f1 + f2);
        resabs += wgk[j] * (f1.abs() + f2.abs());
        if j % 2 == 1 {
            resg += wg[j / 2] * (f1 + f2);
        }
    }
    let reskh = resk * 0.5;
    let mut resasc = wgk[pairs] * (fc - reskh).abs();
    for (weight, &[f1, f2]) in wgk.iter().zip(values) {
        resasc += weight * ((f1 - reskh).abs() + (f2 - reskh).abs());
    }
    // Halving a subnormal width rounds before the panel's sum is applied.
    let scale_to_interval = |sum: f64| {
        if width.abs() < f64::MIN_POSITIVE {
            (sum * 0.5) * width
        } else {
            sum * half
        }
    };
    let result = scale_to_interval(resk);
    resabs = scale_to_interval(resabs).abs();
    resasc = scale_to_interval(resasc).abs();
    let mut err = scale_to_interval(resk - resg).abs();
    if resasc != 0.0 && err != 0.0 {
        err = resasc * 1.0_f64.min((200.0 * err / resasc).powf(1.5));
    }
    if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        err = err.max(50.0 * f64::EPSILON * resabs);
    }
    (result, err, count)
}

#[derive(Clone, Copy, Default)]
struct Estimate {
    value: f64,
    error: f64,
}

impl Estimate {
    fn merge(
        self,
        right: Self,
    ) -> Self {
        Self {
            value: self.value + right.value,
            error: self.error + right.error,
        }
    }

    fn integral(
        self,
        evaluations: usize,
        tolerance: f64,
    ) -> Integral {
        Integral {
            value: self.value,
            error: self.error,
            evaluations,
            converged: self.error <= tolerance,
        }
    }
}

/// Keeps a pairwise reduction of the panel estimates in insertion order.
/// Revising a panel touches only its ancestors, without subtracting an old
/// contribution from a rounded total.
struct Totals {
    nodes: Vec<Estimate>,
    leaves: usize,
}

impl Totals {
    fn new() -> Self {
        Self {
            nodes: vec![Estimate::default(); 2],
            leaves: 0,
        }
    }

    fn combined(&self) -> Estimate {
        self.nodes[1]
    }

    fn set(
        &mut self,
        leaf: usize,
        estimate: Estimate,
    ) {
        let mut index = self.nodes.len() / 2 + leaf;
        self.nodes[index] = estimate;
        while index > 1 {
            index /= 2;
            self.nodes[index] = self.nodes[2 * index].merge(self.nodes[2 * index + 1]);
        }
    }

    fn push(
        &mut self,
        estimate: Estimate,
    ) {
        let capacity = self.nodes.len() / 2;
        if self.leaves == capacity {
            // The previous tree becomes the entire left subtree. Copying by
            // level preserves its rounded partial sums exactly.
            self.nodes.resize(self.nodes.len() * 2, Estimate::default());
            let mut width = capacity;
            while width > 0 {
                self.nodes.copy_within(width..2 * width, 2 * width);
                self.nodes[3 * width..4 * width].fill(Estimate::default());
                width /= 2;
            }
        }
        self.set(self.leaves, estimate);
        self.leaves += 1;
    }
}

/// Endpoints plus the lookup index and error that determine refinement order.
struct Panel {
    lo: f64,
    hi: f64,
    leaf: usize,
    error: f64,
}

impl PartialEq for Panel {
    fn eq(
        &self,
        other: &Self,
    ) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for Panel {}

impl PartialOrd for Panel {
    fn partial_cmp(
        &self,
        other: &Self,
    ) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Panel {
    fn cmp(
        &self,
        other: &Self,
    ) -> Ordering {
        self.error
            .total_cmp(&other.error)
            .then_with(|| other.leaf.cmp(&self.leaf))
    }
}

/// Bound the size and cost of the epsilon table even with a large panel budget.
const EXTRAPOLATION_WINDOW: usize = 32;

pub(super) fn integrate(
    gk: impl Fn(f64, f64) -> (f64, f64, usize),
    initial: impl IntoIterator<Item = (f64, f64)>,
    abs_tol: f64,
    rel_tol: f64,
    max_panels: usize,
) -> Result<Integral, (GkError, Integral)> {
    let mut evaluations = 0;
    let tolerance = |value: f64| abs_tol.max(rel_tol * value.abs());
    let measure = |lo, hi, leaf, evaluations: &mut usize| {
        if crate::cancellation::is_cancelled() {
            return Err((GkError::Cancelled, *evaluations));
        }
        let (value, error, count) = gk(lo, hi);
        *evaluations += count;
        if crate::cancellation::is_cancelled() {
            return Err((GkError::Cancelled, *evaluations));
        }
        if !(value.is_finite() && error.is_finite()) {
            return Err((GkError::NonFiniteEstimate, *evaluations));
        }
        Ok((Panel { lo, hi, leaf, error }, Estimate { value, error }))
    };
    let invalid = |(error, evaluations)| {
        (
            error,
            Integral {
                value: f64::NAN,
                error: f64::INFINITY,
                evaluations,
                converged: false,
            },
        )
    };
    let mut initial = initial.into_iter().peekable();
    let Some((lo, hi)) = initial.next() else {
        return Ok(Estimate::default().integral(0, 0.0));
    };
    let (panel, estimate) = measure(lo, hi, 0, &mut evaluations).map_err(invalid)?;
    let first_tolerance = tolerance(estimate.value);
    if initial.peek().is_none() && (estimate.error <= first_tolerance || max_panels <= 1) {
        return Ok(estimate.integral(evaluations, first_tolerance));
    }
    let mut heap = BinaryHeap::from([panel]);
    let mut totals = Totals::new();
    totals.push(estimate);
    for (lo, hi) in initial {
        let (panel, estimate) =
            measure(lo, hi, totals.leaves, &mut evaluations).map_err(invalid)?;
        totals.push(estimate);
        heap.push(panel);
    }
    let mut total = totals.combined();
    let mut recent = VecDeque::from([total.value]);
    while total.value.is_finite()
        && total.error.is_finite()
        && total.error > tolerance(total.value)
        && totals.leaves < max_panels
    {
        let Some(panel) = heap.pop() else {
            break;
        };
        let mid = f64::midpoint(panel.lo, panel.hi);
        if mid <= panel.lo.min(panel.hi) || mid >= panel.lo.max(panel.hi) {
            continue;
        }
        let (left, left_estimate) =
            measure(panel.lo, mid, panel.leaf, &mut evaluations).map_err(invalid)?;
        let (right, right_estimate) =
            measure(mid, panel.hi, totals.leaves, &mut evaluations).map_err(invalid)?;
        totals.set(left.leaf, left_estimate);
        totals.push(right_estimate);
        heap.push(left);
        heap.push(right);
        total = totals.combined();
        if recent.len() == EXTRAPOLATION_WINDOW {
            recent.pop_front();
        }
        recent.push_back(total.value);
    }
    if crate::cancellation::is_cancelled() {
        return Err(invalid((GkError::Cancelled, evaluations)));
    }
    if !(total.value.is_finite() && total.error.is_finite()) {
        return Err(invalid((GkError::NonFiniteEstimate, evaluations)));
    }
    let mut result = total.integral(evaluations, tolerance(total.value));
    if !result.converged && recent.len() > 8 {
        let accelerated = crate::kernels::convergence::wynn_epsilon(recent.make_contiguous());
        if let Some(&last) = accelerated.get(1..).and_then(<[f64]>::last) {
            let estimate = (last - total.value).abs().max(f64::EPSILON * last.abs());
            if last.is_finite() && estimate < total.error {
                result.value = last;
                result.error = estimate;
                result.converged = estimate <= tolerance(last);
            }
        }
    }
    Ok(result)
}
