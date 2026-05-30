use rand::prelude::*;
use rand::rngs::StdRng;
use rssn::numerical::stats::variance_with_type;

/// Running stats calculator for computing variance in a single pass.
/// <https://www.johndcook.com/skewness_kurtosis.html>
#[derive(Default)]
struct RunningStats(f64, f64, f64);

impl RunningStats {
    fn push(
        &mut self,
        x: f64,
    ) {
        let n1 = self.0;
        self.0 += 1.0;
        let delta = x - self.1;
        let delta_n = delta / self.0;
        let term1 = delta * delta_n * n1;
        self.1 += delta_n;
        self.2 += term1;
    }

    fn variance(&self) -> Option<f64> {
        if self.0 < 2.0 {
            None
        } else {
            Some(self.2 / (self.0 - 1.0))
        }
    }
}

impl<T: Iterator<Item = f64>> From<T> for RunningStats {
    fn from(xs: T) -> Self {
        let mut r = Self::default();
        for x in xs {
            r.push(x);
        }
        r
    }
}

fn main() {
    let mut rng = StdRng::seed_from_u64(30_613_700); // seed from random.org
    let mut xs = vec![0f64; 1_000_000];
    // John D. Cook's second example: samples ~ U([0, 1]), shifted by 10^12
    for x in &mut xs {
        *x = 1e12 + rng.r#gen::<f64>();
    }
    println!(
        "rssn:         {}\nRunningStats: {}",
        variance_with_type(&xs, rssn::numerical::stats::VarianceType::Sample).unwrap(),
        RunningStats::from(xs.into_iter())
            .variance()
            .expect("nonempty")
    );
}
