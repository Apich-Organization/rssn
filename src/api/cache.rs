//! # The compute cache
//!
//! A [`Session`](crate::api::Session) remembers what it has already computed,
//! at two granularities.
//!
//! * **Answers** (request level, on by default). The store is hash-consed, so
//!   a request *is* its store node: the same term built twice is the same
//!   [`NodeId`]. An answer is keyed by that node, the template it ran on, the
//!   template generation and a canonical fingerprint of every [`Config`]
//!   field that can change the result (target, bindings, assumptions,
//!   dimension assumptions and mode, base dimensions, the budget). A hit
//!   returns a value identical to the one a fresh run would produce, report
//!   included, so enabling the cache never changes an answer. Runs that hit
//!   the wall-clock limit ([`Stop::Time`]) are not stored: their outcome
//!   depends on the machine, not on the request.
//! * **Lemmas** (sub-term level, opt-in). After a run, every *heavy*
//!   sub-request of the original term (an `integral`, a `diff`, a `solve`, …)
//!   whose class reached a closed form is recorded as the proven identity
//!   `heavy = closed`. A later request that contains the same heavy sub-term
//!   under the same context starts with that equality already asserted, so
//!   `diff(integral(f, x), x)` after `integral(f, x)` does not integrate
//!   twice. The context of a lemma is the answer fingerprint *without* the
//!   budget: an identity proven under a small budget is still an identity
//!   under a large one. Because a seeded equality can let extraction reach a
//!   different (equal) closed form than a cold run would, lemma reuse is off
//!   by default; turn it on with [`CachePolicy::lemmas`].
//!
//! Both tables are cleared when the session's templates are rebuilt (a
//! plugin or unit is registered, or the policy disables a table), so a stale
//! entry can never be returned. A [`fork`](crate::api::Session::fork) copies
//! the store and therefore the cache; a merge keeps the receiving session's
//! node ids, so its cache stays valid.
//!
//! Capacity is a bound on the number of answers and on the number of lemmas,
//! enforced first-in first-out; `capacity = 0` disables storing altogether.

use std::collections::HashMap;
use std::collections::VecDeque;
use std::fmt::Write as _;

use crate::api::ComputeError;
use crate::api::Config;
use crate::graph::DimError;
use crate::graph::Dimension;
use crate::graph::NodeId;
use crate::graph::Report;
use crate::graph::Stop;

/// What a session caches; see the [module documentation](self).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CachePolicy {
    /// Memoise whole answers (exact: a hit equals a fresh run).
    pub answers: bool,
    /// Record and reuse closed forms of heavy sub-requests.
    pub lemmas: bool,
    /// Maximum number of entries in each table.
    pub capacity: usize,
}

impl Default for CachePolicy {
    fn default() -> Self {
        Self {
            answers: true,
            lemmas: false,
            capacity: 4096,
        }
    }
}

impl CachePolicy {
    /// No caching at all.
    #[must_use]
    pub const fn off() -> Self {
        Self {
            answers: false,
            lemmas: false,
            capacity: 0,
        }
    }

    /// Answers and lemmas.
    #[must_use]
    pub const fn full() -> Self {
        Self {
            answers: true,
            lemmas: true,
            capacity: 4096,
        }
    }
}

/// Counters of a session's cache, for tests and profiling.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub struct CacheStats {
    /// Requests answered from the answer table.
    pub answer_hits: u64,
    /// Requests that ran the engine.
    pub answer_misses: u64,
    /// Answers currently stored.
    pub answers: usize,
    /// Lemmas currently stored.
    pub lemmas: usize,
    /// Lemma equalities asserted into a job before its run.
    pub lemma_hits: u64,
    /// Entries dropped to respect the capacity.
    pub evictions: u64,
    /// Times the tables were cleared because the templates changed.
    pub invalidations: u64,
}

/// An answer as the store holds it, without the session borrow.
#[derive(Clone, Debug)]
pub struct StoredAnswer {
    pub node: NodeId,
    pub value: Option<f64>,
    pub error: Option<f64>,
    pub reduced: bool,
    pub dimension: Option<Dimension>,
    pub diagnostics: Vec<DimError>,
    pub report: Report,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct AnswerKey {
    template: usize,
    root: NodeId,
    config: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct LemmaKey {
    template: usize,
    context: String,
    heavy: NodeId,
}

/// The two tables and their counters.
#[derive(Clone, Debug, Default)]
pub struct ComputeCache {
    pub policy: CachePolicy,
    answers: HashMap<AnswerKey, Result<StoredAnswer, ComputeError>>,
    answer_order: VecDeque<AnswerKey>,
    lemmas: HashMap<LemmaKey, NodeId>,
    lemma_order: VecDeque<LemmaKey>,
    stats: CacheStats,
}

/// The canonical text of every field of `config` that can change an answer,
/// except the rule sets (they select the template, which is keyed apart).
///
/// `Debug` of an `f64` is the shortest text that reads back to the same
/// bits, so two configurations with equal fingerprints are equal.
pub fn fingerprint(
    config: &Config,
    with_budget: bool,
) -> String {
    let mut text = String::new();
    let _ = write!(
        text,
        "t={:?};m={:?};by={};",
        config.target, config.dimension_mode, config.bypass_dimensions
    );
    for (name, value) in &config.bindings {
        let _ = write!(text, "b:{name}={:?};", value.to_bits());
    }
    for (name, facts) in &config.assumptions {
        let _ = write!(text, "a:{name}={facts:?};");
    }
    for (name, dim) in &config.dimension_assumptions {
        let _ = write!(text, "d:{name}={dim};");
    }
    for (name, dim) in &config.dimension_terms {
        let _ = write!(text, "D:{name}={dim:?};");
    }
    for name in &config.base_dimensions {
        let _ = write!(text, "B:{name};");
    }
    if with_budget {
        let _ = write!(text, "u={:?};", config.budget);
    }
    text
}

impl ComputeCache {
    /// Forgets everything; called when the templates are rebuilt.
    pub fn invalidate(&mut self) {
        if !self.answers.is_empty() || !self.lemmas.is_empty() {
            self.stats.invalidations = self.stats.invalidations.saturating_add(1);
        }
        self.answers.clear();
        self.answer_order.clear();
        self.lemmas.clear();
        self.lemma_order.clear();
    }

    pub fn stats(&self) -> CacheStats {
        CacheStats {
            answers: self.answers.len(),
            lemmas: self.lemmas.len(),
            ..self.stats
        }
    }

    /// The stored answer for a request, counting a hit or a miss.
    pub fn answer(
        &mut self,
        template: usize,
        root: NodeId,
        config: &str,
    ) -> Option<Result<StoredAnswer, ComputeError>> {
        if !self.policy.answers {
            return None;
        }
        let key = AnswerKey {
            template,
            root,
            config: config.to_owned(),
        };
        let found = self.answers.get(&key).cloned();
        if found.is_some() {
            self.stats.answer_hits = self.stats.answer_hits.saturating_add(1);
        } else {
            self.stats.answer_misses = self.stats.answer_misses.saturating_add(1);
        }
        found
    }

    /// Stores an answer unless the run's outcome depended on the clock.
    pub fn store_answer(
        &mut self,
        template: usize,
        root: NodeId,
        config: String,
        answer: Result<StoredAnswer, ComputeError>,
    ) {
        if !self.policy.answers || self.policy.capacity == 0 {
            return;
        }
        if let Ok(stored) = &answer
            && stored.report.stop == Stop::Time
        {
            return;
        }
        let key = AnswerKey {
            template,
            root,
            config,
        };
        if self.answers.insert(key.clone(), answer).is_none() {
            self.answer_order.push_back(key);
        }
        while self.answers.len() > self.policy.capacity {
            let Some(oldest) = self.answer_order.pop_front() else {
                break;
            };
            if self.answers.remove(&oldest).is_some() {
                self.stats.evictions = self.stats.evictions.saturating_add(1);
            }
        }
    }

    /// The closed form recorded for the heavy store node `heavy`.
    pub fn lemma(
        &self,
        template: usize,
        context: &str,
        heavy: NodeId,
    ) -> Option<NodeId> {
        if !self.policy.lemmas {
            return None;
        }
        self.lemmas
            .get(&LemmaKey {
                template,
                context: context.to_owned(),
                heavy,
            })
            .copied()
    }

    pub const fn note_lemma_hits(
        &mut self,
        count: usize,
    ) {
        self.stats.lemma_hits = self.stats.lemma_hits.saturating_add(count as u64);
    }

    /// Records `heavy = closed` (both store nodes).
    pub fn store_lemma(
        &mut self,
        template: usize,
        context: &str,
        heavy: NodeId,
        closed: NodeId,
    ) {
        if !self.policy.lemmas || self.policy.capacity == 0 || heavy == closed {
            return;
        }
        let key = LemmaKey {
            template,
            context: context.to_owned(),
            heavy,
        };
        if self.lemmas.insert(key.clone(), closed).is_none() {
            self.lemma_order.push_back(key);
        }
        while self.lemmas.len() > self.policy.capacity {
            let Some(oldest) = self.lemma_order.pop_front() else {
                break;
            };
            if self.lemmas.remove(&oldest).is_some() {
                self.stats.evictions = self.stats.evictions.saturating_add(1);
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use crate::api::CachePolicy;
    use crate::api::Config;
    use crate::api::Session;

    #[test]
    fn a_repeated_request_is_answered_from_the_cache() {
        let s = Session::new();
        let request = s.parse("diff(sin(x) * cos(x), x)").unwrap();
        let first = s.compute(request, &Config::new()).unwrap();
        let second = s.compute(request, &Config::new()).unwrap();
        assert_eq!(first.term.to_string(), second.term.to_string());
        assert_eq!(first.report.iterations, second.report.iterations);
        let stats = s.cache_stats();
        assert_eq!(stats.answer_hits, 1);
        assert_eq!(stats.answer_misses, 1);
        assert_eq!(stats.answers, 1);
        // The same term built again is the same store node: still a hit.
        let rebuilt = s.parse("diff(sin(x) * cos(x), x)").unwrap();
        let _ = s.compute(rebuilt, &Config::new()).unwrap();
        assert_eq!(s.cache_stats().answer_hits, 2);
    }

    #[test]
    fn a_different_configuration_misses() {
        let s = Session::new();
        let request = s.parse("diff(x^3, x)").unwrap();
        let symbolic = s.compute(request, &Config::new()).unwrap();
        let at_two = s
            .compute(request, &Config::new().numeric(1e-12).bind("x", 2.0))
            .unwrap();
        let at_three = s
            .compute(request, &Config::new().numeric(1e-12).bind("x", 3.0))
            .unwrap();
        assert_eq!(symbolic.term.to_string(), "3*x^2");
        assert_eq!(at_two.value, Some(12.0));
        assert_eq!(at_three.value, Some(27.0));
        assert_eq!(s.cache_stats().answer_hits, 0);
        assert_eq!(s.cache_stats().answers, 3);
    }

    #[test]
    fn the_policy_can_switch_the_cache_off() {
        let s = Session::new();
        s.set_cache_policy(CachePolicy::off());
        let request = s.parse("diff(x^3, x)").unwrap();
        let _ = s.compute(request, &Config::new()).unwrap();
        let _ = s.compute(request, &Config::new()).unwrap();
        let stats = s.cache_stats();
        assert_eq!(stats.answer_hits, 0);
        assert_eq!(stats.answers, 0);
    }

    #[test]
    fn capacity_bounds_the_table() {
        let s = Session::new();
        s.set_cache_policy(CachePolicy {
            capacity: 2,
            ..CachePolicy::default()
        });
        for k in 1..=4 {
            let request = s.parse(&format!("diff(x^{k}, x)")).unwrap();
            let _ = s.compute(request, &Config::new()).unwrap();
        }
        let stats = s.cache_stats();
        assert_eq!(stats.answers, 2);
        assert_eq!(stats.evictions, 2);
    }

    #[test]
    fn compute_all_uses_and_fills_the_cache() {
        let s = Session::new();
        let a = s.parse("diff(x^2, x)").unwrap();
        let b = s.parse("diff(x^5, x)").unwrap();
        let _ = s.compute(a, &Config::new()).unwrap();
        let all = s.compute_all(&[a, b, a], &Config::new());
        let texts: Vec<String> = all
            .into_iter()
            .map(|answer| answer.unwrap().term.to_string())
            .collect();
        assert_eq!(texts, ["2*x", "5*x^4", "2*x"]);
        // `a` twice from the cache, `b` once from the engine.
        assert_eq!(s.cache_stats().answer_hits, 2);
        let again = s.compute(b, &Config::new()).unwrap();
        assert_eq!(again.term.to_string(), "5*x^4");
        assert_eq!(s.cache_stats().answer_hits, 3);
    }

    #[test]
    fn registering_a_unit_invalidates_the_cache() {
        let mut s = Session::new();
        let request = s.parse("diff(x^3, x)").unwrap();
        let _ = s.compute(request, &Config::new()).unwrap();
        assert_eq!(s.cache_stats().answers, 1);
        let unit = crate::dim::UnitDef::linear(
            "stoneish",
            crate::kernel::signature::DimExpr::base("Mass"),
            crate::dim::Scale::one(),
        );
        if s.register_unit(unit).is_ok() {
            assert_eq!(s.cache_stats().answers, 0);
            assert_eq!(s.cache_stats().invalidations, 1);
        }
    }

    #[test]
    fn a_lemma_is_reused_by_a_larger_request() {
        let s = Session::new();
        s.set_cache_policy(CachePolicy::full());
        let inner = s.parse("integral(x^2 * exp(x), x)").unwrap();
        let first = s.compute(inner, &Config::new()).unwrap();
        assert!(first.reduced, "{}", first.term);
        assert!(s.cache_stats().lemmas >= 1);
        let outer = s
            .parse("diff(integral(x^2 * exp(x), x), x) + integral(x^2 * exp(x), x)")
            .unwrap();
        let warm = s.compute(outer, &Config::new()).unwrap();
        assert!(warm.reduced, "{}", warm.term);
        assert!(s.cache_stats().lemma_hits >= 1);
        // The warm answer is equal to the cold one (compare numerically).
        let cold_session = Session::new();
        let cold_outer = cold_session
            .parse("diff(integral(x^2 * exp(x), x), x) + integral(x^2 * exp(x), x)")
            .unwrap();
        let cold = cold_session.compute(cold_outer, &Config::new()).unwrap();
        for x in [0.3_f64, 1.7] {
            let w = s
                .compute(warm.term, &Config::new().numeric(1e-12).bind("x", x))
                .unwrap()
                .value
                .unwrap();
            let c = cold_session
                .compute(cold.term, &Config::new().numeric(1e-12).bind("x", x))
                .unwrap()
                .value
                .unwrap();
            // Antiderivatives may differ by a constant; the derivative part
            // pins everything else.
            let wd = s
                .compute(
                    s.parse(&format!("diff({}, x)", warm.term)).unwrap(),
                    &Config::new().numeric(1e-12).bind("x", x),
                )
                .unwrap()
                .value
                .unwrap();
            let cd = cold_session
                .compute(
                    cold_session
                        .parse(&format!("diff({}, x)", cold.term))
                        .unwrap(),
                    &Config::new().numeric(1e-12).bind("x", x),
                )
                .unwrap()
                .value
                .unwrap();
            assert!((wd - cd).abs() < 1e-9, "{wd} vs {cd}");
            assert!(w.is_finite() && c.is_finite());
        }
    }

    #[test]
    fn lemmas_respect_their_context() {
        let s = Session::new();
        s.set_cache_policy(CachePolicy::full());
        let request = s.parse("integral(x^2, x)").unwrap();
        let _ = s.compute(request, &Config::new()).unwrap();
        let lemmas = s.cache_stats().lemmas;
        assert!(lemmas >= 1);
        // A different assumption set is a different context: no reuse.
        let other = s.parse("integral(x^2, x) + 1").unwrap();
        let _ = s
            .compute(
                other,
                &Config::new().assume("x", crate::graph::Facts::POSITIVE),
            )
            .unwrap();
        assert_eq!(s.cache_stats().lemma_hits, 0);
    }
}
