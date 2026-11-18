//! Rustdoc of the public functions behind operators, harvested at build time.
//!
//! `build.rs` (through `buildsrc/rustdocs.rs`) scans `src/**/*.rs` and writes
//! every documented `pub fn` into one text table that this module embeds. The
//! manifest then fills an operator's `doc` from it when the registration did
//! not say enough, and records *where* the text came from
//! ([`DocSource`](crate::kernel::manifest::DocSource)), so a host can tell a
//! description the rule set wrote from one it borrowed.
//!
//! # How an operator finds its function
//!
//! 1. **`rust_path`.** An operator that declares the public function it
//!    exports (`crate::kernels::special::erf`) is matched by that exact path.
//!    This is a statement of the registration, not a guess.
//! 2. **Name.** Otherwise a documented free function with the operator's own
//!    name is used, if there is exactly one, or exactly one in a module whose
//!    path names the operator's rule set or category. Two candidates that
//!    cannot be told apart are left alone: a wrong description is worse than
//!    none.
//!
//! Nothing here changes what an operator does, and nothing is read at run
//! time: the table is a `&'static str` and the index is built on first use.

use std::collections::HashMap;
use std::sync::OnceLock;

use crate::graph::OpDescriptor;
use crate::kernel::manifest::DocSource;
use crate::kernel::signature::Example;

static TABLE: &str = include_str!(concat!(env!("OUT_DIR"), "/rustdocs.txt"));

/// One documented public function.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct RustDoc {
    /// Module path, `crate::kernels::special::erf`.
    pub path: &'static str,
    /// The function's own name.
    pub name: &'static str,
    /// The doc comment, as markdown.
    pub doc: &'static str,
}

struct Index {
    by_path: HashMap<&'static str, RustDoc>,
    by_name: HashMap<&'static str, Vec<RustDoc>>,
}

fn index() -> &'static Index {
    static INDEX: OnceLock<Index> = OnceLock::new();
    INDEX.get_or_init(|| {
        let mut by_path = HashMap::new();
        let mut by_name: HashMap<&'static str, Vec<RustDoc>> = HashMap::new();
        for record in TABLE.split('\u{1e}') {
            let mut fields = record.splitn(3, '\u{1f}');
            let (Some(path), Some(name), Some(doc)) =
                (fields.next(), fields.next(), fields.next())
            else {
                continue;
            };
            let entry = RustDoc { path, name, doc };
            by_path.insert(path, entry);
            by_name.entry(name).or_default().push(entry);
        }
        Index { by_path, by_name }
    })
}

/// How many documented functions the build harvested.
#[must_use]
pub fn count() -> usize {
    index().by_path.len()
}

/// The documentation of the function at `path` (`crate::a::b::f`).
#[must_use]
pub fn by_path(path: &str) -> Option<RustDoc> {
    index().by_path.get(path).copied()
}

/// Whether `path` is a free function (not a method, not in a test module).
fn is_free_function(path: &str) -> bool {
    let mut segments = path.rsplit("::");
    let _function = segments.next();
    let owner = segments.next().unwrap_or("");
    !owner.starts_with(|c: char| c.is_ascii_uppercase())
        && !path.contains("::tests::")
        && !path.contains("::test::")
}

/// The one documented free function called `name`, narrowed by `hints`
/// (module path fragments such as the rule set and the category) when several
/// share the name. `None` when it is absent or cannot be decided.
#[must_use]
pub fn by_name(
    name: &str,
    hints: &[&str],
) -> Option<RustDoc> {
    let candidates: Vec<RustDoc> = index()
        .by_name
        .get(name)?
        .iter()
        .copied()
        .filter(|entry| is_free_function(entry.path))
        .collect();
    match candidates.as_slice() {
        | [] => None,
        | [one] => Some(*one),
        | many => {
            let narrowed: Vec<RustDoc> = many
                .iter()
                .copied()
                .filter(|entry| {
                    entry.path.split("::").any(|segment| {
                        hints
                            .iter()
                            .any(|hint| !hint.is_empty() && segment == *hint)
                    })
                })
                .collect();
            match narrowed.as_slice() {
                | [one] => Some(*one),
                | _ => None,
            }
        },
    }
}

/// What the manifest should say about an operator's documentation.
pub struct Enriched {
    /// The markdown to publish.
    pub doc: String,
    /// Where it came from.
    pub source: DocSource,
    /// The function whose rustdoc was used, when one was.
    pub path: Option<String>,
    /// Examples recovered from the rustdoc, when the registration gave none.
    pub examples: Vec<Example>,
}

/// Combines what the registration declared with the rustdoc of its function.
#[must_use]
pub fn enrich(desc: &OpDescriptor) -> Enriched {
    let declared = desc.doc.trim();
    let found = desc
        .rust_path
        .as_deref()
        .and_then(by_path)
        .map(|doc| (doc, DocSource::RustPath))
        .or_else(|| {
            let rule_set: &str = &desc.rule_set;
            let category: &str = &desc.category;
            by_name(&desc.name, &[rule_set, category]).map(|doc| (doc, DocSource::RustName))
        });
    let paragraphs = declared
        .lines()
        .filter(|line| !line.trim().is_empty())
        .count();
    let Some((rust, source)) = found else {
        return Enriched {
            doc: declared.to_owned(),
            source: if declared.is_empty() {
                DocSource::Missing
            } else {
                DocSource::Declared
            },
            path: None,
            examples: Vec::new(),
        };
    };
    let examples = if desc.examples.is_empty() {
        examples_from(rust.doc, &desc.name, rust.name)
    } else {
        Vec::new()
    };
    let doc = if declared.is_empty() {
        rust.doc.to_owned()
    } else if paragraphs >= 3 {
        // The registration already wrote a real description.
        return Enriched {
            doc: declared.to_owned(),
            source: DocSource::Declared,
            path: None,
            examples,
        };
    } else if rust.doc.starts_with(declared) {
        rust.doc.to_owned()
    } else if declared.starts_with(rust.doc.lines().next().unwrap_or("")) {
        declared.to_owned()
    } else {
        format!("{declared}\n\n{}", rust.doc)
    };
    Enriched {
        doc,
        source,
        path: Some(rust.path.to_owned()),
        examples,
    }
}

/// `assert_eq!(f(2), 4)` lines of the rustdoc, as worked examples.
///
/// Only a single-line assertion whose left side calls the function is taken;
/// anything else in a doctest is setup, not an example of the call.
fn examples_from(
    doc: &str,
    op: &str,
    function: &str,
) -> Vec<Example> {
    let mut out = Vec::new();
    let mut fenced = false;
    for line in doc.lines() {
        let line = line.trim();
        if line.starts_with("```") {
            fenced = !fenced;
            continue;
        }
        if !fenced {
            continue;
        }
        let Some(body) = line
            .strip_prefix("assert_eq!(")
            .and_then(|rest| rest.strip_suffix(");"))
        else {
            continue;
        };
        let Some((left, right)) = split_top_level(body) else {
            continue;
        };
        let calls = left.contains(&format!("{function}(")) || left.contains(&format!("{op}("));
        if !calls || right.contains("assert") {
            continue;
        }
        out.push(Example {
            input: left.trim().to_owned(),
            output: right.trim().to_owned(),
            note: String::from("from the Rust documentation"),
        });
        if out.len() == 3 {
            break;
        }
    }
    out
}

/// Splits `a, b` at the first comma outside brackets and strings.
fn split_top_level(text: &str) -> Option<(&str, &str)> {
    let mut depth = 0_i32;
    let mut quoted = false;
    let mut first: Option<usize> = None;
    for (index, ch) in text.char_indices() {
        match ch {
            | '"' => quoted = !quoted,
            | _ if quoted => {},
            | '(' | '[' | '{' => depth += 1,
            | ')' | ']' | '}' => depth -= 1,
            | ',' if depth == 0 => {
                if first.is_some() {
                    // A third argument is an assertion message.
                    return first.map(|at| (&text[..at], &text[at + 1..index]));
                }
                first = Some(index);
            },
            | _ => {},
        }
    }
    first.map(|at| (&text[..at], &text[at + 1..]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_assertions_at_the_top_level() {
        assert_eq!(split_top_level("f(1, 2), 3"), Some(("f(1, 2)", " 3")));
        assert_eq!(split_top_level("f(\"a,b\"), 3, \"m\""), Some(("f(\"a,b\")", " 3")));
        assert_eq!(split_top_level("nothing"), None);
    }

    #[test]
    fn the_table_is_populated_and_keyed_by_module_path() {
        assert!(count() > 100, "the build harvested {} docs", count());
        let (path, _) = index()
            .by_path
            .iter()
            .find(|(path, _)| path.starts_with("crate::"))
            .expect("at least one entry");
        assert!(by_path(path).is_some());
    }

    #[test]
    fn examples_come_only_from_calls_to_the_function() {
        let doc = "Square.\n\n```\nlet x = 2;\nassert_eq!(square(x), 4);\nassert_eq!(other(1), 2);\n```";
        let examples = examples_from(doc, "square", "square");
        assert_eq!(examples.len(), 1);
        assert_eq!(examples[0].input, "square(x)");
        assert_eq!(examples[0].output, "4");
    }
}
