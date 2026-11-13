//! Self-describing scenarios: the typed catalog behind [`scenario::run`].
//!
//! [`scenario::run`](super::scenario::run) is a JSON-in, JSON-out entry point.
//! That keeps the C interface small, but a typed host language (Aevia) cannot
//! check a JSON string, and writing a typed wrapper by hand for each of the
//! scenarios is how the two drift apart. This module makes every scenario
//! describe itself, so that wrappers are *generated*:
//!
//! * **Parameters** come from the scenario's own request type: the engine
//!   scenarios serialise their `Default` request ([`scenario_engines`]), the
//!   fixed-geometry ones carry a small runnable example. Names, types and
//!   defaults are *inferred* from that JSON ([`infer`]) — no second
//!   description to keep in sync.
//! * **Results** are *witnessed*: [`witness`] runs each scenario on its
//!   defaults and infers the shape of what came back. The witnessed shapes are
//!   checked in (`scenario_results.json`, like `signatures.witnessed.json` for
//!   operators) and a test fails when a scenario's output changes shape.
//! * **What inference cannot see** is a short annotation table ([`NOTES`]):
//!   the accepted values of a string field (a method name becomes an
//!   enumeration), a physical unit, a one-line description, and the type of
//!   a field whose default is `null` or an empty list. A test runs every
//!   annotated choice and checks that an unannotated value is refused, so the
//!   annotations cannot rot either.
//!
//! [`catalog`] assembles the three; [`catalog_json`] is the form hosts read.
//!
//! [`scenario_engines`]: super::scenario_engines

use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;

/// The type of a parameter or result field.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum FieldType {
    /// `true` / `false`.
    Bool,
    /// An integer.
    Int,
    /// A floating-point number.
    Float,
    /// Free text.
    Str,
    /// One of a fixed set of strings (a method or model name).
    Choice {
        /// The accepted values, in documentation order.
        values: Vec<String>,
    },
    /// A rectangular array of numbers: a list for rank 1, a list of equal
    /// rows for rank 2, ... (also an `ndarray` `{v, dim, data}` object).
    Tensor {
        /// Number of axes.
        rank: usize,
        /// Element type (`Int`, `Float` or `Bool`).
        elem: Box<Self>,
    },
    /// A list whose elements have one type but no common length.
    List {
        /// The element type.
        elem: Box<Self>,
    },
    /// A fixed-length list of differently typed parts.
    Tuple {
        /// The parts.
        items: Vec<Self>,
    },
    /// A JSON object with named fields.
    Record {
        /// The fields.
        fields: Vec<Field>,
    },
    /// A value that may be absent (`null`).
    Optional {
        /// The type when present.
        inner: Box<Self>,
    },
    /// Nothing could be inferred (an empty list, `null` without a note).
    Any,
}

/// A named field of a parameter or result record.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Field {
    /// The JSON key.
    pub name: String,
    /// Its type.
    #[serde(rename = "type")]
    pub ty: FieldType,
    /// The value used when the field is omitted (parameters of engine
    /// scenarios), as JSON text; `None` when the field is required.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub default: Option<String>,
    /// The physical unit the number is in, as unit text (`"m"`, `"K"`).
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub unit: Option<String>,
    /// A one-line description.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub doc: Option<String>,
}

/// The description of one scenario.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ScenarioSchema {
    /// The name [`scenario::run`](super::scenario::run) accepts.
    pub name: String,
    /// What the scenario simulates.
    pub summary: String,
    /// Whether every parameter has a default (`{}` is a valid request).
    pub all_defaults: bool,
    /// The parameters.
    pub params: Vec<Field>,
    /// The fields of the result object.
    pub result: Vec<Field>,
    /// A request that runs (the defaults, or the built-in example).
    pub example: String,
}

/// An annotation of one field: what inference cannot see.
#[derive(Clone, Copy, Debug)]
pub struct Note {
    /// The scenario.
    pub scenario: &'static str,
    /// The parameter (or, prefixed with `->`, the result field).
    pub field: &'static str,
    /// Accepted values of a string field; empty otherwise.
    pub choices: &'static [&'static str],
    /// Unit text; empty when dimensionless or unknown.
    pub unit: &'static str,
    /// Description; may be empty.
    pub doc: &'static str,
    /// A type in the catalog's own notation (`"int?"`, `"float[]"`), for a
    /// field inference cannot type; empty otherwise.
    pub ty: &'static str,
}

const fn note(
    scenario: &'static str,
    field: &'static str,
    choices: &'static [&'static str],
    unit: &'static str,
    doc: &'static str,
) -> Note {
    Note {
        scenario,
        field,
        choices,
        unit,
        doc,
        ty: "",
    }
}

const fn typed(
    scenario: &'static str,
    field: &'static str,
    ty: &'static str,
    doc: &'static str,
) -> Note {
    Note {
        scenario,
        field,
        choices: &[],
        unit: "",
        doc,
        ty,
    }
}

const BC: &[&str] = &["dirichlet", "neumann", "periodic"];

/// The annotation table (see the module documentation).
pub const NOTES: &[Note] = &[
    note("poisson_multigrid", "cycle", &["v", "w", "fmg"], "", "multigrid cycle"),
    note("poisson_multigrid", "dim", &[], "", "space dimension, 1 to 3"),
    note("heat_nd", "bc", BC, "", "boundary condition on every face"),
    note("heat_nd", "dim", &[], "", "space dimension, 1 to 3"),
    note(
        "euler_fv",
        "recon",
        &[
            "weno5",
            "first_order",
            "muscl_minmod",
            "muscl_vanleer",
            "muscl_mc",
            "muscl_superbee",
        ],
        "",
        "reconstruction of interface states",
    ),
    note("euler_fv", "riemann", &["hllc", "roe", "rusanov"], "", "approximate Riemann solver"),
    note("dg_conservation", "law", &["advection", "burgers", "euler"], "", "conservation law"),
    note("navier_stokes_mac", "case", &["cavity", "taylor_green"], "", "test problem"),
    note("split_step", "mode", &["ground_state", "evolve"], "", "imaginary-time ground state or real-time evolution"),
    note("bem_laplace", "problem", &["interior", "exterior"], "", "boundary-value problem"),
    note("bem_laplace", "method", &["aca", "p1"], "", "adaptive cross approximation (P0) or dense P1 elements"),
    note("md_lj", "thermostat", &["nve", "langevin"], "", "ensemble"),
    note("spin_model", "model", &["ising", "xy", "heisenberg"], "", "spin model"),
    note(
        "spin_model",
        "update",
        &["wolff", "metropolis", "heat_bath", "swendsen_wang"],
        "",
        "Monte Carlo update",
    ),
    note("sod_shock_tube", "limiter", &["van_leer", "minmod", "none"], "", "slope limiter"),
    typed("sod_shock_tube", "limiter", "choice?", "slope limiter (default van_leer)"),
    typed("ising", "seed", "int?", "random seed (a fixed default when absent)"),
    typed("ising_wolff", "seed", "int?", "random seed (a fixed default when absent)"),
    typed("kepler_symplectic", "order", "int?", "integrator order, 2 or 4 (default 4)"),
    typed("heat_nd", "adaptive_tol", "float?", "local error tolerance of adaptive time steps (fixed steps when absent)"),
    typed("navier_stokes_mac", "->taylor_green_decay", "float?", "measured decay rate (Taylor-Green case only)"),
    note("ising", "temperature", &[], "", "temperature in units of J/k_B"),
    note("geodesic", "black_hole_mass", &[], "", "mass in geometrized units"),
];

/// Infers the type of a JSON value.
#[must_use]
pub fn infer(value: &Value) -> FieldType {
    match value {
        | Value::Null => FieldType::Optional {
            inner: Box::new(FieldType::Any),
        },
        | Value::Bool(_) => FieldType::Bool,
        | Value::Number(n) => {
            if n.is_f64() {
                FieldType::Float
            } else {
                FieldType::Int
            }
        },
        | Value::String(_) => FieldType::Str,
        | Value::Array(items) => infer_array(items),
        | Value::Object(map) => {
            // `ndarray`'s encoding: {"v": 1, "dim": [..], "data": [..]}.
            if let (Some(Value::Array(dim)), Some(Value::Array(data))) = (map.get("dim"), map.get("data"))
                && map.contains_key("v")
            {
                let elem = scalar_join(data.iter().map(infer)).unwrap_or(FieldType::Any);
                return match elem {
                    // Complex entries serialise as [re, im] pairs.
                    | FieldType::Tensor { rank, elem } => FieldType::Tensor {
                        rank: dim.len() + rank,
                        elem,
                    },
                    | elem => FieldType::Tensor {
                        rank: dim.len(),
                        elem: Box::new(elem),
                    },
                };
            }
            FieldType::Record {
                fields: map
                    .iter()
                    .map(|(k, v)| {
                        Field {
                            name: k.clone(),
                            ty: infer(v),
                            default: None,
                            unit: None,
                            doc: None,
                        }
                    })
                    .collect(),
            }
        },
    }
}

/// The common scalar type of several values (`Int` and `Float` join to
/// `Float`), or a common tensor type of equal shape.
fn scalar_join(mut types: impl Iterator<Item = FieldType>) -> Option<FieldType> {
    let first = types.next()?;
    types.try_fold(first, |acc, t| {
        match (acc, t) {
            | (a, b) if a == b => Some(a),
            | (FieldType::Int, FieldType::Float) | (FieldType::Float, FieldType::Int) => {
                Some(FieldType::Float)
            },
            | _ => None,
        }
    })
}

fn infer_array(items: &[Value]) -> FieldType {
    if items.is_empty() {
        return FieldType::List {
            elem: Box::new(FieldType::Any),
        };
    }
    let types: Vec<FieldType> = items.iter().map(infer).collect();
    // A rectangular array of numbers (or booleans) is a tensor: the rows
    // must agree on their length as well as their type.
    let lengths: Vec<Option<usize>> = items
        .iter()
        .map(|v| v.as_array().map(Vec::len))
        .collect();
    let rectangular = lengths.windows(2).all(|w| w[0] == w[1]);
    if let Some(joined) = scalar_join(types.iter().cloned()) {
        return match joined {
            | FieldType::Int | FieldType::Float | FieldType::Bool => FieldType::Tensor {
                rank: 1,
                elem: Box::new(joined),
            },
            | FieldType::Tensor { rank, elem } if rectangular => FieldType::Tensor {
                rank: rank + 1,
                elem,
            },
            | other => FieldType::List {
                elem: Box::new(other),
            },
        };
    }
    // Rows of different lengths, or different element types: a ragged list
    // when the element types agree up to length, a tuple otherwise.
    let ragged = types.iter().all(|t| matches!(t, FieldType::Tensor { .. }))
        && types.windows(2).all(|w| {
            matches!((&w[0], &w[1]), (FieldType::Tensor { rank: a, elem: x }, FieldType::Tensor { rank: b, elem: y }) if a == b && x == y)
        });
    if ragged {
        return FieldType::List {
            elem: Box::new(types[0].clone()),
        };
    }
    FieldType::Tuple { items: types }
}

/// Parses the catalog's type notation of [`Note::ty`]: `bool`, `int`,
/// `float`, `str`, `choice`, a trailing `[]` per tensor axis, a trailing `?`
/// for an optional value.
#[must_use]
pub fn parse_type(text: &str) -> Option<FieldType> {
    if let Some(inner) = text.strip_suffix('?') {
        return Some(FieldType::Optional {
            inner: Box::new(parse_type(inner)?),
        });
    }
    let mut rank = 0;
    let mut base = text;
    while let Some(rest) = base.strip_suffix("[]") {
        rank += 1;
        base = rest;
    }
    let scalar = match base {
        | "bool" => FieldType::Bool,
        | "int" => FieldType::Int,
        | "float" => FieldType::Float,
        | "str" => FieldType::Str,
        | "choice" => FieldType::Choice { values: Vec::new() },
        | _ => return None,
    };
    Some(if rank == 0 {
        scalar
    } else {
        FieldType::Tensor {
            rank,
            elem: Box::new(scalar),
        }
    })
}

/// One-line summaries, in the order of [`SCENARIOS`](super::scenario::SCENARIOS).
const SUMMARIES: &[(&str, &str)] = &[
    ("fdtd", "2D finite-difference time-domain electromagnetics with a point source"),
    ("geodesic", "Timelike geodesic around a Schwarzschild black hole"),
    ("gpe", "Gross-Pitaevskii ground state in a harmonic trap (imaginary time)"),
    ("ising", "2D Ising model, Metropolis Monte Carlo"),
    ("elasticity", "Plane linear elasticity on a quadrilateral mesh"),
    ("lid_driven_cavity", "Incompressible lid-driven cavity (projection method)"),
    ("channel_flow", "Incompressible channel flow past an obstacle mask"),
    ("schrodinger", "Time-dependent 2D Schroedinger equation (split operator)"),
    ("ising_wolff", "2D Ising model, Wolff cluster updates, with statistics"),
    ("sod_shock_tube", "Sod shock tube: HLLC/MUSCL against the exact Riemann solution"),
    ("kepler_symplectic", "Kepler orbit with a symplectic integrator"),
    ("poisson_multigrid", "Poisson equation on a unit box by geometric multigrid, 1D-3D"),
    ("heat_nd", "Heat equation on a unit box, 1D-3D"),
    ("wave_nd", "Wave equation on a unit box, 1D-3D (leapfrog)"),
    ("euler_fv", "Compressible Euler equations, finite volumes"),
    ("dg_conservation", "Discontinuous Galerkin for a 1D conservation law"),
    ("acoustics_pml", "Acoustic waves with a perfectly matched layer"),
    ("yee_fdtd", "Maxwell's equations on a Yee grid"),
    ("navier_stokes_mac", "Incompressible Navier-Stokes on a MAC grid"),
    ("split_step", "Nonlinear Schroedinger / GPE by split-step Fourier"),
    ("fem_simplex", "Poisson equation by P1 finite elements on simplices"),
    ("elasticity_fem", "Linear elasticity by finite elements"),
    ("bem_laplace", "Laplace equation by the boundary element method"),
    ("sph", "Smoothed-particle hydrodynamics (dam break)"),
    ("dpd", "Dissipative particle dynamics"),
    ("md_lj", "Lennard-Jones molecular dynamics with cell lists"),
    ("spin_model", "Ising / XY / Heisenberg lattice spin models, 1D-3D"),
    ("nbody_gravity", "Gravitational N-body problem"),
    ("kerr_geodesic", "Geodesics in the Kerr spacetime"),
];

/// Runnable examples of the scenarios without defaults (small, so that
/// witnessing their results is quick).
fn example(name: &str) -> Option<Value> {
    use serde_json::json;
    Some(match name {
        | "fdtd" => json!({"width": 16, "height": 16, "time_steps": 8, "source_pos": [8, 8], "source_freq": 0.1}),
        | "geodesic" => json!({"black_hole_mass": 1.0, "initial_state": [10.0, 0.0, 0.0, 0.035], "proper_time_end": 5.0, "initial_dt": 0.01}),
        | "gpe" => json!({"nx": 16, "ny": 16, "lx": 10.0, "ly": 10.0, "d_tau": 0.01, "time_steps": 5, "g": 1.0, "trap_strength": 1.0}),
        | "ising" => json!({"width": 8, "height": 8, "temperature": 2.0, "mc_steps": 10, "seed": 1}),
        | "elasticity" => json!({
            "nodes": [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
            "elements": [[0, 1, 2, 3]],
            "youngs_modulus": 1000.0,
            "poissons_ratio": 0.3,
            "fixed_nodes": [0, 3],
            "loads": [[1, 1.0, 0.0], [2, 1.0, 0.0]],
        }),
        | "lid_driven_cavity" => json!({"nx": 8, "ny": 8, "re": 100.0, "dt": 0.01, "n_iter": 5, "lid_velocity": 1.0}),
        | "channel_flow" => json!({"nx": 8, "ny": 8, "re": 50.0, "dt": 0.01, "n_iter": 5,
            "obstacle_mask": {"v": 1, "dim": [8, 8], "data": vec![false; 64]}}),
        | "schrodinger" => {
            let psi: Vec<[f64; 2]> = (0..64).map(|i| [if i == 27 { 1.0 } else { 0.0 }, 0.0]).collect();
            json!({"params": {"nx": 8, "ny": 8, "lx": 1.0, "ly": 1.0, "dt": 0.001, "time_steps": 3,
                "hbar": 1.0, "mass": 1.0, "potential": vec![0.0; 64]}, "initial_psi": psi})
        },
        | "ising_wolff" => json!({"width": 8, "height": 8, "temperature": 1.5, "thermalisation": 5, "sweeps": 10, "seed": 3}),
        | "sod_shock_tube" => json!({"cells": 32, "t": 0.1, "limiter": "van_leer"}),
        | "kepler_symplectic" => json!({"q": [1.0, 0.0], "p": [0.0, 1.0], "dt": 0.01, "steps": 10, "order": 4}),
        | "nbody_gravity" => json!({"positions": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            "masses": [1.0, 0.5, 0.25], "g": 1.0, "softening": 1e-3, "theta": 0.5}),
        | _ => return None,
    })
}

/// The request a scenario is described and witnessed with, and whether it
/// is the scenario's own defaults.
fn request(name: &str) -> Option<(Value, bool)> {
    if let Some(defaults) = super::scenario_engines::defaults(name) {
        return Some((defaults, true));
    }
    example(name).map(|e| (e, false))
}

fn notes_for<'a>(
    scenario: &'a str,
    field: &'a str,
) -> impl Iterator<Item = &'static Note> + 'a {
    NOTES
        .iter()
        .filter(move |n| n.scenario == scenario && n.field == field)
}

fn annotate(
    scenario: &str,
    key: &str,
    field: &mut Field,
) {
    for n in notes_for(scenario, key) {
        if !n.ty.is_empty()
            && let Some(ty) = parse_type(n.ty)
        {
            field.ty = ty;
        }
        if !n.choices.is_empty() {
            let values: Vec<String> = n.choices.iter().map(|c| (*c).to_owned()).collect();
            field.ty = match std::mem::replace(&mut field.ty, FieldType::Any) {
                | FieldType::Optional { .. } => FieldType::Optional {
                    inner: Box::new(FieldType::Choice { values }),
                },
                | _ => FieldType::Choice { values },
            };
        }
        if !n.unit.is_empty() {
            field.unit = Some(n.unit.to_owned());
        }
        if !n.doc.is_empty() {
            field.doc = Some(n.doc.to_owned());
        }
    }
}

fn fields_of(
    scenario: &str,
    value: &Value,
    defaults: bool,
    prefix: &str,
) -> Vec<Field> {
    let Value::Object(map) = value else {
        return Vec::new();
    };
    map.iter()
        .map(|(key, v)| {
            let mut field = Field {
                name: key.clone(),
                ty: infer(v),
                default: defaults.then(|| v.to_string()),
                unit: None,
                doc: None,
            };
            annotate(scenario, &format!("{prefix}{key}"), &mut field);
            // Choices also apply to the choice inside an optional field.
            if let FieldType::Optional { inner } = &mut field.ty
                && let FieldType::Choice { values } = inner.as_mut()
                && values.is_empty()
                && let Some(n) = notes_for(scenario, key).find(|n| !n.choices.is_empty())
            {
                *values = n.choices.iter().map(|c| (*c).to_owned()).collect();
            }
            field
        })
        .collect()
}

/// The witnessed result shapes, keyed by scenario: `{name: {field: type}}`.
const WITNESSED: &str = include_str!("scenario_results.json");

/// Runs every scenario on its request and infers the shape of its result
/// (the content of `scenario_results.json`).
///
/// # Errors
/// A scenario failed on its own defaults or example.
pub fn witness() -> Result<Value, String> {
    let mut out = serde_json::Map::new();
    for &name in super::scenario::SCENARIOS {
        let (params, _) = request(name).ok_or_else(|| format!("no request for `{name}`"))?;
        let text = super::scenario::run(name, &params.to_string())
            .map_err(|e| format!("`{name}` failed on its request: {e}"))?;
        let result: Value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
        let fields = fields_of(name, &result, false, "->");
        out.insert(
            name.to_owned(),
            serde_json::to_value(fields).map_err(|e| e.to_string())?,
        );
    }
    Ok(Value::Object(out))
}

/// The catalog of every scenario (see the module documentation).
#[must_use]
pub fn catalog() -> Vec<ScenarioSchema> {
    let witnessed: serde_json::Map<String, Value> =
        serde_json::from_str(WITNESSED).unwrap_or_default();
    super::scenario::SCENARIOS
        .iter()
        .filter_map(|&name| {
            let (params, all_defaults) = request(name)?;
            let mut result: Vec<Field> = witnessed
                .get(name)
                .and_then(|v| serde_json::from_value(v.clone()).ok())
                .unwrap_or_default();
            for field in &mut result {
                annotate(name, &format!("->{}", field.name), field);
            }
            Some(ScenarioSchema {
                name: name.to_owned(),
                summary: SUMMARIES
                    .iter()
                    .find(|(n, _)| *n == name)
                    .map_or_else(String::new, |(_, s)| (*s).to_owned()),
                all_defaults,
                params: fields_of(name, &params, all_defaults, ""),
                result,
                example: params.to_string(),
            })
        })
        .collect()
}

/// [`catalog`] as JSON text, the form a host reads.
#[must_use]
pub fn catalog_json() -> String {
    serde_json::to_string(&catalog()).unwrap_or_else(|_| "[]".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inference_recognises_tensors_tuples_and_ndarrays() {
        use serde_json::json;
        assert_eq!(
            infer(&json!([[1.0, 2.0], [3, 4.5]])),
            FieldType::Tensor {
                rank: 2,
                elem: Box::new(FieldType::Float)
            }
        );
        assert_eq!(
            infer(&json!([[1, 2], [3]])),
            FieldType::List {
                elem: Box::new(FieldType::Tensor {
                    rank: 1,
                    elem: Box::new(FieldType::Int)
                })
            }
        );
        assert_eq!(
            infer(&json!([1, 1.5, "x"])),
            FieldType::Tuple {
                items: vec![FieldType::Int, FieldType::Float, FieldType::Str]
            }
        );
        assert_eq!(
            infer(&json!({"v": 1, "dim": [2, 3], "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]})),
            FieldType::Tensor {
                rank: 2,
                elem: Box::new(FieldType::Float)
            }
        );
        assert_eq!(parse_type("float[][]?"), Some(FieldType::Optional {
            inner: Box::new(FieldType::Tensor {
                rank: 2,
                elem: Box::new(FieldType::Float)
            })
        }));
    }

    #[test]
    fn every_scenario_is_described() {
        let catalog = catalog();
        assert_eq!(catalog.len(), super::super::scenario::SCENARIOS.len());
        for s in &catalog {
            assert!(!s.summary.is_empty(), "{} has no summary", s.name);
            assert!(!s.params.is_empty(), "{} has no parameters", s.name);
            assert!(!s.result.is_empty(), "{} has no witnessed result", s.name);
        }
        let poisson = catalog
            .iter()
            .find(|s| s.name == "poisson_multigrid")
            .unwrap_or_else(|| panic!("listed"));
        assert!(poisson.all_defaults);
        let cycle = poisson
            .params
            .iter()
            .find(|f| f.name == "cycle")
            .unwrap_or_else(|| panic!("has a cycle"));
        assert_eq!(cycle.default.as_deref(), Some("\"v\""));
        assert!(matches!(&cycle.ty, FieldType::Choice { values } if values.len() == 3));
    }

    /// Regenerate with `RSSN_BLESS=1 cargo test --lib sim::schema`.
    #[test]
    fn witnessed_results_are_current() {
        let fresh = witness().unwrap_or_else(|e| panic!("{e}"));
        let text = serde_json::to_string_pretty(&fresh).unwrap_or_default();
        if std::env::var_os("RSSN_BLESS").is_some() {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/sim/scenario_results.json");
            std::fs::write(path, format!("{text}\n")).unwrap_or_else(|e| panic!("{e}"));
            return;
        }
        let stored: Value = serde_json::from_str(WITNESSED).unwrap_or(Value::Null);
        assert_eq!(stored, fresh, "a scenario's result changed shape; re-bless scenario_results.json");
    }

    #[test]
    fn annotated_choices_run_and_others_are_refused() {
        for n in NOTES.iter().filter(|n| !n.choices.is_empty()) {
            let Some((base, _)) = request(n.scenario) else {
                panic!("{}: no request", n.scenario)
            };
            for &choice in n.choices {
                let mut params = base.clone();
                params[n.field] = Value::String(choice.to_owned());
                // The value must be *recognised*: it may still fail on
                // defaults chosen for another method (the Euler law on the
                // advection defaults), but not as an unknown name.
                let outcome = super::super::scenario::run(n.scenario, &params.to_string());
                if let Err(e) = &outcome {
                    assert!(
                        !e.contains("unknown") && !e.contains(choice),
                        "{}.{} = {choice} is not recognised: {e}",
                        n.scenario,
                        n.field
                    );
                }
            }
            let mut params = base.clone();
            params[n.field] = Value::String("not-a-choice".to_owned());
            assert!(
                super::super::scenario::run(n.scenario, &params.to_string()).is_err(),
                "{}.{} accepts an unlisted value",
                n.scenario,
                n.field
            );
        }
    }
}
