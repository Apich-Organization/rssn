//! Real JSON scenario dispatch plus the native schema API. JSON numbers retain
//! the units documented by each native scenario schema; measured language values
//! use the separate physical adapters instead of losing their unit witnesses.
use super::codecs::Arg;
use super::native_object_api::Dynamic;
use super::native_object_api::Input;
use super::native_object_api::object_binding;
use super::native_object_api::{
    self as api,
};
use crate::graph::RuleError;
use crate::graph::rule::Installer;
use crate::kernel::CallCtx;
use crate::kernel::OpError;
use crate::kernel::Value;
use crate::kernel::signature::TypeExpr;
use crate::kernel::value::TensorValue;
use crate::kernel::value::VariantValue;
use crate::sim::scenario;
use crate::sim::schema;
use serde_json::Value as Json;
fn json(text: &str) -> Result<Json, OpError> {
    if text.len() > 1_048_576 {
        return Err(api::error("native JSON text budget exceeded"));
    }
    let value: Json = serde_json::from_str(text).map_err(api::error)?;
    let mut stack = vec![(&value, 0usize)];
    let mut count = 0usize;
    while let Some((v, depth)) = stack.pop() {
        count += 1;
        if count > 262144 || depth > 32 {
            return Err(api::error("native JSON node/depth budget exceeded"));
        }
        match v {
            | Json::Array(values) => {
                if values.len() > 65536 {
                    return Err(api::error("native JSON array allocation bound exceeded"));
                }
                stack.extend(values.iter().map(|v| (v, depth + 1)));
            },
            | Json::Object(values) => stack.extend(values.values().map(|v| (v, depth + 1))),
            | Json::Number(n) => {
                if n.as_f64().is_none_or(|v| !v.is_finite()) {
                    return Err(api::error("native JSON number is not finite"));
                }
            },
            | _ => {},
        }
    }
    Ok(value)
}
fn text(value: String) -> Result<String, OpError> {
    if value.len() > 1_048_576 {
        return Err(api::error("native JSON output text budget exceeded"));
    }
    Ok(value)
}
fn integer(
    v: &Json,
    name: &str,
) -> Result<Option<usize>, OpError> {
    v.get(name)
        .map(|v| {
            v.as_u64()
                .and_then(|n| usize::try_from(n).ok())
                .ok_or_else(|| {
                    api::error(format!(
                        "native scenario {name} must be an unsigned machine integer"
                    ))
                })
        })
        .transpose()
}
fn product(
    a: usize,
    b: usize,
) -> Result<usize, OpError> {
    a.checked_mul(b)
        .filter(|n| *n <= 50_000_000)
        .ok_or_else(|| api::error("native scenario allocation/work product exceeded"))
}
fn validate(
    name: &str,
    parameters: &str,
    budget: usize,
) -> Result<usize, OpError> {
    if budget == 0 || budget > 50_000_000 {
        return Err(api::error("scenario work budget must be1..50million"));
    }
    let request = json(parameters)?;
    let Json::Object(user) = &request else {
        return Err(api::error(
            "native scenario parameters must be a JSON object",
        ));
    };
    let descriptor = schema::catalog()
        .into_iter()
        .find(|entry| entry.name == name)
        .ok_or_else(|| api::error("unknown native scenario"))?;
    let mut expanded = if descriptor.all_defaults {
        json(&descriptor.example)?
    } else {
        Json::Object(serde_json::Map::new())
    };
    let Json::Object(map) = &mut expanded else {
        return Err(api::error("native default request is not an object"));
    };
    for (key, value) in user {
        if !descriptor.params.iter().any(|p| p.name == *key) {
            return Err(api::error(format!(
                "unknown native scenario parameter {key}"
            )));
        }
        map.insert(key.clone(), value.clone());
    }
    let p = expanded.get("params").unwrap_or(&expanded);
    let d = integer(p, "dim")?.unwrap_or(2);
    if !(1..=3).contains(&d) {
        return Err(api::error("scenario dimension must be1..3"));
    }
    let mut cells = 1usize;
    for key in ["nx", "ny", "nz", "width", "height"] {
        if let Some(n) = integer(p, key)? {
            if n == 0 || n > 4096 {
                return Err(api::error("scenario grid side must be1..4096"));
            }
            cells = product(cells, n)?;
        }
    }
    if let Some(n) = integer(p, "n")? {
        if n == 0 || n > 4096 {
            return Err(api::error("scenario n must be1..4096"));
        }
        cells = product(
            cells,
            n.checked_pow(d as u32)
                .ok_or_else(|| api::error("scenario shape overflows"))?,
        )?;
    }
    if let Some(n) = integer(p, "cells")? {
        cells = product(cells, n)?;
    }
    if let Some(positions) = p.get("positions").and_then(Json::as_array) {
        if positions.is_empty() || positions.len() > 512 {
            return Err(api::error("scenario source count must be1..512"));
        }
        cells = cells.max(product(positions.len(), positions.len())?);
        if p.get("masses")
            .and_then(Json::as_array)
            .is_none_or(|m| m.len() != positions.len())
        {
            return Err(api::error("scenario mass/source counts differ"));
        }
    }
    if name == "bem_laplace" {
        let subdivisions = integer(p, "subdivisions")?.unwrap_or(0);
        if subdivisions > 2 {
            return Err(api::error("BEM scenario subdivision budget is2"));
        }
        cells = product(20, 4usize.pow(subdivisions as u32))?;
        cells = product(cells, cells)?;
    }
    if matches!(name, "fem_simplex" | "elasticity_fem") {
        let degree = integer(p, "degree")?.unwrap_or(1);
        if !(1..=2).contains(&degree) {
            return Err(api::error("FEM scenario degree must be1or2"));
        }
        cells = product(cells, 6 * degree.pow(d as u32))?;
        if cells > 8192 {
            return Err(api::error("FEM scenario DOF estimate exceeds8192"));
        }
    }
    if cells > 262144 {
        return Err(api::error(
            "scenario grid allocation estimate exceeds262144",
        ));
    }
    let mut rounds = 1usize;
    for key in [
        "steps",
        "time_steps",
        "n_iter",
        "mc_steps",
        "thermalisation",
        "sweeps",
    ] {
        if let Some(n) = integer(p, key)? {
            if n > 100000 {
                return Err(api::error("scenario iteration count exceeds100000"));
            }
            rounds = rounds
                .checked_add(n)
                .ok_or_else(|| api::error("scenario iteration count overflows"))?;
        }
    }
    for key in ["dt", "d_tau", "initial_dt", "cfl"] {
        if let Some(v) = p.get(key) {
            let x = v
                .as_f64()
                .ok_or_else(|| api::error("native scenario control must be numeric"))?;
            if x <= 0.0 {
                return Err(api::error(format!(
                    "native scenario {key} must be positive"
                )));
            }
        }
    }
    if let (Some(end), Some(dt)) = (
        p.get("proper_time_end").and_then(Json::as_f64),
        p.get("initial_dt").and_then(Json::as_f64),
    ) {
        if end < 0.0 || end / dt > 100000.0 {
            return Err(api::error("geodesic scenario nominal step budget exceeded"));
        }
        rounds = rounds.max((end / dt).ceil() as usize);
    }
    if let Some(t) = p.get("temperature").and_then(Json::as_f64)
        && t < 0.0
    {
        return Err(api::error("scenario temperature must be nonnegative"));
    }
    if name == "elasticity"
        && let Some(nodes) = p.get("nodes").and_then(Json::as_array)
    {
        if nodes.is_empty() || nodes.len() > 512 {
            return Err(api::error("elasticity scenario node budget is512"));
        }
        cells = product(nodes.len(), nodes.len())?;
    }
    if name == "schrodinger" {
        let count = integer(p, "nx")?
            .and_then(|nx| {
                integer(p, "ny")
                    .ok()
                    .flatten()
                    .and_then(|ny| nx.checked_mul(ny))
            })
            .ok_or_else(|| api::error("Schrodinger grid shape missing/overflow"))?;
        if p.get("potential")
            .and_then(Json::as_array)
            .is_none_or(|v| v.len() != count)
            || expanded
                .get("initial_psi")
                .and_then(Json::as_array)
                .is_none_or(|v| v.len() != count)
        {
            return Err(api::error(
                "Schrodinger potential/psi arrays must match grid",
            ));
        }
    }
    for key in ["t", "t_end", "proper_time_end"] {
        if let Some(v) = p.get(key).and_then(Json::as_f64)
            && !(0.0..=10000.0).contains(&v)
        {
            return Err(api::error("scenario duration must be0..10000nativeunits"));
        }
    }
    if name == "sod_shock_tube"
        && (integer(p, "cells")?.unwrap_or(0) > 4096
            || p.get("t").and_then(Json::as_f64).is_none_or(|t| t > 10.0))
    {
        return Err(api::error(
            "Sod scenario needs<=4096cells and nativefinaltime<=10",
        ));
    }
    if let Some(tol) = p.get("adaptive_tol").and_then(Json::as_f64)
        && tol < 1e-12
    {
        return Err(api::error("scenario adaptive tolerance below1e-12budget"));
    }
    let work = product(cells, rounds)?;
    if work > budget {
        return Err(api::error(
            "native scenario request exceeds explicit work budget",
        ));
    }
    Ok(work)
}
object_binding!(infer,"num_native_schema_infer","crate::sim::schema::infer","Infer the actual native scenario field type from a bounded JSON value; result is the native serialized FieldType enum.",true,(json_value:String)->String,|_ctx:&CallCtx<'_>,json_value:String|{let value=json(&json_value)?;text(serde_json::to_string(&schema::infer(&value)).map_err(api::error)?)});
object_binding!(parse_type,"num_native_schema_parse_type","crate::sim::schema::parse_type","Parse native schema type syntax and return its serialized native FieldType enum if valid.",true,(name:String)->Option<String>,|_ctx:&CallCtx<'_>,name:String|{if name.len()>4096{return Err(api::error("schema type text budget exceeded"));}schema::parse_type(&name).map(|v|serde_json::to_string(&v).map_err(api::error)).transpose()});
object_binding!(catalog,"num_native_schema_catalog","crate::sim::schema::catalog","Serialize the actual native scenario catalog, including physical unit notes and defaults; this discovery operator does not substitute for scenario execution.",true,()->String,|_ctx:&CallCtx<'_>|text(serde_json::to_string(&schema::catalog()).map_err(api::error)?));
object_binding!(catalog_json,"num_native_schema_catalog_json","crate::sim::schema::catalog_json","Read the native scenario catalog JSON including parameter defaults and unit notes.",true,()->String,|_ctx:&CallCtx<'_>|text(schema::catalog_json()));
object_binding!(witness,"num_native_schema_witness","crate::sim::schema::witness","Execute every actual native scenario on its declared example/default request and infer output schemas. All examples are preflighted against one explicit aggregate work budget before execution.",true,(work_budget:usize)->String,|_ctx:&CallCtx<'_>,work_budget:usize|{let mut remaining=work_budget;for entry in schema::catalog(){let cost=validate(&entry.name,&entry.example,remaining)?;remaining=remaining.checked_sub(cost).ok_or_else(||api::error("schema witness aggregate work budget exceeded"))?;}text(serde_json::to_string(&crate::sim::models_geodesic_relativity::with_native_budget(work_budget.min(100000),schema::witness).map_err(api::error)?.map_err(api::error)?).map_err(api::error)?)});
object_binding!(run,"num_native_scenario_run","crate::sim::scenario::run","Execute the actual native named scenario with bounded JSON parameters. The contract must be native-schema: JSON numbers use the units documented by that scenario's schema and preserve native defaults. Typed physical quantity APIs are separately exposed for each numerical engine.",true,(name:String,parameters:String,parameter_contract:String,work_budget:usize)->String,|_ctx:&CallCtx<'_>,name:String,parameters:String,parameter_contract:String,work_budget:usize|{if parameter_contract!="native-schema"{return Err(api::error("JSON scenario parameters require the explicit native-schema unit contract"));}validate(&name,&parameters,work_budget)?;let out=crate::sim::models_geodesic_relativity::with_native_budget(work_budget.min(100000),||scenario::run(&name,&parameters)).map_err(api::error)?.map_err(api::error)?;text(out)});
object_binding!(run_many,"num_native_scenario_run_many","crate::sim::scenario::run_many","Execute actual native parallel run_many over bounded named JSON requests. Input records have name and parameters strings; per-request native Ok/Err results and request order are preserved, with invalid requests rejected before native allocation.",true,(requests:Input,parameter_contract:String,work_budget:usize)->Dynamic,|_ctx:&CallCtx<'_>,requests:Input,parameter_contract:String,work_budget:usize|{if parameter_contract!="native-schema"{return Err(api::error("JSON scenarios require the explicit native-schema unit contract"));}let Value::Tensor(t)=&requests.0 else{return Err(api::error("scenario request vector required"));};
if t.shape!=[t.data.len()]||t.data.len()>32{return Err(api::error("scenario batch shape/count exceeds32"));}let mut accepted=vec![];let mut slots:Vec<Option<Result<String,String>>>=vec![None;t.data.len()];let mut remaining=work_budget;for(i,v)in t.data.iter().enumerate(){let name=<String as Arg>::decode(api::field(v,"name")?,0)?;let params=<String as Arg>::decode(api::field(v,"parameters")?,0)?;match validate(&name,&params,remaining){Ok(cost)=>{remaining-=cost;accepted.push((i,name,params));},Err(e)=>slots[i]=Some(Err(e.to_string()))}}let borrowed=accepted.iter().map(|(_,n,p)|(n.as_str(),p.as_str())).collect::<Vec<_>>();let outputs=crate::sim::models_geodesic_relativity::with_native_budget(work_budget.min(100000),||scenario::run_many(&borrowed)).map_err(api::error)?;for((index,_,_),value)in accepted.into_iter().zip(outputs){slots[index]=Some(value);}let data=slots.into_iter().map(|v|{let result=v.ok_or_else(||api::error("native scenario batch lost a result"))?;let result=match result{Ok(s)=>text(s).map_err(|e|e.to_string()),Err(e)=>Err(e)};Ok(Value::Variant(match result{Ok(v)=>VariantValue{tag:"Ok".into(),payload:vec![Value::String(v)]},Err(e)=>VariantValue{tag:"Err".into(),payload:vec![Value::String(e)]}}))}).collect::<Result<Vec<_>,OpError>>()?;let elem=TypeExpr::Variant(vec![("Ok".into(),vec![TypeExpr::String]),("Err".into(),vec![TypeExpr::String])]);Ok(Dynamic(Value::Tensor(TensorValue{shape:vec![data.len()],data,elem})))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &infer::BINDING,
            &parse_type::BINDING,
            &catalog::BINDING,
            &catalog_json::BINDING,
            &witness::BINDING,
            &run::BINDING,
            &run_many::BINDING,
        ],
        "numerical_vector",
    )
}
