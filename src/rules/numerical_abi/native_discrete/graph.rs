#![allow(clippy::wildcard_imports)] // Shared checked graph method declarations.
use super::*;
use crate::kernels::graph::CanonicalLabelling;
use crate::kernels::graph::CellCertificate;
use crate::kernels::graph::ColourGraph;
use crate::kernels::graph::Graph;
use crate::kernels::graph::Permutation;
use crate::kernels::graph::Refinement;
use crate::kernels::graph::RefinementCertificate;
use crate::kernels::graph::SearchBudget;
impl NativeType for Graph {
    const KIND: &'static str = "rssn/native/discrete/graph-graph/v1";
    const MAX_BYTES: usize = 2097152;

    fn validate(&self) -> Result<(), OpError> {
        if self.num_nodes() > 256 {
            return Err(error("weighted graph vertex budget exceeded"));
        }
        let mut edges = 0usize;
        for u in 0..self.num_nodes() {
            edges += self.adj(u).len();
            if self
                .adj(u)
                .iter()
                .any(|(v, w)| *v >= self.num_nodes() || !w.is_finite())
            {
                return Err(error("invalid weighted graph edge"));
            }
        }
        limit(edges, 65536)?;
        Ok(())
    }
}
impl NativeType for ColourGraph {
    const KIND: &'static str = "rssn/native/discrete/graph-colourgraph/v1";
    const MAX_BYTES: usize = 262144;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.num_nodes(), 64)?;
        Ok(())
    }
}
impl NativeType for Permutation {
    const KIND: &'static str = "rssn/native/discrete/graph-permutation/v1";
    const MAX_BYTES: usize = 4096;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.len(), 256)?;
        Ok(())
    }
}
impl NativeType for SearchBudget {
    const KIND: &'static str = "rssn/native/discrete/graph-searchbudget/v1";
    const MAX_BYTES: usize = 64;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.max_nodes(), 1_000_000)?;
        Ok(())
    }
}
impl NativeType for CanonicalLabelling {
    const KIND: &'static str = "rssn/native/discrete/graph-canonicallabelling/v1";
    const MAX_BYTES: usize = 4194304;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.num_nodes(), 64)?;
        Ok(())
    }
}
impl NativeType for Refinement {
    const KIND: &'static str = "rssn/native/discrete/graph-refinement/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.cell_of().len(), 64)?;
        Ok(())
    }
}
impl NativeType for RefinementCertificate {
    const KIND: &'static str = "rssn/native/discrete/graph-refinementcertificate/v1";
    const MAX_BYTES: usize = 1048576;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.num_nodes(), 64)?;
        Ok(())
    }
}
impl NativeType for CellCertificate {
    const KIND: &'static str = "rssn/native/discrete/graph-cellcertificate/v1";
    const MAX_BYTES: usize = 65536;

    fn validate(&self) -> Result<(), OpError> {
        limit(self.size(), 64)?;
        Ok(())
    }
}
fn index(
    value: usize,
    n: usize,
) -> Result<usize, OpError> {
    if value >= n {
        Err(error("graph vertex lies outside its domain"))
    } else {
        Ok(value)
    }
}
fn weighted_vector(
    data: Vec<f64>,
    dimension: DimExpr,
) -> Result<Dynamic, OpError> {
    let values = data
        .into_iter()
        .map(|v| extended_quantity(v, dimension.clone()).map(|v| v.0))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Dynamic(Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values,
        elem: TypeExpr::Float.with_dim(dimension),
    })))
}
fn vertex_pairs(data: &[(usize, u32)]) -> Vec<[usize; 2]> {
    data.iter().map(|(a, b)| [*a, *b as usize]).collect()
}
fn cell_counts(data: &[(usize, u32, usize)]) -> Vec<[usize; 3]> {
    data.iter().map(|(a, b, c)| [*a, *b as usize, *c]).collect()
}
fn graph_search(
    g: &ColourGraph,
    budget: SearchBudget,
) -> Result<SearchBudget, OpError> {
    if g.num_nodes()
        .checked_mul(g.num_nodes())
        .and_then(|n| n.checked_mul(budget.max_nodes()))
        .is_none_or(|n| n > 50_000_000)
    {
        return Err(error(
            "graph canonical search exceeds vertex/node work budget",
        ));
    }
    Ok(budget)
}
object_binding!(graph_new,"num_native_graph_new","crate::kernels::graph::Graph::new","Construct a directed weighted graph with an explicit positive reference quantity declaring edge-weight units. Edge values and all shortest-path distances use coherent units.",false,(n:usize,weight_unit:Input)->NewObject<Graph>,|_ctx:&CallCtx<'_>,n:usize,weight_unit:Input|{limit(n,256)?;let(_,d)=api::reference(&weight_unit.0)?;Ok(NewObject::new(Graph::new(n),vec![d]))});
object_binding!(graph_num_nodes,"num_native_graph_num_nodes","crate::kernels::graph::Graph::num_nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>)->usize,|ctx:&CallCtx<'_>,g:Handle<Graph>|{g.with(ctx,|g|Ok(g.num_nodes()))});
object_binding!(graph_add_edge,"num_native_graph_add_edge","crate::kernels::graph::Graph::add_edge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>,u:usize,v:usize,weight:Input)->bool,|ctx:&CallCtx<'_>,g:Handle<Graph>,u:usize,v:usize,weight:Input|{let(w,d)=scalar(&weight)?;api::same_dimension(&d,&object_dimension(&g,ctx)?)?;g.with_mut(ctx,|g|{index(u,g.num_nodes())?;index(v,g.num_nodes())?;if w.abs()>f64::MAX/g.num_nodes().max(1) as f64{return Err(error("edge magnitude exceeds finite path representation range"));}let edges=(0..g.num_nodes()).map(|v|g.adj(v).len()).sum::<usize>();if edges>=65536{return Err(error("weighted graph edge budget exhausted"));}g.add_edge(u,v,w);Ok(true)})});
object_binding!(graph_adj,"num_native_graph_adj","crate::kernels::graph::Graph::adj","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>,u:usize)->Dynamic,|ctx:&CallCtx<'_>,g:Handle<Graph>,u:usize|{let d=object_dimension(&g,ctx)?;g.with(ctx,|g|{index(u,g.num_nodes())?;let(indices,weights):(Vec<_>,Vec<_>)=g.adj(u).iter().copied().unzip();Ok(Dynamic(api::record(vec![("indices",Output::encode(indices)?),("weights",weighted_vector(weights,d)?.0)])))})});
object_binding!(graph_bfs,"num_native_graph_bfs","crate::kernels::graph::bfs","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>,start:usize)->Vec<usize>,|ctx:&CallCtx<'_>,g:Handle<Graph>,start:usize|{g.with(ctx,|g|{index(start,g.num_nodes())?;Ok(crate::kernels::graph::bfs(g,start))})});
object_binding!(graph_connected_components,"num_native_graph_connected_components","crate::kernels::graph::connected_components","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>)->Vec<usize>,|ctx:&CallCtx<'_>,g:Handle<Graph>|{g.with(ctx,|g|Ok(crate::kernels::graph::connected_components(g)))});
object_binding!(graph_minimum_spanning_tree,"num_native_graph_minimum_spanning_tree","crate::kernels::graph::minimum_spanning_tree","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>)->NewObject<Graph>,|ctx:&CallCtx<'_>,g:Handle<Graph>|{let d=object_dimension(&g,ctx)?;g.with(ctx,|g|Ok(NewObject::new(crate::kernels::graph::minimum_spanning_tree(g),vec![d])))});
object_binding!(graph_dijkstra,"num_native_graph_dijkstra","crate::kernels::graph::dijkstra","Native Dijkstra with nonnegative coherent edge weights; infinity denotes unreachable vertices, and predecessor entries are native Options. Finite edge/path representation and allocation bounds are checked.",false,(g:Handle<Graph>,start:usize)->Dynamic,|ctx:&CallCtx<'_>,g:Handle<Graph>,start:usize|{let d=object_dimension(&g,ctx)?;g.with(ctx,|g|{index(start,g.num_nodes())?;if(0..g.num_nodes()).any(|u|g.adj(u).iter().any(|(_,w)|*w<0.0)){return Err(error("Dijkstra requires nonnegative edge weights"));}let(dist,previous)=crate::kernels::graph::dijkstra(g,start);Ok(Dynamic(api::record(vec![("distances",weighted_vector(dist,d)?.0),("previous",NativeOutput::encode(Many(previous),ctx)?)])))})});
object_binding!(graph_floyd_warshall,"num_native_graph_floyd_warshall","crate::kernels::graph::floyd_warshall","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",false,(g:Handle<Graph>)->Dynamic,|ctx:&CallCtx<'_>,g:Handle<Graph>|{let d=object_dimension(&g,ctx)?;g.with(ctx,|g|{let n=g.num_nodes();let distances=crate::kernels::graph::floyd_warshall(g);if(0..n).any(|i|distances[i*n+i]<0.0){return Err(error("all-pairs distances contain a negative cycle"));}let mut result=weighted_vector(distances,d)?.0;if let Value::Tensor(t)=&mut result{t.shape=vec![n,n];}Ok(Dynamic(result))})});
object_binding!(graph_page_rank,"num_native_graph_page_rank","crate::kernels::graph::page_rank","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<Graph>,damping:f64,tolerance:f64,iterations:usize)->Vec<f64>,|ctx:&CallCtx<'_>,g:Handle<Graph>,damping:f64,tolerance:f64,iterations:usize|{if !(0.0..1.0).contains(&damping)||tolerance<=0.0{return Err(error("PageRank damping/tolerance domain invalid"));}g.with(ctx,|g|{let edges=(0..g.num_nodes()).map(|u|g.adj(u).len()).sum::<usize>();if iterations.checked_mul(edges.max(g.num_nodes())).is_none_or(|n|n>50_000_000){return Err(error("PageRank iteration work budget exceeded"));}Ok(crate::kernels::graph::page_rank(g,damping,tolerance,iterations))})});
object_binding!(colour_graph_new,"num_native_colour_graph_new","crate::kernels::graph::ColourGraph::new","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize,directed:bool)->NewObject<ColourGraph>,|_ctx:&CallCtx<'_>,n:usize,directed:bool|{limit(n,64)?;Ok(NewObject::normalized(ColourGraph::new(n,directed)))});
object_binding!(colour_graph_with_colours,"num_native_colour_graph_with_colours","crate::kernels::graph::ColourGraph::with_colours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize,directed:bool,colours:Vec<u32>)->Maybe<NewObject<ColourGraph>>,|_ctx:&CallCtx<'_>,n:usize,directed:bool,colours:Vec<u32>|{limit(n,64)?;Ok(Maybe(ColourGraph::with_colours(n,directed,&colours).map(NewObject::normalized)))});
object_binding!(colour_graph_num_nodes,"num_native_colour_graph_num_nodes","crate::kernels::graph::ColourGraph::num_nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>)->usize,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>|{g.with(ctx,|g|Ok(g.num_nodes()))});
object_binding!(colour_graph_is_directed,"num_native_colour_graph_is_directed","crate::kernels::graph::ColourGraph::is_directed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>)->bool,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>|{g.with(ctx,|g|Ok(g.is_directed()))});
object_binding!(colour_graph_matrix,"num_native_colour_graph_matrix","crate::kernels::graph::ColourGraph::matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>)->Vec<Vec<u32>>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>|{g.with(ctx,|g|Ok(g.matrix().to_vec()))});
object_binding!(colour_graph_colours,"num_native_colour_graph_colours","crate::kernels::graph::ColourGraph::colours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>)->Vec<u32>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>|{g.with(ctx,|g|Ok(g.colours().to_vec()))});
object_binding!(colour_graph_set_colours,"num_native_colour_graph_set_colours","crate::kernels::graph::ColourGraph::set_colours","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,colours:Vec<u32>)->bool,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,colours:Vec<u32>|{g.with_mut(ctx,|g|Ok(g.set_colours(&colours)))});
object_binding!(colour_graph_set_edge,"num_native_colour_graph_set_edge","crate::kernels::graph::ColourGraph::set_edge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize,v:usize,label:u32)->bool,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize,v:usize,label:u32|{g.with_mut(ctx,|g|Ok(g.set_edge(u,v,label)))});
object_binding!(colour_graph_edge_label,"num_native_colour_graph_edge_label","crate::kernels::graph::ColourGraph::edge_label","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize,v:usize)->u32,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize,v:usize|{g.with(ctx,|g|{index(u,g.num_nodes())?;index(v,g.num_nodes())?;Ok(g.edge_label(u,v))})});
object_binding!(colour_graph_has_edge,"num_native_colour_graph_has_edge","crate::kernels::graph::ColourGraph::has_edge","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize,v:usize)->bool,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize,v:usize|{g.with(ctx,|g|{index(u,g.num_nodes())?;index(v,g.num_nodes())?;Ok(g.has_edge(u,v))})});
object_binding!(colour_graph_degree,"num_native_colour_graph_degree","crate::kernels::graph::ColourGraph::degree","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize)->usize,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize|{g.with(ctx,|g|{index(u,g.num_nodes())?;Ok(g.degree(u))})});
object_binding!(colour_graph_in_arcs,"num_native_colour_graph_in_arcs","crate::kernels::graph::ColourGraph::in_arcs","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize)->Vec<[usize;2]>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize|{g.with(ctx,|g|{index(u,g.num_nodes())?;Ok(vertex_pairs(g.in_arcs(u)))})});
object_binding!(colour_graph_out_arcs,"num_native_colour_graph_out_arcs","crate::kernels::graph::ColourGraph::out_arcs","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,u:usize)->Vec<[usize;2]>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,u:usize|{g.with(ctx,|g|{index(u,g.num_nodes())?;Ok(vertex_pairs(g.out_arcs(u)))})});
object_binding!(graph_permutation_identity,"num_native_graph_permutation_identity","crate::kernels::graph::Permutation::identity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize)->NewObject<Permutation>,|_ctx:&CallCtx<'_>,n:usize|{limit(n,256)?;Ok(NewObject::normalized(Permutation::identity(n)))});
object_binding!(graph_permutation_from_images,"num_native_graph_permutation_from_images","crate::kernels::graph::Permutation::from_images","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(image:Vec<usize>)->Maybe<NewObject<Permutation>>,|_ctx:&CallCtx<'_>,image:Vec<usize>|{limit(image.len(),256)?;Ok(Maybe(Permutation::from_images(image).map(NewObject::normalized)))});
object_binding!(graph_permutation_len,"num_native_graph_permutation_len","crate::kernels::graph::Permutation::len","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->usize,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(o.len()))});
object_binding!(graph_permutation_is_empty,"num_native_graph_permutation_is_empty","crate::kernels::graph::Permutation::is_empty","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->bool,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(o.is_empty()))});
object_binding!(graph_permutation_is_identity,"num_native_graph_permutation_is_identity","crate::kernels::graph::Permutation::is_identity","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->bool,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(o.is_identity()))});
object_binding!(graph_permutation_image,"num_native_graph_permutation_image","crate::kernels::graph::Permutation::image","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(o.image().to_vec()))});
object_binding!(graph_permutation_apply,"num_native_graph_permutation_apply","crate::kernels::graph::Permutation::apply","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>,v:usize)->Option<usize>,|ctx:&CallCtx<'_>,a:Handle<Permutation>,v:usize|{a.with(ctx,|o|Ok(o.apply(v)))});
object_binding!(graph_permutation_inverse,"num_native_graph_permutation_inverse","crate::kernels::graph::Permutation::inverse","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->NewObject<Permutation>,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(NewObject::normalized(o.inverse())))});
object_binding!(graph_permutation_compose,"num_native_graph_permutation_compose","crate::kernels::graph::Permutation::compose","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>,b:Handle<Permutation>)->Maybe<NewObject<Permutation>>,|ctx:&CallCtx<'_>,a:Handle<Permutation>,b:Handle<Permutation>|{a.with(ctx,|o|b.with(ctx,|p|Ok(Maybe(o.compose(p).map(NewObject::normalized)))))});
object_binding!(graph_permutation_cycles,"num_native_graph_permutation_cycles","crate::kernels::graph::Permutation::cycles","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Permutation>)->Vec<Vec<usize>>,|ctx:&CallCtx<'_>,a:Handle<Permutation>|{a.with(ctx,|o|Ok(o.cycles()))});
object_binding!(graph_search_budget_nodes,"num_native_graph_search_budget_nodes","crate::kernels::graph::SearchBudget::nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(nodes:usize)->NewObject<SearchBudget>,|_ctx:&CallCtx<'_>,nodes:usize|{limit(nodes,1_000_000)?;Ok(NewObject::normalized(SearchBudget::nodes(nodes)))});
object_binding!(graph_search_budget_for_size,"num_native_graph_search_budget_for_size","crate::kernels::graph::SearchBudget::for_size","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(n:usize)->NewObject<SearchBudget>,|_ctx:&CallCtx<'_>,n:usize|{limit(n,64)?;Ok(NewObject::normalized(SearchBudget::for_size(n)))});
object_binding!(graph_search_budget_max_nodes,"num_native_graph_search_budget_max_nodes","crate::kernels::graph::SearchBudget::max_nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<SearchBudget>)->usize,|ctx:&CallCtx<'_>,a:Handle<SearchBudget>|{a.with(ctx,|o|Ok(o.max_nodes()))});
object_binding!(graph_refine,"num_native_graph_refine","crate::kernels::graph::refine","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>)->NewObject<Refinement>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>|{g.with(ctx,|g|Ok(NewObject::normalized(crate::kernels::graph::refine(g))))});
object_binding!(graph_canonical_labelling,"num_native_graph_canonical_labelling","crate::kernels::graph::canonical_labelling","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(g:Handle<ColourGraph>,budget:Handle<SearchBudget>)->NewObject<CanonicalLabelling>,|ctx:&CallCtx<'_>,g:Handle<ColourGraph>,budget:Handle<SearchBudget>|{g.with(ctx,|g|budget.with(ctx,|b|Ok(NewObject::normalized(crate::kernels::graph::canonical_labelling(g,graph_search(g,*b)?)))))});
object_binding!(graph_isomorphic,"num_native_graph_isomorphic","crate::kernels::graph::isomorphic","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<ColourGraph>,b:Handle<ColourGraph>,budget:Handle<SearchBudget>)->Option<bool>,|ctx:&CallCtx<'_>,a:Handle<ColourGraph>,b:Handle<ColourGraph>,budget:Handle<SearchBudget>|{a.with(ctx,|a|b.with(ctx,|b|budget.with(ctx,|budget|{graph_search(a,*budget)?;graph_search(b,*budget)?;Ok(crate::kernels::graph::isomorphic(a,b,*budget))})))});
object_binding!(graph_wl_equivalent,"num_native_graph_wl_equivalent","crate::kernels::graph::wl_equivalent","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<ColourGraph>,b:Handle<ColourGraph>)->bool,|ctx:&CallCtx<'_>,a:Handle<ColourGraph>,b:Handle<ColourGraph>|{a.with(ctx,|a|b.with(ctx,|b|Ok(crate::kernels::graph::wl_equivalent(a,b))))});
object_binding!(graph_refinement_cell_of,"num_native_graph_refinement_cell_of","crate::kernels::graph::Refinement::cell_of","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.cell_of().to_vec()))});
object_binding!(graph_refinement_cell_sizes,"num_native_graph_refinement_cell_sizes","crate::kernels::graph::Refinement::cell_sizes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.cell_sizes()))});
object_binding!(graph_refinement_cells,"num_native_graph_refinement_cells","crate::kernels::graph::Refinement::cells","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->Vec<Vec<usize>>,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.cells().to_vec()))});
object_binding!(graph_refinement_is_directed,"num_native_graph_refinement_is_directed","crate::kernels::graph::Refinement::is_directed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->bool,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.is_directed()))});
object_binding!(graph_refinement_is_discrete,"num_native_graph_refinement_is_discrete","crate::kernels::graph::Refinement::is_discrete","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->bool,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.is_discrete()))});
object_binding!(graph_refinement_num_cells,"num_native_graph_refinement_num_cells","crate::kernels::graph::Refinement::num_cells","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->usize,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.num_cells()))});
object_binding!(graph_refinement_rounds,"num_native_graph_refinement_rounds","crate::kernels::graph::Refinement::rounds","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->usize,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(o.rounds()))});
object_binding!(graph_certificate_cell_sizes,"num_native_graph_certificate_cell_sizes","crate::kernels::graph::RefinementCertificate::cell_sizes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(o.cell_sizes()))});
object_binding!(graph_certificate_is_directed,"num_native_graph_certificate_is_directed","crate::kernels::graph::RefinementCertificate::is_directed","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->bool,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(o.is_directed()))});
object_binding!(graph_certificate_is_discrete,"num_native_graph_certificate_is_discrete","crate::kernels::graph::RefinementCertificate::is_discrete","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->bool,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(o.is_discrete()))});
object_binding!(graph_certificate_num_cells,"num_native_graph_certificate_num_cells","crate::kernels::graph::RefinementCertificate::num_cells","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->usize,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(o.num_cells()))});
object_binding!(graph_certificate_num_nodes,"num_native_graph_certificate_num_nodes","crate::kernels::graph::RefinementCertificate::num_nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->usize,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(o.num_nodes()))});
object_binding!(graph_cell_certificate_colour,"num_native_graph_cell_certificate_colour","crate::kernels::graph::CellCertificate::colour","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CellCertificate>)->u32,|ctx:&CallCtx<'_>,a:Handle<CellCertificate>|{a.with(ctx,|o|Ok(o.colour()))});
object_binding!(graph_cell_certificate_size,"num_native_graph_cell_certificate_size","crate::kernels::graph::CellCertificate::size","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CellCertificate>)->usize,|ctx:&CallCtx<'_>,a:Handle<CellCertificate>|{a.with(ctx,|o|Ok(o.size()))});
object_binding!(graph_canonical_code,"num_native_graph_canonical_code","crate::kernels::graph::CanonicalLabelling::code","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->Vec<u32>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.code().to_vec()))});
object_binding!(graph_canonical_group_order,"num_native_graph_canonical_group_order","crate::kernels::graph::CanonicalLabelling::group_order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->num_bigint::BigInt,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.group_order().clone()))});
object_binding!(graph_canonical_is_exact,"num_native_graph_canonical_is_exact","crate::kernels::graph::CanonicalLabelling::is_exact","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->bool,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.is_exact()))});
object_binding!(graph_canonical_matrix,"num_native_graph_canonical_matrix","crate::kernels::graph::CanonicalLabelling::matrix","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->Vec<Vec<u32>>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.matrix()))});
object_binding!(graph_canonical_num_nodes,"num_native_graph_canonical_num_nodes","crate::kernels::graph::CanonicalLabelling::num_nodes","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->usize,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.num_nodes()))});
object_binding!(graph_canonical_num_nodes_visited,"num_native_graph_canonical_num_nodes_visited","crate::kernels::graph::CanonicalLabelling::num_nodes_visited","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->usize,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.num_nodes_visited()))});
object_binding!(graph_canonical_orbits,"num_native_graph_canonical_orbits","crate::kernels::graph::CanonicalLabelling::orbits","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->Vec<Vec<usize>>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.orbits().to_vec()))});
object_binding!(graph_canonical_order,"num_native_graph_canonical_order","crate::kernels::graph::CanonicalLabelling::order","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->Vec<usize>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(o.order().to_vec()))});
object_binding!(graph_refinement_certificate,"num_native_graph_refinement_certificate","crate::kernels::graph::Refinement::certificate","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>)->NewObject<RefinementCertificate>,|ctx:&CallCtx<'_>,a:Handle<Refinement>|{a.with(ctx,|o|Ok(NewObject::normalized(o.certificate().clone())))});
object_binding!(graph_refinement_verify,"num_native_graph_refinement_verify","crate::kernels::graph::Refinement::verify","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<Refinement>,g:Handle<ColourGraph>)->bool,|ctx:&CallCtx<'_>,a:Handle<Refinement>,g:Handle<ColourGraph>|{a.with(ctx,|o|g.with(ctx,|g|Ok(o.verify(g))))});
object_binding!(graph_certificate_cells,"num_native_graph_certificate_cells","crate::kernels::graph::RefinementCertificate::cells","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<RefinementCertificate>)->Many<NewObject<CellCertificate>>,|ctx:&CallCtx<'_>,a:Handle<RefinementCertificate>|{a.with(ctx,|o|Ok(Many(o.cells().iter().cloned().map(NewObject::normalized).collect())))});
object_binding!(graph_cell_certificate_in_counts,"num_native_graph_cell_certificate_in_counts","crate::kernels::graph::CellCertificate::in_counts","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CellCertificate>)->Vec<[usize;3]>,|ctx:&CallCtx<'_>,a:Handle<CellCertificate>|{a.with(ctx,|o|Ok(cell_counts(o.in_counts())))});
object_binding!(graph_cell_certificate_out_counts,"num_native_graph_cell_certificate_out_counts","crate::kernels::graph::CellCertificate::out_counts","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CellCertificate>)->Vec<[usize;3]>,|ctx:&CallCtx<'_>,a:Handle<CellCertificate>|{a.with(ctx,|o|Ok(cell_counts(o.out_counts())))});
object_binding!(graph_canonical_generators,"num_native_graph_canonical_generators","crate::kernels::graph::CanonicalLabelling::generators","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>)->Many<NewObject<Permutation>>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>|{a.with(ctx,|o|Ok(Many(o.generators().iter().cloned().map(NewObject::normalized).collect())))});
object_binding!(graph_canonical_position,"num_native_graph_canonical_position","crate::kernels::graph::CanonicalLabelling::position","Checked native method over explicitly normalized values. Session/kind, shape and work/allocation bounds are enforced.",true,(a:Handle<CanonicalLabelling>,vertex:usize)->Option<usize>,|ctx:&CallCtx<'_>,a:Handle<CanonicalLabelling>,vertex:usize|{a.with(ctx,|o|Ok(o.position(vertex)))});
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install(
        i,
        &[
            &graph_new::BINDING,
            &graph_num_nodes::BINDING,
            &graph_add_edge::BINDING,
            &graph_adj::BINDING,
            &graph_bfs::BINDING,
            &graph_connected_components::BINDING,
            &graph_minimum_spanning_tree::BINDING,
            &graph_dijkstra::BINDING,
            &graph_floyd_warshall::BINDING,
            &graph_page_rank::BINDING,
            &colour_graph_new::BINDING,
            &colour_graph_with_colours::BINDING,
            &colour_graph_num_nodes::BINDING,
            &colour_graph_is_directed::BINDING,
            &colour_graph_matrix::BINDING,
            &colour_graph_colours::BINDING,
            &colour_graph_set_colours::BINDING,
            &colour_graph_set_edge::BINDING,
            &colour_graph_edge_label::BINDING,
            &colour_graph_has_edge::BINDING,
            &colour_graph_degree::BINDING,
            &colour_graph_in_arcs::BINDING,
            &colour_graph_out_arcs::BINDING,
            &graph_permutation_identity::BINDING,
            &graph_permutation_from_images::BINDING,
            &graph_permutation_len::BINDING,
            &graph_permutation_is_empty::BINDING,
            &graph_permutation_is_identity::BINDING,
            &graph_permutation_image::BINDING,
            &graph_permutation_apply::BINDING,
            &graph_permutation_inverse::BINDING,
            &graph_permutation_compose::BINDING,
            &graph_permutation_cycles::BINDING,
            &graph_search_budget_nodes::BINDING,
            &graph_search_budget_for_size::BINDING,
            &graph_search_budget_max_nodes::BINDING,
            &graph_refine::BINDING,
            &graph_canonical_labelling::BINDING,
            &graph_isomorphic::BINDING,
            &graph_wl_equivalent::BINDING,
            &graph_refinement_cell_of::BINDING,
            &graph_refinement_cell_sizes::BINDING,
            &graph_refinement_cells::BINDING,
            &graph_refinement_is_directed::BINDING,
            &graph_refinement_is_discrete::BINDING,
            &graph_refinement_num_cells::BINDING,
            &graph_refinement_rounds::BINDING,
            &graph_certificate_cell_sizes::BINDING,
            &graph_certificate_is_directed::BINDING,
            &graph_certificate_is_discrete::BINDING,
            &graph_certificate_num_cells::BINDING,
            &graph_certificate_num_nodes::BINDING,
            &graph_cell_certificate_colour::BINDING,
            &graph_cell_certificate_size::BINDING,
            &graph_canonical_code::BINDING,
            &graph_canonical_group_order::BINDING,
            &graph_canonical_is_exact::BINDING,
            &graph_canonical_matrix::BINDING,
            &graph_canonical_num_nodes::BINDING,
            &graph_canonical_num_nodes_visited::BINDING,
            &graph_canonical_orbits::BINDING,
            &graph_canonical_order::BINDING,
            &graph_refinement_certificate::BINDING,
            &graph_refinement_verify::BINDING,
            &graph_certificate_cells::BINDING,
            &graph_cell_certificate_in_counts::BINDING,
            &graph_cell_certificate_out_counts::BINDING,
            &graph_canonical_generators::BINDING,
            &graph_canonical_position::BINDING,
        ],
        "numerical_vector",
    )
}
