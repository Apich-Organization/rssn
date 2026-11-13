//! A planner that fuses chains of kernel launches.
//!
//! Launching a kernel costs a dispatch and a round trip through memory for
//! every column it reads and writes. A pipeline such as
//! `y = f(x); z = g(y, w); s = sum(h(z))` written as three kernels moves `y`
//! and `z` through device memory twice for nothing. A [`Plan`] describes the
//! pipeline as a graph of kernels whose inputs are external columns or the
//! outputs of earlier kernels; [`Plan::fuse`] merges every producer into the
//! consumer that reads it, so the intermediates become registers, and
//! [`Plan::run`] executes what is left with one launch per remaining kernel.
//!
//! # What fuses
//!
//! Producer `P` fuses into consumer `C` when
//!
//! * `C` reads the columns that come from `P` only with `Input` (an addressed
//!   `Load`, a `Len` or a `Gather` needs the whole column in memory, so the
//!   edge stays a launch boundary);
//! * at most one of the two has a [`Fold`](super::kernel::Fold) (a kernel has
//!   one loop);
//! * the merged program validates ([`KernelIr::analyse`]).
//!
//! If other kernels or the results still need an output of `P`, the merged
//! kernel keeps it as an additional output, so nothing is recomputed.
//!
//! Every column of a plan has `rows` elements (the launch size); a program
//! whose columns have other lengths is a single [`Device::launch_cols`] and
//! needs no planner.
//!
//! Fusion changes only where intermediate values live: each value is
//! computed by the same instruction as before, so a fused plan gives the
//! same bits as the unfused one on a deterministic device.

use super::Buffer;
use super::Device;
use super::DeviceError;
use super::KernelIr;
use super::KernelOp;
use super::kernel::Fold;
use super::preferred_dtype;
use crate::backend::tape::Reg;

/// A kernel of a [`Plan`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(pub usize);

/// Where a kernel input column comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Source {
    /// Column `i` of the data handed to [`Plan::run`].
    External(usize),
    /// Output `k` of an earlier kernel.
    Output(NodeId, usize),
}

#[derive(Clone, Debug)]
struct Node {
    ir: KernelIr,
    inputs: Vec<Source>,
    params: Vec<f64>,
}

/// A graph of kernels over columns of equal length.
#[derive(Clone, Debug, Default)]
pub struct Plan {
    nodes: Vec<Option<Node>>,
    externals: usize,
    results: Vec<Source>,
}

fn invalid(why: impl Into<String>) -> DeviceError {
    DeviceError::InvalidKernel(why.into())
}

/// `op` with its registers, columns and parameters renamed.
fn rename(
    op: &KernelOp,
    reg: &dyn Fn(Reg) -> Reg,
    col: &dyn Fn(u32) -> u32,
    param: u32,
) -> KernelOp {
    use KernelOp as K;
    match *op {
        | K::Const(c) => K::Const(c),
        | K::Input(i) => K::Input(col(i)),
        | K::Param(i) => K::Param(i + param),
        | K::Add(a, b) => K::Add(reg(a), reg(b)),
        | K::Sub(a, b) => K::Sub(reg(a), reg(b)),
        | K::Mul(a, b) => K::Mul(reg(a), reg(b)),
        | K::Div(a, b) => K::Div(reg(a), reg(b)),
        | K::Neg(a) => K::Neg(reg(a)),
        | K::PowI(a, n) => K::PowI(reg(a), n),
        | K::PowF(a, b) => K::PowF(reg(a), reg(b)),
        | K::Unary(f, a) => K::Unary(f, reg(a)),
        | K::Atan2(a, b) => K::Atan2(reg(a), reg(b)),
        | K::Cmp(c, a, b) => K::Cmp(c, reg(a), reg(b)),
        | K::Select(c, a, b) => K::Select(reg(c), reg(a), reg(b)),
        | K::Gather(base, a) => K::Gather(base, reg(a)),
        | K::Row => K::Row,
        | K::Rows => K::Rows,
        | K::Len(i) => K::Len(col(i)),
        | K::Inner => K::Inner,
        | K::IConst(c) => K::IConst(c),
        | K::IAdd(a, b) => K::IAdd(reg(a), reg(b)),
        | K::ISub(a, b) => K::ISub(reg(a), reg(b)),
        | K::IMul(a, b) => K::IMul(reg(a), reg(b)),
        | K::IDiv(a, b) => K::IDiv(reg(a), reg(b)),
        | K::IMod(a, b) => K::IMod(reg(a), reg(b)),
        | K::IMin(a, b) => K::IMin(reg(a), reg(b)),
        | K::IMax(a, b) => K::IMax(reg(a), reg(b)),
        | K::ICmp(c, a, b) => K::ICmp(c, reg(a), reg(b)),
        | K::ToNum(a) => K::ToNum(reg(a)),
        | K::ToIdx(a) => K::ToIdx(reg(a)),
        | K::Load(i, a) => K::Load(col(i), reg(a)),
        | K::Acc(k) => K::Acc(k),
    }
}

/// The merged kernel of producer `p` into consumer `c`. `from_p[k]` is
/// `Some(j)` when input column `k` of `c` is output `j` of `p`, and
/// `column[k]` the column of the merged kernel for every other `k`. Returns
/// the kernel and the register of each output of `p` in it.
fn merge(
    p: &KernelIr,
    c: &KernelIr,
    from_p: &[Option<usize>],
    column: &[u32],
) -> Option<(KernelIr, Vec<Reg>)> {
    if p.fold.is_some() && c.fold.is_some() {
        return None;
    }
    for op in &c.ops {
        match *op {
            // A gather addresses a range of columns by position.
            | KernelOp::Gather(..) => return None,
            | KernelOp::Len(i) | KernelOp::Load(i, _) if from_p[i as usize].is_some() => {
                return None;
            },
            | _ => {},
        }
    }
    let mut ops = p.ops.clone();
    let mut map: Vec<Reg> = Vec::with_capacity(c.ops.len());
    for op in &c.ops {
        if let KernelOp::Input(i) = *op
            && let Some(j) = from_p[i as usize]
        {
            map.push(*p.outputs.get(j)?);
            continue;
        }
        let renamed = rename(
            op,
            &|r| map.get(r as usize).copied().unwrap_or(r),
            &|i| column[i as usize],
            u32::try_from(p.params).ok()?,
        );
        ops.push(renamed);
        map.push(Reg::try_from(ops.len() - 1).ok()?);
    }
    let fold = match (&p.fold, &c.fold) {
        | (Some(f), None) => Some(f.clone()),
        | (None, Some(f)) => {
            Some(Fold {
                count: map[f.count as usize],
                accs: f
                    .accs
                    .iter()
                    .map(|&(op, r)| (op, map[r as usize]))
                    .collect(),
            })
        },
        | _ => None,
    };
    let columns = column
        .iter()
        .zip(from_p)
        .filter(|(_, from)| from.is_none())
        .map(|(&k, _)| k as usize + 1)
        .max()
        .unwrap_or(0)
        .max(p.inputs);
    let merged = KernelIr {
        ops,
        inputs: columns,
        params: p.params + c.params,
        outputs: c.outputs.iter().map(|&r| map[r as usize]).collect(),
        fold,
    };
    Some((merged, p.outputs.clone()))
}

impl Plan {
    /// A plan over `externals` input columns.
    #[must_use]
    pub const fn new(externals: usize) -> Self {
        Self {
            nodes: Vec::new(),
            externals,
            results: Vec::new(),
        }
    }

    /// Adds a kernel reading `inputs` (one source per input column of `ir`)
    /// with the uniform `params`.
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] when the kernel does not validate, the
    /// counts do not match, or a source names a later kernel, a missing
    /// output or a missing external column.
    pub fn add(
        &mut self,
        ir: KernelIr,
        inputs: Vec<Source>,
        params: Vec<f64>,
    ) -> Result<NodeId, DeviceError> {
        ir.analyse()?;
        if inputs.len() != ir.inputs || params.len() != ir.params {
            return Err(invalid(format!(
                "kernel takes {} columns and {} parameters, {} and {} given",
                ir.inputs,
                ir.params,
                inputs.len(),
                params.len()
            )));
        }
        for &s in &inputs {
            self.check(s)?;
        }
        self.nodes.push(Some(Node { ir, inputs, params }));
        Ok(NodeId(self.nodes.len() - 1))
    }

    fn check(
        &self,
        s: Source,
    ) -> Result<(), DeviceError> {
        match s {
            | Source::External(i) if i < self.externals => Ok(()),
            | Source::External(i) => Err(invalid(format!("external column {i} does not exist"))),
            | Source::Output(NodeId(n), k) => {
                match self.nodes.get(n) {
                    | Some(Some(node)) if k < node.ir.outputs.len() => Ok(()),
                    | _ => Err(invalid(format!("kernel {n} has no output {k}"))),
                }
            },
        }
    }

    /// Asks for `source` to be returned by [`Plan::run`] (results come back
    /// in the order they were asked for).
    ///
    /// # Errors
    /// [`DeviceError::InvalidKernel`] for a source that does not exist.
    pub fn result(
        &mut self,
        source: Source,
    ) -> Result<(), DeviceError> {
        self.check(source)?;
        self.results.push(source);
        Ok(())
    }

    /// The number of launches [`Plan::run`] performs.
    #[must_use]
    pub fn launches(&self) -> usize {
        self.nodes.iter().flatten().count()
    }

    /// Merges every producer into a consumer that reads it, as far as the
    /// rules in the module documentation allow. Returns the number of
    /// launches removed.
    pub fn fuse(&mut self) -> usize {
        let before = self.launches();
        'again: loop {
            for c in 0..self.nodes.len() {
                let Some(consumer) = &self.nodes[c] else {
                    continue;
                };
                let producers: Vec<usize> = consumer
                    .inputs
                    .iter()
                    .filter_map(|s| {
                        match *s {
                            | Source::Output(NodeId(p), _) => Some(p),
                            | Source::External(_) => None,
                        }
                    })
                    .collect();
                for p in producers {
                    if self.try_fuse(p, c) {
                        continue 'again;
                    }
                }
            }
            break;
        }
        before - self.launches()
    }

    /// Fuses node `p` into node `c` when possible.
    fn try_fuse(
        &mut self,
        p: usize,
        c: usize,
    ) -> bool {
        let (Some(producer), Some(consumer)) = (&self.nodes[p], &self.nodes[c]) else {
            return false;
        };
        // The merged kernel's columns: the producer's, then the consumer's
        // other sources, without repeats.
        let mut sources = producer.inputs.clone();
        let mut from_p = Vec::with_capacity(consumer.inputs.len());
        let mut column = Vec::with_capacity(consumer.inputs.len());
        for &s in &consumer.inputs {
            if let Source::Output(NodeId(n), j) = s
                && n == p
            {
                from_p.push(Some(j));
                column.push(0);
                continue;
            }
            from_p.push(None);
            let at = sources.iter().position(|&t| t == s).unwrap_or_else(|| {
                sources.push(s);
                sources.len() - 1
            });
            column.push(u32::try_from(at).unwrap_or(u32::MAX));
        }
        let Some((mut merged, p_regs)) = merge(&producer.ir, &consumer.ir, &from_p, &column) else {
            return false;
        };
        merged.inputs = sources.len();
        // Outputs of `p` that something other than `c` still reads are kept
        // as extra outputs of the merged kernel.
        let mut redirect: Vec<(usize, usize)> = Vec::new();
        for (j, &reg) in p_regs.iter().enumerate() {
            let needed_elsewhere = self.results.contains(&Source::Output(NodeId(p), j))
                || self.nodes.iter().enumerate().any(|(n, node)| {
                    n != c
                        && node
                            .as_ref()
                            .is_some_and(|node| node.inputs.contains(&Source::Output(NodeId(p), j)))
                });
            if needed_elsewhere {
                merged.outputs.push(reg);
                redirect.push((j, merged.outputs.len() - 1));
            }
        }
        if merged.analyse().is_err() {
            return false;
        }
        let mut params = producer.params.clone();
        params.extend_from_slice(&consumer.params);
        let node = Node {
            ir: merged,
            inputs: sources,
            params,
        };
        self.nodes[p] = None;
        self.nodes[c] = Some(node);
        let moved = |s: &mut Source| {
            if let Source::Output(NodeId(n), j) = *s
                && n == p
                && let Some(&(_, k)) = redirect.iter().find(|&&(from, _)| from == j)
            {
                *s = Source::Output(NodeId(c), k);
            }
        };
        for node in self.nodes.iter_mut().flatten() {
            node.inputs.iter_mut().for_each(moved);
        }
        self.results.iter_mut().for_each(moved);
        true
    }

    /// Runs the plan on `device`: `externals[i]` is external column `i`
    /// (`rows` values each). Returns the requested results.
    ///
    /// # Errors
    /// [`DeviceError::Mismatch`] when a column has the wrong length or the
    /// count of externals is wrong; any error of the device.
    pub fn run(
        &self,
        device: &dyn Device,
        rows: usize,
        externals: &[&[f64]],
    ) -> Result<Vec<Vec<f64>>, DeviceError> {
        if externals.len() != self.externals || externals.iter().any(|c| c.len() != rows) {
            return Err(DeviceError::Mismatch(format!(
                "the plan reads {} external columns of {rows} rows",
                self.externals
            )));
        }
        let dtype = preferred_dtype(device.info());
        let mut flat = Vec::with_capacity(self.externals * rows);
        for c in externals {
            flat.extend_from_slice(c);
        }
        let external = device.upload(&flat, dtype)?;
        let mut outputs: Vec<Option<Buffer>> = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            let Some(node) = node else {
                outputs.push(None);
                continue;
            };
            let mut input = device.alloc(node.inputs.len() * rows, dtype)?;
            for (k, &s) in node.inputs.iter().enumerate() {
                match s {
                    | Source::External(i) => {
                        device.copy(&external, i * rows, &mut input, k * rows, rows)?;
                    },
                    | Source::Output(NodeId(n), j) => {
                        let produced = outputs
                            .get(n)
                            .and_then(Option::as_ref)
                            .ok_or_else(|| invalid("a kernel reads a fused-away output"))?;
                        device.copy(produced, j * rows, &mut input, k * rows, rows)?;
                    },
                }
            }
            let kernel = device.compile(&node.ir)?;
            let mut out = device.alloc(node.ir.outputs.len() * rows, dtype)?;
            device.launch(&*kernel, rows, &input, &node.params, &mut out)?;
            outputs.push(Some(out));
        }
        let mut answers = Vec::with_capacity(self.results.len());
        for &s in &self.results {
            let mut host = vec![0.0; rows];
            match s {
                | Source::External(i) => host.copy_from_slice(externals[i]),
                | Source::Output(NodeId(n), j) => {
                    let produced = outputs
                        .get(n)
                        .and_then(Option::as_ref)
                        .ok_or_else(|| invalid("a result names a fused-away output"))?;
                    let mut column = device.alloc(rows, dtype)?;
                    device.copy(produced, j * rows, &mut column, 0, rows)?;
                    device.download(&column, &mut host)?;
                },
            }
            answers.push(host);
        }
        Ok(answers)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::backend::Intrinsic;
    use crate::backend::device::CpuDevice;
    use crate::backend::tape::Inst;
    use crate::backend::tape::Tape;

    fn kernel(
        insts: Vec<Inst>,
        inputs: usize,
        params: usize,
        outputs: Vec<u32>,
    ) -> KernelIr {
        KernelIr::from_tape(&Tape {
            insts,
            args: Vec::new(),
            inputs,
            params,
            outputs,
        })
        .unwrap()
    }

    /// `a * x + 1` with `a` a parameter.
    fn affine() -> KernelIr {
        kernel(
            vec![
                Inst::Input(0),
                Inst::Param(0),
                Inst::Mul(1, 0),
                Inst::Const(1.0),
                Inst::Add(2, 3),
            ],
            1,
            1,
            vec![4],
        )
    }

    /// `sin(y) * w`.
    fn modulate() -> KernelIr {
        kernel(
            vec![
                Inst::Input(0),
                Inst::Input(1),
                Inst::Unary(Intrinsic::Sin, 0),
                Inst::Mul(2, 1),
            ],
            2,
            0,
            vec![3],
        )
    }

    fn data() -> (Vec<f64>, Vec<f64>) {
        let x: Vec<f64> = (0..100).map(|i| f64::from(i) * 0.05 - 2.0).collect();
        let w: Vec<f64> = (0..100).map(|i| 1.0 + f64::from(i % 7)).collect();
        (x, w)
    }

    #[test]
    fn a_chain_fuses_into_one_launch_with_the_same_bits() {
        let (x, w) = data();
        let build = || {
            let mut plan = Plan::new(2);
            let y = plan
                .add(affine(), vec![Source::External(0)], vec![3.0])
                .unwrap();
            let z = plan
                .add(
                    modulate(),
                    vec![Source::Output(y, 0), Source::External(1)],
                    vec![],
                )
                .unwrap();
            plan.result(Source::Output(z, 0)).unwrap();
            plan
        };
        let device = CpuDevice::sequential();
        let plain = build();
        assert_eq!(plain.launches(), 2);
        let want = plain.run(&device, 100, &[&x, &w]).unwrap();
        let mut fused = build();
        assert_eq!(fused.fuse(), 1);
        assert_eq!(fused.launches(), 1);
        let got = fused.run(&device, 100, &[&x, &w]).unwrap();
        assert_eq!(got, want);
        for (i, v) in got[0].iter().enumerate() {
            assert_eq!(*v, (3.0 * x[i] + 1.0).sin() * w[i]);
        }
    }

    #[test]
    fn an_output_needed_elsewhere_is_kept_and_not_recomputed() {
        let (x, w) = data();
        let mut plan = Plan::new(2);
        let y = plan
            .add(affine(), vec![Source::External(0)], vec![2.0])
            .unwrap();
        let z = plan
            .add(
                modulate(),
                vec![Source::Output(y, 0), Source::External(1)],
                vec![],
            )
            .unwrap();
        // Both the intermediate and the final value are results.
        plan.result(Source::Output(y, 0)).unwrap();
        plan.result(Source::Output(z, 0)).unwrap();
        let device = CpuDevice::sequential();
        let want = plan.run(&device, 100, &[&x, &w]).unwrap();
        assert_eq!(plan.fuse(), 1);
        assert_eq!(plan.launches(), 1);
        assert_eq!(plan.run(&device, 100, &[&x, &w]).unwrap(), want);
        assert_eq!(want[0][10], 2.0 * x[10] + 1.0);
    }

    #[test]
    fn parameters_of_the_consumer_are_renumbered_after_the_producers() {
        let (x, _) = data();
        let scale = kernel(
            vec![Inst::Input(0), Inst::Param(0), Inst::Mul(0, 1)],
            1,
            1,
            vec![2],
        );
        let mut plan = Plan::new(1);
        let a = plan
            .add(affine(), vec![Source::External(0)], vec![4.0])
            .unwrap();
        let b = plan
            .add(scale, vec![Source::Output(a, 0)], vec![0.5])
            .unwrap();
        plan.result(Source::Output(b, 0)).unwrap();
        let device = CpuDevice::sequential();
        let want = plan.run(&device, 100, &[&x]).unwrap();
        assert_eq!(plan.fuse(), 1);
        let got = plan.run(&device, 100, &[&x]).unwrap();
        assert_eq!(got, want);
        for (i, v) in got[0].iter().enumerate() {
            assert_eq!(*v, f64::midpoint(4.0 * x[i], 1.0));
        }
    }

    #[test]
    fn addressed_reads_and_two_folds_stay_launch_boundaries() {
        let (x, _) = data();
        // The consumer loads an arbitrary element of the produced column.
        let loader = KernelIr {
            ops: vec![KernelOp::IConst(0), KernelOp::Load(0, 0)],
            inputs: 1,
            params: 0,
            outputs: vec![1],
            fold: None,
        };
        let mut plan = Plan::new(1);
        let y = plan
            .add(affine(), vec![Source::External(0)], vec![1.0])
            .unwrap();
        let z = plan
            .add(loader, vec![Source::Output(y, 0)], vec![])
            .unwrap();
        // Nor does a gather.
        let gather = kernel(vec![Inst::Const(0.0), Inst::Gather(0, 0)], 1, 0, vec![1]);
        let g = plan
            .add(gather, vec![Source::Output(y, 0)], vec![])
            .unwrap();
        plan.result(Source::Output(g, 0)).unwrap();
        plan.result(Source::Output(z, 0)).unwrap();
        let device = CpuDevice::sequential();
        let want = plan.run(&device, 100, &[&x]).unwrap();
        assert_eq!(plan.fuse(), 0);
        assert_eq!(plan.launches(), 3);
        assert_eq!(plan.run(&device, 100, &[&x]).unwrap(), want);
        // The gather reads its own row; the load reads element 0.
        assert_eq!(want[0][7], x[7] + 1.0);
        assert!(want[1].iter().all(|&v| v.total_cmp(&(x[0] + 1.0)).is_eq()));
    }

    #[test]
    fn independent_kernels_and_shared_externals_share_columns() {
        let (x, w) = data();
        let mut plan = Plan::new(2);
        let a = plan
            .add(affine(), vec![Source::External(0)], vec![1.0])
            .unwrap();
        // The consumer also reads external 0 directly: the merged kernel
        // reads it once.
        let both = kernel(
            vec![
                Inst::Input(0),
                Inst::Input(1),
                Inst::Input(2),
                Inst::Add(0, 1),
                Inst::Add(3, 2),
            ],
            3,
            0,
            vec![4],
        );
        let c = plan
            .add(
                both,
                vec![
                    Source::Output(a, 0),
                    Source::External(0),
                    Source::External(1),
                ],
                vec![],
            )
            .unwrap();
        plan.result(Source::Output(c, 0)).unwrap();
        let device = CpuDevice::sequential();
        let want = plan.run(&device, 100, &[&x, &w]).unwrap();
        assert_eq!(plan.fuse(), 1);
        assert_eq!(plan.nodes.iter().flatten().next().unwrap().inputs.len(), 2);
        assert_eq!(plan.run(&device, 100, &[&x, &w]).unwrap(), want);
    }

    #[test]
    fn bad_plans_are_rejected() {
        let mut plan = Plan::new(1);
        assert!(plan.add(affine(), vec![], vec![1.0]).is_err());
        assert!(
            plan.add(affine(), vec![Source::External(5)], vec![1.0])
                .is_err()
        );
        assert!(
            plan.add(affine(), vec![Source::External(0)], vec![])
                .is_err()
        );
        assert!(
            plan.add(affine(), vec![Source::Output(NodeId(3), 0)], vec![1.0])
                .is_err()
        );
        assert!(plan.result(Source::External(2)).is_err());
        let device = CpuDevice::sequential();
        assert!(plan.run(&device, 3, &[]).is_err());
    }
}
