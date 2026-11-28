//! Exact dense tensor transport through the public C ABI.
use num_bigint::BigInt;
use rssn::ffi::*;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use std::ffi::CStr;
use std::ptr;

struct Session(*mut RssnSession);

impl Session {
    fn new() -> Self {
        let raw = rssn_session_new();
        assert!(!raw.is_null());
        Self(raw)
    }

    fn stored(
        &self,
        op: &CStr,
        args: &[u32],
    ) -> u32 {
        let mut output = 0;
        // SAFETY: self retains the live session. Arguments and output are valid.
        assert_eq!(
            unsafe {
                rssn_value_call(
                    self.0,
                    op.as_ptr(),
                    args.as_ptr(),
                    args.len(),
                    &raw mut output,
                )
            },
            RssnStatus::Ok,
        );
        output
    }

    fn integer(
        &self,
        digits: &CStr,
    ) -> u32 {
        let mut output = 0;
        // SAFETY: self retains the session. Digits are terminated and output is writable.
        assert_eq!(
            unsafe { rssn_value_new_int_decimal(self.0, digits.as_ptr(), &raw mut output) },
            RssnStatus::Ok,
        );
        output
    }

    fn value(
        &self,
        cell: u32,
    ) -> Value {
        let mut raw = ptr::null_mut();
        // SAFETY: self retains the session, cell is live, and the output is writable.
        assert_eq!(
            unsafe { rssn_value_json(self.0, cell, &raw mut raw) },
            RssnStatus::Ok
        );
        assert!(!raw.is_null());
        // SAFETY: RSSN returned an owned NUL-terminated JSON string.
        let result = unsafe { CStr::from_ptr(raw) }.to_bytes();
        let decoded = serde_json::from_slice(result).expect("valid RSSN value JSON");
        // SAFETY: the JSON borrow ended before the owned RSSN string is freed.
        unsafe { rssn_string_free(raw) };
        decoded
    }

    fn tensor(
        &self,
        cells: &[u32],
        shape: &[usize],
    ) -> u32 {
        self.tensor_with_exemplar(cells, shape, 0)
    }

    fn tensor_with_exemplar(
        &self,
        cells: &[u32],
        shape: &[usize],
        exemplar: u32,
    ) -> u32 {
        let mut output = 0;
        // SAFETY: input slices hold the declared numbers of values.
        assert_eq!(
            unsafe {
                rssn_value_new_tensor_cells(
                    self.0,
                    cells.as_ptr(),
                    cells.len(),
                    shape.as_ptr(),
                    shape.len(),
                    exemplar,
                    &raw mut output,
                )
            },
            RssnStatus::Ok,
        );
        output
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        // SAFETY: Session owns this live pointer and releases it once.
        unsafe { rssn_session_free(self.0) }
    }
}

fn tensor(value: Value) -> TensorValue {
    let Value::Tensor(tensor) = value else {
        panic!("expected a tensor value");
    };
    tensor
}

#[test]
fn copies_exact_integers_and_rationals_in_row_major_order() {
    let root = Session::new();
    let a = root.integer(c"18446744073709551617");
    let b = root.integer(c"2");
    let c = root.integer(c"3");
    let d = root.integer(c"4");
    let one = root.integer(c"1");
    let matrix = root.tensor(&[a, b, c, d], &[2, 2]);
    let fraction = root.stored(c"divide", &[a, c]);
    let mixed = root.tensor(&[fraction, b, c, d], &[2, 2]);
    let fraction_value = root.value(fraction);
    for id in [a, b, c, d, fraction] {
        // SAFETY: the constructor copied these source-owned cells.
        unsafe { rssn_value_free(root.0, id) };
    }
    let wide = BigInt::from(18_446_744_073_709_551_617_u128);
    let packed = tensor(root.value(matrix));
    assert_eq!(packed.shape, [2, 2]);
    assert_eq!(packed.elem, TypeExpr::Int);
    assert_eq!(
        packed.data,
        vec![
            Value::int(wide.clone()),
            Value::int(2),
            Value::int(3),
            Value::int(4)
        ],
    );
    let mixed = tensor(root.value(mixed));
    assert_eq!(mixed.elem, TypeExpr::Real);
    assert_eq!(mixed.data.first(), Some(&fraction_value));
    let normalized = root.stored(c"num_native_qmatrix_new", &[matrix, one]);
    let determinant = root.stored(c"num_native_qlinalg_det", &[normalized]);
    assert_eq!(root.value(determinant), Value::int(wide * 4_i32 - 6_i32));
}

#[test]
fn homogeneous_measured_cells_carry_their_dimension_and_representation() {
    let s = Session::new();
    let mut source = 0;
    let mut measured = 0;
    let dim = br#"{"base":"Length"}"#;
    // SAFETY: all outputs are writable, the session is live and dimension JSON is a byte slice.
    unsafe {
        assert_eq!(
            rssn_value_new_float(s.0, -0.0, &raw mut source),
            RssnStatus::Ok
        );
        assert_eq!(
            rssn_value_with_coherent_dimension(
                s.0,
                source,
                dim.as_ptr(),
                dim.len(),
                &raw mut measured,
            ),
            RssnStatus::Ok,
        );
    }
    let result = tensor(s.value(s.tensor(&[measured], &[1, 1])));
    assert_eq!(result.shape, [1, 1]);
    assert_eq!(
        result.elem,
        TypeExpr::Float.with_dim(DimExpr::base("Length"))
    );
    let [Value::Dyn(value)] = result.data.as_slice() else {
        panic!("expected a measured Float");
    };
    assert!(matches!(*value.value, Value::Float(v) if v.to_bits() == (-0.0_f64).to_bits()));
}

#[test]
fn mixed_exact_quantities_share_their_dimension() {
    let s = Session::new();
    let mut unit = 0;
    // SAFETY: the session is live and the symbol name is NUL terminated.
    let metre = unsafe { rssn_term_sym(s.0, c"m".as_ptr()) };
    // SAFETY: metre belongs to s and the output is writable.
    assert_eq!(
        unsafe { rssn_value_new_term(s.0, metre, &raw mut unit) },
        RssnStatus::Ok
    );
    let one = s.integer(c"1");
    let wide = s.integer(c"18446744073709551617");
    let three = s.integer(c"3");
    let fraction = s.stored(c"divide", &[s.integer(c"2"), three]);
    let measured =
        [wide, fraction, three, s.integer(c"4")].map(|cell| s.stored(c"quantity", &[cell, unit]));
    let matrix = s.tensor(&measured, &[2, 2]);
    let packed = tensor(s.value(matrix));
    assert_eq!(
        packed.elem,
        TypeExpr::Real.with_dim(DimExpr::base("Length"))
    );
    let reference = s.stored(c"quantity", &[one, unit]);
    let normalized = s.stored(c"num_native_qmatrix_new", &[matrix, reference]);
    let determinant = s.stored(c"num_native_qlinalg_det", &[normalized]);
    assert_eq!(
        s.value(determinant),
        Value::int(73_786_976_294_838_206_466_u128)
    );
    let empty = tensor(s.value(s.tensor_with_exemplar(&[], &[0, 2], measured[1])));
    assert_eq!(empty.elem, packed.elem);
    assert_eq!(empty.shape, [0, 2]);
    assert!(empty.data.is_empty());
}

#[test]
fn validates_extents_and_cell_liveness_without_mutating_the_bank_on_failure() {
    let s = Session::new();
    let foreign = Session::new();
    let foreign_id = foreign.integer(c"6");
    let stale = s.integer(c"11");
    // SAFETY: s owns this live cell and no longer needs it.
    unsafe { rssn_value_free(s.0, stale) };
    let cell = s.integer(c"7");
    assert!(tensor(s.value(s.tensor(&[cell], &[]))).shape.is_empty());
    let empty = tensor(s.value(s.tensor_with_exemplar(&[], &[0, 9, 2], cell)));
    assert_eq!(empty.elem, TypeExpr::Int);
    assert_eq!(empty.shape, [0, 9, 2]);
    assert!(empty.data.is_empty() && empty.is_consistent());
    let zero_extent = tensor(s.value(s.tensor(&[], &[usize::MAX, 2, 0])));
    assert!(zero_extent.is_consistent());
    assert_eq!(zero_extent.extent(), 0);

    // SAFETY: sessions live throughout; each non-null input points at exactly the claimed count.
    unsafe {
        let baseline = rssn_value_count(s.0);
        let singleton = [1usize];
        for (input, count, shape, rank, status) in [
            (&raw const cell, 1, ptr::null(), 1, RssnStatus::NullArgument),
            (
                ptr::null(),
                1,
                singleton.as_ptr(),
                1,
                RssnStatus::NullArgument,
            ),
        ] {
            let mut output = u32::MAX;
            assert_eq!(
                rssn_value_new_tensor_cells(s.0, input, count, shape, rank, 0, &raw mut output),
                status,
            );
            assert_eq!(output, u32::MAX);
        }
        for (cells, shape, status) in [
            (&[cell][..], &[2, 1][..], RssnStatus::InvalidArgument),
            (&[][..], &[usize::MAX, 2][..], RssnStatus::InvalidArgument),
            (&[][..], &[][..], RssnStatus::InvalidArgument),
            (&[stale][..], &[1][..], RssnStatus::NoValue),
            // A low-level cell ID is interpreted in the supplied session's bank.
            (&[foreign_id][..], &[1][..], RssnStatus::NoValue),
        ] {
            let mut output = u32::MAX;
            assert_eq!(
                rssn_value_new_tensor_cells(
                    s.0,
                    cells.as_ptr(),
                    cells.len(),
                    shape.as_ptr(),
                    shape.len(),
                    0,
                    &raw mut output,
                ),
                status,
            );
            assert_eq!(output, u32::MAX);
        }
        assert_eq!(
            rssn_value_new_tensor_cells(
                s.0,
                &raw const cell,
                1,
                ptr::null(),
                0,
                0,
                ptr::null_mut()
            ),
            RssnStatus::NullArgument,
        );
        let empty_shape = [0usize, 2];
        for (exemplar, status) in [
            (stale, RssnStatus::NoValue),
            (foreign_id, RssnStatus::NoValue),
        ] {
            let mut output = u32::MAX;
            assert_eq!(
                rssn_value_new_tensor_cells(
                    s.0,
                    ptr::null(),
                    0,
                    empty_shape.as_ptr(),
                    empty_shape.len(),
                    exemplar,
                    &raw mut output,
                ),
                status,
            );
            assert_eq!(output, u32::MAX);
        }
        let mut output = u32::MAX;
        assert_eq!(
            rssn_value_new_tensor_cells(
                s.0,
                &raw const cell,
                1,
                ptr::null(),
                0,
                cell,
                &raw mut output,
            ),
            RssnStatus::InvalidArgument,
        );
        assert_eq!(output, u32::MAX);
        assert_eq!(rssn_value_count(s.0), baseline);
    }
}
