"""Unit tests for the Python binding.

    python -m unittest discover -s bindings/python/tests -v

They need a built rssn library (see `bindings/README.md`); the conformance table
in `bindings/tests/cases.tsv` is the cross-language half of the suite.
"""

from __future__ import annotations

import ctypes
import json
import math
import os
import pathlib
import sys
import tempfile
import threading
import unittest
from fractions import Fraction
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import rssn  # noqa: E402


class LibraryTests(unittest.TestCase):
    def test_source_checkout_searches_bazel_before_cargo(self) -> None:
        from rssn import _ffi

        with tempfile.TemporaryDirectory() as temp:
            root = pathlib.Path(temp)
            package = root / "bindings" / "python" / "rssn"
            (root / "rssn.h").touch()
            (root / "MODULE.bazel").touch()
            with patch.object(_ffi, "_HERE", str(package)), patch.object(_ffi, "_REPO", str(root)):
                candidates = list(map(pathlib.Path, _ffi._candidate_dirs()))
            targets = ("ffi_distribution", "ffi", "ffi_distribution_portable", "ffi_dev")
            indices = [candidates.index(root / "bazel-bin" / target / "ffi") for target in targets]
            self.assertEqual(indices, sorted(indices))
            self.assertLess(indices[-1], candidates.index(root / "target" / "debug"))

    def test_platform_loader_is_tried_after_local_directories(self) -> None:
        # Mock only the system loader and ABI declarations. Discovery should
        # still resolve a platform library when the checkout has no copy.
        from rssn import _ffi

        with patch.dict(os.environ, {"RSSN_LIB": ""}), patch.object(
            _ffi, "_candidate_dirs", return_value=[]
        ), patch.object(_ffi, "_bind"), patch.object(_ffi.ctypes, "CDLL") as loader, patch.object(
            _ffi, "_lib", None
        ), patch.object(_ffi, "_lib_path", None):
            self.assertIs(rssn.load(), loader.return_value)
            loader.assert_called_once_with(_ffi._LIB_NAMES[0])
            self.assertEqual(rssn.library_path(), _ffi._LIB_NAMES[0])

    def test_loading_a_directory_and_then_a_bad_library(self) -> None:
        from rssn import _ffi

        with tempfile.TemporaryDirectory() as directory, patch.object(
            _ffi, "_lib", None
        ), patch.object(_ffi, "_lib_path", None):
            candidate = pathlib.Path(directory) / _ffi._LIB_NAMES[0]
            candidate.touch()
            with patch.object(_ffi, "_bind"), patch.object(_ffi.ctypes, "CDLL") as loader:
                bound = rssn.load(directory)
                loader.assert_called_once_with(str(candidate))
                self.assertEqual(rssn.library_path(), str(candidate))
            with self.assertRaisesRegex(OSError, "cannot load"):
                rssn.load(candidate)
            self.assertIs(_ffi.library(), bound)


class SessionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.session = rssn.Session()
        self.x = self.session.sym("x")

    def tearDown(self) -> None:
        self.session.close()

    def test_version_and_manifest(self) -> None:
        self.assertTrue(rssn.version().startswith("3."))
        self.assertNotEqual(rssn.manifest_hash(), 0)
        manifest = rssn.manifest()
        self.assertIn("ops", manifest)
        self.assertEqual(len(manifest["ops"]), len(rssn.OPS))

    def test_terms_render(self) -> None:
        expression = self.x * self.x + 2 * self.x + 1
        # A term is a request: `x*x` stays a product until it is computed.
        self.assertEqual(str(expression), "x*x + 2*x + 1")
        self.assertEqual(self.session.compute(expression).text, "x^2 + 2*x + 1")
        self.assertIn("x", expression.latex)
        self.assertTrue(expression.typst)
        self.assertTrue(expression.pretty)

    def test_arithmetic_sugar(self) -> None:
        # The engine has no `sub`, `neg` or arithmetic `div` operator, so the
        # sugar builds `add(a, mul(-1, b))` and `mul(a, pow(b, -1))`.
        x = self.x
        self.assertEqual(str(x - 1), "x - 1")
        self.assertEqual(str(1 - x), "1 - x")
        self.assertEqual(str(-x), "-x")
        self.assertEqual(str(x / 2), "x/2")
        self.assertEqual(str(2 / x), "2/x")
        self.assertEqual(str(x**2), "x^2")
        self.assertEqual(str(2**x), "2^x")
        self.assertEqual(str(2 * x + 1), "2*x + 1")
        self.assertEqual(self.session.compute((x**2 - 1) / (x - 1)).text, "x + 1")

    def test_named_arithmetic_operators(self) -> None:
        # `sub`, `neg` and `divide` are engine operators now, and the sugar
        # goes through them.
        self.assertEqual(self.session.apply("sub", self.x, 1).text, "x - 1")
        self.assertEqual(self.session.apply("neg", self.x).text, "-x")
        self.assertEqual(self.session.apply("divide", self.x, 2).text, "x/2")
        self.assertEqual(self.session.compute(self.session.apply("sub", 5, 3)).text, "2")
        self.assertEqual(self.session.compute(self.session.apply("divide", 6, 3)).text, "2")
        self.assertEqual(
            self.session.compute((self.x**2 - 1) / (self.x - 1)).text, "x + 1"
        )

    def test_parse_and_literal_values(self) -> None:
        self.assertEqual(self.session.parse("1/3").text, "1/3")
        self.assertEqual(self.session.rational(1, 3).text, "1/3")
        self.assertEqual(self.session.int(7).value, 7.0)
        self.assertIsNone(self.x.value)
        self.assertEqual(self.session.term(Fraction(2, 5)).text, "2/5")
        self.assertEqual(self.session.term("x + 1").text, "x + 1")
        self.assertEqual(self.session.term(2).text, "2")

    def test_compute_symbolic(self) -> None:
        answer = self.session.compute(rssn.diff(rssn.sin(self.x) * rssn.cos(self.x), self.x))
        self.assertEqual(answer.text, "cos(2*x)")
        self.assertTrue(answer.reduced)
        self.assertIsNone(answer.value)
        self.assertEqual(str(answer), "cos(2*x)")

    def test_compute_numeric(self) -> None:
        config = rssn.Config().numeric(1e-12).bind("x", 0.0)
        answer = self.session.compute(rssn.diff(rssn.sin(self.x) * rssn.cos(self.x), self.x), config)
        self.assertAlmostEqual(answer.value, 1.0, places=12)
        self.assertAlmostEqual(answer.float(), 1.0, places=12)

    def test_float_of_term(self) -> None:
        self.assertAlmostEqual(float(self.session.parse("2 + 3")), 5.0)

    def test_tensor(self) -> None:
        matrix = self.session.parse(
            "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))"
        )
        entries = self.session.compute(matrix).tensor()
        shape = getattr(entries, "shape", None)
        if shape is not None:  # numpy installed
            self.assertEqual(tuple(shape), (2, 2))
            self.assertEqual(float(entries[0][1]), 1.0)
        else:
            self.assertEqual(entries, [[2.0, 1.0], [4.0, 3.0]])

    def test_compute_all(self) -> None:
        answers = self.session.compute_all(["diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"])
        self.assertEqual(len(answers), 3)
        self.assertTrue(all(isinstance(answer, rssn.Answer) for answer in answers))

    def test_batch_dimension_findings_follow_the_request(self) -> None:
        with rssn.Session() as session:
            requests = ["quantity(3, m) + quantity(4, s)", "quantity(2, m)"]
            individual = [session.compute(request) for request in requests]
            batch = session.compute_all(requests)
            self.assertIn("mismatch", individual[0].diagnostic_codes)
            self.assertEqual(
                [(answer.dimension, answer.diagnostic_codes) for answer in batch],
                [(answer.dimension, answer.diagnostic_codes) for answer in individual],
            )
            off = rssn.Config().dimension_mode("off")
            self.assertTrue(all(answer.dimension is None and not answer.diagnostic_codes
                                for answer in session.compute_all(requests, off)))

    def test_compute_all_collects_failures(self) -> None:
        # A numeric target on a symbolic answer is a per-request failure.
        config = rssn.Config().numeric(1e-12)
        answers = self.session.compute_all(["1 + 1", "diff(x^2, x)"], config, errors="collect")
        self.assertIsInstance(answers[0], rssn.Answer)
        self.assertIsInstance(answers[1], rssn.RssnError)
        self.assertEqual(answers[1].status, rssn.Status.NOT_NUMERIC)

    def test_compute_all_rejects_an_unknown_error_policy(self) -> None:
        with self.assertRaisesRegex(ValueError, "errors must be"):
            self.session.compute_all([self.x], errors="typo")

    def test_compute_all_raises(self) -> None:
        config = rssn.Config().numeric(1e-12)
        with self.assertRaises(rssn.RssnError):
            self.session.compute_all(["diff(x^2, x)"], config)

    def test_compute_json_record(self) -> None:
        record = self.session.compute_json("diff(sin(x), x)")
        for field in ("term", "value", "error", "reduced", "iterations", "dimension", "diagnostics"):
            self.assertIn(field, record)
        self.assertEqual(record["term"], "cos(x)")

    def test_compute_bincode_is_bytes(self) -> None:
        blob = self.session.compute_bincode("diff(sin(x), x)")
        self.assertIsInstance(blob, bytes)
        self.assertGreater(len(blob), 4)
        batch = self.session.compute_all_bincode(["diff(x^2, x)", "gcd(2, 4)"])
        self.assertIsInstance(batch, bytes)

    def test_eval_without_rules(self) -> None:
        self.assertAlmostEqual(self.session.eval("sin(x)", {"x": 0.0}), 0.0)

    def test_dimension(self) -> None:
        config = rssn.Config().dimension("x", "m")
        self.assertTrue(self.session.dimension("diff(x^2, x)", config))
        self.assertEqual(self.session.dimension("2 + 2"), "1")

    def test_strict_dimension_conflict(self) -> None:
        config = rssn.Config().dimension("x", "m").dimension("y", "s").dimension_mode("strict")
        with self.assertRaises(rssn.RssnError) as caught:
            self.session.compute("x + y", config)
        self.assertEqual(caught.exception.status, rssn.Status.DIMENSION)

    def test_custom_base_dimension(self) -> None:
        config = rssn.Config().define_base_dimension("Currency").dimension("price", "Currency")
        self.assertEqual(self.session.dimension("price", config), "Currency")

    def test_budget(self) -> None:
        config = rssn.Config().budget(max_iterations=1, max_nodes=64)
        self.session.compute("diff(sin(x), x)", config)  # must not raise

    def test_assumptions(self) -> None:
        config = rssn.Config().assume("n", rssn.Assumption.INTEGER)
        self.session.compute("diff(x^n, x)", config)

    def test_closed_session_is_an_error(self) -> None:
        session = rssn.Session()
        session.close()
        with self.assertRaises(RuntimeError):
            _ = session.ptr


class ErrorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.session = rssn.Session()

    def tearDown(self) -> None:
        self.session.close()

    def test_parse_error(self) -> None:
        with self.assertRaises(rssn.RssnError) as caught:
            self.session.parse("sin(")
        self.assertEqual(caught.exception.status, rssn.Status.PARSE)

    def test_unknown_operator(self) -> None:
        with self.assertRaises(rssn.RssnError) as caught:
            self.session.apply("definitely_not_an_operator", 1)
        self.assertEqual(caught.exception.status, rssn.Status.UNKNOWN_OPERATOR)

    def test_arity_is_checked_before_the_call(self) -> None:
        with self.assertRaises(rssn.RssnError) as caught:
            self.session.apply("gcd", 1)
        self.assertIn("arguments", caught.exception.message)

    def test_not_numeric(self) -> None:
        config = rssn.Config().numeric(1e-12)
        with self.assertRaises(rssn.RssnError) as caught:
            self.session.compute("diff(x^2, x)", config)
        self.assertEqual(caught.exception.status, rssn.Status.NOT_NUMERIC)

    def test_zero_denominator(self) -> None:
        with self.assertRaises(rssn.RssnError):
            self.session.rational(1, 0)


class FacadeTests(unittest.TestCase):
    def test_keyword_names_are_mangled(self) -> None:
        self.assertTrue(callable(rssn.and_))
        self.assertEqual(rssn.NAMES["and"], "and_")
        self.assertIn("and", rssn.OPS)

    def test_zero_arity_operator(self) -> None:
        session = rssn.default_session()
        term = rssn.E(session=session)
        self.assertIsInstance(term, rssn.Term)

    def test_variadic_operator(self) -> None:
        session = rssn.Session()
        try:
            total = rssn.add(session.int(1), session.int(2), session.int(3))
            self.assertEqual(session.compute(total).text, "6")
        finally:
            session.close()

    def test_alias_forwards(self) -> None:
        self.assertEqual(rssn.NAMES.get("airy_ai", "airyai"), "airyai")

    def test_groups(self) -> None:
        self.assertIn("diff", rssn.BY_RULE_SET["calculus"])
        self.assertTrue(rssn.BY_CATEGORY)

    def test_default_session_is_per_thread(self) -> None:
        main_session = rssn.default_session()
        seen: list[int] = []

        def worker() -> None:
            seen.append(id(rssn.default_session()))

        thread = threading.Thread(target=worker)
        thread.start()
        thread.join()
        self.assertNotEqual(seen[0], id(main_session))

    def test_use_session(self) -> None:
        session = rssn.Session()
        try:
            with rssn.use_session(session):
                self.assertIs(rssn.default_session(), session)
            self.assertIsNot(rssn.default_session(), session)
        finally:
            session.close()


class CheckedFacadeTests(unittest.TestCase):
    def test_dynamic_facade_accepts_native_arguments_and_preserves_term_calls(self) -> None:
        with rssn.Session() as session:
            request = rssn.num_special_ext_lerch(-1234.0, 2.0, 1.0, session=session)
            self.assertIsInstance(request, rssn.Term)
            self.assertTrue(str(request).startswith("num_special_ext_lerch("))
            with rssn.num_special_ext_lerch(-1234.0, 2.0, 1.0, session=session, checked=True) as lerch:
                self.assertAlmostEqual(lerch.float(), 0.021861589304464496, places=15)
            with rssn.num_special_ext_polylog(2.0, -1234.0, session=session, checked=True) as polylog:
                self.assertAlmostEqual(polylog.float() / -1234.0, 0.021861589304464496, places=15)
            with rssn.num_uhp_to_disk(1j, session=session, checked=True) as point:
                self.assertEqual(point.json()["complex"]["float"], {"re": 0.0, "im": 0.0})
                with rssn.num_disk_to_uhp(point, checked=True) as recovered:
                    self.assertIs(recovered.session, session)
                    self.assertEqual(recovered.json()["complex"]["float"], {"re": 0.0, "im": 1.0})
            with self.assertRaises(rssn.RssnError) as error:
                rssn.num_disk_to_uhp(1.0 + 0j, session=session, checked=True)
            self.assertEqual(error.exception.status, rssn.Status.INVALID_ARGUMENT)


class ReferenceTests(unittest.TestCase):
    """The operator reference: the manifest's documentation, searchable."""

    def test_op_help_prints_the_operator(self) -> None:
        text = rssn.op_help("diff")
        self.assertIn("rssn.diff(f, x", text)
        self.assertIn("calculus", text)
        self.assertIn("derivative", text)

    def test_op_help_refuses_an_unknown_name(self) -> None:
        with self.assertRaises(KeyError):
            rssn.op_help("definitely_not_an_operator")

    def test_op_search(self) -> None:
        matches = rssn.op_search("brent", limit=3)
        self.assertTrue(matches)
        self.assertTrue(all(isinstance(pair, tuple) for pair in matches))
        self.assertIn("num_brent", [name for name, _ in rssn.op_search("brent", limit=50)])

    def test_op_search_filters(self) -> None:
        calculus = rssn.op_search(category="calculus", limit=5)
        self.assertTrue(calculus)
        for name, _ in calculus:
            self.assertIn(rssn.NAMES[name], rssn.BY_CATEGORY["calculus"])
        self.assertEqual(rssn.op_search(rule_set="no_such_rule_set"), [])

    def test_the_tables(self) -> None:
        self.assertCountEqual(rssn.op_names(), [op["name"] for op in rssn.manifest()["ops"]])
        self.assertEqual(len(rssn.op_categories()), 40)
        self.assertEqual(len(rssn.op_rule_sets()), 40)


class BackendTests(unittest.TestCase):
    def test_current_backend_queries_retire_transient_handles(self) -> None:
        with rssn.Session() as session:
            def objects():
                with session.call_value("num_native_usage") as usage, usage["objects"] as count:
                    return int(count)

            before = objects()
            self.assertTrue(session.current_backend())
            self.assertEqual(objects(), before)
            self.assertEqual(session.value_count(), 0)

    def test_process_selection_and_explicit_handles(self) -> None:
        with rssn.Session() as session, rssn.Session() as other:
            self.assertIn("interpreter", session.available_backends())
            with session.open_backend("interpreter") as handle:
                try:
                    with session.call_value("num_backend_name", handle) as name:
                        self.assertEqual(name.text(), "interpreter")
                finally:
                    session.release_backend(handle)
                with self.assertRaises(rssn.RssnError):
                    session.call_value("num_backend_name", handle)
            try:
                self.assertEqual(session.set_default_backend("interpreter"), "interpreter")
                self.assertEqual(other.current_backend(), "interpreter")
                with self.assertRaises(rssn.RssnError):
                    other.set_default_backend("not-an-rssn-backend")
                self.assertEqual(other.current_backend(), "interpreter")
                with self.assertRaises(ValueError):
                    other.set_default_backend("interpreter\0ignored")
            finally:
                self.assertEqual(session.set_default_backend("auto"), other.current_backend())

    def test_all_usable_backends_compile_a_measured_program(self) -> None:
        with rssn.Session() as session:
            for selector in session.available_backends():
                with self.subTest(selector=selector), session.open_backend(selector) as handle:
                    with session.call_value("num_backend_reference", 1.0) as reference:
                        with session.call_value("num_backend_value", 3.0) as sample:
                            try:
                                with session.call_value(
                                    "num_backend_compile", handle, "x*x+1",
                                    ["x"], [reference], reference, 16,
                                ) as program:
                                    try:
                                        with session.call_value("num_backend_call", program, [sample]) as result:
                                            self.assertAlmostEqual(result.float(), 10.0, places=5)
                                    finally:
                                        session.release_backend(program)
                            finally:
                                session.release_backend(sample)
                                session.release_backend(reference)
                                session.release_backend(handle)


class NpyTests(unittest.TestCase):
    def test_rank_and_ieee_samples_round_trip(self) -> None:
        if not rssn.manifest()["features"].get("npy"):
            self.skipTest("this build lacks the NumPy codec")
        with rssn.Session() as session:
            for shape, data in [
                ((2, 3), [-0.0, 1.0, 2.0, 3.0, float("inf"), float("nan")]),
                ((), [-0.0]),
                ((2, 0, 3), []),
            ]:
                with self.subTest(shape=shape):
                    payload = session.encode_npy(shape, data)
                    self.assertTrue(payload.startswith(b"\x93NUMPY"))
                    with session.decode_npy(memoryview(payload)) as array:
                        extents, flat = array._tensor_parts()
                        self.assertEqual(extents, list(shape))
                        self.assertEqual(len(flat), len(data))
                        for actual, expected in zip(flat, data):
                            if math.isnan(expected):
                                self.assertTrue(math.isnan(actual))
                            else:
                                self.assertEqual(actual, expected)
                                self.assertEqual(math.copysign(1.0, actual),
                                                 math.copysign(1.0, expected))
            with self.assertRaises(rssn.RssnError):
                session.decode_npy(b"truncated")

    def test_complex_samples_are_rejected(self) -> None:
        if not rssn.manifest()["features"].get("npy"):
            self.skipTest("this build lacks the NumPy codec")
        with rssn.Session() as session:
            with self.assertRaisesRegex(TypeError, "requires real samples"):
                session.encode_npy((1,), [1 + 2j])
            self.assertEqual(session.value_count(), 0)

    def test_numpy_scalar_types_do_not_lose_imaginary_parts(self) -> None:
        if not rssn.manifest()["features"].get("npy"):
            self.skipTest("this build lacks the NumPy codec")
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is not installed")

        with rssn.Session() as session:
            for sample in (np.complex64(1 + 2j), np.complex128(1 + 0j)):
                with self.subTest(dtype=sample.dtype), self.assertRaisesRegex(
                    TypeError, "requires real samples"
                ):
                    session.encode_npy((1,), iter([sample]))
            values = np.array([1.25, -0.0], dtype=np.float32)
            with session.decode_npy(session.encode_npy(values.shape, values)) as decoded:
                np.testing.assert_array_equal(decoded.tensor(), values)
                np.testing.assert_array_equal(np.signbit(decoded.tensor()), np.signbit(values))
            self.assertEqual(session.value_count(), 0)

    def test_read_fortran_order_float32(self) -> None:
        if not rssn.manifest()["features"].get("npy"):
            self.skipTest("this build lacks the NumPy codec")
        import struct

        header = b"{'descr': '<f4', 'fortran_order': True, 'shape': (2, 3), }"
        header += b" " * ((-10 - len(header) - 1) % 64) + b"\n"
        raw = b"\x93NUMPY\x01\x00" + struct.pack("<H", len(header)) + header
        raw += struct.pack("<6f", 1, 4, 2, 5, 3, 6)
        with rssn.Session() as session, session.decode_npy(raw) as array:
            self.assertEqual(array._tensor_parts(), ([2, 3], [1, 2, 3, 4, 5, 6]))

    def test_exchange_npy_files_with_numpy(self) -> None:
        if not rssn.manifest()["features"].get("npy"):
            self.skipTest("this build lacks the NumPy codec")
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is not installed")
        from io import BytesIO

        samples = (
            np.array(-0.0, dtype=">f8"),
            np.empty((2, 0, 3), dtype="<f4"),
            np.array(np.arange(60).reshape(3, 4, 5), dtype=">f4", order="F"),
        )
        with rssn.Session() as session:
            for sample in samples:
                with self.subTest(shape=sample.shape, dtype=sample.dtype):
                    numpy_file = BytesIO()
                    np.save(numpy_file, sample, allow_pickle=False)
                    with session.decode_npy(numpy_file.getbuffer()) as value:
                        decoded = value.tensor()
                    np.testing.assert_array_equal(decoded, sample)
                    payload = session.encode_npy(decoded.shape, decoded.flat)
                    restored = np.load(BytesIO(payload), allow_pickle=False)
                    self.assertEqual(restored.dtype, np.dtype("<f8"))
                    np.testing.assert_array_equal(restored, sample)
                    np.testing.assert_array_equal(np.signbit(restored), np.signbit(sample))


class ValueTests(unittest.TestCase):
    """The value ABI: operators whose arguments have no term form."""

    def setUp(self) -> None:
        self.session = rssn.Session()

    def tearDown(self) -> None:
        self.session.close()

    def test_numpy_arrays_keep_rank_and_empty_extents(self) -> None:
        try:
            import numpy as np
        except ImportError:
            self.skipTest("NumPy is not installed")

        arrays = (np.array(3.5), np.empty((0, 2, 3)), np.arange(24.0).reshape(2, 3, 4).T)
        for array in arrays:
            with rssn.Value.of(self.session, array) as stored:
                actual = stored.tensor()
                self.assertEqual(actual.shape, array.shape)
                np.testing.assert_array_equal(actual, array)

    def test_integer_cells_keep_values_outside_i64(self) -> None:
        for original in (2**100 + 7, -(2**100 + 7)):
            with rssn.Value.of(self.session, original) as exact:
                self.assertEqual(exact.kind(), "int")
                self.assertEqual(int(exact.text()), original)
                self.assertEqual(int(exact), original)
            self.assertEqual(int(self.session.int(original).text), original)

    def test_integer_cells_round_trip_at_32768_bits(self) -> None:
        wide = (1 << 32767) + 123
        digit_limit = getattr(sys, "get_int_max_str_digits", lambda: None)
        before = digit_limit()
        for original in (wide, -wide):
            with rssn.Value.of(self.session, original) as exact:
                self.assertEqual(exact.kind(), "int")
                self.assertEqual(int(exact), original)
                self.assertEqual(self.session.int(original).text, exact.text())
        self.assertEqual(digit_limit(), before)

    def test_backend_handle(self) -> None:
        backend = self.session.call_value("num_backend_current")
        self.assertEqual(backend.kind(), "opaque")
        name = self.session.call_value("num_backend_name", backend)
        self.assertEqual(name.kind(), "string")
        self.assertEqual(name.text(), self.session.current_backend())
        # The handle has no term form.
        with self.assertRaises(rssn.RssnError):
            backend.to_term()

    def test_ad_tape(self) -> None:
        tape = self.session.call_value("num_ad_tape_new", 16)
        x = self.session.call_value("num_ad_tape_var", tape, 1.5)
        z = self.session.call_value(
            "num_ad_var_mul", self.session.call_value("num_ad_var_sin", x), x
        )
        self.assertAlmostEqual(self.session.call_value("num_ad_var_value", z).float(),
                               math.sin(1.5) * 1.5, places=12)
        gradient = self.session.call_value("num_ad_var_backward", z)
        self.assertAlmostEqual(
            self.session.call_value("num_ad_var_grad_in", x, gradient).float(),
            math.cos(1.5) * 1.5 + math.sin(1.5), places=12,
        )
        self.assertEqual(self.session.call_value("num_ad_tape_len", tape).float(), 3.0)

    def test_value_json_and_lists(self) -> None:
        backend = self.session.call_value("num_backend_current")
        self.assertIsInstance(backend.json(), dict)
        # A list argument keeps its element kinds.
        values = rssn.Value.of(self.session, [1, 2.5, 3])
        self.assertEqual(values.kind(), "tensor")
        self.assertEqual(len(values), 3)
        self.assertAlmostEqual(values[1].float(), 2.5)
        with self.assertRaises(TypeError):
            values.field_names()
        with self.assertRaises(TypeError):
            values["a"]

    def test_a_dict_becomes_a_record(self) -> None:
        from rssn.callback import to_python

        record = rssn.Value.of(self.session, {"a": 1.0, "b": "two"})
        self.assertEqual(record.kind(), "record")
        self.assertEqual(len(record), 2)
        self.assertEqual(record[0].text(), "1")
        self.assertEqual(record.json()["record"]["fields"][1][0], "b")
        self.assertEqual(to_python(self.session, record.cell), {"a": 1.0, "b": "two"})

    def test_rolling_moments_record_supports_field_names(self) -> None:
        with self.session.call_value(
            "num_timeseries_rolling_moments", [1.0, 3.0, 5.0], 2, 0
        ) as result:
            self.assertEqual(result.field_names(), ("mean", "variance"))
            with result["mean"] as mean, result["variance"] as variance:
                self.assertEqual(list(mean.tensor()), [2.0, 4.0])
                self.assertEqual(list(variance.tensor()), [1.0, 1.0])
            with result[0] as first:
                self.assertEqual(list(first.tensor()), [2.0, 4.0])
            with self.assertRaises(KeyError) as missing:
                result["absent"]
            self.assertEqual(missing.exception.args, ("absent",))

    def test_eval_term(self) -> None:
        term = self.session.parse("sin(x) + 1")
        self.assertAlmostEqual(self.session.eval("sin(x) + 1", {"x": 0.0}), 1.0)
        self.assertAlmostEqual(float(self.session.parse("2 + 3")), 5.0)
        self.assertIsNotNone(term)

    def test_typed_handle_is_checked(self) -> None:
        tape = self.session.call_value("num_ad_tape_new", 16)
        var = self.session.call_value("num_ad_tape_var", tape, 1.0)
        with self.assertRaises(rssn.RssnError):
            # A variable is not a gradient.
            self.session.call_value("num_ad_var_grad_in", var, var)

    def test_term_round_trip(self) -> None:
        value = self.session.call_value("num_backend_name",
                                        self.session.call_value("num_backend_current"))
        self.assertEqual(value.kind(), "string")
        cell = rssn.Value.of(self.session, 2.5)
        self.assertEqual(cell.kind(), "float")
        self.assertAlmostEqual(cell.float(), 2.5)

    def test_cells_are_released(self) -> None:
        session = rssn.Session()
        self.assertEqual(session.value_count(), 0)
        value = session.call_value("num_backend_current")
        self.assertEqual(session.value_count(), 1)
        value.close()
        self.assertEqual(session.value_count(), 0)
        session.close()


class CallbackTests(unittest.TestCase):
    """Host callbacks: a Python callable the engine calls."""

    def setUp(self) -> None:
        self.session = rssn.Session()

    def tearDown(self) -> None:
        self.session.close()

    @staticmethod
    def residual(x: float) -> float:
        """The residual whose positive root is sqrt(2)."""
        return x * x - 2.0

    def test_a_callable_argument_becomes_a_callback(self) -> None:
        answer = self.session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, self.residual)
        self.assertAlmostEqual(float(answer[0]), math.sqrt(2.0), places=11)
        self.assertGreater(float(answer[1]), 0.0)

    def test_the_generated_facade_takes_callables(self) -> None:
        answer = rssn.num_brent(1.0, 2.0, 1e-12, 64, self.residual, session=self.session)
        self.assertAlmostEqual(float(answer[0]), math.sqrt(2.0), places=11)

    def test_a_callback_receives_native_complex_arguments(self) -> None:
        center = 0.3 + 0.4j
        sampled = []

        def square(argument):
            sampled.append(argument)
            return argument * argument

        with rssn.num_complex_derivative_scaled(center, 1.0, 1.0, 64, square, session=self.session) as answer:
            parts = answer.json()["complex"]["float"]
            self.assertAlmostEqual(complex(parts["re"], parts["im"]), 2 * center, places=8)
        self.assertTrue(sampled)
        self.assertTrue(all(isinstance(z, complex) and abs(z - center) < 1e-4 for z in sampled))
        self.assertEqual(self.session.callback_count(), 0)

    def test_callback_complex_reader_checks_wire_parts(self) -> None:
        from rssn.callback import to_python

        with rssn.Value.of(self.session, 1j) as argument:
            for tag in ("float", "Float"):
                with patch.object(rssn.Value, "json", return_value={"complex": {tag: {"re": 1.0, "im": -2.0}}}):
                    self.assertEqual(to_python(self.session, argument.cell), 1.0 - 2.0j)
            for tag in ("exact", "Exact"):
                payload = {"re": {"Int": 3}, "im": {"Rat": [1, 4]}}
                with patch.object(rssn.Value, "json", return_value={"complex": {tag: payload}}):
                    self.assertEqual(to_python(self.session, argument.cell), 3.0 + 0.25j)
            with patch.object(rssn.Value, "json", return_value={"complex": {"float": {"re": 1.0}}}):
                with self.assertRaisesRegex(KeyError, "im"):
                    to_python(self.session, argument.cell)

    def test_a_callback_can_read_an_exact_large_integer(self) -> None:
        from rssn.callback import to_python

        integer = (1 << 32768) - 1
        with rssn.Value.of(self.session, integer) as argument:
            self.assertEqual(to_python(self.session, argument.cell), integer)

    def test_callback_can_return_a_retained_value(self) -> None:
        with rssn.Value.of(self.session, 0.0) as zero:
            with self.session.callback(lambda _: zero, signature="(float) -> float") as residual:
                answer = rssn.num_brent(1.0, 2.0, 1e-12, 64, residual)
                self.assertIs(answer.session, self.session)
                self.assertIn(answer[0].float(), (1.0, 2.0))
                self.assertEqual(zero.float(), 0.0)
                with zero.clone() as duplicate:
                    self.assertNotEqual(duplicate.cell, zero.cell)
                    self.assertEqual(duplicate.float(), 0.0)

    def test_an_explicit_callback_outlives_one_call(self) -> None:
        with self.session.callback(self.residual, signature="(float) -> float") as residual:
            first = self.session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, residual)
            second = self.session.call_value("num_brent", 1.0, 2.0, 1e-10, 32, residual)
        self.assertAlmostEqual(float(first[0]), math.sqrt(2.0), places=11)
        self.assertAlmostEqual(float(second[0]), math.sqrt(2.0), places=8)
        self.assertEqual(self.session.callback_count(), 0)

    def test_signatures_are_listed(self) -> None:
        listed = rssn.callback_signatures(self.session, "num_brent", 4)
        self.assertEqual(len(listed), 1)
        self.assertEqual(listed[0]["params"], ["float"])
        self.assertEqual(listed[0]["ret"], "float")
        # An argument that is not a function lists nothing.
        self.assertEqual(rssn.callback_signatures(self.session, "num_brent", 0), [])

    def test_parse_signature(self) -> None:
        self.assertEqual(
            rssn.parse_signature("(real, real) -> real"),
            {"params": ["real", "real"], "ret": "real"},
        )
        self.assertEqual(rssn.parse_signature("() -> bool")["params"], [])
        vector = rssn.parse_signature("float[] -> float")
        self.assertEqual(vector["params"][0]["tensor"]["rank"], {"lit": 1})
        self.assertEqual(vector["params"][0]["tensor"]["elem"], "float")
        # JSON passes through.
        payload = {"params": ["int"], "ret": "bool"}
        self.assertEqual(rssn.parse_signature(json.dumps(payload)), payload)

    def test_a_callback_that_raises_aborts_the_operator(self) -> None:
        def refuses(_x: float) -> float:
            raise ValueError("no thanks")

        with self.assertRaises(rssn.RssnError) as raised:
            self.session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, refuses)
        self.assertEqual(raised.exception.status, rssn.Status.INVALID_ARGUMENT)
        self.assertIn("no thanks", raised.exception.message)
        # The implicit callback is retired even though the call failed.
        self.assertEqual(self.session.callback_count(), 0)

    def test_a_callback_receives_native_objects(self) -> None:
        seen: list[object] = []

        def objective(point):
            seen.append(type(point).__name__)
            return float(sum((x - 1.0) ** 2 for x in point))

        def gradient(point):
            return [2.0 * (x - 1.0) for x in point]

        answer = self.session.call_value(
            "num_bfgs", [3.0, -2.0], 1e-10, 100, 1000, objective, gradient
        )
        point = [float(value) for value in answer[0]]
        self.assertAlmostEqual(point[0], 1.0, places=6)
        self.assertAlmostEqual(point[1], 1.0, places=6)
        self.assertLess(float(answer[1]), 1e-8)
        self.assertTrue(seen)
        self.assertIn(seen[0], {"ndarray", "list"})

    def test_a_callable_is_refused_by_the_term_abi(self) -> None:
        with self.assertRaises(rssn.RssnError) as raised:
            self.session.apply("num_brent", 1.0, 2.0, 1e-12, 64, self.residual)
        self.assertEqual(raised.exception.status, rssn.Status.INVALID_ARGUMENT)


class DataTests(unittest.TestCase):
    def test_units(self) -> None:
        self.assertEqual(len(rssn.units.UNITS), 82)
        metre = rssn.units.get("m")
        self.assertEqual(metre.dimension, "Length")
        self.assertEqual(metre.scale, Fraction(1))
        self.assertEqual(rssn.units.get("metre").name, "m")
        foot = rssn.units.get("ft")
        self.assertEqual(foot.scale, Fraction(381, 1250))
        self.assertFalse(foot.prefixable)
        self.assertTrue(rssn.units.get("m").prefixable)

    def test_constants(self) -> None:
        light = rssn.constants.get("c")
        self.assertEqual(light.symbol, "c")
        self.assertEqual(int(light.numer), 299792458)
        self.assertTrue(light.exact)
        self.assertGreater(rssn.constants.get("G").uncertainty_value, 0)
        self.assertAlmostEqual(light.value(), 299792458.0)

    def test_scenarios(self) -> None:
        self.assertEqual(len(rssn.scenarios.SCENARIOS), 29)
        self.assertIn("ising", rssn.sim.names())
        self.assertIn("dim", rssn.sim.params("poisson_multigrid"))
        self.assertEqual(rssn.sim.defaults("poisson_multigrid")["dim"], 2)


class SimulationTests(unittest.TestCase):
    def test_ising(self) -> None:
        output = rssn.sim.run(
            "ising", {"width": 4, "height": 4, "temperature": 2.0, "mc_steps": 5, "seed": 7}
        )
        self.assertIn("spins", output)
        self.assertIn("magnetization", output)

    def test_defaults(self) -> None:
        output = rssn.sim.run("poisson_multigrid", {"dim": 2, "n": 9, "tol": 1e-8})
        self.assertIsInstance(output, dict)
        json.dumps(output)  # must be JSON-serialisable

    def test_unknown_scenario(self) -> None:
        with self.assertRaises(rssn.RssnError):
            rssn.sim.run("no_such_scenario", {})


class StructTests(unittest.TestCase):
    def test_answer_layout(self) -> None:
        # The ctypes mirror must agree with the C struct in rssn.h.
        answer = rssn._ffi.RssnAnswer()
        self.assertEqual(ctypes.sizeof(answer), 40)
        self.assertEqual(ctypes.sizeof(rssn._ffi.RssnTensor), 32)


if __name__ == "__main__":
    unittest.main()
