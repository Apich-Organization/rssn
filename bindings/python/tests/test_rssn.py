"""Unit tests for the Python binding.

    python -m unittest discover -s bindings/python/tests -v

They need a built rssn library (see `bindings/README.md`); the conformance table
in `bindings/tests/cases.tsv` is the cross-language half of the suite.
"""

from __future__ import annotations

import ctypes
import json
import math
import pathlib
import sys
import threading
import unittest
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import rssn  # noqa: E402


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
            self.assertEqual(entries, [[0.0, 1.0], [3.0, 2.0]])

    def test_compute_all(self) -> None:
        answers = self.session.compute_all(["diff(x^3, x)", "integral(x^2, x)", "gcd(4, 6)"])
        self.assertEqual(len(answers), 3)
        self.assertTrue(all(isinstance(answer, rssn.Answer) for answer in answers))

    def test_compute_all_collects_failures(self) -> None:
        # A numeric target on a symbolic answer is a per-request failure.
        config = rssn.Config().numeric(1e-12)
        answers = self.session.compute_all(["1 + 1", "diff(x^2, x)"], config, errors="collect")
        self.assertIsInstance(answers[0], rssn.Answer)
        self.assertIsInstance(answers[1], rssn.RssnError)
        self.assertEqual(answers[1].status, rssn.Status.NOT_NUMERIC)

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


class ValueTests(unittest.TestCase):
    """The value ABI: operators whose arguments have no term form."""

    def setUp(self) -> None:
        self.session = rssn.Session()

    def tearDown(self) -> None:
        self.session.close()

    def test_backend_handle(self) -> None:
        backend = self.session.call_value("num_backend_current")
        self.assertEqual(backend.kind(), "opaque")
        name = self.session.call_value("num_backend_name", backend)
        self.assertEqual(name.kind(), "string")
        self.assertEqual(name.text(), "interpreter")
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
