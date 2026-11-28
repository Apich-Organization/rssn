// A tour of the C++ client: build it against the rssn shared library.
//
//   cmake -S bindings/cpp -B build/cpp -DRSSN_ROOT=<repo>
//   cmake --build build/cpp
//   ./build/cpp/rssn_smoke
//
// Or by hand, with the header directory and the C library on the command line:
//
//   g++ -std=c++17 -I bindings/cpp/include -I . tests/smoke.cpp -L<dir> -lrssn
//   -o smoke

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "rssn/client.hpp"
#include "rssn/constants.hpp"
#include "rssn/ops.hpp"
#include "rssn/scenarios.hpp"
#include "rssn/units.hpp"
#include "rssn/value.hpp"

namespace {

int failures = 0;

void expect(bool condition, const char *what) {
  if (!condition) {
    std::printf("FAIL %s\n", what);
    ++failures;
  } else {
    std::printf("ok   %s\n", what);
  }
}

} // namespace

int main() {
  using namespace rssn::client;

  std::printf("rssn %s\n", Session::version());

  Session session;
  const Term x = session.sym("x");

  const Term mixed_units = session.parse("quantity(3, m) + quantity(4, s)");
  const Term length = session.parse("quantity(2, m)");
  const Answer mixed_units_single = session.compute(mixed_units);
  const auto unit_batch = session.compute_all({mixed_units, length});
  expect(!mixed_units_single.diagnostic_codes.empty(),
         "single answer reports the mismatch");
  expect(unit_batch[0].dimension == mixed_units_single.dimension &&
             unit_batch[0].diagnostic_codes ==
                 mixed_units_single.diagnostic_codes &&
             unit_batch[1].dimension == "Length" &&
             unit_batch[1].diagnostic_codes.empty(),
         "batch dimensions and diagnostics stay with their requests");

  // A symbolic request, reduced to a closed form.
  const Answer derivative = session.compute(diff(sin(x) * cos(x), x));
  std::printf("d/dx sin(x)cos(x) = %s\n", derivative.text.c_str());
  expect(derivative.text == "cos(2*x)", "derivative is cos(2*x)");
  expect(!derivative.latex().empty(), "latex is non-empty");

  // The same request, evaluated.
  Config numeric;
  numeric.numeric(1e-12).bind("x", 0.0);
  const Answer at_zero = session.compute(diff(sin(x) * cos(x), x), &numeric);
  expect(at_zero.value.has_value() && std::abs(*at_zero.value - 1.0) < 1e-12,
         "derivative at x = 0 is 1");

  // Text is parsed as infix; numbers become exact literals.
  const Answer factored = session.compute(session.parse("factorint(360)"));
  std::printf("factorint(360) = %s\n", factored.text.c_str());

  // A matrix comes back as a flat row-major array with its shape.
  const Term matrix = session.parse(
      "matmul(list(list(1, 2), list(3, 4)), list(list(0, 1), list(1, 0)))");
  std::vector<std::size_t> shape;
  const std::vector<double> entries =
      session.compute(matrix).term.tensor(&shape);
  std::printf("matrix %zux%zu:", shape.empty() ? 0 : shape[0],
              shape.size() < 2 ? 0 : shape[1]);
  for (double entry : entries) {
    std::printf(" %g", entry);
  }
  std::printf("\n");
  expect(shape.size() == 2 && shape[0] == 2 && shape[1] == 2, "matrix is 2x2");

  // Arithmetic operators build requests.
  const Term polynomial = x * x + x * 2.0 + 1.0;
  std::printf("x^2 + 2x + 1 = %s\n", session.compute(polynomial).text.c_str());

  // The sugar: subtraction, division and negation are `add`/`mul`/`pow`
  // built the way the engine builds them.
  expect((x - 1).text() == "x - 1", "x - 1 is a subtraction");
  expect((1 - x).text() == "1 - x", "1 - x is a subtraction");
  expect((-x).text() == "-x", "negation is a product with -1");
  expect((x / 2).text() == "x/2", "x / 2 is a quotient");
  expect((2 / x).text() == "2/x", "2 / x is a quotient");
  expect(session.compute((x * x - 1) / (x - 1)).text == "x + 1",
         "ratio simplifies");

  // Batch computation.
  const std::vector<Term> requests = {
      session.parse("diff(x^3, x)"),
      session.parse("integral(x^2, x)"),
      session.parse("solve(x^2 - 4, x)"),
  };
  const std::vector<Answer> answers = session.compute_all(requests);
  for (const Answer &answer : answers) {
    std::printf("  batch: %s\n", answer.text.c_str());
  }
  expect(answers.size() == 3, "three batch answers");
  Session other;
  bool foreign_batch_rejected = false;
  try {
    (void)session.compute_all({x, other.sym("foreign")});
  } catch (const Error &error) {
    foreign_batch_rejected = error.status() == Status::InvalidArgument;
  }
  expect(foreign_batch_rejected,
         "a batch rejects term handles from other sessions");

  // Failures carry a status code.
  bool raised = false;
  try {
    (void)session.parse("sin(");
  } catch (const Error &error) {
    raised = true;
    std::printf("parse error: %s: %s\n", status_name(error.status()),
                error.what());
    expect(error.status() == Status::Parse, "parse error has the Parse status");
  }
  expect(raised, "a bad parse raises");

  // The value ABI: operators whose arguments are native objects.
  {
    const Value backend = session.call_value("num_backend_current");
    expect(backend.kind() == "opaque", "a backend handle is an opaque value");
    const Value name = session.call_value("num_backend_name", backend);
    expect(name.text() == session.current_backend(), "the active backend name");
    bool raised = false;
    try {
      (void)backend.to_term();
    } catch (const Error &) {
      raised = true;
    }
    expect(raised, "a native handle has no term form");

    const Value tape = session.call_value("num_ad_tape_new", 16LL);
    const Value x0 = session.call_value("num_ad_tape_var", tape, 1.5);
    const Value product = session.call_value(
        "num_ad_var_mul", session.call_value("num_ad_var_sin", x0), x0);
    const Value number = session.call_value("num_ad_var_value", product);
    expect(std::abs(number.as_f64() - std::sin(1.5) * 1.5) < 1e-12,
           "the AD tape evaluates sin(x)*x");
  }

  // Read record fields by name and retain positional access.
  {
    const Value moments =
        session.call_value("num_timeseries_rolling_moments",
                           std::vector<double>{1.0, 3.0, 5.0}, 2, 0);
    expect(moments.field_name(0) == "mean" &&
               moments.field_name(1) == "variance",
           "rolling record exposes its field names");
    expect(moments.field("mean").tensor() == std::vector<double>{2.0, 4.0} &&
               moments.field("variance").tensor() == moments.at(1).tensor(),
           "named record access reads the requested tensor");

    bool missing = false;
    try {
      (void)moments.field("missing");
    } catch (const Error &error) {
      missing = error.status() == Status::InvalidArgument;
    }
    expect(missing, "named record access rejects an unknown field");
  }

  // The manifest is readable.
  expect(Session::manifest_hash() != 0, "manifest hash is non-zero");
  expect(!Session::manifest_json().empty(), "manifest JSON is non-empty");

  // Units and constants travel with the manifest.
  expect(lookup_unit("m") != nullptr, "metre is in the unit table");
  expect(std::string(lookup_unit("metre")->name) == "m", "metre resolves to m");
  expect(lookup_constant("c") != nullptr, "the speed of light is in the table");
  expect(std::abs(constant_value(*lookup_constant("c")) - 299792458.0) < 1e-6,
         "the speed of light is exact");

  // The operator table answers a lookup by name.
  const OpSpec *spec = op_spec("diff");
  expect(spec != nullptr && spec->min_arity == 2 && spec->max_arity == 2,
         "diff takes two arguments");
  expect(op_spec("definitely_not_an_operator") == nullptr,
         "an unknown name is not registered");

  // A simulation, and its catalog.
  expect(lookup_scenario("ising") != nullptr,
         "ising is in the scenario catalog");
  const std::string simulation = Session::run_scenario(
      "poisson_multigrid", R"({"dim": 2, "n": 9, "tol": 1e-8})");
  expect(simulation.find("converged") != std::string::npos,
         "the multigrid run converged");
  bool simulation_raised = false;
  try {
    (void)Session::run_scenario("no_such_scenario");
  } catch (const Error &error) {
    simulation_raised = error.status() == Status::Simulation;
  }
  expect(simulation_raised, "an unknown scenario is a Simulation error");

  std::printf("%s\n", failures == 0 ? "all checks passed" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
