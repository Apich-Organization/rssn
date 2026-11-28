// The workflows `bindings/cpp/README.md` shows, as one program that runs them.
//
//     cmake --build build/cpp --target rssn_workflow &&
//     ./build/cpp/rssn_workflow
//
// A README snippet is a fragment of a function, not a translation unit, so the
// other bindings' examples are extracted and executed by
// `bindings/tests/run_readme_examples.py` while the C++ ones are compiled here.
// Every block of the README appears below in the same order, with a check that
// it produced what the comment says.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "rssn/callback.hpp"
#include "rssn/client.hpp"
#include "rssn/constants.hpp"
#include "rssn/ops.hpp"
#include "rssn/scenarios.hpp"
#include "rssn/units.hpp"

namespace {

int failures = 0;

void expect(bool condition, const std::string &what) {
  if (!condition) {
    std::printf("FAIL %s\n", what.c_str());
    ++failures;
  } else {
    std::printf("ok   %s\n", what.c_str());
  }
}

} // namespace

int main() {
  using namespace rssn::client;

  // -- symbolic calculus -------------------------------------------------
  Session session;
  const auto x = session.sym("x");

  const auto request = diff(sin(x) * cos(x), x);
  expect(request.text() == "diff(cos(x)*sin(x), x)",
         "a request is not a result");
  expect(session.compute(request).text == "cos(2*x)", "the derivative");

  expect(session.compute(integral(x * x, x)).text == "1/3*x^3",
         "an indefinite integral");
  expect(session.compute("defint(x^2, x, 0, 1)").text == "1/3",
         "a definite integral");
  expect(session.compute("taylor(sin(x), x, 0, 5)").text ==
             "x - 1/6*x^3 + 1/120*x^5",
         "a Taylor series");
  expect(session.compute("limit(sin(x)/x, x, 0)").text == "1", "a limit");

  // A `Term` overloads the arithmetic operators.
  const Term polynomial = x * x + 2 * x + 1;
  expect(session.compute(polynomial).text == "x^2 + 2*x + 1",
         "arithmetic builds a request");
  expect(session.compute("diff(x^2, x)").text == "2*x",
         "text is parsed where it is computed");

  // -- a numeric target --------------------------------------------------
  Config cfg;
  cfg.numeric(1e-12).bind("x", 0.5);
  const auto answer = session.compute(diff(sin(x), x), &cfg);
  expect(answer.value.has_value() &&
             std::abs(*answer.value - 0.877582561890373) < 1e-12,
         "a numeric target");
  expect(answer.error.has_value(), "the estimated error");

  const std::vector<Answer> answers = session.compute_all(
      {diff(x * x * x, x), integral(x * x, x), session.parse("gcd(4, 6)")});
  expect(answers.size() == 3, "a batch");

  // -- linear algebra and tensors ----------------------------------------
  const auto matrix = session.parse(
      "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))");
  std::vector<std::size_t> shape;
  const std::vector<double> entries =
      session.compute(matrix).term.tensor(&shape);
  expect(shape.size() == 2 && shape[0] == 2 && shape[1] == 2, "the extents");
  expect(entries.size() == 4 && entries[0] == 2.0 && entries[1] == 1.0,
         "row-major reading order");

  // -- a simulation ------------------------------------------------------
  const std::string simulation = Session::run_scenario(
      "ising",
      R"({"width":8,"height":8,"temperature":2.0,"mc_steps":10,"seed":3})");
  expect(simulation.find("magnetization") != std::string::npos, "an Ising run");
  const ScenarioDef *scenario = lookup_scenario("poisson_multigrid");
  expect(scenario != nullptr, "the scenario catalog");
  expect(std::string(scenario->params).find("dim:int=2") != std::string::npos,
         "the parameters a scenario reads");

  // -- units and constants ----------------------------------------------
  const UnitDef *foot = lookup_unit("ft");
  expect(foot != nullptr && std::string(foot->dimension) == "Length",
         "a unit's dimension");
  const ConstantDef *light = lookup_constant("c");
  expect(light != nullptr &&
             std::abs(constant_value(*light) - 299792458.0) < 1e-6,
         "the speed of light");

  const auto dimension = R"({"base":"Length"})";
  const auto differences = session.value(std::vector<double>{2.0, 3.0})
                               .with_coherent_dimension(dimension);
  const auto initial = session.value(std::vector<double>{0.0})
                           .with_coherent_dimension(dimension);
  const auto positions =
      session.call_value("num_timeseries_integrate", differences, initial);
  expect(positions.tensor() == std::vector<double>{0.0, 2.0, 5.0},
         "integrating coherent Length differences");
  expect(positions.at(0).json().find("Length") != std::string::npos,
         "the measured result retains its wire dimension");

  // -- checked dynamic calls ---------------------------------------------
  const auto builder = &num_special_ext_lerch;
  expect(builder(session.parse("-1234"), session.parse("2"), session.parse("1"))
                 .text()
                 .find("num_special_ext_lerch(") == 0,
         "the original dynamic request builder");
  const auto lerch_value =
      checked::num_special_ext_lerch(session, -1234.0, 2.0, 1.0);
  expect(std::abs(lerch_value.as_f64() - 0.021861589304464496) < 1e-15,
         "a checked native Lerch call");
  const auto polylog_value =
      checked::num_special_ext_polylog(session, 2.0, -1234.0);
  expect(std::abs(polylog_value.as_f64() / -1234.0 - lerch_value.as_f64()) <
             1e-15,
         "checked polylog agrees");
  const auto point =
      checked::num_uhp_to_disk(session, std::complex<double>(0.0, 1.0));
  expect(point.kind() == "complex" &&
             point.json().find("\"im\":0.0") != std::string::npos,
         "a checked native complex call");
  const auto recovered = checked::num_disk_to_uhp(session, point);
  expect(recovered.json().find("\"im\":1.0") != std::string::npos,
         "pass a checked value to another call");
  try {
    (void)checked::num_disk_to_uhp(session, std::complex<double>(1.0, 0.0));
    expect(false, "checked facade preserves domain errors");
  } catch (const Error &error) {
    expect(error.status() == Status::InvalidArgument,
           "checked facade preserves domain errors");
  }

  // -- finding an operator ----------------------------------------------
  const OpSpec *spec = op_spec("diff");
  expect(spec != nullptr && spec->min_arity == 2 && spec->max_arity == 2,
         "the operator table");
  expect(spec != nullptr && spec->ffi,
         "the operator is reachable through the C ABI");

  // -- dynamic values ----------------------------------------------------
  const Value backend = session.call_value("num_backend_current");
  expect(backend.kind() == "opaque", "an opaque handle");
  expect(session.call_value("num_backend_name", backend).text() ==
             session.current_backend(),
         "the active backend's name");
  session.release_backend(backend);

  // -- backend selection --------------------------------------------------
  const auto available = session.available_backends();
  expect(std::find(available.begin(), available.end(), "interpreter") !=
             available.end(),
         "the interpreter is available");
  const auto interpreter = session.open_backend("interpreter");
  expect(session.call_value("num_backend_name", interpreter).text() ==
             "interpreter",
         "an explicitly selected backend");
  session.release_backend(interpreter);
  const Session other;
  expect(session.set_default_backend("interpreter") == "interpreter" &&
             other.current_backend() == "interpreter",
         "sessions share the process default");
  try {
    (void)other.set_default_backend("not-an-rssn-backend");
    expect(false, "unknown backends report a typed error");
  } catch (const Error &error) {
    expect(error.status() == Status::InvalidArgument &&
               other.current_backend() == "interpreter",
           "a failed selection preserves the default");
  }
  const auto automatic = session.set_default_backend("auto");
  expect(other.current_backend() == automatic, "restore automatic selection");

  const Value zero = session.value(0.0);
  auto residual = Callback::for_argument(
      session.raw(), "num_brent", 4,
      [&](SessionHandle, const std::vector<Value> &) {
        return Value::borrowed(zero.session(), zero.cell());
      });
  const Value endpoint = num_brent(1.0, 2.0, 1e-12, 64, residual);
  expect(endpoint.session() == session.raw() && zero.as_f64() == 0.0,
         "the facade infers the callback session and preserves a retained "
         "return value");
  expect(session.value(true).kind() == "bool" &&
             session.value(std::complex<double>(1.0, 2.0)).kind() ==
                 "complex" &&
             session.value("abc").text() == "abc",
         "host scalar value kinds");
  const Value copy = zero.clone();
  expect(copy.cell() != zero.cell() && copy.as_f64() == 0.0,
         "independent value cells");

  const std::string digits = "-1267650600228229401496703205383";
  const auto exact = value_integer_decimal(session.raw(), digits);
  expect(exact.kind() == "int" && exact.text() == digits,
         "decimal integer cell");

  // -- NumPy byte payloads ------------------------------------------------
  const std::string manifest = Session::manifest_json();
  constexpr std::string_view npy_flag = "\"npy\":";
  const auto flag = manifest.find(npy_flag);
  const auto enabled =
      flag == std::string::npos
          ? flag
          : manifest.find_first_not_of(" \t\r\n", flag + npy_flag.size());
  if (enabled != std::string::npos &&
      manifest.compare(enabled, 4, "true") == 0) {
    const std::vector<double> samples{-0.0, 1.0, 2.0, 3.0, INFINITY, NAN};
    const auto payload = session.encode_npy({2, 3}, samples);
    expect(payload.size() > 6 && payload[0] == 0x93 && payload[1] == 'N',
           "the NumPy file signature");
    std::vector<std::size_t> extents;
    const auto recovered = session.decode_npy(payload).tensor(&extents);
    expect(extents == std::vector<std::size_t>{2, 3} &&
               recovered.size() == samples.size() &&
               std::signbit(recovered[0]) && recovered[3] == 3.0 &&
               std::isinf(recovered[4]) && std::isnan(recovered[5]),
           "NumPy shape, row order and IEEE samples");
    const auto empty =
        session.decode_npy(session.encode_npy({2, 0, 3}, {})).tensor(&extents);
    expect(empty.empty() && extents == std::vector<std::size_t>{2, 0, 3},
           "empty NumPy extents");
    const auto scalar =
        session.decode_npy(session.encode_npy({}, {-0.0})).tensor(&extents);
    expect(extents.empty() && scalar.size() == 1 && std::signbit(scalar[0]),
           "a NumPy scalar");
    try {
      (void)session.decode_npy({0x93, 'N'});
      expect(false, "truncated NumPy bytes fail");
    } catch (const Error &error) {
      expect(error.status() == Status::InvalidArgument,
             "truncated NumPy bytes fail");
    }
  }

  std::printf("%s\n", failures == 0 ? "every README example ran" : "FAILURES");
  return failures == 0 ? 0 : 1;
}
