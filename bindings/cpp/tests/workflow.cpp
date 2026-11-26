// The workflows `bindings/cpp/README.md` shows, as one program that runs them.
//
//     cmake --build build/cpp --target rssn_workflow && ./build/cpp/rssn_workflow
//
// A README snippet is a fragment of a function, not a translation unit, so the
// other bindings' examples are extracted and executed by
// `bindings/tests/run_readme_examples.py` while the C++ ones are compiled here.
// Every block of the README appears below in the same order, with a check that
// it produced what the comment says.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "rssn/client.hpp"
#include "rssn/constants.hpp"
#include "rssn/ops.hpp"
#include "rssn/scenarios.hpp"
#include "rssn/units.hpp"

namespace {

int failures = 0;

void expect(bool condition, const std::string& what) {
    if (!condition) {
        std::printf("FAIL %s\n", what.c_str());
        ++failures;
    } else {
        std::printf("ok   %s\n", what.c_str());
    }
}

}  // namespace

int main() {
    using namespace rssn::client;

    // -- symbolic calculus -------------------------------------------------
    Session session;
    const auto x = session.sym("x");

    const auto request = diff(sin(x) * cos(x), x);
    expect(request.text() == "diff(cos(x)*sin(x), x)", "a request is not a result");
    expect(session.compute(request).text == "cos(2*x)", "the derivative");

    expect(session.compute(integral(x * x, x)).text == "1/3*x^3", "an indefinite integral");
    expect(session.compute("defint(x^2, x, 0, 1)").text == "1/3", "a definite integral");
    expect(session.compute("taylor(sin(x), x, 0, 5)").text == "x - 1/6*x^3 + 1/120*x^5",
           "a Taylor series");
    expect(session.compute("limit(sin(x)/x, x, 0)").text == "1", "a limit");

    // A `Term` overloads the arithmetic operators.
    const Term polynomial = x * x + 2 * x + 1;
    expect(session.compute(polynomial).text == "x^2 + 2*x + 1", "arithmetic builds a request");
    expect(session.compute("diff(x^2, x)").text == "2*x", "text is parsed where it is computed");

    // -- a numeric target --------------------------------------------------
    Config cfg;
    cfg.numeric(1e-12).bind("x", 0.5);
    const auto answer = session.compute(diff(sin(x), x), &cfg);
    expect(answer.value.has_value() && std::abs(*answer.value - 0.877582561890373) < 1e-12,
           "a numeric target");
    expect(answer.error.has_value(), "the estimated error");

    const std::vector<Answer> answers =
        session.compute_all({diff(x * x * x, x), integral(x * x, x), gcd(4, 6)});
    expect(answers.size() == 3, "a batch");

    // -- linear algebra and tensors ----------------------------------------
    const auto matrix = session.parse(
        "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))");
    std::vector<std::size_t> shape;
    const std::vector<double> entries = session.compute(matrix).term.tensor(&shape);
    expect(shape.size() == 2 && shape[0] == 2 && shape[1] == 2, "the extents");
    expect(entries.size() == 4 && entries[0] == 2.0 && entries[1] == 1.0,
           "row-major reading order");

    // -- a simulation ------------------------------------------------------
    const std::string simulation = Session::run_scenario(
        "ising", R"({"width":8,"height":8,"temperature":2.0,"mc_steps":10,"seed":3})");
    expect(simulation.find("magnetization") != std::string::npos, "an Ising run");
    const ScenarioDef* scenario = lookup_scenario("poisson_multigrid");
    expect(scenario != nullptr, "the scenario catalog");
    expect(std::string(scenario->params).find("dim:int=2") != std::string::npos,
           "the parameters a scenario reads");

    // -- units and constants ----------------------------------------------
    const UnitDef* foot = lookup_unit("ft");
    expect(foot != nullptr && std::string(foot->dimension) == "Length", "a unit's dimension");
    const ConstantDef* light = lookup_constant("c");
    expect(light != nullptr && std::abs(constant_value(*light) - 299792458.0) < 1e-6,
           "the speed of light");

    // -- finding an operator ----------------------------------------------
    const OpSpec* spec = op_spec("diff");
    expect(spec != nullptr && spec->min_arity == 2 && spec->max_arity == 2,
           "the operator table");
    expect(spec != nullptr && spec->ffi, "the operator is reachable through the C ABI");

    // -- dynamic values ----------------------------------------------------
    const Value backend = session.call_value("num_backend_current");
    expect(backend.kind() == "opaque", "an opaque handle");
    expect(session.call_value("num_backend_name", backend).text() == "interpreter",
           "the backend's name");

    std::printf("%s\n", failures == 0 ? "every README example ran" : "FAILURES");
    return failures == 0 ? 0 : 1;
}
