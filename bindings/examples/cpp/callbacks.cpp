// Higher-order operators: C++ functions the engine calls.
//
//     cmake -S bindings/cpp -B build/cpp -DRSSN_ROOT=$PWD
//     cmake --build build/cpp --target rssn_callbacks
//     ./build/cpp/rssn_callbacks
//
// Some operators take a function and call it — root finders, quadrature,
// optimisers, bootstrap. A `Callback` wraps a `std::function` in the C function
// pointer the engine calls.

#include <cmath>
#include <cstdio>
#include <vector>

#include "rssn/callback.hpp"
#include "rssn/client.hpp"
#include "rssn/ops.hpp"

namespace {

using rssn::client::Callback;
using rssn::client::SessionHandle;
using rssn::client::Value;
using rssn::client::value_number;

/// `x -> x*x - 2`, the residual whose positive root is `sqrt(2)`.
Value residual(SessionHandle session, const std::vector<Value> &arguments) {
  const double x = arguments[0].as_f64();
  return value_number(session, x * x - 2.0);
}

void root_of_two() {
  // The generated facade goes through the thread's default session, so the
  // callback has to be registered on that one.
  rssn::client::Session &session = rssn::client::default_session();
  auto function =
      Callback::for_argument(session.raw(), "num_brent", 4, residual);
  const Value answer = rssn::client::num_brent(1.0, 2.0, 1e-12, 64, function);
  std::printf("num_brent:        sqrt(2) = %.15g\n", answer.at(0).as_f64());
  std::printf("                  evaluations = %.0f\n", answer.at(1).as_f64());
}

void explicit_callback() {
  rssn::client::Session session;
  auto function = Callback::create(session.raw(), residual,
                                   R"({"params": ["float"], "ret": "float"})");
  const Value answer =
      session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, function);
  std::printf("explicit Callback: %.15g\n", answer.at(0).as_f64());
}

void minimise() {
  rssn::client::Session session;
  auto objective = Callback::for_argument(
      session.raw(), "num_bfgs", 4,
      [](SessionHandle handle, const std::vector<Value> &arguments) {
        double total = 0.0;
        const Value &point = arguments[0];
        for (std::size_t index = 0; index < point.size(); ++index) {
          const double difference = point.at(index).as_f64() - 1.0;
          total += difference * difference;
        }
        return value_number(handle, total);
      });
  auto gradient = Callback::for_argument(
      session.raw(), "num_bfgs", 5,
      [](SessionHandle handle, const std::vector<Value> &arguments) {
        std::vector<double> out;
        const Value &point = arguments[0];
        for (std::size_t index = 0; index < point.size(); ++index) {
          out.push_back(2.0 * (point.at(index).as_f64() - 1.0));
        }
        return rssn::client::value_tensor(handle, out, {out.size()});
      });
  const Value answer =
      session.call_value("num_bfgs", std::vector<double>{3.0, -2.0}, 1e-10, 100,
                         1000, objective, gradient);
  const std::vector<double> point = answer.at(0).tensor();
  std::printf("num_bfgs:         minimum at (%.6g, %.6g), value %.6g\n",
              point[0], point[1], answer.at(1).as_f64());
}

void failures_cross() {
  rssn::client::Session session;
  auto function = Callback::for_argument(
      session.raw(), "num_brent", 4,
      [](SessionHandle, const std::vector<Value> &) -> Value {
        throw std::runtime_error("no thanks");
      });
  try {
    session.call_value("num_brent", 1.0, 2.0, 1e-12, 64, function);
  } catch (const rssn::client::Error &error) {
    std::printf("callback failure: %s\n", error.what());
  }
}

} // namespace

int main() {
  root_of_two();
  explicit_callback();
  minimise();
  failures_cross();
  return 0;
}
