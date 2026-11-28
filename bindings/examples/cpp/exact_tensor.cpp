// Build an exact rational matrix from row-major value cells.
#include "rssn/value.hpp"

#include <iostream>
#include <vector>

int main() {
  using namespace rssn::client;
  Session session;
  Session foreign;
  std::vector<Value> entries;
  entries.push_back(
      value_integer_decimal(session.raw(), "18446744073709551617"));
  entries.push_back(session.call_value("divide", 2, 3));
  entries.push_back(value_integer(session.raw(), 3));
  entries.push_back(value_integer(session.raw(), 4));

  auto matrix = value_tensor_cells(session.raw(), entries, {2, 2});
  entries.clear();
  auto normalized = session.call_value("num_native_qmatrix_new", matrix, 1);
  auto determinant = session.call_value("num_native_qlinalg_det", normalized);
  if (determinant.text() != "73786976294838206466")
    return 1;
  std::cout << "Exact determinant: " << determinant.text() << '\n';

  entries.push_back(value_integer(session.raw(), 7));
  auto scalar = value_tensor_cells(session.raw(), entries, {});
  auto empty =
      value_tensor_cells(session.raw(), {}, {2, 0, 3}, &entries.front());
  if (scalar.at(0).text() != "7" || empty.size() != 0 ||
      empty.json().find("\"elem\":\"int\"") == std::string::npos)
    return 2;

  try {
    (void)value_tensor_cells(session.raw(), entries, {2, 2});
    return 3;
  } catch (const Error &error) {
    if (error.status() != Status::InvalidArgument)
      return 4;
  }

  entries.clear();
  entries.push_back(value_integer(foreign.raw(), 7));
  try {
    (void)value_tensor_cells(session.raw(), entries, {1});
    return 5;
  } catch (const Error &error) {
    if (error.status() != Status::InvalidArgument)
      return 6;
  }
  entries.clear();
  return 0;
}
