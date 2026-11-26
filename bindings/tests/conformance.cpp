// The cross-language conformance runner for the C++ client.
//
//     rssn_conformance bindings/tests/cases.tsv
//
// It reads the same table as the Python and Julia runners and checks the same
// answers, so the three wrappers are held to one set of expectations.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "rssn/callback.hpp"
#include "rssn/client.hpp"
#include "rssn/constants.hpp"
#include "rssn/ops.hpp"
#include "rssn/units.hpp"

namespace {

struct Case {
    std::string name;
    std::string kind;
    std::string payload;
    std::string target;
    std::string tolerance;
    std::string bindings;
    std::string expect_text;
    std::string expect_value;
    std::string expect_status;
};

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream stream(line);
    while (std::getline(stream, field, '\t')) {
        if (!field.empty() && field.back() == '\r') {
            field.pop_back();
        }
        fields.push_back(field);
    }
    if (!line.empty() && line.back() == '\t') {
        fields.emplace_back();
    }
    return fields;
}

std::vector<Case> read_cases(const char* path) {
    std::ifstream input(path);
    if (!input) {
        std::fprintf(stderr, "cannot open %s\n", path);
        std::exit(2);
    }
    std::vector<Case> cases;
    bool seen_header = false;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        const std::vector<std::string> fields = split_tabs(line);
        if (!seen_header) {
            seen_header = true;
            continue;
        }
        Case entry;
        std::string* slots[] = {&entry.name,        &entry.kind,        &entry.payload,
                                &entry.target,      &entry.tolerance,   &entry.bindings,
                                &entry.expect_text, &entry.expect_value, &entry.expect_status};
        for (std::size_t index = 0; index < 9; ++index) {
            if (index < fields.size()) {
                *slots[index] = fields[index];
            }
        }
        cases.push_back(entry);
    }
    return cases;
}

std::string format_double(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

bool close_enough(double left, double right) {
    const double scale = std::max(1.0, std::abs(right));
    return std::abs(left - right) <= 1e-9 * scale;
}

/// Every comma-separated key of `keys` appears as a string in `json`.
bool has_keys(const std::string& json, const std::string& keys) {
    std::size_t start = 0;
    while (start <= keys.size()) {
        const std::size_t comma = keys.find(',', start);
        const std::string key =
            keys.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!key.empty() && json.find("\"" + key + "\"") == std::string::npos) {
            return false;
        }
        if (comma == std::string::npos) {
            return true;
        }
        start = comma + 1;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace rssn::client;
    const char* path = argc > 1 ? argv[1] : "bindings/tests/cases.tsv";
    const std::vector<Case> cases = read_cases(path);

    int failures = 0;
    int passed = 0;
    for (const Case& entry : cases) {
        const std::string wanted = entry.expect_status.empty() ? "ok" : entry.expect_status;
        std::string status = "ok";
        std::string text;
        std::string value;
        try {
            if (entry.kind == "manifest") {
                value = std::to_string(Session::manifest_hash());
            } else if (entry.kind == "units") {
                value = std::to_string(std::size(kUnits));
            } else if (entry.kind == "constants") {
                value = std::to_string(std::size(kConstants));
            } else if (entry.kind == "sim") {
                const std::size_t space = entry.payload.find(' ');
                const std::string name = entry.payload.substr(0, space);
                const std::string parameters =
                    space == std::string::npos ? "{}" : entry.payload.substr(space + 1);
                char* out = nullptr;
                const rssn::RssnStatus code = rssn::rssn_sim_run(name.c_str(), parameters.c_str(), &out);
                if (code != rssn::RSSN_STATUS_OK) {
                    status = status_name(static_cast<Status>(code));
                } else {
                    const std::string json = detail::take_string(out);
                    const bool complete =
                        !json.empty() && json.front() == '{' && has_keys(json, entry.expect_text);
                    text = complete ? entry.expect_text : json.substr(0, 40);
                }
            } else {
                Session session;
                Config config;
                bool has_config = false;
                if (entry.target == "numeric") {
                    config.numeric(entry.tolerance.empty() ? 1e-12 : std::stod(entry.tolerance));
                    has_config = true;
                }
                std::string remaining = entry.bindings;
                while (!remaining.empty()) {
                    const std::size_t semicolon = remaining.find(';');
                    const std::string pair = remaining.substr(0, semicolon);
                    if (!pair.empty()) {
                        const std::size_t equals = pair.find('=');
                        config.bind(pair.substr(0, equals), std::stod(pair.substr(equals + 1)));
                        has_config = true;
                    }
                    if (semicolon == std::string::npos) {
                        break;
                    }
                    remaining = remaining.substr(semicolon + 1);
                }
                if (entry.kind == "apply") {
                    std::vector<Term> arguments;
                    std::string remaining_arguments = entry.bindings;
                    while (!remaining_arguments.empty()) {
                        const std::size_t comma = remaining_arguments.find(',');
                        const std::string literal =
                            remaining_arguments.substr(0, comma);
                        if (!literal.empty()) {
                            arguments.push_back(Term(std::stod(literal)));
                        }
                        if (comma == std::string::npos) {
                            break;
                        }
                        remaining_arguments = remaining_arguments.substr(comma + 1);
                    }
                    const Term built =
                        session.apply(entry.payload.c_str(), arguments.data(), arguments.size());
                    text = built.text();
                    if (text != entry.expect_text) {
                        std::printf("FAIL %-32s text %s != %s\n", entry.name.c_str(),
                                    text.c_str(), entry.expect_text.c_str());
                        ++failures;
                    } else {
                        ++passed;
                    }
                    continue;
                }
                if (entry.kind == "value" || entry.kind == "callback2") {
                    // The value ABI: an operator called with native objects.
                    // `value` passes the numbers in the payload; `callback2`
                    // starts from the point (3, -2) every runner uses and
                    // appends the objective and the gradient `num_bfgs`
                    // declares.
                    const std::size_t space = entry.payload.find(' ');
                    const std::string name = entry.payload.substr(0, space);
                    std::vector<Value> values;
                    if (entry.kind == "callback2") {
                        values.push_back(session.value(std::vector<double>{3.0, -2.0}));
                    }
                    std::string remaining = space == std::string::npos
                                                ? std::string()
                                                : entry.payload.substr(space + 1);
                    while (!remaining.empty()) {
                        const std::size_t comma = remaining.find(',');
                        const std::string literal = remaining.substr(0, comma);
                        if (!literal.empty()) {
                            const bool integral =
                                literal.find('.') == std::string::npos &&
                                literal.find('e') == std::string::npos &&
                                literal.find('E') == std::string::npos;
                            if (integral) {
                                values.push_back(
                                    session.value(static_cast<long long>(std::stoll(literal))));
                            } else {
                                values.push_back(session.value(std::stod(literal)));
                            }
                        }
                        if (comma == std::string::npos) {
                            break;
                        }
                        remaining = remaining.substr(comma + 1);
                    }
                    // The callbacks must outlive the call: a `Callback` retires
                    // its registration when it is destroyed, and the cell it
                    // handed over would no longer be live.
                    Callback objective;
                    Callback gradient;
                    if (entry.kind == "callback2") {
                        const std::size_t objective_index = values.size();
                        objective = Callback::for_argument(
                            session.raw(), name.c_str(), objective_index,
                            [](SessionHandle handle, const std::vector<Value>& arguments) {
                                double total = 0.0;
                                for (std::size_t index = 0; index < arguments[0].size(); ++index) {
                                    const double delta = arguments[0].at(index).as_f64() - 1.0;
                                    total += delta * delta;
                                }
                                return value_number(handle, total);
                            });
                        gradient = Callback::for_argument(
                            session.raw(), name.c_str(), objective_index + 1,
                            [](SessionHandle handle, const std::vector<Value>& arguments) {
                                std::vector<double> out;
                                for (std::size_t index = 0; index < arguments[0].size(); ++index) {
                                    out.push_back(2.0 * (arguments[0].at(index).as_f64() - 1.0));
                                }
                                return value_numbers(handle, out);
                            });
                        values.push_back(objective.value());
                        values.push_back(gradient.value());
                    }
                    const Value answer =
                        session.call_values(name.c_str(), values.data(), values.size());
                    std::ostringstream formatted;
                    if (entry.kind == "callback2") {
                        formatted << std::setprecision(10) << answer.at(0).at(0).as_f64();
                        text = formatted.str();
                    } else {
                        text = answer.text();
                    }
                    if (text != entry.expect_text) {
                        std::printf("FAIL %-32s text %s != %s\n", entry.name.c_str(),
                                    text.c_str(), entry.expect_text.c_str());
                        ++failures;
                    } else {
                        ++passed;
                    }
                    continue;
                }
                if (entry.kind == "callback") {
                    // A higher-order operator: `<op> <args>`, with the same
                    // residual `f(x) = x*x - 2` every runner uses. The first
                    // field of the answer record is the root.
                    const std::size_t space = entry.payload.find(' ');
                    const std::string name = entry.payload.substr(0, space);
                    std::vector<double> numbers;
                    std::vector<Value> values;
                    std::string remaining = space == std::string::npos
                                                ? std::string()
                                                : entry.payload.substr(space + 1);
                    while (!remaining.empty()) {
                        const std::size_t comma = remaining.find(',');
                        const std::string literal = remaining.substr(0, comma);
                        if (!literal.empty()) {
                            // An operator with an `Int` parameter refuses a
                            // `Float` argument, so keep integer literals
                            // integral.
                            if (literal.find('.') == std::string::npos &&
                                literal.find('e') == std::string::npos &&
                                literal.find('E') == std::string::npos) {
                                values.push_back(session.value(
                                    static_cast<long long>(std::stoll(literal))));
                            } else {
                                numbers.push_back(std::stod(literal));
                                values.push_back(session.value(numbers.back()));
                            }
                        }
                        if (comma == std::string::npos) {
                            break;
                        }
                        remaining = remaining.substr(comma + 1);
                    }
                    auto residual = Callback::for_argument(
                        session.raw(), name.c_str(), values.size(),
                        [](SessionHandle handle, const std::vector<Value>& arguments) {
                            const double x = arguments[0].as_f64();
                            return value_number(handle, x * x - 2.0);
                        });
                    values.push_back(residual.value());
                    const Value answer =
                        session.call_values(name.c_str(), values.data(), values.size());
                    std::ostringstream formatted;
                    formatted << std::setprecision(10) << answer.at(0).as_f64();
                    text = formatted.str();
                    if (text != entry.expect_text) {
                        std::printf("FAIL %-32s text %s != %s\n", entry.name.c_str(),
                                    text.c_str(), entry.expect_text.c_str());
                        ++failures;
                    } else {
                        ++passed;
                    }
                    continue;
                }
                const Answer answer =
                    session.compute(session.parse(entry.payload), has_config ? &config : nullptr);
                text = answer.text;
                if (answer.value && entry.target == "numeric") {
                    value = format_double(*answer.value);
                }
            }
        } catch (const Error& error) {
            status = status_name(error.status());
        }

        if (status != wanted) {
            std::printf("FAIL %-32s status %s != %s\n", entry.name.c_str(), status.c_str(),
                        wanted.c_str());
            ++failures;
            continue;
        }
        if (wanted == "ok") {
            if (text != entry.expect_text) {
                std::printf("FAIL %-32s text %s != %s\n", entry.name.c_str(), text.c_str(),
                            entry.expect_text.c_str());
                ++failures;
                continue;
            }
            if (!entry.expect_value.empty() && !value.empty()) {
                if (!close_enough(std::stod(value), std::stod(entry.expect_value))) {
                    std::printf("FAIL %-32s value %s != %s\n", entry.name.c_str(), value.c_str(),
                                entry.expect_value.c_str());
                    ++failures;
                    continue;
                }
            }
        }
        ++passed;
    }

    std::printf("%d/%zu cases pass\n", passed, cases.size());
    return failures == 0 ? 0 : 1;
}
