#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "rssn/client.hpp"
#include "rssn/value.hpp"

namespace {

using rssn::client::Error;
using rssn::client::Session;
using rssn::client::Status;
using rssn::client::Value;
using rssn::client::value_tensor_cells;

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool contains_tag(const Value& values, const char* tag) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (values.at(i).text() == tag) {
            return true;
        }
    }
    return false;
}

std::string name(Session& session, const Value& backend) {
    return session.call_value("num_backend_name", backend).text();
}

Value compile(Session& session, const Value* backend) {
    const auto metre = session.value(1.0).with_coherent_dimension(R"({"base":"Length"})");
    const auto reference = session.call_value("num_backend_reference", metre);
    std::vector<Value> names;
    names.push_back(session.value("x"));
    std::vector<Value> references;
    references.push_back(reference.clone());
    const auto inputs = value_tensor_cells(session.raw(), names, {names.size()});
    const auto scales = value_tensor_cells(session.raw(), references, {references.size()});
    if (backend != nullptr) {
        return session.call_value("num_backend_compile", *backend, "x+x", inputs, scales, reference, 4096);
    }
    return session.call_value("num_backend_compile_current", "x+x", inputs, scales, reference, 4096);
}

void check_value(Session& session, const Value& program) {
    const auto measured = session.value(3.25).with_coherent_dimension(R"({"base":"Length"})");
    std::vector<Value> values;
    values.push_back(session.call_value("num_backend_value", measured));
    const auto inputs = value_tensor_cells(session.raw(), values, {values.size()});
    const auto result = session.call_value("num_backend_call", program, inputs);
    require(result.as_f64() == 6.5 && result.json().find(R"("base":"Length")") != std::string::npos,
            "measured compiled program must return 6.5 Length");
}

void test_backends() {
    Session session;
    const auto names = session.call_value("num_backend_available");
    require(contains_tag(names, "interpreter") && contains_tag(names, "cranelift"),
            "release library must advertise the interpreter and Cranelift");
    const auto current = session.call_value("num_backend_current_name").text();
    const auto automatic = session.call_value("num_backend_named", "auto");
    require(current == name(session, automatic), "fresh sessions must use automatic backend selection");
    require(current == "tiered" || current == "interpreter", "unexpected automatic backend");
    check_value(session, compile(session, nullptr));

    const auto interpreter = session.call_value("num_backend_named", "interpreter");
    require(name(session, interpreter) == "interpreter", "explicit interpreter selection failed");
    check_value(session, compile(session, &interpreter));

    const auto native = session.call_value("num_backend_named", "cranelift");
    require(name(session, native) == "cranelift", "explicit Cranelift selection failed");
    Value native_program;
    try {
        native_program = compile(session, &native);
    } catch (const Error& error) {
        // Respect executable-memory policy. Other JIT errors are regressions.
        const std::string message = error.what();
        if (error.status() != Status::InvalidArgument ||
            message.find("code generation failed:") == std::string::npos ||
            message.find("unable to make memory readable+executable") == std::string::npos) {
            throw;
        }
        std::puts("Cranelift returned a typed executable-memory denial. Auto and interpreter passed.");
        return;
    }
    check_value(session, native_program);
    require(current == "tiered", "auto must select tiered when native JIT can run");
    std::puts("Auto, explicit interpreter and native Cranelift compiled and evaluated a measured program.");
}

void check_npy(const Value& record, const std::vector<double>& expected) {
    const auto shape = record.field("shape");
    require(shape.size() == 2 && shape.at(0).text() == "2" && shape.at(1).text() == "3",
            "NumPy record did not retain its exact 2x3 shape");
    const auto numbers = record.field("data").tensor();
    require(numbers == expected && std::signbit(numbers[0]), "NumPy data or the sign of zero changed");
}

void test_numpy() {
    Session session;
    // Written independently by numpy.save with dtype '<f8'.
    const std::string fixture =
        "934e554d5059010076007b276465736372273a20273c6638272c2027666f727472616e5f6f72646572273a204661"
        "6c73652c20277368617065273a2028322c2033292c207d2020202020202020202020202020202020202020202020"
        "20202020202020202020202020202020202020202020202020202020202020202020200a00000000000000800000"
        "00000000f43f00000000000004400000000000000e4000000000000012c00000000000001840";
    const std::vector<double> expected{-0.0, 1.25, 2.5, 3.75, -4.5, 6.0};
    const auto decoded = session.call_value("num_npy_decode", fixture);
    check_npy(decoded, expected);

    const auto encoded = session.call_value("num_npy_encode", decoded.field("shape"), session.value(expected)).text();
    require(encoded.rfind("934e554d5059", 0) == 0, "NumPy encoder did not emit an NPY magic prefix");
    // Two hex characters for each of the eight bytes of a binary64 cell.
    const auto payload = fixture.substr(fixture.size() - expected.size() * 16);
    require(encoded.size() >= payload.size() &&
            encoded.compare(encoded.size() - payload.size(), payload.size(), payload) == 0,
            "NumPy encoder must write the independent fixture's row-major binary64 bytes");
    check_npy(session.call_value("num_npy_decode", encoded), expected);
    std::puts("Read an independent NumPy fixture and round-tripped a 2x3 array.");
}

void print_manifest() {
    char* manifest = rssn::rssn_manifest_json();
    require(manifest != nullptr, "release manifest missing");
    const int status = std::fputs(manifest, stdout);
    rssn::rssn_string_free(manifest);
    require(status >= 0, "failed to write release manifest");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected backend, npy or manifest");
        const std::string mode = argv[1];
        if (mode == "backend") {
            test_backends();
        } else if (mode == "npy") {
            test_numpy();
        } else if (mode == "manifest") {
            print_manifest();
        } else {
            throw std::runtime_error("unknown feature test");
        }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
    return 0;
}
