// rssn C++ client: an RAII layer over the C ABI in `rssn.h`.
//
// The C interface is a handle-and-status API: `uint32_t` terms in one session,
// an `RssnStatus` code, and a thread-local message. This header turns that into
// values with destructors, exceptions carrying the code, and operators.
//
//     #include "rssn/client.hpp"
//     #include "rssn/ops.hpp"      // one function per registered operator
//
//     rssn::client::Session session;
//     auto x = session.sym("x");
//     auto answer = session.compute(rssn::client::diff(rssn::client::sin(x) * x, x));
//     std::cout << answer.text << "\n";
//
// Lifetime: a `Term` is a handle owned by its session, exactly as in C. It stays
// valid until the session is destroyed; destroying the session invalidates every
// term made through it. `Term` is copyable and cheap (a pointer and a `uint32_t`).
//
// Threads: a session is single-threaded (its store sits behind a `RefCell` in
// Rust). Use one session per thread; `default_session()` is thread-local.

#ifndef RSSN_CLIENT_HPP
#define RSSN_CLIENT_HPP

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "rssn.h"

namespace rssn {
namespace client {

/// The status codes of `rssn.h`, as a scoped enum.
enum class Status : int {
    Ok = RSSN_STATUS_OK,
    NullArgument = RSSN_STATUS_NULL_ARGUMENT,
    InvalidUtf8 = RSSN_STATUS_INVALID_UTF8,
    Parse = RSSN_STATUS_PARSE,
    UnknownOperator = RSSN_STATUS_UNKNOWN_OPERATOR,
    NotNumeric = RSSN_STATUS_NOT_NUMERIC,
    Inconsistent = RSSN_STATUS_INCONSISTENT,
    Rule = RSSN_STATUS_RULE,
    NoSuchTerm = RSSN_STATUS_NO_SUCH_TERM,
    NoValue = RSSN_STATUS_NO_VALUE,
    Simulation = RSSN_STATUS_SIMULATION,
    Panic = RSSN_STATUS_PANIC,
    Dimension = RSSN_STATUS_DIMENSION,
    InvalidArgument = RSSN_STATUS_INVALID_ARGUMENT,
    Cancelled = RSSN_STATUS_CANCELLED,
};

/// The name of a status code, for messages.
inline const char* status_name(Status status) {
    switch (status) {
        case Status::Ok: return "ok";
        case Status::NullArgument: return "null_argument";
        case Status::InvalidUtf8: return "invalid_utf8";
        case Status::Parse: return "parse";
        case Status::UnknownOperator: return "unknown_operator";
        case Status::NotNumeric: return "not_numeric";
        case Status::Inconsistent: return "inconsistent";
        case Status::Rule: return "rule";
        case Status::NoSuchTerm: return "no_such_term";
        case Status::NoValue: return "no_value";
        case Status::Simulation: return "simulation";
        case Status::Panic: return "panic";
        case Status::Dimension: return "dimension";
        case Status::InvalidArgument: return "invalid_argument";
        case Status::Cancelled: return "cancelled";
    }
    return "unknown";
}

class Session;
class Config;
class Term;
class Value;
struct Answer;
class Error;

/// The canonical name, rule set and signature digest from a manifest entry.
struct OpKey {
    std::string rule_set;
    std::string name;
    std::array<std::uint8_t, 16> signature_hash{};
};

namespace detail {

/// The message of the last failed call on this thread (may be empty).
inline std::string last_error() {
    const char* message = rssn_last_error();
    return message == nullptr ? std::string() : std::string(message);
}

/// Reads a caller-owned `char *` and releases it with `rssn_string_free`.
inline std::string take_string(char* pointer) {
    if (pointer == nullptr) {
        return std::string();
    }
    std::string text(pointer);
    rssn_string_free(pointer);
    return text;
}

}  // namespace detail

/// A failed rssn call: a status code plus the message rssn reported.
class Error : public std::runtime_error {
public:
    Error(Status status, std::string message)
        : std::runtime_error(message.empty() ? std::string(status_name(status)) : message),
          status_(status) {}

    /// The status code of the failure.
    Status status() const noexcept { return status_; }

private:
    Status status_;
};

namespace detail {

/// Throws `Error` unless `status` is `RSSN_STATUS_OK`.
inline void check(RssnStatus code) {
    if (code != RSSN_STATUS_OK) {
        throw Error(static_cast<Status>(code), last_error());
    }
}

/// Throws `Error` with `hint` and the thread's message.
[[noreturn]] inline void fail(Status hint, std::string message = std::string()) {
    if (message.empty()) {
        message = last_error();
    }
    throw Error(hint, std::move(message));
}

}  // namespace detail

/// The session that owns a term: the raw C pointer, so a `Term` stays usable
/// while the `Session` object is moved.
using SessionHandle = RssnSession*;

/// The calling thread's session; defined after `Session` below.
Session& default_session();

/// A node of a session's graph: a `uint32_t` handle and the session that owns it.
///
/// Operators are requests, not results: `diff(sin(x) * x, x)` builds the
/// unevaluated derivative and `Session::compute` reduces it.
class Term {
public:
    /// An invalid term (no session, `RSSN_NO_TERM`).
    Term() = default;

    /// A handle `session` minted. Prefer the session's constructors.
    Term(SessionHandle session, std::uint32_t handle) : session_(session), handle_(handle) {}

    /// A number in the default session. Implicit, so `sin(1.0)` reads naturally.
    template <class Number,
              std::enable_if_t<std::is_arithmetic_v<Number> &&
                                   !std::is_same_v<std::decay_t<Number>, bool>,
                               int> = 0>
    Term(Number value);  // NOLINT(google-explicit-constructor)

    /// Infix text parsed in the default session, so `diff("x^2", "x")` works.
    Term(const char* text);  // NOLINT(google-explicit-constructor)

    /// Infix text parsed in the default session.
    Term(const std::string& text) : Term(text.c_str()) {}  // NOLINT

    /// Whether this term can be used: it has a session and a handle.
    bool valid() const noexcept { return session_ != nullptr && handle_ != RSSN_NO_TERM; }

    /// The raw `uint32_t` handle.
    std::uint32_t handle() const noexcept { return handle_; }

    /// The session that owns this term.
    SessionHandle session() const noexcept { return session_; }

    /// The term as infix text, parseable by `Session::parse`.
    std::string text() const {
        return detail::take_string(rssn_term_to_string(session_, handle_));
    }

    /// The term as LaTeX math.
    std::string latex() const {
        return detail::take_string(rssn_term_to_latex(session_, handle_));
    }

    /// The term as Typst math.
    std::string typst() const {
        return detail::take_string(rssn_term_to_typst(session_, handle_));
    }

    /// The term as a multi-line Unicode drawing.
    std::string pretty() const {
        return detail::take_string(rssn_term_to_pretty(session_, handle_));
    }

    /// The value of a literal number; `std::nullopt` for a symbolic term.
    std::optional<double> value() const {
        double out = 0.0;
        const RssnStatus code = rssn_term_to_f64(session_, handle_, &out);
        if (code == RSSN_STATUS_OK) {
            return out;
        }
        if (code == RSSN_STATUS_NO_VALUE) {
            return std::nullopt;
        }
        detail::fail(static_cast<Status>(code));
    }

    /// A dense `f64` extraction of a list term, row-major; `shape` receives the
    /// extents (empty for a scalar).
    std::vector<double> tensor(std::vector<std::size_t>* shape = nullptr) const {
        RssnTensor tensor{};
        detail::check(rssn_term_tensor_data(session_, handle_, &tensor));
        std::vector<double> data(tensor.data, tensor.data + tensor.len);
        if (shape != nullptr) {
            shape->assign(tensor.shape, tensor.shape + tensor.rank);
        }
        rssn_tensor_free(&tensor);
        return data;
    }

    /// Reduce this term to the target of `config` (symbolic by default).
    Answer compute(const Config* config = nullptr) const;

    /// Evaluate numerically without running any rule.
    double eval(const std::vector<std::pair<std::string, double>>& bindings = {}) const {
        std::vector<const char*> names;
        std::vector<double> values;
        names.reserve(bindings.size());
        values.reserve(bindings.size());
        for (const auto& binding : bindings) {
            names.push_back(binding.first.c_str());
            values.push_back(binding.second);
        }
        double out = 0.0;
        detail::check(rssn_term_eval(session_, handle_, names.empty() ? nullptr : names.data(),
                                     values.empty() ? nullptr : values.data(), names.size(), &out));
        return out;
    }

    /// The numeric value, computing it when the term is not a literal.
    double to_double() const;

    /// The term as a LaTeX-free `std::string`, for logging.
    friend bool operator==(const Term& left, const Term& right) {
        return left.session_ == right.session_ && left.handle_ == right.handle_;
    }

    friend bool operator!=(const Term& left, const Term& right) { return !(left == right); }

private:
    SessionHandle session_ = nullptr;
    std::uint32_t handle_ = RSSN_NO_TERM;
};

/// The answer of a computation: a reduced term plus the engine's report.
struct Answer {
    /// The answer term.
    Term term;
    /// The answer as infix text.
    std::string text;
    /// The numeric value, when one is known.
    std::optional<double> value;
    /// The estimated absolute error of `value`.
    std::optional<double> error;
    /// Whether the request was fully reduced.
    bool reduced = false;
    /// Outer iterations the engine performed.
    std::size_t iterations = 0;
    /// The inferred dimension as text, when one was inferred.
    std::optional<std::string> dimension;
    /// Stable codes of the advisory dimensional findings (`"mismatch"`, ...).
    std::vector<std::string> diagnostic_codes;

    /// The answer term as LaTeX math.
    std::string latex() const { return term.latex(); }
    /// The answer term as Typst math.
    std::string typst() const { return term.typst(); }
    /// The answer term as a Unicode drawing.
    std::string pretty() const { return term.pretty(); }
    /// The answer as a dense row-major array, when it is a rectangular list.
    std::vector<double> tensor(std::vector<std::size_t>* shape = nullptr) const {
        return term.tensor(shape);
    }
};

/// An owned `RssnConfig`: target, bindings, assumptions and dimension mode.
class Config {
public:
    /// A configuration with every standard rule set and a symbolic target.
    Config() {
        pointer_ = rssn_config_new();
        if (pointer_ == nullptr) {
            detail::fail(Status::Rule);
        }
    }

    ~Config() { release(); }

    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;

    Config(Config&& other) noexcept : pointer_(other.pointer_) { other.pointer_ = nullptr; }

    Config& operator=(Config&& other) noexcept {
        if (this != &other) {
            release();
            pointer_ = other.pointer_;
            other.pointer_ = nullptr;
        }
        return *this;
    }

    /// The raw `RssnConfig *`.
    RssnConfig* raw() const noexcept { return pointer_; }

    /// Ask for a closed form (the default).
    Config& symbolic() {
        detail::check(rssn_config_symbolic(pointer_));
        return *this;
    }

    /// Ask for a numeric answer within `tolerance`.
    Config& numeric(double tolerance) {
        detail::check(rssn_config_numeric(pointer_, tolerance));
        return *this;
    }

    /// Bind the symbol `name` to `value` for numeric evaluation.
    Config& bind(const std::string& name, double value) {
        detail::check(rssn_config_bind(pointer_, name.c_str(), value));
        return *this;
    }

    /// Declare facts about `name`; `facts` is a bit set (1 real, 2 nonzero,
    /// 4 nonnegative, 8 nonpositive, 16 integer, so 7 is positive).
    Config& assume(const std::string& name, std::uint8_t facts) {
        detail::check(rssn_config_assume(pointer_, name.c_str(), facts));
        return *this;
    }

    /// Limit the search; zero keeps the current value.
    Config& budget(std::size_t max_iterations, std::size_t max_nodes) {
        detail::check(rssn_config_budget(pointer_, max_iterations, max_nodes));
        return *this;
    }

    /// Set the dimension mode: 0 off, 1 advisory (the default), 2 strict.
    Config& dimension_mode(std::uint8_t mode) {
        detail::check(rssn_config_dimension_mode(pointer_, mode));
        return *this;
    }

    /// Disable hard dimensional errors even in strict mode.
    Config& bypass_dimensions() {
        detail::check(rssn_config_bypass_dimensions(pointer_));
        return *this;
    }

    /// Declare the dimension of `symbol` as text (`"m/s"`, `"Length/Time"`, `"any"`).
    Config& dimension(const std::string& symbol, const std::string& text) {
        detail::check(rssn_config_dimension(pointer_, symbol.c_str(), text.c_str()));
        return *this;
    }

    /// Register a base dimension the engine does not have (`"Currency"`).
    Config& define_base_dimension(const std::string& name) {
        detail::check(rssn_config_define_base_dimension(pointer_, name.c_str()));
        return *this;
    }

private:
    void release() noexcept {
        if (pointer_ != nullptr) {
            rssn_config_free(pointer_);
            pointer_ = nullptr;
        }
    }

    RssnConfig* pointer_ = nullptr;
};

/// An owned `RssnSession`: the owner of every term made through it.
class Session {
public:
    /// A session with every standard rule set.
    Session() {
        pointer_ = rssn_session_new();
        if (pointer_ == nullptr) {
            detail::fail(Status::Rule);
        }
    }

    ~Session() { release(); }

    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;

    Session(Session&& other) noexcept : pointer_(other.pointer_) { other.pointer_ = nullptr; }

    Session& operator=(Session&& other) noexcept {
        if (this != &other) {
            release();
            pointer_ = other.pointer_;
            other.pointer_ = nullptr;
        }
        return *this;
    }

    /// The raw `RssnSession *`.
    RssnSession* raw() const noexcept { return pointer_; }

    /// Whether the session is alive.
    bool valid() const noexcept { return pointer_ != nullptr; }

    /// The symbol `name`.
    Term sym(std::string_view name) const {
        const std::string text(name);
        return wrap(rssn_term_sym(pointer_, text.c_str()), Status::UnknownOperator);
    }

    /// An exact integer.
    Term integer(long long value) const {
        return wrap(rssn_term_int(pointer_, value), Status::InvalidArgument);
    }

    /// A floating point number.
    Term real(double value) const {
        return wrap(rssn_term_float(pointer_, value), Status::InvalidArgument);
    }

    /// The exact fraction `numerator / denominator`.
    Term rational(long long numerator, long long denominator) const {
        return wrap(rssn_term_rational(pointer_, numerator, denominator), Status::InvalidArgument);
    }

    /// Parse infix text such as `"diff(sin(x)^2, x)"`.
    Term parse(std::string_view source) const {
        const std::string text(source);
        return wrap(rssn_term_parse(pointer_, text.c_str()), Status::Parse);
    }

    /// Apply the operator `op` to `args`, as an unevaluated request.
    ///
    /// An argument from another session (a literal built in the default
    /// session, say) is copied in through its text form, exactly as the free
    /// `apply` does.
    Term apply(const char* op, const Term* args, std::size_t count) const;

    /// Apply the operator `op` to a list of terms.
    Term apply(const char* op, std::initializer_list<Term> args) const {
        return apply(op, args.begin(), args.size());
    }

    /// Reduce `term` to the target of `config` (symbolic by default).
    Answer compute(const Term& term, const Config* config = nullptr) const;

    /// Reduce infix text, parsed in **this** session.
    ///
    /// A bare `Term("x^2")` is parsed in the thread's default session, so
    /// `session.compute("x^2")` would otherwise mix two sessions; this overload
    /// is the one a string argument reaches, and it parses where it computes.
    Answer compute(const char* text, const Config* config = nullptr) const {
        return compute(parse(text), config);
    }

    /// Reduce every term in `requests`; the first failure raises.
    std::vector<Answer> compute_all(const std::vector<Term>& requests,
                                    const Config* config = nullptr) const;

    /// Reduce every term in `requests`, keeping a failed request's status in
    /// place instead of raising.
    std::vector<std::pair<Status, Answer>> compute_all_each(const std::vector<Term>& requests,
                                                            const Config* config = nullptr) const;

    /// The compute record (`RssnComputeRecord`) as JSON text.
    std::string compute_json(const Term& term, const Config* config = nullptr) const {
        char* out = nullptr;
        detail::check(rssn_compute_json(pointer_, term.handle(), config == nullptr ? nullptr : config->raw(), &out));
        return detail::take_string(out);
    }

    /// The compute record in the deterministic `bincode-next` wire format.
    std::vector<std::uint8_t> compute_bincode(const Term& term,
                                              const Config* config = nullptr) const {
        std::uint8_t* data = nullptr;
        std::size_t length = 0;
        detail::check(rssn_compute_bincode(pointer_, term.handle(),
                                           config == nullptr ? nullptr : config->raw(), &data,
                                           &length));
        std::vector<std::uint8_t> bytes(data, data + length);
        rssn_bytes_free(data, length);
        return bytes;
    }

    /// The inferred dimension of `term` as text (`""` when none).
    std::string dimension(const Term& term, const Config* config = nullptr) const {
        char* out = nullptr;
        detail::check(rssn_dimension(pointer_, term.handle(),
                                     config == nullptr ? nullptr : config->raw(), &out));
        return detail::take_string(out);
    }

    /// Evaluate `term` numerically without running any rule.
    double eval(const Term& term,
                const std::vector<std::pair<std::string, double>>& bindings = {}) const {
        return term.eval(bindings);
    }

    /// Store `term` as a dynamic value (see `rssn/value.hpp`).
    Value value(const Term& term) const;
    /// Store a number as a dynamic value.
    Value value(double number) const;
    /// Store an integer as a dynamic value.
    Value value(long long number) const;
    /// Store an `int` as a dynamic value (an exact overload of the two above).
    Value value(int number) const;
    template <class Integer, std::enable_if_t<std::is_integral_v<Integer>, int> = 0>
    Value value(Integer number) const;
    /// Preserve the Boolean or complex scalar kind.
    Value value(bool flag) const;
    Value value(std::complex<double> number) const;
    /// Store text as a dynamic value.
    Value value(const std::string& text) const;
    Value value(const char* text) const;
    /// Store a dense vector as a rank-1 tensor value.
    Value value(const std::vector<double>& numbers) const;
    /// Call `op` through the value ABI (see `rssn/value.hpp`).
    Value call_value(const char* op) const;

    /// The selectors usable on this host with the loaded library.
    std::vector<std::string> available_backends() const;
    /// The effective backend without changing process selection.
    std::string current_backend() const;
    /// An explicit backend handle. Release it with release_backend.
    Value open_backend(const std::string& selector) const;
    /// Release the native handle and invalidate every copy of it.
    void release_backend(const Value& backend) const;
    /// Select a process-wide default. "auto" restores the native policy.
    std::string set_default_backend(const std::string& selector) const;

    /// Encode flat row-major binary64 samples as NumPy bytes.
    std::vector<std::uint8_t> encode_npy(const std::vector<std::size_t>& shape,
                                        const std::vector<double>& samples) const;
    /// Decode a NumPy Float32 or Float64 payload as an owned binary64 tensor.
    Value decode_npy(const std::vector<std::uint8_t>& payload) const;
    /// Call `op` with `count` value arguments.
    Value call_values(const char* op, const Value* args, std::size_t count) const;
    /// Call `op` with arguments of any kind `Session::value` accepts.
    template <class... Args>
    Value call_value(const char* op, Args&&... args) const;

    /// Call the exact overload identified by a key from the live manifest.
    template <class... Args>
    Value call_key(const OpKey& key, Args&&... args) const;

    /// The message of the last failed call on this thread.
    static std::string last_error() { return detail::last_error(); }

    /// The first 64 bits of the manifest hash.
    static std::uint64_t manifest_hash() { return rssn_manifest_hash(); }

    /// The standard kernel manifest as JSON text.
    static std::string manifest_json() { return detail::take_string(rssn_manifest_json()); }

    /// One live operator manifest entry, with overloads, keys and dimensions.
    static std::string op_info_json(const std::string& name) {
        if (name.empty() || name.find('\0') != std::string::npos) {
            detail::fail(Status::InvalidArgument, "provide one operator name");
        }
        char* json = nullptr;
        detail::check(rssn_op_info_json(name.c_str(), &json));
        return detail::take_string(json);
    }

    /// Read one overload key at its zero-based position in the live operator contract.
    static OpKey op_key(const std::string& name, std::size_t index = 0) {
        if (name.empty() || name.find('\0') != std::string::npos) {
            detail::fail(Status::InvalidArgument, "provide one operator name");
        }
        OpKey key;
        char* rule_set = nullptr;
        char* canonical_name = nullptr;
        detail::check(rssn_op_key_at(name.c_str(), index, &rule_set, &canonical_name,
                                    key.signature_hash.data()));
        const std::unique_ptr<char, decltype(&rssn_string_free)> rule_set_owner(rule_set, rssn_string_free);
        const std::unique_ptr<char, decltype(&rssn_string_free)> name_owner(canonical_name, rssn_string_free);
        key.rule_set = rule_set_owner.get();
        key.name = name_owner.get();
        return key;
    }

    /// The crate version.
    static const char* version() { return rssn_version(); }

    /// Run the simulation scenario `name` and answer its JSON result.
    ///
    /// `parameters` is a JSON object; `{}` selects the defaults of the
    /// scenarios that have them. The catalog is `rssn::client::kScenarios`
    /// (see `rssn/scenarios.hpp`), and a scenario that does not exist is an
    /// `Error` with `Status::Simulation`.
    static std::string run_scenario(const char* name, const char* parameters = "{}") {
        char* out = nullptr;
        detail::check(rssn_sim_run(name, parameters, &out));
        return detail::take_string(out);
    }

private:
    Term wrap(std::uint32_t handle, Status hint) const {
        if (handle == RSSN_NO_TERM) {
            detail::fail(hint);
        }
        return Term(pointer_, handle);
    }

    void release() noexcept {
        if (pointer_ != nullptr) {
            rssn_session_free(pointer_);
            pointer_ = nullptr;
        }
    }

    RssnSession* pointer_ = nullptr;
};

/// The calling thread's session, created on first use.
///
/// The implicit numeric and text constructors of `Term` and the zero-argument
/// operator functions need a session; this is the one they use.
inline Session& default_session() {
    thread_local Session session;
    return session;
}

// ---------------------------------------------------------------------------
// Term: constructors that need a session
// ---------------------------------------------------------------------------

template <class Number,
          std::enable_if_t<std::is_arithmetic_v<Number> &&
                               !std::is_same_v<std::decay_t<Number>, bool>,
                           int>>
Term::Term(Number value) : session_(default_session().raw()) {
    if constexpr (std::is_integral_v<Number>) {
        handle_ = rssn_term_int(session_, static_cast<long long>(value));
    } else {
        handle_ = rssn_term_float(session_, static_cast<double>(value));
    }
    if (handle_ == RSSN_NO_TERM) {
        detail::fail(Status::InvalidArgument);
    }
}

inline Term::Term(const char* text) : session_(default_session().raw()) {
    handle_ = rssn_term_parse(session_, text);
    if (handle_ == RSSN_NO_TERM) {
        detail::fail(Status::Parse);
    }
}

namespace detail {

/// Apply `op` in `session`.
///
/// A handle belongs to the session that minted it, so an argument from another
/// session is copied into this one through its text form — exact, because the
/// engine's text output round-trips. That is what makes `x * 2.0` work: the
/// literal was built in the default session and is re-created in `x`'s.
inline Term apply_in(SessionHandle session, const char* op, const Term* args, std::size_t count) {
    std::vector<std::uint32_t> handles(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (!args[index].valid()) {
            detail::fail(Status::NoSuchTerm, "argument is not a valid term");
        }
        if (args[index].session() == session) {
            handles[index] = args[index].handle();
            continue;
        }
        const std::string text = args[index].text();
        const std::uint32_t copied = rssn_term_parse(session, text.c_str());
        if (copied == RSSN_NO_TERM) {
            detail::fail(Status::NoSuchTerm,
                         "cannot copy a term from another session: " + text);
        }
        handles[index] = copied;
    }
    const std::uint32_t handle =
        rssn_term_apply(session, op, handles.empty() ? nullptr : handles.data(), count);
    if (handle == RSSN_NO_TERM) {
        detail::fail(Status::UnknownOperator);
    }
    return Term(session, handle);
}

}  // namespace detail

/// Apply `op` in the first valid term's session, or in the default session.
inline Term apply(const char* op, const Term* args, std::size_t count) {
    SessionHandle session = nullptr;
    for (std::size_t index = 0; index < count; ++index) {
        if (args[index].valid()) {
            session = args[index].session();
            break;
        }
    }
    if (session == nullptr) {
        session = default_session().raw();
    }
    return detail::apply_in(session, op, args, count);
}

/// Apply `op` to a list of terms.
inline Term apply(const char* op, std::initializer_list<Term> args) {
    return apply(op, args.begin(), args.size());
}

inline Term Session::apply(const char* op, const Term* args, std::size_t count) const {
    return detail::apply_in(pointer_, op, args, count);
}

// ---------------------------------------------------------------------------
// Answers
// ---------------------------------------------------------------------------

namespace detail {

class MetadataStrings {
public:
    MetadataStrings(RssnAnswerMetadata* rows, std::size_t count) : rows_(rows), count_(count) {}
    MetadataStrings(const MetadataStrings&) = delete;
    MetadataStrings& operator=(const MetadataStrings&) = delete;
    ~MetadataStrings() {
        for (std::size_t index = 0; index < count_; ++index) {
            rssn_string_free(rows_[index].dimension);
            rssn_string_free(rows_[index].diagnostic_codes);
            rssn_string_free(rows_[index].message);
        }
    }

private:
    RssnAnswerMetadata* rows_;
    std::size_t count_;
};

inline Answer answer_from(SessionHandle session, const RssnAnswer& raw,
                          const RssnAnswerMetadata& metadata);

}  // namespace detail

inline Answer Term::compute(const Config* config) const {
    RssnAnswer raw{};
    RssnStatus status = RSSN_STATUS_OK;
    RssnAnswerMetadata metadata{};
    const detail::MetadataStrings owned(&metadata, 1);
    detail::check(rssn_compute_all_with_metadata(session_, &handle_, 1,
                                                config == nullptr ? nullptr : config->raw(),
                                                &raw, &status, &metadata));
    if (status != RSSN_STATUS_OK) {
        detail::fail(static_cast<Status>(status), metadata.message == nullptr ? "" : metadata.message);
    }
    return detail::answer_from(session_, raw, metadata);
}

inline double Term::to_double() const {
    if (const std::optional<double> literal = value()) {
        return *literal;
    }
    Config config;
    config.numeric(1e-12);
    const Answer answer = compute(&config);
    if (!answer.value) {
        detail::fail(Status::NotNumeric, "term has no numeric value");
    }
    return *answer.value;
}

inline Answer Session::compute(const Term& term, const Config* config) const {
    return term.compute(config);
}

inline std::vector<Answer> Session::compute_all(const std::vector<Term>& requests,
                                                const Config* config) const {
    std::vector<Answer> answers;
    answers.reserve(requests.size());
    for (const auto& entry : compute_all_each(requests, config)) {
        if (entry.first != Status::Ok) {
            throw Error(entry.first, detail::last_error());
        }
        answers.push_back(entry.second);
    }
    return answers;
}

inline std::vector<std::pair<Status, Answer>> Session::compute_all_each(
    const std::vector<Term>& requests, const Config* config) const {
    std::vector<std::pair<Status, Answer>> out;
    const std::size_t count = requests.size();
    if (count == 0) {
        return out;
    }
    std::vector<std::uint32_t> handles(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (requests[index].session() != pointer_) {
            detail::fail(Status::InvalidArgument, "a batch request belongs to another session");
        }
        handles[index] = requests[index].handle();
    }
    std::vector<RssnAnswer> answers(count);
    std::vector<RssnStatus> statuses(count);
    std::vector<RssnAnswerMetadata> metadata(count);
    const detail::MetadataStrings owned(metadata.data(), count);
    detail::check(rssn_compute_all_with_metadata(pointer_, handles.data(), count,
                                                config == nullptr ? nullptr : config->raw(),
                                                answers.data(), statuses.data(), metadata.data()));
    out.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const Status status = static_cast<Status>(statuses[index]);
        if (status == Status::Ok) {
            out.emplace_back(status, detail::answer_from(pointer_, answers[index], metadata[index]));
        } else {
            out.emplace_back(status, Answer{});
        }
    }
    return out;
}

namespace detail {

inline Answer answer_from(SessionHandle session, const RssnAnswer& raw,
                          const RssnAnswerMetadata& metadata) {
    Answer answer;
    answer.term = Term(session, raw.term);
    answer.text = detail::take_string(rssn_term_to_string(session, raw.term));
    if (raw.has_value) {
        answer.value = raw.value;
        answer.error = raw.error;
    }
    answer.reduced = raw.reduced;
    answer.iterations = raw.iterations;
    if (metadata.dimension != nullptr) {
        answer.dimension = metadata.dimension;
    }
    if (metadata.diagnostic_codes != nullptr) {
        std::string_view codes(metadata.diagnostic_codes);
        while (!codes.empty()) {
            const auto delimiter = codes.find('\n');
            answer.diagnostic_codes.emplace_back(codes.substr(0, delimiter));
            if (delimiter == std::string_view::npos) break;
            codes.remove_prefix(delimiter + 1);
        }
    }
    return answer;
}

}  // namespace detail

/// The operators that read as arithmetic, as free functions.
// `sub`, `neg` and `divide` are the engine's named forms of subtraction,
// negation and quotient: it rewrites them to the canonical `add` / `mul` /
// `pow` shape and prints both spellings the same way. The arithmetic quotient
// is `divide` because `div` is the divergence of a vector field.
inline Term operator+(const Term& left, const Term& right) { return apply("add", {left, right}); }
inline Term operator*(const Term& left, const Term& right) { return apply("mul", {left, right}); }
inline Term operator-(const Term& left, const Term& right) { return apply("sub", {left, right}); }
inline Term operator/(const Term& left, const Term& right) { return apply("divide", {left, right}); }
inline Term operator-(const Term& value) { return apply("neg", {value}); }

// `pow(base, exponent)` is the generated operator in `rssn/ops.hpp`.

}  // namespace client
}  // namespace rssn

#endif  // RSSN_CLIENT_HPP
