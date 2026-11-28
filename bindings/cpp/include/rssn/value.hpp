// rssn C++ client: dynamic values.
//
// The term API carries only values that are nodes of the expression graph. An
// operator that returns a native object — an AD tape, a compiled backend, a
// solver handle — returns a `Value` instead, and an operator that takes one can
// only be called with a `Value`. Those are the operators `op_spec(name)->ffi`
// reports as unreachable through the term ABI; `Session::call_value` is how they
// are called.
//
//     rssn::client::Session session;
//     auto backend = session.call_value("num_backend_current");
//     std::cout << session.call_value("num_backend_name", backend).text() << "\n";
//
// A `Value` owns a cell in its session and frees it when destroyed. It is
// move-only: two copies of one cell would free it twice.

#ifndef RSSN_CLIENT_VALUE_HPP
#define RSSN_CLIENT_VALUE_HPP

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "rssn/client.hpp"

namespace rssn {
namespace client {

class Session;

/// A dynamic value owned by a session.
///
/// Construct one with `Session::value(...)` or receive it from
/// `Session::call_value`; the session keeps the value alive and releases the
/// rest with itself.
class Value {
public:
    /// A value that names no cell.
    Value() = default;

    /// A cell `session` minted.
    Value(SessionHandle session, std::uint32_t cell, bool owned = true)
        : session_(session), cell_(cell), owned_(owned) {}

    /// A view of a cell the caller does not own — what a host callback receives
    /// for its arguments. Destroying it does not free the cell.
    static Value borrowed(SessionHandle session, std::uint32_t cell) {
        return Value(session, cell, false);
    }

    ~Value() { release(); }

    Value(const Value&) = delete;
    Value& operator=(const Value&) = delete;

    Value(Value&& other) noexcept
        : session_(other.session_), cell_(other.cell_), owned_(other.owned_) {
        other.session_ = nullptr;
        other.cell_ = 0;
    }

    Value& operator=(Value&& other) noexcept {
        if (this != &other) {
            release();
            session_ = other.session_;
            cell_ = other.cell_;
            owned_ = other.owned_;
            other.session_ = nullptr;
            other.cell_ = 0;
        }
        return *this;
    }

    /// Whether this value names a live cell.
    bool valid() const noexcept { return session_ != nullptr && cell_ != 0; }

    /// Whether this wrapper owns the cell.
    bool owns_cell() const noexcept { return owned_ && valid(); }

    /// A separate cell which can be released or transferred independently.
    Value clone() const {
        // Both entry points are available in earlier shared libraries.
        std::uint32_t holder_cell = 0;
        detail::check(rssn_value_new_list(session_, &cell_, 1, &holder_cell));
        Value holder(session_, holder_cell);
        return holder.at(0);
    }

    /// The raw cell id.
    std::uint32_t cell() const noexcept { return cell_; }

    /// The session that owns the cell.
    SessionHandle session() const noexcept { return session_; }

    /// The value's kind: `"float"`, `"tensor"`, `"opaque"`, ...
    std::string kind() const;

    /// The numeric value.
    double as_f64() const;

    /// The text of a string or symbol, or any value's display form.
    std::string text() const;

    /// The value in its `serde` JSON form.
    std::string json() const;

    /// Attach a concrete dimension to a Float already in coherent units.
    /// Pass DimExpr JSON such as `{"base":"Length"}`. Returns a fresh cell.
    Value with_coherent_dimension(std::string_view dimension_json) const;

    /// The value as a term, when it has a graph form.
    Term to_term() const;

    /// A dense extraction of a tensor value; `shape` receives the extents.
    std::vector<double> tensor(std::vector<std::size_t>* shape = nullptr) const;

    /// The number of elements (a tensor's length, a record's field count).
    std::size_t size() const;

    /// Element `index` as a new value.
    Value at(std::size_t index) const;

    /// The name of record field `index`.
    std::string field_name(std::size_t index) const;

    /// Read one record field by name. Missing fields report the InvalidArgument status.
    Value field(const std::string& name) const;

    /// Release the cell now rather than at destruction.
    void close() { release(); }

    /// Hand the cell to the caller: the value stops owning it.
    std::uint32_t detach() noexcept {
        const std::uint32_t cell = cell_;
        cell_ = 0;
        return cell;
    }

private:
    void release() noexcept {
        if (valid() && owned_) {
            rssn_value_free(session_, cell_);
        }
        session_ = nullptr;
        cell_ = 0;
    }

    SessionHandle session_ = nullptr;
    std::uint32_t cell_ = 0;
    bool owned_ = true;
};

namespace detail {

/// Stores a cell of the kind `make` builds, mapping failure to an exception.
template <class Make>
inline Value value_of(SessionHandle session, Make make) {
    std::uint32_t cell = 0;
    detail::check(make(session, &cell));
    return Value(session, cell);
}

}  // namespace detail

/// A float value, for a host callback to answer with.
inline Value value_number(SessionHandle session, double number) {
    return detail::value_of(session, [number](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_float(s, number, out);
    });
}

/// An exact integer value.
inline Value value_integer(SessionHandle session, long long number) {
    return detail::value_of(session, [number](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_int(s, number, out);
    });
}

/// An exact integer from signed decimal digits.
inline Value value_integer_decimal(SessionHandle session, const std::string& digits) {
    if (digits.find('\0') != std::string::npos) {
        detail::fail(Status::InvalidArgument, "decimal integer contains an embedded NUL");
    }
    return detail::value_of(session, [&](SessionHandle owner, std::uint32_t* out) {
        return rssn_value_new_int_decimal(owner, digits.c_str(), out);
    });
}

/// A truth value. A `bool` has no numeric form an operator would accept in its
/// place, so this is the only way to answer a `Bool` callback parameter.
inline Value value_bool(SessionHandle session, bool flag) {
    return detail::value_of(session, [flag](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_bool(s, flag ? 1 : 0, out);
    });
}

/// A complex value.
inline Value value_complex(SessionHandle session, double re, double im) {
    return detail::value_of(session, [re, im](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_complex(s, re, im, out);
    });
}

/// A string value.
inline Value value_text(SessionHandle session, const std::string& text) {
    return detail::value_of(session, [&text](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_text(s, text.c_str(), out);
    });
}

/// A dense `f64` tensor, row-major.
inline Value value_tensor(SessionHandle session, const std::vector<double>& data,
                          const std::vector<std::size_t>& shape) {
    return detail::value_of(session, [&](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_tensor(s, data.empty() ? nullptr : data.data(), data.size(),
                                     shape.empty() ? nullptr : shape.data(), shape.size(), out);
    });
}

/// A row-major tensor made from existing cells, preserving exact values.
/// Rank zero takes one cell. An exemplar supplies the element type of an empty tensor.
/// Each cell must belong to the supplied session.
inline Value value_tensor_cells(SessionHandle session, const std::vector<Value>& values,
                                const std::vector<std::size_t>& shape,
                                const Value* exemplar = nullptr) {
    if (exemplar != nullptr && exemplar->session() != session) {
        detail::fail(Status::InvalidArgument, "tensor exemplar belongs to another session");
    }
    std::vector<std::uint32_t> cells;
    cells.reserve(values.size());
    for (const auto& value : values) {
        if (value.session() != session) {
            detail::fail(Status::InvalidArgument, "tensor cell belongs to another session");
        }
        cells.push_back(value.cell());
    }
    return detail::value_of(session, [&](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_tensor_cells(s, cells.data(), cells.size(),
                                           shape.data(), shape.size(),
                                           exemplar == nullptr ? 0 : exemplar->cell(), out);
    });
}

/// A rank-1 tensor of floats, which is what a callback returning a vector
/// answers (`num_bfgs`'s gradient, say).
inline Value value_numbers(SessionHandle session, const std::vector<double>& numbers) {
    return value_tensor(session, numbers, {numbers.size()});
}

/// A rank-1 vector of exact unsigned 64-bit integers.
inline Value value_integer_vector(SessionHandle session, const std::vector<std::uint64_t>& numbers) {
    std::vector<Value> values;
    std::vector<std::uint32_t> cells;
    values.reserve(numbers.size());
    cells.reserve(numbers.size());
    for (const auto number : numbers) {
        if (number <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
            values.push_back(value_integer(session, static_cast<std::int64_t>(number)));
        } else {
            values.push_back(value_integer_decimal(session, std::to_string(number)));
        }
        cells.push_back(values.back().cell());
    }
    return detail::value_of(session, [&](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_list(s, cells.empty() ? nullptr : cells.data(), cells.size(), out);
    });
}

/// A record of the given fields, in order.
inline Value value_record(SessionHandle session, const std::vector<std::string>& names,
                          const std::vector<Value>& fields) {
    if (names.size() != fields.size()) {
        detail::fail(Status::InvalidArgument, "a record needs one name per field");
    }
    std::vector<const char*> labels;
    std::vector<std::uint32_t> cells;
    labels.reserve(names.size());
    cells.reserve(fields.size());
    for (const auto& name : names) {
        labels.push_back(name.c_str());
    }
    for (const auto& field : fields) {
        if (!field.valid() || field.session() != session) {
            detail::fail(Status::NoValue, "record field is not a live value of this session");
        }
        cells.push_back(field.cell());
    }
    return detail::value_of(session, [&](SessionHandle s, std::uint32_t* out) {
        return rssn_value_new_record(s, labels.empty() ? nullptr : labels.data(),
                                     cells.empty() ? nullptr : cells.data(), cells.size(), out);
    });
}

/// The kind name of a value-kind code (`rssn_value_kind`).
inline const char* value_kind_name(std::uint8_t code) {
    static const char* names[] = {
        "int", "real", "float", "complex", "bool", "string", "symbol", "unit",
        "quantity", "tensor", "sparse", "record", "variant", "term", "function",
        "dyn", "opaque",
    };
    return code < std::size(names) ? names[code] : "unknown";
}

inline std::string Value::kind() const {
    std::uint8_t code = 0;
    detail::check(rssn_value_kind(session_, cell_, &code));
    return value_kind_name(code);
}

inline double Value::as_f64() const {
    double out = 0.0;
    detail::check(rssn_value_to_f64(session_, cell_, &out));
    return out;
}

inline std::string Value::text() const {
    char* out = nullptr;
    detail::check(rssn_value_to_text(session_, cell_, &out));
    return detail::take_string(out);
}

inline std::string Value::json() const {
    char* out = nullptr;
    detail::check(rssn_value_json(session_, cell_, &out));
    return detail::take_string(out);
}

inline Value Value::with_coherent_dimension(std::string_view dimension_json) const {
    std::uint32_t out = 0;
    const char* bytes = dimension_json.empty() ? "" : dimension_json.data();
    detail::check(rssn_value_with_coherent_dimension(session_, cell_,
        reinterpret_cast<const std::uint8_t*>(bytes), dimension_json.size(), &out));
    return Value(session_, out);
}

inline Term Value::to_term() const {
    std::uint32_t handle = 0;
    detail::check(rssn_value_to_term(session_, cell_, &handle));
    return Term(session_, handle);
}

inline std::vector<double> Value::tensor(std::vector<std::size_t>* shape) const {
    RssnTensor raw{};
    detail::check(rssn_value_to_tensor(session_, cell_, &raw));
    const std::unique_ptr<RssnTensor, decltype(&rssn_tensor_free)> release(&raw, rssn_tensor_free);
    std::vector<double> data;
    if (raw.len != 0) {
        data.assign(raw.data, raw.data + raw.len);
    }
    if (shape != nullptr) {
        if (raw.rank != 0) {
            shape->assign(raw.shape, raw.shape + raw.rank);
        } else {
            shape->clear();
        }
    }
    return data;
}

inline std::size_t Value::size() const {
    std::size_t out = 0;
    detail::check(rssn_value_len(session_, cell_, &out));
    return out;
}

inline Value Value::at(std::size_t index) const {
    std::uint32_t cell = 0;
    detail::check(rssn_value_at(session_, cell_, index, &cell));
    return Value(session_, cell);
}

inline std::string Value::field_name(std::size_t index) const {
    char* name = nullptr;
    detail::check(rssn_value_field_name(session_, cell_, index, &name));
    return detail::take_string(name);
}

inline Value Value::field(const std::string& name) const {
    if (kind() != "record") {
        detail::fail(Status::InvalidArgument, "named access requires a record value");
    }
    const std::size_t fields = size();
    for (std::size_t index = 0; index < fields; ++index) {
        if (field_name(index) == name) {
            return at(index);
        }
    }
    detail::fail(Status::InvalidArgument, "record field not found: " + name);
}

// ---------------------------------------------------------------------------
// Session: the value ABI
// ---------------------------------------------------------------------------

namespace detail {

inline Value as_value(SessionHandle session, const Term& term) {
    if (!term.valid() || term.session() != session) {
        fail(Status::InvalidArgument, "term is not a live term of this session");
    }
    return value_of(session, [&](SessionHandle owner, std::uint32_t* out) {
        return rssn_value_new_term(owner, term.handle(), out);
    });
}

inline Value as_value(SessionHandle session, double number) { return value_number(session, number); }
inline Value as_value(SessionHandle session, bool flag) { return value_bool(session, flag); }
inline Value as_value(SessionHandle session, std::complex<double> number) {
    return value_complex(session, number.real(), number.imag());
}
template <class Integer, std::enable_if_t<std::is_integral_v<Integer>, int> = 0>
inline Value as_value(SessionHandle session, Integer number) {
    if constexpr (std::numeric_limits<Integer>::digits > std::numeric_limits<std::int64_t>::digits) {
        if (number > static_cast<Integer>(std::numeric_limits<std::int64_t>::max())) {
            fail(Status::InvalidArgument, "integer exceeds the value ABI's signed 64-bit range");
        }
    }
    return value_integer(session, static_cast<long long>(number));
}
inline Value as_value(SessionHandle session, const std::string& text) { return value_text(session, text); }
inline Value as_value(SessionHandle session, const char* text) {
    if (text == nullptr) fail(Status::InvalidArgument, "null text argument");
    return value_text(session, text);
}
inline Value as_value(SessionHandle session, const std::vector<double>& numbers) {
    return value_numbers(session, numbers);
}

}  // namespace detail

inline Value Session::value(const Term& term) const { return detail::as_value(pointer_, term); }
inline Value Session::value(double number) const { return detail::as_value(pointer_, number); }
inline Value Session::value(bool flag) const { return detail::as_value(pointer_, flag); }
inline Value Session::value(std::complex<double> number) const { return detail::as_value(pointer_, number); }
inline Value Session::value(long long number) const { return detail::as_value(pointer_, number); }
inline Value Session::value(int number) const { return detail::as_value(pointer_, number); }
template <class Integer, std::enable_if_t<std::is_integral_v<Integer>, int>>
inline Value Session::value(Integer number) const { return detail::as_value(pointer_, number); }
inline Value Session::value(const std::string& text) const { return detail::as_value(pointer_, text); }
inline Value Session::value(const char* text) const { return detail::as_value(pointer_, text); }
inline Value Session::value(const std::vector<double>& numbers) const {
    return detail::as_value(pointer_, numbers);
}

/// Call `op` through the value ABI with `count` prepared values.
inline Value Session::call_values(const char* op, const Value* args, std::size_t count) const {
    std::vector<std::uint32_t> cells(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (!args[index].valid() || args[index].session() != pointer_) {
            detail::fail(Status::NoValue, "argument is not a live value of this session");
        }
        cells[index] = args[index].cell();
    }
    std::uint32_t cell = 0;
    detail::check(rssn_value_call(pointer_, op, cells.empty() ? nullptr : cells.data(), count,
                                  &cell));
    return Value(pointer_, cell);
}

/// Call `op` with no arguments.
inline Value Session::call_value(const char* op) const { return call_values(op, nullptr, 0); }

namespace detail {

template <class Invoke, class... Args>
inline Value call_with_cells(SessionHandle session, Invoke&& invoke, Args&&... args) {
    std::vector<Value> owned;
    std::vector<std::uint32_t> cells;
    owned.reserve(sizeof...(Args));
    cells.reserve(sizeof...(Args));
    (
        [&] {
            if constexpr (std::is_same_v<std::decay_t<Args>, Value>) {
                if (!args.valid() || args.session() != session) {
                    detail::fail(Status::NoValue, "argument is not a live value of this session");
                }
                cells.push_back(args.cell());
            } else if constexpr (std::is_convertible_v<Args, Value>) {
                // A `Callback` converts to a borrowed value that names the
                // callback cell; the callback object keeps it alive.
                Value converted = args;
                if (!converted.valid() || converted.session() != session) {
                    detail::fail(Status::NoValue, "argument is not a live value of this session");
                }
                cells.push_back(converted.cell());
                owned.push_back(std::move(converted));
            } else {
                owned.push_back(as_value(session, std::forward<Args>(args)));
                cells.push_back(owned.back().cell());
            }
        }(),
        ...);
    std::uint32_t cell = 0;
    detail::check(invoke(cells.empty() ? nullptr : cells.data(), cells.size(), &cell));
    return Value(session, cell);
}

template <class... Args>
inline Value call_value_in(SessionHandle session, const char* op, Args&&... args) {
    return call_with_cells(session, [&](const std::uint32_t* cells, std::size_t count, std::uint32_t* out) {
        return rssn_value_call(session, op, cells, count, out);
    }, std::forward<Args>(args)...);
}

// A callback converts to a borrowed Value. Plain host arguments have no session.
template <class T>
inline SessionHandle value_owner(const T& argument) {
    if constexpr (std::is_same_v<T, Value> || std::is_same_v<T, Term>) {
        return argument.session();
    } else if constexpr (std::is_convertible_v<const T&, Value>) {
        return static_cast<Value>(argument).session();
    } else {
        return nullptr;
    }
}

template <class... Args>
inline Value call_inferred(const char* op, Args&&... args) {
    SessionHandle session = nullptr;
    const auto use_owner = [&](auto const& argument) {
        if (session == nullptr) session = value_owner(argument);
    };
    (use_owner(args), ...);
    if (session == nullptr) session = default_session().raw();
    return call_value_in(session, op, std::forward<Args>(args)...);
}

}  // namespace detail

template <class... Args>
inline Value Session::call_value(const char* op, Args&&... args) const {
    return detail::call_value_in(pointer_, op, std::forward<Args>(args)...);
}

template <class... Args>
inline Value Session::call_key(const OpKey& key, Args&&... args) const {
    if (key.rule_set.find('\0') != std::string::npos || key.name.find('\0') != std::string::npos) {
        detail::fail(Status::InvalidArgument, "operator key contains an embedded NUL");
    }
    return detail::call_with_cells(pointer_, [&](const std::uint32_t* cells, std::size_t count, std::uint32_t* out) {
        return rssn_value_call_key(pointer_, key.rule_set.c_str(), key.name.c_str(),
                                   key.signature_hash.data(), cells, count, out);
    }, std::forward<Args>(args)...);
}

inline std::vector<std::string> Session::available_backends() const {
    const Value choices = call_value("num_backend_available");
    const std::size_t count = choices.size();
    std::vector<std::string> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        result.push_back(choices.at(index).text());
    }
    return result;
}

inline std::string Session::current_backend() const {
    try {
        return call_value("num_backend_current_name").text();
    } catch (const Error& error) {
        if (error.status() != Status::UnknownOperator) throw;
    }

    const Value backend = call_value("num_backend_current");
    std::string name;
    try {
        name = call_value("num_backend_name", backend).text();
    } catch (...) {
        try {
            release_backend(backend);
        } catch (...) {
            // Preserve the original diagnostic.
        }
        throw;
    }
    release_backend(backend);
    return name;
}

namespace detail {

inline void check_backend_selector(const std::string& selector) {
    if (selector.find('\0') != std::string::npos) {
        fail(Status::InvalidArgument, "the backend selector contains NUL");
    }
}

}  // namespace detail

inline Value Session::open_backend(const std::string& selector) const {
    detail::check_backend_selector(selector);
    return call_value("num_backend_named", selector);
}

inline void Session::release_backend(const Value& backend) const {
    (void)call_value("num_backend_release", backend);
}

inline std::string Session::set_default_backend(const std::string& selector) const {
    detail::check_backend_selector(selector);
    return call_value("num_backend_select_default", selector).text();
}

inline std::vector<std::uint8_t> Session::encode_npy(
    const std::vector<std::size_t>& shape, const std::vector<double>& samples) const {
    std::vector<Value> extents;
    extents.reserve(shape.size());
    for (const std::size_t extent : shape) {
        extents.push_back(value_integer_decimal(pointer_, std::to_string(extent)));
    }
    const Value exemplar = extents.empty() ? value(0) : Value();
    const Value dimensions = value_tensor_cells(pointer_, extents, {extents.size()},
                                                extents.empty() ? &exemplar : nullptr);
    const std::string hex = call_value("num_npy_encode", dimensions, samples).text();
    if (hex.size() % 2 != 0) {
        detail::fail(Status::InvalidArgument, "invalid RSSN NumPy byte payload");
    }
    std::vector<std::uint8_t> bytes;
    bytes.reserve(hex.size() / 2);
    for (std::size_t pos = 0; pos < hex.size(); pos += 2) {
        unsigned int byte = 0;
        const char* first = hex.data() + pos;
        const char* last = first + 2;
        const auto parsed = std::from_chars(first, last, byte, 16);
        if (parsed.ec != std::errc() || parsed.ptr != last) {
            detail::fail(Status::InvalidArgument, "invalid RSSN NumPy byte payload");
        }
        bytes.push_back(static_cast<std::uint8_t>(byte));
    }
    return bytes;
}

inline Value Session::decode_npy(const std::vector<std::uint8_t>& payload) const {
    constexpr char digits[] = "0123456789abcdef";
    std::string hex;
    if (payload.size() > hex.max_size() / 2) {
        detail::fail(Status::InvalidArgument, "NumPy payload is too large");
    }
    hex.reserve(payload.size() * 2);
    for (const std::uint8_t byte : payload) {
        hex.push_back(digits[byte >> 4]);
        hex.push_back(digits[byte & 15]);
    }
    const Value record = call_value("num_npy_decode", hex);
    const Value dimensions = record.field("shape");
    const std::size_t rank = dimensions.size();
    std::vector<std::size_t> shape;
    shape.reserve(rank);
    for (std::size_t index = 0; index < rank; ++index) {
        const std::string text = dimensions.at(index).text();
        std::size_t extent = 0;
        const char* last = text.data() + text.size();
        const auto parsed = std::from_chars(text.data(), last, extent);
        if (parsed.ec != std::errc() || parsed.ptr != last) {
            detail::fail(Status::InvalidArgument, "invalid RSSN NumPy extent");
        }
        shape.push_back(extent);
    }
    return value_tensor(pointer_, record.field("data").tensor(), shape);
}

}  // namespace client
}  // namespace rssn

#endif  // RSSN_CLIENT_VALUE_HPP
