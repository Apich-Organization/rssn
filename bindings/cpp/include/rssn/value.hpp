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

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <string>
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

    /// The value as a term, when it has a graph form.
    Term to_term() const;

    /// A dense extraction of a tensor value; `shape` receives the extents.
    std::vector<double> tensor(std::vector<std::size_t>* shape = nullptr) const;

    /// The number of elements (a tensor's length, a record's field count).
    std::size_t size() const;

    /// Element `index` as a new value.
    Value at(std::size_t index) const;

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

inline Term Value::to_term() const {
    std::uint32_t handle = 0;
    detail::check(rssn_value_to_term(session_, cell_, &handle));
    return Term(session_, handle);
}

inline std::vector<double> Value::tensor(std::vector<std::size_t>* shape) const {
    RssnTensor raw{};
    detail::check(rssn_value_to_tensor(session_, cell_, &raw));
    std::vector<double> data(raw.data, raw.data + raw.len);
    if (shape != nullptr) {
        shape->assign(raw.shape, raw.shape + raw.rank);
    }
    rssn_tensor_free(&raw);
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

// ---------------------------------------------------------------------------
// Session: the value ABI
// ---------------------------------------------------------------------------

/// Store `value` as a cell.
inline Value Session::value(const Term& term) const {
    return detail::value_of(pointer_, [&](SessionHandle session, std::uint32_t* out) {
        return rssn_value_new_term(session, term.handle(), out);
    });
}

/// Store a number as a cell.
inline Value Session::value(double number) const {
    return detail::value_of(pointer_, [&](SessionHandle session, std::uint32_t* out) {
        return rssn_value_new_float(session, number, out);
    });
}

/// Store an integer as a cell.
inline Value Session::value(long long number) const {
    return detail::value_of(pointer_, [&](SessionHandle session, std::uint32_t* out) {
        return rssn_value_new_int(session, number, out);
    });
}

/// Store an `int` as a cell.
inline Value Session::value(int number) const {
    return value(static_cast<long long>(number));
}

/// Store text as a cell.
inline Value Session::value(const std::string& text) const {
    return detail::value_of(pointer_, [&](SessionHandle session, std::uint32_t* out) {
        return rssn_value_new_text(session, text.c_str(), out);
    });
}

/// Store a dense vector as a rank-1 tensor cell.
inline Value Session::value(const std::vector<double>& numbers) const {
    return value_tensor(pointer_, numbers, {numbers.size()});
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

/// Call `op` with arguments of any of the kinds `Session::value` accepts.
///
/// A `Value` argument is used as it is; a `Term`, a number or a string becomes
/// a temporary cell for the duration of the call.
template <class... Args>
inline Value Session::call_value(const char* op, Args&&... args) const {
    std::vector<Value> owned;
    std::vector<std::uint32_t> cells;
    owned.reserve(sizeof...(Args));
    cells.reserve(sizeof...(Args));
    (
        [&] {
            if constexpr (std::is_same_v<std::decay_t<Args>, Value>) {
                if (!args.valid() || args.session() != pointer_) {
                    detail::fail(Status::NoValue, "argument is not a live value of this session");
                }
                cells.push_back(args.cell());
            } else if constexpr (std::is_convertible_v<Args, Value>) {
                // A `Callback` converts to a borrowed value that names the
                // callback cell; the callback object keeps it alive.
                Value converted = args;
                if (!converted.valid() || converted.session() != pointer_) {
                    detail::fail(Status::NoValue, "argument is not a live value of this session");
                }
                cells.push_back(converted.cell());
                owned.push_back(std::move(converted));
            } else {
                owned.push_back(this->value(std::forward<Args>(args)));
                cells.push_back(owned.back().cell());
            }
        }(),
        ...);
    std::uint32_t cell = 0;
    detail::check(rssn_value_call(pointer_, op, cells.empty() ? nullptr : cells.data(),
                                  cells.size(), &cell));
    return Value(pointer_, cell);
}

}  // namespace client
}  // namespace rssn

#endif  // RSSN_CLIENT_VALUE_HPP
