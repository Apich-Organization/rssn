// rssn C++ client: host callbacks.
//
// Some operators are higher-order: `num_bfgs`, `num_brent`,
// `num_analysis_integrate_adaptive_scaled`, `num_bootstrap` and about 180 others
// take a function and call it. A `Callback` is how a C++ host supplies one: it
// wraps a `std::function` in the C function pointer the engine calls.
//
//     #include "rssn/client.hpp"
//     #include "rssn/callback.hpp"
//     #include "rssn/ops.hpp"
//
//     rssn::client::Session session;
//     auto f = rssn::client::Callback::for_argument(
//         session.raw(), "num_brent", 4,
//         [](rssn::client::SessionHandle s, const std::vector<rssn::client::Value>& args) {
//             const double x = args[0].as_f64();
//             return rssn::client::value_number(s, x * x - 2.0);
//         });
//     auto answer = num_brent(1.0, 2.0, 1e-12, 64, f);
//
// The body receives the session it was registered on (so it can build its
// result) and the arguments as borrowed `Value`s. It may throw: the exception's
// message becomes the failure the enclosing operator reports.

#ifndef RSSN_CLIENT_CALLBACK_HPP
#define RSSN_CLIENT_CALLBACK_HPP

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "rssn/client.hpp"
#include "rssn/value.hpp"

namespace rssn {
namespace client {

/// The body of a host callback: the session it was registered on, and the
/// arguments the engine passed.
using CallbackBody = std::function<Value(SessionHandle, const std::vector<Value>&)>;

namespace detail {

struct CallbackState {
    CallbackBody body;
};

/// The C function the engine calls. It cannot capture, so the body travels
/// through `user`.
extern "C" inline RssnStatus callback_thunk(void* user, const RssnSession* session,
                                            const std::uint32_t* args, std::size_t count,
                                            std::uint32_t* out) {
    auto* state = static_cast<CallbackState*>(user);
    if (state == nullptr || out == nullptr) {
        return RSSN_STATUS_NULL_ARGUMENT;
    }
    try {
        std::vector<Value> values;
        values.reserve(count);
        for (std::size_t index = 0; index < count; ++index) {
            values.push_back(Value::borrowed(const_cast<RssnSession*>(session), args[index]));
        }
        Value result = state->body(const_cast<RssnSession*>(session), values);
        if (!result.valid()) {
            rssn_callback_fail("the callback returned no value");
            return RSSN_STATUS_INVALID_ARGUMENT;
        }
        // The engine reads the cell and frees it; the callback stops owning it.
        *out = result.detach();
        return RSSN_STATUS_OK;
    } catch (const std::exception& error) {
        rssn_callback_fail(error.what());
        return RSSN_STATUS_INVALID_ARGUMENT;
    } catch (...) {
        rssn_callback_fail("the callback threw an unknown exception");
        return RSSN_STATUS_INVALID_ARGUMENT;
    }
}

}  // namespace detail

/// A C++ function registered with a session as a host callback.
///
/// Move-only: it owns one callback registration and releases it when destroyed.
class Callback {
public:
    /// A callback that names nothing.
    Callback() = default;

    ~Callback() { release(); }

    Callback(const Callback&) = delete;
    Callback& operator=(const Callback&) = delete;

    Callback(Callback&& other) noexcept
        : session_(other.session_), cell_(other.cell_), state_(std::move(other.state_)) {
        other.session_ = nullptr;
        other.cell_ = 0;
    }

    Callback& operator=(Callback&& other) noexcept {
        if (this != &other) {
            release();
            session_ = other.session_;
            cell_ = other.cell_;
            state_ = std::move(other.state_);
            other.session_ = nullptr;
            other.cell_ = 0;
        }
        return *this;
    }

    /// Registers `body` with an explicit signature: `"(float) -> float"`, or the
    /// JSON object `{"params": [...], "ret": ...}`.
    static Callback create(SessionHandle session, CallbackBody body, const char* signature) {
        Callback callback;
        callback.state_ = std::make_shared<detail::CallbackState>();
        callback.state_->body = std::move(body);
        std::uint32_t cell = 0;
        detail::check(rssn_callback_new(session, signature, &detail::callback_thunk,
                                        callback.state_.get(), &cell));
        callback.session_ = session;
        callback.cell_ = cell;
        return callback;
    }

    /// Registers `body` as the callback for argument `index` of operator `op`,
    /// taking the signature from the operator's own declaration.
    static Callback for_argument(SessionHandle session, const char* op, std::size_t index,
                                 CallbackBody body) {
        Callback callback;
        callback.state_ = std::make_shared<detail::CallbackState>();
        callback.state_->body = std::move(body);
        std::uint32_t cell = 0;
        detail::check(rssn_callback_new_for(session, op, index, &detail::callback_thunk,
                                            callback.state_.get(), &cell));
        callback.session_ = session;
        callback.cell_ = cell;
        return callback;
    }

    /// Whether this callback names a live registration.
    bool valid() const noexcept { return session_ != nullptr && cell_ != 0; }

    /// The raw cell id.
    std::uint32_t cell() const noexcept { return cell_; }

    /// The callback as a value, to pass to `Session::call_value`.
    Value value() const { return Value::borrowed(session_, cell_); }

    /// Implicit conversion, so a `Callback` can be passed where a value is
    /// expected.
    operator Value() const { return value(); }

    /// Retire the callback now rather than at destruction.
    void close() { release(); }

private:
    void release() noexcept {
        if (session_ != nullptr && cell_ != 0) {
            rssn_callback_free(session_, cell_);
        }
        session_ = nullptr;
        cell_ = 0;
        state_.reset();
    }

    SessionHandle session_ = nullptr;
    std::uint32_t cell_ = 0;
    std::shared_ptr<detail::CallbackState> state_;
};

}  // namespace client
}  // namespace rssn

#endif  // RSSN_CLIENT_CALLBACK_HPP
