# Dynamic values: the cell ABI behind operators that take native objects.
#
# The term API carries only values that are nodes of the expression graph. An
# operator that returns a native object — an AD tape, a compiled backend, a
# solver handle — returns a `Value` instead, and an operator that takes one can
# only be called with a `Value`. Those are the operators `OPS` marks with
# `ffi = false`; they are reached with [`call_value`](@ref).
#
#     session = Session()
#     backend = call_value(session, "num_backend_current")
#     text(call_value(session, "num_backend_name", backend))   # "interpreter"
#
#     tape = call_value(session, "num_ad_tape_new", 16)
#     x = call_value(session, "num_ad_tape_var", tape, 1.5)
#     z = call_value(session, "num_ad_var_mul", call_value(session, "num_ad_var_sin", x), x)
#     value(call_value(session, "num_ad_var_value", z))        # 1.4962424799060816

#: The kind codes `rssn_value_kind` reports, in declaration order.
const VALUE_KINDS = (
    "int", "real", "float", "complex", "bool", "string", "symbol", "unit",
    "quantity", "tensor", "sparse", "record", "variant", "term", "function",
    "dyn", "opaque",
)

"""
    Value

A dynamic value owned by a [`Session`](@ref): a handle into the session's cell
table. It is released by [`close`](@ref) or when it is collected, and only the
session that minted it can read it.
"""
mutable struct Value
    session::Session
    cell::UInt32
    owned::Bool

    function Value(session::Session, cell::Integer; owned::Bool = true)
        cell == 0 && throw(RssnError(Int32(9), last_error()))
        value = new(session, UInt32(cell), owned)
        finalizer(finalize_value, value)
        return value
    end
end


function finalize_value(value::Value)
    if value.cell != 0 && value.owned
        value_free(value.session.ptr, value.cell)
    end
    value.cell = 0
    return nothing
end

"""
    close(value)

Release the cell. The session releases the rest with itself. Safe to call twice.
"""
function Base.close(value::Value)
    finalize_value(value)
    return nothing
end

Base.show(io::IO, value::Value) =
    print(io, "Value(", value.cell == 0 ? "closed" : kind(value), ")")

"""
    value(session, x) -> Value

`x` as a cell of `session`: a `Value`, a `Term`, a number, a string, or a
vector of those.
"""
function value(session::Session, x::Value)
    x.session === session || throw(ArgumentError("value belongs to another session"))
    return x
end

function value(session::Session, x::Term)
    x.session === session || throw(ArgumentError("term belongs to another session"))
    out = Ref{UInt32}(0)
    _check(value_new_term(session.ptr, x.handle, out))
    return Value(session, out[])
end

function value(session::Session, x::Bool)
    return value(session, Int(x))
end

function value(session::Session, x::Integer)
    out = Ref{UInt32}(0)
    _check(value_new_int(session.ptr, x, out))
    return Value(session, out[])
end

function value(session::Session, x::AbstractFloat)
    out = Ref{UInt32}(0)
    _check(value_new_float(session.ptr, x, out))
    return Value(session, out[])
end

function value(session::Session, x::AbstractString)
    out = Ref{UInt32}(0)
    _check(value_new_text(session.ptr, x, out))
    return Value(session, out[])
end

"""
    value(session, x::Union{AbstractArray, Tuple}) -> Value

A rank-1 numeric vector becomes a list of its elements, which keeps each
element's own kind (an integer vector stays `tensor<int>`); a matrix or a
higher-rank numeric array becomes a dense tensor. Julia arrays are
column-major and the ABI is row-major, so the data is transposed on the way
out and [`tensor`](@ref) transposes it back.
"""
function value(session::Session, x::Union{AbstractArray, Tuple})
    if x isa AbstractArray && ndims(x) >= 2 && eltype(x) <: Number
        shape = Csize_t[size(x, dimension) for dimension in 1:ndims(x)]
        data = Cdouble[vec(permutedims(Array(x)));]
        out = Ref{UInt32}(0)
        _check(value_new_tensor(session.ptr, data, shape, out))
        return Value(session, out[])
    end
    cells = UInt32[value(session, element).cell for element in x]
    out = Ref{UInt32}(0)
    _check(value_new_list(session.ptr, cells, out))
    return Value(session, out[])
end

"""Release ownership of the cell, which the caller then owns."""
function detach(value::Value)
    cell = value.cell
    value.cell = 0
    return cell
end

"""
    call_value(session, op, args...) -> Value

Call `op` through the value ABI. This is the path to the operators whose
arguments have no term form: `OPS[op][3]` is `false` for exactly those.

Arguments may be `Value`s, `Term`s, numbers, strings, vectors, arrays, dicts
(records) — and functions. A `Function` argument becomes a host callback with
the signature the operator declares for that position:

```julia
call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64, x -> x * x - 2)
```

The `do` form puts the function last, which is where the callback sits in most
of these operators:

```julia
call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x
    x * x - 2
end
```
"""
function call_value(session::Session, op::AbstractString, args...)
    if any(argument -> argument isa Function, args)
        return _call_with_callbacks(session, op, args)
    end
    cells = UInt32[value(session, arg).cell for arg in args]
    out = Ref{UInt32}(0)
    _check(value_call(session.ptr, op, cells, out))
    return Value(session, out[])
end

"""
    call_value(f::Function, session, op, args...)

The `do`-block form: `call_value(session, op, args...) do ... end` reaches this
method, with the function first and the remaining arguments in order. The
function is appended as the last argument of `op`.
"""
call_value(function_::Function, session::Session, op::AbstractString, args...) =
    call_value(session, op, args..., function_)

"""
Run `op` with every `Function` in `args` registered as a host callback.

A few operators declare more than one function type at one position (a residual
that answers a number or a vector, say); the declared ones are tried in order
and the first the engine accepts wins.
"""
function _call_with_callbacks(session::Session, op::AbstractString, args)
    positions = findall(argument -> argument isa Function, args)
    options = [callback_signatures(session, op, index) for index in positions]
    for option in options
        isempty(option) && (push!(option, nothing))
    end
    first_error = nothing
    attempts = 0
    for combination in Iterators.product(options...)
        attempts += 1
        attempts > 8 && break
        callbacks = Callback[]
        try
            cells = UInt32[]
            for (position, argument) in enumerate(args)
                if argument isa Function
                    slot = findfirst(==(position), positions)
                    declared = combination[slot]
                    callback = declared === nothing ?
                               Callback(session, argument; op = op, index = position) :
                               Callback(session, argument; signature = declared)
                    push!(callbacks, callback)
                    push!(cells, callback.value.cell)
                else
                    push!(cells, value(session, argument).cell)
                end
            end
            out = Ref{UInt32}(0)
            try
                _check(value_call(session.ptr, op, cells, out))
                return Value(session, out[])
            catch error
                first_error === nothing && (first_error = error)
                continue
            end
        finally
            for callback in callbacks
                close(callback)
            end
        end
    end
    first_error === nothing && throw(ArgumentError("no callback signature to try for $op"))
    throw(first_error)
end

"""The kind of a value: `"float"`, `"tensor"`, `"opaque"`, ..."""
function kind(value::Value)
    out = Ref{UInt8}(0)
    _check(value_kind(value.session.ptr, value.cell, out))
    code = Int(out[])
    return code < length(VALUE_KINDS) ? VALUE_KINDS[code + 1] : "kind$code"
end

"""The numeric value of a value."""
function value(number::Value)
    out = Ref{Cdouble}(0.0)
    _check(value_to_f64(number.session.ptr, number.cell, out))
    return out[]
end

"""The text of a string or symbol value, or any value's display form."""
function text(value::Value)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(value_to_text(value.session.ptr, value.cell, out))
    return take_string(out[])
end

"""The value in its `serde` JSON form, decoded."""
function json(value::Value)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(value_json(value.session.ptr, value.cell, out))
    return json_decode(take_string(out[]))
end

"""The value as a term, when it has a graph form."""
function term(value::Value)
    out = Ref{UInt32}(0)
    _check(value_to_term(value.session.ptr, value.cell, out))
    return Term(value.session, out[])
end

"""A dense extraction of a tensor value, reshaped to its extents."""
function tensor(value::Value)
    out = Ref{CTensor}()
    _check(value_to_tensor(value.session.ptr, value.cell, out))
    raw = out[]
    shape = [Int(unsafe_load(raw.shape, index)) for index in 1:Int(raw.rank)]
    data = [unsafe_load(raw.data, index) for index in 1:Int(raw.len)]
    tensor_free(out)
    if isempty(shape)
        return fill(isempty(data) ? 0.0 : data[1])
    end
    dims = Tuple(shape)
    filled = reshape(data, reverse(dims))
    return length(dims) == 1 ? filled : permutedims(filled, length(dims):-1:1)
end

function Base.length(value::Value)
    out = Ref{Csize_t}(0)
    _check(value_len(value.session.ptr, value.cell, out))
    return Int(out[])
end

function Base.getindex(value::Value, index::Integer)
    # The C side counts from zero; a Julia index counts from one.
    1 <= index <= length(value) || throw(BoundsError(value, index))
    out = Ref{UInt32}(0)
    _check(value_at(value.session.ptr, value.cell, index - 1, out))
    return Value(value.session, out[])
end

Base.firstindex(value::Value) = 1
Base.lastindex(value::Value) = length(value)

"""
    iterate(value::Value, state = 1)

Elements of a tensor or record value, so `for element in value` works. Each
element is a new value owned by the caller.
"""
function Base.iterate(value::Value, state::Int = 1)
    state > length(value) && return nothing
    return (value[state], state + 1)
end

"""The number of live value cells of `session`; a diagnostic for leak hunting."""
value_count(session::Session) = Int(value_count(session.ptr))
