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
    out = Ref{UInt32}(0)
    _check(value_new_bool(session.ptr, x, out))
    return Value(session, out[])
end

function value(session::Session, x::Integer)
    out = Ref{UInt32}(0)
    if typemin(Int64) <= x <= typemax(Int64)
        _check(value_new_int(session.ptr, x, out))
    else
        _check(value_new_int_decimal(session.ptr, string(x), out))
    end
    return Value(session, out[])
end

function value(session::Session, x::Complex)
    out = Ref{UInt32}(0)
    _check(value_new_complex(session.ptr, real(x), imag(x), out))
    return Value(session, out[])
end

function value(session::Session, x::Real)
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
    value(session, x::Union{AbstractDict, NamedTuple}) -> Value

Creates a record using the input's iteration order. Use a named tuple when
field order matters. Both forms can serve as records for the value ABI.
"""
function value(session::Session, x::Union{AbstractDict, NamedTuple})
    names = String[String(key) for key in keys(x)]
    # The members own their cells: a comprehension of bare cells would let the
    # temporary `Value`s be finalized -- and their cells freed -- before the
    # record is built.
    members = Value[value(session, element) for element in values(x)]
    cells = UInt32[member.cell for member in members]
    out = Ref{UInt32}(0)
    GC.@preserve members _check(value_new_record(session.ptr, names, cells, out))
    return Value(session, out[])
end

"""
    value(session, x::Union{AbstractArray, Tuple}) -> Value

A floating-point array becomes a dense Float tensor, preserving rank and
empty extents. A vector of integers stays `tensor<int>`. Higher-rank real
arrays use the Float tensor constructor. Julia arrays are
column-major and the ABI is row-major, so the data is transposed on the way
out and [`tensor`](@ref) transposes it back.
"""
function value(session::Session, x::Union{AbstractArray, Tuple})
    if x isa AbstractArray && (eltype(x) <: AbstractFloat ||
            (ndims(x) >= 2 && eltype(x) <: Number))
        rank = ndims(x)
        shape = Csize_t[size(x, axis) for axis in 1:rank]
        storage = rank < 2 ? x : permutedims(x, rank:-1:1)
        data = Cdouble[vec(storage);]
        out = Ref{UInt32}(0)
        _check(value_new_tensor(session.ptr, data, shape, out))
        return Value(session, out[])
    end
    members = Value[value(session, element) for element in x]
    cells = UInt32[member.cell for member in members]
    out = Ref{UInt32}(0)
    GC.@preserve members _check(value_new_list(session.ptr, cells, out))
    return Value(session, out[])
end

"""
    tensor_cells(session::Session, x::AbstractArray; exemplar=nothing) -> Value

Copy an array into a tensor with its original rank and extents. Integers and
rationals retain their exact values. Supply an exemplar to choose the element
type of an empty array. Existing `Value` cells must belong to `session`.
Julia's column-major layout is converted to the ABI's row-major layout.
"""
function tensor_cells(session::Session, x::AbstractArray; exemplar=nothing)
    if exemplar !== nothing && !isempty(x)
        throw(ArgumentError("an element exemplar requires an empty tensor"))
    end
    as_cell(element) = element isa Rational ?
        call_value(session, "divide", numerator(element), denominator(element)) : value(session, element)
    rank = ndims(x)
    shape = Csize_t[size(x, axis) for axis in 1:rank]
    storage = rank < 2 ? x : permutedims(x, rank:-1:1)
    members = Value[
        if element isa Union{AbstractArray, Tuple}
            throw(ArgumentError("supply a rectangular array or explicit Value cells"))
        else
            as_cell(element)
        end for element in vec(storage)
    ]
    cells = UInt32[member.cell for member in members]
    out = Ref{UInt32}(0)
    proof = exemplar === nothing ? nothing : as_cell(exemplar)
    exemplar_cell = proof === nothing ? UInt32(0) : proof.cell
    GC.@preserve members proof _check(value_new_tensor_cells(session.ptr, cells, shape, exemplar_cell, out))
    return Value(session, out[])
end

"""An independently owned value cell in the same session."""
function Base.copy(source::Value)
    session = source.session
    entry = _optional_sym(:rssn_value_clone)
    if entry == C_NULL
        holder = value(session, (source,))
        try
            return holder[1]
        finally
            close(holder)
        end
    end
    out = Ref{UInt32}(0)
    GC.@preserve source _check(ccall(entry, Cint, (Ptr{Cvoid}, UInt32, Ptr{UInt32}),
                                     session.ptr, source.cell, out))
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
arguments have no term form, including records and native object handles.

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
    members = Value[value(session, arg) for arg in args]
    return _invoke_value(session, op, members)
end

"""
    call_key(session, key, args...) -> Value

Call the exact overload identified by a key from [`op_info`](@ref). Pass a
registered `Callback` as the value for a function parameter.
"""
function call_key(session::Session, key::AbstractDict, args...)
    rule_set = get(key, "rule_set", nothing)
    name = get(key, "name", nothing)
    if !(rule_set isa AbstractString && !isempty(rule_set) &&
         name isa AbstractString && !isempty(name))
        throw(ArgumentError("operator key requires a rule set and name"))
    end
    hash = get(key, "signature_hash", nothing)
    if !(hash isa AbstractVector && length(hash) == 16 &&
         all(byte -> byte isa Integer && 0 <= byte <= 255, hash))
        throw(ArgumentError("operator key requires 16 signature hash bytes"))
    end
    digest = UInt8[hash...]
    members = Value[value(session, argument) for argument in args]
    cells = UInt32[member.cell for member in members]
    out = Ref{UInt32}(0)
    GC.@preserve members _check(value_call_key(session.ptr, rule_set, name, digest, cells, out))
    return Value(session, out[])
end

function _invoke_value(session::Session, op::AbstractString, members::Vector{Value})
    cells = UInt32[member.cell for member in members]
    out = Ref{UInt32}(0)
    GC.@preserve members _check(value_call(session.ptr, op, cells, out))
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
            members = Value[]
            sizehint!(members, length(args))
            for (position, argument) in enumerate(args)
                if argument isa Function
                    slot = findfirst(==(position), positions)
                    declared = combination[slot]
                    callback = declared === nothing ?
                               Callback(session, argument; op = op, index = position) :
                               Callback(session, argument; signature = declared)
                    push!(callbacks, callback)
                    push!(members, callback.value)
                else
                    push!(members, value(session, argument))
                end
            end
            try
                return _invoke_value(session, op, members)
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

"""Attach a concrete dimension to a Float value already in coherent units.

Pass a dimension object such as `Dict("base" => "Length")`, or `"dimless"`.
The result owns a separate cell and keeps the tensor shape.
"""
function with_coherent_dimension(source::Value, dimension)
    out = Ref{UInt32}(0)
    payload = json_encode(dimension)
    GC.@preserve source _check(value_with_coherent_dimension(source.session.ptr, source.cell, payload, out))
    return Value(source.session, out[])
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
    shape, data = try
        extents = [Int(unsafe_load(raw.shape, index)) for index in 1:Int(raw.rank)]
        samples = [unsafe_load(raw.data, index) for index in 1:Int(raw.len)]
        extents, samples
    finally
        tensor_free(out)
    end
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

function _value_at(value::Value, index::Integer)
    # The C side counts from zero; a Julia index counts from one.
    out = Ref{UInt32}(0)
    _check(value_at(value.session.ptr, value.cell, index - 1, out))
    return Value(value.session, out[])
end

function Base.getindex(value::Value, index::Integer)
    1 <= index <= length(value) || throw(BoundsError(value, index))
    return _value_at(value, index)
end

function _value_field_name(record::Value, index::Integer)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(value_field_name(record.session.ptr, record.cell, index - 1, out))
    return take_string(out[])
end

"""Read one field of a record by its exact name."""
function Base.getindex(record::Value, name::AbstractString)
    kind(record) == "record" || throw(ArgumentError("named indexing requires an RSSN record"))
    for index in 1:length(record)
        _value_field_name(record, index) == name && return _value_at(record, index)
    end
    throw(KeyError(name))
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
