# Host callbacks: Julia functions the engine can call.
#
# Some operators are higher-order: `num_bfgs`, `num_brent`,
# `num_analysis_integrate_adaptive_scaled`, `num_bootstrap` and about 180 others
# take a function and call it. In Julia any function works, because this file
# turns one into the C function pointer the engine calls.
#
#     session = Session()
#     answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64,
#                         x -> x * x - 2)
#     answer["root"]                     # 1.4142135623730951
#
# The `do` form puts the function last, which is where the callback sits in most
# of these operators:
#
#     answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x
#         x * x - 2
#     end
#
# A callback receives native Julia objects — `Float64`, `Int`, `Bool`,
# `Complex`, `String`, arrays and `Dict`s for records — and whatever it returns
# is converted back into a value.

#: The status a callback returns when the Julia function raised.
#: `RssnStatus::InvalidArgument`, the same code `RssnError` reports.
const CALLBACK_FAILED = Cint(13)

#: Registered callbacks, keyed by the integer the engine carries as `user`.
#: The value keeps both the session and the Julia function alive.
const _CALLBACKS = Dict{UInt, Tuple{Session, Any}}()

const _NEXT_CALLBACK = Ref(UInt(0))
const _CALLBACK_FN = Ref{Ptr{Cvoid}}(C_NULL)

#: Scalar spellings the compact signature grammar accepts.
const _SIGNATURE_SCALARS = Dict(
    "int" => "int", "integer" => "int",
    "real" => "real", "rational" => "real",
    "float" => "float", "double" => "float",
    "complex" => "complex",
    "bool" => "bool", "boolean" => "bool",
    "string" => "string", "str" => "string", "text" => "string",
    "symbol" => "symbol",
)

"""
    parse_signature(text) -> String

`"(float, float) -> float"` as the JSON object `rssn_callback_new` reads.

The arrow may be `->` or `=>`; `[]` appends a tensor rank (`"float[]"` is a
vector, `"float[][]"` a matrix). A string that is already JSON is returned
unchanged.
"""
function parse_signature(text::AbstractString)
    stripped = strip(text)
    (startswith(stripped, "{") || startswith(stripped, "[")) && return String(stripped)
    arrow = occursin("->", stripped) ? "->" : "=>"
    occursin(arrow, stripped) ||
        throw(ArgumentError("a callback signature needs `->`: $(repr(text))"))
    left, right = split(stripped, arrow; limit = 2)
    left = strip(left)
    if startswith(left, "(") && endswith(left, ")")
        left = strip(left[2:(end - 1)])
    end
    params = isempty(left) ? String[] : [_type_json(part) for part in _split_top(left)]
    return string("{\"params\":[", join(params, ","), "],\"ret\":", _type_json(right), "}")
end

"""Split on commas that are not inside brackets."""
function _split_top(text::AbstractString)
    parts = String[]
    depth = 0
    current = IOBuffer()
    for character in text
        if character in "([{"
            depth += 1
        elseif character in ")]}"
            depth -= 1
        end
        if character == ',' && depth == 0
            push!(parts, strip(String(take!(current))))
        else
            print(current, character)
        end
    end
    tail = strip(String(take!(current)))
    isempty(tail) || push!(parts, tail)
    return parts
end

"""One compact type (`real`, `float[]`) as a `TypeExpr` JSON fragment."""
function _type_json(token::AbstractString)
    name = strip(token)
    rank = 0
    while endswith(name, "[]")
        rank += 1
        name = strip(name[1:(end - 2)])
    end
    base = get(_SIGNATURE_SCALARS, lowercase(name), nothing)
    base === nothing && throw(ArgumentError(
        "unknown type $(repr(name)) in a callback signature; use one of " *
        join(sort(unique(collect(values(_SIGNATURE_SCALARS)))), ", ") *
        " or append [] for a tensor"))
    rank == 0 && return string('"', base, '"')
    return string("{\"tensor\":{\"rank\":{\"lit\":", rank, "},\"shape\":null,\"elem\":\"",
                  base, "\"}}")
end

"""
    Callback

A Julia function registered with a [`Session`](@ref) as a host callback.

```julia
f = Callback(session, x -> x * x; signature = "(float) -> float")
g = Callback(session, x -> x * x; op = "num_brent", index = 4)
```

A callback lives until [`close`](@ref) or until the session is released.
"""
mutable struct Callback
    session::Session
    value::Value
    function_::Any
    id::UInt
end

"""
    Callback(session, f; signature = nothing, op = nothing, index = 0)

Register `f` as a host callback. Give either `signature` (a compact
`"(float) -> float"` string or the JSON object) or the operator and argument
position the callback is for, in which case the engine supplies the type from
that operator's declaration.
"""
function Callback(session::Session, function_::Function; signature = nothing,
                  op::Union{Nothing, AbstractString} = nothing, index::Integer = 0)
    id = (_NEXT_CALLBACK[] += 1)
    _CALLBACKS[id] = (session, function_)
    out = Ref{UInt32}(0)
    pointer = _callback_pointer()
    try
        if signature !== nothing
            text = signature isa AbstractString ? parse_signature(signature) :
                   json_encode(signature)
            _check(callback_new(session.ptr, text, pointer, Ptr{Cvoid}(id), out))
        elseif op !== nothing
            _check(callback_new_for(session.ptr, op, index, pointer, Ptr{Cvoid}(id), out))
        else
            throw(ArgumentError("a Callback needs a signature, or an operator and an index"))
        end
    catch
        delete!(_CALLBACKS, id)
        rethrow()
    end
    return Callback(session, Value(session, out[]), function_, id)
end

"""
The trampoline the engine calls. It is a plain C function, so it cannot capture
anything: the Julia function is found through the integer the engine carries as
`user`.
"""
function _callback_thunk(user::Ptr{Cvoid}, _session::Ptr{Cvoid}, args::Ptr{UInt32},
                         count::Csize_t, out::Ptr{UInt32})::Cint
    entry = get(_CALLBACKS, UInt(user), nothing)
    entry === nothing && return CALLBACK_FAILED
    session, function_ = entry
    try
        values = Any[to_julia(session, unsafe_load(args, index)) for index in 1:Int(count)]
        cell = from_julia(session, function_(values...))
        unsafe_store!(out, UInt32(cell))
        return Cint(0)
    catch error
        try
            callback_fail(string(typeof(error), ": ", sprint(showerror, error)))
        catch
        end
        return CALLBACK_FAILED
    end
end

function _callback_pointer()
    if _CALLBACK_FN[] == C_NULL
        _CALLBACK_FN[] = @cfunction(_callback_thunk, Cint,
                                    (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{UInt32}))
    end
    return _CALLBACK_FN[]
end



"""
    close(callback)

Retire the callback and release its cell. Safe to call twice.
"""
function Base.close(callback::Callback)
    if callback.value.cell != 0
        callback_free(callback.session.ptr, callback.value.cell)
        callback.value.cell = 0
    end
    delete!(_CALLBACKS, callback.id)
    return nothing
end

"""
    value(session, callback::Callback) -> Value

The callback as a value, to pass to [`call_value`](@ref).
"""
function value(session::Session, callback::Callback)
    callback.session === session ||
        throw(ArgumentError("callback belongs to another session"))
    return callback.value
end

Base.show(io::IO, callback::Callback) =
    print(io, "Callback(", callback.value.cell == 0 ? "closed" : string(callback.value.cell), ")")

"""
    callback_signatures(session, op, index) -> Vector

The signatures `op` declares for a callback at argument `index` (Julia counts
from one). Each element is the JSON object [`Callback`](@ref) accepts.
"""
function callback_signatures(session::Session, op::AbstractString, index::Integer)
    index >= 1 || throw(ArgumentError("Julia argument positions start at 1"))
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(callback_signatures(session.ptr, op, index - 1, out))
    listed = json_decode(take_string(out[]))
    return listed isa AbstractVector ? listed : Any[]
end

# ---------------------------------------------------------------------------
# Julia <-> value conversion
# ---------------------------------------------------------------------------

"""
    to_julia(session, cell)

The Julia object a value cell holds. This is what a callback receives for each
argument: numbers, `Bool`, `Complex`, `String`, arrays, and `Dict`s for records.
"""
function to_julia(session::Session, cell::Integer)
    wrapped = Value(session, cell; owned = false)
    kind_ = kind(wrapped)
    kind_ == "int" && return Int(round(value(wrapped)))
    (kind_ == "real" || kind_ == "float") && return value(wrapped)
    kind_ == "bool" && return Bool(json(wrapped)["bool"])
    if kind_ == "complex"
        body = json(wrapped)["complex"]
        if haskey(body, "Exact")
            exact = body["Exact"]
            return Complex(_number(exact["re"]), _number(exact["im"]))
        end
        floating = body["Float"]
        return Complex(Float64(floating["re"]), Float64(floating["im"]))
    end
    (kind_ == "string" || kind_ == "symbol") && return text(wrapped)
    if kind_ == "tensor"
        try
            return tensor(wrapped)
        catch error
            error isa RssnError || rethrow()
            return Any[to_julia(session, element.cell) for element in wrapped]
        end
    end
    if kind_ == "record"
        fields = Dict{String, Any}()
        for index in 1:length(wrapped)
            out = Ref{Ptr{UInt8}}(C_NULL)
            _check(value_field_name(session.ptr, wrapped.cell, index - 1, out))
            fields[take_string(out[])] = to_julia(session, wrapped[index].cell)
        end
        return fields
    end
    if kind_ == "variant"
        out = Ref{Ptr{UInt8}}(C_NULL)
        _check(value_tag(session.ptr, wrapped.cell, out))
        tag = take_string(out[])
        return (tag, Any[to_julia(session, element.cell) for element in wrapped])
    end
    return json(wrapped)
end

"""A `Number` JSON (`Dict("Float" => 1.0)`, `Dict("Int" => 3)`, ...) as a float."""
function _number(payload)
    payload isa Number && return Float64(payload)
    if payload isa AbstractDict && length(payload) == 1
        tag, body = first(payload)
        tag in ("Float", "float") && return Float64(body)
        tag in ("Int", "int") && return Float64(body)
        if tag in ("Rat", "rat")
            body isa AbstractDict && return Float64(body["numer"]) / Float64(body["denom"])
            return Float64(body[1]) / Float64(body[2])
        end
    end
    throw(ArgumentError("not a number: $(repr(payload))"))
end

"""
    from_julia(session, x) -> UInt32

A cell holding `x`, which the caller then owns. This is how a callback's return
value crosses back.
"""
function from_julia(session::Session, x)
    out = Ref{UInt32}(0)
    x isa Value && return x.cell
    x isa Callback && return x.value.cell
    if x isa Term
        x.session === session || throw(ArgumentError("term belongs to another session"))
        _check(value_new_term(session.ptr, x.handle, out))
        return out[]
    end
    if x isa Bool
        _check(value_new_bool(session.ptr, x, out))
        return out[]
    end
    if x isa Integer
        _check(value_new_int(session.ptr, x, out))
        return out[]
    end
    if x isa Real
        _check(value_new_float(session.ptr, Float64(x), out))
        return out[]
    end
    if x isa Complex
        _check(value_new_complex(session.ptr, real(x), imag(x), out))
        return out[]
    end
    if x isa AbstractString
        _check(value_new_text(session.ptr, x, out))
        return out[]
    end
    if x isa AbstractDict
        names = String[String(key) for key in keys(x)]
        cells = UInt32[from_julia(session, element) for element in values(x)]
        _check(value_new_record(session.ptr, names, cells, out))
        return out[]
    end
    if x isa AbstractArray
        if x isa AbstractVector && all(element -> element isa Union{Bool, Integer, AbstractFloat, AbstractString}, x)
            # A vector of scalars: keep each element's own kind, exactly as
            # `value` does for an argument.
            return detach(value(session, x))
        end
        if eltype(x) <: Number && ndims(x) >= 2
            shape = Csize_t[size(x, dimension) for dimension in 1:ndims(x)]
            data = Cdouble[vec(permutedims(Array(x)));]
            _check(value_new_tensor(session.ptr, data, shape, out))
            return out[]
        end
        cells = UInt32[from_julia(session, element) for element in x]
        _check(value_new_list(session.ptr, cells, out))
        return out[]
    end
    throw(ArgumentError("cannot use $(typeof(x)) as a callback result"))
end
