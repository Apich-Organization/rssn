# Terms and sessions: the owner of every handle, and the values it holds.

"""
    RssnError(status, message)

A failed rssn call. `status` is the `RssnStatus` code (`0` is success).
"""
struct RssnError <: Exception
    status::Int32
    message::String
end

function Base.showerror(io::IO, error::RssnError)
    print(io, "RssnError(", _STATUS_NAMES[error.status + 1], "): ", error.message)
end

const _STATUS_NAMES = [
    "ok", "null_argument", "invalid_utf8", "parse", "unknown_operator", "not_numeric",
    "inconsistent", "rule", "no_such_term", "no_value", "simulation", "panic",
    "dimension", "invalid_argument", "cancelled",
]

@inline function _check(code::Integer)
    code == 0 && return nothing
    throw(RssnError(Int32(code), last_error()))
end

"""
    Session()

An owned `RssnSession`: the owner of every term made through it.

A session is single-threaded — its store sits behind a `RefCell` in Rust, so
two tasks using one session at once is a borrow conflict. Use one session per
task, or [`compute_all`](@ref) for parallel requests inside one session.

```julia
session = Session()
x = sym(session, "x")
answer = compute(session, diff(sin(x) * cos(x), x))
answer.text            # "cos(2*x)"
close(session)
```
"""
mutable struct Session
    ptr::Ptr{Cvoid}

    function Session()
        pointer = session_new()
        pointer == C_NULL && throw(RssnError(Int32(7), last_error()))
        session = new(pointer)
        finalizer(finalize_session, session)
        return session
    end
end

function finalize_session(session::Session)
    if session.ptr != C_NULL
        session_free(session.ptr)
        session.ptr = C_NULL
    end
    return nothing
end

"""
    close(session)

Release the session and every term made through it. Safe to call twice.

This extends `Base.close`, so no import is needed.
"""
function Base.close(session::Session)
    finalize_session(session)
    return nothing
end

Base.show(io::IO, session::Session) =
    print(io, "Session(", session.ptr == C_NULL ? "closed" : string(session.ptr), ")")

"""
    Term

A node of a [`Session`](@ref)'s graph: a `UInt32` handle and the session that
owns it. Terms are immutable values; `x * x` makes a new node.

Operators are requests, not results: `diff(sin(x) * x, x)` builds the
unevaluated derivative and [`compute`](@ref) reduces it.
"""
struct Term
    session::Session
    handle::UInt32
end

Base.show(io::IO, term::Term) = print(io, "Term(", repr(text(term)), ")")

function Base.:(==)(left::Term, right::Term)
    return left.session === right.session && left.handle == right.handle
end
Base.hash(term::Term, seed::UInt) = hash(term.handle, hash(objectid(term.session), seed))

"""The raw `UInt32` handle of a term."""
handle(term::Term) = term.handle

"""
    text(term) -> String

The term as infix text, parseable by [`parse_term`](@ref).
"""
text(term::Term) = take_string(term_to_string(term.session.ptr, term.handle))

"""The term as LaTeX math."""
latex(term::Term) = take_string(term_to_latex(term.session.ptr, term.handle))

"""The term as Typst math."""
typst(term::Term) = take_string(term_to_typst(term.session.ptr, term.handle))

"""The term as a multi-line Unicode drawing."""
pretty(term::Term) = take_string(term_to_pretty(term.session.ptr, term.handle))

"""
    value(term) -> Union{Float64, Nothing}

The value of a literal number, or `nothing` for a symbolic term.
"""
function value(term::Term)
    out = Ref{Cdouble}(0.0)
    code = term_to_f64(term.session.ptr, term.handle, out)
    code == 0 && return out[]
    code == 9 && return nothing    # RSSN_STATUS_NO_VALUE
    throw(RssnError(Int32(code), last_error()))
end

"""
    tensor(term) -> Array{Float64}

A dense `f64` extraction of a list term, reshaped to its extents. A scalar
term gives a zero-dimensional array.
"""
function tensor(term::Term)
    out = Ref{CTensor}()
    _check(term_tensor(term.session.ptr, term.handle, out))
    raw = out[]
    shape = [Int(unsafe_load(raw.shape, index)) for index in 1:Int(raw.rank)]
    data = [unsafe_load(raw.data, index) for index in 1:Int(raw.len)]
    tensor_free(out)
    if isempty(shape)
        return fill(isempty(data) ? 0.0 : data[1])
    end
    # The C side is row-major; Julia fills column-major, so reshape with the
    # reversed extents and permute the axes back.
    dims = Tuple(shape)
    filled = reshape(data, reverse(dims))
    return length(dims) == 1 ? filled : permutedims(filled, length(dims):-1:1)
end

"""
    to_float(term) -> Float64

The numeric value of the term, computing it when the term is not a literal.
"""
function to_float(term::Term)
    literal = value(term)
    literal !== nothing && return literal
    config = numeric!(Config(), 1e-12)
    answer = compute(term, config)
    close(config)
    answer.value === nothing && throw(RssnError(Int32(5), "term has no numeric value"))
    return answer.value
end

Base.Float64(term::Term) = to_float(term)

# -- construction -----------------------------------------------------------

"""The symbol `name`."""
sym(session::Session, name::AbstractString) =
    _wrap(session, term_sym(session.ptr, name), Int32(4))

"""An exact integer."""
function integer(session::Session, x::Integer)
    if typemin(Int64) <= x <= typemax(Int64)
        return _wrap(session, term_int(session.ptr, x), Int32(13))
    end
    exact = value(session, x)
    try
        return term(exact)
    finally
        close(exact)
    end
end

"""A floating point number."""
number(session::Session, value::Real) =
    _wrap(session, term_float(session.ptr, value), Int32(13))

"""The exact fraction `numerator / denominator`."""
function rational(session::Session, numerator::Integer, denominator::Integer = 1)
    handle = term_rational(session.ptr, numerator, denominator)
    handle == NO_TERM && throw(RssnError(Int32(13), last_error()))
    return Term(session, handle)
end

"""
    parse_term(session, source) -> Term

Parse infix text such as `"diff(sin(x)^2, x)"`.
"""
function parse_term(session::Session, source::AbstractString)
    handle = term_parse(session.ptr, source)
    handle == NO_TERM && throw(RssnError(Int32(3), last_error()))
    return Term(session, handle)
end

"""
    term(session, value) -> Term

`value` as a term: a `Term`, a number, a `Rational`, or infix text.
"""
function term(session::Session, value::Term)
    value.session === session ||
        throw(ArgumentError("term belongs to another session"))
    return value
end

function term(session::Session, value::Bool)
    return integer(session, Int(value))
end

term(session::Session, value::Integer) = integer(session, value)
term(session::Session, value::Rational) = rational(session, numerator(value), denominator(value))
term(session::Session, value::AbstractFloat) = number(session, value)
term(session::Session, value::AbstractString) = parse_term(session, value)
term(session::Session, value) =
    throw(ArgumentError("cannot use $(typeof(value)) as a term"))

function _wrap(session::Session, handle::UInt32, status::Int32)
    handle == NO_TERM && throw(RssnError(status, last_error()))
    return Term(session, handle)
end

"""
    apply(session, op, args...) -> Term

Apply the operator `op` to `args`, as an unevaluated request.
"""
function apply(session::Session, op::AbstractString, args...)
    handles = UInt32[term(session, value).handle for value in args]
    handle = term_apply(session.ptr, op, handles)
    handle == NO_TERM && throw(RssnError(Int32(4), last_error()))
    return Term(session, handle)
end

"""
    _apply(op, args, session) -> Term

The operator facade's entry point: use `session` when it is given, else the
first term's session, else the default session.
"""
function _apply(op::AbstractString, args::Tuple, session)
    if session === nothing
        for value in args
            if value isa Term
                session = value.session
                break
            end
        end
        session === nothing && (session = default_session())
    end
    return apply(session, op, args...)
end

"""
    _call(op, args, session) -> Value

The facade's entry point for a higher-order operator: the same session
resolution as [`_apply`](@ref), but through the value ABI, which is the one that
can carry a host callback. A `Function` in `args` is registered with the
signature `op` declares for that position.
"""
function _call(op::AbstractString, args::Tuple, session)
    if session === nothing
        for value in args
            if value isa Union{Term, Value, Callback}
                session = value.session
                break
            end
        end
        session === nothing && (session = default_session())
    end
    return call_value(session, op, args...)
end

# -- arithmetic, as requests ------------------------------------------------

#: What an operator accepts besides a `Term`: numbers and infix text.
const Coercible = Union{Number, AbstractString}

# `sub`, `neg` and `divide` are the engine's named forms of subtraction,
# negation and quotient: it rewrites them to the canonical `add` / `mul` /
# `pow` shape and prints both spellings the same way. The arithmetic quotient
# is `divide` because `div` is the divergence of a vector field.
_difference(session, left, right) = apply(session, "sub", left, right)
_quotient(session, left, right) = apply(session, "divide", left, right)

Base.:+(left::Term, right::Term) = _apply("add", (left, right), nothing)
Base.:+(left::Term, right::Coercible) = _apply("add", (left, right), nothing)
Base.:+(left::Coercible, right::Term) = _apply("add", (left, right), nothing)

Base.:*(left::Term, right::Term) = _apply("mul", (left, right), nothing)
Base.:*(left::Term, right::Coercible) = _apply("mul", (left, right), nothing)
Base.:*(left::Coercible, right::Term) = _apply("mul", (left, right), nothing)

Base.:-(left::Term, right::Term) = _difference(left.session, left, right)
Base.:-(left::Term, right::Coercible) = _difference(left.session, left, right)
Base.:-(left::Coercible, right::Term) = _difference(right.session, left, right)

Base.:/(left::Term, right::Term) = _quotient(left.session, left, right)
Base.:/(left::Term, right::Coercible) = _quotient(left.session, left, right)
Base.:/(left::Coercible, right::Term) = _quotient(right.session, left, right)

Base.:-(value::Term) = apply(value.session, "neg", value)
Base.:+(value::Term) = value

Base.:^(base::Term, exponent::Term) = _apply("pow", (base, exponent), nothing)
Base.:^(base::Term, exponent::Coercible) = _apply("pow", (base, exponent), nothing)
Base.:^(base::Coercible, exponent::Term) = _apply("pow", (base, exponent), nothing)
Base.abs(term::Term) = _apply("abs", (term,), nothing)

# ---------------------------------------------------------------------------
# answers
# ---------------------------------------------------------------------------

"""
    Answer

The result of a computation: the reduced term plus everything the engine
reported about it.

Fields: `term`, `text`, `value`, `error`, `reduced`, `iterations`, `dimension`,
`diagnostic_codes` (the stable codes of the advisory dimensional findings).
"""
struct Answer
    term::Term
    text::String
    value::Union{Float64, Nothing}
    error::Union{Float64, Nothing}
    reduced::Bool
    iterations::Int
    dimension::Union{String, Nothing}
    diagnostic_codes::Vector{String}
end

function Base.show(io::IO, answer::Answer)
    print(io, "Answer(", repr(answer.text))
    answer.value === nothing || print(io, ", value = ", answer.value)
    answer.reduced || print(io, ", reduced = false")
    answer.dimension === nothing || print(io, ", dimension = ", repr(answer.dimension))
    isempty(answer.diagnostic_codes) || print(io, ", codes = ", answer.diagnostic_codes)
    print(io, ")")
end

"""`answer.latex`, `answer.typst`, `answer.pretty` and `answer.tensor`."""
function Base.getproperty(answer::Answer, name::Symbol)
    name === :latex && return latex(getfield(answer, :term))
    name === :typst && return typst(getfield(answer, :term))
    name === :pretty && return pretty(getfield(answer, :term))
    name === :tensor && return tensor(getfield(answer, :term))
    return getfield(answer, name)
end

Base.:(==)(left::Answer, right::Answer) = left.text == right.text
Base.hash(answer::Answer, seed::UInt) = hash(answer.text, seed)

"""
    compute(term, config = nothing) -> Answer

Reduce `term` to the target of `config` (symbolic by default).
"""
function compute(term::Term, config::Union{Config, Nothing} = nothing)
    return first(compute_all(term.session, (term,); config = config))::Answer
end

compute(session::Session, request, config::Union{Config, Nothing} = nothing) =
    compute(term(session, request), config)

"""
    simplify(session, source) -> Answer

Parse and reduce `source` with the default symbolic configuration.
"""
function simplify(session::Session, source::AbstractString)
    handle = simplify(session.ptr, source)
    handle == NO_TERM && throw(RssnError(Int32(3), last_error()))
    return compute(Term(session, handle))
end

"""
    compute_all(session, requests; config = nothing, errors = :raise)

Compute several requests in parallel inside the library.

`errors` decides what happens to a request that failed: `:raise` raises the
first failure, `:collect` puts the [`RssnError`](@ref) in its place, and
`:ignore` puts `nothing` there.
"""
function compute_all(session::Session, requests;
                     config::Union{Config, Nothing} = nothing,
                     errors::Symbol = :raise)
    errors in (:raise, :collect, :ignore) ||
        throw(ArgumentError("errors must be :raise, :collect, or :ignore"))
    terms = Term[term(session, request) for request in requests]
    isempty(terms) && return Any[]
    handles = UInt32[entry.handle for entry in terms]
    answers = Vector{CAnswer}(undef, length(terms))
    statuses = Vector{Cint}(undef, length(terms))
    symbol = _batch_metadata_symbol()
    metadata = symbol == C_NULL ? nothing : fill(_EMPTY_ANSWER_METADATA, length(terms))
    config_ptr = config === nothing ? C_NULL : config.ptr
    out = Vector{Any}(undef, length(terms))
    first_error = nothing
    try
        if metadata === nothing
            _check(compute_all(session.ptr, handles, config_ptr, answers, statuses))
        else
            _check(compute_all_with_metadata(symbol, session.ptr, handles, config_ptr,
                                             answers, statuses, metadata))
        end
        batch_error = last_error()
        for index in eachindex(terms)
            status = Int32(statuses[index])
            dim, codes, message = metadata === nothing ?
                                  ("", "", batch_error) : _read_metadata!(metadata, index)
            if status == 0
                if metadata === nothing
                    dim = dimension(terms[index], config)
                    diagnostic_codes = _dimension_codes()
                else
                    diagnostic_codes = String.(split(codes, '\n'; keepempty = false))
                end
                out[index] = _answer(session, answers[index], dim, diagnostic_codes)
            else
                error = RssnError(status, message)
                first_error === nothing && (first_error = error)
                out[index] = errors === :collect ? error : nothing
            end
        end
    finally
        metadata === nothing || foreach(_free_metadata, metadata)
    end
    errors === :raise && first_error !== nothing && throw(first_error)
    return out
end

"""The compute record (`RssnComputeRecord`) as a JSON string."""
function compute_json(term::Term, config::Union{Config, Nothing} = nothing)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(compute_json(term.session.ptr, term.handle,
                        config === nothing ? C_NULL : config.ptr, out))
    return take_string(out[])
end

"""The compute record in the deterministic `bincode-next` wire format."""
function compute_bincode(term::Term, config::Union{Config, Nothing} = nothing)
    out = Ref{Ptr{UInt8}}(C_NULL)
    length = Ref{Csize_t}(0)
    _check(compute_bincode(term.session.ptr, term.handle,
                           config === nothing ? C_NULL : config.ptr, out, length))
    return take_bytes(out[], length[])
end

"""The batch record (`RssnComputeBatchRecord`) in the `bincode-next` format."""
function compute_all_bincode(session::Session, requests,
                             config::Union{Config, Nothing} = nothing)
    handles = UInt32[term(session, request).handle for request in requests]
    out = Ref{Ptr{UInt8}}(C_NULL)
    length = Ref{Csize_t}(0)
    _check(compute_all_bincode(session.ptr, handles,
                               config === nothing ? C_NULL : config.ptr, out, length))
    return take_bytes(out[], length[])
end

"""
    eval_term(term, bindings = Dict()) -> Float64
    eval_term(term; x = 0.5, ...) -> Float64

Evaluate a term numerically without running any rule.

```julia
eval_term(parse_term(session, "sin(x) + 1"), Dict("x" => 0.0))   # 1.0
eval_term(parse_term(session, "sin(x) + 1"); x = 0.0)            # 1.0
```
"""
function eval_term(term::Term; bindings...)
    table = Dict{String, Float64}(string(key) => Float64(value)
                                  for (key, value) in pairs(bindings))
    return eval_term(term, table)
end

function eval_term(term::Term,
                   bindings::AbstractDict{<:AbstractString, <:Real})
    names = String[String(name) for name in keys(bindings)]
    values = Float64[Float64(bindings[name]) for name in keys(bindings)]
    out = Ref{Cdouble}(0.0)
    _check(term_eval(term.session.ptr, term.handle, names, values, out))
    return out[]
end

"""The inferred dimension of a term as text (`""` when none)."""
function dimension(term::Term, config::Union{Config, Nothing} = nothing)
    out = Ref{Ptr{UInt8}}(C_NULL)
    _check(dimension(term.session.ptr, term.handle,
                     config === nothing ? C_NULL : config.ptr, out))
    return take_string(out[])
end

function _dimension_codes()
    # `dimension_count()` is a `Csize_t`; convert before subtracting so an empty
    # report does not wrap around to a huge range.
    findings = Int(dimension_count())
    codes = Vector{String}(undef, findings)
    for index in 0:(findings - 1)
        pointer = dimension_code(index)
        codes[index + 1] = pointer == C_NULL ? "" : unsafe_string(pointer)
    end
    return codes
end

function _answer(session::Session, raw::CAnswer, dimension_text::String, codes::Vector{String})
    answer_term = Term(session, raw.term)
    return Answer(
        answer_term,
        text(answer_term),
        raw.has_value ? raw.value : nothing,
        raw.has_value ? raw.error : nothing,
        raw.reduced,
        Int(raw.iterations),
        isempty(dimension_text) ? nothing : dimension_text,
        codes,
    )
end

# ---------------------------------------------------------------------------
# a default session per task, for the operator facade
# ---------------------------------------------------------------------------

const _DEFAULT_KEY = :rssn_default_session

"""
    default_session() -> Session

The current task's session, created on first use.
"""
function default_session()
    storage = task_local_storage()
    session = get(storage, _DEFAULT_KEY, nothing)
    if session === nothing || session.ptr == C_NULL
        session = Session()
        storage[_DEFAULT_KEY] = session
    end
    return session
end

"""
    with_session(f, session)

Run `f()` with `session` as the default session:

```julia
with_session(session) do
    diff(sin(x) * x, x)
end
```
"""
function with_session(f, session::Session)
    storage = task_local_storage()
    previous = get(storage, _DEFAULT_KEY, nothing)
    storage[_DEFAULT_KEY] = session
    try
        return f()
    finally
        storage[_DEFAULT_KEY] = previous
    end
end

"""The standard kernel manifest as a JSON string."""
manifest_json_string() = take_string(manifest_json())
