# The C ABI of `rssn.h`, one `ccall` per entry point.

using Libdl
#
# Handles are `UInt32` terms in one session, failures are `RssnStatus` codes
# plus a thread-local message. `session.jl`, `term.jl` and `config.jl` build the
# ergonomic layer on top of this file.

# ---------------------------------------------------------------------------
# library resolution
# ---------------------------------------------------------------------------

"""The library file name for this platform."""
libname() = Sys.iswindows() ? "rssn.dll" :
            Sys.isapple() ? "librssn.dylib" : "librssn.so"

const _LIB = Ref{String}("")
const _HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const _SYMBOLS = Dict{Symbol, Ptr{Cvoid}}()
const _LIB_LOCK = ReentrantLock()

function _library_location(path::AbstractString)
    if isdir(path)
        candidate = joinpath(path, libname())
        isfile(candidate) || error("no rssn shared library in $(repr(path)); expected $(libname())")
        return abspath(candidate)
    end
    return isfile(path) ? abspath(path) : String(path)
end

function _candidate_paths()
    explicit = get(ENV, "RSSN_LIB", "")
    isempty(explicit) || return [_library_location(explicit)]

    here = @__DIR__
    root = normpath(joinpath(here, "..", "..", ".."))
    paths = [joinpath(here, "..", "lib", libname())]
    if here == joinpath(root, "bindings", "julia", "src") && isfile(joinpath(root, "rssn.h"))
        push!(paths, joinpath(root, libname()))
        if isfile(joinpath(root, "MODULE.bazel"))
            for target in ("ffi_distribution", "ffi", "ffi_distribution_portable",
                           "ffi_dev", "ffi_distribution_dev")
                push!(paths, joinpath(root, "bazel-bin", target, "ffi", libname()))
            end
        end
    end
    for profile in ("debug", "ffi", "release")
        for directory in ("target", "target-bindings", "target-ffi", joinpath("dist", "rssn"))
            push!(paths, joinpath(root, directory, profile, libname()))
        end
    end
    filter!(isfile, paths)
    push!(paths, libname())
    return paths
end

"""
    libpath() -> String

The path or platform loader name of the bound RSSN shared library. On first
use, resolves `RSSN_LIB`, the package or extracted distribution, a source
checkout's build directories, then the platform loader path.
"""
function libpath()
    lock(_LIB_LOCK) do
        if isempty(_LIB[])
            requested = get(ENV, "RSSN_LIB", "")
            for path in _candidate_paths()
                handle = Libdl.dlopen_e(path, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
                handle == C_NULL && continue
                if Libdl.dlsym_e(handle, :rssn_version) == C_NULL
                    Libdl.dlclose(handle)
                    isempty(requested) || error("the library at $(repr(path)) does not expose the RSSN C ABI")
                    continue
                end
                _HANDLE[] = handle
                _LIB[] = path
                break
            end
            if isempty(_LIB[])
                isempty(requested) || error("cannot load the RSSN shared library at $(repr(requested))")
                error("cannot load the RSSN shared library. Set RSSN_LIB to the library file or its " *
                      "directory, or build it with `bazel build //:ffi_dev`.")
            end
        end
        return _LIB[]
    end
end

"""
    set_library!(path)

Bind RSSN to a shared-library file or directory before creating sessions.
An empty path restores automatic discovery.
"""
function set_library!(path::AbstractString)
    lock(_LIB_LOCK) do
        resolved = isempty(path) ? "" : _library_location(path)
        handle = isempty(resolved) ? C_NULL :
                 Libdl.dlopen(resolved, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)
        if handle != C_NULL && Libdl.dlsym_e(handle, :rssn_version) == C_NULL
            Libdl.dlclose(handle)
            error("the library at $(repr(resolved)) does not expose the RSSN C ABI")
        end
        if _HANDLE[] != C_NULL
            Libdl.dlclose(_HANDLE[])
        end
        _HANDLE[] = handle
        _LIB[] = resolved
        empty!(_SYMBOLS)
        return resolved
    end
end

# Sessions can live on different host threads. Protect the shared symbol cache
# while discovering the library and on its first use from each entry point.
function _sym(name::Symbol)
    lock(_LIB_LOCK) do
        get!(_SYMBOLS, name) do
            libpath()
            Libdl.dlsym(_HANDLE[], name)
        end
    end
end

function _optional_sym(name::Symbol)
    lock(_LIB_LOCK) do
        get!(_SYMBOLS, name) do
            libpath()
            Libdl.dlsym_e(_HANDLE[], name)
        end
    end
end

_batch_metadata_symbol() = _optional_sym(:rssn_compute_all_with_metadata)

# ---------------------------------------------------------------------------
# C types
# ---------------------------------------------------------------------------

"""The handle value that means "no term" (a failed constructor)."""
const NO_TERM = typemax(UInt32)

"""`struct RssnAnswer`: the result of `rssn_compute`."""
struct CAnswer
    term::UInt32
    has_value::Bool
    value::Cdouble
    error::Cdouble
    reduced::Bool
    iterations::Csize_t
end

"""Owned metadata strings from `rssn_compute_all_with_metadata`."""
struct CAnswerMetadata
    dimension::Ptr{UInt8}
    diagnostic_codes::Ptr{UInt8}
    message::Ptr{UInt8}
end

const _EMPTY_ANSWER_METADATA = CAnswerMetadata(C_NULL, C_NULL, C_NULL)

"""`struct RssnTensor`: a dense row-major `f64` array extracted from a term."""
struct CTensor
    rank::Csize_t
    shape::Ptr{Csize_t}
    len::Csize_t
    data::Ptr{Cdouble}
end

# ---------------------------------------------------------------------------
# ownership helpers
# ---------------------------------------------------------------------------

"""Read a caller-owned `char *` and release it with `rssn_string_free`."""
function take_string(ptr::Ptr{UInt8})
    ptr == C_NULL && return ""
    text = unsafe_string(ptr)
    ccall(_sym(:rssn_string_free), Cvoid, (Ptr{UInt8},), ptr)
    return text
end

function _free_metadata(metadata::CAnswerMetadata)
    for pointer in (metadata.dimension, metadata.diagnostic_codes, metadata.message)
        pointer == C_NULL && continue
        ccall(_sym(:rssn_string_free), Cvoid, (Ptr{UInt8},), pointer)
    end
end

function _read_metadata!(rows, index)
    metadata = rows[index]
    rows[index] = _EMPTY_ANSWER_METADATA
    try
        read_string(p) = p == C_NULL ? "" : unsafe_string(p)
        return (read_string(metadata.dimension), read_string(metadata.diagnostic_codes),
                read_string(metadata.message))
    finally
        _free_metadata(metadata)
    end
end

"""Read a caller-owned byte buffer and release it with `rssn_bytes_free`."""
function take_bytes(ptr::Ptr{UInt8}, length::Integer)
    ptr == C_NULL && return UInt8[]
    bytes = Vector{UInt8}(undef, length)
    unsafe_copyto!(ptr, Base.pointer(bytes), length)
    ccall(_sym(:rssn_bytes_free), Cvoid, (Ptr{UInt8}, Csize_t), ptr, length)
    return bytes
end

"""The message of the last failed call on this thread (may be empty)."""
function last_error()
    message = ccall(_sym(:rssn_last_error), Ptr{UInt8}, ())
    return message == C_NULL ? "" : unsafe_string(message)
end

# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------

version() = unsafe_string(ccall(_sym(:rssn_version), Ptr{UInt8}, ()))

session_new() = ccall(_sym(:rssn_session_new), Ptr{Cvoid}, ())
session_free(ptr::Ptr{Cvoid}) =
    ccall(_sym(:rssn_session_free), Cvoid, (Ptr{Cvoid},), ptr)

config_new() = ccall(_sym(:rssn_config_new), Ptr{Cvoid}, ())
config_free(ptr::Ptr{Cvoid}) =
    ccall(_sym(:rssn_config_free), Cvoid, (Ptr{Cvoid},), ptr)
config_symbolic(ptr::Ptr{Cvoid}) =
    ccall(_sym(:rssn_config_symbolic), Cint, (Ptr{Cvoid},), ptr)
config_numeric(ptr::Ptr{Cvoid}, tolerance::Real) =
    ccall(_sym(:rssn_config_numeric), Cint, (Ptr{Cvoid}, Cdouble), ptr, tolerance)
config_bind(ptr::Ptr{Cvoid}, name::AbstractString, value::Real) =
    ccall(_sym(:rssn_config_bind), Cint, (Ptr{Cvoid}, Cstring, Cdouble),
          ptr, name, value)
config_assume(ptr::Ptr{Cvoid}, name::AbstractString, facts::Integer) =
    ccall(_sym(:rssn_config_assume), Cint, (Ptr{Cvoid}, Cstring, UInt8),
          ptr, name, facts)
config_budget(ptr::Ptr{Cvoid}, iterations::Integer, nodes::Integer) =
    ccall(_sym(:rssn_config_budget), Cint, (Ptr{Cvoid}, Csize_t, Csize_t),
          ptr, iterations, nodes)
config_dimension_mode(ptr::Ptr{Cvoid}, mode::Integer) =
    ccall(_sym(:rssn_config_dimension_mode), Cint, (Ptr{Cvoid}, UInt8), ptr, mode)
config_bypass_dimensions(ptr::Ptr{Cvoid}) =
    ccall(_sym(:rssn_config_bypass_dimensions), Cint, (Ptr{Cvoid},), ptr)
config_dimension(ptr::Ptr{Cvoid}, symbol::AbstractString, text::AbstractString) =
    ccall(_sym(:rssn_config_dimension), Cint, (Ptr{Cvoid}, Cstring, Cstring),
          ptr, symbol, text)
config_define_base_dimension(ptr::Ptr{Cvoid}, name::AbstractString) =
    ccall(_sym(:rssn_config_define_base_dimension), Cint, (Ptr{Cvoid}, Cstring),
          ptr, name)

term_parse(ptr::Ptr{Cvoid}, source::AbstractString) =
    ccall(_sym(:rssn_term_parse), UInt32, (Ptr{Cvoid}, Cstring), ptr, source)
term_sym(ptr::Ptr{Cvoid}, name::AbstractString) =
    ccall(_sym(:rssn_term_sym), UInt32, (Ptr{Cvoid}, Cstring), ptr, name)
term_int(ptr::Ptr{Cvoid}, value::Integer) =
    ccall(_sym(:rssn_term_int), UInt32, (Ptr{Cvoid}, Int64), ptr, value)
term_float(ptr::Ptr{Cvoid}, value::Real) =
    ccall(_sym(:rssn_term_float), UInt32, (Ptr{Cvoid}, Cdouble), ptr, value)
term_rational(ptr::Ptr{Cvoid}, numerator::Integer, denominator::Integer) =
    ccall(_sym(:rssn_term_rational), UInt32, (Ptr{Cvoid}, Int64, Int64),
          ptr, numerator, denominator)
term_apply(ptr::Ptr{Cvoid}, op::AbstractString, args::Vector{UInt32}) =
    ccall(_sym(:rssn_term_apply), UInt32, (Ptr{Cvoid}, Cstring, Ptr{UInt32}, Csize_t),
          ptr, op, isempty(args) ? C_NULL : Base.pointer(args), length(args))
# `pointer()` on a `Vector{String}` points at the `String` objects, not at their
# bytes, so the `const char *const *` arguments go through `cstring_vector`: each
# element is converted to a `Cstring`, and the caller keeps the strings alive
# across the call with `GC.@preserve`.
function cstring_vector(strings::Vector{String})
    labels = Vector{Cstring}(undef, length(strings))
    for index in eachindex(strings)
        labels[index] = Base.unsafe_convert(Cstring, Base.cconvert(Cstring, strings[index]))
    end
    return labels
end

function term_eval(ptr::Ptr{Cvoid}, term::UInt32, names::Vector{String},
                   values::Vector{Float64}, out::Base.RefValue{Cdouble})
    GC.@preserve names begin
        labels = cstring_vector(names)
        return ccall(_sym(:rssn_term_eval), Cint,
                     (Ptr{Cvoid}, UInt32, Ptr{Cstring}, Ptr{Cdouble}, Csize_t,
                      Ptr{Cdouble}),
                     ptr, term, isempty(labels) ? C_NULL : Base.pointer(labels),
                     isempty(values) ? C_NULL : Base.pointer(values), length(labels), out)
    end
end
term_to_f64(ptr::Ptr{Cvoid}, term::UInt32, out::Base.RefValue{Cdouble}) =
    ccall(_sym(:rssn_term_to_f64), Cint, (Ptr{Cvoid}, UInt32, Ptr{Cdouble}),
          ptr, term, out)
term_tensor(ptr::Ptr{Cvoid}, term::UInt32, out::Base.RefValue{CTensor}) =
    ccall(_sym(:rssn_term_tensor_data), Cint, (Ptr{Cvoid}, UInt32, Ptr{CTensor}),
          ptr, term, out)
tensor_free(tensor::Base.RefValue{CTensor}) =
    ccall(_sym(:rssn_tensor_free), Cvoid, (Ptr{CTensor},), tensor)
term_to_string(ptr::Ptr{Cvoid}, term::UInt32) =
    ccall(_sym(:rssn_term_to_string), Ptr{UInt8}, (Ptr{Cvoid}, UInt32), ptr, term)
term_to_latex(ptr::Ptr{Cvoid}, term::UInt32) =
    ccall(_sym(:rssn_term_to_latex), Ptr{UInt8}, (Ptr{Cvoid}, UInt32), ptr, term)
term_to_typst(ptr::Ptr{Cvoid}, term::UInt32) =
    ccall(_sym(:rssn_term_to_typst), Ptr{UInt8}, (Ptr{Cvoid}, UInt32), ptr, term)
term_to_pretty(ptr::Ptr{Cvoid}, term::UInt32) =
    ccall(_sym(:rssn_term_to_pretty), Ptr{UInt8}, (Ptr{Cvoid}, UInt32), ptr, term)

compute(ptr::Ptr{Cvoid}, term::UInt32, config::Ptr{Cvoid}, out::Base.RefValue{CAnswer}) =
    ccall(_sym(:rssn_compute), Cint,
          (Ptr{Cvoid}, UInt32, Ptr{Cvoid}, Ptr{CAnswer}), ptr, term, config, out)
compute_all(ptr::Ptr{Cvoid}, terms::Vector{UInt32}, config::Ptr{Cvoid},
            answers::Vector{CAnswer}, statuses::Vector{Cint}) =
    ccall(_sym(:rssn_compute_all), Cint,
          (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cvoid}, Ptr{CAnswer}, Ptr{Cint}),
          ptr, Base.pointer(terms), length(terms), config, Base.pointer(answers), Base.pointer(statuses))
compute_all_with_metadata(symbol::Ptr{Cvoid}, ptr::Ptr{Cvoid}, terms::Vector{UInt32},
                          config::Ptr{Cvoid}, answers::Vector{CAnswer},
                          statuses::Vector{Cint}, metadata::Vector{CAnswerMetadata}) =
    ccall(symbol, Cint,
          (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cvoid}, Ptr{CAnswer}, Ptr{Cint}, Ptr{CAnswerMetadata}),
          ptr, terms, length(terms), config, answers, statuses, metadata)
compute_json(ptr::Ptr{Cvoid}, term::UInt32, config::Ptr{Cvoid}, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_compute_json), Cint,
          (Ptr{Cvoid}, UInt32, Ptr{Cvoid}, Ptr{Ptr{UInt8}}), ptr, term, config, out)
compute_bincode(ptr::Ptr{Cvoid}, term::UInt32, config::Ptr{Cvoid},
                out::Base.RefValue{Ptr{UInt8}}, length::Base.RefValue{Csize_t}) =
    ccall(_sym(:rssn_compute_bincode), Cint,
          (Ptr{Cvoid}, UInt32, Ptr{Cvoid}, Ptr{Ptr{UInt8}}, Ptr{Csize_t}),
          ptr, term, config, out, length)
compute_all_bincode(ptr::Ptr{Cvoid}, terms::Vector{UInt32}, config::Ptr{Cvoid},
                    out::Base.RefValue{Ptr{UInt8}}, length::Base.RefValue{Csize_t}) =
    ccall(_sym(:rssn_compute_all_bincode), Cint,
          (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cvoid}, Ptr{Ptr{UInt8}}, Ptr{Csize_t}),
          ptr, Base.pointer(terms), length(terms), config, out, length)
simplify(ptr::Ptr{Cvoid}, source::AbstractString) =
    ccall(_sym(:rssn_simplify), UInt32, (Ptr{Cvoid}, Cstring), ptr, source)
dimension(ptr::Ptr{Cvoid}, term::UInt32, config::Ptr{Cvoid}, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_dimension), Cint,
          (Ptr{Cvoid}, UInt32, Ptr{Cvoid}, Ptr{Ptr{UInt8}}), ptr, term, config, out)
dimension_count() = ccall(_sym(:rssn_dimension_count), Csize_t, ())
dimension_code(index::Integer) =
    ccall(_sym(:rssn_dimension_code), Ptr{UInt8}, (Csize_t,), index)

manifest_json() = ccall(_sym(:rssn_manifest_json), Ptr{UInt8}, ())
manifest_hash() = ccall(_sym(:rssn_manifest_hash), UInt64, ())
op_info_json(name::AbstractString, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_op_info_json), Cint, (Cstring, Ptr{Ptr{UInt8}}), name, out)
sim_run(name::AbstractString, params::AbstractString, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_sim_run), Cint, (Cstring, Cstring, Ptr{Ptr{UInt8}}),
          name, params, out)

# ---------------------------------------------------------------------------
# dynamic values (the value ABI)
# ---------------------------------------------------------------------------

value_new_int(ptr::Ptr{Cvoid}, value::Integer, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_int), Cint, (Ptr{Cvoid}, Int64, Ptr{UInt32}), ptr, value, out)
function value_new_int_decimal(ptr::Ptr{Cvoid}, digits::AbstractString, out::Base.RefValue{UInt32})
    entry = _optional_sym(:rssn_value_new_int_decimal)
    entry == C_NULL && throw(ArgumentError("this RSSN library only supports signed 64-bit integer cells"))
    return ccall(entry, Cint, (Ptr{Cvoid}, Cstring, Ptr{UInt32}), ptr, digits, out)
end
function value_with_coherent_dimension(ptr::Ptr{Cvoid}, cell::UInt32, payload::String,
                                       out::Base.RefValue{UInt32})
    entry = _optional_sym(:rssn_value_with_coherent_dimension)
    entry == C_NULL && throw(ArgumentError("this RSSN library lacks coherent dimension cells"))
    GC.@preserve payload return ccall(entry, Cint,
        (Ptr{Cvoid}, UInt32, Ptr{UInt8}, Csize_t, Ptr{UInt32}),
        ptr, cell, pointer(payload), ncodeunits(payload), out)
end
value_new_float(ptr::Ptr{Cvoid}, value::Real, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_float), Cint, (Ptr{Cvoid}, Cdouble, Ptr{UInt32}), ptr, value, out)
value_new_text(ptr::Ptr{Cvoid}, text::AbstractString, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_text), Cint, (Ptr{Cvoid}, Cstring, Ptr{UInt32}), ptr, text, out)
value_new_symbol(ptr::Ptr{Cvoid}, name::AbstractString, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_symbol), Cint, (Ptr{Cvoid}, Cstring, Ptr{UInt32}), ptr, name, out)
value_new_term(ptr::Ptr{Cvoid}, term::UInt32, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_term), Cint, (Ptr{Cvoid}, UInt32, Ptr{UInt32}), ptr, term, out)
value_new_list(ptr::Ptr{Cvoid}, cells::Vector{UInt32}, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_list), Cint, (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{UInt32}),
          ptr, isempty(cells) ? C_NULL : Base.pointer(cells), length(cells), out)
value_call(ptr::Ptr{Cvoid}, op::AbstractString, cells::Vector{UInt32},
           out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_call), Cint, (Ptr{Cvoid}, Cstring, Ptr{UInt32}, Csize_t, Ptr{UInt32}),
          ptr, op, isempty(cells) ? C_NULL : Base.pointer(cells), length(cells), out)
value_call_key(ptr::Ptr{Cvoid}, rule_set::AbstractString, name::AbstractString,
               hash::Vector{UInt8}, cells::Vector{UInt32}, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_call_key), Cint,
          (Ptr{Cvoid}, Cstring, Cstring, Ptr{UInt8}, Ptr{UInt32}, Csize_t, Ptr{UInt32}),
          ptr, rule_set, name, hash, cells, length(cells), out)
value_kind(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{UInt8}) =
    ccall(_sym(:rssn_value_kind), Cint, (Ptr{Cvoid}, UInt32, Ptr{UInt8}), ptr, cell, out)
value_to_f64(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{Cdouble}) =
    ccall(_sym(:rssn_value_to_f64), Cint, (Ptr{Cvoid}, UInt32, Ptr{Cdouble}), ptr, cell, out)
value_to_text(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_value_to_text), Cint, (Ptr{Cvoid}, UInt32, Ptr{Ptr{UInt8}}), ptr, cell, out)
value_to_term(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_to_term), Cint, (Ptr{Cvoid}, UInt32, Ptr{UInt32}), ptr, cell, out)
value_to_tensor(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{CTensor}) =
    ccall(_sym(:rssn_value_to_tensor), Cint, (Ptr{Cvoid}, UInt32, Ptr{CTensor}), ptr, cell, out)
value_len(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{Csize_t}) =
    ccall(_sym(:rssn_value_len), Cint, (Ptr{Cvoid}, UInt32, Ptr{Csize_t}), ptr, cell, out)
value_at(ptr::Ptr{Cvoid}, cell::UInt32, index::Integer, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_at), Cint, (Ptr{Cvoid}, UInt32, Csize_t, Ptr{UInt32}),
          ptr, cell, index, out)
value_json(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_value_json), Cint, (Ptr{Cvoid}, UInt32, Ptr{Ptr{UInt8}}), ptr, cell, out)
value_free(ptr::Ptr{Cvoid}, cell::UInt32) =
    ccall(_sym(:rssn_value_free), Cvoid, (Ptr{Cvoid}, UInt32), ptr, cell)
value_count(ptr::Ptr{Cvoid}) = ccall(_sym(:rssn_value_count), Csize_t, (Ptr{Cvoid},), ptr)

# ---------------------------------------------------------------------------
# host callbacks (the value ABI's function arguments)
# ---------------------------------------------------------------------------

"""The C signature of a host callback: `(user, session, args, count, out)`."""
const CALLBACK_SIGNATURE = (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{UInt32})

value_new_bool(ptr::Ptr{Cvoid}, value::Bool, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_bool), Cint, (Ptr{Cvoid}, Cint, Ptr{UInt32}),
          ptr, value ? 1 : 0, out)
value_new_complex(ptr::Ptr{Cvoid}, re::Real, im::Real, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_complex), Cint, (Ptr{Cvoid}, Cdouble, Cdouble, Ptr{UInt32}),
          ptr, re, im, out)
function value_new_record(ptr::Ptr{Cvoid}, names::Vector{String}, cells::Vector{UInt32},
                          out::Base.RefValue{UInt32})
    GC.@preserve names begin
        labels = cstring_vector(names)
        return ccall(_sym(:rssn_value_new_record), Cint,
                     (Ptr{Cvoid}, Ptr{Cstring}, Ptr{UInt32}, Csize_t, Ptr{UInt32}),
                     ptr, isempty(labels) ? C_NULL : Base.pointer(labels),
                     isempty(cells) ? C_NULL : Base.pointer(cells), length(labels), out)
    end
end
value_field_name(ptr::Ptr{Cvoid}, cell::UInt32, index::Integer, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_value_field_name), Cint, (Ptr{Cvoid}, UInt32, Csize_t, Ptr{Ptr{UInt8}}),
          ptr, cell, index, out)
value_tag(ptr::Ptr{Cvoid}, cell::UInt32, out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_value_tag), Cint, (Ptr{Cvoid}, UInt32, Ptr{Ptr{UInt8}}), ptr, cell, out)

callback_new(ptr::Ptr{Cvoid}, signature::AbstractString, function_::Ptr{Cvoid},
             user::Ptr{Cvoid}, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_callback_new), Cint,
          (Ptr{Cvoid}, Cstring, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt32}),
          ptr, signature, function_, user, out)
callback_new_for(ptr::Ptr{Cvoid}, op::AbstractString, index::Integer, function_::Ptr{Cvoid},
                 user::Ptr{Cvoid}, out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_callback_new_for), Cint,
          (Ptr{Cvoid}, Cstring, Csize_t, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{UInt32}),
          ptr, op, index, function_, user, out)
callback_signatures(ptr::Ptr{Cvoid}, op::AbstractString, index::Integer,
                    out::Base.RefValue{Ptr{UInt8}}) =
    ccall(_sym(:rssn_callback_signatures), Cint,
          (Ptr{Cvoid}, Cstring, Csize_t, Ptr{Ptr{UInt8}}), ptr, op, index, out)
callback_free(ptr::Ptr{Cvoid}, cell::UInt32) =
    ccall(_sym(:rssn_callback_free), Cvoid, (Ptr{Cvoid}, UInt32), ptr, cell)
callback_count(ptr::Ptr{Cvoid}) = ccall(_sym(:rssn_callback_count), Csize_t, (Ptr{Cvoid},), ptr)
callback_fail(message::AbstractString) =
    ccall(_sym(:rssn_callback_fail), Cvoid, (Cstring,), message)

function value_new_tensor_cells(ptr::Ptr{Cvoid}, cells::Vector{UInt32},
                                shape::Vector{Csize_t}, exemplar::UInt32,
                                out::Base.RefValue{UInt32})
    entry = _optional_sym(:rssn_value_new_tensor_cells)
    entry == C_NULL && throw(ArgumentError("this RSSN library lacks tensor cell construction"))
    return ccall(entry, Cint,
                 (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Csize_t}, Csize_t, UInt32, Ptr{UInt32}),
                 ptr, cells, length(cells), shape, length(shape), exemplar, out)
end

value_new_tensor(ptr::Ptr{Cvoid}, data::Vector{Cdouble}, shape::Vector{Csize_t},
                 out::Base.RefValue{UInt32}) =
    ccall(_sym(:rssn_value_new_tensor), Cint,
          (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t, Ptr{Csize_t}, Csize_t, Ptr{UInt32}),
          ptr, data, length(data), shape, length(shape), out)
