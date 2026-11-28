# A minimal JSON codec, so the package has no dependencies.
#
# It covers what the simulation entry point needs: objects, arrays, strings,
# numbers, booleans and null. Encoding accepts `AbstractDict`, `AbstractVector`,
# `Tuple`, `AbstractString`, `Integer`, `AbstractFloat`, `Bool` and `Nothing`;
# decoding returns `Dict{String, Any}`, `Vector{Any}`, `String`, `Int64`,
# `Float64`, `Bool` and `nothing`.

const _JSON_ESCAPES = Dict(
    '"' => "\\\"", '\\' => "\\\\", '\b' => "\\b", '\f' => "\\f",
    '\n' => "\\n", '\r' => "\\r", '\t' => "\\t",
)

function _json_escape(text::AbstractString)
    out = IOBuffer()
    for character in text
        escape = get(_JSON_ESCAPES, character, nothing)
        if escape !== nothing
            print(out, escape)
        elseif character < ' '
            print(out, "\\u", lpad(string(UInt16(character); base = 16), 4, '0'))
        else
            print(out, character)
        end
    end
    return String(take!(out))
end

"""Encode `value` as JSON text."""
function json_encode(value)
    out = IOBuffer()
    _json_write(out, value)
    return String(take!(out))
end

function _json_write(out::IO, value)
    if value === nothing
        print(out, "null")
    elseif value isa Bool
        print(out, value ? "true" : "false")
    elseif value isa Integer
        print(out, value)
    elseif value isa AbstractFloat
        isfinite(value) || throw(ArgumentError("JSON cannot represent $value"))
        print(out, repr(Float64(value)))
    elseif value isa AbstractString
        print(out, '"', _json_escape(value), '"')
    elseif value isa AbstractDict
        print(out, '{')
        first = true
        for (key, entry) in value
            first || print(out, ',')
            first = false
            print(out, '"', _json_escape(string(key)), "\":")
            _json_write(out, entry)
        end
        print(out, '}')
    elseif value isa Union{AbstractVector, Tuple}
        print(out, '[')
        first = true
        for entry in value
            first || print(out, ',')
            first = false
            _json_write(out, entry)
        end
        print(out, ']')
    else
        throw(ArgumentError("cannot encode $(typeof(value)) as JSON"))
    end
    return nothing
end

"""A cursor over JSON text."""
mutable struct _JsonState
    text::String
    index::Int
end

"""Decode JSON text into Julia values."""
function json_decode(text::AbstractString)
    state = _JsonState(String(text), 1)
    value = _json_parse(state)
    _json_skip!(state)
    state.index <= lastindex(state.text) &&
        throw(ArgumentError("trailing data at byte $(state.index)"))
    return value
end

function _json_skip!(state)
    while state.index <= lastindex(state.text) && isspace(state.text[state.index])
        state.index = nextind(state.text, state.index)
    end
    return nothing
end

function _json_parse(state)
    _json_skip!(state)
    state.index > lastindex(state.text) && throw(ArgumentError("unexpected end of JSON"))
    character = state.text[state.index]
    if character == '{'
        return _json_object(state)
    elseif character == '['
        return _json_array(state)
    elseif character == '"'
        return _json_string(state)
    elseif startswith(SubString(state.text, state.index), "true")
        state.index += 4
        return true
    elseif startswith(SubString(state.text, state.index), "false")
        state.index += 5
        return false
    elseif startswith(SubString(state.text, state.index), "null")
        state.index += 4
        return nothing
    else
        return _json_number(state)
    end
end

function _json_expect!(state, character)
    _json_skip!(state)
    state.index <= lastindex(state.text) && state.text[state.index] == character ||
        throw(ArgumentError("expected '$character' at byte $(state.index)"))
    state.index = nextind(state.text, state.index)
    return nothing
end

function _json_object(state)
    dict = Dict{String, Any}()
    _json_expect!(state, '{')
    _json_skip!(state)
    if state.index <= lastindex(state.text) && state.text[state.index] == '}'
        state.index = nextind(state.text, state.index)
        return dict
    end
    while true
        _json_skip!(state)
        key = _json_string(state)
        _json_expect!(state, ':')
        dict[key] = _json_parse(state)
        _json_skip!(state)
        if state.index <= lastindex(state.text) && state.text[state.index] == ','
            state.index = nextind(state.text, state.index)
            continue
        end
        _json_expect!(state, '}')
        return dict
    end
end

function _json_array(state)
    entries = Any[]
    _json_expect!(state, '[')
    _json_skip!(state)
    if state.index <= lastindex(state.text) && state.text[state.index] == ']'
        state.index = nextind(state.text, state.index)
        return entries
    end
    while true
        push!(entries, _json_parse(state))
        _json_skip!(state)
        if state.index <= lastindex(state.text) && state.text[state.index] == ','
            state.index = nextind(state.text, state.index)
            continue
        end
        _json_expect!(state, ']')
        return entries
    end
end

function _json_string(state)
    _json_expect!(state, '"')
    out = IOBuffer()
    while true
        state.index > lastindex(state.text) && throw(ArgumentError("unterminated string"))
        character = state.text[state.index]
        state.index = nextind(state.text, state.index)
        if character == '"'
            return String(take!(out))
        elseif character == '\\'
            escape = state.text[state.index]
            state.index = nextind(state.text, state.index)
            if escape == 'u'
                code = parse(UInt16, state.text[state.index:(state.index + 3)]; base = 16)
                state.index += 4
                if 0xD800 <= code <= 0xDBFF
                    state.text[state.index] == '\\' || throw(ArgumentError("lone surrogate"))
                    state.index += 2
                    low = parse(UInt16, state.text[state.index:(state.index + 3)]; base = 16)
                    state.index += 4
                    print(out, Char(0x10000 + (UInt32(code - 0xD800) << 10) + UInt32(low - 0xDC00)))
                else
                    print(out, Char(code))
                end
            else
                print(out, get(Dict(
                    '"' => '"', '\\' => '\\', '/' => '/', 'b' => '\b', 'f' => '\f',
                    'n' => '\n', 'r' => '\r', 't' => '\t',
                ), escape) do
                    throw(ArgumentError("unknown escape \\\\$escape"))
                end)
            end
        else
            print(out, character)
        end
    end
end

function _json_number(state)
    start = state.index
    while state.index <= lastindex(state.text)
        character = state.text[state.index]
        if isdigit(character) || character in ('-', '+', '.', 'e', 'E')
            state.index = nextind(state.text, state.index)
        else
            break
        end
    end
    literal = state.text[start:(state.index - 1)]
    isempty(literal) && throw(ArgumentError("expected a value at byte $start"))
    if occursin('.', literal) || occursin('e', literal) || occursin('E', literal)
        return parse(Float64, literal)
    end
    return parse(Int64, literal)
end
