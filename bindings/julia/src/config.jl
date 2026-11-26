# The computation configuration: target, bindings, assumptions, dimensions.

"""The fact bits `rssn_config_assume` accepts."""
module Assumption
const REAL = UInt8(1)
const NONZERO = UInt8(2)
const NONNEGATIVE = UInt8(4)
const NONPOSITIVE = UInt8(8)
const INTEGER = UInt8(16)
const POSITIVE = REAL | NONZERO | NONNEGATIVE   # 7
const NEGATIVE = REAL | NONZERO | NONPOSITIVE   # 11
end

"""The dimension checking mode."""
module DimensionMode
const OFF = UInt8(0)
const ADVISORY = UInt8(1)
const STRICT = UInt8(2)
end

"""
    Config()

An owned `RssnConfig`: the target of a computation, symbol bindings,
assumptions and the dimension mode.

Every setter mutates the configuration and returns it, so a target reads in one
expression:

```julia
cfg = Config()
cfg = numeric!(cfg, 1e-12)
cfg = bind!(cfg, "x", 0.5)
```

Release it with `close(cfg)`; a `do` block or a `finalizer` also releases it.
"""
mutable struct Config
    ptr::Ptr{Cvoid}

    function Config()
        pointer = config_new()
        pointer == C_NULL && throw(RssnError(Int32(7), last_error()))
        config = new(pointer)
        finalizer(finalize_config, config)
        return config
    end
end

function finalize_config(config::Config)
    if config.ptr != C_NULL
        config_free(config.ptr)
        config.ptr = C_NULL
    end
    return nothing
end

"""
    close(config)

Release the configuration. Safe to call twice.

This extends `Base.close`, so no import is needed.
"""
function Base.close(config::Config)
    finalize_config(config)
    return nothing
end

Base.show(io::IO, config::Config) =
    print(io, "Config(", config.ptr == C_NULL ? "closed" : string(config.ptr), ")")

"""Ask for a closed form (the default)."""
function symbolic!(config::Config)
    _check(config_symbolic(config.ptr))
    return config
end

"""Ask for a numeric answer within `tolerance`."""
function numeric!(config::Config, tolerance::Real)
    _check(config_numeric(config.ptr, tolerance))
    return config
end

"""Bind the symbol `name` to `value` for numeric evaluation."""
function bind!(config::Config, name::AbstractString, value::Real)
    _check(config_bind(config.ptr, name, value))
    return config
end

"""Declare facts about `name`: an `Assumption` bit set or a name such as `"positive"`."""
function assume!(config::Config, name::AbstractString, facts)
    _check(config_assume(config.ptr, name, _facts(facts)))
    return config
end

"""Limit the search; zero keeps the current value."""
function budget!(config::Config, iterations::Integer = 0, nodes::Integer = 0)
    _check(config_budget(config.ptr, iterations, nodes))
    return config
end

const _MODE_WORDS = Dict(
    "off" => DimensionMode.OFF,
    "advisory" => DimensionMode.ADVISORY,
    "strict" => DimensionMode.STRICT,
)

"""Set the dimension mode (`DimensionMode.OFF`, `.ADVISORY`, `.STRICT`, or a name)."""
function dimension_mode!(config::Config, mode)
    value = if mode isa AbstractString
        _MODE_WORDS[lowercase(mode)]
    else
        UInt8(mode)
    end
    _check(config_dimension_mode(config.ptr, value))
    return config
end

"""Disable hard dimensional errors even in strict mode."""
function bypass_dimensions!(config::Config)
    _check(config_bypass_dimensions(config.ptr))
    return config
end

"""Declare the dimension of `symbol` as text (`"m/s"`, `"Length/Time"`, `"any"`)."""
function dimension!(config::Config, symbol::AbstractString, text::AbstractString)
    _check(config_dimension(config.ptr, symbol, text))
    return config
end

"""Register a base dimension the engine does not have (for example `"Currency"`)."""
function define_base_dimension!(config::Config, name::AbstractString)
    _check(config_define_base_dimension(config.ptr, name))
    return config
end

const _ASSUMPTION_WORDS = Dict(    "real" => Assumption.REAL,
    "nonzero" => Assumption.NONZERO,
    "nonnegative" => Assumption.NONNEGATIVE,
    "nonpositive" => Assumption.NONPOSITIVE,
    "integer" => Assumption.INTEGER,
    "positive" => Assumption.POSITIVE,
    "negative" => Assumption.NEGATIVE,
)

function _facts(value)
    value isa Integer && return UInt8(value)
    value isa AbstractString && return _ASSUMPTION_WORDS[lowercase(value)]
    total = UInt8(0)
    for item in value
        total |= item isa AbstractString ? _ASSUMPTION_WORDS[lowercase(item)] : UInt8(item)
    end
    return total
end
