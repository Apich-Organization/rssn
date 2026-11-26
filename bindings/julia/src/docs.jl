# The operator reference: the engine's own documentation, searchable.
#
# Every operator's documentation travels with the package: it is generated from
# the kernel manifest into the docstring above each function, so `?diff` and
# `Docs.doc` already show it. This file adds the two things a REPL user wants on
# top of that -- one operator by name, and a search over all of them.
#
#     julia> print(op_help("diff"))
#     julia> op_search("brent")            # [(name, first line of the doc), ...]
#     julia> op_search(category = "calculus")[1:5]

"""
    op_help(name) -> String

The documentation of one operator, as text: the binding, the arity, the
category and the engine's own description.

```julia
print(op_help("diff"))
print(op_help("num_brent"))
```
"""
function op_help(name::AbstractString)
    haskey(OPS, name) || throw(ArgumentError("'$name' is not a registered operator"))
    binding = NAMES[name]
    minimum, maximum, _ = OPS[name]
    arity = maximum === nothing ? "$(minimum)+ arguments" :
            minimum == maximum ? (minimum == 1 ? "1 argument" : "$(minimum) arguments") :
            "$(minimum) to $(maximum) arguments"
    category = something(_category_of(binding), "uncategorised")
    lines = [
        "RSSN.$binding",
        "",
        "  operator   : $name",
        "  arity      : $arity",
        "  category   : $category",
        "",
        strip(get(OPERATOR_DOCS, name, "(no documentation)")),
    ]
    return join(lines, "\n")
end

# The category whose member list holds this binding name.
function _category_of(binding::AbstractString)
    for (category, members) in BY_CATEGORY
        binding in members && return category
    end
    return nothing
end

"""
    op_names() -> Vector{String}

Every registered operator name, sorted.
"""
op_names() = sort!(collect(keys(OPS)))

"""
    op_categories() -> Vector{String}

The manifest's categories, sorted.
"""
op_categories() = sort!(collect(keys(BY_CATEGORY)))

"""
    op_rule_sets() -> Vector{String}

The engine's rule sets, sorted.
"""
op_rule_sets() = sort!(collect(keys(BY_RULE_SET)))

"""
    op_search(pattern = nothing; category = nothing, rule_set = nothing, limit = 50)

Operators whose name or documentation matches `pattern`, a regular expression
matched case-insensitively.

Returns `(name, first line of the documentation)` pairs in name order.

```julia
op_search("brent")
op_search(category = "calculus")
op_search(rule_set = "numerical_matrix", limit = 10)
```
"""
function op_search(pattern = nothing; category = nothing, rule_set = nothing, limit = 50)
    names = sort!(collect(keys(OPS)))
    if category !== nothing
        members = Set(get(BY_CATEGORY, category, String[]))
        names = [name for name in names if NAMES[name] in members]
    end
    if rule_set !== nothing
        members = Set(get(BY_RULE_SET, rule_set, String[]))
        names = [name for name in names if NAMES[name] in members]
    end
    expression = pattern === nothing ? nothing : Regex(pattern, "i")
    results = Tuple{String, String}[]
    for name in names
        doc = strip(get(OPERATOR_DOCS, name, ""))
        summary = isempty(doc) ? "" : first(split(doc, '\n'))
        if length(summary) > 120
            summary = string(first(summary, 119), "…")
        end
        if expression === nothing || occursin(expression, name) || occursin(expression, doc)
            push!(results, (name, summary))
            length(results) >= limit && break
        end
    end
    return results
end
