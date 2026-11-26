# Tests for the Julia binding.
#
#   julia --project=bindings/julia bindings/julia/test/runtests.jl
#
# The conformance cases are the same file the Python and C++ runners use, so a
# wrapper that disagrees with the engine — or with another wrapper — fails here.

using Test
using Printf

const HERE = @__DIR__
push!(LOAD_PATH, normpath(joinpath(HERE, "..", "src")))

using RSSN

const CASES = normpath(joinpath(HERE, "..", "..", "tests", "cases.tsv"))

struct Case
    name::String
    kind::String
    payload::String
    target::String
    tolerance::String
    bindings::String
    expect_text::String
    expect_value::String
    expect_status::String
end

function read_cases(path)
    cases = Case[]
    header = nothing
    for raw in eachline(path)
        line = chomp(raw)
        (isempty(strip(line)) || startswith(lstrip(line), "#")) && continue
        fields = split(line, '\t')
        if header === nothing
            header = fields
            continue
        end
        while length(fields) < length(header)
            push!(fields, "")
        end
        cases = push!(cases, Case(fields...))
    end
    return cases
end

function config_for(case)
    config = nothing
    if case.target == "numeric" || !isempty(case.bindings)
        config = Config()
        if case.target == "numeric"
            numeric!(config, isempty(case.tolerance) ? 1e-12 : parse(Float64, case.tolerance))
        end
        for pair in split(case.bindings, ';')
            isempty(strip(pair)) && continue
            parts = split(pair, '='; limit = 2)
            length(parts) == 2 || continue
            bind!(config, strip(parts[1]), parse(Float64, strip(parts[2])))
        end
    end
    return config
end

status_name(error::RssnError) = RSSN._STATUS_NAMES[error.status + 1]
"""A numeric literal from the table, keeping integers integral: an operator
with an `Int` parameter refuses a `Float` argument."""
function numeric_literal(text)
    literal = strip(text)
    (occursin('.', literal) || occursin('e', literal) || occursin('E', literal)) &&
        return parse(Float64, literal)
    return parse(Int, literal)
end

function run_case(case)
    kind = case.kind
    if kind == "manifest"
        digest = manifest_hash()
        @test digest != 0
        return (text = "", value = string(digest), status = "ok")
    elseif kind == "units"
        return (text = "", value = string(length(UNITS)), status = "ok")
    elseif kind == "constants"
        return (text = "", value = string(length(CONSTANTS)), status = "ok")
    elseif kind == "value"
        # The second ABI: an operator called with native objects. The display
        # form is the engine's own, so it is stable across languages.
        parts = split(case.payload, ' '; limit = 2)
        arguments = length(parts) == 2 ?
                    [numeric_literal(entry) for entry in split(parts[2], ',')
                     if !isempty(strip(entry))] : Any[]
        answer = call_value(default_session(), parts[1], arguments...)
        return (text = text(answer), value = "", status = "ok")
    elseif kind == "sim"
        parts = split(case.payload, ' ', limit = 2)
        parameters = length(parts) == 2 ? RSSN.json_decode(parts[2]) : Dict{String, Any}()
        output = run_scenario(parts[1], parameters)
        @test output isa AbstractDict
        names = sort!(collect(String.(keys(output))))
        return (text = join(names, ","), value = "", status = "ok")
    end

    session = Session()
    try
        if kind == "apply"
            arguments = [parse(Float64, strip(entry)) for entry in split(case.bindings, ',')
                         if !isempty(strip(entry))]
            built = apply(session, case.payload, arguments...)
            return (text = text(built), value = "", status = "ok")
        end
        if kind == "callback2"
            # Two callbacks in one call: the objective and its gradient that
            # `num_bfgs` declares, from the same start point every runner uses.
            parts = split(case.payload, ' '; limit = 2)
            arguments = isempty(strip(parts[2])) ? Any[] :
                        [numeric_literal(entry) for entry in split(parts[2], ',')
                         if !isempty(strip(entry))]
            objective(point) = sum((x - 1.0)^2 for x in point)
            gradient(point) = [2.0 * (x - 1.0) for x in point]
            answer = call_value(session, parts[1], [3.0, -2.0], arguments...,
                                objective, gradient)
            coordinate = RSSN.value(answer[1][1])
            return (text = @sprintf("%.10g", coordinate), value = "", status = "ok")
        end
        if kind == "callback"
            # A higher-order operator: `<op> <args>`, with the same residual
            # `f(x) = x*x - 2` every runner uses. The first field of the answer
            # record is the root.
            parts = split(case.payload, ' '; limit = 2)
            arguments = isempty(strip(parts[2])) ? Any[] :
                        [numeric_literal(entry) for entry in split(parts[2], ',')
                         if !isempty(strip(entry))]
            answer = call_value(session, parts[1], arguments..., x -> x * x - 2)
            root = RSSN.value(answer[1])
            return (text = @sprintf("%.10g", root), value = "", status = "ok")
        end
        answer = kind == "simplify" ? simplify(session, case.payload) :
                 compute(term(session, case.payload), config_for(case))
        value = ""
        if answer.value !== nothing && case.target == "numeric"
            value = repr(answer.value)
        end
        return (text = answer.text, value = value, status = "ok")
    finally
        close(session)
    end
end

@testset "conformance" begin
    for case in read_cases(CASES)
        wanted = isempty(case.expect_status) ? "ok" : case.expect_status
        observed = try
            run_case(case)
        catch error
            error isa RssnError ? (text = "", value = "", status = status_name(error)) : rethrow()
        end
        @test observed.status == wanted
        if wanted == "ok"
            @test observed.text == case.expect_text
            if !isempty(case.expect_value) && !isempty(observed.value)
                @test isapprox(parse(Float64, observed.value), parse(Float64, case.expect_value);
                               rtol = 1e-9, atol = 1e-12)
            end
        end
    end
end

@testset "terms and sessions" begin
    session = Session()
    x = sym(session, "x")

    @test text(x) == "x"
    @test text(x * x) == "x*x"   # a term is a request
    @test compute(x * x).text == "x^2"
    @test text(sym(session, "y") + 1) == "y + 1"
    @test text(2 * x) == "2*x"
    @test text(x - 1) == "x - 1"
    @test text(x^2) == "x^2"
    @test text(x / 2) == "x/2"
    @test text(2 / x) == "2/x"
    @test text(-x) == "-x"
    @test compute((x^2 - 1) / (x - 1)).text == "x + 1"
    @test handle(x) isa UInt32

    answer = compute(diff(sin(x) * cos(x), x))
    @test answer.text == "cos(2*x)"
    @test occursin("cos", answer.latex)   # property, not a field
    @test answer.reduced

    # Base extensions: sin, diff and friends accept a Term.
    @test text(sin(x)) == "sin(x)"
    @test compute(diff(sin(x), x)).text == "cos(x)"
    @test text(sqrt(x)) == "sqrt(x)"

    # Implicit coercion of numbers, rationals and text.
    @test text(term(session, 2)) == "2"
    @test text(term(session, 1 // 3)) == "1/3"
    @test text(term(session, "1/4")) == "1/4"

    # Tensor extraction.
    matrix = parse_term(session, "matmul(list(list(1,2),list(3,4)), list(list(0,1),list(1,0)))")
    entries = tensor(compute(matrix).term)
    @test size(entries) == (2, 2)
    @test entries[1, 2] == 1.0
    @test entries[2, 1] == 4.0

    # Batch computation.
    answers = compute_all(session, ["diff(x^3, x)", "integral(x^2, x)", "solve(x^2 - 4, x)"])
    @test length(answers) == 3
    @test all(answer -> answer isa Answer, answers)

    close(session)
end

@testset "errors" begin
    session = Session()
    @test_throws RssnError parse_term(session, "sin(")
    error = try
        parse_term(session, "sin(")
    catch caught
        caught
    end
    @test status_name(error) == "parse"
    close(session)
end

@testset "values" begin
    session = Session()
    backend = call_value(session, "num_backend_current")
    @test kind(backend) == "opaque"
    @test text(call_value(session, "num_backend_name", backend)) == "interpreter"

    tape = call_value(session, "num_ad_tape_new", 16)
    x = call_value(session, "num_ad_tape_var", tape, 1.5)
    z = call_value(session, "num_ad_var_mul", call_value(session, "num_ad_var_sin", x), x)
    @test isapprox(value(call_value(session, "num_ad_var_value", z)), sin(1.5) * 1.5; rtol = 1e-12)
    gradient = call_value(session, "num_ad_var_backward", z)
    @test isapprox(value(call_value(session, "num_ad_var_grad_in", x, gradient)),
                   cos(1.5) * 1.5 + sin(1.5); rtol = 1e-12)
    @test json(call_value(session, "num_ad_tape_len", tape)) isa AbstractDict

    # A list argument keeps its element kinds.
    list = value(session, Any[1, 2.5, 3])
    @test kind(list) == "tensor"
    @test length(list) == 3
    @test isapprox(value(list[2]), 2.5)

    # A record value: the ABI's `record` kind, which a dict builds.
    record = value(session, Dict("a" => 1.0, "b" => "two"))
    @test kind(record) == "record"
    @test length(record) == 2
    @test text(record[1]) == "1"
    @test json(record)["record"]["fields"][2][1] == "b"

    # Evaluation without any rule, in both spellings.
    @test eval_term(parse_term(session, "sin(x) + 1"); x = 0.0) == 1.0
    @test eval_term(parse_term(session, "sin(x) + 1"), Dict("x" => 0.0)) == 1.0
    @test eval_term(parse_term(session, "2 + 3")) == 5.0

    # The handles are typed: a variable is not a gradient.
    @test_throws RssnError call_value(session, "num_ad_var_grad_in", x, x)
    close(session)
end

@testset "json codec" begin
    @test RSSN.json_decode("{}") == Dict{String, Any}()
    @test RSSN.json_decode("[]") == Any[]
    @test RSSN.json_decode("1") == 1
    @test RSSN.json_decode("1.5") == 1.5
    @test RSSN.json_decode("true") === true
    @test RSSN.json_decode("null") === nothing
    @test RSSN.json_decode("\"a\\nb\"") == "a\nb"
    @test RSSN.json_decode("{\"a\": [1, 2.5, \"x\"]}")["a"][3] == "x"
    @test RSSN.json_encode(Dict("a" => [1, 2])) == "{\"a\":[1,2]}"
end

@testset "units and constants" begin
    @test length(UNITS) == 82
    @test unit("m").dimension == "Length"
    @test unit("metre").name == "m"
    @test unit_names()[1] isa String
    @test value(constant("c")) > 2.9e8
    @test constant("G_0").dimension == "Length^-2*Mass^-1*Time^3*Current^2"
end

@testset "the operator reference" begin
    @test length(op_names()) == 5411
    @test "diff" in op_names()
    @test length(op_categories()) == 40
    @test "calculus" in op_categories()
    @test "numerical_matrix" in op_rule_sets()

    # One operator, with its arity, category and the engine's own documentation.
    text = op_help("diff")
    @test occursin("diff", text)
    @test occursin("calculus", text)
    @test occursin("derivative", text)
    @test_throws ArgumentError op_help("definitely_not_an_operator")

    # A search over names and documentation.
    matches = op_search("brent", limit = 3)
    @test !isempty(matches)
    @test all(pair -> pair isa Tuple{String, String}, matches)
    @test any(pair -> pair[1] == "num_brent", op_search("brent", limit = 50))
    @test all(pair -> NAMES[pair[1]] in BY_CATEGORY["calculus"],
              op_search(category = "calculus", limit = 5))
    @test isempty(op_search(rule_set = "no_such_rule_set"))

    # The documentation table is the one the docstrings interpolate.
    @test occursin("derivative", OPERATOR_DOCS["diff"])
    @test occursin("host callbacks", OPERATOR_DOCS["num_brent"])
end

@testset "scenarios" begin
    @test "ising" in SCENARIO_NAMES
    @test haskey(SCENARIO_PARAMS, "poisson_multigrid")
    @test length(SCENARIO_NAMES) == 29
end

@testset "host callbacks" begin
    session = Session()

    # A `Function` argument becomes a host callback, with the signature the
    # operator declares for that position.
    answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64, x -> x * x - 2)
    @test isapprox(value(answer[1]), sqrt(2); rtol = 1e-10)

    # The `do` form puts the function last.
    answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x
        x * x - 2
    end
    @test isapprox(value(answer[1]), sqrt(2); rtol = 1e-10)

    # An explicit callback, with its signature spelled out.
    callback = Callback(session, x -> x * x - 2; signature = "(float) -> float")
    answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64, callback)
    @test isapprox(value(answer[1]), sqrt(2); rtol = 1e-10)
    close(callback)

    # The declared signatures are listed, and the parser reads the compact form.
    listed = callback_signatures(session, "num_brent", 5)
    @test length(listed) == 1
    @test listed[1]["params"] == ["float"]
    @test RSSN.parse_signature("(float) -> float") ==
          "{\"params\":[\"float\"],\"ret\":\"float\"}"
    @test isempty(callback_signatures(session, "num_brent", 1))

    # An objective and its gradient: two callbacks in one call.
    objective(xs) = sum((x - 1.0)^2 for x in xs)
    gradient(xs) = [2.0 * (x - 1.0) for x in xs]
    answer = call_value(session, "num_bfgs", [3.0, -2.0], 1e-10, 100, 1000,
                        objective, gradient)
    point = [value(x) for x in answer[1]]
    @test isapprox(point[1], 1.0; atol = 1e-6)
    @test isapprox(point[2], 1.0; atol = 1e-6)

    # A callback that throws aborts the operator with its message.
    failure = try
        call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do _
            Base.error("no thanks")
        end
        nothing
    catch caught
        caught
    end
    @test failure isa RssnError
    @test occursin("no thanks", failure.message)

    close(session)
end
