# Higher-order operators: Julia functions the engine calls.
#
#     RSSN_LIB=/path/to/rssn.dll julia --project=bindings/julia bindings/examples/julia/callbacks.jl
#
# Some operators take a function and call it — root finders, quadrature,
# optimisers, bootstrap. In Julia any function works: the binding turns it into
# the C function pointer the engine calls and converts the arguments and the
# result for you.

using RSSN

function root_of_two()
    session = Session()
    answer = call_value(session, "num_brent",
                        1.0, 2.0,      # a, b: a bracketing interval
                        1e-12,         # tolerance
                        64,            # iteration budget
                        x -> x * x - 2)
    println("num_brent:        sqrt(2) = ", value(answer[1]))
    println("                  evaluations = ", value(answer[2]))
    close(session)
end

function root_of_two_do()
    session = Session()
    # `do` puts the function last, which is where the callback sits in most of
    # these operators.
    answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do x
        x * x - 2
    end
    println("do-block:         sqrt(2) = ", value(answer[1]))
    close(session)
end

function minimise()
    session = Session()
    objective(xs) = sum((x - 1.0)^2 for x in xs)
    gradient(xs) = [2.0 * (x - 1.0) for x in xs]
    answer = call_value(session, "num_bfgs", [3.0, -2.0], 1e-10, 100, 1000,
                        objective, gradient)
    println("num_bfgs:         minimum at ", [value(x) for x in answer[1]])
    println("                  value = ", value(answer[2]))
    close(session)
end

function explicit_callback()
    session = Session()
    println("declared at num_brent/5: ", callback_signatures(session, "num_brent", 5))
    callback = Callback(session, x -> x * x - 2; signature = "(float) -> float")
    answer = call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64, callback)
    println("explicit Callback: ", value(answer[1]))
    close(callback)
    close(session)
end

function failures_cross()
    session = Session()
    try
        call_value(session, "num_brent", 1.0, 2.0, 1e-12, 64) do _
            error("no thanks")
        end
    catch error
        error isa RssnError || rethrow()
        println("callback failure: ", error)
    finally
        close(session)
    end
end

root_of_two()
root_of_two_do()
minimise()
explicit_callback()
failures_cross()
