# Build an exact rational matrix while preserving Julia array dimensions.
using RSSN

session = Session()
try
    wide = (big(1) << 64) + 1
    matrix = tensor_cells(session, Any[wide 2//3; 3 4])
    normalized = call_value(session, "num_native_qmatrix_new", matrix, 1)
    determinant = call_value(session, "num_native_qlinalg_det", normalized)
    @assert text(determinant) == string(4 * wide - 2)
    println("Exact determinant: ", text(determinant))

    scalar = tensor_cells(session, fill(1//3))
    @assert isempty(json(scalar)["tensor"]["shape"])
    @assert text(scalar[1]) == "1/3"
    empty = tensor_cells(session, Array{Int}(undef, 2, 0, 3); exemplar=0)
    @assert json(empty)["tensor"]["shape"] == [2, 0, 3]
    @assert json(empty)["tensor"]["elem"] == "int"

    try
        tensor_cells(session, [[1], [2, 3]])
        error("ragged tensor input was accepted")
    catch failure
        @assert failure isa ArgumentError
    end
    foreign = Session()
    try
        cell = value(foreign, wide)
        try
            tensor_cells(session, [cell])
            error("a foreign tensor cell was accepted")
        catch failure
            @assert failure isa ArgumentError
        end
    finally
        close(foreign)
    end
finally
    close(session)
end
