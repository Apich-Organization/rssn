"""
    encode_npy(session, samples) -> Vector{UInt8}

Encode a numeric array as little-endian Float64 NumPy data in memory. The
array's original rank and extents are preserved.
"""
function encode_npy(session::Session, samples::AbstractArray{<:Real})
    dimensions = collect(size(samples))
    dimension_cell = tensor_cells(session, dimensions; exemplar=isempty(dimensions) ? 0 : nothing)
    try
        rank = ndims(samples)
        row_major = rank < 2 ? samples : permutedims(samples, rank:-1:1)
        data = value(session, Float64.(vec(row_major)))
        try
            payload = call_value(session, "num_npy_encode", dimension_cell, data)
            try
                return hex2bytes(text(payload))
            finally
                close(payload)
            end
        finally
            close(data)
        end
    finally
        close(dimension_cell)
    end
end

"""
    decode_npy(session, payload) -> Value

Decode NumPy Float32 or Float64 bytes as an owned binary64 tensor. Use
[`tensor`](@ref) to read its samples as a Julia array.
"""
function decode_npy(session::Session, payload::AbstractVector{UInt8})
    record = call_value(session, "num_npy_decode", bytes2hex(payload))
    try
        dimensions = record["shape"]
        shape = Csize_t[]
        try
            for index in 1:length(dimensions)
                extent = dimensions[index]
                try
                    size = parse(Csize_t, text(extent))
                    size <= typemax(Int) || throw(ArgumentError("NumPy extent exceeds Julia's array index range"))
                    push!(shape, size)
                finally
                    close(extent)
                end
            end
        finally
            close(dimensions)
        end
        data = record["data"]
        samples = try
            convert(Vector{Cdouble}, vec(tensor(data)))
        finally
            close(data)
        end
        out = Ref{UInt32}(0)
        _check(value_new_tensor(session.ptr, samples, shape, out))
        return Value(session, out[])
    finally
        close(record)
    end
end
