"""The backend selectors the loaded RSSN library can use on this host."""
function available_backends(session::Session=default_session())
    choices = call_value(session, "num_backend_available")
    try
        names = String[]
        for index in 1:length(choices)
            choice = choices[index]
            try
                push!(names, text(choice))
            finally
                close(choice)
            end
        end
        return names
    finally
        close(choices)
    end
end

"""Return the active backend's name without changing the selection."""
function current_backend(session::Session=default_session())
    name = try
        call_value(session, "num_backend_current_name")
    catch error
        # C ABI status 4 means that an older library lacks this operator.
        error isa RssnError && error.status == 4 || rethrow()
        nothing
    end
    if name !== nothing
        try
            return text(name)
        finally
            close(name)
        end
    end

    backend = call_value(session, "num_backend_current")
    try
        name = call_value(session, "num_backend_name", backend)
        try
            return text(name)
        finally
            close(name)
        end
    finally
        try
            release_backend!(backend)
        finally
            close(backend)
        end
    end
end

"""
    open_backend(session, name)
    open_backend(name)

Open a named backend for explicit compilation in this session. Release its
native handle with [`release_backend!`](@ref) when finished.
"""
function open_backend(session::Session, name::AbstractString)
    occursin('\0', name) && throw(ArgumentError("the backend selector contains NUL"))
    return call_value(session, "num_backend_named", name)
end

open_backend(name::AbstractString) = open_backend(default_session(), name)

"""Release a native backend handle and invalidate its copies."""
function release_backend!(backend::Value)
    close(call_value(backend.session, "num_backend_release", backend))
    return nothing
end

"""
    set_default_backend!(session, name)
    set_default_backend!(name)

Select the process-wide backend and return its effective name. Use `"auto"` to
restore the library's automatic selection policy.
"""
function set_default_backend!(session::Session, name::AbstractString)
    occursin('\0', name) && throw(ArgumentError("the backend selector contains NUL"))
    selected = call_value(session, "num_backend_select_default", name)
    try
        return text(selected)
    finally
        close(selected)
    end
end

set_default_backend!(name::AbstractString) = set_default_backend!(default_session(), name)
