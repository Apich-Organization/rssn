# Package state and the library search path.
#
# The shared library is loaded at run time, never linked: any build of
# `rssn.dll` / `librssn.so` / `librssn.dylib` works with this package. The
# search order is
#
#   1. `rssn_load(path)`, called explicitly by the user;
#   2. `options(rssn.library = ...)`;
#   3. the `RSSN_LIB` environment variable (a file or a directory);
#   4. the package's own `libs/` directory (an installed copy shipped with it);
#   5. the build directories of a source checkout, relative to the working
#      directory (`target/debug`, `target-bindings/debug`, ...);
#   6. the bare library name, which lets the platform loader search its own
#      path.
#
# `library(rssn)` does not fail when the library is missing: the error is
# deferred to the first call that needs it, so the help pages stay readable.

.rssn_state <- new.env(parent = emptyenv())
.rssn_state$loaded <- FALSE
.rssn_state$path <- ""
.rssn_state$error <- NULL
.rssn_state$default <- NULL
.rssn_state$manifest <- NULL

#: The file name of the shared library on this platform.
.rssn_library_names <- function() {
    if (.Platform$OS.type == "windows") {
        c("rssn.dll", "librssn.dll")
    } else if (identical(Sys.info()[["sysname"]], "Darwin")) {
        c("librssn.dylib", "rssn.dylib", "librssn.so")
    } else {
        c("librssn.so", "librssn.so.0", "rssn.so")
    }
}

#' Every path `rssn` looks for the shared library in, in order
#'
#' @return A character vector of candidate files and directories.
#' @export
rssn_library_candidates <- function() {
    candidates <- character()
    option <- getOption("rssn.library")
    if (is.character(option) && length(option) > 0L) {
        candidates <- c(candidates, option)
    }
    environment <- Sys.getenv("RSSN_LIB", "")
    if (nzchar(environment)) {
        candidates <- c(candidates, environment)
    }
    installed <- system.file("libs", package = "rssn")
    if (nzchar(installed)) {
        candidates <- c(candidates, installed, system.file("lib", package = "rssn"))
    }
    # A source checkout: the working directory, then the parents of the package
    # when it is being developed in place.
    roots <- c(".", "..", "../..", "../../..")
    profiles <- c("debug", "ffi", "release")
    for (root in roots) {
        for (directory in c("target", "target-bindings", "target-ffi",
                            file.path("dist", "rssn"))) {
            for (profile in profiles) {
                candidates <- c(candidates, file.path(root, directory, profile))
            }
        }
    }
    c(candidates, .rssn_library_names()[1])
}

#' Point the package at a specific rssn shared library
#'
#' @param path A file, a directory holding the library, or `NULL` to restore the
#'   automatic search.
#' @return The path of the library that was loaded, invisibly.
#' @export
rssn_load <- function(path = NULL) {
    candidates <- if (is.null(path)) rssn_library_candidates() else as.character(path)
    loaded <- .Call(C_rssn_load, candidates)
    .rssn_state$loaded <- TRUE
    .rssn_state$path <- loaded
    .rssn_state$error <- NULL
    invisible(loaded)
}

#' The path of the rssn shared library in use
#'
#' @return A string; `""` when nothing has been loaded yet.
#' @export
rssn_library <- function() {
    if (!.rssn_state$loaded) {
        .rssn_require()
    }
    .Call(C_rssn_libpath)
}

# Called by every entry point. The cost is one logical test in the common case.
.rssn_require <- function() {
    if (.rssn_state$loaded) {
        return(invisible(NULL))
    }
    error <- tryCatch({
        rssn_load()
        NULL
    }, error = function(condition) condition)
    if (!is.null(error)) {
        .rssn_state$error <- error
        stop(error)
    }
    invisible(NULL)
}

#' The rssn version string
#'
#' @return A string such as `"3.0.0"`.
#' @export
rssn_version <- function() {
    .rssn_require()
    .Call(C_rssn_version)
}

.onLoad <- function(libname, pkgname) {
    tryCatch(rssn_load(), error = function(condition) {
        .rssn_state$error <- condition
        .rssn_state$loaded <- FALSE
    })
    invisible(NULL)
}

.onUnload <- function(libpath) {
    .rssn_state$loaded <- FALSE
    invisible(NULL)
}
