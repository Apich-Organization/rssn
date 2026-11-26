# Metadata: the kernel manifest, the operator reference, the unit and constant
# tables.
#
# The manifest is the engine's own description of itself: the operator table
# (name, arity, category, rule set, documentation), the units with their exact
# scales, the physical constants and the feature flags. The bindings are
# generated from it, so `OPS` and `OPERATORS` cannot drift from the engine
# without the generator's `--check` failing.
#
# `op_help("diff")` prints one operator's documentation; `op_table()` searches
# the whole table.

#' The kernel manifest
#'
#' The engine's own description of itself: the versions and feature flags, the
#' rule sets, the base dimensions, the value kinds, the operator table, the
#' units, the constants and the simulation catalog.
#'
#' `ops` is the reference data frame (`OPERATORS`) and `units`, `constants` and
#' `scenarios` are the generated tables, so this is assembled without parsing
#' the engine's JSON: it is the same data the package was built from. Use
#' [manifest_json()] for the raw text.
#'
#' @return A named list, assembled once per session.
#' @examples
#' manifest()$rssn_version
#' nrow(manifest()$ops)
#' @export
manifest <- function() {
    if (is.null(.rssn_state$manifest)) {
        .rssn_state$manifest <- c(
            MANIFEST_META,
            list(ops = OPERATORS, units = UNITS, constants = CONSTANTS,
                 scenarios = SCENARIOS)
        )
    }
    .rssn_state$manifest
}

#' The kernel manifest as the JSON text the engine emits
#'
#' The engine's own serialisation, about 18 MB.  [manifest()] is the fast path:
#' it assembles the same information from the tables the package was generated
#' from, and returns them as R data frames.
#'
#' @return A length-one character vector.
#' @export
manifest_json <- function() {
    .rssn_require()
    .Call(C_rssn_manifest_json)
}

#' The first 64 bits of the manifest hash
#'
#' Returned as a decimal string: the value is a `uint64_t` and R's double cannot
#' hold it exactly.
#'
#' @return A length-one character vector.
#' @export
manifest_hash <- function() {
    .rssn_require()
    .Call(C_rssn_manifest_hash)
}

# ---------------------------------------------------------------------------
# the operator reference
# ---------------------------------------------------------------------------

#' Look one operator up
#'
#' @param name The operator name, as the manifest spells it (`"diff"`).
#' @return A one-row data frame of class `rssn_operator`; an error when the name
#'   is not registered.
#' @examples
#' op_help("diff")
#' op_help("num_brent")
#' @export
op_help <- function(name) {
    name <- as.character(name)[[1L]]
    row <- OPERATORS[OPERATORS$name == name, , drop = FALSE]
    if (nrow(row) == 0L) {
        stop(sprintf("'%s' is not a registered operator", name))
    }
    structure(row, class = c("rssn_operator", "data.frame"))
}

#' @export
print.rssn_operator <- function(x, ...) {
    cat(x$name[[1L]], " (", x$signature[[1L]], ")\n", sep = "")
    cat("  binding   : rssn_", x$name[[1L]], "()\n", sep = "")
    cat("  rule set  : ", x$rule_set[[1L]], "\n", sep = "")
    cat("  category  : ", x$category[[1L]], "\n", sep = "")
    cat("  returns   : ", x$returns[[1L]], "\n", sep = "")
    if (nzchar(x$since[[1L]])) cat("  since     : ", x$since[[1L]], "\n", sep = "")
    if (nzchar(x$stability[[1L]])) cat("  stability : ", x$stability[[1L]], "\n", sep = "")
    cat("\n")
    cat(strwrap(x$doc[[1L]], width = 78L, prefix = "  "), sep = "\n")
    cat("\n")
    invisible(x)
}

#' Search the operator table
#'
#' @param pattern A regular expression matched against the operator name and its
#'   documentation; `NULL` keeps every operator.
#' @param category Keep one manifest category (`"calculus"`, `"numerical_matrix"`,
#'   ...); `NULL` keeps all of them.
#' @param rule_set Keep one rule set; `NULL` keeps all of them.
#' @param limit The maximum number of rows to return.
#' @return A data frame of the matching operators.
#' @examples
#' nrow(op_table(category = "calculus"))
#' head(op_table(pattern = "brent"))
#' @export
op_table <- function(pattern = NULL, category = NULL, rule_set = NULL, limit = 200L) {
    table <- OPERATORS
    if (!is.null(category)) {
        table <- table[table$category == category, , drop = FALSE]
    }
    if (!is.null(rule_set)) {
        table <- table[table$rule_set == rule_set, , drop = FALSE]
    }
    if (!is.null(pattern)) {
        matches <- grepl(pattern, table$name, ignore.case = TRUE) |
            grepl(pattern, table$doc, ignore.case = TRUE)
        table <- table[matches, , drop = FALSE]
    }
    if (!is.null(limit) && nrow(table) > limit) {
        table <- table[seq_len(limit), , drop = FALSE]
    }
    rownames(table) <- NULL
    table
}

#' Every operator name
#'
#' @return A sorted character vector of the 5411 registered operator names.
#' @export
op_names <- function() OPERATORS$name

# ---------------------------------------------------------------------------
# units and constants
# ---------------------------------------------------------------------------

#' Look a unit up by name or alias
#'
#' @param name A unit name or alias: `"m"`, `"metre"`, `"ft"`, ...
#' @return A unit record: a list with `name`, `display`, `dimension`, the exact
#'   `numer`/`denom`/`pi` scale, the affine offset and the aliases.
#' @examples
#' unit_get("ft")
#' @export
unit_get <- function(name) {
    unit <- UNIT_BY_ALIAS[[as.character(name)[[1L]]]]
    if (is.null(unit)) stop(sprintf("unknown unit '%s'", name))
    unit
}

#' Every canonical unit name
#' @return A sorted character vector.
#' @export
unit_names <- function() sort(names(UNITS))

#' The unit table as a data frame
#'
#' @return A data frame with one row per unit: name, display, dimension, the
#'   exact scale as text, the `pi` power, the affine offset, the aliases and the
#'   category.
#' @examples
#' head(units_table())
#' @export
units_table <- function() {
    rows <- lapply(UNITS, function(unit) {
        data.frame(
            name = unit$name, display = unit$display, dimension = unit$dimension,
            numer = unit$numer, denom = unit$denom, pi = unit$pi,
            affine = unit$affine, prefixable = unit$prefixable,
            category = unit$category, aliases = paste(unit$aliases, collapse = " "),
            stringsAsFactors = FALSE
        )
    })
    table <- do.call(rbind, rows)
    rownames(table) <- NULL
    table
}

#' The numeric scale of a unit relative to its base unit
#'
#' The scale is the exact rational `numer / denom * pi^pi`.
#'
#' @param unit A unit record from [unit_get()].
#' @return A number.
#' @examples
#' unit_scale(unit_get("ft"))    # 0.3048
#' @export
unit_scale <- function(unit) {
    as.numeric(unit$numer) / as.numeric(unit$denom) * pi^unit$pi
}

#' Look a physical constant up by display name or symbol
#'
#' @param name A display name (`"speed of light"`) or a symbol (`"c"`).
#' @return A constant record: a list with `name`, `symbol`, the exact
#'   `numer`/`denom`/`pi` value, the dimension, the unit and the uncertainty.
#' @examples
#' constant_get("c")
#' @export
constant_get <- function(name) {
    name <- as.character(name)[[1L]]
    constant <- CONSTANTS[[name]]
    if (is.null(constant)) constant <- CONSTANT_BY_SYMBOL[[name]]
    if (is.null(constant)) stop(sprintf("unknown constant '%s'", name))
    constant
}

#' The constant table as a data frame
#'
#' @return A data frame with one row per constant: name, symbol, the exact value
#'   as text, the `pi` power, the dimension, the unit, whether it is exact, and
#'   the standard uncertainty.
#' @export
constants_table <- function() {
    rows <- lapply(CONSTANTS, function(constant) {
        data.frame(
            name = constant$name, symbol = constant$symbol, numer = constant$numer,
            denom = constant$denom, pi = constant$pi, dimension = constant$dimension,
            unit = constant$unit, exact = constant$exact,
            uncertainty = constant$uncertainty %||% NA_character_,
            source = constant$source, stringsAsFactors = FALSE
        )
    })
    table <- do.call(rbind, rows)
    rownames(table) <- NULL
    table
}

#' The numeric value of a physical constant
#'
#' The value is `numer / denom * pi^pi`.
#'
#' @param constant A constant record from [constant_get()].
#' @return A number.
#' @examples
#' constant_value(constant_get("c"))    # 299792458
#' @export
constant_value <- function(constant) {
    as.numeric(constant$numer) / as.numeric(constant$denom) * pi^constant$pi
}
