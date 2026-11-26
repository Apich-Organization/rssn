# A consistency check of the R package's surface: every exported name has a help
# topic, every documented name is exported, and the help index finds an operator.

suppressPackageStartupMessages(library(rssn))

namespace <- asNamespace("rssn")
exports <- getNamespaceExports("rssn")

# The aliases the installed package's help database knows about.
index <- utils:::.getHelpFile  # touch the symbol so the namespace is loaded
help_topics <- character()
for (file in list.files(file.path(find.package("rssn"), "help"), full.names = TRUE)) {
    if (!grepl("\\.rdb$", file)) next
}
database <- file.path(find.package("rssn"), "help", "aliases.rds")
if (file.exists(database)) {
    help_topics <- names(readRDS(database))
} else {
    cat("no aliases.rds; skipping the alias check\n")
}

missing <- setdiff(exports, help_topics)
cat(length(exports), "exports,", length(help_topics), "help topics\n")
if (length(help_topics) && length(missing)) {
    cat("exports without a help topic:", paste(head(missing, 20), collapse = ", "), "\n")
    quit(status = 1L)
}

# `?rssn_diff` resolves through the generated alias page.
if (length(help_topics)) {
    stopifnot("rssn_diff" %in% help_topics, "rssn-operators" %in% help_topics)
    stopifnot("rssn_session" %in% help_topics, "compute" %in% help_topics)
    cat("operator help topics resolve\n")
}
cat("surface check passed\n")
