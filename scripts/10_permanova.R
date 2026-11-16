#!/usr/bin/env Rscript
## PERMANOVA PLWH vs HC per mode (exploratory, prereg section 12)

suppressPackageStartupMessages({
    library("data.table")
    library("vegan")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
res_dir <- file.path(root, "results")
if (length(args) && args[1] == "--root") res_dir <- file.path(args[2], "results")

modes <- c(pos = "positive", neg = "negative")
rows <- list()
for (tg in names(modes)) {
    tag <- modes[[tg]]
    s <- fread(file.path(res_dir, sprintf("%s_samples.tsv", tg)))[excluded == FALSE]
    d <- fread(file.path(res_dir, sprintf("%s_matrix_primary.tsv.gz", tg)))
    m <- t(as.matrix(d[, s$sample_name, with = FALSE]))
    ## dist() needs complete data: half minimum here, and only here
    m <- apply(m, 2, function(x) { x[is.na(x)] <- min(x, na.rm = TRUE) - 1; x })
    m <- m[, apply(m, 2, sd) > 0, drop = FALSE]
    set.seed(20261109)
    fit <- adonis2(dist(scale(m)) ~ group, data = data.frame(group = s$group), permutations = 9999)
    rows[[tag]] <- data.table(mode = tag, n_plwh = sum(s$group == "PLWH"), n_hc = sum(s$group == "HC"), features = ncol(m),
                              R2 = fit$R2[1], F = fit$F[1], p = fit$`Pr(>F)`[1], permutations = 9999,
                              missing_values = "half minimum, this analysis only")
    message(sprintf("%s: R2 %.3f, F %.2f, p %.4f", tag, fit$R2[1], fit$F[1], fit$`Pr(>F)`[1]))
}
fwrite(rbindlist(rows), file.path(res_dir, "permanova.tsv"), sep = "\t")
