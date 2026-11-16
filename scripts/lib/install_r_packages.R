#!/usr/bin/env Rscript
## called by 01_install.sh; records what got installed rather than pinning, so a rebuild reports its own versions

args <- commandArgs(trailingOnly = TRUE)
out_tsv <- args[1]
out_session <- args[2]

## R.utils because fread won't open .gz without it
pkgs <- c("limma", "ropls", "vegan", "data.table", "ggplot2", "ggrepel",
          "patchwork", "jsonlite", "statmod", "R.utils")

repos <- "https://cloud.r-project.org"
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
bioc_version <- as.character(BiocManager::version())
cat("R:", R.version.string, "\nBioconductor:", bioc_version, "\n")

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) {
    ## binaries only: conda R on Windows has no compiler
    options(install.packages.compile.from.source = "never")
    BiocManager::install(missing, update = FALSE, ask = FALSE, type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))
}
still_missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(still_missing)) stop("packages failed to install: ", paste(still_missing, collapse = ", "))

ip <- installed.packages()[, c("Package", "Version")]
ip <- ip[ip[, "Package"] %in% c(pkgs, "BiocManager"), , drop = FALSE]
ip <- ip[order(ip[, "Package"]), , drop = FALSE]
write.table(data.frame(package = ip[, "Package"], version = ip[, "Version"],
                       r_version = R.version.string, bioconductor = bioc_version,
                       recorded = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
            out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

## a broken install should fail here, not halfway through stage 09
suppressPackageStartupMessages(for (p in pkgs) library(p, character.only = TRUE))
writeLines(capture.output(sessionInfo()), out_session)
cat("wrote", out_tsv, "and", out_session, "\n")
