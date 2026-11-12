## Stage 09: the pre-registered test of PLWH against HC, its sensitivity
## analyses, the exploratory INR against IR contrast, and power
## (preregistration.md, sections 10, 11 and 13).
##
## Reads the log2 matrices of the representative features written by
## stage 08 (results/<mode>_matrix_<variant>.tsv.gz) and the sample
## sheets (results/<mode>_samples.tsv).
##
## Usage: scripts/lib/in_env.sh Rscript scripts/09_statistics.R [--force]
##
## Outputs:
##   results/<mode>_limma.tsv.gz           primary result per feature, with annotation
##   results/<mode>_sensitivity.tsv.gz     every sensitivity analysis, per feature
##   results/sensitivity_summary.tsv       primary hits kept and sign agreement, per analysis
##   results/<mode>_limma_inr_ir.tsv.gz    exploratory INR against IR
##   results/test_counts.tsv               features tested and called, per mode and analysis; n_reference
##                                         is HC (IR for the exploratory contrast), n_comparison PLWH (INR)
##   results/power.tsv                     minimum detectable standardised difference

suppressPackageStartupMessages({
    library("data.table")
    library("limma")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
res_dir <- file.path(root, "results")
if (length(args) && args[1] == "--root") res_dir <- file.path(args[2], "results")
set.seed(20261109)

modes <- c(pos = "positive", neg = "negative")
min_per_group <- 3L

## ---- Helpers -------------------------------------------------------------

read_matrix <- function(tag, variant) {
    f <- file.path(res_dir, sprintf("%s_matrix_%s.tsv.gz", tag, variant))
    if (!file.exists(f)) return(NULL)
    d <- fread(f)
    m <- as.matrix(d[, -1])
    rownames(m) <- d$feature_ID
    m
}

## Features with at least 3 non-missing values in each group are tested
## (section 10); the same rule for every analysis, so they stay comparable.
testable <- function(m, g) {
    ok <- rowSums(!is.na(m[, g == "HC", drop = FALSE])) >= min_per_group &
        rowSums(!is.na(m[, g == "PLWH", drop = FALSE])) >= min_per_group
    m[ok, , drop = FALSE]
}

limma_fit <- function(m, g, extra = NULL, ref = "HC") {
    group <- relevel(factor(g), ref = ref)
    design <- if (is.null(extra)) model.matrix(~ group) else model.matrix(~ group + extra)
    fit <- eBayes(lmFit(m, design))
    coef <- grep("^group", colnames(design), value = TRUE)[1]
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none", confint = 0.95)
    data.table(feature_ID = as.integer(rownames(tt)), logFC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R,
               AveExpr = tt$AveExpr, t = tt$t, p = tt$P.Value, q = p.adjust(tt$P.Value, "BH"))
}

per_feature_test <- function(m, g, fun) {
    p <- apply(m, 1, function(x) {
        a <- x[g == "PLWH"]; b <- x[g == "HC"]
        fun(a[!is.na(a)], b[!is.na(b)])$p.value
    })
    lfc <- apply(m, 1, function(x) mean(x[g == "PLWH"], na.rm = TRUE) - mean(x[g == "HC"], na.rm = TRUE))
    data.table(feature_ID = as.integer(rownames(m)), logFC = lfc, p = p, q = p.adjust(p, "BH"))
}

## Half of the feature's smallest observed value; on the log2 scale that is
## the minimum minus one.
half_min <- function(m) {
    t(apply(m, 1, function(x) { x[is.na(x)] <- min(x, na.rm = TRUE) - 1; x }))
}

## Exact two-sample t power by integrating the normal tails over the
## chi-square distribution of the variance (section 11).
power_t <- function(d, n1, n2, alpha) {
    df <- n1 + n2 - 2; k <- sqrt(n1 * n2 / (n1 + n2)); tc <- qt(1 - alpha / 2, df)
    integrate(function(v) (pnorm(tc * sqrt(v / df) - d * k, lower.tail = FALSE) + pnorm(-tc * sqrt(v / df) - d * k)) * dchisq(v, df),
              0, Inf, rel.tol = 1e-10)$value
}
min_d <- function(n1, n2, alpha) uniroot(function(d) power_t(d, n1, n2, alpha) - 0.8, c(0.05, 5), tol = 1e-8)$root

## ---- Analyses ----------------------------------------------------------------

summary_rows <- list()
counts <- list()
power_rows <- list()

for (tag in names(modes)) {
    samples <- fread(file.path(res_dir, sprintf("%s_samples.tsv", tag)))
    ann <- fread(file.path(res_dir, sprintf("%s_representatives.tsv", tag)))
    primary_m <- read_matrix(tag, "primary")
    incl <- samples[excluded == FALSE]
    g <- incl$group[match(colnames(primary_m), incl$sample_name)]
    m <- testable(primary_m, g)

    prim <- limma_fit(m, g)
    prim[, `:=`(n_hc = rowSums(!is.na(m[, g == "HC", drop = FALSE])), n_plwh = rowSums(!is.na(m[, g == "PLWH", drop = FALSE])))]
    out <- merge(prim, ann[, .(feature_ID, group_ID, `m/z`, RT, adduct, level, annotation, InChIKey, formula,
                               similarity, matched_peak_number, precursor_ppm)], by = "feature_ID", all.x = TRUE)
    fwrite(out[order(p)], file.path(res_dir, sprintf("%s_limma.tsv.gz", tag)), sep = "\t")
    hits <- prim[q < 0.05, feature_ID]
    counts[[length(counts) + 1]] <- data.table(mode = modes[[tag]], analysis = "primary limma", tested = nrow(prim),
                                               called = length(hits), higher_in_plwh = prim[q < 0.05 & logFC > 0, .N],
                                               lower_in_plwh = prim[q < 0.05 & logFC < 0, .N],
                                               n_reference = sum(g == "HC"), n_comparison = sum(g == "PLWH"))

    ## Sensitivity analyses: each returns feature-level results comparable with the primary.
    sens <- list()
    sens[["welch"]] <- per_feature_test(m, g, function(a, b) t.test(a, b, var.equal = FALSE))
    sens[["wilcoxon"]] <- per_feature_test(m, g, function(a, b) wilcox.test(a, b, exact = FALSE))
    sens[["half_minimum_imputation"]] <- limma_fit(half_min(m), g)
    ord <- incl$injection_order[match(colnames(m), incl$sample_name)]
    ord_s <- as.numeric(scale(ord))
    sens[["injection_order_covariate"]] <- limma_fit(m, g, extra = cbind(ord_s, ord_s^2))
    for (v in c("all_injections", "area", "rsd20", "with_lowess", "without_lowess", "without_qc11")) {
        mv <- read_matrix(tag, v)
        if (is.null(mv)) next
        sv <- if (v == "all_injections") samples else incl
        gv <- sv$group[match(colnames(mv), sv$sample_name)]
        sens[[v]] <- limma_fit(testable(mv, gv), gv)
    }
    long <- rbindlist(lapply(names(sens), function(n) sens[[n]][, .(analysis = n, feature_ID, logFC, p, q)]))
    fwrite(long, file.path(res_dir, sprintf("%s_sensitivity.tsv.gz", tag)), sep = "\t")
    for (n in names(sens)) {
        s <- sens[[n]]
        kept <- s[feature_ID %in% hits & q < 0.05, feature_ID]
        agree <- merge(prim[feature_ID %in% hits, .(feature_ID, l0 = logFC)], s[, .(feature_ID, l1 = logFC)], by = "feature_ID")
        summary_rows[[length(summary_rows) + 1]] <- data.table(
            mode = modes[[tag]], analysis = n, tested = nrow(s), called = s[q < 0.05, .N],
            primary_hits = length(hits), primary_hits_kept = length(kept),
            primary_hits_same_sign = agree[sign(l0) == sign(l1), .N], primary_hits_compared = nrow(agree))
        counts[[length(counts) + 1]] <- data.table(mode = modes[[tag]], analysis = n, tested = nrow(s), called = s[q < 0.05, .N],
                                                   higher_in_plwh = s[q < 0.05 & logFC > 0, .N], lower_in_plwh = s[q < 0.05 & logFC < 0, .N],
                                                   n_reference = NA_integer_, n_comparison = NA_integer_)
    }

    ## Exploratory: INR against IR among PLWH, same model and rules.
    sub <- incl[group == "PLWH"]
    ms <- primary_m[, sub$sample_name, drop = FALSE]
    gs <- ifelse(sub$subgroup == "INR", "PLWH", "HC")  # testable() reads the two labels; INR plays PLWH, IR plays HC
    ms <- testable(ms, gs)
    inr <- limma_fit(ms, ifelse(gs == "PLWH", "INR", "IR"), ref = "IR")
    fwrite(inr[order(p)], file.path(res_dir, sprintf("%s_limma_inr_ir.tsv.gz", tag)), sep = "\t")
    counts[[length(counts) + 1]] <- data.table(mode = modes[[tag]], analysis = "exploratory INR against IR", tested = nrow(inr),
                                               called = inr[q < 0.05, .N], higher_in_plwh = NA_integer_, lower_in_plwh = NA_integer_,
                                               n_reference = sum(sub$subgroup == "IR"), n_comparison = sum(sub$subgroup == "INR"))

    n1 <- sum(g == "PLWH"); n2 <- sum(g == "HC")
    power_rows[[length(power_rows) + 1]] <- data.table(
        mode = modes[[tag]], n_plwh = n1, n_hc = n2, features_tested = nrow(prim),
        d_alpha_0.05 = min_d(n1, n2, 0.05), alpha_bonferroni = 0.05 / nrow(prim),
        d_alpha_bonferroni = min_d(n1, n2, 0.05 / nrow(prim)))
    message(sprintf("%s: %d features tested, %d called at q < 0.05 (%d higher in PLWH, %d lower)", modes[[tag]], nrow(prim),
                    length(hits), prim[q < 0.05 & logFC > 0, .N], prim[q < 0.05 & logFC < 0, .N]))
}

fwrite(rbindlist(summary_rows), file.path(res_dir, "sensitivity_summary.tsv"), sep = "\t")
fwrite(rbindlist(counts), file.path(res_dir, "test_counts.tsv"), sep = "\t")
fwrite(rbindlist(power_rows), file.path(res_dir, "power.tsv"), sep = "\t")
