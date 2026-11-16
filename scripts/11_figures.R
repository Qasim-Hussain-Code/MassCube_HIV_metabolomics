#!/usr/bin/env Rscript
## Figures from results/; only the QC RSDs in fig 2 are computed here. PNG at 300 dpi plus PDF.

suppressPackageStartupMessages({
    library("data.table")
    library("ggplot2")
    library("patchwork")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) && args[1] == "--root") args[2] else root
res <- function(f) file.path(base, "results", f)
fig_dir <- file.path(base, "figures")
dir.create(file.path(fig_dir, "supplementary"), recursive = TRUE, showWarnings = FALSE)
set.seed(20261109)
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank()))

## Okabe-Ito; a group keeps its colour in every figure
okabe <- c(HC = "#0072B2", PLWH = "#D55E00", QC = "#E69F00", INR = "#CC79A7", IR = "#009E73")
level_col <- c(`2` = "#009E73", `3` = "#CC79A7", `4` = "#BBBBBB")
modes <- c(pos = "positive", neg = "negative")

save_fig <- function(p, name, w, h) {
    for (ext in c("png", "pdf")) ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300, bg = "white")
    message("wrote figures/", name, ".png and .pdf")
}
n_label <- function(g) {
    n <- table(g)
    setNames(sprintf("%s (n = %d)", names(n), as.integer(n)), names(n))
}

## ---- Figure 1: QC overview ---------------------------------------------------

panels <- list()
for (tg in names(modes)) {
    inj <- fread(res(sprintf("%s_injection_metrics.tsv", tg)))
    ord <- fread(res(sprintf("%s_run_order.tsv", tg)))
    qcv <- fread(res("injection_qc.tsv"))[mode == modes[[tg]], .(sample_name, excluded)]
    d <- merge(merge(inj, ord[, .(sample_name, injection_order, group)], by = "sample_name"), qcv, by = "sample_name")
    lab <- n_label(d$group)
    panels[[paste0(tg, "a")]] <- ggplot(d, aes(injection_order, summed_height / 1e9, colour = group, shape = excluded)) +
        geom_point(size = 1.6) + scale_colour_manual(values = okabe, labels = lab) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c(`FALSE` = "kept", `TRUE` = "excluded")) +
        labs(x = "injection order", y = "summed feature height (x 10^9, counts)", colour = NULL, shape = NULL,
             title = sprintf("%s mode: total intensity by run order", modes[[tg]]))
    pc <- fread(res(sprintf("%s_pca.tsv", tg)))
    panels[[paste0(tg, "b")]] <- ggplot(pc, aes(PC1, PC2, colour = group)) + geom_point(size = 1.6, alpha = 0.85) +
        scale_colour_manual(values = okabe, labels = n_label(pc$group)) +
        labs(x = sprintf("PC1 (%.1f%%)", 100 * pc$var_PC1[1]), y = sprintf("PC2 (%.1f%%)", 100 * pc$var_PC2[1]), colour = NULL,
             title = sprintf("%s mode: PCA of %d complete features, QCs included", modes[[tg]], pc$features_used[1]))
}
save_fig(wrap_plots(panels, ncol = 2), "fig1_qc_overview", 9, 7)

## ---- Figure 2: QC RSD with and without drift correction -----------------------------

rsd_rows <- list()
for (tg in names(modes)) {
    ord <- fread(res(sprintf("%s_run_order.tsv", tg)))
    qcs <- ord[is_qc == TRUE & conditioning_qc == FALSE, sample_name]
    keep <- fread(res(sprintf("%s_feature_qc.tsv.gz", tg)))[detected_qc > 0.7, feature_ID]
    for (f in list.files(file.path(base, "results"), pattern = sprintf("^%s_height_normalised(_with_lowess|_without_lowess)?\\.tsv\\.gz$", tg))) {
        h <- fread(file.path(base, "results", f))[feature_ID %in% keep]
        v <- as.matrix(h[, qcs, with = FALSE])
        rsd <- apply(v, 1, function(x) 100 * sd(x) / mean(x))
        lowess <- if (grepl("with_lowess", f)) TRUE else if (grepl("without_lowess", f)) FALSE else
            fread(res(sprintf("%s_drift_rule.tsv", tg)))$lowess
        rsd_rows[[f]] <- data.table(mode = modes[[tg]], lowess = ifelse(lowess, "LOWESS on", "LOWESS off"), rsd = rsd)
    }
}
r <- rbindlist(rsd_rows)
med <- r[, .(median = median(rsd, na.rm = TRUE), n = .N), by = .(mode, lowess)]
p2 <- ggplot(r, aes(rsd, colour = lowess)) + stat_ecdf() + geom_vline(xintercept = 30, linetype = 2, colour = "grey40") +
    facet_wrap(~ mode) + coord_cartesian(xlim = c(0, 100)) +
    scale_colour_manual(values = c(`LOWESS on` = "#009E73", `LOWESS off` = "#000000")) +
    labs(x = "QC relative standard deviation (%)", y = "cumulative fraction of features", colour = NULL,
         caption = paste(sprintf("%s, %s: median %.1f%% over %d features", med$mode, med$lowess, med$median, med$n), collapse = "; "))
save_fig(p2, "fig2_qc_rsd", 8, 3.6)

## ---- Figure 3: volcano plots --------------------------------------------------------

vols <- list()
for (tg in names(modes)) {
    d <- fread(res(sprintf("%s_limma.tsv.gz", tg)))
    d[, lv := factor(level, levels = c(2, 3, 4))]
    cut <- if (d[q < 0.05, .N]) d[q < 0.05, max(p)] else NA_real_
    vols[[tg]] <- ggplot(d[order(-level)], aes(logFC, -log10(p), colour = lv)) + geom_point(size = 0.9, alpha = 0.8) +
        { if (!is.na(cut)) geom_hline(yintercept = -log10(cut), linetype = 2, colour = "grey40") } +
        scale_colour_manual(values = level_col, drop = FALSE, labels = c(`2` = "level 2", `3` = "level 3", `4` = "level 4")) +
        labs(x = "log2 fold change, PLWH over HC", y = "-log10 p (limma)", colour = "annotation",
             title = sprintf("%s mode: %d features tested, %d at q < 0.05", modes[[tg]], nrow(d), d[q < 0.05, .N]),
             caption = sprintf("HC n = %d, PLWH n = %d; dashed line: q = 0.05", d$n_hc[1], d$n_plwh[1]))
}
save_fig(wrap_plots(vols, ncol = 2), "fig3_volcano", 9, 4)

## ---- Figure 4: antiretroviral positive control ------------------------------------------

ac <- fread(res("antiretroviral_control.tsv"))
if (nrow(ac)) {
    long <- rbind(ac[, .(mode, drug, group = "PLWH", frac = plwh_detected / plwh_n, n = plwh_n)],
                  ac[, .(mode, drug, group = "HC", frac = hc_detected / hc_n, n = hc_n)])
    p4 <- ggplot(long, aes(drug, 100 * frac, fill = group)) + geom_col(position = position_dodge(0.8), width = 0.75) +
        facet_wrap(~ mode) + scale_fill_manual(values = okabe, labels = c(HC = sprintf("HC (n = %d)", ac$hc_n[1]), PLWH = sprintf("PLWH (n = %d)", ac$plwh_n[1]))) +
        labs(x = NULL, y = "detected before gap filling (% of group)", fill = NULL) +
        theme(axis.text.x = element_text(angle = 40, hjust = 1))
    save_fig(p4, "fig4_antiretrovirals", 8, 3.8)
}

## ---- Figure 5: pathways ----------------------------------------------------------------

ora <- fread(res("pathway_ora_hsa.tsv"))
if (nrow(ora)) {
    top <- head(ora[order(p)], 15)
    top[, label := sprintf("%s (%d of %d)", name, significant_members, background_members)]
    top[, label := factor(label, levels = rev(label))]
    p5 <- ggplot(top, aes(-log10(q), label)) + geom_col(fill = "#56B4E9") + geom_vline(xintercept = -log10(0.05), linetype = 2) +
        labs(x = "-log10 q (Fisher, one-sided, BH across pathways)", y = NULL,
             caption = sprintf("background: %d endogenous level 2 compounds with a KEGG identifier, %d of them called",
                               ora$background_size[1], ora$significant_size[1]))
    save_fig(p5, "fig5_pathways", 8, 4.5)
}

## ---- Supplementary: MS2 mirror plots ------------------------------------------------------

parse_ms2 <- function(s) {
    if (is.na(s) || !nzchar(s)) return(data.table(mz = numeric(), int = numeric()))
    v <- do.call(rbind, strsplit(strsplit(s, "|", fixed = TRUE)[[1]], ";", fixed = TRUE))
    d <- data.table(mz = as.numeric(v[, 1]), int = as.numeric(v[, 2]))
    d[, int := 100 * int / max(int)]
}
named <- rbindlist(lapply(names(modes), function(tg) fread(res(sprintf("%s_feature_qc.tsv.gz", tg)))[drug_level == 2, .(tag = tg, feature_ID)]))
for (tg in unique(named$tag)) {
    f <- file.path(base, "work", tg, "aligned_feature_table.txt")
    if (!file.exists(f)) next
    t <- fread(f, select = c("feature_ID", "annotation", "MS2", "matched_MS2", "similarity"))
    for (id in named[tag == tg, unique(feature_ID)]) {
        row <- t[feature_ID == id]
        if (!nrow(row) || is.na(row$matched_MS2)) next
        a <- parse_ms2(row$MS2); b <- parse_ms2(row$matched_MS2)
        p <- ggplot() + geom_segment(data = a, aes(mz, 0, xend = mz, yend = int), colour = "#000000") +
            geom_segment(data = b, aes(mz, 0, xend = mz, yend = -int), colour = "#D55E00") + geom_hline(yintercept = 0) +
            labs(x = "m/z", y = "relative intensity (%): measured above, library below",
                 title = sprintf("%s mode feature %d: %s, entropy similarity %.2f", modes[[tg]], id, row$annotation, row$similarity))
        for (ext in c("png", "pdf")) ggsave(file.path(fig_dir, "supplementary", sprintf("ms2_%s_%d.%s", tg, id, ext)), p, width = 6, height = 3.5, dpi = 300, bg = "white")
    }
}
