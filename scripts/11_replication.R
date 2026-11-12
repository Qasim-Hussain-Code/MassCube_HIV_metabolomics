## Stage 11: the published claims of Tian et al. 2024, judged by the rules
## fixed in preregistration.md, section 12, and the paper's own rule rerun
## on these data.
##
## The original rule: Student t-test with equal variances at unadjusted
## p < 0.05, and OPLS-DA VIP > 1, on the annotated compounds of both modes
## combined. OPLS-DA in ropls: one predictive component, orthogonal
## components chosen by ropls, unit-variance scaling, 7-fold
## cross-validation, 1,000 permutations. VIP is descriptive.
##
## The authors' metabolite assignment files are opened here, after the
## primary analysis, as the pre-registration requires, and only to report
## whether each named metabolite appears in them.
##
## Usage: scripts/lib/in_env.sh Rscript scripts/11_replication.R
##
## Outputs:
##   results/combined_compounds.tsv     level 2 and 3 compounds of both modes, one per InChIKey first block
##   results/original_rule.tsv          per compound: Student t p, VIP, called by the original rule
##   results/oplsda_summary.tsv         R2X, R2Y, Q2 and both permutation p values
##   results/maf_summary.tsv            what the authors' assignment files contain for each named metabolite
##   results/replication_tian2024.tsv   one row per claim, with verdict and reason

suppressPackageStartupMessages({
    library("data.table")
    library("ropls")
})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), winslash = "/")
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) && args[1] == "--root") args[2] else root
res_dir <- file.path(base, "results")
cfg_dir <- file.path(base, "config")
meta_dir <- file.path(base, "data", "metadata")
modes <- c(pos = "positive", neg = "negative")
set.seed(20261109)

first_block <- function(k) sub("-.*$", "", k)

## ---- Compounds of both modes, one per InChIKey first block --------------------

lim <- rbindlist(lapply(names(modes), function(tg) {
    d <- fread(file.path(res_dir, sprintf("%s_limma.tsv.gz", tg)))
    d[, `:=`(mode = modes[[tg]], tag = tg)]
}), fill = TRUE)
comp <- lim[level <= 3 & !is.na(InChIKey) & InChIKey != ""]
comp[, block := first_block(InChIKey)]
setorder(comp, block, level, q)
comp <- comp[, .SD[1], by = block]
fwrite(comp[, .(block, InChIKey, annotation, formula, level, mode, feature_ID, logFC, q)], file.path(res_dir, "combined_compounds.tsv"), sep = "\t")

## Participant by compound matrix, joining the modes on the participant.
wide <- lapply(names(modes), function(tg) {
    s <- fread(file.path(res_dir, sprintf("%s_samples.tsv", tg)))[excluded == FALSE]
    m <- fread(file.path(res_dir, sprintf("%s_matrix_primary.tsv.gz", tg)))
    ids <- comp[tag == tg, feature_ID]
    m <- m[feature_ID %in% ids]
    x <- t(as.matrix(m[, s$sample_name, with = FALSE]))
    colnames(x) <- comp[tag == tg][match(m$feature_ID, feature_ID), block]
    rownames(x) <- s$participant
    list(x = x, group = setNames(s$group, s$participant))
})
participants <- sort(unique(unlist(lapply(wide, function(w) rownames(w$x)))))
x <- do.call(cbind, lapply(wide, function(w) {
    out <- matrix(NA_real_, length(participants), ncol(w$x), dimnames = list(participants, colnames(w$x)))
    out[rownames(w$x), ] <- w$x
    out
}))
group <- unlist(lapply(wide, function(w) w$group))
group <- group[!duplicated(names(group))][participants]

## ---- The original rule ----------------------------------------------------------

p_t <- apply(x, 2, function(v) {
    a <- v[group == "PLWH"]; b <- v[group == "HC"]
    if (sum(!is.na(a)) < 2 || sum(!is.na(b)) < 2) return(NA_real_)
    t.test(a, b, var.equal = TRUE)$p.value
})
fit <- opls(x, factor(group), predI = 1, orthoI = NA, crossvalI = 7, permI = 1000, scaleC = "standard",
            fig.pdfC = "none", info.txtC = "none")
ortho_note <- "orthogonal components chosen by ropls"
if (length(fit@vipVn) == 0) {
    ## With orthoI = NA, ropls 1.44 returns an empty object, with a message
    ## and no error, when its automatic selection judges the first predictive
    ## component not significant (it did so on simulated data with Q2 = 0.70
    ## at one orthogonal component). One orthogonal component is then fixed,
    ## and the summary says so.
    fit <- opls(x, factor(group), predI = 1, orthoI = 1, crossvalI = 7, permI = 1000, scaleC = "standard",
                fig.pdfC = "none", info.txtC = "none")
    ortho_note <- "ropls built no model with automatic selection; one orthogonal component fixed"
}
vip <- getVipVn(fit)
rule <- data.table(block = colnames(x), p_student = p_t, vip = vip[colnames(x)])
rule[, called := !is.na(p_student) & p_student < 0.05 & vip > 1]
rule <- merge(comp[, .(block, annotation, level, mode)], rule, by = "block")
fwrite(rule, file.path(res_dir, "original_rule.tsv"), sep = "\t")
sm <- fit@summaryDF
fwrite(data.table(participants = nrow(x), compounds = ncol(x), predictive = sm$pre, orthogonal = sm$ort,
                  R2X = sm$`R2X(cum)`, R2Y = sm$`R2Y(cum)`, Q2 = sm$`Q2(cum)`, pR2Y = sm$pR2Y, pQ2 = sm$pQ2,
                  permutations = 1000, cross_validation_folds = 7, note = ortho_note),
       file.path(res_dir, "oplsda_summary.tsv"), sep = "\t")

## ---- The authors' assignment files, opened only now ----------------------------

claims <- fread(file.path(cfg_dir, "tian2024_claims.tsv"))
mafs <- list.files(meta_dir, pattern = "^m_.*maf\\.tsv$", full.names = TRUE)
maf <- rbindlist(lapply(mafs, function(f) fread(f, colClasses = "character")[, maf_file := basename(f)]), fill = TRUE)
idcol <- intersect(c("metabolite_identification", "Metabolite Identification"), names(maf))[1]
maf_rows <- claims[claim_type == "metabolite", {
    hit <- if (is.na(idcol)) integer(0) else which(tolower(maf[[idcol]]) == tolower(claim))
    .(in_authors_files = length(hit) > 0, rows = length(hit),
      mz = if (length(hit) && "mass_to_charge" %in% names(maf)) paste(unique(maf$mass_to_charge[hit]), collapse = ";") else "",
      rt = if (length(hit) && "retention_time" %in% names(maf)) paste(unique(maf$retention_time[hit]), collapse = ";") else "")
}, by = claim]
fwrite(cbind(maf_rows, maf_rows_total = nrow(maf)), file.path(res_dir, "maf_summary.tsv"), sep = "\t")

## ---- Verdicts, section 12 -----------------------------------------------------------

ann_all <- rbindlist(lapply(names(modes), function(tg) {
    a <- fread(file.path(res_dir, sprintf("%s_annotation.tsv.gz", tg)))
    q <- fread(file.path(res_dir, sprintf("%s_feature_qc.tsv.gz", tg)))[, .(feature_ID, passes)]
    merge(a, q, by = "feature_ID")[, mode := modes[[tg]]]
}))
drugs <- fread(file.path(cfg_dir, "antiretrovirals.tsv"))
out <- list()
for (i in seq_len(nrow(claims))) {
    cl <- claims[i]
    row <- data.table(claim_type = cl$claim_type, claim = cl$claim, published = as.character(cl$published), detected = FALSE, feature_ID = NA_integer_,
                      best_level = NA_integer_, mode = NA_character_, direction = NA_character_, logFC = NA_real_,
                      CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, q = NA_real_, verdict = "not testable", reason = "")
    if (cl$claim_type == "metabolite") {
        hit <- lim[level <= 3 & first_block(InChIKey) == first_block(cl$inchikey)]
        setorder(hit, level, q)
        pub_up <- grepl("^higher", cl$published)
        in_maf <- maf_rows[claim == cl$claim, in_authors_files]
        if (nrow(hit)) {
            h <- hit[1]
            same <- (h$logFC > 0) == pub_up
            row[, `:=`(detected = TRUE, feature_ID = h$feature_ID, best_level = h$level, mode = h$mode, direction = ifelse(h$logFC > 0, "higher in PLWH", "lower in PLWH"),
                       logFC = h$logFC, CI_low = h$CI.L, CI_high = h$CI.R, p = h$p, q = h$q)]
            row[, verdict := fcase(h$level == 2 & h$q < 0.05 & same, "held",
                                   same & ((h$level == 3 & h$q < 0.05) | (h$p < 0.05 & h$q >= 0.05)), "partly held",
                                   default = "not supported")]
            drug_like <- !is.na(h$formula) && h$formula %in% drugs$formula
            row[, reason := sprintf("level %d in %s mode, %s, q = %.3g%s", h$level, h$mode, ifelse(same, "same direction", "opposite direction"),
                                    h$q, ifelse(drug_like, "; formula matches an antiretroviral", ""))]
        } else {
            gone <- ann_all[level <= 3 & first_block(InChIKey) == first_block(cl$inchikey)]
            iso <- ann_all[level <= 3 & formula == cl$formula]
            row[, reason := if (nrow(gone)) sprintf("annotated at level %d but removed by the filters", min(gone$level))
                            else if (nrow(iso)) sprintf("not annotated; %d level 2 or 3 features with the same formula are annotated as other isomers (%s)",
                                                       nrow(iso), paste(head(unique(iso$annotation), 3), collapse = "; "))
                            else "not annotated at level 3 or better in either mode"]
        }
        row[, reason := paste0(reason, ifelse(isTRUE(in_maf), "; listed in the authors' assignment file", "; not in the authors' assignment file"))]
    } else if (cl$claim_type == "pathway") {
        ora <- fread(file.path(res_dir, "pathway_ora_hsa.tsv"))
        pw <- ora[pathway == sub("^hsa", "map", cl$kegg_id)]
        if (nrow(pw)) {
            row[, `:=`(detected = TRUE, p = pw$p, q = pw$q,
                       verdict = fcase(pw$q < 0.05, "held", pw$p < 0.05, "partly held", default = "not supported"),
                       reason = sprintf("%d of %d background compounds in the pathway called", pw$significant_members, pw$background_members))]
        } else {
            row[, reason := "fewer than 3 background compounds in the pathway"]
        }
    } else {
        here <- if (cl$claim == "metabolites analysed") nrow(comp) else rule[called == TRUE, .N]
        ratio <- here / as.numeric(cl$published)
        row[, `:=`(detected = TRUE, verdict = fcase(abs(ratio - 1) <= 0.25, "held", abs(ratio - 1) <= 0.5, "partly held", default = "not supported"),
                   reason = sprintf("%d here against %s published (%.0f per cent)", here, cl$published, 100 * ratio))]
    }
    out[[i]] <- row
}
rep <- rbindlist(out)
fwrite(rep, file.path(res_dir, "replication_tian2024.tsv"), sep = "\t")
print(rep[, .(claim, verdict)])
