## Print the numbers the README quotes, each from the results file that
## holds it, so that a reader can regenerate them and check the README
## against them.
##
## Usage: scripts/lib/in_env.sh Rscript scripts/lib/readme_numbers.R

suppressPackageStartupMessages(library("data.table"))
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."), winslash = "/")
res <- function(f) file.path(root, "results", f)
show <- function(f, cols = NULL) {
    cat("\n==== results/", f, " ====\n", sep = "")
    if (!file.exists(res(f))) { cat("(not written yet)\n"); return(invisible()) }
    d <- fread(res(f))
    if (!is.null(cols)) d <- d[, intersect(cols, names(d)), with = FALSE]
    print(d, nrows = 200, class = FALSE)
}

show("mtbls9264_study.tsv", c("raw_files", "raw_bytes", "positive_files", "positive_bytes", "negative_files", "negative_bytes",
                              "unassigned_files", "unassigned_bytes", "study_size_bytes_reported"))
show("sample_counts.tsv")
show("isa_assay_consistency.tsv")
for (f in c("pos_mz_rule.tsv", "neg_mz_rule.tsv", "pos_pilot_projection.tsv", "neg_pilot_projection.tsv",
            "pos_rt_rule.tsv", "neg_rt_rule.tsv", "pos_drift_rule.tsv", "neg_drift_rule.tsv",
            "pos_alignment_summary.tsv", "neg_alignment_summary.tsv")) show(f)
show("injection_qc.tsv", c("mode", "sample_name", "excluded", "reasons"))
show("run_order_test.tsv")
for (f in c("pos_filter_log.tsv", "neg_filter_log.tsv")) show(f, c("step", "features_before", "features_after", "removed"))
show("database_version.tsv"); show("database_load.tsv"); show("drug_library_subset.tsv")
show("antiretroviral_control.tsv"); show("antiretroviral_control_verdict.tsv")
show("test_counts.tsv"); show("sensitivity_summary.tsv"); show("power.tsv")
for (f in c("pos_pca_order.tsv", "neg_pca_order.tsv")) show(f)
show("permanova.tsv"); show("classifier.tsv")
show("pathway_ora_hsa.tsv", c("pathway", "name", "background_members", "significant_members", "p", "q"))
show("oplsda_summary.tsv"); show("replication_tian2024.tsv"); show("maf_summary.tsv")

cat("\n==== conversion_manifest.tsv ====\n")
if (file.exists(res("conversion_manifest.tsv"))) {
    m <- fread(res("conversion_manifest.tsv"))
    print(m[, .(files = .N, raw_gb = sum(raw_bytes) / 1024^3, mzml_gb = sum(mzml_bytes) / 1024^3,
                median_last_scan_min = median(last_scan_min), median_features = median(features),
                instrument = paste(unique(instrument_model), collapse = "; ")), by = mode])
}
cat("\n==== logs/*.resources.tsv ====\n")
r <- rbindlist(lapply(list.files(file.path(root, "logs"), pattern = "resources\\.tsv$", full.names = TRUE), fread), fill = TRUE)
print(r[, .(stage, elapsed_s, peak_rss_mb, disk_before_gb, disk_after_gb)])
