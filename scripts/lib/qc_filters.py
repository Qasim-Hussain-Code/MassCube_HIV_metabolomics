"""Injection QC, run order, feature filters, the positive control and PCA.

preregistration.md, sections 5, 6, 8 and 9, applied to the tables the
annotation stage wrote. Nothing here compares the two groups' feature
values; that is the statistics stage's job.

Writes, all to results/:
  injection_qc.tsv                one row per injection, both modes, with the exclusion verdict
  run_order_test.tsv              injection order of PLWH against HC, per mode
  <mode>_filter_log.tsv           each filter with its threshold, the number removed and why
  <mode>_feature_qc.tsv.gz        per feature: detection by group, QC detection, RSD, D-ratio, flags
  <mode>_representatives.tsv      the feature tested for each MassCube feature group
  antiretroviral_control.tsv      detection of each drug by group; in positive mode also the counts
                                  under the Sample Name link that section 3 rejected
  antiretroviral_control_verdict.tsv  whether the control passes, by the rule of section 9
  <mode>_pca.tsv, <mode>_pca_order.tsv
  <mode>_samples.tsv              study samples with group, order and inclusion
  <mode>_matrix_<variant>.tsv.gz  log2 values of the representatives for each analysis variant

Usage: python qc_filters.py --root ROOT
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_sample_tables import read_isa  # noqa: E402
from mc_common import MODES  # noqa: E402

DETECT_GROUP = 0.5
DETECT_QC = 0.7
RSD_MAX = 0.30
RSD_SENSITIVITY = 0.20
DRATIO_MAX = 0.50


def injection_verdicts(res: str, tag: str) -> pd.DataFrame:
    """Section 6: each criterion computed within the mode."""
    m = pd.read_csv(os.path.join(res, f"{tag}_injection_metrics.tsv"), sep="\t")
    m["mode"] = MODES[tag]
    third = lambda c: m[c] < m[c].median() / 3
    m["masscube_check"] = third("max_height")
    m["low_feature_count"] = third("features")
    m["low_summed_height"] = third("summed_height")
    m["truncated_run"] = m["last_scan_min"] < 0.9 * m["last_scan_min"].median()
    is_rate = m["is_detected"].astype(bool).mean()
    # Dropped for the mode when the standard is found in under 90 per cent of
    # files, since it would then measure the standard rather than the injection.
    m["internal_standard_used"] = is_rate >= 0.9
    m["internal_standard_missing"] = (~m["is_detected"].astype(bool)) & m["internal_standard_used"]
    crit = ["masscube_check", "low_feature_count", "low_summed_height", "truncated_run", "internal_standard_missing"]
    m["excluded"] = m[crit].any(axis=1)
    m["reasons"] = m[crit].apply(lambda r: ";".join(c for c in crit if r[c]), axis=1)
    return m


def rsd(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.std(x, ddof=1) / np.mean(x)) if len(x) > 2 and np.mean(x) > 0 else np.nan


def log2_na(a: np.ndarray) -> np.ndarray:
    a = a.astype(float)
    a[a <= 0] = np.nan
    return np.log2(a)


def filter_mode(res: str, tag: str, samples: pd.DataFrame, qcs: list[str], height: pd.DataFrame,
                detected: pd.DataFrame, rsd_max: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sections 8.1 to 8.4, in order. Returns the per-feature table and the log."""
    hc = list(samples.loc[samples["group"] == "HC", "sample_name"])
    plwh = list(samples.loc[samples["group"] == "PLWH", "sample_name"])
    study = hc + plwh
    q = pd.DataFrame({"feature_ID": height["feature_ID"]})
    det = detected.set_index("feature_ID").loc[height["feature_ID"]]
    q["detected_hc"] = det[hc].mean(axis=1).to_numpy()
    q["detected_plwh"] = det[plwh].mean(axis=1).to_numpy()
    q["detected_qc"] = det[qcs].mean(axis=1).to_numpy() if qcs else np.nan
    h = height.set_index("feature_ID")
    q["qc_rsd"] = [rsd(r) for r in h[qcs].to_numpy(float)] if len(qcs) > 2 else np.nan
    lq, ls = log2_na(h[qcs].to_numpy()), log2_na(h[study].to_numpy())
    with np.errstate(invalid="ignore", divide="ignore"):
        q["d_ratio"] = (np.nanstd(lq, axis=1, ddof=1) / np.nanstd(ls, axis=1, ddof=1)) if len(qcs) > 2 else np.nan

    steps = [
        ("detection by group", (q["detected_hc"] >= DETECT_GROUP) | (q["detected_plwh"] >= DETECT_GROUP),
         f"detected before gap filling in at least {DETECT_GROUP:.0%} of HC or of PLWH",
         "a feature seen in too few samples of either group cannot be compared; independent of the group difference"),
        ("QC detection", q["detected_qc"] > DETECT_QC, f"detected before gap filling in more than {DETECT_QC:.0%} of usable QCs",
         "precision cannot be measured for a feature the pooled QC does not show (Broadhurst et al. 2018)"),
        ("QC RSD", q["qc_rsd"] <= rsd_max, f"relative standard deviation of normalised QC heights at most {rsd_max:.0%}",
         "technical imprecision too large for the feature to be trusted"),
        ("D-ratio", q["d_ratio"] <= DRATIO_MAX, f"SD of log2 QC values over SD of log2 study values at most {DRATIO_MAX:.0%}",
         "technical variance too close to the total variance to carry biological information"),
    ]
    if not qcs or len(qcs) <= 2:
        steps = steps[:1]  # no QCs: RSD, D-ratio and QC detection are dropped (section 7)
    keep = pd.Series(True, index=q.index)
    log = []
    for name, ok, crit, why in steps:
        before = int(keep.sum())
        keep &= ok.fillna(False)
        log.append({"step": name, "features_before": before, "features_after": int(keep.sum()),
                    "removed": before - int(keep.sum()), "criterion": crit, "why": why})
    log.append({"step": "blank filter", "features_before": int(keep.sum()), "features_after": int(keep.sum()), "removed": 0,
                "criterion": "not applied", "why": "no blanks were deposited"})
    q["passes"] = keep.to_numpy()
    return q, pd.DataFrame(log)


def representatives(ann: pd.DataFrame, q: pd.DataFrame, height: pd.DataFrame, study: list[str]) -> pd.DataFrame:
    """One feature per MassCube feature group (section 8)."""
    a = ann.merge(q[["feature_ID", "passes"]], on="feature_ID")
    a = a[a["passes"] & ~a["is_isotope"].astype(bool) & ~a["is_in_source_fragment"].astype(bool)].copy()
    med = height.set_index("feature_ID")[study].replace(0, np.nan).median(axis=1)
    a["median_height"] = a["feature_ID"].map(med)
    a["group_key"] = a["group_ID"].where(a["group_ID"].notna(), -a["feature_ID"])
    rep = a.sort_values(["group_key", "median_height"], ascending=[True, False]).groupby("group_key").head(1)
    return rep.drop(columns=["group_key"])


def matrix(height: pd.DataFrame, feats: pd.Series, cols: list[str]) -> pd.DataFrame:
    h = height.set_index("feature_ID").loc[feats, cols]
    m = pd.DataFrame(log2_na(h.to_numpy()), index=h.index, columns=cols)
    return m.reset_index()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = a.root
    res = os.path.join(root, "results")
    drugs = pd.read_csv(os.path.join(root, "config", "antiretrovirals.tsv"), sep="\t")
    drug_first = {str(k).split("-")[0]: d for k, d in zip(drugs["inchikey"], drugs["drug"])}

    inj = pd.concat([injection_verdicts(res, t) for t in MODES])
    inj.to_csv(os.path.join(res, "injection_qc.tsv"), sep="\t", index=False)

    order_tests, control = [], []
    for tag, mode in MODES.items():
        order = pd.read_csv(os.path.join(res, f"{tag}_run_order.tsv"), sep="\t")
        verdict = inj[inj["mode"] == mode].set_index("sample_name")
        order["excluded"] = order["sample_name"].map(verdict["excluded"]).fillna(False).astype(bool)
        order["exclusion_reasons"] = order["sample_name"].map(verdict["reasons"]).fillna("")
        study = order[~order["is_qc"].astype(bool)]
        order[~order["is_qc"].astype(bool)][["sample_name", "participant", "group", "subgroup", "injection_order", "batch_id",
                                             "excluded", "exclusion_reasons"]].to_csv(os.path.join(res, f"{tag}_samples.tsv"), sep="\t", index=False)
        p_o, h_o = study.loc[study["group"] == "PLWH", "injection_order"], study.loc[study["group"] == "HC", "injection_order"]
        w = stats.mannwhitneyu(p_o, h_o, alternative="two-sided")
        order_tests.append({"mode": mode, "n_plwh": len(p_o), "n_hc": len(h_o), "median_order_plwh": p_o.median(),
                            "median_order_hc": h_o.median(), "U": w.statistic, "p": w.pvalue, "associated": w.pvalue < 0.05})

        height = pd.read_csv(os.path.join(res, f"{tag}_height_normalised.tsv.gz"), sep="\t")
        detected = pd.read_csv(os.path.join(res, f"{tag}_detected.tsv.gz"), sep="\t")
        ann = pd.read_csv(os.path.join(res, f"{tag}_annotation.tsv.gz"), sep="\t")
        usable_qc = list(order.loc[order["is_qc"].astype(bool) & ~order["conditioning_qc"].astype(bool), "sample_name"])
        inc = study[~study["excluded"]]

        q, log = filter_mode(res, tag, inc, usable_qc, height, detected, RSD_MAX)
        log.to_csv(os.path.join(res, f"{tag}_filter_log.tsv"), sep="\t", index=False)
        q20, _ = filter_mode(res, tag, inc, usable_qc, height, detected, RSD_SENSITIVITY)
        q["passes_rsd20"] = q20["passes"].to_numpy()
        # A feature is drug-derived when its main-database match is one of the
        # antiretrovirals, or when the GNPS drug spectra matched it; it then
        # carries the level of that drug match, not of the main annotation.
        hits = pd.read_csv(os.path.join(res, f"{tag}_drug_hits.tsv"), sep="\t").set_index("feature_ID")
        drug_by_db = ann["InChIKey"].astype(str).str.split("-").str[0].map(drug_first)
        by_lib = ann["feature_ID"].map(hits["drug"]) if len(hits) else pd.Series(np.nan, index=ann.index)
        lib_level = ann["feature_ID"].map(hits["level"]) if len(hits) else pd.Series(np.nan, index=ann.index)
        ann["drug"] = drug_by_db.where(drug_by_db.notna(), by_lib)
        ann["drug_level"] = np.where(drug_by_db.notna(), ann["level"], lib_level)
        q = q.merge(ann[["feature_ID", "level", "drug", "drug_level"]], on="feature_ID", how="left")
        q.to_csv(os.path.join(res, f"{tag}_feature_qc.tsv.gz"), sep="\t", index=False)

        all_study = list(study["sample_name"])
        inc_study = list(inc["sample_name"])
        rep = representatives(ann, q, height, inc_study)
        rep.to_csv(os.path.join(res, f"{tag}_representatives.tsv"), sep="\t", index=False)
        variants = {"primary": (height, rep, inc_study), "all_injections": (height, rep, all_study)}
        area = pd.read_csv(os.path.join(res, f"{tag}_area_normalised.tsv.gz"), sep="\t")
        variants["area"] = (area, rep, inc_study)
        rep20 = representatives(ann, q.assign(passes=q["passes_rsd20"]), height, inc_study)
        variants["rsd20"] = (height, rep20, inc_study)
        for alt in ("with", "without"):
            p = os.path.join(res, f"{tag}_height_normalised_{alt}_lowess.tsv.gz")
            if os.path.exists(p):
                variants[f"{alt}_lowess"] = (pd.read_csv(p, sep="\t"), rep, inc_study)
        if tag == "pos" and "QC11" in usable_qc:
            q11, _ = filter_mode(res, tag, inc, [c for c in usable_qc if c != "QC11"], height, detected, RSD_MAX)
            variants["without_qc11"] = (height, representatives(ann, q11, height, inc_study), inc_study)
        for name, (tab, r, cols) in variants.items():
            matrix(tab, r["feature_ID"], cols).to_csv(os.path.join(res, f"{tag}_matrix_{name}.tsv.gz"), sep="\t", index=False)

        # Positive control, section 9: detection before gap filling, by group.
        det = detected.set_index("feature_ID")
        drug_feats = ann[ann["drug"].notna() & (ann["drug_level"] <= 3)]
        for drug, g in drug_feats.groupby("drug"):
            best_level = int(g["drug_level"].min())
            hit = det.loc[g.loc[g["drug_level"] == best_level, "feature_ID"]].max(axis=0)
            control.append({"mode": mode, "drug": drug, "best_level": best_level, "features": len(g),
                            "plwh_detected": int(hit[list(study.loc[study["group"] == "PLWH", "sample_name"])].sum()),
                            "plwh_n": int((study["group"] == "PLWH").sum()),
                            "hc_detected": int(hit[list(study.loc[study["group"] == "HC", "sample_name"])].sum()),
                            "hc_n": int((study["group"] == "HC").sum())})
            if tag == "pos":
                a_pos = [f for f in os.listdir(os.path.join(root, "data", "metadata")) if f.startswith("a_") and "_positive_" in f][0]
                s = read_isa(os.path.join(root, "data", "metadata", "s_MTBLS9264.txt"))
                lab = {r["Sample Name"]: r["Factor Value[disease]"] for r in s}
                wrong = {os.path.splitext(os.path.basename(r["Raw Spectral Data File"]))[0]: lab.get(r["Sample Name"], "")
                         for r in read_isa(os.path.join(root, "data", "metadata", a_pos))}
                grp = lambda f: "HC" if wrong.get(f) == "Healthy Control" else "PLWH" if wrong.get(f, "").startswith("Immune") else "QC"
                cols = [c for c in hit.index if c in wrong]
                control[-1]["plwh_detected_sample_name_link"] = int(sum(hit[c] for c in cols if grp(c) == "PLWH"))
                control[-1]["hc_detected_sample_name_link"] = int(sum(hit[c] for c in cols if grp(c) == "HC"))

        # PCA on log2, autoscaled data with QCs, features complete across the included injections.
        cols = usable_qc + inc_study
        x = log2_na(height.set_index("feature_ID").loc[rep["feature_ID"], cols].to_numpy()).T
        x = x[:, np.all(np.isfinite(x), axis=0)]
        if x.shape[1] >= 2:
            x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
            u, sv, vt = np.linalg.svd(x - x.mean(axis=0), full_matrices=False)
            var = sv**2 / np.sum(sv**2)
            pca = pd.DataFrame({"sample_name": cols, "PC1": u[:, 0] * sv[0], "PC2": u[:, 1] * sv[1]})
            pca = pca.merge(order[["sample_name", "group", "injection_order"]], on="sample_name")
            pca["var_PC1"], pca["var_PC2"], pca["features_used"] = var[0], var[1], x.shape[1]
            pca.to_csv(os.path.join(res, f"{tag}_pca.tsv"), sep="\t", index=False)
            st = pca[pca["group"] != "QC"]
            pd.DataFrame([{"mode": mode, "pc": pc, "spearman_rho": stats.spearmanr(st[pc], st["injection_order"])[0],
                           "p": stats.spearmanr(st[pc], st["injection_order"])[1]} for pc in ("PC1", "PC2")]).to_csv(
                os.path.join(res, f"{tag}_pca_order.tsv"), sep="\t", index=False)
        print(f"{mode}: {int(inj[(inj['mode'] == mode)]['excluded'].sum())} injections excluded; "
              f"{int(q['passes'].sum())} of {len(q)} features pass the filters; {len(rep)} representatives")

    pd.DataFrame(order_tests).to_csv(os.path.join(res, "run_order_test.tsv"), sep="\t", index=False)
    c = pd.DataFrame(control)
    if len(c):
        best = c.groupby("drug").apply(lambda g: g.loc[g["plwh_detected"].idxmax()], include_groups=False).reset_index()
        passes = bool(((best["best_level"] == 2) & (best["plwh_detected"] >= 0.5 * best["plwh_n"]) & (best["hc_detected"] <= 1)).any())
    else:
        passes = False
    c.to_csv(os.path.join(res, "antiretroviral_control.tsv"), sep="\t", index=False)
    pd.DataFrame([{"control_passes": passes, "rule": "a level 2 antiretroviral detected in at least 50% of PLWH and at most one HC"}]).to_csv(
        os.path.join(res, "antiretroviral_control_verdict.tsv"), sep="\t", index=False)
    print(f"positive control {'passes' if passes else 'fails'}")


if __name__ == "__main__":
    main()
