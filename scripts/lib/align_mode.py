#!/usr/bin/env python3
# One mode, all files processed: run order, RT and drift rules (prereg 5 to 7), MassCube alignment, then clean-up.
# --test aligns whatever is there into work/<mode>/test/ and touches nothing else.

from __future__ import annotations

import argparse
import gzip
import os
import pickle
import shutil
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disk_usage import tree_bytes  # noqa: E402
from mc_common import (MODES, internal_standard, match_top_features, mode_params, read_manifest,  # noqa: E402
                       read_single, run_order, set_parameter)

RT_LIMIT_MIN = 0.05
DRIFT_MIN_QC = 8
DRIFT_MAX_GAP = 20


def longest_gap(meta: pd.DataFrame, qc: pd.Series) -> int:
    run = best = 0
    for is_study, is_qc in zip(~meta["is_qc"] & ~meta["is_blank"], qc):
        if is_qc:
            run = 0
        elif is_study:
            run += 1
            best = max(best, run)
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODES))
    ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    root, tag, mode = a.root, a.mode, MODES[a.mode]
    work = os.path.join(root, "work", tag)
    out = os.path.join(work, "test") if a.test else os.path.join(root, "results")
    os.makedirs(out, exist_ok=True)

    from masscube.alignment import convert_features_to_df, feature_alignment, output_feature_table
    from masscube.feature_grouping import group_features_after_alignment
    from masscube.utils_functions import convert_signals_to_string

    table = pd.read_csv(os.path.join(root, "config", f"sample_table_{tag}.csv"), dtype=str)
    meta = run_order(root, tag)
    single = os.path.join(work, "single_files")
    missing = [n for n in table["sample_name"] if n not in set(meta["sample_name"])
               or not os.path.exists(os.path.join(single, n + ".txt"))]
    if missing and not a.test:
        sys.exit(f"{mode}: {len(missing)} of {len(table)} files not processed yet ({', '.join(missing[:8])}...)")

    # run order goes back into the sample table too
    man = read_manifest(root, mode).set_index("sample_name")
    order = pd.DataFrame({
        "injection_order": meta["analytical_order"] + 1,
        "sample_name": meta["sample_name"], "group": meta["group"], "subgroup": meta["subgroup"],
        "participant": meta["participant"], "is_qc": meta["is_qc"], "conditioning_qc": meta["conditioning_qc"],
        "run_start": [man.loc[n, "run_start"].isoformat() for n in meta["sample_name"]],
        "hours_from_first": (meta["time"] / 3600).round(4), "batch_id": meta["batch_id"],
        "last_scan_min": meta["last_scan_min"],
    })
    order.to_csv(os.path.join(out, f"{tag}_run_order.tsv"), sep="\t", index=False)
    if not a.test:
        pos = order.set_index("sample_name")
        table["injection_order"] = [str(pos.loc[n, "injection_order"]) for n in table["sample_name"]]
        table["batch"] = [str(pos.loc[n, "batch_id"]) for n in table["sample_name"]]
        table.to_csv(os.path.join(root, "config", f"sample_table_{tag}.csv"), index=False, lineterminator="\n")

    # correct_rt on if the median RT range across QCs > 0.05 min
    qc_names = list(meta.loc[meta["is_qc"] & ~meta["conditioning_qc"], "sample_name"])
    tables = {n: read_single(os.path.join(single, n + ".txt")) for n in meta["sample_name"]}
    if len(qc_names) >= 2:
        ref, m = match_top_features({n: tables[n] for n in qc_names})
        rt_range = m.groupby("feature")["rt"].agg(lambda x: x.max() - x.min())
        rt_median = float(rt_range.median()) if len(rt_range) else float("nan")
    else:
        ref, rt_range, rt_median = "", pd.Series(dtype=float), float("nan")
    correct_rt = bool(rt_median > RT_LIMIT_MIN)
    pd.DataFrame([{"mode": mode, "qcs_used": len(qc_names), "reference_qc": ref, "matched_in_all_qcs": int(rt_range.size),
                   "rt_range_median_min": rt_median, "rt_range_p95_min": float(rt_range.quantile(0.95)) if len(rt_range) else np.nan,
                   "limit_min": RT_LIMIT_MIN, "correct_rt": correct_rt}]).to_csv(os.path.join(out, f"{tag}_rt_rule.tsv"), sep="\t", index=False)

    # LOWESS only with >= 8 usable QCs, a QC either side of the study samples, and no stretch of > 20 samples without one
    usable_qc = meta["is_qc"] & ~meta["conditioning_qc"]
    study_idx = meta.index[~meta["is_qc"] & ~meta["is_blank"]]
    qc_idx = meta.index[usable_qc]
    before = bool(len(study_idx) and len(qc_idx) and qc_idx.min() < study_idx.min())
    after = bool(len(study_idx) and len(qc_idx) and qc_idx.max() > study_idx.max())
    gap = longest_gap(meta, usable_qc)
    drift = bool(len(qc_idx) >= DRIFT_MIN_QC and before and after and gap <= DRIFT_MAX_GAP)
    pd.DataFrame([{"mode": mode, "usable_qcs": len(qc_idx), "min_qcs": DRIFT_MIN_QC, "qc_before_first_sample": before,
                   "qc_after_last_sample": after, "longest_study_run_without_qc": gap, "max_run": DRIFT_MAX_GAP,
                   "lowess": drift}]).to_csv(os.path.join(out, f"{tag}_drift_rule.tsv"), sep="\t", index=False)
    if not a.test:
        set_parameter(root, tag, "correct_rt", "true" if correct_rt else "false")
        set_parameter(root, tag, "signal_normalization", "true" if drift else "false")

    params = mode_params(root, tag)
    if a.test:
        params.correct_rt = correct_rt
        params.project_file_dir = out
    meta_mc = meta.copy()
    meta_mc["VALID"] = True
    meta_mc["ABSOLUTE_PATH"] = None
    params.sample_metadata = meta_mc
    meta_mc.to_csv(os.path.join(params.project_file_dir, "sample_table_with_time.csv"), index=False)
    t0 = time.monotonic()
    features = feature_alignment(params.single_file_dir, params)
    t_align = time.monotonic() - t0
    with open(os.path.join(params.project_file_dir, "sample_metadata.pkl"), "wb") as fh:
        pickle.dump(params.sample_metadata, fh)
    ft = convert_features_to_df(features=features, sample_names=params.sample_metadata.iloc[:, 0], quant_method=params.quant_method)
    output_feature_table(ft, os.path.join(params.project_file_dir, "aligned_feature_table_before_annotation.txt"))

    # Grouping before annotation, unlike MassCube: it reads the .mzpkl caches and quietly
    # leaves features ungrouped without them, and the caches have to go before the
    # database comes down. So a library adduct can't feed in here (deviation 2).
    t0 = time.monotonic()
    group_features_after_alignment(features, params)
    for f in features:
        f.isotope_signals = convert_signals_to_string(f.isotope_signals)
    t_group = time.monotonic() - t0
    # one gzipped copy, after grouping (positive: 1.14 GB to 0.24): the 3.3 GB database
    # has to fit beside both modes' features under the 9.5 GB ceiling
    with gzip.open(os.path.join(params.project_file_dir, "aligned_features_grouped.pkl.gz"), "wb", compresslevel=6) as fh:
        pickle.dump(features, fh)

    rows = []
    for n in meta["sample_name"]:
        t = tables[n]
        rows.append({"sample_name": n, "features": len(t), "summed_height": float(t["peak_height"].sum()),
                     "max_height": float(t["peak_height"].max()),
                     "last_scan_min": float(meta.loc[meta["sample_name"] == n, "last_scan_min"].iloc[0]),
                     **internal_standard(t, mode)})
    pd.DataFrame(rows).to_csv(os.path.join(out, f"{tag}_injection_metrics.tsv"), sep="\t", index=False)

    # caches out, tables gzipped (~28 MB to ~10 each), or the other mode won't fit
    before_bytes = tree_bytes(root)
    if not a.test:
        for n in meta["sample_name"]:
            cache = os.path.join(work, "tmp", n + ".mzpkl")
            if os.path.exists(cache):
                os.remove(cache)
            txt = os.path.join(single, n + ".txt")
            if os.path.exists(txt):
                with open(txt, "rb") as src, gzip.open(txt + ".gz", "wb", compresslevel=6) as dst:
                    shutil.copyfileobj(src, dst)
                os.remove(txt)
    after_bytes = tree_bytes(root)
    n_detected = np.array([np.sum(f.feature_id_arr != -1) for f in features])
    groups = pd.Series([f.feature_group_id for f in features])
    pd.DataFrame([{"mode": mode, "files": len(meta), "aligned_features": len(features),
                   "feature_groups": int(groups.nunique()),
                   "isotopes": int(sum(f.is_isotope for f in features)),
                   "in_source_fragments": int(sum(f.is_in_source_fragment for f in features)),
                   "median_files_detected_before_gap_filling": float(np.median(n_detected)) if len(n_detected) else np.nan,
                   "correct_rt": correct_rt, "align_s": round(t_align, 1), "group_s": round(t_group, 1),
                   "project_gb_before_cleanup": round(before_bytes / 1024**3, 3),
                   "project_gb_after_cleanup": round(after_bytes / 1024**3, 3)}]).to_csv(
        os.path.join(out, f"{tag}_alignment_summary.tsv"), sep="\t", index=False)
    print(f"{mode}: {len(features)} aligned features from {len(meta)} files in {t_align:.0f} s; "
          f"retention time correction {'on' if correct_rt else 'off'}, LOWESS {'on' if drift else 'off'}; "
          f"project {before_bytes / 1024**3:.2f} GB before clean-up, {after_bytes / 1024**3:.2f} GB after")


if __name__ == "__main__":
    main()
