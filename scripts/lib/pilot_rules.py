#!/usr/bin/env python3
# Prereg sections 7 and 16: m/z rule on the pilot (widen to 0.015 Da if the p99 m/z range > 0.01, never narrow)
# and a projection for the whole mode, written before it runs. Exit 10 = widened, 05_pilot_rules.sh redoes the pilot.

from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disk_usage import tree_bytes  # noqa: E402
from mc_common import MODES, internal_standard, match_top_features, read_conf, read_manifest, read_single, set_parameter  # noqa: E402
from process_mode import pilot_files  # noqa: E402

LIMIT_DA = 0.01
WIDENED_DA = 0.015


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODES))
    ap.add_argument("--test", action="store_true", help="whatever pilot files exist, output to work/<mode>/test/")
    ap.add_argument("--reset-pilot", action="store_true",
                    help="remove the pilot's MassCube outputs and manifest rows so it can be processed again")
    a = ap.parse_args()
    root, tag, mode = a.root, a.mode, MODES[a.mode]
    conf = read_conf(root)

    if a.reset_pilot:
        with open(os.path.join(root, "config", f"sample_table_{tag}.csv"), encoding="utf-8") as fh:
            names = pilot_files(list(csv.DictReader(fh)), int(conf["SEED"]))
        for n in names:
            for p in (os.path.join(root, "work", tag, "single_files", n + ".txt"), os.path.join(root, "work", tag, "tmp", n + ".mzpkl")):
                if os.path.exists(p):
                    os.remove(p)
        path = os.path.join(root, "results", "conversion_manifest.tsv")
        with open(path, encoding="utf-8") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        drop = {n + ".raw" for n in names}
        kept = [rows[0]] + [r for r in rows[1:] if not (r[0] in drop and r[1] == mode)]
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(kept)
        print(f"{mode}: pilot outputs and {len(rows) - len(kept)} manifest rows removed for reprocessing")
        return 0
    out_dir = os.path.join(root, "work", tag, "test") if a.test else os.path.join(root, "results")
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(root, "config", f"sample_table_{tag}.csv"), encoding="utf-8") as fh:
        table = list(csv.DictReader(fh))
    names = pilot_files(table, int(conf["SEED"]))
    single = os.path.join(root, "work", tag, "single_files")
    man = read_manifest(root, mode).set_index("sample_name")
    have = [n for n in names if n in man.index and os.path.exists(os.path.join(single, n + ".txt"))]
    if len(have) < len(names) and not a.test:
        sys.exit(f"pilot incomplete: {len(have)} of {len(names)} files processed ({', '.join(sorted(set(names) - set(have)))} missing)")
    if len(have) < 2:
        sys.exit("the m/z rule needs at least two processed pilot files")

    tables = {n: read_single(os.path.join(single, n + ".txt")) for n in have}
    ref, matched = match_top_features(tables, top=200, mz_tol=0.02, rt_tol=0.5)
    spread = matched.groupby("feature")["mz"].agg(lambda x: x.max() - x.min())
    p99 = float(spread.quantile(0.99)) if len(spread) else float("nan")

    params = pd.read_csv(os.path.join(root, "config", f"masscube_parameters_{tag}.csv"), dtype=str).set_index("parameter")["value"]
    current = float(params["mz_tol_ms1"])
    widen = bool(p99 > LIMIT_DA) and current < WIDENED_DA
    decision = ("widen to 0.015 Da and reprocess the pilot" if widen
                else "already widened" if current >= WIDENED_DA else "keep 0.01 Da")
    rule = pd.DataFrame([{
        "mode": mode, "pilot_files": len(have), "reference_file": ref, "reference_features": 200,
        "matched_in_all_files": int(spread.size), "mz_range_median_da": float(spread.median()) if len(spread) else np.nan,
        "mz_range_p99_da": p99, "limit_da": LIMIT_DA, "mz_tol_ms1_before": current,
        "mz_tol_alignment_before": float(params["mz_tol_alignment"]), "decision": decision,
        "evaluated": pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ"),
    }])
    # A widened mode runs the rule twice; the first row is the decision itself, so append, don't overwrite
    path = os.path.join(out_dir, f"{tag}_mz_rule.tsv")
    if os.path.exists(path) and not a.test:
        rule = pd.concat([pd.read_csv(path, sep="\t"), rule], ignore_index=True)
    rule.to_csv(path, sep="\t", index=False)

    rows = []
    for n in have:
        rows.append({"file": n, "subgroup": man.loc[n, "subgroup"], **internal_standard(tables[n], mode)})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, f"{tag}_pilot_internal_standard.tsv"), sep="\t", index=False)

    p = man.loc[have]
    num = lambda c: pd.to_numeric(p[c], errors="coerce")
    n_mode = len(table)
    per_file_out = (num("single_file_bytes") + num("mzpkl_bytes")).mean()
    # download_s is NA for files that were already on disk
    dl = p[pd.to_numeric(p["download_s"], errors="coerce").notna()]
    raw_mb_s = (pd.to_numeric(dl["raw_bytes"]).sum() / pd.to_numeric(dl["download_s"]).clip(lower=1).sum() / 1024**2
                if len(dl) else float("nan"))
    project_now = tree_bytes(root)
    kept_raw = sum(os.path.getsize(os.path.join(root, "data", "raw", tag, n + ".raw"))
                   for n in have if os.path.exists(os.path.join(root, "data", "raw", tag, n + ".raw")))
    proj = {
        "mode": mode, "pilot_files": len(have), "mode_files": n_mode,
        "mean_download_s": num("download_s").mean(), "pilot_download_mb_per_s": raw_mb_s,
        "mean_convert_s": num("convert_s").mean(), "mean_masscube_s": num("masscube_s").mean(),
        "max_peak_rss_mb": num("peak_rss_so_far_mb").max(),
        "mean_single_file_mb": num("single_file_bytes").mean() / 1024**2, "mean_mzpkl_mb": num("mzpkl_bytes").mean() / 1024**2,
        "max_mzml_mb": num("mzml_bytes").max() / 1024**2, "max_raw_mb": num("raw_bytes").max() / 1024**2,
        "projected_processing_h": n_mode * (num("convert_s").mean() + num("masscube_s").mean()) / 3600,
        "projected_download_h_at_pilot_rate": num("raw_bytes").mean() * n_mode / max(raw_mb_s, 1e-9) / 1024**2 / 3600,
        "project_gb_now": project_now / 1024**3,
        # now, minus the pilot's raw files, plus every remaining file's outputs and one raw+mzML in flight
        "projected_peak_gb": (project_now - kept_raw + (n_mode - len(have)) * per_file_out
                              + num("raw_bytes").max() + num("mzml_bytes").max()) / 1024**3,
        "ceiling_gb": float(conf["DISK_CEILING_GB"]),
    }
    pd.DataFrame([proj]).round(4).to_csv(os.path.join(out_dir, f"{tag}_pilot_projection.tsv"), sep="\t", index=False)
    print(f"{mode}: m/z rule on {len(have)} pilot files, {spread.size} of 200 features matched, "
          f"99th percentile range {p99:.4f} Da: {decision}")
    print(f"{mode}: projected {proj['projected_processing_h']:.1f} h of conversion and MassCube, "
          f"peak project size {proj['projected_peak_gb']:.2f} GB against {conf['DISK_CEILING_GB']} GB")

    if a.test:
        return 0
    if widen:
        set_parameter(root, tag, "mz_tol_ms1", str(WIDENED_DA))
        set_parameter(root, tag, "mz_tol_alignment", str(WIDENED_DA))
        return 10
    for n in have:
        raw = os.path.join(root, "data", "raw", tag, n + ".raw")
        if os.path.exists(raw):
            os.remove(raw)
    return 0


if __name__ == "__main__":
    sys.exit(main())
