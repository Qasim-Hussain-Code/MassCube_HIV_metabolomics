# Bits every MassCube stage shares. Params are always built by mode_params() so 04, 06 and 07
# see the same folders and parameter file, the way MassCube's own workflow would.

from __future__ import annotations

import csv
import datetime
import os
import re

import numpy as np
import pandas as pd

MODES = {"pos": "positive", "neg": "negative"}

# Internal standard, 2-amino-3-(2-chlorophenyl)propionic acid (C9H10ClNO2), 4 ppm in every
# extract per the ISA protocol. Worked out from element masses rather than copied from somewhere.
_MASS = {"C": 12.0, "H": 1.00782503207, "Cl": 34.96885268, "N": 14.0030740048, "O": 15.99491461956}
PROTON = 1.007276466812
IS_NEUTRAL = 9 * _MASS["C"] + 10 * _MASS["H"] + _MASS["Cl"] + _MASS["N"] + 2 * _MASS["O"]
IS_ION = {"positive": IS_NEUTRAL + PROTON, "negative": IS_NEUTRAL - PROTON}


def read_conf(root: str) -> dict:
    conf = {}
    with open(os.path.join(root, "project.conf"), encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r'^([A-Z_]+)=("?)(.*)\2\s*$', line.strip())
            if m:
                conf[m.group(1)] = m.group(3)
    return conf


def native(path: str) -> str:
    # project.conf stores MSYS paths (/c/Users/...); Windows Python needs C:/Users/...
    m = re.match(r"^/([a-zA-Z])/(.*)$", path)
    return f"{m.group(1).upper()}:/{m.group(2)}" if (m and os.name == "nt") else path


def mode_params(root: str, tag: str):
    # what _untargeted_metabolomics_workflow_preparation would set up, minus its demand
    # that every mzML of the mode sits in data/ at once (they never do here)
    from masscube.params import Params

    work = os.path.join(root, "work", tag)
    params = Params()
    params.project_dir = work
    params.sample_dir = os.path.join(work, "data")
    params.single_file_dir = os.path.join(work, "single_files")
    params.tmp_file_dir = os.path.join(work, "tmp")
    params.ms2_matching_dir = os.path.join(work, "ms2_matching")
    params.bpc_dir = os.path.join(work, "chromatograms")
    params.project_file_dir = os.path.join(work, "project_files")
    params.statistics_dir = os.path.join(work, "statistical_analysis")
    params.normalization_dir = os.path.join(work, "normalization_results")
    for d in (params.sample_dir, params.single_file_dir, params.tmp_file_dir, params.ms2_matching_dir,
              params.bpc_dir, params.project_file_dir, params.statistics_dir, params.normalization_dir):
        os.makedirs(d, exist_ok=True)
    # copied so work/<mode>/parameters.csv is exactly what the run used
    src = os.path.join(root, "config", f"masscube_parameters_{tag}.csv")
    dst = os.path.join(work, "parameters.csv")
    with open(src, encoding="utf-8") as a, open(dst, "w", encoding="utf-8", newline="\n") as b:
        b.write(a.read())
    params.read_parameters_from_csv(dst)
    params.output_single_file = True
    params.output_ms1_scans = True
    params.output_aligned_file = True
    return params


def set_parameter(root: str, tag: str, name: str, value: str) -> None:
    path = os.path.join(root, "config", f"masscube_parameters_{tag}.csv")
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    hit = [r for r in rows if r and r[0] == name]
    if len(hit) != 1:
        raise KeyError(f"{name} appears {len(hit)} times in {path}")
    hit[0][1] = value
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        csv.writer(fh, lineterminator="\n").writerows(rows)


def read_manifest(root: str, mode: str) -> pd.DataFrame:
    m = pd.read_csv(os.path.join(root, "results", "conversion_manifest.tsv"), sep="\t", dtype=str)
    m = m[m["mode"] == mode].copy()
    m["sample_name"] = m["file"].str.replace(r"\.raw$", "", regex=True)
    m["run_start"] = pd.to_datetime(m["run_start"], utc=True)
    return m


def read_single(path: str) -> pd.DataFrame:
    # 06 gzips these once a mode is aligned
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        path = path + ".gz"
    return pd.read_csv(path, sep="\t", low_memory=False, usecols=["m/z", "RT", "peak_height"])


def internal_standard(table: pd.DataFrame, mode: str, ppm: float = 10.0) -> dict:
    target = IS_ION[mode]
    hit = table[(table["m/z"] - target).abs() <= target * ppm * 1e-6]
    if hit.empty:
        return {"is_detected": False, "is_mz": np.nan, "is_ppm": np.nan, "is_rt": np.nan, "is_height": np.nan}
    best = hit.loc[hit["peak_height"].idxmax()]
    return {"is_detected": True, "is_mz": best["m/z"], "is_ppm": (best["m/z"] - target) / target * 1e6,
            "is_rt": best["RT"], "is_height": best["peak_height"]}


def match_top_features(tables: dict[str, pd.DataFrame], top: int = 200, mz_tol: float = 0.02,
                       rt_tol: float = 0.5) -> tuple[str, pd.DataFrame]:
    # Prereg section 7. Reference = index n // 2 by total height, the same way MassCube picks
    # its RT anchor file; only features found in every file come back.
    names = list(tables)
    totals = np.array([tables[n]["peak_height"].sum() for n in names])
    ref = names[np.argsort(totals)[len(names) // 2]]
    ref_top = tables[ref].nlargest(top, "peak_height").reset_index(drop=True)
    rows = []
    for i, f in ref_top.iterrows():
        found = {ref: (f["m/z"], f["RT"])}
        for n in names:
            if n == ref:
                continue
            t = tables[n]
            c = t[((t["m/z"] - f["m/z"]).abs() <= mz_tol) & ((t["RT"] - f["RT"]).abs() <= rt_tol)]
            if c.empty:
                break
            best = c.iloc[(c["m/z"] - f["m/z"]).abs().argmin()]
            found[n] = (best["m/z"], best["RT"])
        if len(found) == len(names):
            rows += [{"feature": i, "file": n, "mz": mz, "rt": rt} for n, (mz, rt) in found.items()]
    return ref, pd.DataFrame(rows, columns=["feature", "file", "mz", "rt"])


def run_order(root: str, tag: str) -> pd.DataFrame:
    # MassCube's own run ordering, except the times. Its reader gives None for TRFP's seven
    # fractional digits (2022-10-13T20:05:39.4620643Z), so they come from the manifest instead.
    from masscube.params import Params
    from masscube.utils_functions import label_batch_id

    mode = MODES[tag]
    p = Params()
    p.read_sample_metadata(os.path.join(root, "config", f"sample_table_{tag}.csv"))
    meta = p.sample_metadata
    man = read_manifest(root, mode).set_index("sample_name")
    meta["time"] = [man["run_start"].get(n, pd.NaT) for n in meta.iloc[:, 0]]
    meta["last_scan_min"] = [float(man["last_scan_min"].get(n, "nan")) for n in meta.iloc[:, 0]]
    meta = meta[meta["time"].notna()].sort_values(by=["time"]).reset_index(drop=True)
    meta["analytical_order"] = np.arange(len(meta))
    # label_batch_id works on numbers; seconds since the first file keep its ratio test intact.
    meta["time"] = (meta["time"] - meta["time"].min()).dt.total_seconds()
    meta = label_batch_id(meta)
    meta.index = np.arange(len(meta))
    study = ~meta["is_qc"] & ~meta["is_blank"]
    first_study = meta.index[study].min() if study.any() else len(meta)
    early_qc = meta.index[meta["is_qc"] & (meta.index < first_study)]
    # Broadhurst 2018: the QCs before the first sample are conditioning runs; all but the last stay out of RSD and drift
    meta["conditioning_qc"] = meta.index.isin(early_qc[:-1])
    return meta
