#!/usr/bin/env python3
# Both modes against the MassCube MS/MS database and the GNPS antiretroviral spectra, then
# LOWESS (if 06 switched it on) and PQN. Prereg sections 8 and 9.
# --test --library x.msp: work/pos/test/ against a local library, no downloads.

from __future__ import annotations

import argparse
import copy
import gzip
import os
import pickle
import shutil
import sys
import time

import numpy as np
import pandas as pd
import psutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disk_usage import tree_bytes  # noqa: E402
from fetch import download  # noqa: E402
from mc_common import MODES, mode_params, read_conf  # noqa: E402

GNPS_RECORD = "17232042"
GNPS_SPECTRA = "GNPS_Drug_Library_Spectra_Drugs_and_Metabolites.mgf"
SLIM_DROP = ["MS2", "matched_MS2", "peak_shape", "isotopes"]


def check_space(root: str, conf: dict, incoming: int) -> None:
    used, free = tree_bytes(root), shutil.disk_usage(root).free
    if used + incoming > float(conf["DISK_CEILING_GB"]) * 1024**3 or free - incoming < float(conf["FREE_FLOOR_GB"]) * 1024**3:
        sys.exit(f"disk budget: project {used / 1024**3:.2f} GB, {free / 1024**3:.2f} GB free, "
                 f"next download {incoming / 1024**3:.2f} GB; stopping")


def fetch_database(root: str, conf: dict) -> str:
    from record_versions import newest_database
    db = newest_database()
    recorded = pd.read_csv(os.path.join(root, "config", "tool_versions.tsv"), sep="\t").set_index("tool").loc["MassCube MS/MS database"]
    pd.DataFrame([{**db, "recorded_at_install": recorded["checksum"],
                   "same_as_recorded": recorded["checksum"] == "md5:" + db["md5"]}]).to_csv(
        os.path.join(root, "results", "database_version.tsv"), sep="\t", index=False)
    dest = os.path.join(root, "data", "db", db["file"])
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if not os.path.exists(dest):
        check_space(root, conf, int(db["size"]))
    download(f"https://zenodo.org/api/records/{db['record']}/files/{db['file']}/content", dest,
             size=int(db["size"]), md5=db["md5"], connections=4)
    return dest


def drug_library(root: str, conf: dict) -> tuple[str, pd.DataFrame]:
    # Written as MSP: MassCube 1.2.22 has no index_msp_to_pkl, and load_ms2_db indexes an MSP
    # in memory with the same unweighted entropy as the main database (deviation 3).
    # TODO(qasim): the MGF field names below (inchikey, ionmode, pepmass, adduct, spectrumid)
    # were never checked against the real file; look at its first record before trusting the subset.
    from pyteomics import mgf
    from record_versions import fetch_json

    drugs = pd.read_csv(os.path.join(root, "config", "antiretrovirals.tsv"), sep="\t")
    first = {k.split("-")[0]: d for k, d in zip(drugs["inchikey"], drugs["drug"]) if isinstance(k, str) and k != "NA"}
    rec = fetch_json(f"https://zenodo.org/api/records/{GNPS_RECORD}")
    f = [x for x in rec["files"] if x["key"] == GNPS_SPECTRA][0]
    src = os.path.join(root, "data", "db", GNPS_SPECTRA)
    if not os.path.exists(src):
        check_space(root, conf, int(f["size"]))
    download(f"https://zenodo.org/api/records/{GNPS_RECORD}/files/{GNPS_SPECTRA}/content", src,
             size=int(f["size"]), md5=f["checksum"].split(":", 1)[1], connections=4)

    out = os.path.join(root, "data", "db", "antiretrovirals_gnps.msp")
    kept = []
    with mgf.read(src, use_index=False, convert_arrays=1) as reader, open(out, "w", encoding="utf-8", newline="\n") as msp:
        for s in reader:
            p = {k.lower(): v for k, v in s["params"].items()}
            key = str(p.get("inchikey", "")).split("-")[0].upper()
            if key not in first:
                continue
            mode = str(p.get("ionmode", "")).capitalize()
            prec = p.get("pepmass", (None,))[0] if isinstance(p.get("pepmass"), tuple) else p.get("pepmass")
            msp.write(f"NAME: {first[key]}\nPRECURSORMZ: {prec}\nPRECURSORTYPE: {p.get('adduct', '')}\n"
                      f"IONMODE: {mode}\nINCHIKEY: {p.get('inchikey', '')}\nFORMULA: {p.get('formula', '')}\n"
                      f"COMMENT: GNPS Drug Library {p.get('spectrumid', p.get('title', ''))}\n"
                      f"Num Peaks: {len(s['m/z array'])}\n")
            for mz, it in zip(s["m/z array"], s["intensity array"]):
                msp.write(f"{mz:.5f}\t{it:.1f}\n")
            msp.write("\n")
            kept.append({"drug": first[key], "ion_mode": mode, "adduct": p.get("adduct", ""), "precursor_mz": prec})
    os.remove(src)  # 338 MB, and only the subset is needed
    summary = (pd.DataFrame(kept, columns=["drug", "ion_mode", "adduct", "precursor_mz"])
               .groupby(["drug", "ion_mode"]).size().rename("spectra").reset_index())
    summary = drugs[["drug", "class"]].merge(summary, on="drug", how="left").fillna({"spectra": 0})
    summary.to_csv(os.path.join(root, "results", "drug_library_subset.tsv"), sep="\t", index=False)
    return out, summary


def level_of(f) -> int:
    # prereg section 9: level 2 needs identity search, similarity > 0.7, >= 3 peaks and <= 10 ppm
    mode = f.search_mode or ""
    if mode == "identity_search" and f.similarity is not None and f.similarity > 0.7 \
            and (f.matched_peak_number or 0) >= 3 and f.matched_precursor_mz:
        ppm = abs(f.mz - float(f.matched_precursor_mz)) / float(f.matched_precursor_mz) * 1e6
        if ppm <= 10:
            return 2
    if mode in ("identity_search", "identity_search_with_rt", "fuzzy_search"):
        return 3
    return 4


def annotate_mode(root: str, tag: str, library: str, drug_msp: str | None, out: str, work_files: str) -> None:
    from masscube.alignment import convert_features_to_df, output_feature_table, output_feature_to_msp
    from masscube.annotation import annotate_aligned_features
    from masscube.normalization import sample_normalization, signal_normalization
    from masscube.visualization import plot_ms2_matching_from_feature_table

    mode = MODES[tag]
    params = mode_params(root, tag)
    params.project_file_dir = work_files
    with gzip.open(os.path.join(work_files, "aligned_features_grouped.pkl.gz"), "rb") as fh:
        features = pickle.load(fh)
    with open(os.path.join(work_files, "sample_metadata.pkl"), "rb") as fh:
        meta = pickle.load(fh).reset_index(drop=True)
    params.sample_metadata = meta
    names = list(meta.iloc[:, 0])

    # drugs on copies, identity search only, so the main annotation below starts clean
    drug_rows = []
    if drug_msp:
        p_drug = copy.copy(params)
        p_drug.ms2_library_path = drug_msp
        p_drug.fuzzy_search = False
        with_ms2 = [copy.deepcopy(f) for f in features if len(f.ms2_seq)]
        for f in annotate_aligned_features(with_ms2, p_drug):
            if f.annotation is not None:
                ppm = (f.mz - float(f.matched_precursor_mz)) / float(f.matched_precursor_mz) * 1e6 if f.matched_precursor_mz else np.nan
                drug_rows.append({"feature_ID": f.id, "m/z": round(f.mz, 4), "RT": round(f.rt, 3), "drug": f.annotation,
                                  "similarity": round(float(f.similarity), 4), "matched_peak_number": f.matched_peak_number,
                                  "precursor_ppm": round(ppm, 2), "adduct": f.matched_adduct_type, "level": level_of(f)})
    pd.DataFrame(drug_rows, columns=["feature_ID", "m/z", "RT", "drug", "similarity", "matched_peak_number",
                                     "precursor_ppm", "adduct", "level"]).to_csv(os.path.join(out, f"{tag}_drug_hits.tsv"), sep="\t", index=False)

    params.ms2_library_path = library
    t0 = time.monotonic()
    features = annotate_aligned_features(features, params)
    t_annot = time.monotonic() - t0
    with gzip.open(os.path.join(work_files, "aligned_features_annotated.pkl.gz"), "wb", compresslevel=6) as fh:
        pickle.dump(features, fh)

    ft = convert_features_to_df(features=features, sample_names=names, quant_method="peak_height")
    output_feature_table(ft, os.path.join(params.project_dir, "aligned_feature_table.txt"))
    output_feature_to_msp(ft, os.path.join(work_files, "features.msp"))
    if params.plot_ms2:
        plot_ms2_matching_from_feature_table(ft, params)

    levels = {f.id: level_of(f) for f in features}
    ann = pd.DataFrame([{
        "feature_ID": f.id, "group_ID": f.feature_group_id, "m/z": round(f.mz, 4), "RT": round(f.rt, 3),
        "adduct": f.adduct_type, "is_isotope": f.is_isotope, "is_in_source_fragment": f.is_in_source_fragment,
        "search_mode": f.search_mode, "annotation": f.annotation, "formula": f.formula, "InChIKey": f.inchikey,
        "similarity": None if f.similarity is None else round(float(f.similarity), 4),
        "matched_peak_number": f.matched_peak_number, "matched_precursor_mz": f.matched_precursor_mz,
        "precursor_ppm": (round((f.mz - float(f.matched_precursor_mz)) / float(f.matched_precursor_mz) * 1e6, 2)
                          if f.matched_precursor_mz else None),
        "level": levels[f.id], "database": f.database} for f in features])
    ann.to_csv(os.path.join(out, f"{tag}_annotation.tsv.gz"), sep="\t", index=False)

    # detected before gap filling, which is what 08's detection filters count
    det = pd.DataFrame(np.array([f.feature_id_arr != -1 for f in features], dtype=np.int8), columns=names)
    det.insert(0, "feature_ID", [f.id for f in features])
    det.to_csv(os.path.join(out, f"{tag}_detected.tsv.gz"), sep="\t", index=False)

    ft.drop(columns=[c for c in SLIM_DROP if c in ft.columns] + names).to_csv(os.path.join(out, f"{tag}_features.tsv.gz"), sep="\t", index=False)
    msp = os.path.join(work_files, "features.msp")
    with open(msp, "rb") as src, gzip.open(os.path.join(out, f"{tag}_features.msp.gz"), "wb") as dst:
        shutil.copyfileobj(src, dst)
    if os.path.getsize(os.path.join(out, f"{tag}_features.msp.gz")) > 20 * 1024**2:
        os.remove(os.path.join(out, f"{tag}_features.msp.gz"))  # too big to track; stays in work/

    # conditioning QCs out of the normalisation
    keep = ~meta["conditioning_qc"].astype(bool)
    meta_n = meta[keep].reset_index(drop=True)
    dropped = list(meta.loc[~keep, meta.columns[0]])
    def normalise(t: pd.DataFrame, lowess: bool) -> pd.DataFrame:
        t = t.drop(columns=dropped)
        if lowess:
            t = signal_normalization(t, meta_n, params.signal_norm_method)
        if params.sample_normalization:
            t = sample_normalization(t, meta_n, params.sample_norm_method)
        return t[["feature_ID"] + list(meta_n.iloc[:, 0])]

    for quant in ("peak_height", "peak_area"):
        t = convert_features_to_df(features=features, sample_names=names, quant_method=quant)
        t[["feature_ID"] + names].to_csv(os.path.join(out, f"{tag}_{quant.replace('peak_', '')}_raw.tsv.gz"), sep="\t", index=False)
        normalise(t, params.signal_normalization).to_csv(
            os.path.join(out, f"{tag}_{quant.replace('peak_', '')}_normalised.tsv.gz"), sep="\t", index=False)
    # the other LOWESS choice too, for the sensitivity analysis; LOWESS won't run on fewer than 3 QCs
    if int(meta_n["is_qc"].sum()) >= 3:
        t = convert_features_to_df(features=features, sample_names=names, quant_method="peak_height")
        normalise(t, not params.signal_normalization).to_csv(
            os.path.join(out, f"{tag}_height_normalised_{'without' if params.signal_normalization else 'with'}_lowess.tsv.gz"),
            sep="\t", index=False)
    print(f"{mode}: {len(features)} features annotated in {t_annot:.0f} s; level 2: {sum(v == 2 for v in levels.values())}, "
          f"level 3: {sum(v == 3 for v in levels.values())}; drug hits: {len(drug_rows)}; "
          f"LOWESS {'on' if params.signal_normalization else 'off'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--library", help="local MSP or pickle used instead of the downloads (--test)")
    a = ap.parse_args()
    root = a.root
    conf = read_conf(root)

    import masscube.annotation as mca

    if a.test:
        library, drug_msp, tags = a.library, a.library, ["pos"]
    else:
        for tag in MODES:
            if not os.path.exists(os.path.join(root, "logs", f"06_align_{tag}.done")):
                sys.exit(f"{MODES[tag]} mode is not aligned yet; the database is fetched only after both modes are")
        library = fetch_database(root, conf)
        drug_msp, _ = drug_library(root, conf)
        tags = list(MODES)

    me = psutil.Process()
    rss0, t0 = me.memory_info().rss, time.monotonic()
    loaded = {library: mca.load_ms2_db(library)}
    load_s, rss1 = time.monotonic() - t0, me.memory_info().rss
    if not a.test:
        # 3.3 GB, and load_once below serves both modes from memory, so it goes before any
        # output is written; a rerun fetches it again
        os.remove(library)
    es = loaded[library]
    from importlib.metadata import version
    pd.DataFrame([{"database": os.path.basename(library), "spectra": len(es.precursor_mz_array), "load_s": round(load_s, 1),
                   "rss_before_mb": round(rss0 / 1024**2), "rss_after_mb": round(rss1 / 1024**2),
                   "rss_increase_mb": round((rss1 - rss0) / 1024**2), "ms_entropy": version("ms_entropy")}]).to_csv(
        os.path.join(root, "work", "pos", "test", "database_load.tsv") if a.test else os.path.join(root, "results", "database_load.tsv"),
        sep="\t", index=False)

    original = mca.load_ms2_db

    def load_once(path):
        # annotate_aligned_features reloads the library on every call; one load serves both modes
        if path not in loaded:
            loaded[path] = original(path)
        return loaded[path]

    mca.load_ms2_db = load_once
    for tag in tags:
        work_files = os.path.join(root, "work", tag, "test" if a.test else "project_files")
        out = work_files if a.test else os.path.join(root, "results")
        annotate_mode(root, tag, library, drug_msp, out, work_files)


if __name__ == "__main__":
    main()
