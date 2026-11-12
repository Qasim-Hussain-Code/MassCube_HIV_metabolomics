"""Download, convert and process one ion mode, one raw file at a time.

For each file, in turn: check the disk budget, download the raw file and
verify its published sha256, convert it with ThermoRawFileParser, check
the mzML, run MassCube's single-file processing on it, record its start
time, and delete the raw file and the mzML. What stays is MassCube's
per-file feature table (single_files/<name>.txt), its scan cache
(tmp/<name>.mzpkl, which gap filling reads in place of the mzML) and a
row of results/conversion_manifest.tsv. A file whose outputs already
exist is skipped, so an interrupted run resumes at the next file.

MassCube's own workflow function expects every mzML of a mode in its data
folder at once, which the disk budget does not allow (preregistration.md,
sections 4 and 15). This driver prepares Params the way that function's
preparation step does and calls process_single_file with a copy of them
for each file, as the workflow's joblib workers receive a pickled copy.

Usage:
  python process_mode.py --root ROOT --mode pos|neg [--pilot] [--limit N] [--connections K]

--pilot processes the six pilot files of preregistration.md section 16
and writes results/<mode>_pilot.tsv; --limit N stops after N files.
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time

import psutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disk_usage import tree_bytes  # noqa: E402
from fetch import download, file_hash  # noqa: E402
from mc_common import MODES, mode_params, native, read_conf  # noqa: E402

BASE = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS9264/FILES"
MANIFEST_FIELDS = [
    "file", "mode", "subgroup", "raw_bytes", "raw_sha256", "mzml_bytes", "mzml_sha256",
    "ms1_scans", "ms2_scans", "all_centroid", "polarity", "instrument_model", "instrument_serial",
    "mass_resolution", "expected_runtime", "method_file", "vial", "creation_date",
    "run_start", "first_scan_min", "last_scan_min", "masscube_ms_type", "masscube_ion_mode",
    "masscube_start_time", "features", "single_file_bytes", "mzpkl_bytes",
    "download_s", "convert_s", "masscube_s", "rss_after_mb", "peak_rss_so_far_mb",
    "project_bytes_after", "processed",
]




def check_space(root: str, conf: dict, incoming: int) -> None:
    used = tree_bytes(root)
    free = shutil.disk_usage(root).free
    ceiling = float(conf["DISK_CEILING_GB"]) * 1024**3
    floor = float(conf["FREE_FLOOR_GB"]) * 1024**3
    if used + incoming > ceiling:
        sys.exit(f"disk budget: project is {used / 1024**3:.2f} GB and the next file needs "
                 f"{incoming / 1024**3:.2f} GB, above the {conf['DISK_CEILING_GB']} GB ceiling; stopping")
    if free - incoming < floor:
        sys.exit(f"disk budget: {free / 1024**3:.2f} GB free on the drive and the next file needs "
                 f"{incoming / 1024**3:.2f} GB, leaving less than {conf['FREE_FLOOR_GB']} GB; stopping")


def inspect_mzml(path: str) -> dict:
    """Scan counts, centroid flags, polarity and times, read line by line.

    ThermoRawFileParser writes one element per line, so the cvParams that
    matter can be read without decoding any binary array.
    """
    out = {"ms1": 0, "ms2": 0, "centroid": 0, "profile": 0, "pos": 0, "neg": 0,
           "first": None, "last": None, "start": None, "spectra": 0}
    # Flags are counted once per spectrum and only inside <spectrum>, so a
    # cvParam in a header group or on a chromatogram cannot inflate them.
    inside = False
    seen: set[str] = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "<binary>" in line:
                continue
            if out["start"] is None and "startTimeStamp=" in line:
                m = re.search(r'startTimeStamp="([^"]+)"', line)
                out["start"] = m.group(1) if m else None
            if "<spectrum " in line:
                out["spectra"] += 1
                inside = True
                seen = set()
                continue
            if "</spectrum>" in line:
                inside = False
                continue
            if not inside:
                continue
            for acc, key in (("MS:1000127", "centroid"), ("MS:1000128", "profile"),
                             ("MS:1000130", "pos"), ("MS:1000129", "neg")):
                if f'accession="{acc}"' in line and key not in seen:
                    out[key] += 1
                    seen.add(key)
            if 'accession="MS:1000511"' in line:  # ms level
                level = int(re.search(r'value="(\d+)"', line).group(1))
                out["ms1" if level == 1 else "ms2"] += 1
            elif 'accession="MS:1000016"' in line:  # scan start time
                v = float(re.search(r'value="([^"]+)"', line).group(1))
                if 'unitName="second"' in line:
                    v /= 60
                out["first"] = v if out["first"] is None else min(out["first"], v)
                out["last"] = v if out["last"] is None else max(out["last"], v)
    return out


def parse_start(stamp: str) -> datetime.datetime:
    """The mzML run's startTimeStamp as an aware datetime.

    ThermoRawFileParser writes seven fractional digits
    (2022-10-13T20:05:39.4620643Z). Python's %f takes at most six, which is
    why MassCube's own get_start_time returns None for these files and its
    workflow would fall back to sample-table order; the seventh digit is
    dropped here instead.
    """
    m = re.fullmatch(r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d+))?(Z|[+-]\d{2}:?\d{2})?", stamp)
    if not m:
        raise ValueError(f"unreadable startTimeStamp {stamp!r}")
    base, frac, tz = m.groups()
    t = datetime.datetime.strptime(base, "%Y-%m-%dT%H:%M:%S")
    if frac:
        t = t.replace(microsecond=int(frac[:6].ljust(6, "0")))
    if tz in (None, "Z"):
        return t.replace(tzinfo=datetime.timezone.utc)
    offset = datetime.datetime.strptime(tz.replace(":", ""), "%z").tzinfo
    return t.replace(tzinfo=offset)


TRFP_TERMS = {
    # accession or name in ThermoRawFileParser 1.4.5's metadata JSON -> manifest column
    "MS:1000494": "instrument_model",     # Thermo Scientific instrument model
    "MS:1000529": "instrument_serial",    # instrument serial number
    "MS:1000011": "mass_resolution",      # mass resolution
    "NCIT:C69199": "creation_date",       # content creation date
    "expected runtime": "expected_runtime",
    "NCIT:C41275": "vial",
    "AFR:0002045": "method_file",         # device acquisition method
}


def trfp_metadata(path: str) -> dict:
    """The instrument and acquisition fields of ThermoRawFileParser's JSON.

    Each section is a list of CV terms; keys are matched case-insensitively
    because the serialiser's capitalisation is not part of its documentation.
    """
    with open(path, encoding="utf-8") as fh:
        meta = json.load(fh)
    out = {col: "" for col in TRFP_TERMS.values()}
    for section in meta.values():
        if not isinstance(section, list):
            continue
        for term in section:
            t = {str(k).lower(): v for k, v in term.items()}
            for key in (t.get("accession") or "", str(t.get("name") or "").lower()):
                if key in TRFP_TERMS and not out[TRFP_TERMS[key]]:
                    out[TRFP_TERMS[key]] = str(t.get("value", "")).replace("\t", " ").strip()
    return out


def pilot_files(table: list[dict], seed: int) -> list[str]:
    qcs = sorted((r for r in table if r["is_qc"] == "yes"), key=lambda r: int(re.sub(r"\D", "", r["sample_name"])))
    rng = random.Random(seed)
    pick = [qcs[0]["sample_name"], qcs[-1]["sample_name"]]
    for sub, k in (("HC", 2), ("INR", 1), ("IR", 1)):
        pool = sorted(r["sample_name"] for r in table if r["subgroup"] == sub)
        pick += rng.sample(pool, k)
    return pick


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODES))
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--connections", type=int, default=8)
    a = ap.parse_args()

    from masscube.raw_data_utils import find_raw_data_info
    from masscube.utils_functions import get_start_time
    from masscube.workflows import process_single_file

    root = a.root
    conf = read_conf(root)
    mode = MODES[a.mode]
    trfp = native(conf["TRFP"])
    work = os.path.join(root, "work", a.mode)
    raw_dir = os.path.join(root, "data", "raw", a.mode)
    mzml_dir = os.path.join(work, "data")
    meta_dir = os.path.join(work, "trfp_metadata")
    for d in (raw_dir, mzml_dir, meta_dir):
        os.makedirs(d, exist_ok=True)

    with open(os.path.join(root, "config", f"sample_table_{a.mode}.csv"), encoding="utf-8") as fh:
        table = list(csv.DictReader(fh))
    with open(os.path.join(root, "results", "mtbls9264_files.tsv"), encoding="utf-8") as fh:
        listing = {os.path.splitext(r["file"])[0]: r for r in csv.DictReader(fh, delimiter="\t")}

    params = mode_params(root, a.mode)

    names = [r["sample_name"] for r in table]
    if a.pilot:
        names = pilot_files(table, int(conf["SEED"]))
    sub_of = {r["sample_name"]: r["subgroup"] for r in table}

    manifest = os.path.join(root, "results", "conversion_manifest.tsv")
    done = set()
    if os.path.exists(manifest):
        with open(manifest, encoding="utf-8") as fh:
            done = {r["file"] for r in csv.DictReader(fh, delimiter="\t") if r["mode"] == mode}
    else:
        with open(manifest, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(MANIFEST_FIELDS) + "\n")

    me = psutil.Process()
    processed = 0
    for name in names:
        raw_name = name + ".raw"
        txt = os.path.join(params.single_file_dir, name + ".txt")
        mzpkl = os.path.join(params.tmp_file_dir, name + ".mzpkl")
        if raw_name in done and os.path.exists(txt) and os.path.exists(mzpkl):
            print(f"{name}: already processed; skipping")
            continue
        if a.limit and processed >= a.limit:
            break
        entry = listing[name]
        size = int(entry["size_bytes"])
        # The raw file, its mzML (measured on the pilot at up to a few times the
        # raw size) and MassCube's outputs must all fit before the download starts.
        check_space(root, conf, size * 4)

        raw = os.path.join(raw_dir, raw_name)
        t0 = time.monotonic()
        got = download(f"{BASE}/{raw_name}", raw, size=size, sha256=entry["sha256_published"], connections=a.connections)
        t_download = time.monotonic() - t0
        # A raw file already on disk from an interrupted run is only re-verified;
        # its time is not a download time and is recorded as NA.
        verb = "already present" if got.get("resumed") else f"downloaded in {t_download:.0f} s"
        print(f"{name}: {size / 1024**2:.0f} MB {verb}, sha256 verified")

        t0 = time.monotonic()
        cmd = [trfp, "-i", raw, "-o", mzml_dir, "-f", "1", "-m", "0", "-z", "-l", "3"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        t_convert = time.monotonic() - t0
        mzml = os.path.join(mzml_dir, name + ".mzML")
        if res.returncode != 0 or not os.path.exists(mzml):
            sys.exit(f"{name}: ThermoRawFileParser failed (exit {res.returncode}): {res.stdout[-2000:]} {res.stderr[-2000:]}")
        meta_json = os.path.join(mzml_dir, name + "-metadata.json")
        if not os.path.exists(meta_json):
            sys.exit(f"{name}: ThermoRawFileParser wrote no metadata JSON")
        acq = trfp_metadata(meta_json)
        shutil.move(meta_json, os.path.join(meta_dir, name + "-metadata.json"))

        info = inspect_mzml(mzml)
        ms_type, mc_mode, mc_centroid, _ = find_raw_data_info(mzml)
        mc_start = get_start_time(mzml)
        problems = []
        if info["ms1"] == 0:
            problems.append("no MS1 scans")
        if info["profile"] or info["centroid"] != info["spectra"]:
            problems.append(f"{info['centroid']} of {info['spectra']} spectra flagged centroid")
        want_pos = mode == "positive"
        if (info["pos"] if want_pos else info["neg"]) != info["spectra"]:
            problems.append(f"polarity: {info['pos']} positive and {info['neg']} negative of {info['spectra']} spectra")
        if mc_mode != mode:
            problems.append(f"MassCube reads the ion mode as {mc_mode}")
        if not mc_centroid:
            problems.append("MassCube does not read the file as centroid")
        if ms_type != "orbitrap":
            problems.append(f"MassCube reads the instrument type as {ms_type}")
        if info["start"] is None:
            problems.append("no startTimeStamp")
        else:
            try:
                run_start = parse_start(info["start"])
            except ValueError as e:
                problems.append(str(e))
        if problems:
            sys.exit(f"{name}: mzML check failed: " + "; ".join(problems))

        t0 = time.monotonic()
        process_single_file(mzml, copy.deepcopy(params), return_data=False)
        t_masscube = time.monotonic() - t0
        # process_single_file reports its own failures by printing and returning
        # None, so its outputs are the only reliable sign that it worked.
        if not (os.path.exists(txt) and os.path.exists(mzpkl)):
            sys.exit(f"{name}: MassCube single-file processing left no output (see the console log)")
        with open(txt, encoding="utf-8") as fh:
            n_features = sum(1 for _ in fh) - 1

        row = {
            "file": raw_name, "mode": mode, "subgroup": sub_of[name], "raw_bytes": got["size"], "raw_sha256": entry["sha256_published"],
            "mzml_bytes": os.path.getsize(mzml), "mzml_sha256": file_hash(mzml),
            "ms1_scans": info["ms1"], "ms2_scans": info["ms2"], "all_centroid": "yes", "polarity": mode,
            **acq,
            "run_start": run_start.isoformat(), "first_scan_min": f"{info['first']:.4f}", "last_scan_min": f"{info['last']:.4f}",
            # Kept to show what MassCube itself reads from the same header (None for these files).
            "masscube_ms_type": ms_type, "masscube_ion_mode": mc_mode, "masscube_start_time": mc_start.isoformat() if mc_start else "None",
            "features": n_features, "single_file_bytes": os.path.getsize(txt), "mzpkl_bytes": os.path.getsize(mzpkl),
            "download_s": "NA" if got.get("resumed") else f"{t_download:.0f}", "convert_s": f"{t_convert:.0f}", "masscube_s": f"{t_masscube:.0f}",
            "rss_after_mb": f"{me.memory_info().rss / 1024**2:.0f}",
            "peak_rss_so_far_mb": f"{getattr(me.memory_info(), 'peak_wset', me.memory_info().rss) / 1024**2:.0f}",
            "project_bytes_after": "", "processed": datetime.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z"),
        }
        # The pilot keeps its raw files until the m/z tolerance rule has been
        # applied (preregistration.md, section 7), so reprocessing needs no
        # second download; every other raw file goes as soon as it is done.
        if not a.pilot:
            os.remove(raw)
        os.remove(mzml)
        row["project_bytes_after"] = tree_bytes(root)
        with open(manifest, "a", encoding="utf-8", newline="\n") as fh:
            fh.write("\t".join(str(row[k]) for k in MANIFEST_FIELDS) + "\n")
        processed += 1
        print(f"{name}: {info['ms1']} MS1 and {info['ms2']} MS2 scans to {info['last']:.2f} min, {acq['instrument_model']}, "
              f"{n_features} features; mzML {row['mzml_bytes'] / 1024**2:.0f} MB, mzpkl {row['mzpkl_bytes'] / 1024**2:.0f} MB; "
              f"convert {t_convert:.0f} s, MassCube {t_masscube:.0f} s")

    if a.pilot:
        # Rebuilt from the manifest, so a pilot that was interrupted and
        # resumed still lists all of its files, in pilot order.
        with open(manifest, encoding="utf-8") as fh:
            rows = {r["file"]: r for r in csv.DictReader(fh, delimiter="\t") if r["mode"] == mode}
        pilot = [rows[n + ".raw"] for n in names if n + ".raw" in rows]
        out = os.path.join(root, "results", f"{a.mode}_pilot.tsv")
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(pilot)
        print(f"wrote {out}: {len(pilot)} of {len(names)} pilot files processed")


if __name__ == "__main__":
    main()
