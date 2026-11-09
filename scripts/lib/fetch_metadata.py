"""Fetch the MTBLS9264 ISA-Tab, assignment files, checksums and file list.

Everything is read from the public FTP tree. The six metadata files are
checked against the sha256 sums MetaboLights publishes in
HASHES/metadata_sha256.json. The raw files are not downloaded here: each
one is sized with a HEAD request and matched to the published sha256 in
HASHES/data_sha256.json, and its ion mode is taken from the two ISA assay
tables, not guessed from the file name.

Usage: python fetch_metadata.py --root ROOT

Writes data/metadata/ (untracked) and, tracked:
  results/mtbls9264_files.tsv            one row per raw file: mode, bytes, sha256
  results/mtbls9264_metadata_files.tsv   the metadata files and their checksum status
  results/mtbls9264_study.tsv            study record as retrieved
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime
import json
import os
import re
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fetch import download, file_hash, remote_size  # noqa: E402

STUDY = "MTBLS9264"
BASE = f"https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{STUDY}"
WS = f"https://www.ebi.ac.uk/metabolights/ws/studies/public/study/{STUDY}"


def read_isa_table(path: str) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    # ISA headers repeat (several "Term Source REF" columns); keep the first of each name.
    keep = {}
    for i, h in enumerate(header):
        keep.setdefault(h, i)
    return [{h: (r[i] if i < len(r) else "") for h, i in keep.items()} for r in rows[1:] if any(r)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    meta = os.path.join(a.root, "data", "metadata")
    res = os.path.join(a.root, "results")
    os.makedirs(meta, exist_ok=True)
    os.makedirs(res, exist_ok=True)
    now = datetime.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")

    for f in ("metadata_sha256.json", "data_sha256.json"):
        download(f"{BASE}/HASHES/{f}", os.path.join(meta, f))
    published_meta = json.load(open(os.path.join(meta, "metadata_sha256.json"), encoding="utf-8"))
    published_data = json.load(open(os.path.join(meta, "data_sha256.json"), encoding="utf-8"))

    meta_rows = []
    for name, sha in sorted(published_meta.items()):
        r = download(f"{BASE}/{name}", os.path.join(meta, name), sha256=sha)
        meta_rows.append({"file": name, "size_bytes": r["size"], "sha256": sha, "checked_against": "HASHES/metadata_sha256.json", "retrieved": now})

    study_json = os.path.join(meta, "study.json")
    with urllib.request.urlopen(WS, timeout=300) as r, open(study_json, "wb") as out:
        out.write(r.read())
    st = json.load(open(study_json, encoding="utf-8"))["content"]
    meta_rows.append({"file": "study.json (MetaboLights web service)", "size_bytes": os.path.getsize(study_json),
                      "sha256": file_hash(study_json), "checked_against": "none published", "retrieved": now})

    # Mode per raw file, from the assay tables.
    mode_of = {}
    for name in published_meta:
        if not name.startswith("a_"):
            continue
        mode = "positive" if "_positive_" in name else "negative" if "_negative_" in name else None
        if mode is None:
            raise SystemExit(f"cannot tell the ion mode of assay table {name}")
        for row in read_isa_table(os.path.join(meta, name)):
            raw = row.get("Raw Spectral Data File", "")
            if raw:
                if raw in mode_of and mode_of[raw][0] != mode:
                    raise SystemExit(f"{raw} is listed under both ion modes")
                mode_of[raw] = (mode, row.get("Sample Name", ""))

    with urllib.request.urlopen(f"{BASE}/FILES/", timeout=300) as r:
        index = r.read().decode("utf-8", "replace")
    listed = sorted(set(re.findall(r'href="([^"/?]+\.raw)"', index, flags=re.I)))

    def size_of(name: str):
        return name, *remote_size(f"{BASE}/FILES/{name}")

    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        sizes = {n: (s, lm) for n, s, lm in pool.map(size_of, listed)}

    rows = []
    for name in listed:
        key = f"FILES/{name}"
        mode, sample = mode_of.get(key, mode_of.get(name, ("unassigned", "")))
        rows.append({"file": name, "mode": mode, "isa_sample_name": sample, "size_bytes": sizes[name][0],
                     "sha256_published": published_data.get(key, "NA"), "last_modified": sizes[name][1]})
    missing_hash = [r["file"] for r in rows if r["sha256_published"] == "NA"]
    hashed_not_listed = sorted(k for k in published_data if k.startswith("FILES/") and k[6:] not in sizes)
    in_isa_not_listed = sorted(k for k in mode_of if k.replace("FILES/", "") not in sizes)

    def write(path, recs):
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(recs)

    write(os.path.join(res, "mtbls9264_files.tsv"), rows)
    write(os.path.join(res, "mtbls9264_metadata_files.tsv"), meta_rows)
    by_mode = {}
    for r in rows:
        n, b = by_mode.get(r["mode"], (0, 0))
        by_mode[r["mode"]] = (n + 1, b + r["size_bytes"])
    study = {
        "study": STUDY, "title": st.get("title", ""), "revision": st.get("revisionNumber", ""),
        "revision_datetime": st.get("revisionDatetime", ""),
        "public_release_date": datetime.datetime.fromtimestamp(st["studyPublicReleaseDate"] / 1000, datetime.timezone.utc).date().isoformat(),
        "submission_date": datetime.datetime.fromtimestamp(st["studySubmissionDate"] / 1000, datetime.timezone.utc).date().isoformat(),
        "study_size_bytes_reported": st.get("studySize", ""),
        "raw_files": len(rows), "raw_bytes": sum(r["size_bytes"] for r in rows),
        **{f"{m}_files": by_mode[m][0] for m in sorted(by_mode)},
        **{f"{m}_bytes": by_mode[m][1] for m in sorted(by_mode)},
        "raw_without_published_sha256": len(missing_hash),
        "published_sha256_without_file": len(hashed_not_listed),
        "assay_files_not_in_listing": len(in_isa_not_listed),
        "retrieved": now,
    }
    write(os.path.join(res, "mtbls9264_study.tsv"), [study])
    for label, items in (("raw files without a published sha256", missing_hash),
                         ("published sha256 entries without a file", hashed_not_listed),
                         ("assay table entries not in the listing", in_isa_not_listed)):
        if items:
            print(f"{label}: {', '.join(items[:20])}")
    print(f"{len(rows)} raw files, {study['raw_bytes'] / 1024**3:.2f} GB; " +
          ", ".join(f"{m} {n} files {b / 1024**3:.2f} GB" for m, (n, b) in sorted(by_mode.items())))


if __name__ == "__main__":
    main()
