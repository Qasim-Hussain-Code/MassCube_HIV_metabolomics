"""Build the MassCube sample tables from the deposited ISA-Tab.

One table per ion mode, in MassCube's format: the first column is the raw
file stem (MassCube matches it to the file), then is_qc and is_blank as
"yes"/"no" strings. MassCube reads those two columns as booleans only
when they are strings; a column of Python booleans is silently replaced
by False, which would turn every QC into a study sample. The remaining
columns are for this repository: group (PLWH, HC or QC), subgroup (the
deposited disease factor: HC, INR, IR or QC), participant (the ISA
sample name, shared by the two modes), isa_listed, and injection_order
and batch, which stay empty here and are filled from the acquisition
timestamps once every file of a mode has been converted.

Labels come only from s_MTBLS9264.txt and the two assay tables. A file in
the FILES/ listing that no assay table mentions is kept only if its name
says it is a pooled QC and its polarity is then confirmed from its scans
by the conversion stage; it is marked isa_listed = no.

A file is linked to its participant through the assay table's Extract
Name, not its Sample Name. In the positive assay table as deposited
(revision 2), the Sample Name column is out of step with every other
identifier column in all 114 rows: Extract Name, MS Assay Name, Derived
Spectral Data File and Normalization Name all equal the raw file's stem,
while Sample Name holds the same 114 names in another order, the ten QCs
first. Read literally it would assign 40 files to the wrong group. The
negative table is consistent in every column. The script checks those
agreements on every run and writes them to
results/isa_assay_consistency.tsv, so the reason for the choice is
regenerated with the tables.

Usage: python build_sample_tables.py --root ROOT
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import re
import sys

SUBGROUP = {"Healthy Control": "HC", "Immune Non-Responders": "INR", "Immune Responders": "IR"}
GROUP = {"HC": "HC", "INR": "PLWH", "IR": "PLWH"}
EXPECTED = {"HC": 34, "INR": 32, "IR": 38}  # Tian et al. 2024 and the study factor


def read_isa(path: str) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    first = {}
    for i, h in enumerate(rows[0]):
        first.setdefault(h, i)
    return [{h: r[i] if i < len(r) else "" for h, i in first.items()} for r in rows[1:] if any(r)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    meta = os.path.join(a.root, "data", "metadata")

    study = read_isa(os.path.join(meta, "s_MTBLS9264.txt"))
    factor = {}
    for r in study:
        kind = r["Characteristics[Sample type]"]
        if kind == "pooled quality control sample":
            factor[r["Sample Name"]] = "QC"
        elif kind == "experimental sample":
            label = r["Factor Value[disease]"]
            if label not in SUBGROUP:
                sys.exit(f"sample {r['Sample Name']} has no usable group label ({label!r}); stopping")
            factor[r["Sample Name"]] = SUBGROUP[label]
        else:
            sys.exit(f"sample {r['Sample Name']} has an unknown sample type {kind!r}")

    listing = {}
    with open(os.path.join(a.root, "results", "mtbls9264_files.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            listing[r["file"]] = r

    def stem(path: str) -> str:
        return os.path.splitext(os.path.basename(path))[0]

    counts = []
    consistency = []
    for mode in ("positive", "negative"):
        assay = [f for f in os.listdir(meta) if f.startswith("a_") and f"_{mode}_" in f]
        if len(assay) != 1:
            sys.exit(f"expected one {mode} assay table, found {assay}")
        table = read_isa(os.path.join(meta, assay[0]))
        # Negative-mode files carry an "n" prefix that participant names do not.
        prefix = "n" if mode == "negative" else ""

        def matches(r: dict, col: str) -> bool:
            file_id = stem(r["Raw Spectral Data File"]).removeprefix(prefix)
            value = r[col]
            if col in ("MS Assay Name", "Derived Spectral Data File"):
                value = stem(value).removeprefix(prefix)
            return file_id == value

        agree = {col: sum(matches(r, col) for r in table)
                 for col in ("Extract Name", "MS Assay Name", "Derived Spectral Data File", "Normalization Name", "Sample Name")}
        regrouped = sum(1 for r in table if factor.get(r["Sample Name"]) != factor.get(r["Extract Name"]))
        consistency.append({"mode": mode, "assay_table": assay[0], "rows": len(table),
                            **{f"{c.lower().replace(' ', '_')}_matches_file": v for c, v in agree.items()},
                            "files_whose_group_differs_by_sample_name": regrouped,
                            "link_used": "Extract Name"})
        for col in ("Extract Name", "MS Assay Name", "Derived Spectral Data File", "Normalization Name"):
            if agree[col] != len(table):
                sys.exit(f"{mode}: {col} matches the raw file in only {agree[col]} of {len(table)} rows; "
                         "the Extract Name link cannot be trusted, stopping")
        rows = []
        seen = set()
        for r in table:
            raw = r["Raw Spectral Data File"].replace("FILES/", "")
            if raw not in listing:
                sys.exit(f"{raw} is in the {mode} assay table but not on the FTP listing")
            sub = factor.get(r["Extract Name"])
            if sub is None:
                sys.exit(f"{raw}: participant {r['Extract Name']} is not in the study table")
            rows.append((raw, r["Extract Name"], sub, "yes"))
            seen.add(raw)
        # Deposited but unlisted files of this mode: only pooled QCs are admitted.
        for raw, r in listing.items():
            if r["mode"] == "unassigned" and raw not in seen:
                name = stem(raw)
                if re.fullmatch(r"n?QC\d+", name) and name.startswith("n") == (mode == "negative"):
                    rows.append((raw, name.removeprefix("n"), "QC", "no"))
        rows.sort(key=lambda t: (t[2] != "QC", t[2], int(re.sub(r"\D", "", t[0]) or 0)))

        n = collections.Counter(sub for _, _, sub, _ in rows)
        for sub, want in EXPECTED.items():
            if n[sub] != want:
                sys.exit(f"{mode}: {n[sub]} {sub} samples, expected {want}; stopping")
        if n["QC"] == 0:
            print(f"{mode}: no pooled QC deposited", file=sys.stderr)

        tag = "pos" if mode == "positive" else "neg"
        out = os.path.join(a.root, "config", f"sample_table_{tag}.csv")
        with open(out, "w", encoding="utf-8", newline="\n") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["sample_name", "is_qc", "is_blank", "group", "subgroup", "participant", "isa_listed",
                        "injection_order", "batch"])
            for raw, participant, sub, listed in rows:
                w.writerow([os.path.splitext(raw)[0], "yes" if sub == "QC" else "no", "no",
                            "QC" if sub == "QC" else GROUP[sub], sub, participant, listed, "", ""])
        for sub in ("HC", "INR", "IR", "QC"):
            counts.append({"mode": mode, "subgroup": sub, "group": "QC" if sub == "QC" else GROUP[sub],
                           "n": n[sub], "not_in_assay_table": sum(1 for r in rows if r[2] == sub and r[3] == "no")})
        print(f"{mode}: {len(rows)} files, " + ", ".join(f"{k} {n[k]}" for k in ("HC", "INR", "IR", "QC")) + f" -> {out}")

    for name, recs in (("sample_counts.tsv", counts), ("isa_assay_consistency.tsv", consistency)):
        with open(os.path.join(a.root, "results", name), "w", encoding="utf-8", newline="\n") as fh:
            w = csv.DictWriter(fh, fieldnames=list(recs[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(recs)


if __name__ == "__main__":
    main()
