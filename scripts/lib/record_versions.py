"""Write the installed versions to config/.

python_packages.tsv lists every distribution in the environment.
tool_versions.tsv names the pieces a reader needs to reproduce the
processing: MassCube and the wheel PyPI served for it, ms_entropy (the
MassCube database was built with 1.3.4), ThermoRawFileParser, R, Python,
and the newest version of the MassCube MS/MS database on Zenodo with its
md5. The database itself is downloaded much later, by the annotation
stage, which checks it against the md5 recorded here.

Usage: python record_versions.py --root ROOT --trfp-version V --trfp-sha256 S --r-version "R ..."
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import sys
import urllib.request
from importlib.metadata import distributions, version
from pathlib import Path

ZENODO_DB_RECORD = "19121645"  # the version the brief names; newer ones are looked for below


def fetch_json(url: str, tries: int = 4) -> dict:
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except OSError as e:  # slow links here drop connections now and then
            last = e
    raise RuntimeError(f"could not fetch {url}: {last}")


def newest_database() -> dict:
    rec = fetch_json(f"https://zenodo.org/api/records/{ZENODO_DB_RECORD}")
    versions = fetch_json(f"https://zenodo.org/api/records/{ZENODO_DB_RECORD}/versions?sort=version&size=25")
    hits = sorted(versions["hits"]["hits"], key=lambda h: h["metadata"]["publication_date"])
    newest = hits[-1]
    pkl = [f for f in newest["files"] if f["key"].endswith(".pkl")]
    if len(pkl) != 1:
        raise RuntimeError(f"expected one .pkl in Zenodo record {newest['id']}, found {len(pkl)}")
    return {
        "record": str(newest["id"]),
        "concept": str(rec["conceptrecid"]),
        "date": newest["metadata"]["publication_date"],
        "index": f"{len(hits)} of {len(hits)}",
        "file": pkl[0]["key"],
        "size": pkl[0]["size"],
        "md5": pkl[0]["checksum"].split(":", 1)[1],
        "license": newest["metadata"].get("license", {}).get("id", ""),
    }


def masscube_wheel_sha256(v: str) -> str:
    d = fetch_json(f"https://pypi.org/pypi/masscube/{v}/json")
    wheels = [u for u in d["urls"] if u["filename"].endswith(".whl")]
    return wheels[0]["digests"]["sha256"] if wheels else "NA"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--trfp-version", required=True)
    ap.add_argument("--trfp-source", required=True)
    ap.add_argument("--trfp-sha256", required=True)
    ap.add_argument("--r-version", required=True)
    a = ap.parse_args()
    cfg = Path(a.root) / "config"
    now = datetime.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")

    rows = sorted({(d.metadata["Name"].lower(), d.version) for d in distributions()})
    with open(cfg / "python_packages.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("package\tversion\trecorded\n")
        for name, v in rows:
            fh.write(f"{name}\t{v}\t{now}\n")

    mc = version("masscube")
    db = newest_database()
    tools = [
        ("masscube", mc, "PyPI wheel", f"sha256:{masscube_wheel_sha256(mc)}", "CC BY-NC 4.0"),
        ("ms_entropy", version("ms_entropy"), "PyPI, pinned by masscube", "NA", "the database was built with 1.3.4"),
        ("ThermoRawFileParser", a.trfp_version, a.trfp_source, a.trfp_sha256, "vendor peak picking on"),
        ("R", a.r_version, "conda-forge r-base", "NA", ""),
        ("Python", platform.python_version(), "conda-forge python", "NA", sys.platform),
        ("MassCube MS/MS database", db["date"],
         f"Zenodo record {db['record']} (concept {db['concept']}, version {db['index']}), {db['file']}, {db['size']} bytes",
         f"md5:{db['md5']}", db["license"]),
    ]
    with open(cfg / "tool_versions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("tool\tversion\tsource\tchecksum\tnote\trecorded\n")
        for t in tools:
            fh.write("\t".join(t) + f"\t{now}\n")
    print(f"wrote {cfg / 'python_packages.tsv'} ({len(rows)} packages) and {cfg / 'tool_versions.tsv'}")


if __name__ == "__main__":
    main()
