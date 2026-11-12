"""Resolve every DOI the README cites through Crossref, and its PMID through
Europe PMC, and record what each service returned with the date.

Writes results/doi_check.tsv: doi, resolved, title, the first eight authors,
volume, issue, pages and year from Crossref,
pmid from Europe PMC (blank when Europe PMC has none), and the date.

Usage: python check_dois.py --root ROOT
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request


def get_json(url: str, tries: int = 4):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last = e
        except OSError as e:
            last = e
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    text = open(os.path.join(a.root, "README.md"), encoding="utf-8").read()
    # A DOI in the README is written "doi 10.xxxx/..." or "DOI 10.xxxx/...".
    dois = sorted(set(d.rstrip(".,;)") for d in re.findall(r"\b[Dd][Oo][Ii] (10\.\d{4,9}/\S+)", text)))
    today = datetime.date.today().isoformat()
    rows = []
    for d in dois:
        cr = get_json(f"https://api.crossref.org/works/{urllib.parse.quote(d)}")
        msg = (cr or {}).get("message", {})
        year = ""
        for key in ("published-print", "published-online", "issued"):
            if msg.get(key, {}).get("date-parts"):
                year = str(msg[key]["date-parts"][0][0])
                break
        ep = get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&query=" + urllib.parse.quote(f'DOI:"{d}"'))
        hits = (ep or {}).get("resultList", {}).get("result", [])
        authors = "; ".join(f"{x.get('family', '')} {''.join(n[0] for n in x.get('given', '').replace('-', ' ').split())}".strip()
                            for x in msg.get("author", [])[:8])
        rows.append({"doi": d, "resolved": bool(msg), "title": " ".join((msg.get("title") or [""])[0].split()), "authors": authors,
                     "volume": msg.get("volume", ""), "issue": msg.get("issue", ""), "pages": msg.get("page", msg.get("article-number", "")),
                     "container": (msg.get("container-title") or [""])[0], "year": year,
                     "pmid": hits[0].get("pmid", "") if hits else "", "checked": today})
    with open(os.path.join(a.root, "results", "doi_check.tsv"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(rows[0]) + "\n")
        for r in rows:
            fh.write("\t".join(str(v).replace("\t", " ") for v in r.values()) + "\n")
    for r in rows:
        print(r["doi"], r["resolved"], r["year"], r["pmid"], r["title"][:70], sep=" | ")


if __name__ == "__main__":
    main()
