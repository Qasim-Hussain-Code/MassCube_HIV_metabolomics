#!/usr/bin/env python3
# config/antiretrovirals.tsv: the prereg section 9 drug list with PubChem ids; 07 picks GNPS spectra by InChIKey.

from __future__ import annotations

import argparse
import csv
import datetime
import json
import time
import urllib.parse
import urllib.request

# Chinese ART regimens of the study period; azvudine and ainuovirine are China-only approvals
DRUGS = [
    ("tenofovir", "NRTI"), ("tenofovir disoproxil", "NRTI prodrug"), ("tenofovir alafenamide", "NRTI prodrug"),
    ("emtricitabine", "NRTI"), ("lamivudine", "NRTI"), ("zidovudine", "NRTI"), ("abacavir", "NRTI"),
    ("azvudine", "NRTI"), ("efavirenz", "NNRTI"), ("rilpivirine", "NNRTI"), ("nevirapine", "NNRTI"),
    ("ainuovirine", "NNRTI"), ("dolutegravir", "integrase inhibitor"), ("bictegravir", "integrase inhibitor"),
    ("raltegravir", "integrase inhibitor"), ("elvitegravir", "integrase inhibitor"),
    ("lopinavir", "protease inhibitor"), ("ritonavir", "protease inhibitor"), ("darunavir", "protease inhibitor"),
    ("cobicistat", "pharmacokinetic booster"),
]


def pubchem(name: str) -> dict:
    q = urllib.parse.quote(name, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/InChIKey,MolecularFormula,MonoisotopicMass/JSON"
    last = None
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                p = json.load(r)["PropertyTable"]["Properties"][0]
            return {"pubchem_cid": p["CID"], "inchikey": p["InChIKey"], "formula": p["MolecularFormula"],
                    "monoisotopic_mass": p["MonoisotopicMass"]}
        except OSError as e:
            last = e
            time.sleep(3 * (i + 1))
    print(f"{name}: not found in PubChem ({last})")
    return {"pubchem_cid": "NA", "inchikey": "NA", "formula": "NA", "monoisotopic_mass": "NA"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    today = datetime.date.today().isoformat()
    rows = [{"drug": n, "class": c, **pubchem(n), "retrieved": today} for n, c in DRUGS]
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["drug"], r["pubchem_cid"], r["inchikey"], r["formula"], sep="\t")


if __name__ == "__main__":
    main()
