"""Record the published claims of Tian et al. 2024 with their identifiers.

The five named metabolites are looked up in PubChem by the name the paper
uses (CID, InChIKey, formula, monoisotopic mass) and in KEGG (compound
identifier); the four pathways are resolved to KEGG human and reference
map identifiers through the KEGG REST API. The two headline counts are
copied from the paper. The table is written before the authors' own
metabolite assignment files are opened, so that matching claims to
features cannot be tuned to what those files contain.

Usage: python build_claims.py --out config/tian2024_claims.tsv
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import time
import urllib.parse
import urllib.request

METABOLITES = [
    # name as printed in the paper, direction in PLWH, KEGG name used for the KEGG lookup
    ("N-acetyl-L-glutamate 5-semialdehyde", "higher", "N-Acetyl-L-glutamate 5-semialdehyde"),
    ("11b-hydroxyandrost-4-ene-3,17-dione", "higher", "11beta-Hydroxyandrost-4-ene-3,17-dione"),
    ("17alpha,21-dihydroxypregnenolone", "higher", "17alpha,21-Dihydroxypregnenolone"),
    ("N-acetyl-L-phenylalanine", "lower", "N-Acetyl-L-phenylalanine"),
    ("trehalose", "lower", "alpha,alpha-Trehalose"),
]
# PubChem does not know the paper's "11b" shorthand; the spelled-out form is tried after it.
PUBCHEM_ALIASES = {"11b-hydroxyandrost-4-ene-3,17-dione": ["11beta-hydroxyandrost-4-ene-3,17-dione", "11beta-hydroxyandrostenedione"],
                   "trehalose": ["alpha,alpha-trehalose", "trehalose"]}
PATHWAYS = [
    ("phenylalanine metabolism", "00360", 0.005),
    ("arginine and proline metabolism", "00330", 0.007),
    ("phenylalanine, tyrosine and tryptophan biosynthesis", "00400", 0.0256),
    ("arginine biosynthesis", "00220", 0.0457),
]


def get(url: str, tries: int = 5) -> str:
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "MassCube_HIV_metabolomics"}), timeout=120) as r:
                return r.read().decode("utf-8")
        except OSError as e:
            last = e
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def pubchem(name: str) -> dict:
    for n in [name] + PUBCHEM_ALIASES.get(name, []):
        q = urllib.parse.quote(n, safe="")
        try:
            d = json.loads(get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/InChIKey,MolecularFormula,MonoisotopicMass/JSON", tries=2))
        except RuntimeError:
            continue
        p = d["PropertyTable"]["Properties"][0]
        return {"pubchem_query": n, "pubchem_cid": p["CID"], "inchikey": p["InChIKey"],
                "formula": p["MolecularFormula"], "monoisotopic_mass": p["MonoisotopicMass"]}
    return {"pubchem_query": name, "pubchem_cid": "NA", "inchikey": "NA", "formula": "NA", "monoisotopic_mass": "NA"}


def kegg_compound(kegg_name: str) -> str:
    # KEGG's find service rejects commas in the query (HTTP 400); searching on
    # the words and then requiring the exact name among the results is equivalent.
    query = urllib.parse.quote(kegg_name.replace(",", " "))
    hits = []
    for line in get(f"https://rest.kegg.jp/find/compound/{query}").splitlines():
        cid, names = line.split("\t", 1)
        if kegg_name.lower() in [n.strip().lower() for n in names.split(";")]:
            hits.append(cid.replace("cpd:", ""))
    return ";".join(hits) if hits else "NA"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    today = datetime.date.today().isoformat()
    hsa = {l.split("\t")[0].replace("path:", ""): l.split("\t")[1] for l in get("https://rest.kegg.jp/list/pathway/hsa").splitlines()}

    rows = []
    for name, direction, kegg_name in METABOLITES:
        rows.append({"claim_type": "metabolite", "claim": name, "published": f"{direction} in PLWH",
                     "kegg_id": kegg_compound(kegg_name), **pubchem(name), "retrieved": today})
    for name, num, p in PATHWAYS:
        rows.append({"claim_type": "pathway", "claim": name, "published": f"p = {p}",
                     "kegg_id": f"hsa{num}" if f"hsa{num}" in hsa else f"map{num} (no hsa map)",
                     "pubchem_query": "", "pubchem_cid": "", "inchikey": "", "formula": "", "monoisotopic_mass": "",
                     "retrieved": today})
    for claim, value in (("metabolites analysed", 829), ("differential metabolites", 292)):
        rows.append({"claim_type": "count", "claim": claim, "published": value, "kegg_id": "", "pubchem_query": "",
                     "pubchem_cid": "", "inchikey": "", "formula": "", "monoisotopic_mass": "", "retrieved": today})
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r["claim_type"], r["claim"], r["kegg_id"], r["pubchem_cid"], r["inchikey"], r["formula"], sep="\t")


if __name__ == "__main__":
    main()
