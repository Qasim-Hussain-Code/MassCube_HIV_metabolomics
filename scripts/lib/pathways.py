#!/usr/bin/env python3
# KEGG over-representation, exploratory (prereg 12, after Wieder et al. 2021). The background is the
# endogenous level 2 compounds that passed the filters, never all of KEGG. KEGG lookups are cached in config/.

from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mc_common import MODES  # noqa: E402

MIN_MEMBERS = 3


class Refused(Exception):
    pass


def get(url: str, tries: int = 5) -> str:
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            # 404 is "nothing found"; 400 is KEGG refusing a name it can't search, which no retry fixes
            if e.code in (400, 404):
                return ""
            if e.code == 403:
                raise Refused(url) from e
            last = e
        except OSError as e:
            last = e
        time.sleep(3 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1)
    return out


def kegg_by_name(name: str) -> tuple[str, str]:
    # only when exactly one KEGG compound carries the name
    if not isinstance(name, str) or not name.strip():
        return "", "no KEGG identifier found"
    try:
        # safe="": a slash left bare ("Cer 14:1;2O/16:0") reads to KEGG as a find option and gets a 400
        txt = get(f"https://rest.kegg.jp/find/compound/{urllib.parse.quote(name.replace(',', ' '), safe='')}")
    except Refused:
        # KEGG answers some names with 403, e.g. any with "? " in it (the database has "?" where
        # superscript locants were lost). A real block refuses this one too and stops the run.
        get("https://rest.kegg.jp/list/cpd:C00001")
        return "", "KEGG refused the name search (HTTP 403)"
    hits = [l.split("\t")[0].replace("cpd:", "") for l in txt.splitlines()
            if l and name.strip().lower() in [n.strip().lower() for n in l.split("\t", 1)[1].split(";")]]
    return (hits[0], "KEGG exact name match") if len(hits) == 1 else ("", "no KEGG identifier found")


def kegg_ids(pairs: list[tuple[str, str]], cache_path: str) -> pd.DataFrame:
    cache = pd.read_csv(cache_path, sep="\t", dtype=str).fillna("") if os.path.exists(cache_path) else \
        pd.DataFrame(columns=["inchikey", "kegg_id", "source", "retrieved"])
    known = set(cache["inchikey"])
    today = datetime.date.today().isoformat()
    # KEGG's conv/ maps PubChem SIDs, not CIDs, so the C numbers come from PubChem's synonyms (deviation 4)
    new = []
    try:
        for k, name in pairs:
            if k in known:
                continue
            txt = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{urllib.parse.quote(k)}/synonyms/JSON")
            cands = []
            if txt:
                for info in json.loads(txt).get("InformationList", {}).get("Information", []):
                    cands += [s for s in info.get("Synonym", []) if re.fullmatch(r"C\d{5}", s)]
            kid = ""
            for c in dict.fromkeys(cands):
                if get(f"https://rest.kegg.jp/list/cpd:{c}").strip():
                    kid = c
                    break
            source = "PubChem synonyms, confirmed in KEGG" if kid else ""
            if not kid:
                # synonyms are per CID, and the C number can sit on a stereo variant's CID instead
                kid, source = kegg_by_name(name)
            new.append({"inchikey": k, "kegg_id": kid, "source": source, "retrieved": today})
    finally:
        # kept when a lookup fails too, so the rerun starts where this one stopped
        if new:
            cache = pd.concat([cache, pd.DataFrame(new)], ignore_index=True)
            cache.to_csv(cache_path, sep="\t", index=False)
    return cache


def kegg_membership(root: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    links = os.path.join(root, "config", "kegg_pathway_compounds.tsv")
    names = os.path.join(root, "config", "kegg_pathways.tsv")
    if not (os.path.exists(links) and os.path.exists(names)):
        today = datetime.date.today().isoformat()
        rows = [l.split("\t") for l in get("https://rest.kegg.jp/link/compound/pathway").splitlines() if l]
        lk = pd.DataFrame({"pathway": [r[0].replace("path:", "") for r in rows], "kegg_id": [r[1].replace("cpd:", "") for r in rows]})
        lk = lk[lk["pathway"].str.startswith("map")]
        lk["retrieved"] = today
        lk.to_csv(links, sep="\t", index=False)
        ref = [l.split("\t") for l in get("https://rest.kegg.jp/list/pathway").splitlines() if l]
        hsa = {l.split("\t")[0].replace("path:", "").replace("hsa", "map") for l in get("https://rest.kegg.jp/list/pathway/hsa").splitlines() if l}
        nm = pd.DataFrame({"pathway": [r[0].replace("path:", "") for r in ref], "name": [r[1] for r in ref]})
        nm["in_hsa"] = nm["pathway"].isin(hsa)
        nm["retrieved"] = today
        nm.to_csv(names, sep="\t", index=False)
    return pd.read_csv(links, sep="\t"), pd.read_csv(names, sep="\t")


def ora(bg: pd.DataFrame, links: pd.DataFrame, names: pd.DataFrame) -> pd.DataFrame:
    ids = set(bg["kegg_id"])
    sig = set(bg.loc[bg["significant"], "kegg_id"])
    rows = []
    for pw, g in links[links["kegg_id"].isin(ids)].groupby("pathway"):
        members = set(g["kegg_id"])
        if len(members) < MIN_MEMBERS:
            continue
        a = len(members & sig); b = len(sig) - a; c = len(members) - a; d = len(ids) - len(sig) - c
        p = fisher_exact([[a, b], [c, d]], alternative="greater")[1]
        rows.append({"pathway": pw, "background_members": len(members), "significant_members": a,
                     "background_size": len(ids), "significant_size": len(sig), "p": p,
                     "members": ";".join(sorted(members & sig))})
    out = pd.DataFrame(rows, columns=["pathway", "background_members", "significant_members", "background_size", "significant_size", "p", "members"])
    out["q"] = bh(out["p"].to_numpy())
    return out.merge(names[["pathway", "name", "in_hsa"]], on="pathway", how="left").sort_values("p")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root, res = a.root, os.path.join(a.root, "results")
    parts = []
    for tag, mode in MODES.items():
        lim = pd.read_csv(os.path.join(res, f"{tag}_limma.tsv.gz"), sep="\t")
        qc = pd.read_csv(os.path.join(res, f"{tag}_feature_qc.tsv.gz"), sep="\t")[["feature_ID", "drug"]]
        lim = lim.merge(qc, on="feature_ID", how="left")
        lim = lim[(lim["level"] == 2) & lim["drug"].isna() & lim["InChIKey"].notna()].copy()
        lim["mode"] = mode
        parts.append(lim)
    comp = pd.concat(parts)
    # one per InChIKey first block across both modes, lowest q kept
    comp["block"] = comp["InChIKey"].str.split("-").str[0]
    comp = comp.sort_values("q").groupby("block").head(1)
    pairs = list(comp[["InChIKey", "annotation"]].drop_duplicates("InChIKey").itertuples(index=False, name=None))
    ids = kegg_ids(pairs, os.path.join(root, "config", "kegg_compound_ids.tsv")).set_index("inchikey")["kegg_id"]
    comp["kegg_id"] = comp["InChIKey"].map(ids).fillna("")
    comp["significant"] = comp["q"] < 0.05
    bg = comp[comp["kegg_id"] != ""]
    comp[["block", "InChIKey", "annotation", "mode", "feature_ID", "logFC", "q", "significant", "kegg_id"]].to_csv(
        os.path.join(res, "pathway_background.tsv"), sep="\t", index=False)
    links, names = kegg_membership(root)
    # human maps as the main result; all reference maps as the sensitivity version
    hsa_maps = set(names.loc[names["in_hsa"], "pathway"])
    ora(bg, links[links["pathway"].isin(hsa_maps)], names).to_csv(os.path.join(res, "pathway_ora_hsa.tsv"), sep="\t", index=False)
    ora(bg, links, names).to_csv(os.path.join(res, "pathway_ora_reference.tsv"), sep="\t", index=False)
    print(f"{len(comp)} endogenous level 2 compounds, {len(bg)} with a KEGG identifier, "
          f"{int(bg['significant'].sum())} of them called at q < 0.05")


if __name__ == "__main__":
    main()
