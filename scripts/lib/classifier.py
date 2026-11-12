"""Elastic-net classifier of PLWH against HC, exploratory and not diagnostic.

preregistration.md, section 13. Representative features of both modes,
joined by participant; drug-derived features excluded for the main model
and included for the comparison model. Repeated nested cross-validation:
5 outer folds, 5 inner, 10 repeats, stratified. Inside each outer training
fold the missing values are set to half the feature's training minimum (on
the log2 scale, the minimum minus one) and the features autoscaled; the
inner loop, run on that training fold, picks l1_ratio from 0.1, 0.5 and 0.9
and C from 12 values log-spaced from 0.001 to 100, by AUC, and the chosen
model is refitted on the whole training fold. The inner loop is
scikit-learn's LogisticRegressionCV, which searches the same grid over the
same folds as a grid search would but warm-starts along the C path; an
explicit grid search was about 27,000 saga fits for the three models and
had not finished a simulated test in 45 minutes. Each participant's out-of-fold probability is
averaged over the 10 repeats; AUC on those, with a 95 per cent interval
from 2,000 stratified bootstrap resamples. The negative control repeats the
whole procedure with the labels permuted once.

Writes results/classifier.tsv (one row per model) and
results/classifier_predictions.tsv.gz (averaged out-of-fold probabilities).

Usage: python classifier.py --root ROOT [--jobs N]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mc_common import MODES, read_conf  # noqa: E402

L1_RATIOS = [0.1, 0.5, 0.9]
CS = list(np.logspace(-3, 2, 12))


class HalfMin(BaseEstimator, TransformerMixin):
    """Half the training minimum per feature, on log2 values: minimum minus one."""

    def fit(self, x, y=None):
        self.fill_ = np.nanmin(x, axis=0) - 1
        self.fill_ = np.where(np.isfinite(self.fill_), self.fill_, 0.0)
        return self

    def transform(self, x):
        x = np.array(x, dtype=float, copy=True)
        idx = np.where(np.isnan(x))
        x[idx] = self.fill_[idx[1]]
        return x


def participant_matrix(res: str, drugs: bool) -> tuple[pd.DataFrame, pd.Series]:
    blocks, groups = [], {}
    for tag in MODES:
        s = pd.read_csv(os.path.join(res, f"{tag}_samples.tsv"), sep="\t")
        s = s[~s["excluded"].astype(bool)]
        m = pd.read_csv(os.path.join(res, f"{tag}_matrix_primary.tsv.gz"), sep="\t").set_index("feature_ID")
        qc = pd.read_csv(os.path.join(res, f"{tag}_feature_qc.tsv.gz"), sep="\t").set_index("feature_ID")
        if not drugs:
            drug = qc["drug"].notna() & (qc["drug_level"] <= 3)
            m = m.loc[~m.index.isin(qc.index[drug])]
        m = m[s["sample_name"]].T
        m.index = s.set_index("sample_name").loc[m.index, "participant"]
        m.columns = [f"{tag}_{c}" for c in m.columns]
        blocks.append(m)
        groups.update(dict(zip(s["participant"], s["group"])))
    x = pd.concat(blocks, axis=1, join="outer")
    y = pd.Series({p: int(groups[p] == "PLWH") for p in x.index})
    return x, y


def fit_fold(x: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray, r: int, k: int, seed: int) -> tuple:
    inner = StratifiedKFold(5, shuffle=True, random_state=seed + 100 + r)
    pipe = Pipeline([("impute", HalfMin()), ("scale", StandardScaler()),
                     ("clf", LogisticRegressionCV(Cs=CS, l1_ratios=L1_RATIOS, solver="saga", cv=inner, scoring="roc_auc",
                                                  max_iter=5000, n_jobs=1, refit=True, random_state=seed,
                                                  use_legacy_attributes=False))])
    pipe.fit(x[train], y[train])
    clf = pipe.named_steps["clf"]
    return r, test, pipe.predict_proba(x[test])[:, 1], {
        "repeat": r + 1, "fold": k + 1, "C": float(np.ravel(clf.C_)[0]),
        "l1_ratio": float(np.ravel(clf.l1_ratio_)[0]), "nonzero": int(np.sum(clf.coef_ != 0))}


def nested_cv(x: np.ndarray, y: np.ndarray, seed: int, jobs: int) -> tuple[np.ndarray, list[dict]]:
    # The 50 outer folds run in parallel, each with a single-threaded inner
    # search; every fold has a fixed random state, so the result does not
    # depend on how the folds are spread over processes.
    folds = [(train, test, r, k) for r in range(10)
             for k, (train, test) in enumerate(StratifiedKFold(5, shuffle=True, random_state=seed + r).split(x, y))]
    out = Parallel(n_jobs=jobs, backend="loky")(delayed(fit_fold)(x, y, tr, te, r, k, seed) for tr, te, r, k in folds)
    probs = np.zeros((10, len(y)))
    for r, test, p, _ in out:
        probs[r, test] = p
    return probs.mean(axis=0), [c for *_, c in out]


def bootstrap_ci(y: np.ndarray, p: np.ndarray, seed: int, n: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    aucs = [roc_auc_score(y[idx], p[idx]) for idx in
            (np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]) for _ in range(n))]
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--jobs", type=int, default=0)
    a = ap.parse_args()
    res = os.path.join(a.root, "results")
    conf = read_conf(a.root)
    seed = int(conf["SEED"])
    jobs = a.jobs or max(1, int(conf["THREADS"]) - 2)

    rows, preds = [], {}
    for name, drugs, permute in (("drug features excluded", False, False),
                                 ("drug features excluded, labels permuted", False, True),
                                 ("drug features included", True, False)):
        x, y = participant_matrix(res, drugs)
        yv = y.to_numpy()
        if permute:
            yv = np.random.default_rng(seed).permutation(yv)
        t0 = time.monotonic()
        p, chosen = nested_cv(x.to_numpy(float), yv, seed, jobs)
        tuned = pd.DataFrame(chosen)
        auc = roc_auc_score(yv, p)
        lo, hi = bootstrap_ci(yv, p, seed)
        rows.append({"model": name, "participants": len(yv), "plwh": int(yv.sum()), "hc": int(len(yv) - yv.sum()),
                     "features": x.shape[1], "auc": auc, "ci_low": lo, "ci_high": hi,
                     "median_C": tuned["C"].median(), "median_l1_ratio": tuned["l1_ratio"].median(),
                     "median_nonzero_coefficients": tuned["nonzero"].median(), "seconds": round(time.monotonic() - t0)})
        preds[name] = p
        print(f"{name}: AUC {auc:.3f} ({lo:.3f} to {hi:.3f}), {x.shape[1]} features, {len(yv)} participants")
    pd.DataFrame(rows).to_csv(os.path.join(res, "classifier.tsv"), sep="\t", index=False)
    pd.DataFrame({"participant": x.index, **{k: v for k, v in preds.items()}}).to_csv(
        os.path.join(res, "classifier_predictions.tsv.gz"), sep="\t", index=False)


if __name__ == "__main__":
    main()
