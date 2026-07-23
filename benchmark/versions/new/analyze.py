#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ocena poz: symetryczny RMSD (RDKit) docked vs natywna poza krystaliczna.

Roznica wobec wersji legacy: RMSD liczy RDKit (rdMolAlign.CalcRMS), a pozy z
out.pdbqt odtwarza Meeko (RDKitMolCreate) razem z rzedami wiazan zapisanymi
w "REMARK SMILES" - czyli graf czasteczki jest odtworzony, a nie zgadywany
przez perception OpenBabela. CalcRMS uwzglednia symetrie czasteczki i NIE
nasuwa czasteczek na siebie (redocking: obie sa w tym samym ukladzie wspolrzednych).

Metryka pierwszorzedna: RMSD, prog sukcesu <= config.RMSD_SUCCESS, dla top-1 i best-of-N.
Score Vina raportowany tylko diagnostycznie (nieporownywalny miedzy warunkami,
bo bias sztucznie zmienia energie w jednym miejscu siatki).

Wyniki:
    results/results.csv    (per przebieg: cel, warunek, seed, RMSD top-1/best, score)
    results/summary.csv    (per cel x warunek: success@2A, mediany RMSD, Wilcoxon)
    results/plot.png       (success@2A i rozklad RMSD)
"""
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

import config as C

RDLogger.DisableLog("rdApp.warning")


def _neutral_copy(mol):
    """Kopia bez ladunkow formalnych - awaryjne dopasowanie grafu, gdy stany
    protonacji pozy i referencji zapisane sa inaczej."""
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
    m = rw.GetMol()
    Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    return m


def _rmsd_by_matches(probe, conf_id, ref):
    """Symetryczny RMSD po wszystkich dopasowaniach grafu (bez nasuwania)."""
    matches = _neutral_copy(ref).GetSubstructMatches(_neutral_copy(probe),
                                                     uniquify=False, maxMatches=100000)
    if not matches:
        raise RuntimeError("poza nie pasuje grafem do natywnego liganda")
    p = probe.GetConformer(conf_id).GetPositions()
    r = ref.GetConformer().GetPositions()
    return min(float(np.sqrt(((p - r[list(m)]) ** 2).sum(axis=1).mean())) for m in matches)


def rmsd_poses(native_sdf, out_pdbqt):
    """Symetryczny heavy-atom RMSD kazdej pozy vs natywna."""
    ref = Chem.RemoveAllHs(Chem.MolFromMolFile(native_sdf, removeHs=False))
    pmol = PDBQTMolecule.from_file(out_pdbqt, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pmol)
    probe = Chem.RemoveAllHs(mols[0])

    vals = []
    for conf in probe.GetConformers():
        try:
            vals.append(rdMolAlign.CalcRMS(probe, ref, prbId=conf.GetId()))
        except Exception:
            vals.append(_rmsd_by_matches(probe, conf.GetId(), ref))
    return vals


def vina_scores(pdbqt):
    """Lista affinity (kcal/mol) dla kolejnych poz z out.pdbqt."""
    scores = []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split(":")[1].split()[0]))
    return scores


# --------------------------------------------------------------------------
def collect():
    rows = []
    for pdbid in C.TARGETS:
        tdir = os.path.join(C.RESULTS, pdbid)
        native = os.path.join(tdir, "native.sdf")
        for cond in C.CONDITIONS:
            for seed in C.SEEDS:
                out_pdbqt = os.path.join(tdir, cond, "seed_%d" % seed, "out.pdbqt")
                if not os.path.isfile(out_pdbqt):
                    continue
                rmsds = rmsd_poses(native, out_pdbqt)
                scores = vina_scores(out_pdbqt)
                if not rmsds:
                    continue
                top1, best = rmsds[0], min(rmsds)
                rows.append(dict(
                    pdbid=pdbid, condition=cond, seed=seed,
                    top1_rmsd=round(top1, 3), best_rmsd=round(best, 3),
                    best_pose=rmsds.index(best) + 1,
                    top1_score=scores[0] if scores else "",
                    top1_success=int(top1 <= C.RMSD_SUCCESS),
                    best_success=int(best <= C.RMSD_SUCCESS),
                ))
    return pd.DataFrame(rows)


def write_runs(df):
    path = os.path.join(C.RESULTS, "results.csv")
    df.to_csv(path, index=False)
    print("zapisano " + path)


def summarize(df):
    try:
        from scipy.stats import wilcoxon
    except Exception:
        wilcoxon = None

    summary = []
    for pdbid in C.TARGETS:
        # sparowany Wilcoxon: top-1 RMSD conventional vs biased (po seedach)
        paired = {c: df[(df.pdbid == pdbid) & (df.condition == c)].set_index("seed")["top1_rmsd"]
                  for c in C.CONDITIONS}
        common = sorted(set(paired["conventional"].index) & set(paired["biased"].index))
        p_val = ""
        if wilcoxon and common:
            a = paired["conventional"].loc[common].values
            b = paired["biased"].loc[common].values
            if np.any(a != b):
                try:
                    p_val = round(float(wilcoxon(a, b).pvalue), 4)
                except Exception:
                    p_val = ""
        for cond in C.CONDITIONS:
            g = df[(df.pdbid == pdbid) & (df.condition == cond)]
            if g.empty:
                continue
            summary.append(dict(
                pdbid=pdbid, condition=cond, n=len(g),
                top1_success_rate=round(g.top1_success.mean(), 3),
                best_success_rate=round(g.best_success.mean(), 3),
                top1_rmsd_mean=round(g.top1_rmsd.mean(), 3),
                top1_rmsd_median=round(g.top1_rmsd.median(), 3),
                best_rmsd_median=round(g.best_rmsd.median(), 3),
                wilcoxon_p_top1=p_val if cond == "biased" else "",
            ))
    sdf = pd.DataFrame(summary)
    path = os.path.join(C.RESULTS, "summary.csv")
    sdf.to_csv(path, index=False)
    print("zapisano " + path)
    return sdf


def plot(df, sdf):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        print("matplotlib niedostepny - pomijam wykres")
        return
    targets = list(C.TARGETS)
    colors = {"conventional": "#5B8FF9", "biased": "#E8684A"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # lewy: success@2A (top-1) per cel x warunek
    w = 0.38
    for i, cond in enumerate(C.CONDITIONS):
        vals = []
        for t in targets:
            s = sdf[(sdf.pdbid == t) & (sdf.condition == cond)]
            vals.append(float(s.top1_success_rate.iloc[0]) if len(s) else 0.0)
        x = [j + (i - 0.5) * w for j in range(len(targets))]
        ax1.bar(x, vals, w, label=cond, color=colors[cond])
    ax1.set_xticks(range(len(targets)))
    ax1.set_xticklabels(targets)
    ax1.set_ylabel("success@%.1f A (top-1)" % C.RMSD_SUCCESS)
    ax1.set_ylim(0, 1.05)
    ax1.legend(frameon=False)
    ax1.set_title("Odzyskanie pozy (top-1)")

    # prawy: rozklad top-1 RMSD (punkty) per cel x warunek
    for i, cond in enumerate(C.CONDITIONS):
        for j, t in enumerate(targets):
            ys = df[(df.pdbid == t) & (df.condition == cond)].top1_rmsd.values
            ax2.scatter([j + (i - 0.5) * w] * len(ys), ys, s=18, alpha=0.6,
                        color=colors[cond], label=cond if j == 0 else None)
    ax2.axhline(C.RMSD_SUCCESS, ls="--", lw=1, color="gray")
    ax2.set_xticks(range(len(targets)))
    ax2.set_xticklabels(targets)
    ax2.set_ylabel("top-1 RMSD [A]")
    ax2.legend(frameon=False)
    ax2.set_title("Rozklad top-1 RMSD")

    fig.tight_layout()
    path = os.path.join(C.RESULTS, "plot.png")
    fig.savefig(path, dpi=130)
    print("zapisano " + path)


def main():
    df = collect()
    if df.empty:
        print("Brak wynikow - najpierw uruchom dock.py")
        sys.exit(1)
    write_runs(df)
    sdf = summarize(df)
    plot(df, sdf)
    print("\nGotowe (analyze).")


if __name__ == "__main__":
    main()
