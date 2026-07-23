#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit import rdqueries
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

import config as C

# żeby wszystko wyłączyć potzrebne jest rdApp, inaczej nie blokujemy komunikatów z C++
RDLogger.DisableLog("rdApp.warning")


def _topology_query(mol: Chem.Mol) -> Chem.Mol:
    """
    budujemy z mol wzorzecz do GetSubstructMatches, ktory dopasuje atomy tylko po liczbie atomowej
    """
    query = Chem.RWMol(mol)
    for atom in query.GetAtoms():
        query.ReplaceAtomWithQueryAtom(
            atom.GetIdx(),
            rdqueries.AtomNumEqualsQueryAtom(atom.GetAtomicNum())
        )
    return query

#obecnie nieużywane bo jest źle zaprojektowane
def _neutral_copy(mol):
    """
    kopia bez ladunkow formalnych,
    """
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
    """
    Symetryczny RMSD po wszystkich dopasowaniach grafu
    """
    matches = ref.GetSubstructMatches(
        _topology_query(probe),
        uniquify=False,
        maxMatches=100000
    )

    if not matches:
        raise RuntimeError("poza nie pasuje grafem do natywnego liganda")

    p = probe.GetConformer(conf_id).GetPositions()
    r = ref.GetConformer().GetPositions()

    # usunięcie cold golfu, żeby to było poprawnie zrobione
    rmsd_values = []
    for match in matches:
        matched_reference = r[list(match)]
        squared_distances = ((p - matched_reference) ** 2).sum(axis=1)
        rmsd = np.sqrt(squared_distances.mean())
        rmsd_values.append(float(rmsd))

    return min(rmsd_values)


def _rmsd_poses(native_sdf, out_pdbqt):
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


def _vina_scores(pdbqt):
    """
    lista affinity (kcal/mol) dla kolejnych poz z out.pdbqt.
    REMARK VINA RESULT:    -8.4      0.000      0.000
    """

    scores = []

    with open(pdbqt) as f:
        for line in f:
            if not line.startswith("REMARK VINA RESULT:"):
                continue

            _, results_text = line.split("REMARK VINA RESULT:")
            affinity_text, rmsd_lower, rmsd_upper = results_text.split()[:3]

            affinity = float(affinity_text)
            scores.append(affinity)

    return scores


def _collect():
    rows = []
    #przejście na pathlib zamiast os jest wygodniejsze i bardziej czytelne
    results_directory = Path(C.RESULTS)

    for pdbid in C.TARGETS:
        target_directory = results_directory / pdbid
        native_sdf_path = target_directory / "native.sdf"

        for condition in C.CONDITIONS:
            for seed in C.SEEDS:
                output_pdbqt_path = (
                    target_directory / condition / f"seed_{seed}" / "out.pdbqt"
                )
                if not output_pdbqt_path.is_file():
                    continue
                rmsd_values = rmsd_poses(
                    str(native_sdf_path), 
                    str(output_pdbqt_path)
                )
                if not rmsd_values:
                    continue
                # robimy to tutaj a nie na końcu w tworzeniu słownika
                score_values = vina_scores(str(output_pdbqt_path))
                top1_rmsd = rmsd_values[0]
                best_rmsd = min(rmsd_values)
                best_pose_number = rmsd_values.index(best_rmsd) + 1

                if score_values:
                    top1_score = score_values[0]
                else:
                    top1_score = np.nan
                row = {
                    "pdbid": pdbid,
                    "condition": condition,
                    "seed": seed,
                    "top1_rmsd": round(top1_rmsd, 3),
                    "best_rmsd": round(best_rmsd, 3),
                    "best_pose": best_pose_number,
                    "top1_score": top1_score,
                    "top1_success": top1_rmsd <= C.RMSD_SUCCESS,
                    "best_success": best_rmsd <= C.RMSD_SUCCESS,
                }
                rows.append(row)
    return pd.DataFrame(rows)


#helper do zapisywania wyników
def _write_runs(df):
    results_directory = Path(C.RESULTS)
    path = results_directory / "results.csv"
    df.to_csv(path, index=False)
    print("zapisano " + str(path))

# 
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
