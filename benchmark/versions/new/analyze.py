#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdqueries
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate
from scipy.stats import wilcoxon

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
                # Funkcje pomocnicze mają prefiks "_", bo są używane tylko w tym pliku.
                rmsd_values = _rmsd_poses(
                    str(native_sdf_path), 
                    str(output_pdbqt_path)
                )
                if not rmsd_values:
                    continue
                # robimy to tutaj a nie na końcu w tworzeniu słownika
                score_values = _vina_scores(str(output_pdbqt_path))
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

# Tworzymy jeden wiersz podsumowania dla każdego celu i warunku.
def _summarize(df):
    summary_rows = []

    for pdbid in C.TARGETS:

        #wywołanie testu wilcoxona dla pary warunków (conventional vs biased) dla danego celu
        p_value = _calculate_paired_wilcoxon(df, pdbid)
        if np.isnan(p_value):
            print(f"blad testu Wilcoxona dla {pdbid}")

        for condition in C.CONDITIONS:
            group = df[
                (df["pdbid"] == pdbid)
                & (df["condition"] == condition)
            ]

            if group.empty:
                continue

            row = {
                "pdbid": pdbid,
                "condition": condition,
                "n": len(group),
                "top1_success_rate": round(
                    group["top1_success"].mean(),
                    3,
                ),
                "best_success_rate": round(
                    group["best_success"].mean(),
                    3,
                ),
                "top1_rmsd_mean": round(
                    group["top1_rmsd"].mean(),
                    3,
                ),
                "top1_rmsd_median": round(
                    group["top1_rmsd"].median(),
                    3,
                ),
                "best_rmsd_median": round(
                    group["best_rmsd"].median(),
                    3,
                ),
                "wilcoxon_p_top1": (
                    round(p_value, 4)
                    if condition == "biased"
                    and not np.isnan(p_value)
                    else np.nan
                ),
            }

            summary_rows.append(row)
    summary_dataframe = pd.DataFrame(summary_rows)
    output_path = os.path.join(
        C.RESULTS,
        "summary.csv",
    )
    summary_dataframe.to_csv(
        output_path,
        index=False,
    )
    return summary_dataframe


# Test porównuje tylko wyniki mające ten sam seed w obu warunkach.
def _calculate_paired_wilcoxon(
    df,
    pdbid,
    first_condition="conventional",
    second_condition="biased",
):
    target_rows = df[df["pdbid"] == pdbid]

    first_values = (
        target_rows[
            target_rows["condition"] == first_condition
        ]
        .set_index("seed")["top1_rmsd"]
    )

    second_values = (
        target_rows[
            target_rows["condition"] == second_condition
        ]
        .set_index("seed")["top1_rmsd"]
    )

    common_seeds = first_values.index.intersection(
        second_values.index
    )

    if len(common_seeds) == 0:
        return np.nan

    first_paired = first_values.loc[common_seeds].to_numpy()
    second_paired = second_values.loc[common_seeds].to_numpy()

    differences = first_paired - second_paired

    if np.all(differences == 0):
        return np.nan

    result = wilcoxon(first_paired, second_paired)
    return float(result.pvalue)

# Pobieramy wartości sukcesu potrzebne do narysowania jednej serii słupków.
def _get_success_rates(summary_df, targets, condition):
    success_rates = []

    # Kolejność wartości musi odpowiadać kolejności słupków na osi X.
    for target in targets:
        matching_rows = summary_df[
            (summary_df["pdbid"] == target)
            & (summary_df["condition"] == condition)
        ]
        if matching_rows.empty:
            success_rate = 0.0
        else:
            success_rate = float(matching_rows["top1_success_rate"].iloc[0])
        success_rates.append(success_rate)
    return success_rates


# Wyliczamy pozycję warunku wewnątrz grupy danego celu.
def get_condition_position(
    target_index,
    condition_index,
    width,
):
    # Warunki są przesunięte na boki względem środka grupy celu.
    return target_index + (condition_index - 0.5) * width


# Lewy panel pokazuje odsetek udanych wyników dla każdego celu.
def plot_success_rates(
    axis,
    summary_df,
    targets,
    colors,
    bar_width,
):
    for condition_index, condition in enumerate(C.CONDITIONS):
        success_rates = _get_success_rates(
            summary_df,
            targets,
            condition,
        )

        positions = []

        for target_index in range(len(targets)):
            position = get_condition_position(
                target_index,
                condition_index,
                bar_width,
            )
            positions.append(position)

        axis.bar(
            positions,
            success_rates,
            width=bar_width,
            label=condition,
            color=colors[condition],
        )

    axis.set_xticks(range(len(targets)))
    axis.set_xticklabels(targets)
    axis.set_ylabel(
        f"success@{C.RMSD_SUCCESS:.1f} A (top-1)"
    )
    axis.set_ylim(0, 1.05)
    axis.set_title("pozy (top-1)")
    axis.legend(frameon=False)

# Prawy panel pokazuje wszystkie wartości RMSD, bez ich agregowania.
def plot_rmsd_distribution(
    axis,
    results_df,
    targets,
    colors,
    point_offset,
):
    for condition_index, condition in enumerate(C.CONDITIONS):
        for target_index, target in enumerate(targets):
            matching_rows = results_df[
                (results_df["pdbid"] == target)
                & (results_df["condition"] == condition)
            ]

            rmsd_values = matching_rows["top1_rmsd"].to_numpy()

            x_position = get_condition_position(
                target_index,
                condition_index,
                point_offset,
            )

            axis.scatter(
                [x_position] * len(rmsd_values),
                rmsd_values,
                s=18,
                alpha=0.6,
                color=colors[condition],
                label=condition if target_index == 0 else None,
            )

    axis.axhline(
        C.RMSD_SUCCESS,
        linestyle="--",
        linewidth=1,
        color="gray",
    )

    axis.set_xticks(range(len(targets)))
    axis.set_xticklabels(targets)
    axis.set_ylabel("top-1 RMSD (A)")
    axis.set_title("dystrybucja top-1 RMSD")
    axis.legend(frameon=False)

# Składamy oba panele i zapisujemy gotowy wykres do pliku.
def plot(df, sdf):
    targets = list(C.TARGETS)
    colors = {
        "conventional": "#698BD1",
        "biased": "#A1856AFF",
    }
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11, 4.5),
    )
    point_offset = 0.38
    plot_success_rates(
        axis=axes[0],
        summary_df=sdf,
        targets=targets,
        colors=colors,
        bar_width=point_offset,
    )
    plot_rmsd_distribution(
        axis=axes[1],
        results_df=df,
        targets=targets,
        colors=colors,
        point_offset=point_offset,
    )
    figure.tight_layout()
    output_path = os.path.join(
        C.RESULTS,
        "plot.png",
    )
    figure.savefig(output_path, dpi=130)
    plt.close(figure)

def main():
    # Główne kroki analizy są prywatnymi funkcjami tego modułu.
    df = _collect()
    _write_runs(df)
    sdf = _summarize(df)
    plot(df, sdf)

if __name__ == "__main__":
    main()
