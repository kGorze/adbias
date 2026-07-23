#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ocena poz: symetryczny RMSD (OpenBabel obrms) docked vs natywna poza krystaliczna.
Metryka pierwszorzedna: RMSD, prog sukcesu <= config.RMSD_SUCCESS, dla top-1 i best-of-N.
Score Vina raportowany tylko diagnostycznie (nieporownywalny miedzy warunkami).

Wyniki:
    results/results.csv    (per przebieg: cel, warunek, seed, RMSD top-1/best, score)
    results/summary.csv    (per cel x warunek: success@2A, mediany RMSD, Wilcoxon)
    results/plot.png       (success@2A i rozklad RMSD)
"""
import os
import csv
import sys
import subprocess

import config as C


def vina_scores(pdbqt):
    """Lista affinity (kcal/mol) dla kolejnych poz z out.pdbqt."""
    scores = []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                scores.append(float(line.split(":")[1].split()[0]))
    return scores


def rmsd_poses(native_sdf, out_pdbqt, work):
    """Symetryczny heavy-atom RMSD kazdej pozy vs natywna (obabel obrms)."""
    poses_sdf = os.path.join(work, "poses.sdf")
    subprocess.check_call([C.OBABEL, out_pdbqt, "-O", poses_sdf],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = subprocess.check_output([C.OBRMS, "-f", native_sdf, poses_sdf],
                                  stderr=subprocess.DEVNULL).decode()
    vals = []
    for line in out.splitlines():
        # format obrms: "RMSD <ref>:<test> <wartosc>"  -> RMSD = ostatni token
        if line.startswith("RMSD"):
            try:
                vals.append(float(line.split()[-1]))
            except ValueError:
                pass
    return vals


# --------------------------------------------------------------------------
def collect():
    rows = []
    for pdbid in C.TARGETS:
        tdir = os.path.join(C.RESULTS, pdbid)
        native = os.path.join(tdir, "native.sdf")
        for cond in C.CONDITIONS:
            for seed in C.SEEDS:
                seed_dir = os.path.join(tdir, cond, "seed_%d" % seed)
                out_pdbqt = os.path.join(seed_dir, "out.pdbqt")
                if not os.path.isfile(out_pdbqt):
                    continue
                rmsds = rmsd_poses(native, out_pdbqt, seed_dir)
                scores = vina_scores(out_pdbqt)
                if not rmsds:
                    continue
                top1 = rmsds[0]
                best = min(rmsds)
                rows.append(dict(
                    pdbid=pdbid, condition=cond, seed=seed,
                    top1_rmsd=top1, best_rmsd=best,
                    best_pose=rmsds.index(best) + 1,
                    top1_score=scores[0] if scores else "",
                    top1_success=int(top1 <= C.RMSD_SUCCESS),
                    best_success=int(best <= C.RMSD_SUCCESS),
                ))
    return rows


def write_runs(rows):
    path = os.path.join(C.RESULTS, "results.csv")
    cols = ["pdbid", "condition", "seed", "top1_rmsd", "best_rmsd",
            "best_pose", "top1_score", "top1_success", "best_success"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 3) if isinstance(r[k], float) else r[k]) for k in cols})
    print("zapisano " + path)


def summarize(rows):
    try:
        from scipy.stats import wilcoxon
    except Exception:
        wilcoxon = None

    def mean(v):
        return sum(v) / len(v) if v else float("nan")

    def median(v):
        v = sorted(v)
        n = len(v)
        if not n:
            return float("nan")
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0

    path = os.path.join(C.RESULTS, "summary.csv")
    cols = ["pdbid", "condition", "n", "top1_success_rate", "best_success_rate",
            "top1_rmsd_mean", "top1_rmsd_median", "best_rmsd_median", "wilcoxon_p_top1"]
    summary = []
    for pdbid in C.TARGETS:
        # sparowany Wilcoxon: top-1 RMSD conventional vs biased (po seedach)
        paired = {c: {r["seed"]: r["top1_rmsd"] for r in rows
                      if r["pdbid"] == pdbid and r["condition"] == c} for c in C.CONDITIONS}
        common = sorted(set(paired["conventional"]) & set(paired["biased"]))
        p_val = ""
        if wilcoxon and len(common) >= 1:
            a = [paired["conventional"][s] for s in common]
            b = [paired["biased"][s] for s in common]
            if any(x != y for x, y in zip(a, b)):
                try:
                    p_val = round(wilcoxon(a, b).pvalue, 4)
                except Exception:
                    p_val = ""
        for cond in C.CONDITIONS:
            g = [r for r in rows if r["pdbid"] == pdbid and r["condition"] == cond]
            if not g:
                continue
            summary.append(dict(
                pdbid=pdbid, condition=cond, n=len(g),
                top1_success_rate=round(mean([r["top1_success"] for r in g]), 3),
                best_success_rate=round(mean([r["best_success"] for r in g]), 3),
                top1_rmsd_mean=round(mean([r["top1_rmsd"] for r in g]), 3),
                top1_rmsd_median=round(median([r["top1_rmsd"] for r in g]), 3),
                best_rmsd_median=round(median([r["best_rmsd"] for r in g]), 3),
                wilcoxon_p_top1=p_val if cond == "biased" else "",
            ))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for s in summary:
            w.writerow(s)
    print("zapisano " + path)
    return summary


def plot(rows, summary):
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
            s = [r for r in summary if r["pdbid"] == t and r["condition"] == cond]
            vals.append(s[0]["top1_success_rate"] if s else 0)
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
            ys = [r["top1_rmsd"] for r in rows if r["pdbid"] == t and r["condition"] == cond]
            xs = [j + (i - 0.5) * w] * len(ys)
            ax2.scatter(xs, ys, s=18, alpha=0.6, color=colors[cond],
                        label=cond if j == 0 else None)
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
    rows = collect()
    if not rows:
        print("Brak wynikow - najpierw uruchom dock.py")
        sys.exit(1)
    write_runs(rows)
    summary = summarize(rows)
    plot(rows, summary)
    print("\nGotowe (analyze).")


if __name__ == "__main__":
    main()
