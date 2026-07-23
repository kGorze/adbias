#!/usr/bin/env python3
"""
autodocktools a meeko na tych samych celach, warunkach i seedach.
obrms a rdkit

do użycia z interpreterem z vs (rdkit + meeko):
/home/kgorzelanczyk/miniforge3/envs/vs/bin/python compare.py
"""

import os
import csv
import importlib.util
import subprocess
import tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
COMPARISON_RESULTS = os.path.join(HERE, "comparison_results")


def load_config(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LEG = load_config("cfg_legacy", os.path.join(HERE, "legacy", "config.py"))
NEW = load_config("cfg_new", os.path.join(HERE, "new", "config.py"))


def vina_scores(pdbqt):
    out = []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                out.append(float(line.split(":")[1].split()[0]))
    return out


def rmsd_obrms(native_sdf, out_pdbqt):
    """Symetryczny RMSD OpenBabela (metryka wersji legacy)."""
    with tempfile.TemporaryDirectory() as tmp:
        poses = os.path.join(tmp, "poses.sdf")
        subprocess.check_call([LEG.OBABEL, out_pdbqt, "-O", poses],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        txt = subprocess.check_output([LEG.OBRMS, "-f", native_sdf, poses],
                                      stderr=subprocess.DEVNULL).decode()
    vals = []
    for line in txt.splitlines():
        if line.startswith("RMSD"):
            try:
                vals.append(float(line.split()[-1]))
            except ValueError:
                pass
    return vals


def rmsd_rdkit(native_sdf, out_pdbqt):
    """Symetryczny RMSD RDKit (metryka wersji new)."""
    import sys
    sys.path.insert(0, os.path.join(HERE, "new"))
    from analyze import rmsd_poses
    return rmsd_poses(native_sdf, out_pdbqt)


def collect(cfg, version, rmsd_fn, extra_fn=None):
    rows = []
    for pdbid in cfg.TARGETS:
        tdir = os.path.join(cfg.RESULTS, pdbid)
        native = os.path.join(tdir, "native.sdf")
        for cond in cfg.CONDITIONS:
            for seed in cfg.SEEDS:
                out_pdbqt = os.path.join(tdir, cond, "seed_%d" % seed, "out.pdbqt")
                if not os.path.isfile(out_pdbqt):
                    continue
                r = rmsd_fn(native, out_pdbqt)
                s = vina_scores(out_pdbqt)
                if not r:
                    continue
                row = dict(version=version, pdbid=pdbid, condition=cond, seed=seed,
                           top1_rmsd=round(r[0], 3), best_rmsd=round(min(r), 3),
                           top1_score=s[0] if s else "")
                if extra_fn is not None:
                    r2 = extra_fn(native, out_pdbqt)
                    row["top1_rmsd_alt"] = round(r2[0], 3) if r2 else ""
                rows.append(row)
    return rows


def agg(rows, version, pdbid, cond, thr):
    g = [r for r in rows if r["version"] == version and r["pdbid"] == pdbid
         and r["condition"] == cond]
    if not g:
        return None
    top1 = np.array([r["top1_rmsd"] for r in g])
    best = np.array([r["best_rmsd"] for r in g])
    sc = np.array([r["top1_score"] for r in g], dtype=float)
    return dict(n=len(g),
                top1_success=round(float((top1 <= thr).mean()), 3),
                best_success=round(float((best <= thr).mean()), 3),
                top1_median=round(float(np.median(top1)), 3),
                best_median=round(float(np.median(best)), 3),
                score_mean=round(float(sc.mean()), 2))


def main():
    os.makedirs(COMPARISON_RESULTS, exist_ok=True)
    rows = []
    rows += collect(LEG, "legacy", rmsd_obrms)
    rows += collect(NEW, "new", rmsd_rdkit, extra_fn=rmsd_obrms)

    path = os.path.join(COMPARISON_RESULTS, "comparison_runs.csv")
    cols = ["version", "pdbid", "condition", "seed", "top1_rmsd", "best_rmsd",
            "top1_score", "top1_rmsd_alt"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print("zapisano " + path)

    # tabela zbiorcza tylko dla seedow policzonych w OBU wersjach
    common = {}
    for r in rows:
        common.setdefault((r["pdbid"], r["condition"], r["seed"]), set()).add(r["version"])
    keep = {k for k, v in common.items() if v == {"legacy", "new"}}
    rows_c = [r for r in rows if (r["pdbid"], r["condition"], r["seed"]) in keep]

    thr = NEW.RMSD_SUCCESS
    out = []
    print("\nwspolne seedy: %d przebiegow na wersje" % (len(rows_c) // 2))
    hdr = ("%-6s %-13s %-7s %3s  %7s %7s  %8s %8s  %7s"
           % ("cel", "warunek", "wersja", "n", "top1med", "bestmed",
              "succ@%.0fA" % thr, "succ_best", "score"))
    print(hdr)
    print("-" * len(hdr))
    for pdbid in NEW.TARGETS:
        for cond in NEW.CONDITIONS:
            for version in ("legacy", "new"):
                a = agg(rows_c, version, pdbid, cond, thr)
                if a is None:
                    continue
                print("%-6s %-13s %-7s %3d  %7.2f %7.2f  %8.2f %8.2f  %7.2f"
                      % (pdbid, cond, version, a["n"], a["top1_median"], a["best_median"],
                         a["top1_success"], a["best_success"], a["score_mean"]))
                out.append(dict(pdbid=pdbid, condition=cond, version=version, **a))
        print()

    path = os.path.join(COMPARISON_RESULTS, "comparison_summary.csv")
    with open(path, "w", newline="") as f:
        cols = ["pdbid", "condition", "version", "n", "top1_median", "best_median",
                "top1_success", "best_success", "score_mean"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow(r)
    print("zapisano " + path)

    # kontrola metryki: RDKit vs obrms na tych samych pozach wersji "new"
    pairs = [(r["top1_rmsd"], r["top1_rmsd_alt"]) for r in rows
             if r["version"] == "new" and r.get("top1_rmsd_alt") not in ("", None)]
    if pairs:
        d = np.array([abs(a - b) for a, b in pairs])
        print("\nkontrola metryki (te same pozy 'new'): |RDKit - obrms| "
              "mediana %.3f A, max %.3f A, n=%d" % (np.median(d), d.max(), len(d)))


if __name__ == "__main__":
    main()
