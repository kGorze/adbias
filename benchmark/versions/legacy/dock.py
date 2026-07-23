#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Docking Vina (scoring ad4, mapy AD4) dla kazdego celu.
Petla: cel x {conventional, biased} x seed.
Warunki roznia sie WYLACZNIE zestawem map (conventional vs maps_biased/);
te same seedy w obu warunkach (test sparowany).

Uruchamiac po prepare.py.  Argumenty opcjonalne: lista PDB ID.
"""
import os
import sys
import subprocess

import config as C


def run_vina(maps_dir, ligand, out_pdbqt, log_path, seed):
    cmd = [C.VINA,
           "--scoring", "ad4",
           "--maps", "receptor",
           "--ligand", ligand,
           "--out", out_pdbqt,
           "--seed", str(seed),
           "--exhaustiveness", str(C.EXHAUSTIVENESS),
           "--num_modes", str(C.NUM_MODES),
           "--energy_range", str(C.ENERGY_RANGE),
           "--cpu", str(C.CPU)]
    with open(log_path, "w") as log:
        subprocess.check_call(
            cmd, 
            cwd=maps_dir, 
            stdout=log, 
            stderr=subprocess.STDOUT
            )


def dock_target(pdbid):
    print("\n=== %s ===" % pdbid)
    tdir = os.path.join(C.RESULTS, pdbid)
    ligand = os.path.join(tdir, "ligand.pdbqt")
    maps = {"conventional": tdir, "biased": os.path.join(tdir, "maps_biased")}

    for cond in C.CONDITIONS:
        for seed in C.SEEDS:
            seed_dir = os.path.join(tdir, cond, "seed_%d" % seed)
            os.makedirs(seed_dir, exist_ok=True)
            out_pdbqt = os.path.join(seed_dir, "out.pdbqt")
            log_path = os.path.join(seed_dir, "vina.log")
            run_vina(maps[cond], ligand, out_pdbqt, log_path, seed)
        print("  %-12s: %d seedow" % (cond, len(C.SEEDS)))


def main():
    which = sys.argv[1:] or list(C.TARGETS)
    for pdbid in which:
        dock_target(pdbid)
    print("\nGotowe (dock).")


if __name__ == "__main__":
    main()
