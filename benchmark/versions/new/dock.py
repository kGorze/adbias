#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Docking Vina (scoring ad4, mapy AD4) dla kazdego celu - przez python API Vina
(pakiet `vina`), zamiast wolania binarki jak w wersji legacy.

Petla: cel x {conventional, biased} x seed.
Warunki roznia sie WYLACZNIE zestawem map (conventional vs maps_biased/);
te same seedy w obu warunkach (test sparowany).

Uruchamiac po prepare.py.  Argumenty opcjonalne: lista PDB ID.
"""
import os
import sys
import time

from vina import Vina

import config as C


def run_vina(maps_prefix, ligand, out_pdbqt, log_path, seed):
    t0 = time.time()
    v = Vina(sf_name="ad4", cpu=C.CPU, seed=seed, verbosity=0)
    v.load_maps(maps_prefix)
    v.set_ligand_from_file(ligand)
    v.dock(exhaustiveness=C.EXHAUSTIVENESS, n_poses=C.NUM_MODES)
    v.write_poses(out_pdbqt, n_poses=C.NUM_MODES, energy_range=C.ENERGY_RANGE, overwrite=True)
    energies = v.energies(n_poses=C.NUM_MODES)
    dt = time.time() - t0

    with open(log_path, "w") as f:
        f.write("maps        %s\n" % maps_prefix)
        f.write("ligand      %s\n" % ligand)
        f.write("seed        %d\nexhaustiveness %d\nnum_modes   %d\n"
                "energy_range %d\ncpu         %d\n"
                % (seed, C.EXHAUSTIVENESS, C.NUM_MODES, C.ENERGY_RANGE, C.CPU))
        f.write("time_s      %.1f\n\n" % dt)
        f.write("mode |   affinity (kcal/mol)\n-----+---------------------\n")
        for i, e in enumerate(energies, 1):
            f.write("%4d | %10.3f\n" % (i, e[0]))
    return dt


def dock_target(pdbid):
    print("\n=== %s ===" % pdbid)
    tdir = os.path.join(C.RESULTS, pdbid)
    ligand = os.path.join(tdir, "ligand.pdbqt")
    maps = {"conventional": os.path.join(tdir, "receptor"),
            "biased": os.path.join(tdir, "maps_biased", "receptor")}

    for cond in C.CONDITIONS:
        total = 0.0
        for seed in C.SEEDS:
            seed_dir = os.path.join(tdir, cond, "seed_%d" % seed)
            os.makedirs(seed_dir, exist_ok=True)
            total += run_vina(maps[cond], ligand,
                              os.path.join(seed_dir, "out.pdbqt"),
                              os.path.join(seed_dir, "vina.log"), seed)
        print("  %-12s: %d seedow, %.0f s" % (cond, len(C.SEEDS), total))


def main():
    which = sys.argv[1:] or list(C.TARGETS)
    for pdbid in which:
        dock_target(pdbid)
    print("\nGotowe (dock).")


if __name__ == "__main__":
    main()
