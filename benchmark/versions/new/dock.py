#!/usr/bin/env python3
import os
import sys
import time
from vina import Vina

import config as C


def run_vina(maps_prefix, ligand, out_pdbqt, log_path, seed):
    # liczenie czasu
    t0 = time.time()

    #jest to dokowanie z forcefieldem z AD4
    #TODO trzeba zrobić później tak, żeby był eksperyment z innymi forcefieldami, np. vina, a4, vinardo 
    v = Vina(sf_name="ad4", cpu=C.CPU, seed=seed, verbosity=0)
    v.load_maps(maps_prefix)
    v.set_ligand_from_file(ligand)

    #obecnie jest dane 20, żeby było 20 pozycji zamiast 9
    v.dock(exhaustiveness=C.EXHAUSTIVENESS, n_poses=C.NUM_MODES)


    #zapisujemy wszystkie pozycje jakie są dokowane, dodatkowo jest ten parametr zwiększony do 10
    v.write_poses(out_pdbqt, n_poses=C.NUM_MODES, energy_range=C.ENERGY_RANGE, overwrite=True)
    energies = v.energies(n_poses=C.NUM_MODES)

    dt = time.time() - t0
    with open(log_path, "w") as f:
        f.write("maps %s\n"% maps_prefix)
        f.write("ligand %s\n"% ligand)
        f.write("seed %d\nexhaustiveness %d\nnum_modes %d\n"% (seed, C.EXHAUSTIVENESS, C.NUM_MODES))
        f.write("mode affinity (kcal/mol)\n")
        for i, e in enumerate(energies, 1):
            f.write("%4d %10.3f\n" % (i, e[0]))
    return dt
    


def dock_target(pdbid):
    tdir = os.path.join(C.RESULTS, pdbid)
    ligand = os.path.join(tdir, "ligand.pdbqt")

    #tworzymy mapy dla obu warunków, konwencjonalnego i z biasem
    maps = {
        "conventional": os.path.join(tdir, "receptor"),
        "biased":       os.path.join(tdir, "maps_biased", "receptor")}

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
    #pomin nazwe skryptu i parametry a jak nie to z configu
    arguments = sys.argv[1:]
    if arguments:
        which = arguments
    else:
        which = list(C.TARGETS)
    for pdbid in which:
        dock_target(pdbid)

if __name__ == "__main__":
    main()
