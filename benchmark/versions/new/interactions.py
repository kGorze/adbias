#!/usr/bin/env python3
"""
port ideal_interaction_sites.py (adbias/code, Python2 + Bio.PDB) do "new"

uzycie jako: 
- python interactions.py -i 3CS9 -c A -r 239,241,244
- python interactions.py -i /sciezka/do/receptor.pdb -c A -r 22,54,61
"""
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path

import config as C

SITE_NAMES = {"DON": "donor", "ACC": "acceptor", "ARO": "aromatic"}


def resolve_receptor(receptor):
    if receptor in C.TARGETS:
        return Path(C.RESULTS) / receptor / "receptor_prepared.pdb", receptor
    path = Path(receptor)
    if not path.is_file():
        raise SystemExit(
            "receptor '%s' nie jest ani celem z config.TARGETS (%s), ani "
            "istniejacym plikiem" % (receptor, ", ".join(C.TARGETS))
        )
    return path, path.stem


def run_ideal_interaction_sites(receptor_pdb, chain, residues, work_dir):
    command = [C.PY2, C.IDEAL_INTERACTION_SITES, "-i", str(receptor_pdb), "-c", chain, "-r", residues]
    print("  command:", " ".join(command))
    subprocess.run(command, cwd=work_dir, check=True)
    return Path(work_dir) / "interaction_sites.pdb"


def parse_interaction_sites(pdb_path):
    sites = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            sites.append((SITE_NAMES.get(resname, resname), x, y, z))
    return sites


def write_report(sites, out_path, receptor_pdb, chain, residues):
    with open(out_path, "w") as f:
        f.write("# receptor  %s\n" % receptor_pdb)
        f.write("# chain     %s\n" % chain)
        f.write("# residues  %s\n" % residues)
        f.write("# n_sites   %d\n" % len(sites))
        f.write("#\n")
        f.write("#%-9s %-9s %9s %9s %9s\n" % ("index", "type", "x", "y", "z"))
        for index, (site_type, x, y, z) in enumerate(sites, 1):
            f.write("%-10d %-9s %9.3f %9.3f %9.3f\n" % (index, site_type, x, y, z))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--receptor", required=True,
                         help="ID celu z config.TARGETS (np. 3CS9) albo sciezka do PDB z wodorami")
    parser.add_argument("-c", "--chain", required=True, help="ID lancucha (np. A)")
    parser.add_argument("-r", "--residues", required=True,
                         help="numery reszt oddzielone przecinkami (np. 22,54,61)")
    parser.add_argument("-o", "--out",
                         help="plik wyjsciowy .txt; plik .pdb powstaje obok niego z ta sama nazwa bazowa")
    args = parser.parse_args()

    receptor_pdb, label = resolve_receptor(args.receptor)

    if args.out:
        out_path = Path(args.out)
    else:
        residues_tag = args.residues.replace(",", "-")
        out_path = Path(C.RESULTS) / label / "interactions" / ("%s_%s.txt" % (args.chain, residues_tag))
    pdb_out_path = out_path.with_suffix(".pdb")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as work_dir:
        sites_pdb = run_ideal_interaction_sites(receptor_pdb, args.chain, args.residues, work_dir)
        sites = parse_interaction_sites(sites_pdb)
        shutil.copyfile(sites_pdb, pdb_out_path)

    write_report(sites, out_path, receptor_pdb, args.chain, args.residues)
    print("  %d miejsc oddzialywan -> %s" % (len(sites), out_path))
    print("  miejsca oddzialywan PDB -> %s" % pdb_out_path)

#TODO: receptor_prepared.pdb ma histydyny jako HIS (nie HIE/HID/HIP), więc boczny łańcuch His nie generuje donor/acceptor sites — tylko ring aromatyczny działa. Jeśli w miejscu aktywnym chcemy katalityczna His, trzeba by ręcznie przemianować resname w kopii PDB przed uruchomieniem albo zająć się tym w kodzie.

if __name__ == "__main__":
    main()

