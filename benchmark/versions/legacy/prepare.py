#!/usr/bin/env python3
"""
Przygotowanie kazdego celu wg klasycznego protokolu AutoDock4 (wersja "legacy").
Dla kazdego PDB ID z config.TARGETS tworzy w results/<ID>/:
    receptor.pdbqt        (cala jednostka asymetryczna z CIF, bez wod/HETATM poza jonami)
    ligand.pdbqt          (natywny ligand: obabel pH7.4 -> prepare_ligand4.py)
    native.sdf            (natywna poza krystaliczna = referencja do RMSD)
    receptor.gpf + mapy   (prepare_gpf4.py -> autogrid4)
    bias.bpf + bias_report.txt   (miejsce biasu z geometrii kontaktu H-mostka)
    maps_biased/          (zestaw map z wpisanym biasem, do warunku "biased")

z interpreterem config.PY3 (ma numpy do AutoDockTools_py3).
"""
import os
import sys
import shutil
import subprocess

import config as C


def sh(cmd, cwd=None, use_tools_path=False):
    env = os.environ.copy()

    if use_tools_path:
        env["PYTHONPATH"] = C.TOOLS

    print("  $ " + " ".join(cmd))

    #zamiast robić subprocess.run(..., check=True) bo wtedy nie widać stdout/stderr
    subprocess.check_call(cmd, cwd=cwd, env=env)


# --------------------------------------------------------------------------
# 1. CIF -> czysty PDB receptora (cala jednostka asymetryczna, bez wod,
#    bez HETATM poza jonami metali, jeden model, altloc A)
def cif_to_receptor_pdb(cif_path, out_pdb):
    headers = []
    rows = []
    in_loop = False
    with open(cif_path) as f:
        for line in f:
            s = line.strip()
            if s == "loop_":
                in_loop = True
                headers = []
                continue
            if in_loop and s.startswith("_atom_site."):
                headers.append(s.split(".", 1)[1])
                continue
            if headers and (s.startswith("ATOM") or s.startswith("HETATM")):
                rows.append(s.split())
            elif headers and (s == "#" or s.startswith("loop_") or s.startswith("_")):
                # koniec bloku _atom_site
                if rows:
                    break
                headers = []
                in_loop = s == "loop_"

    if not headers or not rows:
        raise RuntimeError("Nie znaleziono bloku _atom_site w %s" % cif_path)

    idx = {name: i for i, name in enumerate(headers)}

    def col(row, *names, default="."):
        for n in names:
            if n in idx:
                return row[idx[n]]
        return default

    first_model = None
    chain_map = {}
    next_chain = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    out = []
    serial = 0
    for row in rows:
        if len(row) < len(headers):
            continue
        model = col(row, "pdbx_PDB_model_num", default="1")
        if first_model is None:
            first_model = model
        if model != first_model:
            continue
        alt = col(row, "label_alt_id", default=".")
        if alt not in (".", "A"):
            continue

        group = row[idx["group_PDB"]]
        resname = col(row, "auth_comp_id", "label_comp_id").strip('"')
        if resname in C.WATER_RESNAMES:
            continue
        if group == "HETATM" and resname not in C.METAL_IONS:
            continue  # zostawiamy tylko jony metali

        element = row[idx["type_symbol"]].strip('"')
        atom = col(row, "auth_atom_id", "label_atom_id").strip('"')
        chain_in = col(row, "auth_asym_id", "label_asym_id", default="A")
        if chain_in not in chain_map:
            chain_map[chain_in] = chain_in if len(chain_in) == 1 else next(next_chain)
        chain = chain_map[chain_in]
        resseq = col(row, "auth_seq_id", "label_seq_id", default="1")
        icode = col(row, "pdbx_PDB_ins_code", default=".")
        icode = "" if icode in (".", "?") else icode
        try:
            resnum = int(resseq)
        except ValueError:
            resnum = 1
        x = float(row[idx["Cartn_x"]])
        y = float(row[idx["Cartn_y"]])
        z = float(row[idx["Cartn_z"]])

        serial += 1
        # kolumny nazwy atomu wg reguly PDB
        if len(atom) >= 4:
            name4 = atom[:4]
        elif len(element) == 2:
            name4 = atom.ljust(4)
        else:
            name4 = (" " + atom).ljust(4)
        out.append("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                   % (group, serial, name4, " ", resname[:3], chain, resnum,
                      (icode or " ")[:1], x, y, z, 1.0, 0.0, element.rjust(2)))
    out.append("END\n")
    with open(out_pdb, "w") as f:
        f.writelines(out)
    return len(out) - 1


# --------------------------------------------------------------------------
def fix_ligand_sdf(src, dst):
    """Normalizuje 2-literowe symbole pierwiastkow w SDF (np. 'BR'->'Br', 'CL'->'Cl'),
    bo obabel inaczej wpisuje atom '*' i prepare_ligand4.py sie wywala."""
    with open(src) as f:
        lines = f.readlines()
    natoms = int(lines[3][:3])                 # linia counts V2000
    for i in range(4, 4 + natoms):
        sym = lines[i][31:34].strip()          # symbol pierwiastka: kolumny 32-34
        if len(sym) == 2:
            lines[i] = lines[i][:31] + (sym[0].upper() + sym[1].lower()).ljust(3) + lines[i][34:]
    with open(dst, "w") as f:
        f.writelines(lines)


def read_pdbqt_atoms(path):
    """Zwraca liste atomow (name, resname, chain, resnum, x, y, z, adtype)."""
    atoms = []
    with open(path) as f:
        for line in f:
            if line[:6] not in ("ATOM  ", "HETATM"):
                continue
            atoms.append(dict(
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain=line[21:22].strip(),
                resnum=line[22:26].strip(),
                x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
                adtype=line.rstrip().split()[-1],
            ))
    return atoms


def dist2(a, b):
    return (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2


# --------------------------------------------------------------------------
# 2. Pudelko: srodek + npts z natywnego liganda + padding
def grid_box(ligand_atoms):
    xs = [a["x"] for a in ligand_atoms]
    ys = [a["y"] for a in ligand_atoms]
    zs = [a["z"] for a in ligand_atoms]
    center = [(min(v) + max(v)) / 2.0 for v in (xs, ys, zs)]
    npts = []
    for v in (xs, ys, zs):
        size = (max(v) - min(v)) + 2 * C.PADDING
        n = int(round(size / C.SPACING))
        n += n % 2                       # autogrid wymaga parzystej liczby
        npts.append(min(n, C.NPTS_MAX))
    return center, npts


# --------------------------------------------------------------------------
# 3. Miejsce biasu: najkrotszy kontakt heteroatom(ligand)-heteroatom(receptor) <= cutoff.
#    Typ 'acc' jesli atom liganda jest akceptorem (OA/NA/SA), inaczej 'don'.
def extract_bias_site(ligand_pdbqt, receptor_pdbqt):
    lig = read_pdbqt_atoms(ligand_pdbqt)
    rec = read_pdbqt_atoms(receptor_pdbqt)
    acc_types = ("OA", "NA", "SA")

    # tylko ciezkie heteroatomy (N/O/S); atomy H (w tym HD) wykluczone
    lig_polar = [a for a in lig if a["adtype"][0] in ("N", "O", "S")]
    rec_polar = [a for a in rec if a["adtype"][0] in ("N", "O")]

    cutoff2 = C.HBOND_CUTOFF ** 2
    best = None
    for la in lig_polar:
        for ra in rec_polar:
            d2 = dist2(la, ra)
            if best is None or d2 < best[0]:
                best = (d2, la, ra)
    if best is None:
        return None, "brak atomow polarnych"

    d2, la, ra = best
    d = d2 ** 0.5
    btype = "acc" if la["adtype"] in acc_types else "don"
    report = ("ligand atom %s (%s)  <->  receptor %s%s/%s (%s)   d=%.2f A   -> bias '%s'%s"
              % (la["name"], la["adtype"], ra["resname"], ra["resnum"], ra["name"],
                 ra["adtype"], d, btype,
                 "" if d <= C.HBOND_CUTOFF else "   [UWAGA: > cutoff %.1f A]" % C.HBOND_CUTOFF))
    site = dict(x=la["x"], y=la["y"], z=la["z"], type=btype)
    return site, report


# --------------------------------------------------------------------------
def prepare_target(pdbid, spec):
    print("\n=== %s ===" % pdbid)
    tdir = os.path.join(C.RESULTS, pdbid)
    if os.path.isdir(tdir):
        shutil.rmtree(tdir)
    os.makedirs(tdir)

    cif = os.path.join(C.DATA, spec["cif"])
    sdf = os.path.join(C.DATA, spec["ligand_sdf"])

    receptor_pdb = os.path.join(tdir, "receptor_clean.pdb")
    receptor_pdbqt = os.path.join(tdir, "receptor.pdbqt")
    ligand_pdb = os.path.join(tdir, "ligand.pdb")
    ligand_pdbqt = os.path.join(tdir, "ligand.pdbqt")
    native_sdf = os.path.join(tdir, "native.sdf")


    #zamiana cif na pdb
    n = cif_to_receptor_pdb(cif, receptor_pdb)

    # odpalenie 
    sh([C.PY3, C.PREPARE_RECEPTOR, "-r", receptor_pdb, "-o", receptor_pdbqt,
        "-A", "checkhydrogens", "-U", "nphs_lps_waters"], cwd=tdir, use_tools_path=True)

    # --- ligand (obabel pH7.4 -> prepare_ligand4.py); native.sdf z pdbqt = ten sam graf co pozy ---
    sdf_fixed = os.path.join(tdir, "ligand_input.sdf")
    fix_ligand_sdf(sdf, sdf_fixed)
    sh([C.OBABEL, sdf_fixed, "-O", ligand_pdb, "-p", "7.4"], cwd=tdir)
    sh([C.PY3, C.PREPARE_LIGAND, "-l", ligand_pdb, "-o", ligand_pdbqt, "-U", "nphs_lps"],
       cwd=tdir, use_tools_path=True)
    sh([C.OBABEL, ligand_pdbqt, "-O", native_sdf], cwd=tdir)

    # --- pudelko + mapy AD4 ---
    lig_atoms = read_pdbqt_atoms(ligand_pdbqt)
    center, npts = grid_box(lig_atoms)
    with open(os.path.join(tdir, "box.txt"), "w") as f:
        f.write("center = %.3f %.3f %.3f\nnpts   = %d %d %d\nspacing= %.3f\n"
                % (center[0], center[1], center[2], npts[0], npts[1], npts[2], C.SPACING))
    print("  box center=%.2f,%.2f,%.2f npts=%s" % (center[0], center[1], center[2], npts))
    sh([C.PY3, C.PREPARE_GPF, "-l", ligand_pdbqt, "-r", receptor_pdbqt, "-o", "receptor.gpf",
        "-p", "gridcenter=%.3f,%.3f,%.3f" % tuple(center),
        "-p", "npts=%d,%d,%d" % tuple(npts),
        "-p", "spacing=%.3f" % C.SPACING], cwd=tdir, use_tools_path=True)
    sh([C.AUTOGRID4, "-p", "receptor.gpf", "-l", "receptor.glg"], cwd=tdir)

    # --- bias: miejsce + mapy biasowane ---
    site, report = extract_bias_site(ligand_pdbqt, receptor_pdbqt)
    with open(os.path.join(tdir, "bias_report.txt"), "w") as f:
        f.write(report + "\n")
    print("  bias: " + report)
    if site is None:
        raise RuntimeError("%s: nie wyznaczono miejsca biasu" % pdbid)

    bpf = os.path.join(tdir, "bias.bpf")
    with open(bpf, "w") as f:
        f.write("x y z Vset r type\n")
        f.write("%.3f %.3f %.3f %.2f %.2f %s\n"
                % (site["x"], site["y"], site["z"], C.BIAS_VSET, C.BIAS_RADIUS, site["type"]))
    # prepare_bias.py (python2) w katalogu z gpf+mapami -> tworzy receptor.*.biased.map
    sh([C.PY2, C.PREPARE_BIAS, "-b", "bias.bpf", "-g", "receptor.gpf"], cwd=tdir)

    # --- zestaw map biased: kopia map + podmiana tych z biasem ---
    mbias = os.path.join(tdir, "maps_biased")
    os.makedirs(mbias)
    for fn in os.listdir(tdir):
        if fn.endswith(".map") and ".biased." not in fn:
            shutil.copy(os.path.join(tdir, fn), os.path.join(mbias, fn))
        elif fn in ("receptor.maps.fld", "receptor.maps.xyz"):
            shutil.copy(os.path.join(tdir, fn), os.path.join(mbias, fn))
    swapped = []
    for fn in os.listdir(tdir):
        if fn.endswith(".biased.map"):
            target = fn.replace(".biased.map", ".map")
            shutil.copy(os.path.join(tdir, fn), os.path.join(mbias, target))
            swapped.append(target)
    print("  mapy biased podmienione: %s" % ", ".join(sorted(swapped)))


# --------------------------------------------------------------------------
def main():
    which = sys.argv[1:] or list(C.TARGETS)
    for pdbid in which:
        prepare_target(pdbid, C.TARGETS[pdbid])
    print("\nGotowe (prepare).")


if __name__ == "__main__":
    main()
