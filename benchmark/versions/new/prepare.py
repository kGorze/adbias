#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Przygotowanie kazdego celu wg nowoczesnego protokolu (wersja "new").
Ten sam eksperyment co ../legacy, ale innym toolchainem:

    legacy                              new (tutaj)
    -----------------------------------------------------------------
    wlasny parser CIF -> PDB            ProDy (parseMMCIF + select)
    prepare_receptor4.py (ADT)          mk_prepare_receptor.py (Meeko)
    obabel -p 7.4 + prepare_ligand4.py  molscrub (pH 7.4) + RDKit + Meeko
    prepare_gpf4.py                     mk_prepare_receptor --write_gpf
    autogrid4                           autogrid4              (wspolne)
    prepare_bias.py (python2)           prepare_bias.py        (wspolne)

Dla kazdego PDB ID z config.TARGETS tworzy w results/<ID>/:
    receptor_clean.pdb    (jednostka asymetryczna z CIF, bez wod/HETATM poza jonami)
    receptor.pdbqt        (Meeko) + receptor_prepared.pdb + receptor_report.txt
    ligand.pdbqt          (Meeko; protonacja pH 7.4 przez molscrub)
    native.sdf            (natywna poza krystaliczna z H = referencja do RMSD)
    receptor.gpf + mapy   (Meeko --write_gpf -> autogrid4)
    bias.bpf + bias_report.txt   (miejsce biasu z geometrii kontaktu H-mostka)
    maps_biased/          (zestaw map z wpisanym biasem, do warunku "biased")

Uruchamiac interpreterem config.PY_NEW.
"""
import os
import re
import sys
import shutil
import subprocess

import numpy as np
import prody
from rdkit import Chem, RDLogger
from rdkit.Geometry import Point3D
from molscrub import Scrub
from meeko import MoleculePreparation, PDBQTWriterLegacy

import config as C

prody.confProDy(verbosity="critical")
RDLogger.DisableLog("rdApp.warning")


# --------------------------------------------------------------------------
def sh(cmd, cwd=None, log=None):
    """Uruchom polecenie, przerwij przy bledzie. Zwraca stdout+stderr."""
    print("  $ " + " ".join(str(c) for c in cmd))
    out = subprocess.run([str(c) for c in cmd], cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, check=True).stdout.decode()
    if log:
        with open(log, "w") as f:
            f.write(out)
    return out


# --------------------------------------------------------------------------
# 1. CIF -> czysty PDB receptora (ProDy zamiast recznego parsera z legacy)
def cif_to_receptor_pdb(cif_path, out_pdb):
    st = prody.parseMMCIF(cif_path, altloc=C.DEFAULT_ALTLOC)
    sel_str = "(protein or nucleic or resname %s) and not water" % " ".join(C.METAL_IONS)
    sel = st.select(sel_str)
    if sel is None:
        raise RuntimeError("pusta selekcja receptora dla %s" % cif_path)
    prody.writePDB(out_pdb, sel)
    return sel.numAtoms()


def residues_in_pdbqt(path):
    """Zbior (chain, resnum) obecnych w pliku PDB/PDBQT."""
    res = set()
    with open(path) as f:
        for line in f:
            if line[:6] in ("ATOM  ", "HETATM"):
                res.add((line[21:22].strip(), line[22:26].strip()))
    return res


def read_pdbqt_atoms(path):
    """Lista atomow (name, resname, chain, resnum, x, y, z, adtype)."""
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
# 2. Ligand: SDF (same ciezkie atomy, wspolrzedne krystaliczne)
#    -> molscrub ustala stan protonacji przy pH 7.4 (odpowiednik "obabel -p 7.4")
#    -> stan protonacji przenosimy na szkielet krystaliczny (wspolrzedne z SDF)
#    -> RDKit AddHs(addCoords=True) dokłada wodory
#    -> Meeko zapisuje PDBQT z typami AD4
#    native.sdf = ta sama czasteczka w pozie krystalicznej (referencja do RMSD)
def _neutral_copy(mol):
    """Kopia bez ladunkow formalnych i jawnych H - do dopasowania grafu ciezkich
    atomow miedzy czasteczka sprotonowana a krystaliczna."""
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(True)
    m = rw.GetMol()
    Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    return m


def prepare_ligand(sdf_path, ligand_pdbqt, native_sdf):
    crystal = Chem.MolFromMolFile(sdf_path)          # ciezkie atomy, 3D z krysztalu
    if crystal is None:
        raise RuntimeError("RDKit nie wczytal %s" % sdf_path)
    crystal_xyz = crystal.GetConformer().GetPositions()   # Scrub potrafi zabrac konformer

    scrub = Scrub(ph_low=C.LIGAND_PH, ph_high=C.LIGAND_PH,
                  skip_tautomers=True, skip_ringfix=True)
    protonated = next(iter(scrub(Chem.Mol(crystal))))     # protonacja + wlasne 3D (odrzucamy)

    skeleton = Chem.RemoveAllHs(protonated)
    match = _neutral_copy(crystal).GetSubstructMatch(_neutral_copy(skeleton))
    if len(match) != skeleton.GetNumAtoms():
        raise RuntimeError("nie udalo sie dopasowac grafu sprotonowanego liganda "
                           "do struktury krystalicznej (%s)" % sdf_path)

    # wspolrzedne krystaliczne na szkielet sprotonowany
    conf = skeleton.GetConformer()
    for i, j in enumerate(match):
        x, y, z = crystal_xyz[j]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    lig = Chem.AddHs(skeleton, addCoords=True)
    Chem.MolToMolFile(lig, native_sdf)

    setups = MoleculePreparation().prepare(lig)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError("Meeko: %s" % err)
    with open(ligand_pdbqt, "w") as f:
        f.write(pdbqt)

    charge = Chem.GetFormalCharge(lig)
    return dict(smiles=Chem.MolToSmiles(Chem.RemoveHs(lig)), charge=charge,
                heavy=skeleton.GetNumAtoms(), atoms=lig.GetNumAtoms())


# --------------------------------------------------------------------------
# 3. Miejsce biasu: najkrotszy kontakt heteroatom(ligand)-heteroatom(receptor).
#    Typ 'acc' jesli atom liganda jest akceptorem (OA/NA/SA), inaczej 'don'.
#    Kryterium identyczne jak w wersji legacy.
def extract_bias_site(ligand_pdbqt, receptor_pdbqt):
    lig = read_pdbqt_atoms(ligand_pdbqt)
    rec = read_pdbqt_atoms(receptor_pdbqt)
    acc_types = ("OA", "NA", "SA")

    lig_polar = [a for a in lig if a["adtype"][0] in ("N", "O", "S")]
    rec_polar = [a for a in rec if a["adtype"][0] in ("N", "O")]

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
    return dict(x=la["x"], y=la["y"], z=la["z"], type=btype), report


# --------------------------------------------------------------------------
# 4. Raport: ktore reszty Meeko wyrzucil (niekompletne lancuchy boczne)
#    i jak daleko sa od natywnego liganda -> odpowiedz na pytanie z plan.md
#    "czy luka w strukturze dotyka kieszeni".
def dropped_residues_report(clean_pdb, receptor_pdbqt, native_sdf, out_txt, prep_log):
    st = prody.parsePDB(clean_pdb)
    lig = Chem.MolFromMolFile(native_sdf, removeHs=False)
    lig_xyz = lig.GetConformer().GetPositions()

    kept = residues_in_pdbqt(receptor_pdbqt)
    lines = []
    seen = set()
    for res in st.iterResidues():
        key = (str(res.getChid()).strip(), str(res.getResnum()))
        if key in kept or key in seen:
            continue
        seen.add(key)
        d = np.linalg.norm(res.getCoords()[:, None, :] - lig_xyz[None, :, :], axis=-1).min()
        lines.append((d, "%s%s (%s)  min. odleglosc od liganda %.1f A"
                      % (res.getResname(), res.getResnum(), res.getChid(), d)))
    lines.sort()

    with open(out_txt, "w") as f:
        f.write("reszty usuniete przez mk_prepare_receptor (--allow_bad_res):\n")
        if not lines:
            f.write("  brak - wszystkie reszty dopasowane do szablonow\n")
        for _, txt in lines:
            f.write("  " + txt + "\n")
        near = [t for d, t in lines if d <= 5.0]
        f.write("\nw promieniu 5 A od natywnego liganda: %d\n" % len(near))
        f.write("\n--- log mk_prepare_receptor (ostrzezenia) ---\n")
        for line in prep_log.splitlines():
            if "Template matching failed" in line or "Ignored due to" in line:
                f.write(line + "\n")
    return lines


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
    ligand_pdbqt = os.path.join(tdir, "ligand.pdbqt")
    native_sdf = os.path.join(tdir, "native.sdf")

    # --- ligand (najpierw, bo z niego liczymy pudelko) ---
    info = prepare_ligand(sdf, ligand_pdbqt, native_sdf)
    print("  ligand: %d atomow (%d ciezkich), ladunek %+d, SMILES %s"
          % (info["atoms"], info["heavy"], info["charge"], info["smiles"]))

    # --- receptor ---
    n = cif_to_receptor_pdb(cif, receptor_pdb)
    print("  receptor_clean.pdb: %d atomow (cala jednostka, bez wod/HETATM poza jonami)" % n)
    cmd = [C.MK_PREPARE_RECEPTOR, "-i", receptor_pdb, "-o", "receptor",
           "-p", "-g", "--write_pdb", "receptor_prepared.pdb",
           "--box_enveloping", "ligand.pdbqt", "--padding", C.PADDING,
           "--default_altloc", C.DEFAULT_ALTLOC]
    if C.ALLOW_BAD_RES:
        cmd.append("--allow_bad_res")
    prep_log = sh(cmd, cwd=tdir, log=os.path.join(tdir, "mk_prepare_receptor.log"))

    dropped = dropped_residues_report(receptor_pdb, receptor_pdbqt, native_sdf,
                                      os.path.join(tdir, "receptor_report.txt"), prep_log)
    print("  reszty usuniete (--allow_bad_res): %d, najblizsza %.1f A od liganda"
          % (len(dropped), dropped[0][0] if dropped else float("nan")))

    # --- mapy AD4 (gpf od Meeko) ---
    gpf = os.path.join(tdir, "receptor.gpf")
    box = dict(re.findall(r"^(npts|gridcenter|spacing)\s+(.*?)\s*$", open(gpf).read(), re.M))
    with open(os.path.join(tdir, "box.txt"), "w") as f:
        f.write("center = %s\nnpts   = %s\nspacing= %s\n"
                % (box.get("gridcenter"), box.get("npts"), box.get("spacing")))
    print("  box center=%s npts=%s" % (box.get("gridcenter"), box.get("npts")))
    sh([C.AUTOGRID4, "-p", "receptor.gpf", "-l", "receptor.glg"], cwd=tdir)

    # --- bias: miejsce + mapy biasowane ---
    site, report = extract_bias_site(ligand_pdbqt, receptor_pdbqt)
    with open(os.path.join(tdir, "bias_report.txt"), "w") as f:
        f.write(report + "\n")
    print("  bias: " + report)
    if site is None:
        raise RuntimeError("%s: nie wyznaczono miejsca biasu" % pdbid)

    with open(os.path.join(tdir, "bias.bpf"), "w") as f:
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
def write_tool_versions():
    import meeko, rdkit, molscrub, vina
    lines = [
        "python        %s" % sys.version.split()[0],
        "meeko         %s" % meeko.__version__,
        "molscrub      %s" % getattr(molscrub, "__version__", "0.2.2"),
        "rdkit         %s" % rdkit.__version__,
        "prody         %s" % prody.__version__,
        "vina (API)    %s" % vina.__version__,
        "autogrid4     %s" % subprocess.run([C.AUTOGRID4, "--version"], stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT).stdout.decode().split("\n")[0],
        "prepare_bias  %s (python2)" % C.PREPARE_BIAS,
    ]
    path = os.path.join(C.RESULTS, "tool_versions.txt")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    which = sys.argv[1:] or list(C.TARGETS)
    os.makedirs(C.RESULTS, exist_ok=True)
    for pdbid in which:
        prepare_target(pdbid, C.TARGETS[pdbid])
    write_tool_versions()
    print("\nGotowe (prepare).")


if __name__ == "__main__":
    main()
