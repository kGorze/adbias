#!/usr/bin/env python3

"""
uzywany toolchain do tego:
    prody(parseMMCIF + select)
    mk_prepare_receptor.py (Meeko)
    molscrub (ph 7.4) + RDKit + Meeko
    mk_prepare_receptor --write_gpf
    autogrid4
    prepare_bias.py

Dla kazdego PDB ID z config.TARGETS w results/<ID>/:
    receptor_clean.pdb - jednostka asymetryczna z CIF, bez wod/HETATM poza jonami)
    receptor.pdbqt - receptor przygotowany przez Meeko
    ligand.pdbqt - (Meeko; protonacja pH 7.4 przez molscrub)
    native.sdf - (natywna poza krystaliczna z H = referencja do RMSD)
    receptor.gpf + mapy - (Meeko --write_gpf -> autogrid4)
    bias.bpf - miejsce biasu z geometrii kontaktu H-mostka
    maps_biased/ - (zestaw map z wpisanym biasem, do warunku "biased")
"""
import sys
import shutil
import subprocess
import meeko, rdkit, molscrub, vina
from pathlib import Path

import prody
from rdkit import Chem, RDLogger
from molscrub import Scrub
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy


from analyze import topology_query
import config as C

prody.confProDy(verbosity="critical")
RDLogger.DisableLog("rdApp.warning")


ACCEPTOR_TYPES = {"OA", "NA", "SA"}

def squared_distance(a, b):
    difference = a["xyz"] - b["xyz"]
    return float(difference @ difference)


def run_command(command, cwd):
    print("  command:", " ".join(str(argument) for argument in command))
    subprocess.run(
        [str(argument) for argument in command],
        cwd=cwd,
        check=True,
    )


# z prody zamist robienia tego ręcznie
def cif_to_receptor_pdb(cif_path, out_pdb):

    #parse mmcif
    st = prody.parseMMCIF(cif_path, altloc=C.DEFAULT_ALTLOC)
    selection_string = "(protein or nucleic or resname %s) and not water" % " ".join(C.METAL_IONS)
    selection = st.select(selection_string)

    #wypisanie pdb
    prody.writePDB(out_pdb, selection)
    return selection.numAtoms()



def prepare_ligand(sdf_path, ligand_pdbqt, native_sdf):
    # Wczytuje ligand z pliku SDF do obiektu RDKit Mol.
    crystal = Chem.MolFromMolFile(sdf_path)

    if crystal is None:
        raise RuntimeError("RDKit nie wczytal %s" % sdf_path)

    # Pobiera współrzędne XYZ wszystkich atomów z konformera struktury krystalicznej.
    crystal_xyz = crystal.GetConformer().GetPositions()
    #TODO: to czy bedziemy robili zakres czy tylko jedno pH jest zalezne od tego co chcemy robić, obecnie jest tylko jedno pH
    scrub = Scrub(
        ph_low=C.LIGAND_PH,
        ph_high=C.LIGAND_PH,
        skip_tautomers=True,
        skip_ringfix=True
    )
    # Protonuje kopię liganda i pobiera pierwszy wariant zwrócony przez MolScrub.
    protonated = next(iter(scrub(Chem.Mol(crystal))))

    # Usuwa wszystkie atomy wodoru, pozostawiając szkielet złożony z atomów ciężkich.
    skeleton = Chem.RemoveAllHs(protonated)

    # Dopasowuje atomy sprotonowanego szkieletu do atomów oryginalnej struktury krystalicznej na podstawie tylko atomów.
    match = crystal.GetSubstructMatch(topology_query(skeleton))
    if len(match) != skeleton.GetNumAtoms():
        raise RuntimeError("nie udalo sie dopasowac grafu sprotonowanego liganda "
                           "do struktury krystalicznej (%s)" % sdf_path)


    skeleton.GetConformer().SetPositions(crystal_xyz[list(match)])

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



# to nie jest wlasciwe miejsce na wyznaczanie biasu, tylko do testowania ze struktura referencyjna.
def extract_bias_site(ligand_pdbqt, receptor_pdbqt):
    ligand_atoms = PDBQTMolecule.from_file(
        ligand_pdbqt,
        skip_typing=True,
    ).atoms()

    receptor_atoms = PDBQTMolecule.from_file(
        receptor_pdbqt,
        skip_typing=True,
    ).atoms()

    ligand_polar_atoms = [
        atom
        for atom in ligand_atoms
        if atom["atom_type"].startswith(("N", "O", "S"))
    ]

    #TODO: sprawdzenie, bo polarne atomy w receptorze mogą być inne
    receptor_polar_atoms = [
        atom
        for atom in receptor_atoms
        if atom["atom_type"].startswith(("N", "O"))
    ]

    closest_pair = _find_closest_atom_pair(
        ligand_polar_atoms,
        receptor_polar_atoms,
    )

    if closest_pair is None:
        return None, "Brak atomów polarnych w ligandzie lub receptorze."

    distance_squared, ligand_atom, receptor_atom = closest_pair
    distance = distance_squared**0.5

    bias_type = (
        "acc"
        if ligand_atom["atom_type"] in ACCEPTOR_TYPES
        else "don"
    )

    warning = ""
    if distance > C.HBOND_CUTOFF:
        warning = (
            f" [UWAGA: odległość przekracza cutoff "
            f"{C.HBOND_CUTOFF:.1f} Å]"
        )

    report = (
        f"ligand atom {ligand_atom['name']} "
        f"({ligand_atom['atom_type']}) <-> "
        f"receptor {receptor_atom['resname']}"
        f"{receptor_atom['resid']}/"
        f"{receptor_atom['name']} "
        f"({receptor_atom['atom_type']}), "
        f"d={distance:.2f} Å -> bias '{bias_type}'"
        f"{warning}"
    )

    x, y, z = map(float, ligand_atom["xyz"])

    bias_site = {
        "x": x,
        "y": y,
        "z": z,
        "type": bias_type,
    }

    return bias_site, report


def _find_closest_atom_pair(ligand_atoms, receptor_atoms):
    closest_pair = None

    for ligand_atom in ligand_atoms:
        for receptor_atom in receptor_atoms:
            distance_squared = squared_distance(
                ligand_atom,
                receptor_atom,
            )

            if (
                closest_pair is None
                or distance_squared < closest_pair[0]
            ):
                closest_pair = (
                    distance_squared,
                    ligand_atom,
                    receptor_atom,
                )
    return closest_pair


# przygotowanie receptora
def _prepare_receptor(cif_path, target_dir):
    receptor_pdb = target_dir / "receptor_clean.pdb"
    atom_count = cif_to_receptor_pdb(str(cif_path), str(receptor_pdb))
    print(f"  receptor_clean.pdb: {atom_count} atomow")

    command = [
        C.MK_PREPARE_RECEPTOR,
        "-i", receptor_pdb,
        "-o", "receptor",
        "-p",
        "-g",
        "--box_enveloping", target_dir / "ligand.pdbqt",
        "--padding", C.PADDING,
        "--default_altloc", C.DEFAULT_ALTLOC,
    ]
    if C.ALLOW_BAD_RES:
        command.append("--allow_bad_res")

    run_command(command, target_dir)
    return target_dir / "receptor.pdbqt"


def _generate_ad4_maps(target_dir):
    run_command(
        [C.AUTOGRID4, "-p", "receptor.gpf", "-l", "receptor.glg"],
        target_dir,
    )


def _prepare_biased_maps(target_dir, ligand_pdbqt, receptor_pdbqt):
    site, description = extract_bias_site(ligand_pdbqt, receptor_pdbqt)
    if site is None:
        raise RuntimeError(description)
    print(f"  bias: {description}")

    bias_file = target_dir / "bias.bpf"
    bias_file.write_text(
        "x y z Vset r type\n"
        f"{site['x']:.3f} {site['y']:.3f} {site['z']:.3f} "
        f"{C.BIAS_VSET:.2f} {C.BIAS_RADIUS:.2f} {site['type']}\n"
    )
    run_command(
        [C.PY2, C.PREPARE_BIAS, "-b", bias_file.name, "-g", "receptor.gpf"],
        target_dir,
    )

    biased_directory = target_dir / "maps_biased"
    biased_directory.mkdir()
    for source in target_dir.glob("*.map"):
        if not source.name.endswith(".biased.map"):
            shutil.copy(source, biased_directory / source.name)
    for filename in ("receptor.maps.fld", "receptor.maps.xyz"):
        shutil.copy(target_dir / filename, biased_directory / filename)

    biased_maps = sorted(target_dir.glob("*.biased.map"))
    if not biased_maps:
        raise RuntimeError("prepare_bias.py nie utworzyl zadnej mapy biasowanej")
    for source in biased_maps:
        target_name = source.name.removesuffix(".biased.map") + ".map"
        source.replace(biased_directory / target_name)

    names = ", ".join(
        source.name.removesuffix(".biased.map") + ".map"
        for source in biased_maps
    )
    print(f"  mapy biased podmienione: {names}")


def prepare_target(pdbid, spec):
    print()
    print(f"=== {pdbid} ===")

    target_dir = Path(C.RESULTS) / pdbid
    if target_dir.exists():
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True)

    ligand_pdbqt = target_dir / "ligand.pdbqt"
    ligand_info = prepare_ligand(
        str(Path(C.DATA) / spec["ligand_sdf"]),
        str(ligand_pdbqt),
        str(target_dir / "native.sdf"),
    )
    print(
        f"  ligand: {ligand_info['atoms']} atomow "
        f"({ligand_info['heavy']} ciezkich), "
        f"ladunek {ligand_info['charge']:+d}, "
        f"SMILES {ligand_info['smiles']}"
    )

    receptor_pdbqt = _prepare_receptor(
        Path(C.DATA) / spec["cif"],
        target_dir,
    )
    _generate_ad4_maps(target_dir)
    _prepare_biased_maps(target_dir, ligand_pdbqt, receptor_pdbqt)


def write_tool_versions():
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
    path = Path(C.RESULTS) / "tool_versions.txt"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    which = sys.argv[1:] or list(C.TARGETS)
    Path(C.RESULTS).mkdir(parents=True, exist_ok=True)
    for pdbid in which:
        prepare_target(pdbid, C.TARGETS[pdbid])
    write_tool_versions()

if __name__ == "__main__":
    main()
