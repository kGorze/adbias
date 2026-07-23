from rdkit import Chem
import pandas as pd
from vina import Vina
import prody
import subprocess
import shutil
from pathlib import Path

receptor = '6JQR.pdb'

data = Chem.SDMolSupplier('molecules.sdf')
molecules = [elem for elem in data]
mol_len = len(molecules)



def find_program(program_name):
    program_path = shutil.which(program_name)
    if program_path is None:
        raise FileNotFoundError(f"{program_name} not found in PATH.")
    return Path(program_path)

"""
dwa aktywne składniki, 
sól z dużym przeciwjonem
kompleks z metalem
ligand jest mały i przeciwjon jest podobnej wielkości
"""

# jest to tylko heurystyka
# tutaj jest taki problem, ze sole, przeciwjony i inne odłączone fragmenty są usuwane w funkcji strip_salts() 
def strip_salts(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) == 1:
        return mol
    
    #można robić przez lambde, ale tak jest lepiej widoczne
    def count_atoms(frag):
        return frag.GetNumAtoms()
    
    return max(frags, key=count_atoms)


def dock(ligand_pdbqt, receptor_pdbqt, center, box_size, output_pdbqt, mol2_dir, exhaustiveness=8, n_poses=9):
    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(output_pdbqt, n_poses=n_poses, overwrite=True)

    mk_export = find_program("mk_export.py")
    output_sdf = output_pdbqt.replace(".pdbqt", ".sdf")
    subprocess.run(
        [mk_export, output_pdbqt, "-s", output_sdf],
        check=True
    )

    obabel = find_program("obabel")
    output_mol2 = f"{mol2_dir}/{Path(output_pdbqt).stem}.mol2"
    subprocess.run(
        [
            obabel,
            output_sdf,
            "-O",
            output_mol2,
        ],
        check=True
    )

    return v.energies(n_poses=1)[0][0]


def save_complex_to_pdb(complex_mol, filename):
    writer = Chem.PDBWriter(filename)
    writer.write(complex_mol)

def save_to_csv(results, filename):
    df = pd.DataFrame(results)
    df.to_csv(filename, index=False)

ligand_id = 0
ligand_smiles = Chem.MolToSmiles(strip_salts(molecules[0]))
pH = 6

ligands_dir = "/home/kgorzelanczyk/vs_dane/docking/ligands"
receptor_dir = "/home/kgorzelanczyk/vs_dane/docking/receptor"
results_dir = "/home/kgorzelanczyk/vs_dane/docking/results"
mol2_dir = "/home/kgorzelanczyk/vs_dane/docking/results/mol2"

args = []
skip_tautomer = True
if skip_tautomer:
    args.append("--skip_tautomers")
skip_acidbase = False
if skip_acidbase:
    args.append("--skip_acidbase")

def prepare_ligand(ligand_smiles, pH, args, id, out_dir):

    ligandPDBQT = f"ligand_{id}.pdbqt"
    ligand_pdbqt = Path(ligandPDBQT)

    ligandName = ligandPDBQT.replace(".pdbqt", "")
    ligandSDF = f"{ligandName}_scrubbed.sdf"

    scrub = find_program("scrub.py")
    mk_prepare_ligand = find_program("mk_prepare_ligand.py")
    

    subprocess.run(
        [
            scrub,
            ligand_smiles,
            "-o",
            f"{out_dir}/{ligandSDF}",
            "--ph",
            str(pH),
            *args,
        ],
        check=True
    )

    subprocess.run(
        [
            mk_prepare_ligand,
            "-i",
            f"{out_dir}/{ligandSDF}",
            "-o",
            f"{out_dir}/{ligand_pdbqt}",
        ],
        check=True
    )


def prepare_receptor(name, pdb, out_dir):
    atoms_from_pdb = prody.parsePDB(pdb)
    receptor_selection = "chain A and not water and not hetero"
    receptor_atoms = atoms_from_pdb.select(receptor_selection)
    prody_receptorPDB = f"{name}_receptor_atoms.pdb"
    prody.writePDB(f"{out_dir}/{prody_receptorPDB}", receptor_atoms)

    ligand_selection = "resname C6F"
    ligand_atoms = atoms_from_pdb.select(ligand_selection)
    center_x, center_y, center_z = prody.calcCenter(ligand_atoms)
    center = [center_x, center_y, center_z]

    prody_ligandPDB = f"{name}_ligand.pdb"
    prody.writePDB(f"{out_dir}/{prody_ligandPDB}", ligand_atoms)

    size_x = 20.0 
    size_y = 20.0 
    size_z = 20.0
    size = [size_x, size_y, size_z]
    prepare_output = f"{name}_receptorFH"

    mk_prepare_receptor = find_program("mk_prepare_receptor.py")

    subprocess.run(
        [
            mk_prepare_receptor,
            "-i",
            f"{out_dir}/{prody_receptorPDB}",
            "-o",
            f"{out_dir}/{prepare_output}",
            "-p",
            "--write_pdb",
            f"{out_dir}/{prepare_output}.pdb",
            "-v",
            "--box_center",
            str(center_x),
            str(center_y),
            str(center_z),
            "--box_size",
            str(size_x),
            str(size_y),
            str(size_z),
            "--allow_bad_res",
            "--default_altloc",
            "A",
        ],
        check=True
    )

    receptor_pdbqt = f"{out_dir}/{prepare_output}.pdbqt"
    return receptor_pdbqt, center, size





Path(receptor_dir).mkdir(parents=True, exist_ok=True)
Path(results_dir).mkdir(parents=True, exist_ok=True)
Path(mol2_dir).mkdir(parents=True, exist_ok=True)
receptor_pdbqt, center, box_size = prepare_receptor("6JQR", receptor, receptor_dir)

results = []
for i, mol in enumerate(molecules):
    # if i == 10:
    #     break
    ligand_smiles = Chem.MolToSmiles(strip_salts(mol))
    ligand_dir = f"{ligands_dir}/ligand_{i}"
    Path(ligand_dir).mkdir(parents=True, exist_ok=True)
    prepare_ligand(ligand_smiles, pH, args, i, ligand_dir)

    ligand_pdbqt = f"{ligand_dir}/ligand_{i}.pdbqt"
    output_pdbqt = f"{results_dir}/ligand_{i}_docked.pdbqt"
    score = dock(ligand_pdbqt, receptor_pdbqt, center, box_size, output_pdbqt, mol2_dir, exhaustiveness=2048, n_poses=9)
    results.append({"ligand_id": i, "smiles": ligand_smiles, "score": score})

save_to_csv(results, f"{results_dir}/results.csv")