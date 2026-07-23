# -*- coding: utf-8 -*-
import os

#jeden słownik parametrów. testujemy tylko meeko, rdkit, prody i vine zamiast legacy tools 

HERE = os.path.dirname(os.path.abspath(__file__))

# Python 3 toolchain: Meeko, molscrub, RDKit, ProDy, Vina, SciPy i pandas.
PY_NEW = "/home/kgorzelanczyk/miniforge3/envs/adbias-py3/bin/python"
MK_PREPARE_RECEPTOR = "/home/kgorzelanczyk/miniforge3/envs/adbias-py3/bin/mk_prepare_receptor.py"

# te dwa sa wspolne dla obu wersji (to jest metoda, ktora testujemy, a nie toolchain)
PY2 = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/python2.7"
AUTOGRID4 = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/autogrid4"
PREPARE_BIAS = os.path.join(HERE, "..", "legacy", "tools", "prepare_bias.py")
PREPARE_BIAS = os.path.abspath(PREPARE_BIAS)

# ideal_interaction_sites.py tez jest Python2 (Bio.PDB + numpy), uzywa tego samego PY2 co prepare_bias.py
IDEAL_INTERACTION_SITES = os.path.join(HERE, "..", "..", "..", "adbias", "code", "ideal_interaction_sites.py")
IDEAL_INTERACTION_SITES = os.path.abspath(IDEAL_INTERACTION_SITES)

DATA    = "/home/kgorzelanczyk/adbias/benchmark/bialka"
RESULTS = os.path.join(HERE, "results")

# cele: PDB ID -> pliki wejsciowe (cif krystaliczny + natywny ligand w sdf)
TARGETS = {
    "6JQR": {"cif": "6JQR/6JQR.cif", "ligand_sdf": "6JQR/6jqr_B_C6F.sdf"},
    "3CS9": {"cif": "3CS9/3CS9.cif", "ligand_sdf": "3CS9/3cs9_E_NIL.sdf"},
    "5N9R": {"cif": "5N9R/5N9R.cif", "ligand_sdf": "5N9R/5n9r_K_8RN.sdf"},
}

#receptor co jest zostawiane
#nie ma w tych celach metali ale do innych warto zostawić
METAL_IONS = [
    "MG", "ZN", "MN", "CA", "NA", "K", "FE", "FE2", "CU", "CU1", "CO", "NI",
    "CD", "HG", "MO", "W", "SR", "BA", "CS", "LI", "RB", "AL", "GA", "PB", "3CO",
]

#usuwanie reszt bez atomów ciężkich(łańcuchy boczne)
ALLOW_BAD_RES = True
DEFAULT_ALTLOC = "A"

#ligand
LIGAND_PH = 7.4

#box
SPACING = 0.375
PADDING = 8.0

#bias
BIAS_VSET = -1.5
BIAS_RADIUS = 1.0
HBOND_CUTOFF = 3.5

#docking
SEEDS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
         11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
EXHAUSTIVENESS = 32
NUM_MODES = 20
ENERGY_RANGE = 10
CPU = 8

#warunek sukcesu redockingu
RMSD_SUCCESS = 2.0

CONDITIONS = ["conventional", "biased"]
