import os

HERE = os.path.dirname(os.path.abspath(__file__))

TOOLS = os.path.join(HERE, "tools")                    # skopiowane AutoDockTools_py3 + prepare_bias.py
PY3   = "/home/kgorzelanczyk/miniforge3/envs/bo-mmpbsa/bin/python"   # numpy+rdkit+obabel, dodatkowo uruchamia skrypty AutoDockTools_py3
PY2   = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/python2.7"   # tylko dla prepare_bias.py, który używa python2

AUTOGRID4 = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/autogrid4"
VINA  = "/home/kgorzelanczyk/AutoDock-Vina/build/linux/release/vina"

OBABEL = "/usr/bin/obabel"
OBRMS  = "/usr/bin/obrms"

# skrypty z ADT py3 (uruchamiane przez PY3 z PYTHONPATH=TOOLS)
#skrót ze zmienną pomocniczną przechowującą ścieżkę do katalogu
U24 = os.path.join(TOOLS, "AutoDockTools", "Utilities24")

#scieżki do skryptów
PREPARE_RECEPTOR = os.path.join(U24, "prepare_receptor4.py")
PREPARE_LIGAND   = os.path.join(U24, "prepare_ligand4.py")
PREPARE_GPF      = os.path.join(U24, "prepare_gpf4.py")
PREPARE_BIAS     = os.path.join(TOOLS, "prepare_bias.py")

DATA    = "/home/kgorzelanczyk/adbias/benchmark/bialka"
RESULTS = os.path.join(HERE, "results")

TARGETS = {
    "6JQR": {"cif": "6JQR/6JQR.cif","ligand_sdf": "6JQR/6jqr_B_C6F.sdf"},
    "3CS9": {"cif": "3CS9/3CS9.cif","ligand_sdf": "3CS9/3cs9_E_NIL.sdf"},
    "5N9R": {"cif": "5N9R/5N9R.cif","ligand_sdf": "5N9R/5n9r_K_8RN.sdf"},
}


WATER_RESNAMES = set(["HOH", "WAT", "DOD", "H2O"])
METAL_IONS = set([
    "MG", "ZN", "MN", "CA", "NA", "K", "FE", "FE2", "CU", "CU1", "CO", "NI",
    "CD", "HG", "MO", "W", "SR", "BA", "CS", "LI", "RB", "AL", "GA", "PB", "3CO",
])

#standardowy spacing jak w autodock
SPACING = 0.375

#padding jako kostka przeszukiwań
PADDING = 8.0

#punkty na osi autogridu
NPTS_MAX = 100

# kcal/mol (glebokosc studni, jak w przykladzie referencyjnym)
BIAS_VSET = -1.5
# A
BIAS_RADIUS = 1.0
#A, kryterium kontaktu heteroatom-heteroatom
HBOND_CUTOFF = 3.5

# --- docking ---------------------------------------------------------------
# wspolna, sparowana pula seedow: te same seedy dla conventional i biased
SEEDS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
         11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
EXHAUSTIVENESS = 16
NUM_MODES = 9
ENERGY_RANGE = 3
CPU = 8

# --- ocena -----------------------------------------------------------------
RMSD_SUCCESS = 2.0    # A, prog sukcesu redocking (top-1 i best-of-N)

CONDITIONS = ["conventional", "biased"]
