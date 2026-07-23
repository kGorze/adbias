# -*- coding: utf-8 -*-
# Konfiguracja benchmarku "new" (nowoczesny protokol: Meeko + RDKit + ProDy + Vina API).
# Jeden slownik parametrow dla wszystkich skryptow (prepare / dock / analyze).
#
# WAZNE: wszystkie parametry naukowe (cele, pudelko, bias, seedy, docking, prog RMSD)
# sa CELOWO identyczne jak w ../legacy/config.py. Jedyna roznica miedzy wersjami
# "legacy" i "new" ma byc toolchain przygotowania (AutoDockTools vs Meeko/RDKit),
# a nie ustawienia eksperymentu.

import os

HERE = os.path.dirname(os.path.abspath(__file__))

# --- narzedzia / interpretery ---------------------------------------------
# srodowisko "vs": meeko 0.7.1, molscrub 0.2.2, rdkit, prody, vina (python API), scipy, pandas
PY_NEW = "/home/kgorzelanczyk/miniforge3/envs/vs/bin/python"
MK_PREPARE_RECEPTOR = "/home/kgorzelanczyk/miniforge3/envs/vs/bin/mk_prepare_receptor.py"

# te dwa sa wspolne dla obu wersji (to jest metoda, ktora testujemy, a nie toolchain)
PY2 = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/python2.7"   # tylko dla prepare_bias.py
AUTOGRID4 = "/home/kgorzelanczyk/miniforge3/envs/adbias/bin/autogrid4"
PREPARE_BIAS = os.path.join(HERE, "..", "legacy", "tools", "prepare_bias.py")
PREPARE_BIAS = os.path.abspath(PREPARE_BIAS)

# --- dane / wyniki ---------------------------------------------------------
DATA    = "/home/kgorzelanczyk/adbias/benchmark/bialka"
RESULTS = os.path.join(HERE, "results")

# cele: PDB ID -> pliki wejsciowe (cif krystaliczny + natywny ligand w sdf)
TARGETS = {
    "6JQR": {"cif": "6JQR/6JQR.cif", "ligand_sdf": "6JQR/6jqr_B_C6F.sdf"},
    "3CS9": {"cif": "3CS9/3CS9.cif", "ligand_sdf": "3CS9/3cs9_E_NIL.sdf"},
    "5N9R": {"cif": "5N9R/5N9R.cif", "ligand_sdf": "5N9R/5n9r_K_8RN.sdf"},
}

# --- receptor: co zostawiamy ----------------------------------------------
# cala jednostka asymetryczna, ale bez wod i bez HETATM poza jonami metali
# (w tych trzech celach jonow metali nie ma - zostaje samo bialko)
METAL_IONS = [
    "MG", "ZN", "MN", "CA", "NA", "K", "FE", "FE2", "CU", "CU1", "CO", "NI",
    "CD", "HG", "MO", "W", "SR", "BA", "CS", "LI", "RB", "AL", "GA", "PB", "3CO",
]

# reszty bez kompletu ciezkich atomow (nieuporzadkowane lancuchy boczne) sa przez
# mk_prepare_receptor usuwane (--allow_bad_res). prepare.py raportuje ktore i jak
# daleko od liganda, zeby bylo widac, czy luka dotyka kieszeni.
ALLOW_BAD_RES = True
DEFAULT_ALTLOC = "A"

# --- ligand ----------------------------------------------------------------
LIGAND_PH = 7.4       # stan protonacji (molscrub), odpowiednik "obabel -p 7.4" z legacy

# --- pudelko / mapy AD4 ----------------------------------------------------
SPACING = 0.375   # A, standard AutoDock
PADDING = 8.0     # A wokol natywnego liganda (z kazdej strony)

# --- bias (identyczny dla wszystkich celow) --------------------------------
BIAS_VSET = -1.5      # kcal/mol (glebokosc studni, jak w przykladzie referencyjnym)
BIAS_RADIUS = 1.0     # A
HBOND_CUTOFF = 3.5    # A, kryterium kontaktu heteroatom-heteroatom

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
