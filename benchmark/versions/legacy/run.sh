#!/usr/bin/env bash
# Caly benchmark legacy: przygotowanie -> docking -> ocena.
# Uzycie:  ./run.sh            (wszystkie cele)
#          ./run.sh 6JQR       (wybrany cel/cele)
set -e
cd "$(dirname "$0")"

# interpreter python3 bierzemy z config.py (jedno zrodlo prawdy)
PY3=$(python3 -c "import config; print(config.PY3)")

echo "########## 1/3 prepare ##########"
"$PY3" prepare.py "$@"
echo "########## 2/3 dock ##########"
"$PY3" dock.py "$@"
echo "########## 3/3 analyze ##########"
"$PY3" analyze.py
