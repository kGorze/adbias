#!/usr/bin/env bash
# Caly benchmark "new": przygotowanie -> docking -> ocena.
# Uzycie:  ./run.sh            (wszystkie cele)
#          ./run.sh 6JQR       (wybrany cel/cele)
set -e
cd "$(dirname "$0")"

# interpreter bierzemy z config.py (jedno zrodlo prawdy)
PY=$(python3 -c "import config; print(config.PY_NEW)")

echo "########## 1/3 prepare ##########"
"$PY" prepare.py "$@"
echo "########## 2/3 dock ##########"
"$PY" dock.py "$@"
echo "########## 3/3 analyze ##########"
"$PY" analyze.py
