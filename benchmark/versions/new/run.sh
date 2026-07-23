#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"

# interpreter bierzemy z config.py
PY=$(python3 -c "import config; print(config.PY_NEW)")

echo "prepare 1"
"$PY" prepare.py "$@"
echo "dock 2"
"$PY" dock.py "$@"
echo "analiza 3"
"$PY" analyze.py
echo "tcl 4"
"$PY" tcl_script.py
