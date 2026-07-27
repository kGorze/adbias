from pathlib import Path

HERE = Path(__file__).resolve().parent

# Python 3 toolchain: Meeko, molscrub, RDKit, ProDy, Vina, SciPy i pandas.
PY_NEW = "/home/kgorzelanczyk/miniforge3/envs/adbias-py3/bin/python"

DATA = "/home/kgorzelanczyk/adbias/benchmark/bialka"
RESULTS = HERE.parent / "results"

SYSTEMS = {
    "3CS9": {
        "path": RESULTS / "3CS9"
    },
    "5N9R": {
        "path": RESULTS / "5N9R"
    },
    "6JQR": {
        "path": RESULTS / "6JQR"
    }
}
