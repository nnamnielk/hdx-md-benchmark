#!/bin/bash
source ~/miniconda3/bin/activate upside-hermes
export UPSIDE_HOME=$HOME/.hermes/sandbox/upside2-md

export OMP_NUM_THREADS=20

cd ~/.hermes/workspace/hdx-md-benchmark
export PYTHONPATH=$PYTHONPATH:$(pwd)/scripts

python -c '
import sys
import upside_replex

sys.argv = [
    "upside_replex",
    "--local",
    "--duration", "1000000",
    "--frame-interval", "100",
    "--replica-interval", "500",
    "--sim-id", "remd_lowtemp"
]
upside_replex.T_LOW = 0.60
upside_replex.T_HIGH = 0.80

upside_replex.run("nug2b")
'
