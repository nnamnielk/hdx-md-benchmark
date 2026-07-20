#!/bin/bash
source ~/miniconda3/bin/activate upside-hermes
export UPSIDE_HOME=$HOME/.hermes/sandbox/upside2-md

# 20 replicas = 20 threads (1 thread per replica)
export OMP_NUM_THREADS=20

cd ~/.hermes/workspace/hdx-md-benchmark
export PYTHONPATH=$PYTHONPATH:$(pwd)/scripts

# Clean start (no --restart) to ensure we don't load corrupted nan frames
python -c '
import sys
import upside_replex
sys.argv = [
    "upside_replex",
    "--local",
    "--duration", "1000000",
    "--frame-interval", "100",
    "--replica-interval", "500"
]
upside_replex.run("nug2b")
'
