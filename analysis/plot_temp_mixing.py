import sys, os, glob
import numpy as np
import tables as tb
import matplotlib.pyplot as plt
import seaborn as sb

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/simulations")
OUT = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/analysis")
h5_dir = os.path.join(RUNS, "nug2b/outputs/REMD")
files = glob.glob(os.path.join(h5_dir, "*.up"))

all_ridx = []
# Just trace replica 0 and replica 19 to see how they walk the temperature ladder
rep0_history = []
rep19_history = []

for f in sorted(files):
    if "nug2b.run.0.up" in f:
        with tb.open_file(f, 'r') as t:
            for g_no, g in enumerate(_output_groups(t)):
                if g_no == 0: rep0_history.extend(g.replica_index[:, 0])
                else: rep0_history.extend(g.replica_index[1:, 0])
    elif "nug2b.run.19.up" in f:
        with tb.open_file(f, 'r') as t:
            for g_no, g in enumerate(_output_groups(t)):
                if g_no == 0: rep19_history.extend(g.replica_index[:, 0])
                else: rep19_history.extend(g.replica_index[1:, 0])

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(rep0_history[:1000], label='Physical Replica 0', alpha=0.8, color='blue', linewidth=1)
ax.plot(rep19_history[:1000], label='Physical Replica 19', alpha=0.8, color='red', linewidth=1)
ax.set_yticks(np.arange(0, 20, 2))
ax.set_ylabel("Temperature Ladder Index")
ax.set_xlabel("Frame (first 1000 shown)")
ax.set_title("REMD Mixing: Temperature Random Walk")
ax.legend()
sb.despine()
plt.tight_layout()
fig.savefig(os.path.join(OUT, "nug2b_remd_mixing_trace.png"), dpi=150)
