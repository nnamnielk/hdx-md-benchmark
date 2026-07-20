import sys, os, glob
import numpy as np
import tables as tb

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/simulations")
h5_dir = os.path.join(RUNS, "nug2b/outputs/REMD")
files = glob.glob(os.path.join(h5_dir, "*.up"))

all_hb = []
all_ridx = []
seq = None

for f in files:
    with tb.open_file(f, 'r') as t:
        if seq is None:
            seq = [s.decode() for s in t.root.input.sequence[:]]
        
        for g_no, g in enumerate(_output_groups(t)):
            if g_no == 0:
                all_hb.append(g.hbond[:])
                all_ridx.append(g.replica_index[:, 0])
            else:
                all_hb.append(g.hbond[1:])
                all_ridx.append(g.replica_index[1:, 0])

hb_all = np.concatenate(all_hb)
ridx_all = np.concatenate(all_ridx)

native_mask = (ridx_all == 0)
hb_native = hb_all[native_mask, :len(seq)]
CRITERION = 0.01

half = len(hb_native) // 2
pf_h1 = (hb_native[:half] > CRITERION).mean(axis=0)
pf_h2 = (hb_native[half:] > CRITERION).mean(axis=0)
pf_full = (hb_native > CRITERION).mean(axis=0)

diff = np.abs(pf_h1 - pf_h2)

print("\n=== Protection Factor Convergence by Residue ===")
print(f"{'Res':>4} {'AA':>4} {'Pf_H1':>7} {'Pf_H2':>7} {'Diff':>7} {'Category':>12}")
print("-" * 50)

categories = {'High (Pf > 0.8)': [], 'Medium (0.2 < Pf <= 0.8)': [], 'Low (Pf <= 0.2)': []}
for i in range(len(seq)):
    if pf_full[i] > 0.8: cat = 'High (Pf > 0.8)'
    elif pf_full[i] > 0.2: cat = 'Medium (0.2 < Pf <= 0.8)'
    else: cat = 'Low (Pf <= 0.2)'
    
    categories[cat].append(diff[i])
    print(f"{i+1:>4} {seq[i]:>4} {pf_h1[i]:>7.3f} {pf_h2[i]:>7.3f} {diff[i]:>7.3f} {cat:>20}")

print("\n=== Average Absolute Difference by Category ===")
for cat, diffs in categories.items():
    if diffs:
        print(f"{cat:<25}: {np.mean(diffs):.4f}  (max: {np.max(diffs):.4f}, n={len(diffs)})")
