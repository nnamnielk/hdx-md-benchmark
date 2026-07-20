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

# Find residues with very high protection (Pf > 0.99)
pf_full = (hb_native > CRITERION).mean(axis=0)
highly_protected_idx = np.where(pf_full > 0.99)[0]

print("\n=== Extremely Rare Opening Events (Pf > 0.99) ===")
if len(highly_protected_idx) == 0:
    print("No residues with Pf > 0.99 found at this temperature.")
else:
    print(f"{'Res':>4} {'AA':>4} {'Pf_Full':>9} {'Openings (frames)':>18} {'Est. P(open)':>15}")
    print("-" * 55)
    for i in highly_protected_idx:
        openings = np.sum(hb_native[:, i] <= CRITERION)
        p_open = openings / len(hb_native)
        print(f"{i+1:>4} {seq[i]:>4} {pf_full[i]:>9.5f} {openings:>17} / 10000 {p_open:>15.5e}")

# Check scaling for the top 5 most protected residues across temperature ladders
print("\n=== Temperature Scaling for Most Protected Residues ===")
top_5_idx = np.argsort(pf_full)[-5:][::-1]
temps = np.arange(20)

for i in top_5_idx:
    print(f"\nResidue {i+1} {seq[i]} (Native Pf: {pf_full[i]:.5f})")
    print("Ladder T | Openings / Frames | P(open)")
    print("-" * 40)
    for t_idx in temps:
        mask = (ridx_all == t_idx)
        n_frames = np.sum(mask)
        if n_frames > 0:
            hb_t = hb_all[mask, i]
            openings = np.sum(hb_t <= CRITERION)
            p_open = openings / n_frames
            print(f"{t_idx:>8} | {openings:>8} / {n_frames:<6} | {p_open:.5e}")
