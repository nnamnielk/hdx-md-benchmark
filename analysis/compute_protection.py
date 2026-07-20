#!/usr/bin/env python3
"""Compute NuG2b protection factors and ΔG_HX from REMD data."""
import sys, os, glob
import numpy as np
import tables as tb
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sb

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups

mpl.rcParams['font.size'] = 13
mpl.rcParams['font.family'] = "sans-serif"

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/simulations")
OUT = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/analysis")

h5_dir = os.path.join(RUNS, "nug2b/outputs/REMD")
files = glob.glob(os.path.join(h5_dir, "*.up"))

all_hb = []
all_ridx = []
seq = None
n_res = 0

for f in files:
    with tb.open_file(f, 'r') as t:
        if seq is None:
            seq = [s.decode() for s in t.root.input.sequence[:]]
            n_res = len(seq)
        
        for g_no, g in enumerate(_output_groups(t)):
            if g_no == 0:
                all_hb.append(g.hbond[:])
                all_ridx.append(g.replica_index[:, 0])
            else:
                all_hb.append(g.hbond[1:])
                all_ridx.append(g.replica_index[1:, 0])

hb_all = np.concatenate(all_hb)
ridx_all = np.concatenate(all_ridx)

print(f"Total frames: {len(hb_all)}")
print(f"Sequence ({n_res} res): {seq[:5]}...{seq[-3:]}")

# Filter native ensemble: replica_index == 0
native_mask = ridx_all == 0
native_frames = native_mask.sum()
print(f"Native ensemble (T=0.80): {native_frames}/{len(hb_all)} frames ({100*native_frames/len(hb_all):.1f}%)")

# Backbone NH H-bond scores (first n_res values)
hb_native = hb_all[native_mask, :n_res]

# Protection: NH is protected if hbond > 0.01
CRITERION = 0.01
protected = hb_native > CRITERION
pf = protected.mean(axis=0)

# ΔG_HX = -RT ln(Pf / (1-Pf))
R = 0.001987
T = 300.0
pf_clipped = np.clip(pf, 0.001, 0.999)
dG_HX = -R * T * np.log(pf_clipped / (1 - pf_clipped))

# Convergence Check: Split native frames into two halves
half = len(hb_native) // 2
pf_h1 = (hb_native[:half] > CRITERION).mean(axis=0)
pf_h2 = (hb_native[half:] > CRITERION).mean(axis=0)

dG_h1 = -R * T * np.log(np.clip(pf_h1, 0.001, 0.999) / (1 - np.clip(pf_h1, 0.001, 0.999)))
dG_h2 = -R * T * np.log(np.clip(pf_h2, 0.001, 0.999) / (1 - np.clip(pf_h2, 0.001, 0.999)))

rmse = np.sqrt(np.mean((dG_h1 - dG_h2)**2))
print(f"Convergence Check (ΔG_HX): RMSE between first half and second half = {rmse:.3f} kcal/mol")
if rmse < 0.5:
    print("-> Sampling looks SUFFICIENT (RMSE < 0.5 kcal/mol).")
else:
    print("-> More sampling might be NEEDED (RMSE >= 0.5 kcal/mol).")

print(f"\n=== Protection Factors & dG_HX ===")
print(f"{'Res':>4} {'AA':>4} {'Pf':>7} {'dG_HX':>8}")
print("-" * 28)
for i in range(n_res):
    print(f"{i+1:>4} {seq[i]:>4} {pf[i]:>6.3f} {dG_HX[i]:>7.2f}")

print(f"\nMean Pf: {pf.mean():.3f}")
print(f"Mean dG_HX: {dG_HX.mean():.2f} kcal/mol")
print(f"Residues with Pf > 0.5: {np.sum(pf > 0.5)}/{n_res}")
print(f"Residues with Pf > 0.9: {np.sum(pf > 0.9)}/{n_res}")

# ========== PLOT: ΔG_HX vs Residue ==========
fig, ax = plt.subplots(figsize=(12, 5))

residues = np.arange(1, n_res + 1)

# Bar plot of ΔG_HX
colors = ['#2196F3' if pf[i] > 0.5 else '#FF9800' for i in range(n_res)]
ax.bar(residues, dG_HX, color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)

# Add residue labels for structured regions
for i in range(n_res):
    if pf[i] > 0.5:
        ax.text(i+1, dG_HX[i] + 0.3, seq[i], ha='center', fontsize=7, rotation=90, color='#1565C0')

ax.set_xlabel("Residue")
ax.set_ylabel("ΔG_HX (kcal/mol)")
ax.set_title("NuG2b — Upside 2.0 REMD Protection Factors\nff_2.1, T=0.80 (native ensemble), 200k steps, 16 replicas", fontweight='bold')

# Horizontal line at ΔG = 0
ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#2196F3', alpha=0.8, label=f'Protected (Pf > 0.5, n={np.sum(pf>0.5)})'),
    Patch(facecolor='#FF9800', alpha=0.8, label=f'Exposed (Pf ≤ 0.5, n={np.sum(pf<=0.5)})'),
]
ax.legend(handles=legend_elements, fontsize=10, loc='upper right')

sb.despine()
plt.tight_layout()
fig.savefig(os.path.join(OUT, "nug2b_dGHX.png"), dpi=150, bbox_inches='tight')
print(f"\nSaved {OUT}/nug2b_dGHX.png")

# ========== PLOT: Protection factor vs residue ==========
fig2, ax2 = plt.subplots(figsize=(12, 5))
ax2.bar(residues, pf, color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)
ax2.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Pf = 0.5')
ax2.set_xlabel("Residue")
ax2.set_ylabel("Protection Factor (Pf)")
ax2.set_title("NuG2b — Protection Factor per Residue", fontweight='bold')
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=10)
sb.despine()
plt.tight_layout()
fig2.savefig(os.path.join(OUT, "nug2b_protection_factor.png"), dpi=150, bbox_inches='tight')
print(f"Saved {OUT}/nug2b_protection_factor.png")
