#!/usr/bin/env python3
"""Skinner-et-al-2014-style protection factor + Rg analysis for all HDX benchmark proteins."""
import sys, os
import numpy as np
import tables as tb
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sb
from collections import defaultdict

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups

# ── Skinner/PNAS style ─────────────────────────────────────────
mpl.rcParams.update({
    'font.size': 9,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-benchmark-desmond-upside/runs")
OUT = os.path.expanduser("~/.hermes/workspace/hdx-benchmark-desmond-upside/analysis")

# Protein definitions
PROTEINS = [
    ("NuG2b (1MI0)", "nug2b", "#2196F3"),
    ("λ D14A", "lambda_d14a", "#FF9800"),
    ("λ YA (1LMB:3)", "lambda_repressor", "#4CAF50"),
    ("Ubiquitin (1UBQ)", "ubiquitin", "#E91E63"),
]

CRITERION = 0.01  # H-bond threshold
R = 0.001987
T_NATIVE = 300.0
N_REPLICAS = 16
T_IDX_NATIVE = 0    # T=0.80
T_IDX_DENATURED = 15  # T=0.95

def collect_hdx_data(h5_path):
    """Collect hbond, replica_index, Rg, and H-bond count across output groups."""
    with tb.open_file(h5_path, 'r') as t:
        seq = [s.decode() for s in t.root.input.sequence[:]]
        n_res = len(seq)
        
        hb_list, ridx_list, temp_list = [], [], []
        
        # Collect from all output groups (handles restarts)
        output_groups = list(_output_groups(t))
        n_groups = len(output_groups)
        
        if n_groups == 0:
            return None
        
        for g_no, g in enumerate(output_groups):
            nf = g.hbond.shape[0]
            sl = slice(0, None) if g_no == 0 else slice(1, None)
            hb_list.append(g.hbond[sl, :n_res])
            ridx_list.append(g.replica_index[sl, 0])
        
        hb_all = np.concatenate(hb_list)
        ridx_all = np.concatenate(ridx_list)
        
    return seq, hb_all, ridx_all

def compute_protection(hb, ridx, t_idx):
    """Compute protection factors and dG_HX for a given T-index."""
    mask = ridx == t_idx
    if mask.sum() < 10:
        return None, None, mask.sum()
    
    hb_sel = hb[mask]
    protected = hb_sel > CRITERION
    pf = protected.mean(axis=0)
    pf_c = np.clip(pf, 0.001, 0.999)
    dG = -R * T_NATIVE * np.log(pf_c / (1 - pf_c))
    return pf, dG, mask.sum()

def compute_per_residue_hbonds(hb, ridx, t_idx):
    """Mean H-bond occupancy per residue at given T-index."""
    mask = ridx == t_idx
    if mask.sum() < 10:
        return None
    return hb[mask].mean(axis=0)

# ── Collect data ────────────────────────────────────────────────
all_data = {}
for label, subdir, color in PROTEINS:
    h5_path = os.path.join(RUNS, subdir, "outputs", f"{subdir}_remd", f"{subdir}.run.0.up")
    if not os.path.exists(h5_path):
        print(f"{label}: no data, skipping")
        continue
    
    try:
        data = collect_hdx_data(h5_path)
        if data is None:
            print(f"{label}: empty, skipping")
            continue
        seq, hb_all, ridx_all = data
        
        pf_native, dG_native, n_native = compute_protection(hb_all, ridx_all, T_IDX_NATIVE)
        pf_denat, dG_denat, n_denat = compute_protection(hb_all, ridx_all, T_IDX_DENATURED)
        
        # Per-residue H-bond scores at native T
        hb_native = compute_per_residue_hbonds(hb_all, ridx_all, T_IDX_NATIVE)
        hb_denat = compute_per_residue_hbonds(hb_all, ridx_all, T_IDX_DENATURED)
        
        all_data[label] = {
            'seq': seq, 'color': color,
            'pf_native': pf_native, 'dG_native': dG_native,
            'pf_denat': pf_denat, 'dG_denat': dG_denat,
            'hb_native': hb_native, 'hb_denat': hb_denat,
            'n_native': n_native, 'n_denat': n_denat,
            'n_res': len(seq),
        }
        n_prot = (pf_native > 0.5).sum() if pf_native is not None else 0
        print(f"{label}: {n_native} native frames, {n_prot}/{len(seq)} protected residues")
        
    except Exception as e:
        print(f"{label}: ERROR - {e}")

# ── FIGURE 1: Protection factor per residue (all proteins) ─────
n_prots = len(all_data)
fig1, axes1 = plt.subplots(n_prots, 1, figsize=(7.2, 2.2 * n_prots), sharex=False)
if n_prots == 1:
    axes1 = [axes1]

for idx, (label, d) in enumerate(all_data.items()):
    ax = axes1[idx]
    residues = np.arange(1, d['n_res'] + 1)
    
    if d['pf_native'] is not None:
        colors_bar = [d['color'] if pf > 0.5 else '#cccccc' for pf in d['pf_native']]
        ax.bar(residues, d['pf_native'], color=colors_bar, edgecolor='white', linewidth=0.3, width=0.9)
    
    ax.axhline(y=0.5, color='black', linestyle=':', alpha=0.4, linewidth=0.7)
    ax.set_ylabel("Protection Factor")
    ax.set_ylim(0, 1.08)
    ax.set_xlim(0.5, d['n_res'] + 0.5)
    ax.set_title(label, fontweight='bold', fontsize=10, loc='left', pad=4)
    
    # Axis styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=3)

axes1[-1].set_xlabel("Residue Number")
fig1.suptitle("Protection Factors — Native Ensemble (T = 0.80)", fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
fig1.savefig(os.path.join(OUT, "skinner_fig_protection.png"), dpi=300, bbox_inches='tight')
print("\nSaved skinner_fig_protection.png")

# ── FIGURE 2: ΔG_HX per residue ────────────────────────────────
fig2, axes2 = plt.subplots(n_prots, 1, figsize=(7.2, 2.2 * n_prots), sharex=False)
if n_prots == 1:
    axes2 = [axes2]

for idx, (label, d) in enumerate(all_data.items()):
    ax = axes2[idx]
    residues = np.arange(1, d['n_res'] + 1)
    
    if d['dG_native'] is not None:
        dG = d['dG_native']
        # Color: blue = protected (dG > 2), gray = marginal, white/red = exposed
        bar_colors = []
        for val, pf in zip(dG, d['pf_native']):
            if pf > 0.8:
                bar_colors.append(d['color'])
            elif pf > 0.5:
                bar_colors.append('#aaaaaa')
            else:
                bar_colors.append('#dddddd')
        ax.bar(residues, dG, color=bar_colors, edgecolor='white', linewidth=0.3, width=0.9)
    
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=0.5)
    ax.axhline(y=2, color='black', linestyle=':', alpha=0.3, linewidth=0.5)
    ax.set_ylabel("ΔG_HX (kcal/mol)")
    ax.set_title(label, fontweight='bold', fontsize=10, loc='left', pad=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=3)

axes2[-1].set_xlabel("Residue Number")
fig2.suptitle("ΔG_HX — Native Ensemble (T = 0.80)", fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
fig2.savefig(os.path.join(OUT, "skinner_fig_dGHX.png"), dpi=300, bbox_inches='tight')
print("Saved skinner_fig_dGHX.png")

# ── FIGURE 3: H-bond score comparison (Native vs Denatured) ────
fig3, axes3 = plt.subplots(n_prots, 1, figsize=(7.2, 2.2 * n_prots), sharex=False)
if n_prots == 1:
    axes3 = [axes3]

for idx, (label, d) in enumerate(all_data.items()):
    ax = axes3[idx]
    residues = np.arange(1, d['n_res'] + 1)
    
    if d['hb_native'] is not None:
        ax.plot(residues, d['hb_native'], '-o', color=d['color'], markersize=3, 
                linewidth=1.2, label='Native (T=0.80)', markerfacecolor='white')
    if d['hb_denat'] is not None:
        ax.plot(residues, d['hb_denat'], '-o', color='red', markersize=3, 
                linewidth=1.2, label='Denatured (T=0.95)', markerfacecolor='white', alpha=0.7)
    
    ax.set_ylabel("H-bond Score")
    ax.set_title(label, fontweight='bold', fontsize=10, loc='left', pad=4)
    ax.legend(fontsize=7, loc='upper right', frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=3)

axes3[-1].set_xlabel("Residue Number")
fig3.suptitle("Per-Residue H-Bond Scores — Native vs Denatured Ensemble", fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
fig3.savefig(os.path.join(OUT, "skinner_fig_hbond_native_vs_denatured.png"), dpi=300, bbox_inches='tight')
print("Saved skinner_fig_hbond_native_vs_denatured.png")

# ── FIGURE 4: Summary table ────────────────────────────────────
print("\n=== SUMMARY ===")
print(f"{'Protein':<25} {'Res':>4} {'Protected':>10} {'Exposed':>8} {'Mean Pf':>8} {'Mean dG_HX':>11}")
print("-" * 70)
for label, d in all_data.items():
    if d['pf_native'] is not None:
        n_prot = (d['pf_native'] > 0.5).sum()
        n_exp = d['n_res'] - n_prot
        print(f"{label:<25} {d['n_res']:>4} {n_prot:>10} {n_exp:>8} {d['pf_native'].mean():>7.3f} {d['dG_native'].mean():>10.2f}")

print("\nDone!")
