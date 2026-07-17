#!/usr/bin/env python3
"""Individual Skinner/PNAS-style charts per protein: Pf, dG_HX, Hbond native vs denat, Rg convergence."""
import sys, os
import numpy as np
import tables as tb
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups, load_upside_traj

mpl.rcParams.update({
    'font.size': 8, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.labelsize': 9, 'axes.titlesize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.7, 'ytick.major.width': 0.7,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'svg.fonttype': 'none',
})

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-benchmark-desmond-upside/runs")
OUT = os.path.expanduser("~/.hermes/workspace/hdx-benchmark-desmond-upside/analysis")

PROTEINS = [
    ("NuG2b (1MI0)",    "nug2b",            "#2166AC"),
    ("Ubiquitin (1UBQ)", "ubiquitin",        "#B2182B"),
    ("λ YA (1LMB:3)",   "lambda_repressor", "#4DAF4A"),
    ("λ D14A",          "lambda_d14a",       "#FF7F00"),
]

CRITERION = 0.01; R = 0.001987; T_K = 300.0

def collect_all(h5_path):
    with tb.open_file(h5_path, 'r') as t:
        seq = [s.decode() for s in t.root.input.sequence[:]]
        n_res = len(seq)
        hb_list, ridx_list = [], []
        for g_no, g in enumerate(_output_groups(t)):
            sl = slice(0, None) if g_no == 0 else slice(1, None)
            hb_list.append(g.hbond[sl, :n_res])
            ridx_list.append(g.replica_index[sl, 0])
        hb_all = np.concatenate(hb_list)
        ridx_all = np.concatenate(ridx_list)
    return seq, hb_all, ridx_all

def compute_dg(hb, ridx, t_idx):
    mask = ridx == t_idx
    if mask.sum() < 10: return None, None, 0
    pf = (hb[mask] > CRITERION).mean(axis=0)
    pf_c = np.clip(pf, 0.001, 0.999)
    dg = -R * T_K * np.log(pf_c / (1 - pf_c))
    return pf, dg, mask.sum()

def compute_rg(h5_path, stride=20):
    traj = load_upside_traj(h5_path, stride=stride)
    masses = np.array([a.element.mass for a in traj.topology.atoms], dtype=np.float64)
    total_mass = np.sum(masses)
    xyz = traj.xyz
    com = np.sum(xyz * masses[None, :, None], axis=1) / total_mass
    diff = xyz - com[:, None, :]
    rg = 10.0 * np.sqrt(np.sum(masses[None, :] * np.sum(diff**2, axis=2), axis=1) / total_mass)
    time_ns = traj.time / 1000.0
    return time_ns, rg, traj.n_residues

# ── Generate standalone charts per protein ─────────────────────
for label, subdir, color in PROTEINS:
    h5_main = os.path.join(RUNS, subdir, "outputs", f"{subdir}_remd", f"{subdir}.run.0.up")
    if not os.path.exists(h5_main):
        print(f"{label}: SKIP")
        continue
    
    seq, hb_all, ridx_all = collect_all(h5_main)
    pf_n, dg_n, n_nat = compute_dg(hb_all, ridx_all, 0)
    pf_d, dg_d, n_den = compute_dg(hb_all, ridx_all, 15)
    n_res = len(seq)
    residues = np.arange(1, n_res + 1)
    
    short = label.split(" (")[0]
    print(f"{short}: {n_nat} native, {n_den} denat frames")
    
    # ── Chart 1: Protection Factor ──
    fig, ax = plt.subplots(figsize=(7, 2.2))
    bar_colors = [color if pf > 0.5 else '#cccccc' for pf in pf_n]
    ax.bar(residues, pf_n, color=bar_colors, edgecolor='white', linewidth=0.3, width=0.85)
    ax.axhline(y=0.5, color='black', linestyle=':', alpha=0.35, linewidth=0.6)
    ax.set_ylabel("Protection Factor"); ax.set_ylim(0, 1.06)
    ax.set_xlim(0.5, n_res + 0.5)
    ax.set_title(f"{label} — Protection Factor (T=0.80, n={n_nat} frames)", fontweight='bold', fontsize=10, loc='left')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{subdir}_protection.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # ── Chart 2: ΔG_HX ──
    fig, ax = plt.subplots(figsize=(7, 2.2))
    colors_dg = [color if pf > 0.8 else ('#aaaaaa' if pf > 0.5 else '#dddddd') for pf in pf_n]
    ax.bar(residues, dg_n, color=colors_dg, edgecolor='white', linewidth=0.3, width=0.85)
    ax.axhline(y=0, color='black', linestyle='-', alpha=0.4, linewidth=0.5)
    ax.axhline(y=2, color='black', linestyle=':', alpha=0.2, linewidth=0.5)
    ax.set_ylabel("ΔG_HX (kcal/mol)")
    ax.set_title(f"{label} — ΔG_HX (T=0.80)", fontweight='bold', fontsize=10, loc='left')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{subdir}_dGHX.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # ── Chart 3: Native vs Denatured H-bond score ──
    fig, ax = plt.subplots(figsize=(7, 2.2))
    hb_nat = hb_all[ridx_all == 0].mean(axis=0)
    hb_den = hb_all[ridx_all == 15].mean(axis=0)
    ax.plot(residues, hb_nat, '-o', color=color, markersize=3, linewidth=1.2,
            markerfacecolor='white', label=f'Native (T=0.80, {n_nat} frames)')
    ax.plot(residues, hb_den, '-o', color='red',   markersize=3, linewidth=1.0,
            markerfacecolor='white', alpha=0.6, label=f'Denatured (T=0.95, {n_den} frames)')
    ax.set_ylabel("Mean H-bond Score"); ax.set_ylim(bottom=-0.02)
    ax.set_xlim(0.5, n_res + 0.5)
    ax.legend(fontsize=7, loc='upper right', frameon=False)
    ax.set_title(f"{label} — H-bond Score: Native vs Denatured", fontweight='bold', fontsize=10, loc='left')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{subdir}_hbond_nat_vs_denat.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # ── Chart 4: Rg convergence ──
    t_ns, rg, n_residues = compute_rg(h5_main, stride=10)
    fig, ax = plt.subplots(figsize=(7, 2.5))
    window = max(1, len(rg) // 20)
    rg_smooth = np.convolve(rg, np.ones(window)/window, mode='valid')
    t_smooth = t_ns[window-1:]
    ax.plot(t_ns, rg, alpha=0.15, color=color, linewidth=0.4)
    ax.plot(t_smooth, rg_smooth, color=color, linewidth=1.5, label='Running avg')
    half = len(rg) // 2
    rmean, rstd = rg[half:].mean(), rg[half:].std()
    ax.axhline(rmean, color='red', linestyle='--', alpha=0.6, linewidth=1, label=f'μ₂ = {rmean:.1f} ± {rstd:.1f} Å')
    ax.axhspan(rmean - rstd, rmean + rstd, alpha=0.06, color='red')
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Rg (Å)")
    ax.legend(fontsize=7, loc='best', frameon=False)
    ax.set_title(f"{label} — Rg Convergence ({n_residues} res, {len(rg)} frames)", fontweight='bold', fontsize=10, loc='left')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{subdir}_rg_convergence.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

print("\nAll charts saved to analysis/")
