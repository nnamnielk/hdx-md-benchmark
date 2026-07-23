import json
import sys

nb_path = "/home/okleinmann/projects/upside/hdx-md-benchmark/notebooks/shaw_style_upside_analysis.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

code = """
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_nug2b_comparisons(protein: str = 'nug2b', save: bool = True) -> plt.Figure:
    if protein != 'nug2b':
        raise ValueError('Only NuG2b supported.')
    rows = [row for row in read_upside_hdx_csv(protein) if np.isfinite(row['skinner_delta_g_hx'])]
    
    exp_g = np.array([row['skinner_delta_g_hx'] for row in rows], dtype=float)
    sim_g = np.array([row['upside_delta_g_hx'] for row in rows], dtype=float)
    sim_m = np.array([row['upside_m_value'] for row in rows], dtype=float)
    
    # Delta G Parity
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(exp_g, sim_g)
    print(f"Delta G_HX Parity: R^2 = {r_value**2:.3f}, Slope = {slope:.3f}")
    
    # M-value comparison
    exp_m_global = SKINNER_NUG2B_GLOBALS['m_hx']['value']
    mean_sim_m = np.mean(sim_m)
    print(f"M_HX Average: Simulated = {mean_sim_m:.3f}, Experimental Global = {exp_m_global:.3f}")
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    
    # Plot 1: Delta G
    ax = axes[0]
    ax.scatter(exp_g, sim_g, color='#d95f02', s=30)
    x_range = np.array([min(exp_g)-0.5, max(exp_g)+0.5])
    ax.plot(x_range, slope*x_range + intercept, 'k--', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')
    ax.plot(x_range, x_range, 'gray', linestyle=':', label='y = x')
    ax.set_xlabel(r'Skinner Experiment $\\Delta G_{HX}$ (kcal/mol)')
    ax.set_ylabel(r'Upside Simulation $\\Delta G_{HX}$ (kcal/mol)')
    ax.set_title(f'NuG2b $\\Delta G_{HX}$ Parity\\n$R = {r_value:.3f}$, $R^2 = {r_value**2:.3f}$')
    ax.grid(ls=':', lw=0.4, alpha=0.45)
    ax.legend(frameon=False)
    
    # Plot 2: M values
    ax = axes[1]
    residues = np.array([row['residue'] for row in rows], dtype=int)
    labels = [f"{row['residue']}{row['skinner_aa']}" for row in rows]
    ax.bar(residues, sim_m, color='#7570b3', alpha=0.8, label='Upside Simulated $m_{HX}$')
    ax.axhline(exp_m_global, color='black', linestyle='--', label=f'Exp Global $m_{HX}$ ({exp_m_global})')
    ax.set_xticks(residues)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_xlabel('NuG2b residue')
    ax.set_ylabel(r'$m_{HX}$ (kcal/mol/M)')
    ax.set_title('NuG2b $m_{HX}$ Comparison')
    ax.grid(axis='y', ls=':', lw=0.4, alpha=0.45)
    ax.legend(frameon=False)
    
    if save:
        fig.savefig(FIG_ROOT / f'{protein}_comparisons.png', bbox_inches='tight')
    return fig

plot_nug2b_comparisons('nug2b')
"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\\n" for line in code.strip().split("\\n")]
}
if new_cell["source"]:
    new_cell["source"][-1] = new_cell["source"][-1].rstrip("\\n")

nb["cells"].append(new_cell)

with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1)

print("Notebook updated.")
