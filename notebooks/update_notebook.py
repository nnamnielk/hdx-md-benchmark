import json
import sys

nb_path = "/home/okleinmann/projects/upside/hdx-md-benchmark/notebooks/shaw_style_upside_analysis.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

code = """
import scipy.stats

def plot_nug2b_delta_g_hx_scatter(protein: str = 'nug2b', save: bool = True) -> plt.Figure:
    if protein != 'nug2b':
        raise ValueError('Skinner Table S1 HDX values are currently encoded only for NuG2b.')
    rows = [row for row in read_upside_hdx_csv(protein) if np.isfinite(row['skinner_delta_g_hx'])]
    
    exp_y = np.array([row['skinner_delta_g_hx'] for row in rows], dtype=float)
    sim_y = np.array([row['upside_delta_g_hx'] for row in rows], dtype=float)
    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(exp_y, sim_y)
    
    fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
    ax.scatter(exp_y, sim_y, color='#d95f02', s=30, label='Residue HDX')
    
    x_range = np.array([min(exp_y)-0.5, max(exp_y)+0.5])
    ax.plot(x_range, slope*x_range + intercept, 'k--', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')
    ax.plot(x_range, x_range, 'gray', linestyle=':', label='y = x')
    
    ax.set_xlabel(r'Skinner Experiment $\\Delta G_{HX}$ (kcal/mol)')
    ax.set_ylabel(r'Upside Simulation $\\Delta G_{HX}$ (kcal/mol)')
    ax.set_title(f'NuG2b HDX Parity\\n$R = {r_value:.3f}$, $R^2 = {r_value**2:.3f}$')
    
    ax.grid(ls=':', lw=0.4, alpha=0.45)
    ax.legend(frameon=False, loc='best')
    
    if save:
        fig.savefig(FIG_ROOT / f'{protein}_delta_g_hx_scatter.png', bbox_inches='tight')
    return fig

plot_nug2b_delta_g_hx_scatter('nug2b')
"""

new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [line + "\\n" for line in code.strip().split("\\n")]
}
# Remove the last newline from the last string
if new_cell["source"]:
    new_cell["source"][-1] = new_cell["source"][-1].rstrip("\\n")

nb["cells"].append(new_cell)

with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1)

print("Notebook updated.")
