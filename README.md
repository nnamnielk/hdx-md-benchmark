# HDX Benchmark: DESMOND vs Upside MD

Compare hydrogen-deuterium exchange (HDX) predictions between DESMOND (all-atom MD) and Upside 2.0 (coarse-grained MD), using the HX benchmarking framework established by Skinner, Yu, Sosnick et al. (2014).

## Structure

```
papers/     — Relevant literature
data/       — Protein cross-reference and experimental data
```

## Papers

- **Skinner, Yu, et al. (2014)** — "Benchmarking all-atom simulations using hydrogen exchange." *PNAS*. The original HX-as-MD-benchmark protocol. Showed Shaw's all-atom NuG2b simulations had overly collapsed, over-H-bonded denatured states.
- **Robustelli, Piana, Shaw (2018)** — "Developing a molecular dynamics force field for both folded and disordered protein states." *PNAS*. Shaw's response: a99SB-disp force field that improves disordered protein accuracy.
- **Sosnick HDX Denaturation** — Additional Sosnick lab reference on HDX and denaturation.

## Benchmark Systems

The active benchmark follows the Skinner-adjacent DESRES systems that are
represented in this repository:

| Protein | Role in benchmark | Simulation set | HDX / validation status |
|---------|-------------------|----------------|-------------------------|
| lambda_d14a | Lambda repressor mutant extension | Lindorff-Larsen / Shaw lambda trajectories | Partial HDX/structural benchmark; verify residue-level HDX data before treating as a full HDX case |
| lambda_repressor | Lambda repressor extension | Lindorff-Larsen / Shaw lambda trajectories | Partial HDX/structural benchmark; verify residue-level HDX data before treating as a full HDX case |
| NuG2b | Primary Skinner HDX anchor | Lindorff-Larsen 2011 / Piana 2013 | Full Skinner HDX benchmark target with dG_HX and m-values |
| ubiquitin | Skinner-adjacent extension | Lindorff-Larsen 2011 / Piana 2013 / Robustelli 2018 | HDX/structural validation target; experimental conditions need per-dataset alignment |

## Approach

1. Replicate HX protection factors from DESMOND trajectories using protocol from Skinner et al. 2014
2. Run equivalent Upside 2.0 simulations on same proteins
3. Compare: H-bond stability, cooperativity, m-values, Rg
4. Assess whether Upside's coarse-grained force field captures the cooperative unfolding seen in experiment

## Upside Replica Exchange

Download available PDB reference structures:

```bash
python scripts/download_pdbs.py
```

The downloader currently prepares canonical PDB references for NuG2b and
ubiquitin. It removes the five-residue `HHHAM` expression tag from 1MI0 and
applies the three NuG2-to-NuG2b substitutions at positions 37, 46, and 47.
Lambda repressor and lambda D14A inputs are already present under
`simulations/<protein>/inputs`; add curated native PDB references for those
systems before using native-structure-dependent observables as publication
metrics.

The active simulations are stored under `simulations/<protein>/outputs/remd`
for `lambda_d14a`, `lambda_repressor`, `nug2b`, and `ubiquitin`. They use 16
replicas spanning reduced temperatures 0.80 to 0.96 with the quadratic spacing
used by the Upside example.

Use `notebooks/shaw_style_upside_analysis.ipynb` to discover the replica
outputs, compute structural observables, and regenerate the plots in
`docs/figures`.
