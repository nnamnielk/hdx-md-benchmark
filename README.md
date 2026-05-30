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

## Proteins with HDX Data (from DESRES simulations)

| Protein | Simulation Set | HDX Data | Source |
|---------|---------------|----------|--------|
| NuG2b | Lindorff-Larsen 2011 / Piana 2013 | dG_HX, m-values | Skinner et al. 2014 |
| HEWL | Robustelli 2018 | dG_HX ~12.2 kcal/mol | Radford et al. 1992 |
| NTL9 | Lindorff-Larsen 2011 | dG_HX ~4.7 kcal/mol | Kuhlman & Raleigh 1998 |
| Ubiquitin | Lindorff-Larsen 2011 / Robustelli 2018 | Yes | Sosnick Lab |

## Approach

1. Replicate HX protection factors from DESMOND trajectories using protocol from Skinner et al. 2014
2. Run equivalent Upside 2.0 simulations on same proteins
3. Compare: H-bond stability, cooperativity, m-values, Rg
4. Assess whether Upside's coarse-grained force field captures the cooperative unfolding seen in experiment
