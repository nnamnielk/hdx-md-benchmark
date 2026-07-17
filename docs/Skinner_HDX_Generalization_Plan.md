# Skinner-HDX Generalization Plan

## Objective

Use the Skinner, Yu, Gichana, Baxa, Hinshaw, Freed, and Sosnick 2014 HDX benchmark as the anchor, then generalize the same logic to Upside and to additional proteins/trajectories where the experimental data are strong enough.

The benchmark is not simply "does the simulation fold?" It asks whether the simulated ensemble has the same rare opening events inferred from denaturant-dependent hydrogen exchange:

- Are the most protected amides opened only by global unfolding?
- Do their Delta G_HX values match the global stability scale?
- Do their m_HX values match m_global, indicating exposure of the full denaturant-sensitive surface?
- When one protected H-bond is broken in simulation, are most/all other native H-bonds also broken, or does the model populate partially folded, overprotected states?

## What The Supplied Papers Establish

### Skinner et al. 2014

This is the reference implementation.

- Primary system: NuG2b, a fast-folding protein G variant.
- Simulation data: four DESRES/DESMOND all-atom NuG2b trajectories from Lindorff-Larsen et al. 2011, run near the simulated melting point at about 350 K.
- Experimental data: denaturant-dependent D-to-H exchange, CD equilibrium denaturation, kinetic chevrons, SAXS, and peptide CD.
- HX experiment: predeuterated NuG2b, exchange in 90% H2O / 10% D2O at 313 K, GdmCl series from 0 to 1.4 M, pH readings approximately 6.5-7.5 depending on denaturant.
- Key experimental conclusion: 12 protected sites exchange with Delta G_HX about 8.1 kcal/mol and m_HX about 1.31 kcal/mol/M, matching global unfolding. This implies cooperative opening into an expanded, mostly unstructured DSE.
- Key simulation conclusion: the DESRES DSE is too compact and too H-bonded. When one of the globally exchanging H-bonds breaks, many other H-bonds remain intact. Simulated m_HX from exposure is much smaller than m_global.

Skinner also includes supporting analyses of four lambda-repressor trajectories and seven ubiquitin trajectories. These are important extension targets, but the paper does not provide the same residue-resolved experimental HDX table for them that it provides for NuG2b.

### Robustelli, Piana, and Shaw 2018

This is not an HDX benchmark paper. It is a force-field benchmark for folded, disordered, and fast-folding proteins.

- Folded proteins include ubiquitin, GB3, HEWL, and BPTI.
- Fast-folding/peptide systems include villin, Trp-cage, GTT, CLN025, and (AAQAA)3.
- Disordered proteins are evaluated mostly against NMR, SAXS/Rg, chemical shifts, RDCs, scalar couplings, order parameters, and melting curves.
- Most production simulations are at 300 K; fast-folding peptides/proteins use simulated tempering from 278 to 400 K.
- These trajectories/data are useful as candidate systems and as context for all-atom force-field behavior, but they do not by themselves supply Skinner-style denaturant-dependent HDX observables.

Implication: HEWL or other Robustelli systems should not be treated as full Skinner-HDX benchmark cases unless we separately obtain residue-level HDX Delta G_HX/m_HX data and suitable folding/unfolding or DSE-sampling trajectories.

### Peng et al. / Sosnick HDX-Upside paper

This is the closest existing Upside analogue of the Skinner logic.

- Systems: EHEE_rd2_0005, HEEH_rd4_0097, mammalian ubiquitin, and ubiquitin L50E.
- Experimental HDX conditions differ by protein:
  - EHEE_rd2_0005: 298 K, pDread 7.1.
  - HEEH_rd4_0097: 278 K, pDread 4.6.
  - Ubiquitin: 273 K, pDread 7.6.
  - L50E: 277 K, pDread 7.5.
- Denaturant treatment is urea-based, not GdmCl-based as in Skinner's NuG2b experiments.
- Upside simulations are compared at temperatures chosen to match experimental stability, not always at the literal experimental temperature.
- Denaturant is modeled by reweighting conformations according to the number of protected backbone NH groups.
- For ubiquitin and L50E, reversible folding was hard to sample, so HDX was computed from native-start trajectories up to unfolding/misfolding cutoff points.

Implication: this paper gives us a practical Upside implementation pattern, but it also shows where rigor depends on sampling quality and temperature/stability matching.

## Are Experimental Condition Differences A Problem?

They are manageable but must be recorded per protein. They become a problem only if we mix raw quantities across incompatible conditions.

Important differences:

- Skinner NuG2b uses GdmCl-dependent HDX at 313 K, while the Upside HDX paper uses urea-dependent HDX at multiple temperatures/pD values.
- Skinner's DESRES simulations are at about 350 K, not 313 K. Skinner handles this by referencing site stabilities to the global stability rather than pretending the conditions are identical.
- Robustelli 2018 mainly compares native-state simulations to NMR/SAXS/melting observables, not denaturant-dependent HDX.
- Ubiquitin all-atom folding simulations from Piana et al. 2013 need their exact trajectory conditions confirmed from the trajectory metadata/original paper, whereas ubiquitin HDX in the Upside paper is near 273 K.

Recommended rule: benchmark each protein against its own experimental Delta G_HX, m_HX, pH/pD, temperature, denaturant, and global stability. Cross-protein comparisons should use normalized/cooperativity metrics such as Delta G_HX - Delta G_global, m_HX / m_global, and H-bond breakage cooperativity, not raw m-values across urea and GdmCl.

## Benchmark Definition

For each protein/trajectory pair, the generalized Skinner-HDX benchmark should produce:

1. Protein manifest
   - Sequence and variant.
   - Native structure/PDB or reference model.
   - Native H-bond list and residue mapping.
   - Simulation source, force field/model, temperature schedule, trajectory count, frame spacing, and total sampling.
   - Experimental source, denaturant, temperature, pH/pD, Delta G_global, m_global, and residue-level Delta G_HX/m_HX if available.

2. Simulation ensemble parsing
   - Classify frames into NSE, DSE, and optionally intermediate/subglobal states.
   - For Skinner reproduction, use NuG2b criteria: backbone RMSD < 4 A, TM-score > 0.6, and native H-bonds > 20 for NSE; all other frames are DSE.
   - For other proteins, define analogous thresholds from population distributions and document them.

3. Protection/opening calculation
   - DESRES/all-atom: reproduce Skinner's H-bond protection definition: O-H distance < 3.5 A and H-N-O angle > 120 deg, or this state accessed within 0.5 ns.
   - Upside: use the Sosnick/Peng adaptation: an NH is protected if it is H-bonded or sufficiently buried; exchange-competent if H-bond broken and exposed. The exact score threshold must be fixed from the Upside implementation.
   - Compute Delta G_HX-like protection from protected/open populations, referenced to global stability where needed.

4. Denaturant/cooperativity calculation
   - If experimental m_HX exists, compare m_HX and m_HX / m_global.
   - For all-atom trajectories, estimate exposure of open states using ASA or equivalent.
   - For Upside, model denaturant by reweighting conformations using protected NH count or an agreed exposure proxy, calibrated per protein to m_global.
   - Classify openings as local, subglobal, or global.

5. Outputs
   - Reproduction of Skinner NuG2b figures/tables as a validation gate.
   - Per-residue Delta G_HX and m_HX comparison plots.
   - Cooperativity plots: number/fraction of remaining H-bonds when each protected site is open.
   - DSE structure summaries: Rg, native/non-native H-bond counts, contact maps, secondary structure frequencies.
   - A per-protein confidence label: full HDX benchmark, partial HDX benchmark, or non-HDX auxiliary validation.

## Implementation Plan

### Phase 1: Reproduce Skinner NuG2b exactly

- Obtain the exact DESRES NuG2b trajectories used by Skinner, ideally from the DESRES/Science 2011 dataset or lab-held copies.
- Build a trajectory reader and manifest.
- Implement Skinner's H-bond, contact, RMSD, TM-score, DSSP, and Rg analyses.
- Reproduce NSE/DSE parsing and Figure 1 / Figure S1-like summaries.
- Reproduce Table S1 comparison logic: Delta G_HX relative to global stability.
- Reproduce the core discrepancy: simulated protected H-bonds exceed global stability and open through partially folded states.

### Phase 2: Generalize the analysis API

- Make protein-specific configuration explicit: native structure, residue mapping, H-bond list, thresholds, experimental condition metadata.
- Separate "protection definition" from "trajectory source" so all-atom DESRES and Upside can use different low-level observables but feed the same benchmark metrics.
- Add a report generator that emits the same metrics for any protein.

### Phase 3: Add Skinner-adjacent DESRES systems

- Add ubiquitin DESRES trajectories from Piana/Lindorff-Larsen/Shaw 2013.
- Add lambda-repressor trajectories from Lindorff-Larsen/Shaw 2011 if available.
- For each, determine whether residue-level denaturant HDX data exist. If not, score structural/cooperativity metrics and global experimental consistency, but mark as partial HDX benchmark.

### Phase 4: Add Upside comparisons

- Run Upside on NuG2b first, because Skinner provides the strongest experimental anchor.
- Match the experimental stability scale, following the Peng/Sosnick precedent, rather than forcing the literal experimental temperature if that breaks comparability.
- Compute Upside protection/opening probabilities and denaturant reweighting.
- Compare Upside vs DESRES vs experiment on normalized metrics.
- Then add ubiquitin/L50E and the designed proteins from the Upside HDX paper as additional Upside-native validation cases.

### Phase 5: Expand only when data qualify

Candidate proteins such as HEWL and NTL9 need a data audit before inclusion as full Skinner-HDX cases. They require:

- residue-level HDX or at least enough protected-site Delta G_HX/m_HX data,
- global stability and m_global under compatible conditions,
- a simulation ensemble that samples native, open, and DSE/subglobal states sufficiently,
- sequence/variant consistency between experiment and simulation.

## Actual Ambiguities To Confirm

1. Exact DESRES data access
   - Do we have lab copies of the NuG2b, lambda, and ubiquitin trajectories, or should we request them from DESRES/original authors?
   - Are the files full all-atom, stripped, or downsampled? The 0.5 ns H-bond memory criterion requires frame spacing fine enough to evaluate short closures/openings.

2. Superset membership
   - The natural Skinner superset is NuG2b plus the lambda-repressor and ubiquitin trajectories already discussed in Skinner.
   - HEWL and NTL9 should not be promoted to full benchmark status until their residue-level HDX data and matching trajectories are verified.

3. Upside temperature policy
   - Use literal experimental temperature, simulated Tm, or the Peng/Sosnick stability-matched temperature? The most defensible default is stability-matched temperature, with experimental-temperature results as sensitivity checks.

4. Upside protection definition
   - We need to confirm the exact Upside FF/version and protection thresholds used for H-bond score and burial level.

5. Denaturant model
   - Skinner uses GdmCl experiments; the Upside HDX paper uses urea. We should not compare raw m-values across denaturants without per-protein calibration.

6. Sampling standard
   - Full equilibrium folding/unfolding is ideal.
   - Native-to-DSE truncation, as used for ubiquitin/L50E in the Upside HDX paper, is acceptable but should be labeled lower confidence.

7. Experimental refitting
   - Decide whether to use published Delta G_HX/m_HX tables directly or refit raw exchange rates with our own intrinsic-rate corrections. Published values are sufficient for a first benchmark; refitting is a later rigor upgrade.
