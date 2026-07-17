#!/usr/bin/env python3
"""Run the Upside HDX example workflow for the local NuG2b REMD data."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


N_REPLICAS = 16
KCAL_PER_MOL_K = 0.00198720425864083
SKINNER_EXPERIMENT_GLOBAL_KEY = "kinetics"
SKINNER_FIG3_EXPERIMENT_GLOBAL_DELTA_G = 7.70
SKINNER_NUG2B_GLOBALS = {
    "CD": {"delta_g": 7.48, "se": 0.23},
    "kinetics": {"delta_g": 7.78, "se": 0.23},
    "m_hx": {"value": 1.29, "se": 0.10},
}
SKINNER_STATE_RMSD_CUTOFF_ANGSTROM = 4.0
SKINNER_STATE_Q_NATIVE_CUTOFF = 0.60
SKINNER_STATE_NATIVE_HBOND_CUTOFF = 20.0
SKINNER_NUG2B_HX_TABLE = [
    (3, "Y", 8.08, 0.11),
    (4, "K", 8.54, 0.11),
    (5, "L", 8.31, 0.11),
    (6, "V", 8.09, 0.11),
    (7, "I", 7.84, 0.11),
    (16, "Y", 7.85, 0.14),
    (18, "T", 7.85, 0.17),
    (26, "A", 8.45, 0.15),
    (27, "E", 7.74, 0.12),
    (30, "F", 7.93, 0.12),
    (31, "K", 7.94, 0.14),
    (51, "T", 8.30, 0.12),
]
N_TERMINAL_TAGS = {
    "nug2b": ("HIS", "HIS", "HIS", "ALA", "MET"),
}
ONE_LETTER_AA = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute NuG2b Delta G_HX using Upside's HDX protection-state example logic."
    )
    parser.add_argument("--protein", default="nug2b")
    parser.add_argument("--sim-id", default="remd")
    parser.add_argument(
        "--run-dir",
        type=Path,
        help=(
            "Directory containing <protein>.run.<replica>.up files. "
            "Defaults to simulations/<protein>/outputs/<sim-id>."
        ),
    )
    parser.add_argument("--start-frame", type=int, default=100)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--n-replicas", type=int, default=N_REPLICAS)
    parser.add_argument("--target-temperature", type=float, default=0.80)
    parser.add_argument("--experimental-temperature-k", type=float, default=313.0)
    parser.add_argument("--denaturant-sensitivity", type=float, default=0.04)
    parser.add_argument("--denaturant-max", type=float, default=15.0)
    parser.add_argument("--denaturant-bins", type=int, default=150)
    parser.add_argument(
        "--m-value-fit-min",
        type=float,
        default=0.0,
        help="Minimum denaturant concentration used for the reported finite-difference m-value.",
    )
    parser.add_argument(
        "--m-value-fit-max",
        type=float,
        default=1.4,
        help="Maximum denaturant concentration used for the reported finite-difference m-value.",
    )
    parser.add_argument("--force-field", default="ff_2.1")
    parser.add_argument(
        "--upside-home",
        default=os.environ.get("UPSIDE_HOME", "/home/okleinmann/projects/upside2-md/dynalab"),
        help="Linux Upside checkout with py/get_protection_state.py and obj/libupside.so.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable to use when invoking Upside helper scripts.",
    )
    return parser.parse_args()


def output_group_names(handle) -> list[str]:
    previous = sorted(
        [name for name in handle.keys() if name.startswith("output_previous_")],
        key=lambda name: int(name.rsplit("_", 1)[1]),
    )
    return previous + (["output"] if "output" in handle else [])


def read_output_series(path: Path, dataset: str, stride: int) -> np.ndarray:
    import h5py

    chunks = []
    with h5py.File(path, "r", locking=False) as handle:
        for group_index, group_name in enumerate(output_group_names(handle)):
            group = handle[group_name]
            if dataset not in group:
                continue
            raw = np.asarray(group[dataset]).squeeze()
            if group_index > 0 and raw.shape[0] > 1:
                raw = raw[1:]
            chunks.append(raw[::stride])
    if not chunks:
        raise KeyError(f"{dataset!r} not found in {path}")
    return np.concatenate(chunks, axis=0)


def read_input_sequence_and_positions(path: Path) -> tuple[list[str], np.ndarray]:
    import h5py

    with h5py.File(path, "r", locking=False) as handle:
        sequence = [
            item.decode("utf-8") if isinstance(item, bytes) else str(item)
            for item in handle["input/sequence"][:]
        ]
        positions = np.asarray(handle["input/pos"]).squeeze().astype(float)
    return sequence, positions


def read_positions(path: Path, stride: int) -> np.ndarray:
    positions = read_output_series(path, "pos", stride).astype(float)
    if positions.ndim == 4 and positions.shape[1] == 1:
        positions = positions[:, 0]
    return positions


def read_hbond_pair_indices(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    import h5py

    with h5py.File(path, "r", locking=False) as handle:
        group = handle["input/potential/protein_hbond"]
        donor_res = np.asarray(group["id1"], dtype=int)
        acceptor_res = np.asarray(group["id2"], dtype=int)
        idx1 = np.asarray(group["index1"], dtype=int)
        idx2 = np.asarray(group["index2"], dtype=int)
    return donor_res, acceptor_res, idx1, idx2


def ca_positions(positions: np.ndarray) -> np.ndarray:
    return positions[:, 1::3, :] if positions.ndim == 3 else positions[1::3, :]


def kabsch_rmsd_series(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    ref = reference - reference.mean(axis=0, keepdims=True)
    rmsd = np.empty(mobile.shape[0], dtype=float)
    for index, coords in enumerate(mobile):
        mob = coords - coords.mean(axis=0, keepdims=True)
        cov = mob.T @ ref
        v, _, wt = np.linalg.svd(cov)
        if np.linalg.det(v @ wt) < 0:
            v[:, -1] *= -1
        aligned = mob @ (v @ wt)
        rmsd[index] = np.sqrt(((aligned - ref) ** 2).sum(axis=1).mean())
    return rmsd


def native_contact_pairs_from_ca(
    reference_ca: np.ndarray,
    min_separation: int = 3,
    cutoff_angstrom: float = 8.0,
) -> np.ndarray:
    pairs = []
    for i in range(len(reference_ca)):
        for j in range(i + min_separation, len(reference_ca)):
            if np.linalg.norm(reference_ca[i] - reference_ca[j]) < cutoff_angstrom:
                pairs.append((i, j))
    return np.asarray(pairs, dtype=int)


def contact_fraction_series(
    ca: np.ndarray,
    pairs: np.ndarray,
    cutoff_angstrom: float = 8.0,
) -> np.ndarray:
    if len(pairs) == 0:
        return np.full(ca.shape[0], np.nan)
    delta = ca[:, pairs[:, 0], :] - ca[:, pairs[:, 1], :]
    distances = np.sqrt((delta**2).sum(axis=2))
    return (distances < cutoff_angstrom).mean(axis=1)


def hbond_pair_scores(path: Path, stride: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hbond = read_output_series(path, "hbond", stride).astype(float)
    donor_res, acceptor_res, idx1, idx2 = read_hbond_pair_indices(path)
    usable = (idx1 < hbond.shape[1]) & (idx2 < hbond.shape[1])
    pair_scores = 0.5 * (hbond[:, idx1[usable]] + hbond[:, idx2[usable]])
    return pair_scores, donor_res[usable], acceptor_res[usable]


def state_mask_for_run(protein: str, run_file: Path, stride: int) -> np.ndarray:
    sequence, reference = read_input_sequence_and_positions(run_file)
    positions = read_positions(run_file, stride)
    pair_scores, donor_res, acceptor_res = hbond_pair_scores(run_file, stride)
    tag_len = n_terminal_tag_length(protein, sequence)
    if tag_len:
        positions = positions[:, tag_len * 3 :, :]
        reference = reference[tag_len * 3 :, :]
        keep_pairs = (donor_res >= tag_len) & (acceptor_res >= tag_len)
        pair_scores = pair_scores[:, keep_pairs]
    ca = ca_positions(positions)
    reference_ca = reference[1::3, :]
    rmsd = kabsch_rmsd_series(ca, reference_ca)
    q_native = contact_fraction_series(ca, native_contact_pairs_from_ca(reference_ca))
    hbond_sites = (pair_scores > 0.5).sum(axis=1)
    n = min(len(rmsd), len(q_native), len(hbond_sites))
    return (
        (rmsd[:n] < SKINNER_STATE_RMSD_CUTOFF_ANGSTROM)
        & (q_native[:n] > SKINNER_STATE_Q_NATIVE_CUTOFF)
        & (hbond_sites[:n] > SKINNER_STATE_NATIVE_HBOND_CUTOFF)
    )


def load_state_masks(args: argparse.Namespace, run_files: list[Path]) -> np.ndarray:
    masks = []
    for replica, run_file in enumerate(run_files):
        print(f"state masks replica {replica}")
        masks.append(state_mask_for_run(args.protein, run_file, args.stride))
    min_frames = min(len(mask) for mask in masks)
    return np.asarray([mask[:min_frames] for mask in masks], dtype=bool)


def n_terminal_tag_length(protein: str, sequence: list[str]) -> int:
    tag = N_TERMINAL_TAGS.get(protein)
    if tag is None or tuple(sequence[: len(tag)]) != tag:
        return 0
    return len(tag)


def write_fasta(path: Path, sequence: list[str], source: Path) -> None:
    one_letter = "".join(ONE_LETTER_AA[item] for item in sequence)
    with path.open("w", encoding="ascii") as handle:
        handle.write(f"> Created from {source}\n")
        for start in range(0, len(one_letter), 80):
            handle.write(one_letter[start : start + 80] + "\n")


def ensure_hdx_top(args: argparse.Namespace, root: Path, input_dir: Path) -> Path:
    upside_home = Path(args.upside_home).expanduser().resolve()
    upside_py = upside_home / "py"
    hdx_dir = input_dir / "hdx"
    hdx_dir.mkdir(parents=True, exist_ok=True)
    hdx_top = hdx_dir / f"{args.protein}-HDX.up"
    base_top = input_dir / f"{args.protein}.up"
    base_sequence, base_positions = read_input_sequence_and_positions(base_top)
    if hdx_top.exists():
        hdx_sequence, _ = read_input_sequence_and_positions(hdx_top)
        if len(hdx_sequence) == len(base_sequence):
            return hdx_top
        hdx_top.unlink()

    prefix = hdx_dir / args.protein
    write_fasta(prefix.with_suffix(".fasta"), base_sequence, base_top)
    np.save(prefix.with_suffix(".initial.npy"), base_positions)
    chain_breaks = prefix.with_suffix(".chain_breaks")
    if chain_breaks.exists():
        chain_breaks.unlink()

    if hdx_top.exists():
        return hdx_top

    sys.path.insert(0, str(upside_py))
    import run_upside as ru

    parameter_dir = upside_home / "parameters"
    common_dir = parameter_dir / "common"
    force_field_dir = parameter_dir / args.force_field
    config_options = dict(
        rama_library=str(common_dir / "rama.dat"),
        rama_sheet_mix_energy=str(force_field_dir / "sheet"),
        reference_state_rama=str(common_dir / "rama_reference.pkl"),
        hbond_energy=str(force_field_dir / "hbond.h5"),
        rotamer_placement=str(force_field_dir / "sidechain.h5"),
        dynamic_rotamer_1body=True,
        rotamer_interaction=str(force_field_dir / "sidechain.h5"),
        environment_potential=str(force_field_dir / "environment.h5"),
        bb_environment_potential=str(force_field_dir / "bb_env.dat"),
        initial_structure=str(prefix.with_suffix(".initial.npy")),
        use_heavy_atom_coverage=True,
    )
    if chain_breaks.exists():
        config_options["chain_break_from_file"] = str(chain_breaks)
    print(ru.upside_config(str(prefix.with_suffix(".fasta")), str(hdx_top), **config_options))
    return hdx_top


def run_protection_state(
    args: argparse.Namespace,
    hdx_top: Path,
    run_file: Path,
    output_npy: Path,
    residue_file: Path,
) -> None:
    if output_npy.exists() and residue_file.exists():
        return
    upside_home = Path(args.upside_home).expanduser().resolve()
    env = os.environ.copy()
    env["UPSIDE_HOME"] = str(upside_home)
    env["HDF5_USE_FILE_LOCKING"] = "FALSE"
    subprocess.run(
        [
            args.python,
            str(upside_home / "py" / "get_protection_state.py"),
            str(hdx_top),
            str(run_file),
            str(output_npy),
            "--residue",
            str(residue_file),
            "--stride",
            str(args.stride),
        ],
        check=True,
        env=env,
    )


def load_arrays(args: argparse.Namespace, run_files: list[Path], result_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    potentials = []
    temperatures = []
    protection = []
    for replica, run_file in enumerate(run_files):
        ps_file = result_dir / f"{args.protein}_{args.sim_id}_stride{args.stride}_{replica}_PS.npy"
        pot = read_output_series(run_file, "potential", args.stride).reshape(-1)
        temp = read_output_series(run_file, "temperature", args.stride).reshape(-1)
        ps = np.load(ps_file)
        n = min(len(pot), len(temp), ps.shape[0])
        potentials.append(pot[:n])
        temperatures.append(temp[:n])
        protection.append(ps[:n])
    return np.asarray(potentials), np.asarray(temperatures), np.asarray(protection)


def discover_run_files(args: argparse.Namespace, run_dir: Path) -> list[Path]:
    run_files = []
    for replica in range(args.n_replicas):
        name = f"{args.protein}.run.{replica}.up"
        matches = sorted(run_dir.rglob(name))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one {name} below {run_dir}, found {len(matches)}: {matches}"
            )
        run_files.append(matches[0])
    return run_files


def trim_to_common_frames(
    potentials: np.ndarray,
    temperatures: np.ndarray,
    protection: np.ndarray,
    state_masks: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    frame_count = min(potentials.shape[1], temperatures.shape[1], protection.shape[1], state_masks.shape[1])
    return (
        potentials[:, :frame_count],
        temperatures[:, :frame_count],
        protection[:, :frame_count],
        state_masks[:, :frame_count],
    )


def mbar_weights(potentials: np.ndarray, temperatures: np.ndarray, start_frame: int, target_temperature: float) -> np.ndarray:
    from pymbar import MBAR

    k_b = 1.0
    replica_temperatures = np.median(temperatures[:, start_frame:], axis=1)
    beta = 1.0 / (k_b * replica_temperatures)
    c_e0 = potentials[:, start_frame:]
    frame_count = c_e0.shape[1]
    counts = np.full(c_e0.shape[0], frame_count, dtype=np.int32)
    reduced = np.zeros((c_e0.shape[0], c_e0.shape[0], frame_count), dtype=np.float64)
    for source in range(c_e0.shape[0]):
        for target in range(c_e0.shape[0]):
            reduced[source, target] = beta[target] * c_e0[source]
    mbar = MBAR(reduced, counts, verbose=False)
    target_u = (c_e0 / (target_temperature * k_b)).reshape(-1)
    try:
        log_weights = mbar._computeUnnormalizedLogWeights(target_u)
    except AttributeError:
        log_weights = mbar._compute_log_weights(target_u)
    weights = np.exp(log_weights - np.max(log_weights))
    return weights / weights.sum()


def build_mbar_context(
    potentials: np.ndarray,
    temperatures: np.ndarray,
    start_frame: int,
):
    from pymbar import MBAR

    k_b = 1.0
    replica_temperatures = np.median(temperatures[:, start_frame:], axis=1)
    beta = 1.0 / (k_b * replica_temperatures)
    c_e0 = potentials[:, start_frame:]
    frame_count = c_e0.shape[1]
    counts = np.full(c_e0.shape[0], frame_count, dtype=np.int32)
    reduced = np.zeros((c_e0.shape[0], c_e0.shape[0], frame_count), dtype=np.float64)
    for source in range(c_e0.shape[0]):
        for target in range(c_e0.shape[0]):
            reduced[source, target] = beta[target] * c_e0[source]
    return MBAR(reduced, counts, verbose=False), c_e0, replica_temperatures


def mbar_weights_from_context(mbar, c_e0: np.ndarray, target_temperature: float) -> np.ndarray:
    target_u = (c_e0 / target_temperature).reshape(-1)
    try:
        log_weights = mbar._computeUnnormalizedLogWeights(target_u)
    except AttributeError:
        log_weights = mbar._compute_log_weights(target_u)
    weights = np.exp(log_weights - np.max(log_weights))
    return weights / weights.sum()


def calculate_hdx(
    args: argparse.Namespace,
    potentials: np.ndarray,
    temperatures: np.ndarray,
    protection: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    start = args.start_frame
    ps = protection[:, start:, :]
    weights = mbar_weights(potentials, temperatures, start, args.target_temperature)
    protected_prob = np.zeros(ps.shape[2], dtype=float)
    flat_ps = ps.reshape(-1, ps.shape[2])
    for residue in range(ps.shape[2]):
        protected_prob[residue] = np.average(flat_ps[:, residue], weights=weights)
    open_prob = 1.0 - protected_prob
    delta_g = KCAL_PER_MOL_K * args.experimental_temperature_k * np.log(
        np.clip(protected_prob, 1e-12, 1.0) / np.clip(open_prob, 1e-12, 1.0)
    )

    pf_frame = ps.sum(axis=2).reshape(-1)
    den = np.linspace(0.0, args.denaturant_max, args.denaturant_bins + 1)
    dghx_den = []
    for denaturant in den:
        den_weights = np.exp((-args.denaturant_sensitivity * pf_frame) * denaturant / args.target_temperature)
        den_weights *= weights
        den_weights /= den_weights.sum()
        prob = den_weights @ flat_ps
        dghx_den.append(
            KCAL_PER_MOL_K
            * args.experimental_temperature_k
            * np.log(np.clip(prob, 1e-12, 1.0) / np.clip(1.0 - prob, 1e-12, 1.0))
        )
    dghx_den = np.asarray(dghx_den)
    m_values = np.diff(dghx_den, axis=0) / (den[0] - den[1])
    return delta_g, dghx_den, m_values


def m_value_window(args: argparse.Namespace) -> np.ndarray:
    if args.m_value_fit_max <= args.m_value_fit_min:
        raise ValueError("--m-value-fit-max must be larger than --m-value-fit-min")
    den = np.linspace(0.0, args.denaturant_max, args.denaturant_bins + 1)
    window = (den[:-1] >= args.m_value_fit_min) & (den[1:] <= args.m_value_fit_max)
    if not np.any(window):
        raise ValueError(
            "No denaturant bins fall inside the requested m-value fit window "
            f"[{args.m_value_fit_min}, {args.m_value_fit_max}]"
        )
    return window


def calculate_delta_g_temperature_sweep(
    args: argparse.Namespace,
    potentials: np.ndarray,
    temperatures: np.ndarray,
    protection: np.ndarray,
    state_masks: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    start = args.start_frame
    ps = protection[:, start:, :]
    nse = state_masks[:, start:].reshape(-1).astype(float)
    flat_ps = ps.reshape(-1, ps.shape[2])
    mbar, c_e0, target_temperatures = build_mbar_context(potentials, temperatures, start)
    delta_g_by_temperature = []
    global_delta_g_by_temperature = []
    nse_probability_by_temperature = []
    for target_temperature in target_temperatures:
        weights = mbar_weights_from_context(mbar, c_e0, float(target_temperature))
        protected_prob = weights @ flat_ps
        open_prob = 1.0 - protected_prob
        delta_g = KCAL_PER_MOL_K * args.experimental_temperature_k * np.log(
            np.clip(protected_prob, 1e-12, 1.0) / np.clip(open_prob, 1e-12, 1.0)
        )
        nse_probability = float(weights @ nse)
        dse_probability = 1.0 - nse_probability
        global_delta_g = KCAL_PER_MOL_K * args.experimental_temperature_k * np.log(
            np.clip(nse_probability, 1e-12, 1.0) / np.clip(dse_probability, 1e-12, 1.0)
        )
        delta_g_by_temperature.append(delta_g)
        global_delta_g_by_temperature.append(global_delta_g)
        nse_probability_by_temperature.append(nse_probability)
    return (
        target_temperatures,
        np.asarray(delta_g_by_temperature),
        np.asarray(global_delta_g_by_temperature),
        np.asarray(nse_probability_by_temperature),
        1.0 - np.asarray(nse_probability_by_temperature),
    )


def tagless_output_indices(protein: str, input_dir: Path, residue_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sequence, _ = read_input_sequence_and_positions(input_dir / f"{protein}.up")
    tag_len = n_terminal_tag_length(protein, sequence)
    residue_ids = residue_ids.astype(int)
    keep = residue_ids >= tag_len
    tagless_residues = residue_ids[keep] - tag_len + 1
    return keep, tagless_residues


def write_outputs(
    args: argparse.Namespace,
    root: Path,
    input_dir: Path,
    residue_ids: np.ndarray,
    delta_g: np.ndarray,
    dghx_den: np.ndarray,
    m_values: np.ndarray,
) -> None:
    data_dir = root / "docs" / "data"
    figure_dir = root / "docs" / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    keep, tagless_residues = tagless_output_indices(args.protein, input_dir, residue_ids)
    delta_g = delta_g[keep]
    m_values = m_values[:, keep]
    exp = {residue: (aa, dg, se) for residue, aa, dg, se in SKINNER_NUG2B_HX_TABLE}
    csv_path = data_dir / f"{args.protein}_upside_hdx.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "residue",
                "skinner_aa",
                "upside_delta_g_hx",
                "upside_m_value",
                "skinner_delta_g_hx",
                "skinner_se",
            ]
        )
        mean_m = np.nanmean(m_values[m_value_window(args)], axis=0)
        for index, residue in enumerate(tagless_residues.astype(int)):
            aa, skinner_dg, skinner_se = exp.get(int(residue), ("", np.nan, np.nan))
            writer.writerow([int(residue), aa, delta_g[index], mean_m[index], skinner_dg, skinner_se])

    import matplotlib.pyplot as plt

    residues = np.array([row[0] for row in SKINNER_NUG2B_HX_TABLE], dtype=int)
    exp_y = np.array([row[2] for row in SKINNER_NUG2B_HX_TABLE], dtype=float)
    exp_se = np.array([row[3] for row in SKINNER_NUG2B_HX_TABLE], dtype=float)
    sim_y = np.array([delta_g[np.where(tagless_residues == residue)[0][0]] for residue in residues])
    labels = [f"{row[0]}{row[1]}" for row in SKINNER_NUG2B_HX_TABLE]
    fig, ax = plt.subplots(figsize=(7.4, 3.7), constrained_layout=True)
    ax.errorbar(residues, exp_y, yerr=exp_se, fmt="s", color="black", capsize=2, label="Skinner experiment")
    ax.plot(residues, sim_y, "o-", color="#d95f02", label="Upside HDX PS + MBAR")
    ax.axhline(7.48, color="#0072b2", ls="--", lw=1.0, label="Skinner CD global")
    ax.axhline(7.78, color="#009e73", ls=":", lw=1.3, label="Skinner kinetic global")
    ax.set_xticks(residues)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_xlabel("NuG2b residue")
    ax.set_ylabel(r"$\Delta G_{HX}$ (kcal/mol)")
    ax.set_title("NuG2b HDX stability: Skinner experiment vs Upside HDX workflow")
    ax.grid(axis="y", ls=":", lw=0.4, alpha=0.45)
    ax.legend(frameon=False)
    fig.savefig(figure_dir / f"{args.protein}_upside_hdx_delta_g_overlay.png", bbox_inches="tight")
    plt.close(fig)

    np.save(data_dir / f"{args.protein}_upside_hdx_dghx_denaturant.npy", dghx_den)


def write_temperature_sweep_outputs(
    args: argparse.Namespace,
    root: Path,
    input_dir: Path,
    residue_ids: np.ndarray,
    target_temperatures: np.ndarray,
    delta_g_by_temperature: np.ndarray,
    global_delta_g_by_temperature: np.ndarray,
    nse_probability_by_temperature: np.ndarray,
    dse_probability_by_temperature: np.ndarray,
) -> None:
    data_dir = root / "docs" / "data"
    figure_dir = root / "docs" / "figures"
    per_temperature_dir = figure_dir / "upside_hdx_by_temperature"
    data_dir.mkdir(parents=True, exist_ok=True)
    per_temperature_dir.mkdir(parents=True, exist_ok=True)

    keep, tagless_residues = tagless_output_indices(args.protein, input_dir, residue_ids)
    delta_g_by_temperature = delta_g_by_temperature[:, keep]
    exp = {residue: (aa, dg, se) for residue, aa, dg, se in SKINNER_NUG2B_HX_TABLE}
    skinner_residues = np.array([row[0] for row in SKINNER_NUG2B_HX_TABLE], dtype=int)
    skinner_labels = [f"{row[0]}{row[1]}" for row in SKINNER_NUG2B_HX_TABLE]
    exp_y = np.array([row[2] for row in SKINNER_NUG2B_HX_TABLE], dtype=float)
    exp_se = np.array([row[3] for row in SKINNER_NUG2B_HX_TABLE], dtype=float)
    experiment_global_delta_g = SKINNER_FIG3_EXPERIMENT_GLOBAL_DELTA_G
    exp_relative = exp_y - experiment_global_delta_g
    exp_mean_relative = float(np.nanmean(exp_relative))
    skinner_indices = [np.where(tagless_residues == residue)[0][0] for residue in skinner_residues]
    sweep_values = delta_g_by_temperature[:, skinner_indices]
    sweep_relative = sweep_values - global_delta_g_by_temperature[:, None]
    sweep_mean_relative = np.nanmean(sweep_relative, axis=1)

    csv_path = data_dir / f"{args.protein}_upside_hdx_by_temperature.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "target_reduced_temperature",
                "residue",
                "skinner_aa",
                "upside_delta_g_hx",
                "upside_delta_g_hx_minus_global",
                "upside_delta_g_global",
                "skinner_delta_g_hx",
                "skinner_delta_g_hx_minus_global",
                "skinner_delta_g_global",
                "skinner_se",
                "at_probability_cap",
            ]
        )
        cap = KCAL_PER_MOL_K * args.experimental_temperature_k * np.log(1.0 / 1e-12)
        for temp_index, target_temperature in enumerate(target_temperatures):
            for index, residue in enumerate(tagless_residues.astype(int)):
                aa, skinner_dg, skinner_se = exp.get(int(residue), ("", np.nan, np.nan))
                writer.writerow(
                    [
                        float(target_temperature),
                        int(residue),
                        aa,
                        delta_g_by_temperature[temp_index, index],
                        delta_g_by_temperature[temp_index, index] - global_delta_g_by_temperature[temp_index],
                        global_delta_g_by_temperature[temp_index],
                        skinner_dg,
                        skinner_dg - experiment_global_delta_g if np.isfinite(skinner_dg) else np.nan,
                        experiment_global_delta_g,
                        skinner_se,
                        bool(delta_g_by_temperature[temp_index, index] >= cap - 1e-6),
                    ]
                )

    summary_path = data_dir / f"{args.protein}_upside_hdx_temperature_sweep_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "target_reduced_temperature",
                "upside_delta_g_global",
                "upside_p_nse",
                "upside_p_dse",
                "upside_mean_delta_g_hx",
                "upside_mean_delta_g_hx_minus_global",
                "skinner_global_reference",
                "skinner_delta_g_global",
                "skinner_mean_delta_g_hx",
                "skinner_mean_delta_g_hx_minus_global",
            ]
        )
        for index, target_temperature in enumerate(target_temperatures):
            writer.writerow(
                [
                    float(target_temperature),
                    global_delta_g_by_temperature[index],
                    nse_probability_by_temperature[index],
                    dse_probability_by_temperature[index],
                    float(np.nanmean(sweep_values[index])),
                    sweep_mean_relative[index],
                    f"{SKINNER_EXPERIMENT_GLOBAL_KEY}_fig3_rounded",
                    experiment_global_delta_g,
                    float(np.nanmean(exp_y)),
                    exp_mean_relative,
                ]
            )

    import matplotlib.pyplot as plt

    def draw_overlay(
        ax,
        target_temperature: float,
        sim_relative: np.ndarray,
        sim_mean_relative: float,
        show_ylabel: bool = False,
    ) -> None:
        ax.errorbar(
            skinner_residues,
            exp_relative,
            yerr=exp_se,
            fmt="s",
            ms=3.0,
            color="black",
            ecolor="0.35",
            elinewidth=0.7,
            capsize=1.5,
        )
        ax.plot(skinner_residues, sim_relative, "o", ms=3.1, color="red")
        ax.axhline(0.0, color="0.25", lw=0.7)
        ax.axhline(exp_mean_relative, color="black", ls=":", lw=0.9)
        ax.axhline(sim_mean_relative, color="red", ls=":", lw=0.9)
        ax.text(
            0.985,
            exp_mean_relative,
            f"{exp_mean_relative:.1f}",
            color="black",
            ha="right",
            va="bottom",
            fontsize=7,
            transform=ax.get_yaxis_transform(),
        )
        ax.text(
            0.985,
            sim_mean_relative,
            f"{sim_mean_relative:.1f}",
            color="red",
            ha="right",
            va="bottom",
            fontsize=7,
            transform=ax.get_yaxis_transform(),
        )
        ax.set_title(f"T={target_temperature:.4f}", fontsize=8)
        y_min = min(-0.4, float(np.nanmin([np.nanmin(exp_relative), np.nanmin(sim_relative), exp_mean_relative, sim_mean_relative])) - 0.4)
        y_max = max(2.8, float(np.nanmax([np.nanmax(exp_relative), np.nanmax(sim_relative), exp_mean_relative, sim_mean_relative])) + 0.4)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks(skinner_residues)
        ax.set_xticklabels(skinner_labels, rotation=60, ha="right", fontsize=6)
        if show_ylabel:
            ax.set_ylabel(r"$\Delta G_i - \Delta G_{global}$ (kcal/mol)")
        ax.grid(axis="y", ls=":", lw=0.3, alpha=0.45)

    ncols = 4
    nrows = int(np.ceil(len(target_temperatures) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.5, 7.2), constrained_layout=True)
    axes_flat = axes.ravel()
    for index, (ax, target_temperature, sim_relative) in enumerate(
        zip(axes_flat, target_temperatures, sweep_relative)
    ):
        draw_overlay(
            ax,
            float(target_temperature),
            sim_relative,
            float(sweep_mean_relative[index]),
            show_ylabel=(index % ncols == 0),
        )
    for ax in axes_flat[len(target_temperatures) :]:
        ax.axis("off")
    fig.suptitle(
        r"NuG2b HDX stabilities referenced to $\Delta G_{global}$ across the Upside REMD temperature ladder",
        y=1.02,
    )
    sweep_path = figure_dir / f"{args.protein}_upside_hdx_delta_g_temperature_sweep.png"
    fig.savefig(sweep_path, bbox_inches="tight")
    plt.close(fig)

    for index, (target_temperature, sim_relative) in enumerate(zip(target_temperatures, sweep_relative)):
        fig, ax = plt.subplots(figsize=(7.4, 3.7), constrained_layout=True)
        draw_overlay(
            ax,
            float(target_temperature),
            sim_relative,
            float(sweep_mean_relative[index]),
            show_ylabel=True,
        )
        ax.set_xlabel("Residue index")
        fig.savefig(
            per_temperature_dir / f"{args.protein}_upside_hdx_delta_g_T{float(target_temperature):.4f}.png",
            bbox_inches="tight",
        )
        plt.close(fig)

    print("wrote", csv_path)
    print("wrote", summary_path)
    print("wrote", sweep_path)
    print("wrote", per_temperature_dir)
    print(
        "Skinner experiment mean DeltaG_HX - DeltaG_global "
        f"({SKINNER_EXPERIMENT_GLOBAL_KEY}, Fig. 3 rounded) = {exp_mean_relative:.3f} kcal/mol"
    )
    for target_temperature, mean_relative, global_delta_g in zip(
        target_temperatures, sweep_mean_relative, global_delta_g_by_temperature
    ):
        print(
            f"Upside T={float(target_temperature):.4f}: "
            f"DeltaG_global={float(global_delta_g):.3f}, "
            f"mean DeltaG_HX - DeltaG_global={float(mean_relative):.3f} kcal/mol"
        )


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    input_dir = root / "simulations" / args.protein / "inputs"
    run_dir = (
        args.run_dir.expanduser().resolve()
        if args.run_dir is not None
        else root / "simulations" / args.protein / "outputs" / args.sim_id
    )
    run_files = discover_run_files(args, run_dir)
    print(f"trajectory directory: {run_dir}")
    print(f"discovered {len(run_files)} replica files")
    result_dir = root / "docs" / "data" / "upside_hdx"
    result_dir.mkdir(parents=True, exist_ok=True)

    hdx_top = ensure_hdx_top(args, root, input_dir)
    residue_file = result_dir / f"{args.protein}_stride{args.stride}.resid"
    for replica, run_file in enumerate(run_files):
        print(f"protection states replica {replica}")
        run_protection_state(
            args,
            hdx_top,
            run_file,
            result_dir / f"{args.protein}_{args.sim_id}_stride{args.stride}_{replica}_PS.npy",
            residue_file,
        )

    potentials, temperatures, protection = load_arrays(args, run_files, result_dir)
    state_masks = load_state_masks(args, run_files)
    potentials, temperatures, protection, state_masks = trim_to_common_frames(
        potentials, temperatures, protection, state_masks
    )
    residue_ids = np.loadtxt(residue_file, dtype=int)
    delta_g, dghx_den, m_values = calculate_hdx(args, potentials, temperatures, protection)
    write_outputs(args, root, input_dir, residue_ids, delta_g, dghx_den, m_values)
    (
        target_temperatures,
        delta_g_by_temperature,
        global_delta_g_by_temperature,
        nse_probability_by_temperature,
        dse_probability_by_temperature,
    ) = calculate_delta_g_temperature_sweep(
        args, potentials, temperatures, protection, state_masks
    )
    write_temperature_sweep_outputs(
        args,
        root,
        input_dir,
        residue_ids,
        target_temperatures,
        delta_g_by_temperature,
        global_delta_g_by_temperature,
        nse_probability_by_temperature,
        dse_probability_by_temperature,
    )
    print("wrote", root / "docs" / "data" / f"{args.protein}_upside_hdx.csv")
    print("wrote", root / "docs" / "figures" / f"{args.protein}_upside_hdx_delta_g_overlay.png")


if __name__ == "__main__":
    main()
