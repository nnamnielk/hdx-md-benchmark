#!/usr/bin/env python3
"""Shared Upside replica-exchange setup and SLURM submission logic."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

N_REPLICAS = 16
T_LOW = 0.80
T_HIGH = 0.96


def parse_args(protein: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Run a 16-replica Upside exchange simulation for {protein}."
    )
    parser.add_argument("--account", default=os.getenv("SLURM_ACCOUNT"))
    parser.add_argument("--partition", default=os.getenv("SLURM_PARTITION"))
    parser.add_argument("--duration", type=float, default=5000)
    parser.add_argument("--frame-interval", type=float, default=50)
    parser.add_argument("--replica-interval", type=float, default=10)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--time", default="36:00:00", dest="run_time")
    parser.add_argument("--force-field", default="ff_2.1")
    parser.add_argument("--sim-id", default="REMD")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument(
        "--local",
        action="store_true",
        help="Run with mpirun instead of submitting through SLURM.",
    )
    return parser.parse_args()


def archive_previous_output(h5_files: list[Path], log_file: Path) -> None:
    try:
        import tables as tb
    except ImportError as exc:
        raise ImportError(
            "Restarting requires PyTables in the active Python environment"
        ) from exc

    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    if log_file.exists():
        shutil.move(log_file, log_file.with_name(f"{log_file.name}.bck_{timestamp}"))

    for filename in h5_files:
        with tb.open_file(filename, "a") as trajectory:
            index = 0
            while f"output_previous_{index}" in trajectory.root:
                index += 1

            if "output" in trajectory.root:
                output = trajectory.root.output
            elif index:
                output = trajectory.get_node(f"/output_previous_{index - 1}")
            else:
                raise RuntimeError(f"{filename} has no trajectory output to restart")

            trajectory.root.input.pos[:, :, 0] = output.pos[-1, 0]
            if "output" in trajectory.root:
                trajectory.root.output._f_rename(f"output_previous_{index}")


def run(protein: str) -> None:
    args = parse_args(protein)
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "Running Upside requires NumPy in the active Python environment"
        ) from exc

    root = Path(__file__).resolve().parents[1]
    pdb_file = root / "pdb" / f"{protein}.pdb"
    work_dir = root / "simulations" / protein
    input_dir = work_dir / "inputs"
    run_dir = work_dir / "outputs" / args.sim_id
    input_dir.mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True, exist_ok=True)

    if not pdb_file.exists():
        raise FileNotFoundError(
            f"{pdb_file} is missing; run scripts/download_pdbs.py first"
        )

    upside_home_value = os.getenv("UPSIDE_HOME")
    if not upside_home_value:
        raise EnvironmentError("UPSIDE_HOME must point to the upside2-md checkout")

    upside_home = Path(upside_home_value).expanduser().resolve()
    upside_utils = upside_home / "py"
    sys.path.insert(0, str(upside_utils))
    import run_upside as ru

    h5_files = [
        run_dir / f"{protein}.run.{replica}.up"
        for replica in range(N_REPLICAS)
    ]
    log_file = run_dir / f"{protein}.run.log"
    config_base = input_dir / f"{protein}.up"
    fasta = input_dir / f"{protein}.fasta"
    chain_breaks = input_dir / f"{protein}.chain_breaks"

    if args.restart:
        missing = [filename for filename in h5_files if not filename.exists()]
        if missing:
            raise FileNotFoundError(
                "Cannot restart; missing replica files: "
                + ", ".join(str(filename) for filename in missing)
            )
    else:
        subprocess.run(
            [
                sys.executable,
                str(upside_utils / "PDB_to_initial_structure.py"),
                str(pdb_file),
                str(input_dir / protein),
                "--record-chain-breaks",
            ],
            check=True,
        )

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
            initial_structure=str(input_dir / f"{protein}.initial.npy"),
        )
        if chain_breaks.exists():
            config_options["chain_break_from_file"] = str(chain_breaks)
        config_stdout = ru.upside_config(
            str(fasta),
            str(config_base),
            **config_options,
        )
        print(config_stdout)
        for filename in h5_files:
            shutil.copyfile(config_base, filename)

    temperatures = np.linspace(np.sqrt(T_LOW), np.sqrt(T_HIGH), N_REPLICAS) ** 2
    temperature_string = ",".join(f"{temperature:.8g}" for temperature in temperatures)
    swap_sets = ru.swap_table2d(1, N_REPLICAS)
    upside_args = [
        "--duration",
        str(args.duration),
        "--frame-interval",
        str(args.frame_interval),
        "--temperature",
        temperature_string,
        "--seed",
        str(args.seed),
        "--replica-interval",
        str(args.replica_interval),
        "--swap-set",
        str(swap_sets[0]),
        "--swap-set",
        str(swap_sets[1]),
    ]

    if args.restart:
        archive_previous_output(h5_files, log_file)

    upside_binary = upside_home / "obj" / "upside"
    simulation_command = [
        str(upside_binary),
        *upside_args,
        *(str(filename) for filename in h5_files),
    ]

    print("Temperatures:", temperature_string)
    if args.local:
        subprocess.run(simulation_command, check=True)
        return

    if not args.account or not args.partition:
        raise ValueError(
            "Set SLURM_ACCOUNT and SLURM_PARTITION or pass --account and --partition"
        )

    sbatch_command = [
        "sbatch",
        f"--account={args.account}",
        f"--partition={args.partition}",
        f"--job-name={protein}_{args.sim_id}",
        f"--output={log_file}",
        f"--time={args.run_time}",
        "--nodes=1",
        f"--ntasks-per-node={N_REPLICAS}",
        "--wrap",
        " ".join(simulation_command),
    ]
    subprocess.run(sbatch_command, check=True)
