#!/usr/bin/env python
"""Run Upside MD for all 3 benchmark proteins."""
import sys, os, shutil
import subprocess as sp
import numpy as np

UPSIDE_HOME = os.environ.get('UPSIDE_HOME', os.path.expanduser('~/.hermes/sandbox/upside2-md'))
upside_utils_dir = os.path.join(UPSIDE_HOME, "py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

# Parameters
FF = 'ff_2.1'
T = 0.8           # reduced temperature (~300K)
DURATION = 5000   # short test; bump for production
FRAME_INTERVAL = 50
SEED = 42

param_dir_base = os.path.join(UPSIDE_HOME, "parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + f'{FF}/'

kwargs = dict(
    rama_library              = param_dir_common + "rama.dat",
    rama_sheet_mix_energy     = param_dir_ff + "sheet",
    reference_state_rama      = param_dir_common + "rama_reference.pkl",
    hbond_energy              = param_dir_ff + "hbond.h5",
    rotamer_placement         = param_dir_ff + "sidechain.h5",
    dynamic_rotamer_1body     = True,
    rotamer_interaction       = param_dir_ff + "sidechain.h5",
    environment_potential     = param_dir_ff + "environment.h5",
    bb_environment_potential  = param_dir_ff + "bb_env.dat",
)

proteins = [
    ('nug2b',           True),
    ('ubiquitin',       True),
    ('lambda_repressor', True),
]

for name, is_native in proteins:
    base_dir = os.path.join(os.path.dirname(__file__), name)
    input_dir  = os.path.join(base_dir, "inputs")
    output_dir = os.path.join(base_dir, "outputs")
    
    sim_id = f"{name}_run"
    run_dir = os.path.join(output_dir, sim_id)
    os.makedirs(run_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"  {name.upper()}")
    print(f"{'='*60}")
    
    # --- Config ---
    fasta = os.path.join(input_dir, f"{name}.fasta")
    config_base = os.path.join(input_dir, f"{name}.up")
    
    prot_kwargs = dict(**kwargs)
    if is_native:
        prot_kwargs['initial_structure'] = os.path.join(input_dir, f"{name}.initial.npy")
    
    # Add chain_breaks if it exists
    cb_file = os.path.join(input_dir, f"{name}.chain_breaks")
    if os.path.exists(cb_file):
        prot_kwargs['chain_break_from_file'] = cb_file
        print(f"  Using chain breaks: {cb_file}")
    
    print(f"  Residues: {len(open(fasta).read().split()[1])} (fasta)")
    print("  Configuring...")
    config_stdout = ru.upside_config(fasta, config_base, **prot_kwargs)
    
    # --- Run ---
    upside_opts = (
        f"--duration {DURATION} "
        f"--frame-interval {FRAME_INTERVAL} "
        f"--temperature {T} "
        f"--seed {SEED} "
    )
    
    h5_file  = os.path.join(run_dir, f"{name}.run.up")
    log_file = os.path.join(run_dir, f"{name}.run.log")
    shutil.copyfile(config_base, h5_file)
    
    cmd = f"{UPSIDE_HOME}/obj/upside {upside_opts} {h5_file}"
    print(f"  Running: {UPSIDE_HOME}/obj/upside (duration={DURATION}, T={T})")
    
    with open(log_file, 'w') as log:
        sp.check_call(cmd, shell=True, stdout=log, stderr=sp.STDOUT)
    
    print(f"  ✓ Done — output: {h5_file}")
    print(f"  Log: {log_file}")

print(f"\n{'='*60}")
print("  ALL DONE")
print(f"{'='*60}")
