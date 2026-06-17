#!/usr/bin/env python
"""Launch 4-protein REMD benchmark matching example 04.HDX settings.
nug2b, ubiquitin, lambda_repressor (YA), lambda_d14a.
16 replicas, sqrt-spaced T=0.80-0.95, swap_table2d neighbor exchange.
"""
import sys, os, shutil, subprocess as sp
import numpy as np
from math import sqrt

UPSIDE_HOME = os.environ.get('UPSIDE_HOME', os.path.expanduser('~/.hermes/sandbox/upside2-md'))
upside_utils_dir = os.path.join(UPSIDE_HOME, "py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

# ── Parameters (mirroring example/04.HDX) ─────────────────────
N_REPLICAS = 16
T_LOW, T_HIGH = 0.80, 0.95
DURATION = 200000          # matches HDX example
FRAME_INTERVAL = 100       # matches HDX example
REPLICA_INTERVAL = 10      # matches HDX example
FF = 'ff_2.1'
BASE_SEED = 42

# Square-root spacing (matching example)
temperatures = np.linspace(sqrt(T_LOW), sqrt(T_HIGH), N_REPLICAS)**2
temps_str = ','.join(f'{t:.4f}' for t in temperatures)
swap_sets = ru.swap_table2d(1, len(temperatures))
print(f"Temperatures ({N_REPLICAS}): {[f'{t:.3f}' for t in temperatures]}")
print(f"Swap sets: {swap_sets}")

param_dir_base = os.path.join(UPSIDE_HOME, "parameters/")
param_dir_common = param_dir_base + "common/"
param_dir_ff = param_dir_base + f'{FF}/'

common_kwargs = dict(
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
    ('nug2b',            True),
    ('ubiquitin',        True),
    ('lambda_repressor', True),
    ('lambda_d14a',      True),
]

procs = []

for idx, (name, is_native) in enumerate(proteins):
    base_dir = os.path.join(os.path.dirname(__file__), name)
    input_dir  = os.path.join(base_dir, "inputs")
    fasta = os.path.join(input_dir, f"{name}.fasta")
    config_base = os.path.join(input_dir, f"{name}.up")
    
    sim_id = f"{name}_remd"
    run_dir = os.path.join(base_dir, "outputs", sim_id)
    os.makedirs(run_dir, exist_ok=True)
    
    n_res = len(open(fasta).read().split('\n')[1].strip())
    print(f"\n[{name}] {n_res} res — configuring...")
    
    prot_kwargs = dict(**common_kwargs)
    if is_native:
        prot_kwargs['initial_structure'] = os.path.join(input_dir, f"{name}.initial.npy")
    cb_file = os.path.join(input_dir, f"{name}.chain_breaks")
    if os.path.exists(cb_file):
        prot_kwargs['chain_break_from_file'] = cb_file
    
    ru.upside_config(fasta, config_base, **prot_kwargs)
    
    # One HDF5 per replica
    h5_files = []
    for rep in range(N_REPLICAS):
        h5_file = os.path.join(run_dir, f"{name}.run.{rep}.up")
        shutil.copyfile(config_base, h5_file)
        h5_files.append(h5_file)
    files_str = ' '.join(h5_files)
    
    log_file = os.path.join(run_dir, f"{name}.run.log")
    seed = BASE_SEED + idx
    
    cmd = (
        f"{UPSIDE_HOME}/obj/upside "
        f"--duration {DURATION} "
        f"--frame-interval {FRAME_INTERVAL} "
        f"--temperature {temps_str} "
        f"--replica-interval {REPLICA_INTERVAL} "
        f"--swap-set {swap_sets[0]} "
        f"--swap-set {swap_sets[1]} "
        f"--seed {seed} "
        f"{files_str}"
    )
    
    print(f"[{name}] Launching REMD ({N_REPLICAS} reps, T=%.2f-%.2f, dur={DURATION})..." % (T_LOW, T_HIGH))
    
    with open(log_file, 'w') as log:
        p = sp.Popen(cmd, shell=True, stdout=log, stderr=sp.STDOUT)
    
    procs.append((name, p, log_file))
    print(f"[{name}] PID {p.pid}")

print(f"\n{'='*60}")
print(f"4 REMD runs, {N_REPLICAS} reps each = {4*N_REPLICAS} total")
print(f"T={T_LOW}-{T_HIGH}, sqrt-spaced, replica_interval={REPLICA_INTERVAL}")
print(f"Waiting...")
print(f"{'='*60}")

for name, p, log_file in procs:
    rc = p.wait()
    status = "✓" if rc == 0 else f"✗ (rc={rc})"
    print(f"  {status} {name}")
    
print("\nALL DONE")
