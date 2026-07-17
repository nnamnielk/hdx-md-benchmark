#!/usr/bin/env python
"""Configure and run a single Upside protein simulation."""
import sys, os, shutil, subprocess as sp
import numpy as np

name = sys.argv[1]
UPSIDE_HOME = os.environ.get('UPSIDE_HOME', os.path.expanduser('~/.hermes/sandbox/upside2-md'))
upside_utils_dir = os.path.join(UPSIDE_HOME, "py")
sys.path.insert(0, upside_utils_dir)
import run_upside as ru

base_dir = os.path.join(os.path.dirname(__file__), name)
input_dir  = os.path.join(base_dir, "inputs")

FF = 'ff_2.1'
T = 0.8
DURATION = 5000
FRAME_INTERVAL = 50
SEED = np.random.randint(0, 100000)

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
    initial_structure         = os.path.join(input_dir, f"{name}.initial.npy"),
)

cb_file = os.path.join(input_dir, f"{name}.chain_breaks")
if os.path.exists(cb_file):
    kwargs['chain_break_from_file'] = cb_file

sim_id = f"{name}_run"
run_dir = os.path.join(base_dir, "outputs", sim_id)
os.makedirs(run_dir, exist_ok=True)

fasta = os.path.join(input_dir, f"{name}.fasta")
config_base = os.path.join(input_dir, f"{name}.up")

print(f"[{name}] Configuring ({len(open(fasta).read().split()[1])} res, T={T}, dur={DURATION})...")
ru.upside_config(fasta, config_base, **kwargs)

h5_file  = os.path.join(run_dir, f"{name}.run.up")
log_file = os.path.join(run_dir, f"{name}.run.log")
shutil.copyfile(config_base, h5_file)

cmd = f"{UPSIDE_HOME}/obj/upside --duration {DURATION} --frame-interval {FRAME_INTERVAL} --temperature {T} --seed {SEED} {h5_file}"
print(f"[{name}] Running...")
with open(log_file, 'w') as log:
    sp.check_call(cmd, shell=True, stdout=log, stderr=sp.STDOUT)
print(f"[{name}] ✓ Done")
