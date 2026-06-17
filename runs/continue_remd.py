#!/usr/bin/env python3
"""Continue nug2b and lambda_d14a REMD simulations from last checkpoint.

Restarts each replica from its last frame, preserving replica exchange state.
Runs 2 proteins in parallel with OMP_NUM_THREADS=16 each (32 cores total).
"""

import sys, os, shutil, subprocess as sp, time
import numpy as np
import tables as tb

UPSIDE_HOME = os.environ.get('UPSIDE_HOME', os.path.expanduser('~/.hermes/sandbox/upside2-md'))
sys.path.insert(0, os.path.join(UPSIDE_HOME, 'py'))
import run_upside as ru

# ── Configuration ──────────────────────────────────────────────────
RUNS_DIR = os.path.expanduser('~/.hermes/workspace/hdx-benchmark-desmond-upside/runs')
PROTEINS = ['nug2b', 'lambda_d14a']  # 2 out of 4
N_REPLICAS = 16
TOTAL_DURATION = 200000
FRAME_INTERVAL = 100
REPLICA_INTERVAL = 10
T_LOW, T_HIGH = 0.80, 0.95
FF = 'ff_2.1'
OMP_THREADS_PER_PROTEIN = 16  # 2 proteins × 16 = 32 cores total

# Temperature schedule (sqrt-spaced, matching original)
temperatures = np.linspace(np.sqrt(T_LOW), np.sqrt(T_HIGH), N_REPLICAS)**2
temps_str = ','.join(f'{t:.4f}' for t in temperatures)
swap_sets = ru.swap_table2d(1, len(temperatures))

print(f"Temperatures: {[f'{t:.3f}' for t in temperatures]}")
print(f"Swap sets: {swap_sets}")
print(f"OMP_NUM_THREADS per protein: {OMP_THREADS_PER_PROTEIN}")

procs = []

for idx, name in enumerate(PROTEINS):
    run_dir = os.path.join(RUNS_DIR, name, 'outputs', f'{name}_remd')
    log_file = os.path.join(run_dir, f'{name}.run.log')
    
    # ── Determine progress ──────────────────────────────────────
    last_step = 0
    try:
        with open(log_file) as f:
            lines = f.readlines()
            # Read backwards to find the last good progress line
            for line in reversed(lines):
                if 'elapsed' in line and '/' in line:
                    try:
                        step = int(line.split()[0])
                        last_step = max(last_step, step)
                    except (ValueError, IndexError):
                        continue
            # Require at least one valid line
            if last_step == 0:
                # Try forward scan
                for line in lines:
                    if 'elapsed' in line and '/' in line:
                        try:
                            step = int(line.split()[0])
                            if step > last_step:
                                last_step = step
                        except (ValueError, IndexError):
                            continue
        if last_step == 0:
            print(f"[{name}] WARNING: Could not parse step from log, assuming fresh start")
    except FileNotFoundError:
        print(f"[{name}] No log file found, assuming fresh start")
    
    remaining = TOTAL_DURATION - last_step
    print(f"\n[{name}] Progress: {last_step}/{TOTAL_DURATION} ({100*last_step/TOTAL_DURATION:.1f}%)")
    print(f"[{name}] Remaining: {remaining} steps ({remaining//FRAME_INTERVAL} frames)")
    
    if remaining <= 0:
        print(f"[{name}] Already complete! Skipping.")
        continue
    
    # ── Prepare HDF5 files for restart ──────────────────────────
    h5_files = []
    for rep in range(N_REPLICAS):
        h5_file = os.path.join(run_dir, f'{name}.run.{rep}.up')
        if not os.path.exists(h5_file):
            print(f"[{name}] ERROR: Missing {h5_file}")
            sys.exit(1)
        
        with tb.open_file(h5_file, 'a') as t:
            # Copy last position to input
            last_pos = t.root.output.pos[-1, 0]  # shape: (n_res, 3)
            t.root.input.pos[:, :, 0] = last_pos
            
            # Rename old output so new run starts fresh
            prev_idx = 0
            while f'output_previous_{prev_idx}' in t.root:
                prev_idx += 1
            t.root.output._f_rename(f'output_previous_{prev_idx}')
        
        h5_files.append(h5_file)
    
    # ── Archive old log ─────────────────────────────────────────
    if os.path.exists(log_file):
        ts = time.strftime('%Y-%m-%d_%H-%M-%S')
        shutil.move(log_file, f'{log_file}.bck_{ts}')
    
    # ── Build command ───────────────────────────────────────────
    seed = 42 + idx
    files_str = ' '.join(h5_files)
    
    cmd = (
        f"OMP_NUM_THREADS={OMP_THREADS_PER_PROTEIN} "
        f"{UPSIDE_HOME}/obj/upside "
        f"--duration {remaining} "
        f"--frame-interval {FRAME_INTERVAL} "
        f"--temperature {temps_str} "
        f"--replica-interval {REPLICA_INTERVAL} "
        f"--swap-set {swap_sets[0]} "
        f"--swap-set {swap_sets[1]} "
        f"--seed {seed} "
        f"{files_str}"
    )
    
    print(f"[{name}] Launching REMD ({N_REPLICAS} reps, {remaining} steps)...")
    
    with open(log_file, 'w') as log:
        p = sp.Popen(cmd, shell=True, stdout=log, stderr=sp.STDOUT)
    
    procs.append((name, p, log_file))
    print(f"[{name}] PID {p.pid} ✓")

print(f"\n{'='*60}")
print(f"2 REMD runs ({N_REPLICAS} replicas each) continuing on 32 cores")
print(f"Remaining: nug2b={TOTAL_DURATION - (TOTAL_DURATION if len(PROTEINS)<2 else 0)} steps, lambda_d14a steps")
print(f"{'='*60}")

if not procs:
    print("Nothing to launch — all done!")
    sys.exit(0)

# ── Monitor ─────────────────────────────────────────────────────
print("Monitoring (Ctrl+C to detach)...")
try:
    for name, p, log_file in procs:
        rc = p.wait()
        status = "✓" if rc == 0 else f"✗ (rc={rc})"
        print(f"  {status} {name}")
except KeyboardInterrupt:
    print("\nDetached. Simulations continue in background.")
    for name, p, _ in procs:
        print(f"  {name}: PID {p.pid}")

print("\nDONE")
