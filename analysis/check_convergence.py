import sys, os, glob
import numpy as np
import tables as tb

sys.path.append(os.path.expanduser("~/.hermes/sandbox/upside2-md/py"))
from mdtraj_upside import _output_groups

RUNS = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/simulations")
out_dir = os.path.join(RUNS, "nug2b/outputs/REMD")
files = glob.glob(os.path.join(out_dir, "*.up"))

print(f"Found {len(files)} trajectory files.")

total_frames = 0
native_frames = 0
rg_native = []
temp_visits = np.zeros(20)

for f in files:
    try:
        with tb.open_file(f, 'r') as t:
            for g_no, g in enumerate(_output_groups(t)):
                r_idx = g.replica_index[:] if g_no == 0 else g.replica_index[1:]
                
                temps = r_idx[:, 0]
                total_frames += len(temps)
                
                for t_idx in temps:
                    temp_visits[t_idx] += 1
                
                native_mask = (temps == 0)
                native_frames += np.sum(native_mask)
    except Exception as e:
        print(f"Error reading {f}: {e}")

print(f"Total frames: {total_frames}")
print(f"Native frames (T index 0): {native_frames}")

print("Temperature visits across all replicas:")
for i, count in enumerate(temp_visits):
    print(f"  T={i}: {count} ({100 * count / total_frames:.1f}%)")
