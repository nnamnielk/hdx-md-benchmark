import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import upside_hdx_nug2b

def main():
    root = Path(__file__).resolve().parents[1]
    run_dir = root / "simulations" / "nug2b" / "outputs" / "remd_lowtemp"
    
    class Args:
        protein = "nug2b"
        stride = 1
        n_replicas = 20
    args = Args()
    
    run_files = upside_hdx_nug2b.discover_run_files(args, run_dir)
    print(f"Found {len(run_files)} replica files in {run_dir}.")
    
    masks = []
    temps = []
    for i, run_file in enumerate(run_files):
        print(f"Processing replica {i}...")
        mask = upside_hdx_nug2b.state_mask_for_run("nug2b", run_file, stride=1)
        temp = upside_hdx_nug2b.read_output_series(run_file, "temperature", stride=1).reshape(-1)
        min_len = min(len(mask), len(temp))
        masks.append(mask[:min_len])
        temps.append(temp[:min_len])
        
    flat_masks = np.concatenate(masks)
    flat_temps = np.concatenate(temps)
    unique_temps = np.unique(np.round(flat_temps, 4))
    
    p_native = []
    valid_temps = []
    for t in unique_temps:
        idx = np.isclose(flat_temps, t, atol=1e-4)
        if np.sum(idx) > 0:
            frac = np.mean(flat_masks[idx])
            p_native.append(frac)
            valid_temps.append(t)
            print(f"T = {t:.4f} : P(native) = {frac:.4f} (Frames: {np.sum(idx)})")
            
    valid_temps = np.array(valid_temps)
    p_native = np.array(p_native)
    
    tm = np.nan
    try:
        for i in range(len(p_native)-1):
            if (p_native[i] >= 0.5 and p_native[i+1] <= 0.5) or (p_native[i] <= 0.5 and p_native[i+1] >= 0.5):
                t1, t2 = valid_temps[i], valid_temps[i+1]
                p1, p2 = p_native[i], p_native[i+1]
                tm = t1 + (0.5 - p1) * (t2 - t1) / (p2 - p1)
                break
    except Exception as e:
        print(f"Could not interpolate Tm: {e}")

    print(f"\nEstimated Melting Temperature (Tm): {tm:.4f}")

    plt.figure(figsize=(8, 5))
    plt.plot(valid_temps, p_native, marker='o', linestyle='-', color='b')
    if not np.isnan(tm):
        plt.axvline(tm, color='r', linestyle='--', label=f'Tm = {tm:.4f}')
        plt.axhline(0.5, color='gray', linestyle=':')
    plt.xlabel('Temperature (Upside Units)')
    plt.ylabel('Probability of Native State')
    plt.title('NuG2b Melting Curve from T-REMD (0.6 - 0.8)')
    plt.legend()
    plt.grid(True)
    out_path = root / "analysis" / "nug2b_melting_curve_lowtemp2.png"
    plt.savefig(out_path)
    print(f"Plot saved to {out_path}")

if __name__ == "__main__":
    main()
