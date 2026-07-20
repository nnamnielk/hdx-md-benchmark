import json
import sys

nb_path = sys.argv[1]
with open(nb_path, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        new_source = []
        for line in cell['source']:
            if 'RUNS = os.path.expanduser("~/.hermes/workspace/hdx-benchmark-desmond-upside/runs")' in line:
                new_source.append('RUNS = os.path.expanduser("~/.hermes/workspace/hdx-md-benchmark/simulations")\n')
            elif 'run_dir = os.path.join(RUNS, subdir, "outputs", f"{subdir}_remd")' in line:
                new_source.append('    run_dir = os.path.join(RUNS, subdir, "outputs", "REMD")\n')
            elif '"λ D14A (80 res)":        "lambda_d14a",' in line:
                pass # skip
            elif '"λ YA (1LMB:3, 80 res)":  "lambda_repressor",' in line:
                pass # skip
            else:
                new_source.append(line)
        cell['source'] = new_source

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)
