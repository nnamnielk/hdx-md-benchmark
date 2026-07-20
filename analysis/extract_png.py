import json
import base64
import sys

nb_path = sys.argv[1]
with open(nb_path, 'r') as f:
    nb = json.load(f)

img_count = 0
for cell in nb['cells']:
    if 'outputs' in cell:
        for output in cell['outputs']:
            if 'data' in output and 'image/png' in output['data']:
                img_data = output['data']['image/png']
                # Sometimes it's a list of strings, sometimes a single string
                if isinstance(img_data, list):
                    img_data = "".join(img_data)
                
                img_bytes = base64.b64decode(img_data)
                out_name = f"/home/okleinmann/.hermes/workspace/hdx-md-benchmark/analysis/notebook_output_{img_count}.png"
                with open(out_name, 'wb') as img_f:
                    img_f.write(img_bytes)
                print(f"Saved {out_name}")
                img_count += 1
