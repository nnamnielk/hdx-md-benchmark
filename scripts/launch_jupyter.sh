#!/bin/bash

# Determine script directory and project root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

# Locate miniconda and activate upside2-env
CONDA_BASE="$HOME/.miniconda3"

if [ -f "$CONDA_BASE/bin/activate" ]; then
    source "$CONDA_BASE/bin/activate" upside2-env
elif [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate upside2-env
else
    echo "Error: Miniconda installation not found at $CONDA_BASE" >&2
    exit 1
fi

# Set PYTHONPATH to include project root and scripts directory
export PYTHONPATH="$PROJECT_ROOT/scripts:$PROJECT_ROOT:${PYTHONPATH}"

echo "=================================================="
echo " Starting Jupyter Server in upside2-env"
echo " Working directory : $PROJECT_ROOT"
echo " Python binary     : $(which python)"
echo " Jupyter binary    : $(which jupyter-notebook || which jupyter-lab)"
echo " Colab Enabled     : Yes (port 8888)"
echo "=================================================="

# Execute jupyter notebook passing Colab parameters and any additional CLI arguments
exec jupyter notebook \
  --no-browser \
  --ServerApp.jpserver_extensions="jupyter_http_over_ws=True" \
  --NotebookApp.allow_origin='https://colab.research.google.com' \
  --ServerApp.allow_origin='https://colab.research.google.com' \
  --port=8888 \
  --NotebookApp.port_retries=0 \
  --ServerApp.port_retries=0 \
  --notebook-dir=/ \
  --YDocExtension.disable_collaboration=True
  --NotebookApp.allow_credentials=True \
  --ServerApp.allow_credentials=True \
  --NotebookApp.disable_check_xsrf=True \
  --ServerApp.disable_check_xsrf=True \
  --ServerApp.websocket_ping_interval=30 \
  --ServerApp.websocket_ping_timeout=30 \
  "$@"
