#!/bin/bash
set -euo pipefail

# Change only this line before running.
CIF_PATH="/Users/jaekwansmac/Desktop/codex/cifplane/cif/LaMnO3.cif"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/results_plane"
CONDA_ENV_NAME="cifplane"
PNG_BOUNDARY_XMAX=3
PNG_BOUNDARY_YMAX=3
FORCE_REEXPORT_VESTA_PNGS="false"

HKL_LIST=(
  "1 0 0"
  "1 1 1"
)

if [[ ! -f "${CIF_PATH}" ]]; then
  echo "ERROR: CIF file not found: ${CIF_PATH}" >&2
  exit 1
fi

if command -v conda >/dev/null 2>&1; then
  CONDA_BIN="$(command -v conda)"
elif [[ -x "/opt/homebrew/anaconda3/bin/conda" ]]; then
  CONDA_BIN="/opt/homebrew/anaconda3/bin/conda"
else
  echo "ERROR: conda executable was not found." >&2
  exit 1
fi

PYTHON_RUN=("${CONDA_BIN}" run -n "${CONDA_ENV_NAME}" python)
CIF_FILENAME="$(basename "${CIF_PATH}")"
CIF_STEM="${CIF_FILENAME%.*}"
TARGET_RESULT_DIR="${OUTPUT_DIR}/${CIF_STEM}"

echo "[1/3] Running termination analysis for ${CIF_FILENAME}"
for hkl_triplet in "${HKL_LIST[@]}"; do
  echo "  - hkl ${hkl_triplet}"
  "${PYTHON_RUN[@]}" "${SCRIPT_DIR}/analyze_planes.py" \
    --cif "${CIF_PATH}" \
    --hkl ${hkl_triplet} \
    --output-dir "${OUTPUT_DIR}" \
    --termination-only \
    --surface-only-export \
    --png-boundary-xmax "${PNG_BOUNDARY_XMAX}" \
    --png-boundary-ymax "${PNG_BOUNDARY_YMAX}"
done

if [[ "${FORCE_REEXPORT_VESTA_PNGS}" == "true" && -d "${TARGET_RESULT_DIR}" ]]; then
  echo "[2/3] Removing existing VESTA PNGs under ${TARGET_RESULT_DIR}"
  find "${TARGET_RESULT_DIR}" -name '*_top_vesta_view.png' -delete
else
  echo "[2/3] Keeping existing VESTA PNGs when present"
fi

echo "[3/3] Exporting VESTA PNGs"
"${PYTHON_RUN[@]}" "${SCRIPT_DIR}/export_vesta_pngs.py" "${TARGET_RESULT_DIR}"

echo ""
echo "Done."
echo "Results: ${TARGET_RESULT_DIR}"
