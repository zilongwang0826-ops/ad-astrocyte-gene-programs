#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
PYTHON="${AD_ASTRO_PYTHON:-$PWD/.venv/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  echo 'Missing project Python environment. Follow environment/README.md.' >&2
  exit 2
fi
command -v Rscript >/dev/null || { echo 'Rscript is required.' >&2; exit 2; }
"$PYTHON" -c 'import numpy,pandas,scipy,pymupdf'

# These suites need the included result tables, but no raw inputs or saved fit.
"$PYTHON" -m unittest discover -s tests -p test_data.py -v
Rscript --vanilla tests/test_statistics.R
Rscript --vanilla tests/test_sensitivity.R --empirical
Rscript --vanilla tests/test_targeted_validation.R --empirical
"$PYTHON" tests/verify_results.py
"$PYTHON" tests/verify_figure_redesign.py
echo 'Distribution checks passed. This check did not rerun the study analysis.'
