#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
download=0
revision_only=0
targeted_validation=0
analysis_only=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --download) download=1 ;;
    --revision-only) revision_only=1 ;;
    --targeted-validation) targeted_validation=1 ;;
    --analysis-only) analysis_only=1 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done
if [[ "$targeted_validation" -eq 1 && ( "$download" -eq 1 || "$revision_only" -eq 1 ) ]]; then
  echo '--targeted-validation uses distributed donor counts; do not combine it with other modes.' >&2
  exit 2
fi
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/.build/matplotlib}"
export AD_ASTRO_INTERNAL_DIR="${AD_ASTRO_INTERNAL_DIR:-$PWD/.build}"
mkdir -p "$AD_ASTRO_INTERNAL_DIR/logs" "$MPLCONFIGDIR"
PYTHON="${AD_ASTRO_PYTHON:-$PWD/.venv/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  echo 'Missing project Python environment. Follow environment/README.md.' >&2
  exit 2
fi
# Document prerequisites apply only when figures and reports are requested.
if [[ "$analysis_only" -eq 0 ]]; then
  "$PYTHON" -c 'from scripts.typography import font_directory; print("Document fonts:", font_directory())'
  "$PYTHON" -c 'import matplotlib,reportlab,pypdf'
fi
command -v Rscript >/dev/null || { echo 'Rscript is required.' >&2; exit 2; }
"$PYTHON" -c 'import numpy,pandas,scipy,h5py'
Rscript -e 'lib<-Sys.getenv("AD_ASTRO_R_LIBRARY", "environment/R-library"); .libPaths(c(lib,.libPaths())); for(p in c("edgeR","limma","statmod","jsonlite")) stopifnot(requireNamespace(p,quietly=TRUE))'
run_stage() {
  local label="$1"; shift
  echo "Running $label"
  "$@" 2>&1 | tee "$AD_ASTRO_INTERNAL_DIR/logs/$label.log"
}
if [[ "$targeted_validation" -eq 1 ]]; then
  run_stage test_targeted_synthetic Rscript --vanilla tests/test_targeted_validation.R
  run_stage 08_targeted_validation Rscript --vanilla scripts/08_targeted_validation.R
  run_stage test_targeted_empirical Rscript --vanilla tests/test_targeted_validation.R --empirical
  if [[ "$analysis_only" -eq 1 ]]; then
    echo 'Targeted numerical outputs and statistical checks complete; figures and reports were not rebuilt.'
    exit 0
  fi
  run_stage 05_figures "$PYTHON" scripts/05_figures.py
  run_stage 06_documents "$PYTHON" scripts/06_documents.py
  run_stage distribution_checks bash check.sh
  echo 'Targeted validation rebuilt from distributed pseudobulk counts. Raw inputs were not reaggregated.'
  exit 0
fi
if [[ "$download" -eq 1 ]]; then
  "$PYTHON" scripts/verify_inputs.py --download
else
  "$PYTHON" scripts/verify_inputs.py
fi
if [[ "$revision_only" -eq 0 ]]; then
  run_stage test_data "$PYTHON" -m unittest discover -s tests -p test_data.py -v
  run_stage test_statistics Rscript tests/test_statistics.R
  run_stage 00_audit "$PYTHON" scripts/00_audit.py
  run_stage 01_pseudobulk "$PYTHON" scripts/01_pseudobulk.py
  run_stage 02_differential_expression Rscript scripts/02_differential_expression.R
  run_stage 03_enrichment Rscript scripts/03_enrichment.R
  run_stage 04_robustness Rscript scripts/04_robustness.R
fi
run_stage test_sensitivity Rscript tests/test_sensitivity.R
run_stage 07_sensitivity Rscript scripts/07_sensitivity.R
run_stage test_sensitivity_empirical Rscript tests/test_sensitivity.R --empirical
run_stage test_targeted_synthetic Rscript --vanilla tests/test_targeted_validation.R
run_stage 08_targeted_validation Rscript --vanilla scripts/08_targeted_validation.R
run_stage test_targeted_empirical Rscript --vanilla tests/test_targeted_validation.R --empirical
run_stage test_empirical Rscript tests/test_statistics.R --empirical
if [[ "$analysis_only" -eq 1 ]]; then
  echo 'Numerical outputs and statistical checks complete; figures and reports were not rebuilt.'
  exit 0
fi
run_stage 05_figures "$PYTHON" scripts/05_figures.py
run_stage 06_documents "$PYTHON" scripts/06_documents.py
run_stage verify_results "$PYTHON" tests/verify_results.py
run_stage verify_figures "$PYTHON" tests/verify_figure_redesign.py
echo 'Analysis and documents rebuilt. Inspect generated PDFs before sharing.'
