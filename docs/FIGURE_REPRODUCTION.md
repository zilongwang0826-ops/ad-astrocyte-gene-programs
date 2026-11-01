# Regenerate figures from included results

Use the locked dependencies and installed Times New Roman/Arial fonts documented in environment/README.md. The historical figure revision used Python 3.12.14. The targeted-validation work copy uses Python 3.12.16 with the same pinned packages; current execution evidence is recorded separately from the historical snapshot. No additional plotting dependency is required.

    export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 MKL_NUM_THREADS=2
    export PYTHONDONTWRITEBYTECODE=1
    export MPLCONFIGDIR="$PWD/.build/matplotlib"
    .venv/bin/python scripts/05_figures.py
    .venv/bin/python scripts/06_documents.py
    .venv/bin/python tests/verify_figure_redesign.py
    .venv/bin/python tests/verify_results.py
    AD_ASTRO_PYTHON="$PWD/.venv/bin/python" ./check.sh

Plot/document commands read included tables and do not download data or refit scientific models. The full run.sh is a different workflow and is unnecessary for a figure rebuild. Existing check.sh includes small synthetic fixtures and saved-result checks.

Colors are in scripts/figure_style.json; dimensions, labels, axes and deterministic ordering are explicit in scripts/05_figures.py. The four main figures retain their original scientific content. Two supplementary figures read the complete targeted-validation outputs: donor MDS coordinates, and the 49-by-18 gene-set omission table with its summary. Both supplementary displays require results/targeted_validation to be populated by the statistical workflow; the plotting script never creates missing statistical results.

No random subset or inferred uncertainty is used. Outputs have stable names in results/figures; plotted values are in results/plot_data and source hashes/artist geometry in results/audit/figure_qa. PDF/SVG retain vector text/geometry; PNG is 400 dpi. Figure 3 uses larger labels and the full document content width, yielding at least 8-point text in both final PDFs. Its full primary-P-ranked ten-set display and all evidence values are retained.

After rebuilding, render both PDFs and inspect every page and standalone figure. The summary is required to remain one page. The document builder resolves the research report, summary and supplementary methods from the claim registry. It preserves the separate CV bullets and application prose; those documents have their own review scope.

The historical BASELINE_NUMERICAL_HASHES.json remains unchanged. The figure verifier protects those numerical and source files, preserves the prior claim registry as PRE_TARGETED_REPORTED_CLAIMS.json, and checks every pre-existing claim while allowing documented new validation claims. Source-to-artist checks cover the two supplementary plots as well as the four main plots. These computational checks complement, rather than replace, visual review of the rendered pages.

Application wording is outside the numerical baseline check. Its earlier hashes remain in the historical baseline, and earlier checker implementations remain in repository history. Current numerical and plotted-coordinate checks continue to run; their new records go under `.build/validation/`, leaving the saved historical figure-validation record unchanged. `AD_ASTRO_INTERNAL_DIR` can redirect current check records outside the project.
