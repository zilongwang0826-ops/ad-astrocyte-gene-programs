# Validation and reproduction scope

The included primary analysis uses donor-level pseudobulk counts. The targeted computation completed on 31 October 2026 uses the included derived counts, donor metadata, retained gene mapping and frozen Hallmark membership. It refits the full-donor model and 18 donor omissions; it does not reconstruct the original single-cell processing or provide external biological replication.

## Included evidence

- [Targeted execution summary](results/targeted_validation/validation_summary.json), [baseline comparison](results/targeted_validation/baseline_agreement.tsv), [input and code hashes](results/targeted_validation/input_code_checksums.tsv), and [R runtime](results/targeted_validation/R_sessionInfo.txt).
- Complete [primary tables](results/tables/), [age/correlation sensitivity](results/sensitivity/), and [gene-set donor omissions](results/targeted_validation/gene_set_lodo.tsv).
- [Input provenance](data/README.md), [environment restoration](environment/README.md), and the prespecified [targeted analysis plan](docs/TARGETED_VALIDATION_PLAN.md).

Historical execution hashes identify the files used for those runs. Current distribution hashes in `PUBLIC_FILES.sha256` identify this snapshot; an editorial change does not imply a new scientific execution.

## Check or recompute

After restoring the declared environments, run `bash check.sh` to exercise synthetic fixtures and verify saved results, report claims and figure coordinates. This command requires no raw-data download and does not refit the study. `bash run.sh --targeted-validation --analysis-only` recomputes the targeted analyses from the included derived counts and runs the statistical checks without requiring document fonts. Omit `--analysis-only` to include rendering. `bash run.sh --download` reconstructs the full workflow from checksummed public inputs. Both analysis entry points use two numerical threads.

Figure and report commands are in [the reproduction guide](docs/FIGURE_REPRODUCTION.md). The legacy `tests/test_statistics.R --empirical` command requires `results/tables/primary_fit.rds`, an intermediate intentionally omitted from the distribution; use the supported `check.sh` route for this snapshot.

Donor omission is an internal influence diagnostic. Direction consistency, full-rank design checks and figure/table agreement do not establish external replication or adequate control of demographic and technical confounding. The report retains these interpretation limits.

Current check records are written under `.build/validation/` (or `AD_ASTRO_INTERNAL_DIR/validation/` when configured). Historical numerical baselines and saved validation records are retained. Application wording is reviewed separately and is not reconstructed from earlier editing decisions by the routine scientific checker.
