# Dataset audit: PASS with interpretation limitations

## Matrix and provenance
Official GEO GSE174367 processed snRNA-seq inputs; URLs, byte sizes and SHA256 hashes are in `data/input_manifest.json`. The matrix is **58,721 genes x 61,770 nuclei**, sparse CSC, with 177,178,648 stored integer entries. Every count block was checked for nonnegative integer values (stored range 1–6767), valid indices and unique coordinates. Stored row indices are not always sorted; sorting restores canonical sparse representation without altering counts. Gene IDs are unique; 1,588 duplicate-symbol entries are retained at gene-ID level and handled explicitly for enrichment.

The HDF5 reference label is `GRCh38.p12.premrna2`. Original Supplementary Note (pp. 11–12) describes Cell Ranger 3.0.2, a GRCh38.p12 pre-mRNA reference and `cellranger aggr` read-depth normalization. These are untransformed integer UMI counts **after the original aggregation/subsampling**, not untouched sequencing reads. Exact annotation release is not stated and is not inferred from the genome label. No new alignment, original-QC reproduction or doublet-removal analysis is claimed.

## Barcode reconciliation
All **61,472** final author metadata barcodes match the HDF5 exactly; no prefix/suffix manipulation is needed. The **298 matrix-only barcodes** are explicitly excluded because final donor/cell-type labels are absent and are listed in `matrix_only_barcodes.tsv`. This count is consistent with the original Supplementary Note reporting removal of 298 snRNA-seq nuclei in small clusters. Their per-barcode exclusion reasons cannot be reconstructed. Metadata-only barcodes: zero. The count-matrix total must not be confused with the astrocyte population.

## Donors, annotation and coverage
The original `Cell.Type = ASC` annotation selects **4,756 astrocyte nuclei**. Supplementary Table 1 supplies **18 unique case-year/case-number pairs**, one for each SampleID, and confirms diagnosis, age, sex, PMI and RIN. Thus there are **11 AD and 7 control donors**, with no duplicated biological replicate. No separate technical-library identifier is available in cell metadata. Region is FC (frontal/prefrontal cortex) in the original study. Original donor-case IDs remain pseudonymous.

Every donor has 95–546 astrocytes; the prespecified minimum of 20 excludes no donor. This is a coverage screen, not a statistical power guarantee. Per-nucleus UMI/detected-gene summaries are derived from the count matrix.

|Donor|Diagnosis|GEO batch|Age|Sex|Astrocytes|
|---|---|---:|---:|---|---:|
|Sample-100|Control|1|79|M|459|
|Sample-17|AD|3|90|F|251|
|Sample-19|AD|3|90|F|317|
|Sample-22|AD|2|86|M|154|
|Sample-27|AD|2|81|M|546|
|Sample-33|AD|3|80|M|187|
|Sample-37|AD|3|87|F|302|
|Sample-43|AD|1|90|F|283|
|Sample-45|AD|1|89|F|156|
|Sample-46|AD|2|90|M|95|
|Sample-47|AD|3|90|M|419|
|Sample-50|AD|2|89|F|286|
|Sample-52|Control|3|83|M|150|
|Sample-58|Control|3|90|M|247|
|Sample-66|Control|2|90|F|189|
|Sample-82|Control|2|79|M|148|
|Sample-90|Control|3|79|F|298|
|Sample-96|Control|1|79|M|269|

## Missingness and confounding
Only PMI is missing (Sample-100); age, sex, RIN and GEO Batch are complete. The literal plaque category `None` means no plaques and is retained as a category, not treated as missing. Age values of 90 are top-coded in the source; they are not interpreted as exact ages above that threshold. AD has 6 female/5 male donors versus control 2 female/5 male. Clinical ranges and medians are in `covariate_summary.tsv`. Both diagnoses occur in every GEO batch (AD/control counts: batch 1=2/2, batch 2=4/2, batch 3=5/3).

Two source discrepancies are retained transparently. Supplementary Table 1 reports 62,669 snRNA-seq nuclei, differing from final cell metadata, and its `Library.Batch` assignment differs from GEO `Batch`. The original snRNA processing code uses a field named `Batch`; the table does not identify whether `Library.Batch` refers to the same assay/stage. It would be unjustified to silently replace one with the other. Both tables are cross-referenced in `supplementary_metadata_crosscheck.tsv`; final counts and author cell identities come from the processed GEO files. Neither candidate batch variable is perfectly confounded with diagnosis, but batch effects are not ruled out.

## Inferential gate and model decision
Reliable donors, unambiguous diagnosis, author astrocyte labels, untransformed integer counts, at least five covered donors per group and no demonstrated complete technical confounding satisfy the gate. Before differential expression, the requested **`~ diagnosis`** model is retained as a **marginal, unadjusted group comparison**. Age/sex imbalance and unresolved technical batch provenance preclude disease-specific or causal interpretation. No covariate is selected by gene-level P values. LODO measures donor influence; it cannot remove systematic confounding. The package must carry these limitations explicitly.
