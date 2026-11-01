# Data provenance and reuse

The analysis uses **GSE174367**, Morabito et al., Nature Genetics (2021), DOI [10.1038/s41588-021-00894-z](https://doi.org/10.1038/s41588-021-00894-z). Only the processed snRNA-seq HDF5 and final author cell metadata were downloaded. No FASTQ, snATAC or bulk RNA-seq data are included.

`input_manifest.json` records the exact GEO and MSigDB URLs, acquisition date, byte sizes and SHA256 hashes. Input files are kept unmodified in ignored `data/raw/`. `scripts/verify_inputs.py --download` fetches absent inputs only and verifies their expected hashes; a conflicting existing file is preserved and causes an error.

The original count matrix is sparse CSC with 58,721 feature rows and 61,770 nucleus columns. The final metadata has 61,472 nuclei. The 298 unmatched HDF5 barcodes are listed explicitly in the dataset audit. ASC annotations identify 4,756 nuclei in 18 independent donors. Counts are untransformed UMI integers after the authors' Cell Ranger 3.0.2 aggregation and read-depth subsampling, using a GRCh38.p12 pre-mRNA reference (`GRCh38.p12.premrna2`). The exact GTF annotation release is not specified in the primary methods.

`reference/supplementary_donor_metadata.tsv` is a compact factual excerpt of original Supplementary Table 1, retaining the 18 rows with `Assay=both`. Its source workbook URL, SHA256 and derivation are in `reference/provenance.json`. Unique case-year/case-number pairs establish the independent donor mapping. Counts and batch labels in that earlier table differ from the final deposited cell metadata; see `results/audit/dataset_audit.md`. These fields are not silently reconciled or substituted.

MSigDB **Hallmark 2026.1.Hs** was retrieved from the official Broad Institute release directory. The collection is credited to MSigDB/Broad Institute and the authors of Liberzon et al. (2015), [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004). The applicable [MSigDB license terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) specify CC BY 4.0 for this release. Gene sets remain separate from this repository's code license.

GEO inputs and original publication materials retain their source terms and attribution; this project does not relicense them. Original processed files remain outside Git. Public compact result tables are derived by this project's code, not copied from the publication. Public pseudonymous donor metadata are used only for the stated research comparison.
