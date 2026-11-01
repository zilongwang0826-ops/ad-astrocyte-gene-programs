# Project descriptions

Project period: August 2026.

## Short description

This project reanalysed 4,756 astrocyte nuclei from 18 donors to compare Alzheimer's disease with controls. Donor-level edgeR analysis found no significant individual genes. Inflammatory and protein-stress gene sets ranked toward AD, but their statistical support changed with correlation assumptions, age adjustment and individual donors.

## Extended description

This reanalysis used public GSE174367 counts and the study's existing astrocyte labels. Counts were combined within 11 AD and 7 control donors, then analysed with edgeR and Hallmark gene-set testing. No gene passed BH FDR < 0.05. Five of 49 gene sets passed at the default correlation of 0.01, but none passed at 0.03 or 0.05. Age adjustment retained three of the original five sets. All five kept their direction in 18 donor-omission fits, while each remained significant in only 9-11 fits. Complete tables and code document these results. Age and sex imbalance, capped ages and unresolved batch information limit biological interpretation.

## Motivation-letter paragraph

This project developed my interest in evaluating computational evidence at the level of independent biological samples. Through a reanalysis of public astrocyte transcriptomes, I examined how donor structure and statistical assumptions affect interpretation of gene-level and gene-set results. I am interested in extending this work in a research setting where computational analyses can be assessed alongside detailed sample metadata and experimental evidence.
