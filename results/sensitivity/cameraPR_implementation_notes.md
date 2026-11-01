# Installed cameraPR implementation check

Inspected limma 3.68.5 cameraPR.default and rankSumTestWithCorrelation, together with the installed camera help topic.

The method default is inter.gene.cor=0.01. With use.ranks=TRUE, cameraPR calls rankSumTestWithCorrelation using infinite degrees of freedom, then forms a two-sided competitive P value from its two directional tails. Direction indicates the relative signed ranking of set genes against non-set genes. BH correction covers the tested collection. This workflow uses an externally fixed correlation scenario and does not estimate correlation from donor expression.

The function default is the reason for the original 0.01 setting; it is not evidence that this cohort's true within-set correlation equals 0.01. Rank-biserial effect and direction remain fixed across the four correlation scenarios because the ranking and membership do not change.

Source documentation: installed limma help(camera); Wu and Smyth (2012), DOI 10.1093/nar/gks461. Full implementation inspection output is retained with internal revision records, not redistributed as original project code.
