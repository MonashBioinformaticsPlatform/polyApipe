
# 2025-02-10

* Import directly from DelayedArray, as this seems to be necessary: rowSums colSums rowMeans colMeans

# 2024-04-09

* No longer import from BiocGenerics (not needed any more?): rowSums colSums rowMeans colMeans

# 2021-06-26

* Add option do_logNormCounts, and try to not error-out on completely empty cells.

# 2021-06-14

* Default to unstranded gaps (i.e. don't extend into gene on opposite strand).
* peak_min_prop added, default 0.01, removes sites with a proportionately low total count from genes.
* computeSumFactors is done using quickCluster, with quickCluster blocked by batch, sizeFactors only computed once (gene level).
* Various changes aiming to reduce memory usage.
