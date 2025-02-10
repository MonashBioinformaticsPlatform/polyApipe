
# 2025-02-10

* Import rowSums,colSums,rowMeans,colMeans from DelayedArray. (DelayedArray seems to turn the base implementation of these into generics, but polyApiper wasn't loading the generics properly.)

# 2021-06-26

* Add option do_logNormCounts, and try to not error-out on completely empty cells.

# 2021-06-14

* Default to unstranded gaps (i.e. don't extend into gene on opposite strand).
* peak_min_prop added, default 0.01, removes sites with a proportionately low total count from genes.
* computeSumFactors is done using quickCluster, with quickCluster blocked by batch, sizeFactors only computed once (gene level).
* Various changes aiming to reduce memory usage.
