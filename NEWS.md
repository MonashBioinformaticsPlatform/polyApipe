
# 2021-06-14

* Default to unstranded gaps (i.e. don't extend into gene on opposite strand).
* peak_min_prop added, default 0.01, removes sites with a proportionately low total count from genes.
* computeSumFactors is done using quickCluster, with quickCluster blocked by batch, sizeFactors only computed once (gene level).
* Various changes aiming to reduce memory usage.
