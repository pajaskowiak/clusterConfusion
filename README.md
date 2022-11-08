# clusterConfusion
  
<!-- PROJECT LOGO -->
<br />
<div align="center">
<h3 align="center">Cluster Confusion</h3>
<p align="center"> Clustering validation that makes sense.  </p>
</div>

<!-- ABOUT THE PROJECT -->
## About The Project

This package provides an implementation of the Area Under the (ROC) Curve for clustering (AUCC), which can be employed to validate clustering results and select the best number of clusters from a pool of partitions. The measure is described in our ECML/PKDD paper, which was published in the Journal track:

 *Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022).* Available at: https://doi.org/10.1007/s10618-022-00829-0

More clustering validation measures and cool stuff coming soon.

<!-- Installation -->
## Installation

```{r}
install.packages("devtools")
devtools::install_github("https://github.com/pajaskowiak/clusterConfusion")
require(clusterConfusion)
```

<!-- USAGE EXAMPLES -->
## Usage

Just open an R Session and check usage examples as below:

```{r}
library(clusterConfusion)
?aucc
```


<!-- CITE -->
## How to Cite

Please, if you use our method cite our paper:

[1] *Jaskowiak, P.A., Costa, I.G. & Campello, R.J.G.B. The area under the ROC curve as a measure of clustering quality. Data Mining and Knowledge Discovery 36, 1219–1245 (2022).* Available at: https://doi.org/10.1007/s10618-022-00829-0

If you also use our implementation, consider citing our git repository as well.

<!-- CONTACT -->
## Contact

Pablo Andretta Jaskowiak - pablo.andretta@ufsc.br

Project Link: [https://github.com/pajaskowiak/clusterConfusion](https://github.com/pajaskowiak/clusterConfusion)
