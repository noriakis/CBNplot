## CBNplot: Bayesian network plot for the enrichment analysis results

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/noriakis/CBNplot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/CBNplot/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
 
Plot bayesian network inferred from expression data based on the enrichment analysis results from libraries including `clusterProfiler` and `ReactomePA` using `bnlearn`.

### Installation

```R
library(devtools)
install_github("noriakis/CBNplot")
```

### Usage
- Documentation: [https://noriakis.github.io/software/CBNplot/](https://noriakis.github.io/software/CBNplot/)
- Web server: [https://cbnplot.com](https://cbnplot.com)

### Plot examples

- The gene-to-gene relationship compared with the reference network.
<img src="https://github.com/noriakis/software/blob/main/images/CBNplot_readme_1.png?raw=true" width="800px">

- The plot is customizable highliting edges and nodes like hub genes.
<img src="https://github.com/noriakis/software/blob/main/images/CBNplot_readme_2.png?raw=true" width="800px">

- Another customized plot.
<img src="https://github.com/noriakis/software/blob/main/images/CBNplot_readme_3.png?raw=true" width="800px">
