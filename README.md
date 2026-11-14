# VIPER

VIPER measures the predictive importance of a compositional component by setting that component to zero, renormalizing the remaining components, and quantifying the resulting change in prediction. This repository provides the estimation functions, simulation scripts, and code for the YuJ_2015 colorectal cancer analysis.

## Repository contents

| File or directory | Description |
| --- | --- |
| [`R2.R`](R2.R) | Shared functions for compositional reparameterization, cross-fitted VIPER estimation and its variance estimation. |
| [`VIPER_partial_linear_root_n.R`](VIPER_partial_linear_root_n.R) | Simulation of the scaled bias and variance of the VIPER estimator across sample sizes. |
| [`Theorem1.R`](Theorem1.R) | Simulation comparing VIPER, log-contrast plug-in, and random forest plug-in estimators in terms of point estimation, confidence interval coverage, and rejection probabilities for a positive importance threshold. |
| [`YuJ_2015_CRC_VIPER_LC_RF.R`](YuJ_2015_CRC_VIPER_LC_RF.R) | Application of the three methods to stool microbiome data from colorectal cancer cases and controls in the YuJ_2015 study. |
| [`figures/`](figures/) | Reference figures for the simulations and real-data analysis. |

## Installation

Install the CRAN packages from an R session.

```r
install.packages(c("MASS", "randomForest", "doParallel", "foreach", "ggplot2"))
```

For the YuJ_2015 analysis, also install the Bioconductor packages.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("curatedMetagenomicData", "SummarizedExperiment"))
```

## Reproducing the figures

Open a terminal in the repository root, where `R2.R` is located. Run the command corresponding to the analysis you want to reproduce.

```bash
Rscript VIPER_partial_linear_root_n.R
Rscript Theorem1.R
Rscript YuJ_2015_CRC_VIPER_LC_RF.R
```

Each script loads `R2.R` and writes its PNG output to the working directory. The reference figures are available in [`figures/`](figures/).

| Analysis | Generated PNG files |
| --- | --- |
| Root-\(n\) simulation | `VIPER_partial_linear_root_n.png` |
| Estimation, coverage, and threshold-test simulation | `VIPER_two_models_estimation_coverage_n=100_q=0.2_c=0.010.png`<br>`VIPER_threshold_rejection_n=100_q=0.2_c=0.010.png` |
| YuJ_2015 real-data analysis | `YuJ_2015_CRC_species_VIPER_LC_RF_comparison.png` |

The default settings run the full Monte Carlo simulations and bootstrap analyses. Set `n_cores` in either simulation script or `ncores` in the real-data script to match the resources available on your machine.

## Data availability

The real-data script retrieves the `2021-03-31` YuJ_2015 relative-abundance resource through Bioconductor's `curatedMetagenomicData` package. ExperimentHub caches the downloaded resource locally. The first run requires an internet connection. The repository does not distribute the input data.
