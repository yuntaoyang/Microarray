# Microarray
* This is an example of downloading and analyzing GEO datasets using R.
* GSE31684: Male vs. Female.
## R scripts
* Check `microarray.R` .
* R requirements: `GEOquery, limma, ggplot2, dplyr`.
## Microarray raw data
* `./data/exp.csv`: the normalized probe expression data.
* `./data/probe.csv`: the probe information data.
* `./data/metadata.csv`: the metadata.
## Analysis output
* `./output/PCA.png`: the PCA plot of samples.
* `./output/differential_analysis.csv`: the result of differential gene expression analysis.
* `./output/gsea.rnk`: the preranked genes for the input of GSEA.
