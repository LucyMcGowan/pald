#' Tissue Data Distance Matrix
#'
#' A [`dist`] object describing distances from a subset of tissue gene
#' expression data from the following papers:
#'  * http://www.ncbi.nlm.nih.gov/pubmed/17906632
#'  * http://www.ncbi.nlm.nih.gov/pubmed/21177656
#'  * http://www.ncbi.nlm.nih.gov/pubmed/24271388
#' obtained from the **tissuesGeneExpression** bioconductor package.
#'
#' The original data frame had 189 rows,  each with a corresponding tissue,
#' such as `colon`, `kidney` or `cerebellum`.
#' There were 22,215 columns corresponding to gene expression data from each of
#' these rows. This was then converted into a distance matrix.
#'
#' @references M. Love and R. Irizarry. tissueGeneExpression. Bioconductor
#' Package
#'
#' @format A [`dist`] object of 189 tissue types
#'
"tissue_dist"
