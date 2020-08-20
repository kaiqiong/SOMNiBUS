#' A simulated methylation dataset based on a real data.
#'
#' This example data include methylation levels on a region with 208 CpGs for 116
#' blood samples.
#'
#' @format A data frame of 6064 rows and 13 columns. Each row represents a CpG site
#' for a sample. Columns include in order
#' \describe{
#'  \item{Meth_Counts}{Number of methylated reads}
#'  \item{Total_Counts}{Total number of reads; read-depth}
#'  \item{Position}{Genomic position (in bp) for the CpG site}
#'  \item{ID}{indicates which sample the CpG site belongs to}
#'  \item{ACPA4}{binary indicator for a biomarker anti-citrullinated protein antibody}
#'  \item{Age}{Age}
#'  \item{Sex}{2-female; 1-male}
#'  \item{Smoking}{1-current or ex-smoker; 0-non-smoker}
#'  \item{Smoking_NA}{1-Smoking info is NA; 0-Smoking info is available}
#'  \item{PC1}{PC1 for the cell type proportions}
#'  \item{PC2}{PC2 for the cell type proportions}
#'  \item{PC3}{PC3 for the cell type proportions}
#'  \item{PC4}{PC4 for the cell type proportions}
#' }
#' @source simulation is based a real data set provided by PI Dr. Sasha Bernatsky
#' (McGill University)
"RAdat2"
