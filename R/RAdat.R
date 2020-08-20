#' Methylation data from a rheumatoid arthritis study
#'
#' A dataset containing methylation levels on one targeted region on chromosome 4
#' near gene BANK1 from cases with rheumatoid arthritis (RA) and controls
#'
#' This example data include methylation levels of cell type separated blood samples
#' of 22 rheumatoid arthritis (RA) patients and 21 healthy individuals. In the data set,
#' 123 CpG sites are measured and there are 25 samples from circulating T cells and 18
#' samples from monocytes.
#'
#' @format A data frame of 5289 rows and 6 columns. Each row represents a CpG site
#' for a sample. Columns include in order
#' \describe{
#'  \item{Meth_Counts}{Number of methylated reads}
#'  \item{Total_Counts}{Total number of reads; read-depth}
#'  \item{Position}{Genomic position (in bp) for the CpG site}
#'  \item{ID}{indicates which sample the CpG site belongs to}
#'  \item{T_cell}{whether a sample is from T cell or monocyte}
#'  \item{RA}{whether a sample is an RA patient or control}
#' }
#' @source Dr. Marie Hudson (McGill University)
"RAdat"
