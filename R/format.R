#' @title Parsing output from the BSseq package
#' @description This function reads and converts a BSseq object into a 
#' \code{list} of \code{data.frame}s (one per chromosome) to a format compatible
#' with SOMNiBUS' main functions \code{runSOMNiBUS} and 
#' \code{binomRegMethModel}.
#' @param bsseq_dat an object of class BSseq.
#' @param verbose logical indicates the level of information provided by the 
#' algorithm during the process. The default value is TRUE.
#' @return This function returns a \code{list} of \code{data.frame}s (one per 
#' chromosome). Each \code{data.frame} contains rows as individual CpGs 
#' appearing in all the samples. The first 4 columns contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths),
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID). 
#' The additional information (such as disease status, sex, age) extracted from
#' the BSseq object are listed in column 5 and onwards and will be considered 
#' as covariate information by SOMNiBUS algorithms. 
#' @author Audrey Lemaçon
#' @family Parsing functions
#' @seealso  \link[bsseq]{BSseq} for the BSseq class.
#' @examples
#' M <- matrix(1:9, 3,3)
#' colnames(M) <- c("A1", "A2", "A3")
#' BStest <- bsseq::BSseq(pos = 1:3, chr = c("chr1", "chr2", "chr1"), 
#' M = M, Cov = M + 2)
#' dat <- formatFromBSseq(BStest, verbose = FALSE)
#' @importFrom bsseq start seqnames getCoverage pData
#' @importFrom reshape2 melt
#'
#' @export
formatFromBSseq <- function(bsseq_dat, verbose = TRUE){
  
  t0 <- Sys.time()
  if(!inherits(x = bsseq_dat, what = "BSseq")){
    stop("'bsseq' should inherit the BSseq class. ?BSseq for more information.")
  }
  
  if(verbose) Message("Converting BSseq into SOMNiBUS")
  
  # extract positions
  chrom <- gsub(x = as.character(bsseq::seqnames(bsseq_dat)), pattern = "chr", 
                replacement = "")
  position <- paste0("chr",chrom,"_",bsseq::start(bsseq_dat))
  
  # (Total_Counts in SOMNiBUS)
  if(verbose) 
    Message("Extracting coverage information",step = NULL)
  
  bsseq_cov <- as.data.frame(bsseq::getCoverage(bsseq_dat))
  bsseq_cov$Position <- position
  bsseq_cov <- reshape2::melt(data = bsseq_cov, id = "Position", 
                              value.name = "Total_Counts", variable.name = "ID")
  
  # (Meth_Counts in SOMNiBUS)
  if(verbose)
    Message("Extracting methylation count", step = NULL) 
  
  bsseq_met <- as.data.frame(bsseq::getCoverage(bsseq_dat, type = "M"))
  bsseq_met$Position <- position
  bsseq_met <- reshape2::melt(data = bsseq_met, id = "Position", 
                              value.name = "Meth_Counts", variable.name = "ID")
  
  # merge counts
  bsseq_df <- base::merge(x = bsseq_cov, y = bsseq_met, by = c("Position","ID"))
  
  # retrieve phenotype data if exists
  pheno_data <- bsseq::pData(bsseq_dat)
  if(!is.null(pheno_data)){
    if(verbose)
      Message("Retrieving phenotype/experimental data", step = NULL)
    pheno_data <- as.data.frame(pheno_data)
    pheno_data$ID <- rownames(pheno_data)
    bsseq_df <- base::merge(x = bsseq_df, y = pheno_data, by = "ID")
  }
  
  # filter out "bad data"
  bad_dat_idx <- which(bsseq_df$Total_Counts == 0)
  if(length(bad_dat_idx) > 0){
    bsseq_df <- bsseq_df[-bad_dat_idx,]
    msg <- paste(" Remove",length(bad_dat_idx),
                 "rows with Total_Counts equal to 0.")
    if(verbose) message(msg)
  }
  
  # split by chromosome
  chr_position <- strsplit(x = bsseq_df$Position, split = "_")
  bsseq_df_chr <- sapply(X = chr_position, FUN = function(i) i[1])
  bsseq_df$Position <- sapply(X = chr_position, FUN = function(i) i[2])
  bsseq_df <- base::split(x = bsseq_df, f = bsseq_df_chr)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(bsseq_df)
  
}

#' @title Parsing output from the Bismark alignment suite
#' @description This function reads and converts Bismark's 
#' \href{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-genome-wide-cytosine-report-optional-is-tab-delimited-in-the-following-format-1-based-coords}{'genome wide cytosine report'} and \href{https://github.com/FelixKrueger/Bismark/tree/master/Docs#the-coverage-output-looks-like-this-tab-delimited-1-based-genomic-coords}{'coverage'} into a 
#' \code{list} of \code{data.frame}s (one per chromosome) to a format compatible
#' with SOMNiBUS' main functions \code{runSOMNiBUS} and 
#' \code{binomRegMethModel}.
#' @param ... parameters from \code{bsseq::read.bismark()} function
#' @param verbose logical indicates the level of information provided by the 
#' algorithm during the process. The default value is TRUE.
#' @return This function returns a \code{list} of \code{data.frame}s (one per 
#' chromosome). Each \code{data.frame} contains rows as individual CpGs 
#' appearing in all the samples. The first 4 columns contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths),
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID). 
#' The additional information (such as disease status, sex, age) extracted from
#' the BSseq object are listed in column 5 and onwards and will be considered 
#' as covariate information by SOMNiBUS algorithms. 
#' @author  Audrey Lemaçon
#' @seealso  \link[bsseq]{read.bismark} for parsing output from the Bismark
#' alignment suite.
#' @family Parsing functions
#' @examples
#' infile <- system.file("extdata/test_data.fastq_bismark.bismark.cov.gz",
#' package = "bsseq")
#' dat <- formatFromBismark(infile, verbose = FALSE)
#' @importFrom bsseq read.bismark
#'
#' @export
formatFromBismark <- function(..., verbose = TRUE){
  
  t0 <- Sys.time()
  if(verbose) Message("Converting Bismark into BSseq")
  
  bsseq_dat <- bsseq::read.bismark(..., verbose = verbose)
  bsseq_df <- formatFromBSseq(bsseq_dat = bsseq_dat, verbose = verbose)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(bsseq_df)
}

#' @title Parsing output from the BSmooth alignment suite
#' @description This function reads and converts BSmooth's output into a 
#' \code{list} of \code{data.frame}s (one per chromosome) to a format compatible
#' with SOMNiBUS' main functions \code{runSOMNiBUS} and 
#' \code{binomRegMethModel}.
#' @param ... parameters from \code{bsseq::read.bsmooth()} function
#' @param verbose logical indicates the level of information provided by the 
#' algorithm during the process. The default value is TRUE.
#' @return This function returns a \code{list} of \code{data.frame}s (one per 
#' chromosome). Each \code{data.frame} contains rows as individual CpGs 
#' appearing in all the samples. The first 4 columns contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths),
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID). 
#' The additional information (such as disease status, sex, age) extracted from
#' the BSseq object are listed in column 5 and onwards and will be considered 
#' as covariate information by SOMNiBUS algorithms. 
#' @author  Audrey Lemaçon
#' @seealso  \link[bsseq]{read.bismark} for parsing output from the Bismark
#' alignment suite.
#' @family Parsing functions
#' @examples
#' indir <- system.file("extdata/ev_bt2_tab", package = "SOMNiBUS")
#' dat <- formatFromBSmooth(indir, verbose = FALSE)
#' @importFrom bsseq read.bsmooth
#'
#' @export
formatFromBSmooth <- function(..., verbose = TRUE){
  
  t0 <- Sys.time()
  if(verbose) Message("Converting BSmooth into BSseq")
  
  bsseq_dat <- bsseq::read.bsmooth(..., verbose = verbose)
  bsseq_df <- formatFromBSseq(bsseq_dat = bsseq_dat, verbose = verbose)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(bsseq_df)
}
