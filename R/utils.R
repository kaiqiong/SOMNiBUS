#' @title hessianComp computes the Hessian matrix
#'
#' @description computes the Hessian matrix for an EM estimate
#' @param w_ij the diagnal values of the weight matrix
#' @param new.par estimate of alpha
#' @param new.lambda estimate of lambda
#' @param X a vector of read depths
#' @param Y a vector of methylated counts
#' @param my.design.matrix design matrix from the final fit
#' @param gam_smoothMat the smooth matrix from the final gam fit
#' @param Z covariate matrix
#' @param pred.pi predicted methylation probability from the final fit
#' @param p0 the probability of observing a methylated read when
#' the underlying true status is unmethylated. \code{p0} is the rate of false
#' methylation calls, i.e. false positive rate.
#' @param p1 the probability of observing a methylated read when
#' the underlying true status is methylated. \code{1-p1} is the rate of false
#' non-methylation calls, i.e. false negative rate.
#' @param disp_est estimated dispersion parameter
#' @param RanEff whether a subject-level Random effect is added or not
#' @param N number of unique samples in the provided dataset
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return a Hessian matrix for all alphas (including the alphas for RE)
#' @references Zhao, Kaiqiong, et al. 'A novel statistical method for
#' modeling  covariate effects in bisulfite sequencing derived measures
#' of DNA methylation.'Biometrics (2020).
#' @author  Kaiqiong Zhao
#' @importFrom Matrix bdiag
#' 
#' @noRd
hessianComp <- function(w_ij, new.par, new.lambda, X, Y, my.design.matrix,
                        gam_smoothMat, Z, pred.pi, p0, p1, disp_est, 
                        RanEff, N, verbose = TRUE) {
  
  t0 <- Sys.time()
  if(verbose) Message("Compute the Hessian matrix for an EM estimate")
  
  ## Q1: the second partial derivative w.r.t alpha^2 Q2: the second
  ## derivative w.r.t alpha & alpha_star
  res <- outer(seq_len(length(new.par)), 
               seq_len(length(new.par)),
               Vectorize(function(l,m) {
                 sum(-X * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
               }))
  smoth.mat <- lapply(as.list(seq_len(ncol(Z) + 1)), function(i) {
    gam_smoothMat[[i]]$S[[1]] * new.lambda[i]
  })  ## extract the penalty matrix
  ## assume the lambda for the constant of the intercept is 0 -- no
  ## penalization
  smoth.mat[[length(smoth.mat) + 1]] <- 0
  if (RanEff) {
    ## !!!! Otherwise, we get very wide CI
    smoth.mat[[length(smoth.mat) + 1]] <- diag(N) * new.lambda[ncol(Z) + 2]
    span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat) - 1,
                                                             (seq_len((length(smoth.mat) - 2))),
                                                             length(smoth.mat))]))
  } else {
    span.penal.matrix <- as.matrix(Matrix::bdiag(smoth.mat[c(length(smoth.mat),
                                                             (seq_len((length(smoth.mat) - 1))))]))
  }
  
  Q1_with_lambda <- res - span.penal.matrix/disp_est
  Q1_no_lambda <- res
  
  Q2 <- outer(seq_len(length(new.par)), seq_len(length(new.par)), 
              Vectorize(function(l, m) {
                term1 <- Y * p1 * p0/(p1 * pred.pi + p0 * (1 - pred.pi))^2 + (X - Y) * 
                  (1 - p1) * (1 - p0)/((1 - p1) * pred.pi + (1 - p0) * (1 - pred.pi))^2
                sum(term1 * w_ij * my.design.matrix[, m] * my.design.matrix[, l])
              }))
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(Q1_with_lambda + Q2)
}

#' @title Split methylation data into regions based on the spacing of CpGs
#' @description This function splits the methylation data into regions 
#' based on the spacing of CpGs.
#' 
#' @param dat a data frame with rows as individual CpGs appearing in all the
#' samples. The first 4 columns should contain the information of `Meth_Counts`
#' (methylated counts), `Total_Counts` (read depths), `Position` (Genomic 
#' position for the CpG site) and `ID` (sample ID). The covariate information, 
#' such as disease status or cell type composition, are listed in column 5 and 
#' onwards.
#' @param gap positive integer defining the gap width 
#' beyond which we consider that two regions are independent. 
#' Odd and decimal values will be rounded to the next even numbers 
#' (e.g. 8.2 and 8.7 become gaps of 8 and 10 respectively). 
#' The default value is \code{1e+6} (1Mb).
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return A named \code{list} of \code{data.frame} containing the data of each 
#' independent region.
#' @author  Audrey Lemaçon
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' results <- splitDataByRegion( dat=RAdat.f, gap = 1e6, min.cpgs = 5, 
#' verbose = FALSE)
#' @importFrom IRanges IRanges reduce gaps start end
#' @export
splitDataByRegion <- function(dat, gap = 1e6, min.cpgs = 50, max.cpgs = 2000,
                              verbose = TRUE){
  
  t0 <- Sys.time()
  msg <- paste("Partitioning the regions based on",
               "the spacing of CpGs")
  if(verbose) Message(msg)
  
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.")
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.")
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.")
    Warning(w_msg)
  }
  
  if(is.null(gap)){ 
    gap <- 1e+6 
  } else {
    if(!is.numeric(gap) || gap < 1){
      gap <- 1e+6
      w_msg <- paste("'gap' should be a positive integer.",
                     "We will use the default value 1e+6 (1Mb).") 
      Warning(w_msg)
    } else {
      ori_gap <- gap
      gap <- round(gap)
      if(gap %% 2 != 0) {
        gap <- gap + 1
      }
    }
  }
  
  # convert dataframe in IRanges
  dat_gr <- IRanges::IRanges(start = dat$Position, end = dat$Position)
  
  # reduce first orders the ranges in x from left to right, 
  # then merges the overlapping or adjacent ones.
  dat_gr <- IRanges::reduce(dat_gr)
  
  # identify gap between regions
  g = IRanges::gaps(dat_gr)
  
  # extend each region with half of the maximal distance d 
  # and attempt reduction to identify overlapping regions
  hemi_d <- ceiling(gap/2)
  regions_gr <- IRanges::reduce(dat_gr + hemi_d) - hemi_d
  
  # save the identified regions
  listOfRegions <- list()
  
  for(i in seq_along(regions_gr)){
    region = regions_gr[i]
    region_dat <- dat[dat$Position >= IRanges::start(region) & 
                        dat$Position <= IRanges::end(region),]
    
    if(length(unique(region_dat$Position)) > max.cpgs){
      msg <- paste0("Dropped region [",IRanges::start(region),"-",
                    IRanges::end(region),"] containing ",
                    length(unique(region_dat$Position)),
                    " measurement points: above maximal ",
                    "region size accepted for analysis (",max.cpgs,").")
      if(verbose) Message(msg)
    } else {
      if(length(unique(region_dat$Position)) < min.cpgs) {
        msg <- paste0("Dropped region [",IRanges::start(region),"-",
                      IRanges::end(region),"] containing ",
                      length(unique(region_dat$Position)),
                      " measurement points: below minimal ",
                      "region size accepted for analysis (",min.cpgs,").")
        if(verbose) Message(msg)
      } else{
        listOfRegions[[length(listOfRegions) + 1]] <- region_dat
      }
    }
  }
  
  # name each independent region
  names(listOfRegions) <- sapply(X = listOfRegions,
                                 FUN = function(region)
                                   paste("region",
                                         min(region$Position),
                                         max(region$Position), sep = "_")
  )
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  # return the list of independent regions
  return(listOfRegions)
}

#' @title Split methylation data into regions based on the density of CpGs
#' @description This function splits the methylation data into regions 
#' based on the density of CpGs.
#' 
#' @param dat a data frame with rows as individual CpGs appearing
#' in all the samples. The first 4 columns should contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths), 
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID). 
#' The covariate information, such as disease status or cell type composition,
#' are listed in column 5 and onwards.
#' @param window.size this positive integer defines the size of the
#' sliding window in bp. Decimal values will be rounded to the nearest integer.
#' The value should be greater than 10. The default value is \code{100} (100 bp)
#' @param by positive integer defines by how many base pairs the
#' window moves at each increment. Decimal values will be rounded to the
#' nearest integer. The default value is \code{1} (1 bp).
#' @param min.density positive integer defines the minimum density
#' threshold for each window. Decimal values will be rounded to the
#' nearest integer. The default value is \code{5} (5 CpGs/\code{window.size}).
#' @param gap positive integer defining the gap width
#' beyond which we consider that two regions are independent.
#' Decimal values will be rounded to the nearest integer.
#' The default value is \code{10} (10bp).
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information. The default value is TRUE.
#' @return A named \code{list} of \code{data.frame} containing the data of each
#' independent region.
#' @author  Audrey Lemaçon
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' results <- splitDataByDensity(dat = RAdat.f, window.size = 100, by = 1, 
#' min.density = 5, gap = 10, min.cpgs = 50, verbose = FALSE)
#' @importFrom IRanges IRanges reduce start end shift
#' 
#' @export
splitDataByDensity <- function(dat, window.size = 100, by = 1, min.density = 5, 
                               gap = 10, min.cpgs = 50, max.cpgs = 2000, 
                               verbose = TRUE){
  
  t0 <- Sys.time()
  msg <- paste0("Partitioning the regions based on ",
                "the density of CpGs")
  if(verbose) Message(msg)
  
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.") 
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.") 
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.")
    Warning(w_msg)
  }
  
  if(is.null(window.size)){
    window.size <- 100
  } else {
    if(!is.numeric(window.size) || window.size < 10) {
      window.size <- 100
      w_msg <- paste("'window.size' should be a positive integer",
                     "greater than 10. We will use the default value 100.") 
      Warning(w_msg)
    } else {
      window.size <- round(window.size)
    }
  }
  
  if(is.null(by)){
    by <- 1
  } else {
    if(!is.numeric(by) || by < 1) {
      by <- 1
      w_msg <- paste("'by' should be a positive integer.",
                     "We will use the default value 1.") 
      Warning(w_msg)
    } else {
      by <- round(by)
    }
  }
  
  if(is.null(min.density)){
    min.density <- 5
  } else {
    if(!is.numeric(min.density) || min.density < 1) {
      min.density <- 5
      w_msg <- paste("'min.density' should be a positive integer.",
                     "We will use the default value 5.") 
      Warning(w_msg)
    } else {
      min.density <- round(min.density)
    }
  }
  
  if(is.null(gap)){
    gap <- 10
  } else {
    if(!is.numeric(gap) || gap < 1) {
      gap <- 10
      w_msg <- paste("'gap' should be a positive integer.",
                     "We will use the default value 10.") 
      Warning(w_msg)
    } else {
      gap <- round(gap)
    }
  }
  
  if(gap >= window.size){
    gap <- round(window.size/2)
    w_msg <- paste0("'gap' should be lesser than the 'window.size' value. ",
                    "We will use half of window size (",gap,").")
    Warning(w_msg)
  }
  
  # convert dataframe in IRanges
  dat_gr <- unique(IRanges::IRanges(start = dat$Position, end = dat$Position))
  
  # save the end index for the while loop
  endOfRegion <- max(IRanges::end(dat_gr))
  
  # create the starting sliding window
  sliding_window <- IRanges::IRanges(start = min(IRanges::start(dat_gr)), 
                                     end = (min(IRanges::start(dat_gr)) + 
                                              window.size - 1))
  cpgs <- NULL
  
  # calculate the density in each sliding window
  while(IRanges::start(sliding_window) <= endOfRegion){
    
    # overlap region to window
    hits <- IRanges::findOverlaps(query = dat_gr, subject = sliding_window)
    
    # save density
    if(is.null(x = cpgs)){
      cpgs <- data.frame(start = IRanges::start(sliding_window), 
                         end = IRanges::end(sliding_window), 
                         density = length(hits))
    } else {
      cpgs <- rbind(cpgs, 
                    data.frame(start = IRanges::start(sliding_window), 
                               end = IRanges::end(sliding_window), 
                               density = length(hits)))
    }
    
    # shift window
    sliding_window <- IRanges::shift(x = sliding_window, shift = by)
  }
  
  # filter on density
  cpgs <- cpgs[cpgs$density >= min.density,]
  
  # transform the windows in IRanges to create regions
  regions_gr <- IRanges::IRanges(start = cpgs$start, end = cpgs$end)
  # each window is resized to identify independent regions
  regions_gr <- IRanges::resize(x = regions_gr, width = (gap + 1))
  regions_gr <- IRanges::reduce(regions_gr)
  
  # save the identified regions
  listOfRegions <- list()
  
  for(i in seq_along(regions_gr)){
    region = regions_gr[i]
    region_dat <- dat[dat$Position >= IRanges::start(region) & 
                        dat$Position <= IRanges::end(region),]
    
    if(length(unique(region_dat$Position)) > max.cpgs){
      msg <- paste0("Dropped region [",IRanges::start(region),"-",
                    IRanges::end(region),"] containing ",
                    length(unique(region_dat$Position)),
                    " measurement points: above maximal ",
                    "region size accepted for analysis (",max.cpgs,").")
      if(verbose) Message(msg)
    } else {
      if(length(unique(region_dat$Position)) < min.cpgs) {
        msg <- paste0("Dropped region [",IRanges::start(region),"-",
                      IRanges::end(region),"] containing ",
                      length(unique(region_dat$Position)),
                      " measurement points: below minimal ",
                      "region size accepted for analysis (",min.cpgs,").")
        if(verbose) Message(msg)
      } else{
        listOfRegions[[length(listOfRegions) + 1]] <- region_dat
      }
    }
  }
  
  # name each independent region
  names(listOfRegions) <- sapply(X = listOfRegions,
                                 FUN = function(region)
                                   paste("region",
                                         min(region$Position),
                                         max(region$Position), sep = "_")
  )
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(listOfRegions)
}


#' @title Split methylation data into regions based on the genomic annotations
#' @description This function splits the methylation data into regions 
#' based on the genomic annotations provided under the form of a GenomicRanges
#' object.
#' 
#' @param dat a data frame with rows as individual CpGs appearing
#' in all the samples. The first 4 columns should contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths), 
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID).
#' The covariate information, such as disease status
#' or cell type composition, are listed in column 5 and onwards.
#' @param chr character vector containing the chromosome information. Its length 
#' should be equal to the number of rows in \code{dat}.
#' @param annots GenomicRanges object containing the annotations
#' @param gap integer defining the maximum gap that is allowed between 
#' two regions to be considered as overlapping. 
#' According to the \code{GenomicRanges::findOverlaps} function, 
#' the gap between 2 ranges is the number of positions that separate them. 
#' The gap between 2 adjacent ranges is 0. By convention when one range has 
#' its start or end strictly inside the other (i.e. non-disjoint ranges), 
#' the gap is considered to be -1.
#' Decimal values will be rounded to the nearest integer. 
#' The default value is \code{-1} (meaning strict overlaping).
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return A named \code{list} of \code{data.frame} containing the data of each 
#' independent region.
#' @author  Audrey Lemaçon
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' annot <- GenomicRanges::GRanges(seqnames = "chr1", IRanges::IRanges(
#' start = c(102711720,102711844,102712006,102712503,102712702), 
#' end = c(102711757,102711909,102712195,102712637,102712712)
#' ))
#' results <- splitDataByGRanges(dat = RAdat.f, 
#' chr = rep(x = "chr1", times = nrow(RAdat.f)), 
#' annots = annot, gap = -1, min.cpgs = 5)
#' @importFrom IRanges IRanges reduce gaps start end shift
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevelsStyle seqinfo
#' @importFrom S4Vectors subjectHits queryHits
#' 
#' @export
splitDataByGRanges <- function(dat, chr, annots, gap = -1, 
                               min.cpgs = 50, max.cpgs = 2000, verbose = TRUE){
  
  t0 <- Sys.time()
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.")
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.")
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.")
    Warning(w_msg)
  }
  
  if(is.null(chr) || length(chr) != nrow(dat)) {
    e_msg <- paste("'chr' should be a character vector which length equals the",
                   "number of rows in 'dat'.") 
    Error(e_msg)
  }
  
  if(!inherits(annots, "GenomicRanges"))
    Error("'annots' should inherit the GRanges class.")
  
  if(is.null(gap)){
    gap <- -1
  } else {
    if(!is.numeric(gap) || gap < -1) {
      gap <- -1
      w_msg <- paste("'gap' should be an integer >= -1.",
                     "We will use the default value -1.") 
      Warning(w_msg)
    } else {
      gap <- round(gap)
    }
  }
  
  # convert dataframe in GRanges
  dat$Chr <- chr
  dat_gr <- GenomicRanges::GRanges(seqnames = dat$Chr, 
                                   ranges = IRanges::IRanges(start = dat$Position, 
                                                             end = dat$Position))
  dat_gr <- unique(dat_gr)
  
  # error handling
  dat_style <- tryCatch({
    # check if style has a compatible entry for the species supported by Seqname
    GenomeInfoDb::seqlevelsStyle(dat_gr)
  }, error=function(e) e)
  
  if(inherits(dat_style, "error")){
    stop_msg <- paste0("Incompatible chromosome name. ",
                       "Please see genomeStyles() for supported styles.")
    Error(stop_msg)
  }
  
  # error handling
  hits <- tryCatch({
    # adapt the style (chromosome nomenclature) of the annotation 
    # to the style of the input data if needed
    if(length(intersect(x = GenomeInfoDb::seqlevelsStyle(annots), 
                        y = dat_style)) == 0)
      suppressWarnings(GenomeInfoDb::seqlevelsStyle(annots) <- dat_style[1])
    
    # find overlaps
    GenomicRanges::findOverlaps(query = dat_gr, subject = annots, maxgap = gap)
  }, error=function(e) e)
  
  
  if(inherits(hits, "error")){
    stop_msg <- paste0("Incompatible formats between ",
                       "the input data and the annotations.")
    Error(stop_msg)
  }
  
  cpgs <- NULL
  
  if(length(hits) > 0){
    for(anno in unique(S4Vectors::subjectHits(hits))){
      hit <- hits[S4Vectors::subjectHits(hits) == anno]
      
      # save cpgs
      if(is.null(x = cpgs)){
        cpgs <- data.frame(chr = GenomeInfoDb::seqlevelsInUse(dat_gr[min(S4Vectors::queryHits(hit))]),
                           start = min(IRanges::start(dat_gr[S4Vectors::queryHits(hit)])), 
                           end = max(IRanges::end(dat_gr[S4Vectors::queryHits(hit)])), 
                           annotation = as.character(annots[anno]))
      } else {
        cpgs <- rbind(cpgs, 
                      data.frame(chr = GenomeInfoDb::seqlevelsInUse(dat_gr[min(S4Vectors::queryHits(hit))]),
                                 start = min(IRanges::start(dat_gr[S4Vectors::queryHits(hit)])), 
                                 end = max(IRanges::end(dat_gr[S4Vectors::queryHits(hit)])), 
                                 annotation = as.character(annots[anno])))
      }
    }
  }
  
  # transform the windows in GRanges to create regions
  regions_gr <- GenomicRanges::GRanges(seqnames = cpgs$chr, 
                                       IRanges::IRanges(start = cpgs$start, 
                                                        end = cpgs$end))
  
  # save the identified regions
  listOfRegions <- list()
  listOfRegions_names <- c()
  
  for(i in seq_along(regions_gr)){
    region = regions_gr[i]
    region_dat <- dat[(dat$Chr == GenomeInfoDb::seqlevelsInUse(region) & 
                         dat$Position >= IRanges::start(region) & 
                         dat$Position <= IRanges::end(region)),]
    
    if(length(unique(region_dat$Position)) > max.cpgs){
      msg <- paste0("Dropped region [",as.character(region),
                    "] containing ",length(unique(region_dat$Position)),
                    " measurement points: above maximal ",
                    "region size accepted for analysis (",max.cpgs,").")
      if(verbose) Message(msg)
    } else {
      if(length(unique(region_dat$Position)) < min.cpgs) {
        msg <- paste0("Dropped region [",as.character(region),
                      "] containing ",length(unique(region_dat$Position)),
                      " measurement points: below minimal ",
                      "region size accepted for analysis (",min.cpgs,").")
        if(verbose) Message(msg)
      } else {
        chrom <- gsub(x = unique(region_dat$Chr),
                      pattern = "chr", 
                      replacement = "")
        
        listOfRegions_names <- c(listOfRegions_names, 
                                 paste(paste0("chr",chrom),
                                       min(region_dat$Position),
                                       max(region_dat$Position), sep = "_"))
        # suppress the temporary 'Chr' column 
        region_dat$Chr <- NULL
        listOfRegions[[length(listOfRegions) + 1]] <- region_dat
      }
    }
  }
  
  # name each independent region
  names(listOfRegions) <- listOfRegions_names
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(listOfRegions)
}

#' @title Split methylation data into regions based on the genomic annotations
#' @description This function splits the methylation data into regions 
#' based on the genomic annotation provided under the form of a 1-based BED file
#' 
#' @param dat a data frame with rows as individual CpGs appearing
#' in all the samples. The first 4 columns should contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths), 
#' `Position` (Genomic position for the CpG site) and `ID` (sample ID). 
#' The covariate information, such as disease status
#' or cell type composition, are listed in column 5 and onwards.
#' @param chr character vector containing the chromosome information. Its length 
#' should be equal to the number of rows in \code{dat}.
#' @param bed character, path to the 1-based BED file containing the annotations
#' @param gap integer defining the maximum gap that is allowed between 
#' two regions to be considered as overlapping. 
#' According to the \code{GenomicRanges::findOverlaps} function, 
#' the gap between 2 ranges is the number of positions that separate them. 
#' The gap between 2 adjacent ranges is 0. By convention when one range has 
#' its start or end strictly inside the other (i.e. non-disjoint ranges), 
#' the gap is considered to be -1.
#' Decimal values will be rounded to the nearest integer. 
#' The default value is \code{-1} (meaning strict overlapping).
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return A named \code{list} of \code{data.frame} containing the data of each 
#' independent region.
#' @author  Audrey Lemaçon
#' 
#' @importFrom IRanges IRanges reduce gaps start end shift
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevelsStyle seqinfo
#' @importFrom S4Vectors subjectHits queryHits
#' 
#' @export
splitDataByBed <- function(dat, chr, bed, gap = -1, 
                           min.cpgs = 50, max.cpgs = 2000, verbose = TRUE){
  
  t0 <- Sys.time()
  msg <- paste0("Partitioning the regions based on ",
                "genomic annotations")
  if(verbose) Message(msg)
  
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.")
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.") 
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.") 
    Warning(w_msg)
  }
  
  if(is.null(chr) || length(chr) != nrow(dat)) {
    e_msg <- paste("'chr' should be a character vector which length equals the",
                   "number of rows in 'dat'.") 
    Error(e_msg)
  }
  
  if(is.null(bed) || !is.character(bed) || !file.exists(bed))
    Error("'bed' should a valid file path.")
  
  # error handling in BED import
  annots <- tryCatch(
    rtracklayer::import(bed),
    error=function(e) e
  )
  
  if(inherits(annots, "error"))
    Error("'bed' should a path to a valid BED file.")
  
  # the BED file is converted automatically by rtracklayer in 0-based,
  # we have to return back to 1-based
  IRanges::start(annots) <- IRanges::start(annots) - 1
  
  if(is.null(gap)){
    gap <- -1
  } else {
    if(!is.numeric(gap) || gap < -1) {
      gap <- -1
      w_msg <- paste("'gap' should be an integer >= -1.",
                     "We will use the default value -1.") 
      Warning(w_msg)
    } else {
      gap <- round(gap)
    }
  }
  
  listOfRegions <- splitDataByGRanges(dat = dat, chr = chr, annots = annots, 
                                      gap = gap, min.cpgs = min.cpgs, 
                                      max.cpgs = max.cpgs, verbose = verbose)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(listOfRegions)
}

#' @title Split methylation data into regions based on the chromatin states 
#' @description This function splits the methylation data into regions 
#' based on the chromatin states predicted by ChromHMM software 
#' (Ernst and Kellis (2012)). 
#' The annotations come from the Bioconductor package `annnotatr`. 
#' Chromatin states determined by chromHMM are 
#' available in hg19 for nine cell lines (Gm12878, H1hesc, Hepg2, Hmec, Hsmm, 
#' Huvec, K562, Nhek, and Nhlf). 
#' 
#' @param dat a data frame with rows as individual CpGs appearing
#' in all the samples. The first 4 columns should contain the information of
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths), 
#' `Position` (Genomic position for the CpG site) and `ID`(sample ID). 
#' The covariate information, such as disease status
#' or cell type composition, are listed in column 5 and onwards.
#' @param chr character vector containing the chromosome information. Its length 
#' should be equal to the number of rows in \code{dat}.
#' @param cell.line character defining the cell line of interest. Nine cell 
#' lines are available: 
#' - \code{"gm12878"}: Lymphoblastoid cells GM12878, 
#' - \code{"h1hesc"}: Embryonic cells H1 hESC, 
#' - \code{"hepg2"}: Liver carcinoma HepG2,  
#' - \code{"hmec"}, Mammary epithelial cells HMEC, 
#' - \code{"hsmm"}, Skeletal muscle myoblasts HSMM, 
#' - \code{"huvec"}: Umbilical vein endothelial HUVEC, 
#' - \code{"k562"}: Myelogenous leukemia K562, 
#' - \code{"nhek"}: Keratinocytes NHEK,
#' - \code{"nhlf"}: Normal human lung fibroblasts NHLF.
#' @param states character vector defining the chromatin states of interest
#' among the following available options:
#' - \code{"ActivePromoter"}: Active Promoter
#' - \code{"WeakPromoter"}: Weak Promoter
#' - \code{"PoisedPromoter"}: Poised Promoter
#' - \code{"StrongEnhancer"}: Strong Enhancer
#' - \code{"WeakEnhancer"}: Weak/poised Enhancer
#' - \code{"Insulator"}: Insulator
#' - \code{"TxnTransition"}: Transcriptional Transition
#' - \code{"TxnElongation"}: Transcriptional Elongation
#' - \code{"WeakTxn"}: Weak Transcribed
#' - \code{"Repressed"}: Polycomb-Repressed
#' - \code{"Heterochrom"}: Heterochromatin; low signal
#' - \code{"RepetitiveCNV"}: Repetitive/Copy Number Variation
#' Use \code{state="all"} to select all the states simultaneously.
#' @param gap this integer defines the maximum gap that is allowed between 
#' two regions to be considered as overlapping. 
#' According to the \code{GenomicRanges::findOverlaps} function, 
#' the gap between 2 ranges is the number of positions that separate them. 
#' The gap between 2 adjacent ranges is 0. By convention when one range has 
#' its start or end strictly inside the other (i.e. non-disjoint ranges), 
#' the gap is considered to be -1.
#' Decimal values will be rounded to the nearest integer. 
#' The default value is \code{-1}.
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return A \code{list} of \code{data.frame} containing the data of each 
#' independent region.
#' @author  Audrey Lemaçon
#' 
#' @importFrom yaml yaml.load_file
#' @importFrom BiocManager install
#' @importFrom annotatr build_annotations
#' @importFrom utils installed.packages
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' results <- splitDataByChromatin(dat = RAdat.f, 
#' cell.line = "huvec", chr = rep(x = "chr4", times = nrow(RAdat.f)),
#' states = "Insulator", verbose = FALSE)
#' 
#' @export
splitDataByChromatin <- function(dat, chr, cell.line, states, 
                                 gap = -1, min.cpgs = 50, max.cpgs = 2000,
                                 verbose = TRUE){
  
  t0 <- Sys.time()
  msg <- paste0("Partitioning the regions based on ",
                "chromatin states")
  if(verbose) Message(msg)
  
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.") 
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.")
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.")
    Warning(w_msg)
  }
  
  if(is.null(chr) || length(chr) != nrow(dat)) {
    e_msg <- paste("'chr' should be a character vector which length equals the",
                   "number of rows in 'dat'.") 
    Error(e_msg)
  }
  
  # load annotation conf
  annots_conf <- yaml.load_file(input = system.file("annotations.yml", 
                                                    package = 'SOMNiBUS'))
  
  # set organism and build
  organism <- "human"
  build <- "hg19"
  
  # check the cell line
  chromatin_annots <- annots_conf[[organism]][[build]][["chromatin"]]
  supported_types <- names(chromatin_annots)
  
  if(is.null(cell.line) || !is.character(cell.line)) {
    stop_msg <- paste0("'cell.line' argument is mandatory. Choose one among ",
                       "these supported cell lines: '",
                       paste(supported_types, collapse = "', '"),"'.")
    Error(stop_msg)
  } else {
    cell.line <- tolower(cell.line)
  }
  
  if(!cell.line %in% supported_types){
    stop_msg <- paste0("'",cell.line,"' cell line is not supported. ",
                       "Choose one among these supported cell lines: '",
                       paste(supported_types, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  # check the chromHMM states
  supported_states <- names(chromatin_annots[[cell.line]])
  
  if(is.null(states) || !is.character(states)) {
    stop_msg <- paste0("'states' argument is mandatory. Choose among ",
                       "these supported states: '",
                       paste(supported_states, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  if(sum(states %in% supported_states) == 0){
    stop_msg <- paste0("None of the requested states are supported. ",
                       "Choose among these supported states: '",
                       paste(supported_states, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  unsupported_states <- unique(states[!states %in% supported_states])
  if(length(unsupported_states) > 0){
    w_msg <- paste0("The following unsupported states will be ignored: '",
                    paste(unsupported_states, collapse = "', '"),"'. ",
                    "Choose among these supported states: '",
                    paste(supported_states, collapse = "', '"),"'.")
    Warning(w_msg)
  }
  
  if(is.null(gap)){
    gap <- -1
  } else {
    if(!is.numeric(gap) || gap < -1) {
      gap <- -1
      w_msg <- paste("'gap' should be an integer >= -1.",
                     "We will use the default value -1.")
      Warning(w_msg)
    } else {
      gap <- round(gap)
    }
  }
  
  # retrieve annotations
  states <- unique(states[!states %in% unsupported_states])
  
  listOfRegions <- lapply(X = states, FUN = function(state){
    annot <- chromatin_annots[[cell.line]][[state]]
    annot <- suppressMessages(annotatr::build_annotations(genome = build, 
                                                          annotations = annot))
    # suppress duplicates
    annot <- unique(annot)
    
    # findoverlap
    annots <- splitDataByGRanges(dat = dat, chr = chr, 
                                 annots = annot, gap = gap, 
                                 min.cpgs = min.cpgs, max.cpgs = max.cpgs,
                                 verbose = verbose)
    # overwrite the name of each independent region
    if(length(annots) > 0) 
      names(annots) <- paste(names(annots), cell.line, state, sep = "_")
    
    return(annots)
  })
  
  # create final listOfRegions
  listOfRegions <- unlist(listOfRegions,recursive=FALSE)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(listOfRegions)
}

#' @title Split methylation data into regions based on the genes annotations
#' @description This function splits the methylation data into regions 
#' based on the genes. The annotations are coming from the Bioconductor
#' package `annnotatr`.
#' 
#' @param dat a data frame with rows as individual CpGs appearing
#' in all the samples. The first 4 columns should contain the information of 
#' `Meth_Counts` (methylated counts), `Total_Counts` (read depths), 
#' `Position` (Genomic position for the CpG site) and `ID`(sample ID). 
#' The covariate information, such as disease status or cell type composition, 
#' are listed in column 5 and onwards.
#' @param chr character vector containing the chromosome information. Its length 
#' should be equal to the number of rows in \code{dat}.
#' @param organism character defining the organism of interest 
#' Only Homo sapiens (\code{"human"}) is available.
#' Additional packages are required for Mus musculus (\code{"mouse"}),
#' Rattus norvegicus (\code{"rat"}) and Drosophila melanogaster (\code{"fly"}). 
#' The matching is case-insensitive. The default value is \code{"human"}.
#' @param build character defining the version of the genome build on which the 
#' methylation data have been mapped. By default, the build is set to 
#' \code{"hg38"}, however the build \code{"hg19"} is also available for
#' Homo sapiens: 
#' Once the additional packages are installed, the following organisms and 
#' builds are available: 
#' - \code{"mm9"} and \code{"mm10"} for Mus musculus; 
#' - \code{"rn4"}, \code{"rn5"} and \code{"rn6"} for Rattus norvegicus; 
#' - \code{"dm3"} and \code{"dm6"} for Drosophila melanogaster; 
#' @param types character vector defining the type of genic annotations
#' to use among the following options:
#' - \code{"upstream"} for the annotations included 1-5Kb upstream of the TSS; 
#' - \code{"promoter"} for the annotations included < 1Kb upstream of the TSS; 
#' - \code{"threeprime"} for the annotations included in 3' UTR;
#' - \code{"fiveprime"} for the annotations included in the 5' UTR;
#' - \code{"exon"} for the annotations included in the exons;
#' - \code{"intron"} for the annotations included in the introns;
#' - \code{"all"} for all the annotations aforementioned. 
#' The default value is \code{"promoter"}.
#' @param gap this integer defines the maximum gap allowed between two regions
#' to be considered as overlapping.
#' According to the \code{GenomicRanges::findOverlaps} function, 
#' the gap between 2 ranges is the number of positions that separate them. 
#' The gap between 2 adjacent ranges is 0. By convention when one range has 
#' its start or end strictly inside the other (i.e. non-disjoint ranges), 
#' the gap is considered to be -1.
#' Decimal values will be rounded to the nearest integer. 
#' The default value is \code{-1}.
#' @param min.cpgs positive integer defining the minimum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 50.
#' @param max.cpgs positive integer defining the maximum number of 
#' CpGs within a region for the algorithm to perform optimally. 
#' The default value is 2000.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information.
#' The default value is TRUE.
#' @return A named \code{list} of \code{data.frame} containing the data of each 
#' independent region.
#' @author  Audrey Lemaçon
#' 
#' @importFrom yaml yaml.load_file
#' @importFrom BiocManager install
#' @importFrom annotatr build_annotations
#' @importFrom utils installed.packages
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' # Add a column containing the chromosome information
#' RAdat$Chr <- "chr4"
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' results <- splitDataByGene(dat = RAdat.f, 
#' chr = rep(x = "chr1", times = nrow(RAdat.f)), verbose = FALSE)
#' 
#' @export
splitDataByGene <- function(dat, chr, organism = "human", 
                            build = "hg38", types = "promoter", gap = -1, 
                            min.cpgs = 50, max.cpgs = 2000, verbose = TRUE){
  
  t0 <- Sys.time()
  msg <- paste0("Partitioning the regions based on ",
                "genic annotations")
  if(verbose) Message(msg)
  
  # test the arguments
  if(is.null(min.cpgs)){
    min.cpgs <- 50
  } else {
    if(!is.numeric(min.cpgs) || min.cpgs < 1) {
      min.cpgs <- 50
      w_msg <- paste("'min.cpgs' should be a positive integer.",
                     "We will use the default value 50.")
      Warning(w_msg)
    } else {
      min.cpgs <- round(min.cpgs)
    }
  }
  
  if(is.null(max.cpgs)){
    max.cpgs <- 2000
  } else {
    if(!is.numeric(max.cpgs) || max.cpgs < 1) {
      max.cpgs <- 2000
      w_msg <- paste("'max.cpgs' should be a positive integer.",
                     "We will use the default value 2000.")
      Warning(w_msg)
    } else {
      max.cpgs <- round(max.cpgs)
    }
  }
  
  if(max.cpgs <= min.cpgs){
    min.cpgs <- 50
    max.cpgs <- 2000
    w_msg <- paste("'max.cpgs' should be a higher than 'min.cpgs'.",
                   "We will use the default values:",
                   "min.cpgs=50 and max.cpgs=2000.")
    Warning(w_msg)
  }
  
  if(is.null(gap)){
    gap <- -1
  } else {
    if(!is.numeric(gap) || gap < -1) {
      gap <- -1
      w_msg <- paste("'gap' should be an integer >= -1.",
                     "We will use the default value -1.")
      Warning(w_msg)
    } else {
      gap <- round(gap)
    }
  }
  
  if(is.null(chr) || length(chr) != nrow(dat)) {
    e_msg <- paste("'chr' should be a character vector which length equals the",
                   "number of rows in 'dat'.") 
    Error(e_msg)
  }
  
  # load annotation conf
  annots_conf <- yaml.load_file(input = system.file("annotations.yml", 
                                                    package = 'SOMNiBUS'))
  
  # check organism
  if(is.null(organism)){
    organism <- "human"
    build <- "hg38"
  } else{
    if(!is.character(organism)) {
      organism <- "human"
      build <- "hg38"
      w_msg <- paste("'organism' should be a non-null character.",
                     "We will use the default value human",
                     "and the build will be set at hg38.")
      Warning(w_msg)
    } else {
      organism <- tolower(organism)
    }
  }
  
  if(!organism %in% names(annots_conf)){
    stop_msg <- paste0("'",organism,"' organism is not supported. ",
                       "Choose one among these supported types: '",
                       paste(names(annots_conf), collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  # check build
  if(is.null(build)){
    build <- annots_conf[[organism]]$default
  } else {
    if(!is.character(build)) {
      build <- annots_conf[[organism]]$default
      warning_msg <- paste0("'build' should be a non-null character.",
                            "We will use the default", 
                            "value for ",organism, ": ", build,".")
      Warning(warning_msg)
    } else {
      build <- tolower(build)
    }
  }
  
  supported_builds <- names(annots_conf[[organism]])
  supported_builds <- supported_builds[!supported_builds %in% "default"]
  if(!build %in% supported_builds){
    stop_msg <- paste0("'",build,"' build is not supported. Choose one among these ",
                       "supported builds: '",
                       paste(supported_builds, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  # check the annotation types
  if(is.null(types)){
    types <- "promoter"
  } else {
    if(!is.character(types)) {
      types <- "promoter"
      w_msg <- paste("'type' should be a non-null character.",
                     "We will use the default value promoter.")
      Warning(w_msg)
    } else {
      types <- tolower(types)
    }
  }
  
  supported_types <- names(annots_conf[[organism]][[build]][["genes"]])
  supported_types <- supported_types[!supported_types %in% "dependencies"]
  
  if(sum(types %in% supported_types) == 0){
    stop_msg <- paste0("None of the requested types are supported. ",
                       "Choose among these supported types: '",
                       paste(supported_types, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  unsupported_types <- unique(types[!types %in% supported_types])
  if(length(unsupported_types) > 0){
    w_msg <- paste0("The following unsupported types will be ignored: '",
                    paste(unsupported_types, collapse = "', '"),"'. ",
                    "Choose among these supported types: '",
                    paste(supported_types, collapse = "', '"),"'.")
    Warning(w_msg)
  }
  
  # retrieve annotations
  types <- unique(types[!types %in% unsupported_types])
  
  # check if the annotation package is installed, if not install stop 
  annotation_packages <- annots_conf[[organism]][[build]]$genes$dependencies
  missing_packages <- annotation_packages[! annotation_packages %in% 
                                            installed.packages()[,"Package"]]
  if (length(missing_packages) > 0){
    stop_msg <- paste0("Please install the following package(s) in order to ",
                       "perform the requested partiniong: '",
                       paste(missing_packages, collapse = "', '"),"'.")
    Error(stop_msg)
  }
  
  listOfRegions <- lapply(X = types, FUN = function(type){
    # retrieve annotations
    annot <- annots_conf[[organism]][[build]]$genes[[type]]
    annot <- suppressMessages(annotatr::build_annotations(genome = build, 
                                                          annotations = annot))
    # suppress duplicates
    annot <- unique(annot)
    
    # findoverlap
    annots <- splitDataByGRanges(dat = dat, chr = chr, annots = annot, 
                                 gap = gap, min.cpgs = min.cpgs, 
                                 max.cpgs = max.cpgs, verbose = verbose)
    # overwrite the name of each independent region
    if(length(annots) > 0) 
      names(annots) <- paste(names(annots), type, sep = "_")
    
    return(annots)
  })
  
  # create final listOfRegions
  listOfRegions <- unlist(listOfRegions,recursive=FALSE)
  
  msg <- paste("Process completed in",
               format(Sys.time() - t0, digits = 2))
  if(verbose) Message(msg, step = "Finished")
  
  return(listOfRegions)
}

