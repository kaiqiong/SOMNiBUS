#' @title Plot the smooth covariate effect
#'
#' @description This function accepts an output object from function
#' \code{binomRegMethModel} and print out a plot of the estimated
#' effect across the region for each test covariate.
#' @param BEM.obj an output object from function \code{binomRegMethModel}
#' @param mfrow A vector of the form c(nr, nc). Subsequent figures will 
#' be drawn in an nr-by-nc array on the device.
#' @param same.range specify whether the plots should be in the same
#' vertical scale
#' @param title the text for the title
#' @param covs a vector of covariate names. The covariates with
#' names in \code{covs} will be included in the plot.
#' When the value is set to NULL all the covariates and the Intercept will be
#' represented. The default value is NULL.
#' @param save file name to create on disk. When the value is set to NULL, 
#' the plot is not saved. The default value is NULL.
#' @param verbose logical indicates if the algorithm should provide progress
#' report information. The default value is TRUE.
#' @return This function prints out a plot of smooth covariate effects
#' and its pointwise confidence intervals
#' @author  Kaiqiong Zhao, Audrey Lemaçon
#' @importFrom data.table as.data.table
#' @importFrom tidyr pivot_longer %>%
#' @importFrom ggplot2 ggplot geom_line theme_bw aes_string theme
#' @importFrom ggplot2 geom_hline geom_rug xlab ggtitle facet_wrap ggsave
#' @examples
#' #------------------------------------------------------------#
#' data(RAdat)
#' head(RAdat)
#' RAdat.f <- na.omit(RAdat[RAdat$Total_Counts != 0, ])
#' out <- binomRegMethModel(
#'   data=RAdat.f, n.k=rep(5, 3), p0=0.003307034, p1=0.9,
#'   epsilon=10^(-6), epsilon.lambda=10^(-3), maxStep=200,
#'   Quasi = FALSE, RanEff = FALSE
#' )
#' binomRegMethModelPlot(out, same.range=FALSE)
#' 
#' 
#' @export
binomRegMethModelPlot <- function(BEM.obj, mfrow = NULL, same.range = FALSE,
                                  title = "Smooth covariate effects", 
                                  covs = NULL, save = NULL,
                                  verbose = TRUE) {
  
  if(verbose) Message("Plot the smooth covariate effect")
  
  # transform data.frame to data.table to ensure the format integrity 
  # even if there is a single cov selected
  BEM.obj$Beta.out <- data.table::as.data.table(BEM.obj$Beta.out)
  BEM.obj$SE.out <- data.table::as.data.table(BEM.obj$SE.out)
  
  if(is.null(mfrow)) mfrow = c(NULL,NULL)
  
  # select and re-order covariates of interest
  if(!is.null(covs)){
    # remove duplicate
    covs <- unique(covs)
    # raise an error if no covs matches the content of the results df
    selected_covs <- covs[covs %in% colnames(BEM.obj$Beta.out)]
    if(length(selected_covs) == 0){
      # some unknown covs
      e_msg <- paste("None of the requested covariates are included in the",
                     "BEM.obj object")
      Error(e_msg)
    }
    
    # check if the covs of interest are present in the results df
    unknown_covs <- covs[!covs %in% colnames(BEM.obj$Beta.out)]
    if(length(unknown_covs) > 0){
      # some unknown covs
      w_msg <- paste0("Some covariates are missing from the BEM.obj object: '",
                      paste(unknown_covs, collapse = ", "),"'")
      Warning(w_msg)
    }
    
    # To avoid being flagged by R CMD check ("no visible binding for global variable ‘..selected_covs’"), 
    # we used the alternative with = FALSE notation.
    # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
    BEM.obj$Beta.out <- BEM.obj$Beta.out[, selected_covs, with = FALSE]
    BEM.obj$SE.out <- BEM.obj$SE.out[, selected_covs, with = FALSE]
  } else {
    selected_covs <- colnames(BEM.obj$Beta.out)
  }
  
  # re-order the lines according to the genomic positions
  Position_order <- order(BEM.obj$uni.pos)
  BEM.obj$Beta.out <- BEM.obj$Beta.out[Position_order,]
  BEM.obj$SE.out <- BEM.obj$SE.out[Position_order,]
  BEM.obj$uni.pos <- BEM.obj$uni.pos[Position_order]
  
  ll <- BEM.obj$Beta.out - 1.96 * BEM.obj$SE.out
  hh <- BEM.obj$Beta.out + 1.96 * BEM.obj$SE.out
  
  estimates <- cbind(Position = BEM.obj$uni.pos, BEM.obj$Beta.out)
  estimates <- as.data.frame(estimates) %>% 
    tidyr::pivot_longer(!Position, 
                        names_to = "Covariate", 
                        values_to = "Estimate")
  
  lower_bounds <- cbind(Position = BEM.obj$uni.pos, ll)
  lower_bounds <- as.data.frame(lower_bounds) %>% 
    tidyr::pivot_longer(!Position, 
                        names_to = "Covariate", 
                        values_to = "Lower")
  
  upper_bounds <- cbind(Position = BEM.obj$uni.pos, hh)
  upper_bounds <- as.data.frame(upper_bounds) %>% 
    tidyr::pivot_longer(!Position, 
                        names_to = "Covariate", 
                        values_to = "Upper")
  
  bounds <- merge(lower_bounds, upper_bounds, by = c("Position", "Covariate"))
  
  estimates <- merge(estimates, bounds, by = c("Position", "Covariate"))
  # Reordering group factor levels
  estimates$Covariate <- factor(estimates$Covariate,      
                                levels = selected_covs)
  
  # creating the plot
  # use aes_string instead of aes to avoid being flagged by R CMD check 
  # "no visible binding for global variable ‘Lower’"
  g <- ggplot(data = estimates, 
              mapping = aes_string(x = "Position", 
                                   y = "Estimate", 
                                   group = "Covariate", 
                                   color = "Covariate")) + geom_line() + 
    geom_line(aes_string(y = "Lower"), linetype = "dashed") +
    geom_line(aes_string(y = "Upper"), linetype = "dashed") + theme_bw() + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_rug(sides = "b", color = "black") + xlab("Genomic position") + 
    ggtitle(title) + theme(legend.position = "none")
  
  if(length(selected_covs) > 1){
    if(same.range){
      g <- g + facet_wrap(~Covariate, 
                          nrow = mfrow[1], ncol = mfrow[2])
    } else {
      g <- g + facet_wrap(~Covariate, scales = "free_y", 
                          nrow = mfrow[1], ncol = mfrow[2])
    }
  }
  
  if(!is.null(save)){
    supported_format <- ".pdf$|png$|tiff$|bmp$|jpeg$|svg$|eps$|tex$|ps$"
    if(length(grep(supported_format, save)) == 0) {
      w_msg <- paste0("Unknown graphics format for '",save,
                      "'. The plot will be saved as a .pdf by default")
      Warning(w_msg)
      ggsave(filename = paste0(save,".pdf"), plot= g, device ="pdf")
    } else {
      ggsave(filename = save, plot = g)
    }
  }
  
  return(g)
}
