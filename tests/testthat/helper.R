# skip function
skip_if_run_bitbucket_pipeline <- function() {
  if (identical(Sys.getenv("BB_PIPELINE"), "")) {
    return(invisible(TRUE))
  }
  
  skip("Not run in bitbucket pipeline")
}
