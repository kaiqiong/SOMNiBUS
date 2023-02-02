# The following code has been copied on december 15, 2021 
# from the repo mjkallen/rlogging 
# https://github.com/mjkallen/rlogging (last commit: d6f2785 on Aug 9, 2016)

#' @title Set log level
#'
#' @description Set log level
#' @param level level to assign
#' @author mjkallen, Audrey Lemaçon
#' @noRd
SetLogLevel <- function(level="INFO") {
  stopifnot(level %in% c("INFO", "WARN", "STOP"))
  assign("loglevel", level, envir=.SOMNiBUSLogOpts)
}

#' @title Get log level
#'
#' @description Get log level
#' @return character, current log level
#' @author mjkallen, Audrey Lemaçon
#' @noRd
GetLogLevel <-function() {
  get("loglevel", envir=.SOMNiBUSLogOpts)
}

#' @title Set time stamp format
#'
#' @description Set time stamp format
#' @param ts.format format to apply to time stamp
#' @author mjkallen, Audrey Lemaçon
#' @noRd
SetTimeStampFormat <- function(ts.format="[%Y-%m-%d %H:%M:%S]") {
  assign("ts.format", ts.format, envir=.SOMNiBUSLogOpts)
}

#' @title Get time stamp format
#'
#' @description Set time stamp format
#' @return current format applied to time stamp
#' @author mjkallen, Audrey Lemaçon
#' @noRd
GetTimeStampFormat <- function() {
  get("ts.format", envir=.SOMNiBUSLogOpts)
}

#' @title Set message suffixes
#'
#' @description Set message suffixes
#' @param file.name.suffixes list of suffixes for each log level
#' @author mjkallen, Audrey Lemaçon
#' @noRd
SetFilenameSuffixes <- function(file.name.suffixes=list(INFO="message",
                                                        WARN="warning",
                                                        STOP="stop")) {
  # Check that a list of length=3 is passed to this function:
  if (!is.list(file.name.suffixes) | length(file.name.suffixes) != 3) {
    stop("argument file.name.suffixes must be a list of length=3.")
  }
  
  # Check if the list of suffixes contains INFO, WARN, and STOP. If not,
  # Return a message indicating the missing elements of file.name.suffixes
  missing.elements <- setdiff(c("INFO", "WARN", "STOP"),
                              names(file.name.suffixes))
  
  if (length(missing.elements)) {
    stop("argument file.name.suffixes is missing element(s): ",
         paste(missing.elements, sep="", collapse=", "))
  }
  
  assign("file.name.suffixes", file.name.suffixes, envir=.SOMNiBUSLogOpts)
}

#' @title Get message suffixes
#'
#' @description Get message suffixes
#' @return list of suffixes for each log level
#' @author mjkallen, Audrey Lemaçon
#' @noRd
GetFilenameSuffixes <- function() {
  get("file.name.suffixes", envir=.SOMNiBUSLogOpts)
}

#' @title Set log file
#'
#' @description Set log file
#' @param base.file log file base name
#' @param folder path to folder in which create the log file
#' @param split.files logical, determines whether the logs should be splitted
#' according to their levels.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
SetLogFile <- function(base.file="somnibus.log", folder=tempdir(),
                       split.files=FALSE) {
  assign("split.files", split.files, envir=.SOMNiBUSLogOpts)
  
  if (is.null(base.file)) {
    assign("logfile.base", NULL, envir=.SOMNiBUSLogOpts)
  } else {
    assign("logfile.base", file.path(folder, base.file),
           envir=.SOMNiBUSLogOpts)
  }
}

#' @title Get log file
#'
#' @description Get log file
#' @param level log level of interest
#' @return log file
#' @author mjkallen, Audrey Lemaçon
#' @noRd
GetLogFile <- function(level) {
  base.file <- get("logfile.base", envir=.SOMNiBUSLogOpts)
  split.files <- get("split.files", envir=.SOMNiBUSLogOpts)
  if (!missing(level)) {
    if (!split.files) {
      warning("level=", level, "provided to GetLogFile(), but log files
              are not split. Ignoring parameter.")
      base.file
    } else {
      file.name.suffix <- get(level, GetFilenameSuffixes())
      replacement <- paste("\\1", "_", file.name.suffix, "\\2", sep="")
      sub("(.+?)(\\.[^.]*$|$)", replacement, base.file)
    }
  } else {
    if (split.files) {
      stop("log files are split, but no level parameter provided to
              GetLogFile().")
    } else {
      base.file
    }
  }
}


#' @title Diagnostic Messages
#'
#' @description Generate a diagnostic message from its arguments.
#' @param ... zero or more objects which can be coerced to character (and which
#' are pasted together with no separator) or (for message only) a single 
#' condition object.
#' @param domain see gettext. If NA, messages will not be translated, see also 
#' the note in stop.
#' @param appendLF logical: should messages given as a character string have a
#' newline appended?
#' @param expr expression to evaluate.
#' @param step text to display before the name of the calling function. 
#' Set to NULL to avoid displaying the name of the calling function.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
Message <- function(..., domain=NULL, appendLF=TRUE, step="Running") {
  
  args <- list(...)
  is.condition <- length(args) == 1L && inherits(args[[1L]], "condition")
  if (is.condition) {
    # bypass the logger if a condition is supplied or if loglevel is set to "NONE"
    base::message(..., domain=domain, appendLF=appendLF)
  } else {
    # if loglevel is set to INFO, then print log message, else do nothing
    loglevel <- GetLogLevel()
    split.files <- get("split.files", envir=.SOMNiBUSLogOpts)
    if (loglevel == "INFO") {
      if (split.files) {
        PrintLogMessage(..., level="INFO")
      } else {
        # if step not NULL, add the name of the calling function in and the 
        # content of step in the log message
        if(!is.null(step)){
          calling_fn <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))[1]
          calling_fn <- gsub(x = calling_fn, pattern = "([A-z])\\(.*", 
                             replacement = "\\1()")
          calling_fn <- paste0(step," ",calling_fn,": ")
        } else {
          calling_fn <- ""
        }
        PrintLogMessage("[INFO] ",calling_fn, ...)
      }
    }
  }
  invisible()
  
}


.SOMNiBUSLogOpts <- new.env()

#' @title Set the environment when loading the package
#'
#' @description Set the environment when loading the package
#' @param libname a character string giving the library directory where the
#' package defining the namespace was found.
#' @param pkgname a character string giving the name of the package.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
.onLoad <- function(libname, pkgname) {
  SetTimeStampFormat()
  SetFilenameSuffixes()
  SetLogFile()
  SetLogLevel()
  lockEnvironment(.SOMNiBUSLogOpts)
}

#' @title Set the environment when loading the package in interactive mode
#'
#' @description Set the environment when loading the package
#' @param libname a character string giving the library directory where the
#' package defining the namespace was found.
#' @param pkgname a character string giving the name of the package.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
.onAttach <- function(libname, pkgname) {
  pkgversion <- read.dcf(system.file("DESCRIPTION", package=pkgname), 
                         fields="Version")
  msg <- paste("Package", pkgname, "version", pkgversion, 
               "\nNOTE: - Logging level is set to", GetLogLevel(), 
               "\n      - Output file is", GetLogFile(), 
               "\n      - See 'package?SOMNiBUS' for help.")
  packageStartupMessage(msg)
}

#' @title Print log message
#'
#' @description Print log message
#' @param ... zero or more objects which can be coerced to character (and which
#' are pasted together with no separator) or (for message only) a single 
#' condition object.
#' @param domain see gettext. If NA, messages will not be translated, see also 
#' the note in stop.
#' @param level log level
#' @author mjkallen, Audrey Lemaçon
#' @noRd
PrintLogMessage <- function(..., domain=NULL, level) {
  
  timestamp <- format(Sys.time(), format=GetTimeStampFormat())
  base::message(timestamp, ..., domain=domain)
  #cat(timestamp, ..., "\n", sep="")

  # print log message to file if this one is not set to NULL
  if (missing(level)) {
    logfile <- GetLogFile()
  } else {
    logfile <- GetLogFile(level=level)
  }
  
  if (!is.null(logfile)) {
    cat(timestamp, ..., "\n", file=logfile, sep="", append=TRUE)
  }
  
}

#' @title Stop Function Execution
#'
#' @description Stops execution of the current expression and executes an error
#' action.
#' @param ... zero or more objects which can be coerced to character (and which
#' are pasted together with no separator) or a single condition object.
#' @param call. logical, indicating if the call should become part of the error
#' message.
#' @param domain see gettext. If NA, messages will not be translated, see also 
#' the note in stop.
#' @param step text to display before the name of the calling function. 
#' Set to NULL to avoid displaying the name of the calling function.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
Error <- function(..., call.=TRUE, domain=NULL, step="In") {
  
  args <- list(...)
  is.condition <- length(args) == 1L && inherits(args[[1L]], "condition")
  if (!is.condition) {
    split.files <- get("split.files", envir=.SOMNiBUSLogOpts)
    # error messages are always printed (i.e. for levels INFO, WARN and STOP)
    if (split.files) {
      PrintLogMessage(..., level="STOP")
    } else {
      # if step not NULL, add the name of the calling function in and the 
      # content of step in the log message
      if(!is.null(step)){
        calling_fn <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))[1]
        calling_fn <- gsub(x = calling_fn, pattern = "([A-z])\\(.*", 
                           replacement = "\\1()")
        calling_fn <- paste0(step," ",calling_fn,": ")
      } else {
        calling_fn <- ""
      }
      PrintLogMessage("[STOP] ",calling_fn, ...)
    }
  }
  base::stop(..., call.=call., domain=domain)
  invisible()
  
}

#' @title Warning Messages
#'
#' @description Generates a warning message that corresponds to its argument(s
#' and (optionally) the expression or function from which it was called.
#' @param ... zero or more objects which can be coerced to character (and which
#' are pasted together with no separator) or a single condition object.
#' @param call. logical, indicating if the call should become part of the error
#' message.
#' @param immediate. logical, indicating if the call should be output 
#' immediately, even if getOption("warn") <= 0.
#' @param noBreaks. logical, indicating as far as possible the message should be
#' output as a single line when options(warn = 1).
#' @param expr expression to evaluate.
#' @param domain see gettext. If NA, messages will not be translated, see also 
#' the note in stop.
#' @param step text to display before the name of the calling function. 
#' Set to NULL to avoid displaying the name of the calling function.
#' @author mjkallen, Audrey Lemaçon
#' @noRd
Warning <- function(..., call.=TRUE, immediate.=FALSE, domain=NULL, 
                    step = "In") {
  
  args <- list(...)
  is.condition <- length(args) == 1L && inherits(args[[1L]], "condition")
  if (!is.condition) {
    # if loglevel is set to INFO or WARN, then print log message
    loglevel <- GetLogLevel()
    split.files <- get("split.files", envir=.SOMNiBUSLogOpts)
    if (loglevel %in% c("INFO", "WARN")) {
      if (split.files) {
        PrintLogMessage(..., level="WARN")
      } else {
        # if step not NULL, add the name of the calling function in and the 
        # content of step in the log message
        if(!is.null(step)){
          calling_fn <- as.character(deparse(sys.calls()[[sys.nframe()-1]]))[1]
          calling_fn <- gsub(x = calling_fn, pattern = "([A-z])\\(.*", 
                             replacement = "\\1()")
          calling_fn <- paste0(step," ",calling_fn,": ")
        } else {
          calling_fn <- ""
        }
        PrintLogMessage("[WARN] ",calling_fn, ...)
      }
    }
    # always collect warnings when printing log messages
    immediate. <- FALSE
  }
  # always call the base warning function to collect warnings
  base::warning(..., call.=call., immediate.=immediate., domain=domain)
  invisible()
  
}




