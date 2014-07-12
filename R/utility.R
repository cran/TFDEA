#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013
# Use granted under BSD license terms
#
# R TFDEA Package
#
# $Author: tshott $
# $Revision: 105 $
# $Date: 2013-08-16 22:47:28 -0700 (Fri, 16 Aug 2013) $
#
# Utilities File
# General internal use only utility functions
#
#******************************************************************************
#
# TODO
# Check for max lp print parameters
# Add debug verbose variables - option
# Log / export info utility - ALL info including orientation, ID string
# * Want to make sure every run is tagged with enough info to reproduce
#
library(lpSolveAPI)

#
# Package Global Variables
#
tfdea.version = "0.9.4"
tfdea.id = paste0("tfdea Version ", tfdea.version)

#
# Options Strings
#
options.orientation.l <- c("input","output")
options.second.l      <- c("min", "max", "none")
options.mode.l        <- c("static", "dynamic")

options.rts.l         <- c("vrs","drs", "crs", "irs")
rts.rhs               <- c(  1,   1,     0,     1)
rts.typ               <- c( "=", "<=",   0,    ">=")

#
# Set global epsilon for numerical checks in package
# Backup value, dea_internal resets when run
epsilon   <- 0.0003162278


# Function: .checkData()
# * Checks for valid types
# * add names if unnamed vars
#
.checkData <- function(x, name){
  # Is value a dataframe, array or vector?
  if ( !is.data.frame(x) && !is.array(x) && (length(x) < 2))
    stop(name," is not a matrix, array (1 or 2 dimensions) or data.frame", call. = FALSE)

  if(length(dim(x)) > 2)
    stop(name," is greater then 2 dimensions", call. = FALSE)

  # If data.frame - convert to array
  if (is.data.frame(x))
    x <- data.matrix(x)

  # If vector convert to a 1 x N array
  if (is.vector(x))
    x <- array(x, c(length(x), 1))

  # Check that all numeric
  if(! is.numeric(x))
    stop(name," must be numeric", call. = FALSE)

  # Add DMU names if empty
  if (length(dim(x)) >= 1 && is.null(rownames(x)))
    rownames(x) <- paste0("DMU", 1:nrow(x))

  # Add Col names if empty
  if (length(dim(x)) >= 2 && is.null(colnames(x)))
    colnames(x) <- paste0(name, 1:ncol(x))

  return(x)
}

# Function .CheckDataGood
# Check input & output data and warn about problems for DEA
.checkDataGood <- function(x, y){

  status <- TRUE

  # Check for any zero's
  # Check for NA's
  # Check for no positive Values on input & output
  if (any(x == 0, na.rm=TRUE)){
    cat("Warning, data has DMU's with inputs that are zero, this may cause numerical problems\n")
    status <- FALSE
  }

  if (any(y == 0, na.rm=TRUE)){
    cat("Warning, data has DMU's with outputs that are zero, this may cause numerical problems\n")
    status <- FALSE
  }

  for (k in (1:nrow(x))){
    if ( all(x[k,] <= 0, na.rm=TRUE) & all(y[k,] <= 0, na.rm=TRUE)){
      cat("Warning, DMU # ",k, "has no positive inputs or outputs, this may cause numerical problems\n")
    status <- FALSE
    }
  }

  return(status)
}

# Function: .checkVector()
# Check input data and convert to dataframe
# * Checks for valid types
# * add names if unnamed vars
# * make sure returned as array
#
.checkVector <- function(x, name){

  # Is value a dataframe, array or vector?
  if ( !is.data.frame(x) && !is.array(x) && !is.vector(x))
    stop(name," is not a vector, array (1 dimensions) or data.frame", call. = FALSE)

  if(is.data.frame(x))
    x <- data.matrix(x)

  if(!is.null(dim(x))){
    if (length(dim(x)) > 2)
      stop(name," is greater then 2 dimensions", call. = FALSE)
    if (length(dim(x)) == 2 & dim(x)[2] > 1)
      stop(name," is greater then 2 dimensions", call. = FALSE)

    x <- as.vector(x)
  }

  if (!is.numeric(x))
    stop(name, " must be numeric ", call. = FALSE)


#   if (!is.null(ncol(x)) && (ncol(x) != 1))
#       stop(name, " must have one col ", call. = FALSE)
#

  return(x)
}


# Function: .checkOption()
# Check that value in legal list options
# If on legal list, return lowercase short form option
#
.checkOption <- function(value, name, options.l=c(TRUE, FALSE)){

  if (length(value) != 1){
    stop("illegal vector for ", name, " option must be single value", call. = FALSE)
  }

  # If options are character
  if (is.character(options.l[1])){
    if (is.character(value)){
      tmp.value <- tolower(value)
      #    if (tmp.value %in% options.l)
      i <- charmatch(tmp.value, options.l, nomatch=-1)
      if (i > 0)
        return(options.l[i])
    }
  } else if (is.logical(options.l[1])){
    # Logical options
    if (is.logical(value)){
      return(value)
    }
  } else {
    # Numeric options
    if (is.numeric(value) && is.finite(value) && (value >= 0 )){
      return(value)
    }
    options.l <- c("numeric >= 0")
  }

  stop("Illegal value=", value, " for ", name, "; legal values are: ",
       paste(options.l, collapse = ", "),
       call. = FALSE)
}

#
# Check that index is OK, comvert to sparse if not
.checkIndex       <- function(index, x, name){

  # Add check for are index values valid?

  if(is.null(index)){                     # index.K defines the subset of DMU's that eff is
    return(c(1:nrow(x)))
  }

  if (is.logical(index)){
    if (length(index) != nrow(x))
      stop("Length ",name,"= ",length(index),"  must equal number of DMU's", call. = FALSE)
    return(which(index))                # Acccept logocal or spase index array - convert to spars
  }

  if (is.numeric(index)){
    if (any(index < 1) | any (index > nrow(x)))
      stop("Value of ",name, "=", paste(index, sep=",", collapse=",")," is < 1 or > number of DMU's", call.=FALSE)

    if(any(duplicated(index))){
      warning("a value of ",name, " is duplicated", call.=TRUE)
    }
  }
  return(index)
}

# Check if value is weekly efficient
# Depends on orientation
# Needs to be smart FP compare
isEfficient <- function(eff, orientation){

  if (!is.numeric(eff))
    stop("Illegal value= eff; legal values are: numeric", call. = FALSE)

  orientation <- .checkOption(orientation,   "orientation",  options.orientation.l)

  if(orientation == "input"){
    return (is.finite(eff) & (eff + epsilon  >= 1))
  } else {
    return (is.finite(eff) & (eff - epsilon  <= 1))
  }
}

# Internal Function, check if value is efficient
# Depends on orientation
# Needs to be smart FP compare
isStdEfficient <- function(eff){

  # May get called with NaN value - if NaN, return false.
  # Need to write as part of vector test

  return( is.finite(eff) & (eff + epsilon) >= 1)
}


#<New Page>
# Function: .lp_debug_file()
#
# Utility function to assist in debugging a linear model.
# Dumps out linear model to a output file.
# Pass a lp_model to print, a log message to use
#
# Use append=FALSE when you want to start a new file, for a new test run for
# example.
#
# NOTE: lp model print function has limit on how big a model it will print.
#       it quietly fails to print the vars if model is too big.
#
# ToDo:
# add a option to control sig digits of values printed
# add a check for if max # vars is being executed and print will fail quietly
# try using another method to dump model out
#
lpDebugFile <- function (lp_model, message="", file_name="lp_model.out", append=TRUE){
  # Check if valid lp object
  if(! (class(lp_model) == "lpExtPtr") )
    stop("not a lp_object", call. = FALSE)

  sink(file_name, type="output", append)  # Send output to file

  cat("\n")
  cat(paste0(message),"\n")
  print(lp_model)
  cat("\n")

  sink(NULL)                              # Resets output to normal location
}

#
# Utility to print the error message from lp solve for non-zero status
#
lp_solve_error_msg <- function(status){

  error_msg <- c(
    "the model is sub-optimal",
    "the model is infeasible",
    "the model is unbounded",
    "the model is degenerate",
    "numerical failure encountered",
    "process aborted",
    "timeout",
    "the model was solved by presolve",
    "the branch and bound routine failed",
    "the branch and bound was stopped because of a break-at-first or break-at-value",
    "a feasible branch and bound solution was found",
    "no feasible branch and bound solution was found")

  if ( !is.numeric(status) || status < 1 || status > length(error_msg) ){
    warning("Status code:", status, " must be integer from 1 - 12", call. = TRUE)
    message <- paste("ERROR: Unknown status message: ", status)
  } else {
    message <- error_msg[status]
  }
  return(message)
}

#
# Error calculations
#

# MAD - Mean Absolute Deviation
# Todo: Add warning messages for NA's
MAD <- function(x, y){
  mad <- mean(abs(x - y), na.rm=TRUE)
  return(mad)
}
