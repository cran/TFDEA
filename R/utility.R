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
# Check for max lp print paramaters
# Add debug verbose variables - option
# Log / export info utility - ALL info including orientation, ID string
# * Want to make sure every run is tagged with enough info to reproduce
#
library(lpSolveAPI) 


#
# Package Global Variables
#
tfdea_version = "0.7.1"
tfdea_id = paste0("tfdea Version ", tfdea_version, 
                  " $Id: utility.R 105 2013-08-17 05:47:28Z tshott $")

# Numeric Versions Match Benchmarking Package dea options numeric codes
orientation_options_l <- c("input","output")

# Numeric Versions Match Benchmarking Package dea options numeric codes
# Benchmarking suppors 0 - 7, open TFDEA supports 1-4
rts_options_l   <- c("vrs","drs", "crs", "irs")
rts_rhs         <- c( 1,    1,     0,     1)
rts_typ         <- c( "=", "<=",  ">=",  ">=")

#
# Set global epsilon for numerical compeasions in package
# WARN: If we ever add an option to allow the user to set lpModel control
# parameters we'll need to move this into the solution loop so we get the right
# value for epislon if user changes the tightness of the solution
#
#lp_model  <- make.lp(1, 1)
#lp_ctl    <- lp.control(lp_model)
#epsilon   <- sqrt(lp_ctl$epsilon["epsint"])
#delete.lp(lp_model)
#remove(lp_ctl, lp_model)
epsilon   <- 0.0003162278


# Function: .check_data()
# Check input data and convert to dataframe with both dimensions named
# * Checks for valid types, can be converted to data frame
# * add names if unnamed vars
# var_type is the Prefix string for un-named col vars. it only controls the names
#
.check_data <- function(x, var_type="X"){
  # Is data numeric?
  if ( !is.numeric(x) && !is.data.frame(x) )
    stop(x," is not a numeric matrix (or data.frame)", call. = FALSE)
  
  # Check if > 2d
  if( length(dim(x)) > 2)
    stop(x," is greater than 2 dimensions", call. = FALSE)
  
  row_names_org <- rownames(x)
  col_names_org <- colnames(x)
  
  # Convert to data frame
  x <- as.data.frame(x)
  
  # If original names were missing add names
  if(is.null(row_names_org)){
    row_names = paste0("DMU", 1:nrow(x))
    rownames(x) <- row_names
  }
  if(is.null(col_names_org)){
    col_names = paste0(var_type, 1:ncol(x))
    colnames(x) <- col_names
  }
  
  return(x)
}

# Function: .check_date()
# Check input data and convert to dataframe
# * Checks for valid types, can be converted to dataframe
# * add names if unnamed vars
# * make sure returned as dataframe
# var_type is the Prefix string for un-named col vars. it only controls the names
#
.check_date <- function(x, var_type="X"){
  # Is data numeric?
  if ( !is.numeric(x) && !is.data.frame(x) )
    stop(x," is not a numeric matrix (or data.frame)", call. = FALSE)
  
  # Checks for dataframe or array
  if(!is.null(dim(x))){
    if( length(dim(x)) > 2)
      stop(x," is greater than 2 dimensions", call. = FALSE)
    if( dim(x)[2] > 1)
      stop(x," is greater than 1 col", call. = FALSE)  
    return(x[,1])     
  }  

  return(x)
}

#<New Page> 
# Function: .check_orientation()
# Check that orientation is legal orientation
# Convert to lower case, accept some of the other variants of orientation parameter
# Return cleaned up value
# Checks against global list of legal orientations
#
.check_orientation <- function(orientation){
  
  if(length(orientation) > 1){
    stop("illegal vector for orientation value - must be single is value", call. = FALSE)
  }

  if(is.numeric(orientation)){
    # If numeric value make sure within bounds and convert to string
    if (orientation < 1 || orientation > length(orientation_options_l)){
      stop(paste0("illegal numeric orientation is < 1 or > ", 
                  length(orientation_options_l)), call. = FALSE)
    } else {
      orientation <- orientation_options_l[orientation]
      return(orientation)
    }
  }
  # if character make sure is in list of legal orientations
  if(is.character(orientation)){
    orientation <- tolower(orientation)
    
    # Check in legal list
    if (orientation %in% orientation_options_l)
      return(orientation)
    
    # Convert Benchmarking Package style orientation (IN, OUT)
    # Check for partial string match using pmatch
    orientation_index <- pmatch(orientation, orientation_options_l, nomatch =-1 )
    if (orientation_index > 0){
      orientation <- orientation_options_l[orientation_index]
      return(orientation)
    }
  }
  
  # Not numeric, not character - not legal
  stop("orientation is not legal orientation", call. = FALSE)
}


# Function: .check_rts()
# Check value of RTS and make sure in right format
# Return cleaned up value
#
.check_rts <- function(rts){
  
  if(length(rts) > 1){
    stop("illegal vector for rts value - must be single is value", call. = FALSE)
  }
    
  if(is.numeric(rts)){
    # If numeric value make sure within bounds and convert to string
    if (rts < 1 || rts > length(rts_options_l)){
      stop(paste0("illegal numeric rts is < 1 or > ", 
                  length(rts_options_l)), call. = FALSE)
    } else {
      rts <- rts_options_l[rts]
      return(rts)
    }
  }
  
  # if non-numeric make sure is in list of legal rts
  if(is.character(rts)){
    rts <- tolower(rts)
    if (rts %in% rts_options_l)
      return(rts)
  }
  stop("character rts is not legal rts", call. = FALSE)
}

# Check if value is efficient
# Depends on orientation
# Needs to be smart FP compare
is.efficient <- function(eff, orientation){
  orientation <- .check_orientation(orientation)

  # May get called with NaN value - if NaN, return false.
  # Need to write as vector test
  if(orientation == "input"){
    return( is.finite(eff) & (eff + epsilon) >= 1)
  } else {
    return( is.finite(eff) & (eff - epsilon) <= 1)
  }
}

# Internal Function, check if value is efficient
# Depends on orientation
# Needs to be smart FP compare
is.stdefficient <- function(eff){

  # May get called with NaN value - if NaN, return false.
  # Need to write as part of vector test

  return( is.finite(eff) & (eff + epsilon) >= 1)
}

normalize_eff <- function(eff, orientation){
  if (orientation == "output") eff <- 1 / eff
  return(eff)
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
#
lp_debug_file <- function (lp_model, message="", file_name="lp_model.out", append=TRUE){
  # Check of valid lp object
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
  
   if ( !is.integer(status) || status < 1 || status > length(error_msg) )
     stop("Status code must be integer from 1 - 13", call. = TRUE)
  
  return(error_msg[status])
}

#
# Error calculations
#

# MAD - Mean Absolute Deviation
MAD <- function(x, y){
  mad <- mean(abs(x - y), na.rm=TRUE)
  return(mad)
}
