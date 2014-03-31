#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013
# Use granted under BSD license terms
#
# R TFDEA Package
# 
# $Author: tshott $
# $Revision: 104 $
# $Date: 2013-08-12 07:32:56 -0700 (Mon, 12 Aug 2013) $
# $Id: tfdea.R 104 2013-08-12 14:32:56Z tshott $
#
# TFDEA function
#
# NOTE: All efficiencies in program are in 0 - 1 form, output orientated values are
# converted to 1/eff on calculation. use "stdeff" for normalized eff vars
#
# Note: Man page - comment that default orientation is output, different from DEA
#
# Global Dimension Variables
#
# Input Arrays
#
# Type          Variable    Size      iteration var     Dimensions
# Inputs        x           m         i                 i, k
# Outputs       y           s         r                 r, k
# DMU's                     n         k
#
# debug option controls how much output the code genrates. 0 - none, 3 is alot
#
#******************************************************************************
#
# TODO
# Also require utilities
# Check use dim(x) - use value instead
# Question  - why not forecast DMU SEAF < 1?

TFDEA <- function (x, y, dmu_date_rel, date_forecast, 
                   rts="vrs", orientation="output", 
                   second="min", mode="static",
                   debug=0){
  #
  # Check Input Values for correctness and consistency in size
  #
  
  # Check input vars
  x <- .check_data(x)
  y <- .check_data(y)
  if (dim(x)[1] != dim(y)[1])
    stop("Number of DMU's in inputs != number of DMU's in outputs")
  
  
  dmu_date_rel <- .check_date(dmu_date_rel)
  if (dim(x)[1] != length((dmu_date_rel)))
    stop("Number of DMU's in date intro != number of DMU's in inputs", call. = FALSE)
  if(!is.numeric(date_forecast) || length(date_forecast) > 1)
    stop("Forecast year must be numeric value, length = 1", call. = FALSE)
  if(date_forecast < min(dmu_date_rel) || date_forecast > max(dmu_date_rel))
    stop("Forecast year must be between min and max intro_date", call. = FALSE)
  
  # Check other options
  rts         <- .check_rts(rts)
  orientation <- .check_orientation(orientation)
  
  if(! second %in% c("min", "max", "none"))
    stop("second option must be none, min or max")
  
  if(! mode %in% c("static","dynamic"))
    stop("mode option must be must be static or dynamic", call. = FALSE)
  
  if(! (is.numeric(debug) && is.finite(debug) && (debug >= 0 ) ))
    stop("debug option mode must be must numeric >= 0", call. = FALSE)
  
  #
  # Setup output values
  #
  # Create arrays for outputs
  n <- nrow(x)  # number of units, firms, DMUs
  dmu_names     <- rownames(x)
  
  #
  # Calculate unique dates for SOA calculations 
  #
  # Complex expresion, extract unique dates
  # We extract all unique dates so can run eff at rel chek for all DMU's
  date_soa_l <- sort(unique(dmu_date_rel))
  
  # Find the largest date, <= to date forecast
  # The user can pick a forecast date that is not in the DMU release dates, in which case we need to
  # deteime which DMU releaes date is the largest and <= date forecast
  date_cur_forecast <- date_soa_l [max(which( date_soa_l <= date_forecast))]
  
  if (debug >= 2) cat("Date Cur Forecast:", date_cur_forecast, "Unique SOA Dates: ", paste0(date_soa_l, ","), "\n")
  
###############################################################################################
#
# Phase 1 - for each technology SOA set determine which DMU's are efficient at intro
#
# Also store last year efficient to assist in debug
#
###############################################################################################
  #  
  # NOTE: All eff's stdeff - normalized to 0 - 1.  Use Normalize option to DEA

  #
  # Values At Release (REL) date product first introduced
  #  
  dmu_stdeff_rel    <- array(NA, c(n),    list(dmu=dmu_names))
  dmu_date_last     <- array(NA, c(n),    list(dmu=dmu_names))
  
  #
  # Current Values (CUR) - values for last date less then forecast date
  #  
  dmu_stdeff_cur    <- array(NA, c(n),    list(dmu=dmu_names))
  dmu_lambda_cur    <- array(NA, c(n,n),  list(dmu=dmu_names, dmu2=dmu_names))
  
  # TODO - add vx, uy
  #
  # Temp values
  #
  dmu_stdeff_tmp    <- array(NA, c(n),    list(dmu=dmu_names))
  dmu_lambda_tmp    <- array(NA, c(n,n),  list(dmu=dmu_names, dmu2=dmu_names))
  
  
  # Loop for each unique date
  for(t in date_soa_l){
    dmu_stdeff_tmp  <- array(NA, c(n),    list(dmu=dmu_names))  # Make sure all values NA to start
    
    # Build Set in logical (_b - binary) vector
    # All DMU's with date <= t are in SOA set for date t
    dmu_in_soa_b  <- (dmu_date_rel <= t)            # Logical vector - true for dates in SOA
    if (debug >= 3) cat("\nEvaluate SOA for date=", t, " SOA DMUs:", dmu_names[dmu_in_soa_b],"\n")
    
    # Calculate Eff for All DMU's in SOA
    # WARNING: need to subset vector on left & right so results line up
    # WARNING: Use normalize option to make all eff 0 - 1, stdeff, even for output orientation
    results <- .dea(x[dmu_in_soa_b,],
                    y[dmu_in_soa_b,],
                    rts, orientation,
                    second=second, z=dmu_date_rel[dmu_in_soa_b],
                    normalize=TRUE)
    dmu_stdeff_tmp[dmu_in_soa_b] <- results$eff
    dmu_lambda_tmp[dmu_in_soa_b, dmu_in_soa_b] <- results$lambda
    
    # Save results from date that matches forecast year (note that date forecast may not be in 
    # set of dates, so we use the date_cur_forecats, which is largest dmu date <= date forecast)
    if (t == date_cur_forecast){
      dmu_stdeff_cur  <- dmu_stdeff_tmp
      dmu_lambda_cur  <- dmu_lambda_tmp
    }    
    
    # Save release eff for DMU's with dmu_stdeff_rel = NA & dmu_stdeff_tmp != NA
    # This way we save only stdeff when initially released
    # In the future, we may not run eff calculation at every date. This check saves all the
    # dates for new DMU's, even if doesn't match curret date exactly
    dmu_eff_update_b <- is.na(dmu_stdeff_rel) & !is.na(dmu_stdeff_tmp)
    dmu_stdeff_rel[dmu_eff_update_b] <- dmu_stdeff_tmp[dmu_eff_update_b]  
    
    # Update last date effecient for DMU's with dmu_stdeff is efficient, so now last date stdeff
    # Not used in forecast claculations, but produces petty tables
    dmu_date_last[is.stdefficient(dmu_stdeff_tmp)] <- t
  }
  
  table <- cbind(dmu_date_rel, dmu_stdeff_rel, dmu_date_last)
  colnames(table) <- c("Date", "Eff_Rel", "Last")
  if (debug >= 2) {
    print("done Phase 1")
    print(table, digits=7)
  }
  
###############################################################################################
#
# Phase 2 - Calculate technology rate of change (ROC) for DMU's <= forecast date
#
# Calculated and saved neef, lambda in Phase 1
#
###############################################################################################
  #
  # NOTE: dmu_date_cur depends on if using static ot dynamic ROC
  # If static, is just set to forecast date. Set all here, reset if dynamic later
  dmu_date_cur      <- array(date_forecast, c(n),  list(dmu=dmu_names))
  dmu_roc           <- array(NA, c(n),             list(dmu=dmu_names))
  
  # DMU's in current set, date release <= forecast date
  dmu_cur_b  <- (dmu_date_rel <= date_cur_forecast)
  if (debug >= 3) cat("Phase 2 - DMU Forecast Sample Set:", dmu_names[dmu_cur_b], "\n")
  
  # Calculate ROC values for all DMU's < forecast, but not at forecast
  dmu_roc_b  <- (dmu_date_rel < date_cur_forecast)
  if (debug >= 3) cat("Phase 2 - DMU Forecast ROC Calc Set:", dmu_names[dmu_roc_b], "\n")
  
  if(mode == "dynamic"){
    # Dynamic Mode - adjust cur year based on weighted lambdas - forecast year based upon peers
    # year is simple weighted average of peers rel year
    for (k in which(dmu_roc_b)){
      dmu_date_cur[k] <- sum(dmu_lambda_cur[k,] * dmu_date_rel, na.rm=TRUE) / 
        sum(dmu_lambda_cur[k,], na.rm=TRUE)
    }
  }
  
  # ToDO - redo calc with <=, check current date, forecast, one check - not two
  # For each DMU which is earlier then forecast year calc ROC
  # Include DMU's efficeint at release and not efficient at current and eff date > 
  for(k in which(dmu_roc_b)){
    if (is.stdefficient(dmu_stdeff_rel[k]) && !is.stdefficient(dmu_stdeff_cur[k])){

      # Check corner case
      # If we calc a dmu_date_cur < then the actual release date do not use in calc
      if ( dmu_date_cur[k] <= dmu_date_rel[k]){
        warning("DMU ", dmu_names[k], " effective current ",
                "date < release date due to dynamic ROC and is being dropped", call. = FALSE)
      } else {
        dmu_roc[k] <- (dmu_stdeff_rel[k] / dmu_stdeff_cur[k]) ^ 
          ( 1 / (dmu_date_cur[k] - dmu_date_rel[k]) )
        
        # Check for problem with numerical instablity in some dynamic ROC calculations
        # In some cases, dynamic ROC will result in very small delta t and a too large ROC
        if (dmu_roc[k] > 10){
          warning("DMU ", dmu_names[k], " has a numerically unstable ROC and is being dropped",
                  call. = FALSE)
          dmu_roc[k] = NA
        }
      }
    }
  }
  average_roc <- mean(dmu_roc, na.rm=TRUE)
  
  
  
  table <- cbind(dmu_date_rel, dmu_stdeff_rel, dmu_date_last, dmu_stdeff_cur, dmu_date_cur, dmu_roc)
  colnames(table) <- c("Date", "Eff_Rel", "Last", "Eff_Cur", "EDate", "ROC")
  if (debug >= 2) {
    print(c("done Phase 2", "Avg ROC=", average_roc), digits=3)
    print(table, digits=7)
  }
  
###############################################################################################
#
# Phase 3 - Forecast (FOR) intro date
#
# Go through DMU's > forecast, calc forecast date based upon average ROC, 
# calculated stdeff (super efficiency) and cur year
#
###############################################################################################
  #
  dmu_stdeff_for    <- array(NA, c(n),    list(dmu=dmu_names))
  dmu_lambda_for    <- array(NA, c(n,n),  list(dmu=dmu_names, dmu2=dmu_names))
  dmu_date_for      <- array(NA, c(n),    list(dmu=dmu_names))
  
  # Check that Average ROC is valid, if not valid ROC skip section. Likely no DMU's to forecast from
  if(is.numeric(average_roc) && !is.nan(average_roc)){
    # Do foreceast for all DMU's Intro Date > Forecast date
    for (k in which(dmu_date_rel > date_forecast)){
      dmu_stdeff_tmp <- rep(NA, n)                        # Make sure all values NA to start
      
      dmu_test_b <- dmu_cur_b                           # build set, current set DMU's plus
      dmu_test_b[k] <- TRUE                             # plus one DMU from forecast set
      if (debug >= 3) cat("Forecast for SOA DMUs:", dmu_names[dmu_test_b],"\n")
      
      results <- .dea(x[dmu_test_b,],
                      y[dmu_test_b,],
                      rts, orientation,
                      normalize=TRUE,
                      second=second, z=dmu_date_rel[dmu_test_b],
                      super=TRUE)

      dmu_stdeff_tmp[dmu_test_b] <- results$eff
      dmu_lambda_tmp[dmu_test_b, dmu_test_b] <- results$lambda
      
      dmu_stdeff_for[k]  <- dmu_stdeff_tmp[k]
      dmu_lambda_for[dmu_test_b, dmu_test_b] <- results$lambda

      # Question - why Dong-Joon code not forecast < 1 stdeff?
      # Problem ToDo - maybe rounding issue here
      if (!is.finite(dmu_stdeff_for[k])) next
      if (dmu_stdeff_for[k] <= 1){
        warning("DMU ", dmu_names[k], " is not super effecient at forecast and ", 
                "will not be forecast", call. = FALSE)
        next
      }

      if(mode == "dynamic")
        dmu_date_cur[k] <-sum(dmu_lambda_for[k,] * dmu_date_rel, na.rm=TRUE) / 
          sum(dmu_lambda_for[k,], na.rm=TRUE)
      
      # Calculate Forecast from ROC, super effeciency & cur date
      dmu_date_for[k]  <- dmu_date_cur[k] + 
        log(dmu_stdeff_tmp[k], exp(1)) / log(average_roc, exp(1))
    }  
  }  
  
  
  table <- cbind(dmu_date_rel, dmu_stdeff_rel, dmu_date_last, dmu_stdeff_cur, dmu_date_cur, 
                 dmu_roc, dmu_stdeff_for, dmu_date_for)
  colnames(table) <- c("Date", "Eff_Rel", "Last", "Eff_Cur", "EDate", 
                       "ROC", "Eff_For", "Date For")
  if (debug >= 1){ 
    print(c("done Phase 3", "Avg ROC=", average_roc), digits=3)
    print(table, digits=7)
  }
  
  #
  # Return results
  # ToDo - Add option to include ux, vy, w
  results <- list(date_soa=date_soa_l,
                  dmu_eff_rel=dmu_stdeff_rel, dmu_eff=dmu_stdeff_cur,
                  dmu_roc=dmu_roc,
                  dmu_date_cur=dmu_date_cur, dmu_eff_for=dmu_stdeff_for, dmu_date_for=dmu_date_for,
                  roc=average_roc)
  return(results)
}
