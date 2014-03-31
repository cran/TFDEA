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
# $Id: SDEA.R 104 2013-08-12 14:32:56Z tshott $
#
# Super Effecient DEA function
#
#******************************************************************************
#
# TODO
#

# Public SDEA interface
# uppercase SDEA to not mask Benchmarking::dea
# Secondary objective not allowed
#
SDEA <- function(x, y, rts="vrs", orientation="input",
                 slack=FALSE,
                 second="none", z=0,
                 round=FALSE,
                 debug=0){
  
  # Forces super TRUE to perform super efficiency calc
  results <- .dea(x, y, rts, orientation,
                  slack=slack,
                  second=second, z=z,
                  round=round,
                  super=TRUE)
  
  return(results)
}
