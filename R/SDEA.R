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
# Super Efficient DEA function
#
#******************************************************************************
#

# Public SDEA interface
# uppercase SDEA to not mask Benchmarking::dea
# Secondary objective not allowed
#
SDEA <- function(x, y, rts="vrs", orientation="input",
                 slack=TRUE,
                 second="none", z=0,
                 round=FALSE, debug=0){

  rts         <- .checkOption(rts,            "rts",              options.rts.l)
  orientation <- .checkOption(orientation,    "orientation",      options.orientation.l)
  slack       <- .checkOption(slack,          "slack",            TRUE)
  second      <- .checkOption(second,         "second",           options.second.l)
  round       <- .checkOption(round,          "round",            TRUE)
  debug       <- .checkOption(debug,          "debug",            0)

  x <- .checkData(x, "x")
  y <- .checkData(y, "y")
  if (nrow(x) != nrow(y))
    stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)

  .checkDataGood(x, y)

  # Check secondary optimization parms
  if (second != "none"){
    z <- .checkVector(z,"z")
    if (nrow(x) != length(z))
      stop("secondary data size (rows) must match number of DMU's", call. = FALSE)
  }

  if (slack && second != "none")
    stop("Can not use both second and slack option; set slack=FALSE", call. = FALSE)

  # Forces super TRUE to perform super efficiency calc
  results <- .dea(x, y, rts, orientation,
                  slack=slack,
                  second=second, z=z,
                  round=round, debug=debug,
                  super=TRUE)

  return(results)
}
