#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014
# Use granted under BSD license terms
#
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
#
# DEA wraper function
#
#******************************************************************************
#

# Public DEA interface
# uppercase dea to not mask Benchmarking::dea
#
DEA <- function(x, y, rts="vrs", orientation="input",
                slack=TRUE, dual=FALSE,
                second="none", z=0,
                round=FALSE, debug=0){

  rts         <- .checkOption(rts,          "rts",          options.rts.l)
  orientation <- .checkOption(orientation,  "orientation",  options.orientation.l)
  slack       <- .checkOption(slack,        "slack",        TRUE)
  dual        <- .checkOption(dual,         "dual",         TRUE)
  second      <- .checkOption(second,       "second",       options.second.l)
  round       <- .checkOption(round,        "round",        TRUE)
  debug       <- .checkOption(debug,        "debug",        0)

  # Check that x & y are legal inputs & convert to standard values
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

  .checkDataGood(x, y)

  if (nrow(x) != nrow(y))
    stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)

  # Check secondary optimization parms
  if (second != "none"){
    z <- .checkVector(z,"z")

    if (nrow(x) != length(z))
      stop("secondary data size (rows) must match number of DMU's", call. = FALSE)
  }

  if (slack && second != "none")
    stop("Can not use both second and slack option; set slack=FALSE", call. = FALSE)

  results <- .dea(x, y, rts=rts, orientation=orientation,
                  slack=slack, dual=dual,
                  second=second, z=z,
                  round=round, debug=debug)

  return(results)
}
