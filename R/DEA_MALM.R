#library(TFDEA)
#-----------------------------------------------------------------------------------#
# DEA Malmquist function
# x,y         = inputs and outputs for each DMU for each time interval
# dates       = date for each DMU for each interval (must be the same length as x and y)
# names       = dmu names (must be the same length as x and y)
# rts         = crs or vrs
# orientation = input or output
# scale       = TRUE to return scale efficiency change
#-----------------------------------------------------------------------------------#
DEA_MALM <- function(x, y, dates, names = NULL, rts = "crs", orientation = "input", scale = FALSE){
  # required for checks, will be removed when incorporated into TFDEA package
  # options.orientation.l <- c("input","output")
  # Not sure about IRS and DRS for Malmquist
  options.rts.l <- c("vrs", "crs")

  # -------Start of checks ------- #

  # check inputs and outputs
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")
  if (nrow(x) != nrow(y))
    stop("Length of inputs != length of outputs", call. = FALSE)

  # check rts and orientation
  rts <- .checkOption(rts, "rts", options.rts.l)
  orientation <- .checkOption(orientation, "orientation", options.orientation.l)

  # checks dates
  dates <- .checkVector(dates, "dmu_malm_dates")
  if (length(dates) != nrow(x))
    stop("Number of dates != length of inputs", call. = FALSE)

  # find unique dates and sort
  dates.l <- sort(unique(dates))
  # determine whether at least one interval of dates exists
  if (length(dates.l) < 2)
    stop("More than one year of data required", call. = FALSE)

  # determine number of dmus for each date
  nd <- length(dates)/length(dates.l)

  # check that each year has the same number of dmus
  dmu.dates <- table(dates)
  if (sum(dmu.dates == nd) != length(dates.l))
    stop("Each year must have inputs, outputs and dates for each DMU", call. = FALSE)

  # if no names were given, assign names
  if (is.null(names))
    names <- rep(paste0("DMU", seq(1:nd)), length(dates.l))

  # check whether names is the same length as inputs
  if (length(names) != nrow(x))
    stop("Length of names != length of inputs", call. = FALSE)

  # check whether DMU names are in the same order for each year
  names.l <- rep(names[1:nd], length(dates.l))
  if (sum(names != names.l) > 0)
    stop("DMU names must be in the same order for each year", call. = FALSE)

  # scale = TRUE can only be used for VRS
  if (scale && (rts == "crs")) {
    scale <- FALSE
    warning("scale = TRUE can only be used for rts = vrs", call. = FALSE)
  }
  # -------End of checks ------- #

  # create results array
  # assign dimension names
  row.names <- paste(names[-1:-nd], dates[-1:-nd], sep="_")
  col.names <- c("date", "eff.tt", "eff.t2t2", "eff.tt2", "eff.t2t",
                 "tec", "tc1", "tc2", "tc", "sec", "mpi")
  # if scale is false, remove the scale efficiency change (sec) column
  if (!scale)
    col.names <- col.names[-which(col.names == "sec")]

  # create results array
  results <- array(NA, dim = c(length(row.names), length(col.names)),
                   dimnames = list(row.names, col.names))

  # assign DMU dates to the respective column
  results[, "date"] <- dates[-1:-nd]

  # loop for each date, except the last. The last loop is not required
  for (j in 1:(length(dates.l) - 1)) {

    # determine which rows of results array to save to
    r.rows <- which(results[, "date"] == dates.l[j + 1])

    # eff.tt is the efficiency of time j dmu's using frontier of time j dmu's
    index.T <- which(dates == dates.l[j])
    index.K <- which(dates == dates.l[j])
    eff.tmp <- .dea(x = x, y = y, rts = rts, orientation = orientation, index.K = index.K,
                           index.T = index.T)$eff
    # Remove results that were not part of index.K
    results[r.rows, "eff.tt"] <- eff.tmp[!is.na(eff.tmp)]

    # eff.t2t2 is the efficiency of time j + 1 dmu's using frontier of time j + 1 dmu's
    index.T <- which(dates == dates.l[j + 1])
    index.K <- which(dates == dates.l[j + 1])
    eff.tmp <- .dea(x = x, y = y, rts = rts, orientation = orientation, index.K = index.K,
                            index.T = index.T)$eff
    # Remove results that were not part of index.K
    results[r.rows, "eff.t2t2"] <- eff.tmp[!is.na(eff.tmp)]

    # eff.tt2 is the efficiency of time j dmu's using frontier of time j + 1 dmu's
    index.T <- which(dates == dates.l[j])
    index.K <- which(dates == dates.l[j + 1])
    eff.tmp <- .dea(x = x, y = y, rts = rts, orientation = orientation, index.K = index.K,
                            index.T = index.T)$eff
    # Remove results that were not part of index.K
    results[r.rows, "eff.tt2"] <- eff.tmp[!is.na(eff.tmp)]

    # eff.t2t is the efficiency of time j + 1 dmu's using frontier of time j dmu's
    index.T <- which(dates == dates.l[j + 1])
    index.K <- which(dates == dates.l[j])
    eff.tmp <- .dea(x = x, y = y, rts = rts, orientation = orientation, index.K = index.K,
                            index.T = index.T)$eff
    # Remove results that were not part of index.K
    results[r.rows, "eff.t2t"] <- eff.tmp[!is.na(eff.tmp)]

    # Calculate scale efficiency change (only if using vrs)
    if(scale){
      # calcuate crs efficiencies as above, required for calculating sec
      eff.tmp <- .dea(x = x, y = y, rts = "crs", orientation = orientation,
                                index.K = which(dates == dates.l[j]),
                                index.T = which(dates == dates.l[j]))$eff
      eff.tt.crs <- eff.tmp[!is.na(eff.tmp)]
      eff.tmp <- .dea(x = x, y = y, rts = "crs", orientation = orientation,
                                 index.K = which(dates == dates.l[j + 1]),
                                 index.T = which(dates == dates.l[j + 1]))$eff
      eff.t2t2.crs <- eff.tmp[!is.na(eff.tmp)]
      eff.tmp <- .dea(x = x, y = y, rts = "crs", orientation = orientation,
                                 index.K = which(dates == dates.l[j]),
                                 index.T = which(dates == dates.l[j + 1]))$eff
      eff.t2t.crs <- eff.tmp[!is.na(eff.tmp)]
      eff.tmp <- .dea(x = x, y = y, rts = "crs", orientation = orientation,
                                 index.K = which(dates == dates.l[j + 1]),
                                 index.T = which(dates == dates.l[j]))$eff
      eff.tt2.crs <- eff.tmp[!is.na(eff.tmp)]
      # calcuate scale efficiency change and save to results array
      sec1 <- (eff.tt2.crs / results[r.rows, "eff.tt2"]) / (eff.tt.crs / results[r.rows, "eff.tt"])
      sec2 <- (eff.t2t2.crs / results[r.rows, "eff.t2t2"]) / (eff.t2t.crs / results[r.rows, "eff.t2t"])
      results[r.rows, "sec"] <- (sec1 * sec2)^0.5
    }
  }

  # Calculate (pure) technical efficiency change
  results[, "tec"] <- results[, "eff.t2t2"] / results[, "eff.tt"]

  # Calculate technology change (technology frontier shift)
  results[, "tc1"] <- results[, "eff.tt2"] / results[, "eff.t2t2"]
  results[, "tc2"] <- results[, "eff.tt"] / results[, "eff.t2t"]
  results[, "tc"] <- (results[, "tc1"] * results[, "tc2"])^0.5

  # Calculate Malmquist Productivity index
  results[, "mpi"] <- results[, "tec"] * results[, "tc"]

  return(results)
}
