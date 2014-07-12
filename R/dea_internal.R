#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013, 2014
# Use granted under BSD license terms
#
# R TFDEA Package
#
# NOTE: See utility.R for description of naming nomenclature, use of .l & .b suffix
# and other coding conventions
#
# Internal DEA function
#
#******************************************************************************
#
library(lpSolveAPI)

#
# Internal function - not public
# Function does both dea & sdea
#
.dea <- function(x, y, rts="vrs", orientation="input", slack=TRUE, dual=FALSE, second="none", z=0,
                 round=FALSE, stdeff=FALSE, super=FALSE, fast=FALSE,
                 index.K=NULL, index.T=NULL, debug=0){

  #
  # Function handle solver errors - inside function so can access local vars
  #
  solverError <- function(status, lpm.model, phase, tmp.eff, tmp.var, tmp.dual){

    # Numerical instablity, set bias, try again, update vars
    # Status 5 = "numerical failure encountered"
    if (status == 5){
      set.basis(lpm.model, default = TRUE)
      status <- solve(lpm.model)

      tmp.eff      <-  get.objective(lpm.model)
      tmp.var      <-  get.variables(lpm.model)[1:lpm.col.n]
      tmp.dual     <-  dual.sign * get.dual.solution(lpm.model)
    }


    # Status 2 = "the model is infeasible"
    # Status 3 = "the model is unbounded"
    if ( (status == 2) || (status == 3)){
      if (debug >= 2){
        msg <- lp_solve_error_msg(status)
        cat("Solver Phase ", phase, "Status: ", status, "DMU k=", k,
            " is infeasible (not in technology set) :", msg, "\n")
      }

      tmp.eff      <- ifelse(orientation.in, Inf, -Inf)
      tmp.var      <- rep(NA, lpm.col.n)
      tmp.dual     <- rep(NA, nx+ny+1)

    } else if (status != 0) {

      msg <- lp_solve_error_msg(status)
      cat("Solver Phase ", phase, "Status: ", status, "DMU k=", k, " solver failed: ", msg, "\n")

      tmp.eff      <- NA
      tmp.var      <- rep(NA, lpm.col.n)
      tmp.dual     <- rep(NA, nx+ny+1)
    }

    return (list(eff=tmp.eff, var=tmp.var, dual=tmp.dual))
  }

  #
  # Functions results variables
  #
  dea.status    <- 0
  dea.msgs      <- NULL

  #
  # Setup results names, results output arrays, build technology set
  #
  # Note on xT, yT, and zT: these are the subet of DMU's that are the technology set compared with

  index.K       <- .checkIndex(index.K, x, "index.K")     # DMU's to calc eff for
  index.T       <- .checkIndex(index.T, x, "index.T")     # Technology Set

  nd            <- nrow(x)                  # number of units, firms, DMUs
  names.dmu     <- rownames(x)
  nT            <- length(index.T)          # Can't use nrow, DMU length 1 has no rows

  nx            <- ncol(x)                  # number of inputs
  names.in      <- colnames(x)
  xT            <- x[index.T, ,drop=FALSE]  # Technology Set - inputs WARNING: Must use drop option
                                            # to prevent dimensions from being reduced

  ny            <- ncol(y)                  # number of outputs
  names.out     <- colnames(y)
  yT            <- y[index.T, ,drop=FALSE]  # Technology Set - Outputs

  names.slack   <- paste("S", c(names.in, names.out), sep="_") # IN's + OUT's

  if (second != "none")
    zT          <- z[index.T, drop=FALSE]   # Technology Set - secndary objective

  efficiency    <- array(NA, c(nd),       list(dmu=names.dmu))

  if(!fast)
    lambda      <- array(NA, c(nd, nd),   list(dmu=names.dmu, dmu2=names.dmu))

  if (dual){
    weight.x    <- array(NA, c(nd, nx),   list(dmu=names.dmu, "in"=names.in))  # Must quote in,
    weight.y    <- array(NA, c(nd, ny),   list(dmu=names.dmu, out=names.out))
    weight.w    <- array(NA, c(nd),       list(dmu=names.dmu))
  }

  if(slack){
    slack.x     <- array(NA, c(nd, nx),   list(dmu=names.dmu, "in"=names.in))  # Must quote in,
    slack.y     <- array(NA, c(nd, ny),   list(dmu=names.dmu, out=names.out))
  }

  #
  # Setup values based upon orientation
  #
  if (orientation == "input"){
    orientation.in  <- TRUE
    dual.sign       <- 1
    slack.sign      <- -1
    lpm.sense       <- "min"
  } else {
    orientation.in  <- FALSE
    dual.sign       <- -1
    slack.sign      <- 1
    lpm.sense       <- "max"
  }

  #
  # Model Size
  #
  # Model Cols
  lpm.col.n     <- nT + 1 + nx + ny                        # Size model col DMUs + Eff + Slacks
  lpm.col.eff   <- nT + 1                                  # Eff in col DMU's + 1

  lpm.model.cname <- c(names.dmu[index.T], "Eff", names.slack) # DMUs + Eff + Slacks

  # Model Rows
  lpm.row.n     <- nx + ny + 3                             # Size Model is X's + Y's + 3
  lpm.row.rts   <- nx + ny + 1                             # RTS is rigth after X's & Y's
  lpm.row.super <- nx + ny + 2                             # Super is next row
  lpm.row.slack <- nx + ny + 3                             # Slack is last row

  lpm.model.rname <- c(names.in, names.out, "Scale",       # INs' + OUT's + Scale + Super + Slack
                       "Super", "Phase2")

  #
  # Create initial LP (Linear Program) model
  # ToDo: add lpm control option input
  lpm.model <- make.lp(lpm.row.n, lpm.col.n, verbose = "neutral")
  name.lp(lpm.model, paste0(getSrcFilename(DEA), orientation, rts, tfdea.id, sep="-"))
  dimnames(lpm.model) <- list(lpm.model.rname, lpm.model.cname)


  # Zero lower bound, except for Eff which is -Inf
  set.bounds (lpm.model, lower = rep(0.0, lpm.col.n))
  set.bounds (lpm.model, lower = c(-Inf), columns = c(lpm.col.eff))

  set.objfn(lpm.model,  1, lpm.col.eff)
  lpm.ctl <- lp.control(lpm.model, sense=lpm.sense, timeout=60)
  if (debug >= 2) lpDebugFile(lpm.model, "Orientation Objective Function Set")

  #
  # Calculate epislon, minimum error
  #
  epsilon   <- sqrt(lpm.ctl$epsilon["epsint"])

  #
  # Setup RTS & load into LP model
  #
  # ToDo: Add rts types fdh, add, maybe irs2
  # options.rts.l, rts.rhs, rts,typ global package vars defined in utility file
  # options.rts.l <- c("vrs","drs", "crs", "irs")
  # rts.rhs       <- c( 1,    1,     0,     1)
  # rts.typ       <- c( "=", "<=",   0,    ">=")

  # lookup rts in list of rts options and return index number
  rts.n = match(rts, options.rts.l)

  # Scale constraints
  set.row(lpm.model, lpm.row.rts, rep(1, nT), c(1:nT))      # Fill in 1's for DMU's
  set.rhs(lpm.model,         rts.rhs[rts.n],  lpm.row.rts)  # 1 or 0
  set.constr.type(lpm.model, rts.typ[rts.n],  lpm.row.rts)  # =, >= or <=

  if (debug >= 2) lpDebugFile(lpm.model, "rts constraints added")

  #
  # Setup constant part LP model - This part is the same for every K
  #
  set.constr.type(lpm.model, rep("=", (nx+ny)),     c(1:(nx+ny)))
  for (i in 1:nx) {                                         # Setup IN's
    set.row(lpm.model, i,      -xT[, i],      c(1:nT))
  }
  for (r in 1:ny) {                                         # Setup OUT's
    set.row(lpm.model, nx+r,    yT[, r],      c(1:nT))
  }

  # Setup slack varables - even if not using them - need to set so not free to be assigned
  for(i in 1:(nx+ny))
    set.mat(lpm.model, i, (lpm.col.eff+i), -1)              # Set Slack Diag to -1

  if (debug >= 2) lpDebugFile(lpm.model, "constraint values added")

  # Setup constant part of super that is same for every K
  if (super){
    set.constr.type(lpm.model, "=",  c(lpm.row.super) )
    set.rhs(lpm.model, 0,            c(lpm.row.super) )

    if (debug >= 2) lpDebugFile(lpm.model, paste0("super efficiency constraint added"))
  }


  #
  # Loop for each DMU - Solve one linear equation and extract eff for each DMU
  # ToDo: Use apply?
  for (k in index.K)  {
    #
    # Setup Variable part LP model that changes for each K
    #

    # NOTE: Have to reset Theta since cleared by set.column; setting row 0 = Set Obj Fun
    if (orientation.in){                                    # Input Orientation
      # 1 for theta, Theta Column Inputs set to K'th input
      # Right Hand Side - Outputs set to K'th output
      set.column(lpm.model, lpm.col.eff, c(1,  x[k,]), c(0, (1:nx)))
      set.rhs(lpm.model, y[k,], c((nx+1):(nx+ny)))
    } else {                                                # Output Orientation
      # 1 for theta, theta column Outputs set to - K'th output
      # Right Hand Side - inputs set to -K'th input
      set.column(lpm.model, lpm.col.eff, c(1, -y[k,]), c(0, (nx+1):(nx+ny)))
      set.rhs(lpm.model, -x[k,], c(1:nx))
    }

    #
    # Super Efficiency Calculation
    # If super efficiency - need to add an extra constraint for K'th value = 0
    # Need to check if K is in index.T, find what col it is in
    if (super && (iT <- match(k, index.T, nomatch=0))){
      set.row(lpm.model,                (lpm.row.super), 1, c(iT))   # Theta forced to 1
      if (k == 1 && debug >= 2) lpDebugFile(lpm.model, paste0("super efficiency constraint added"))
    }

    if (debug >= 2){
      lpDebugFile(lpm.model, paste0("Final constraints added k= ",k)) # print for every K
      cat("K=", k, "\n")
    }
    #
    # Solve LP model - Phase 1
    #
    status <- solve(lpm.model)

    # Save temp value all outputs
    tmp.eff         <-  get.objective(lpm.model)
    tmp.var         <-  get.variables(lpm.model)[1:lpm.col.n]
    tmp.dual        <-  dual.sign * get.dual.solution(lpm.model)

    if (status != 0) {
      tmp <- solverError(status, lpm.model, 1, tmp.eff, tmp.var, tmp.dual)
      tmp.eff     <- tmp$eff
      tmp.var     <- tmp$var
      tmp.dual    <- tmp$dual
      dea.status  <- max(dea.status, status)
    }

    # Extract Phase 1 results for DMU K, Need to do before slack minimization
    phase1.eff    <- tmp.eff
    efficiency[k] <- tmp.eff

    if(!fast)
      lambda[k, index.T]  <- tmp.var[1:nT]

    # If dual option selected, save dual weights
    if (dual){
      weight.w[k]   <-    dual.sign * tmp.dual[      1 ]
      weight.x[k,]  <-                tmp.dual[      2 : (nx+1)  ]
      weight.y[k,]  <-                tmp.dual[ (nx+2) : (nx+ny+1)]
    }

    #
    # Phase 2 - slack maximization or secondary objective
    #

    # Setup for slack or secondary objective is almost the same
    if ( is.finite(phase1.eff) && (slack || ( second != "none"))){

      phase1.eff <- ifelse( abs(phase1.eff-1) < epsilon, 1, phase1.eff)       # Round to 1
      phase1.eff <- ifelse( abs(phase1.eff)   < epsilon, 0, phase1.eff)       # Round to 0

      # Add constraint that efficiency must equal efficiency found in Phase 1
      set.row(lpm.model,                (lpm.row.slack), 1, c(lpm.col.eff))   # Put 1 in for eff
      set.constr.type(lpm.model, "=",  c(lpm.row.slack) )                     # Make = constraint
      set.rhs(lpm.model, phase1.eff,   c(lpm.row.slack) )                     # Make eff = phase 1

      if (slack)                              # Maxamize Slacks  - set theta to 0, slacks to -1
        set.objfn(lpm.model, c(0, rep(slack.sign, nx + ny)), c( (lpm.col.eff) : (lpm.col.n) ) )

      if (second != "none"){
        if ( second == "min" ) sign <- -1 else sign <- 1
        set.objfn(lpm.model, sign * zT, c(1:nT))
      }

      if (k == 1 && debug >= 2) lpDebugFile(lpm.model, "Phase 2 - Pre slack maximizations")

      #
      # Solve Model - Phase 2
      #
      status <- solve(lpm.model)

      # Save temp value all results
      tmp.eff      <-  get.objective(lpm.model)           # Store theta in efficiency
      tmp.var      <-  get.variables(lpm.model)[1:lpm.col.n]
      tmp.dual     <-  dual.sign * get.dual.solution(lpm.model)

      if (status != 0){
        tmp <- solverError(status, lpm.model, 2, tmp.eff, tmp.var, tmp.dual)
        tmp.eff  <- tmp$eff
        tmp.var  <- tmp$var
        tmp.dual <- tmp$dual
      }

      if(!fast)
        lambda[k, index.T]    <- tmp.var[1:nT]

      if (slack){
        slack.x[k,] <- tmp.var[(lpm.col.eff+1) : (lpm.col.eff + nx)]
        slack.y[k,] <- tmp.var[(lpm.col.eff + nx + 1) : (lpm.col.eff + nx + ny)]
      }

    } # End Phase 2

    # Cleanup from slack
    set.constr.type(lpm.model, 0,  c(lpm.row.slack) )       # Make slack constraint free
    set.objfn(lpm.model,  1,   c(lpm.col.eff))              # Obj fun reset

  } # End of DMU Loop

  # Cleanup & return Results
  # WARN: Rounding turned off by default. To match other DEA packages, turn rounding on,
  # but does hide some numerical results.

  if(stdeff && !orientation.in){
    efficiency <- 1 / efficiency
  }

  # Round efficiencies that are within epilson to 1 or zero to 1 or zero.
  if(round){
    # Round Efficiency to 0 and 1
    efficiency[abs(efficiency-1) < epsilon] <- 1
    efficiency[abs(efficiency)   < epsilon] <- 0
  }

  results.l   <- list(status=dea.status, eff=efficiency)

  if(!fast)
    results.l <- c(results.l, list(lambda=lambda))

  if (dual)
    results.l <- c(results.l, list(vx=weight.x, uy=weight.y, w=weight.w))

  if(slack)
    results.l <- c(results.l, list(sx=slack.x, sy=slack.y))

  return(results.l)
}


# End dea function


#
# Wrapper with input checking for internal DEA function
#
.dea_safe     <- function(x, y, rts="vrs", orientation="input", slack=TRUE, dual=FALSE, second="none", z=0,
                          round=FALSE, stdeff=FALSE, super=FALSE, fast=FALSE,
                          index.K=NULL, index.T=NULL, debug=0) {

  rts         <- .checkOption(rts,           "rts",         options.rts.l)
  orientation <- .checkOption(orientation,   "orientation", options.orientation.l)
  slack       <- .checkOption(slack,         "slack",       TRUE)
  dual        <- .checkOption(dual,          "dual",        TRUE)
  super       <- .checkOption(super,         "super",       TRUE)
  second      <- .checkOption(second,        "second",      options.second.l)
  round       <- .checkOption(round,         "round",       TRUE)
  fast        <- .checkOption(fast,          "fast",        TRUE)
  stdeff      <- .checkOption(stdeff,        "stdeff",      TRUE)
  debug       <- .checkOption(debug,         "debug",       0)

  # Check for legal input values
  #
  # Check that x & y are legal inputs & convert to standard values
  x <- .checkData(x, "x")
  y <- .checkData(y, "y")

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

  #
  # Save arguments for print when done - record of what options used
  # Save before clean up?
  arguments <- list(version=tfdea.id,
                    rts=rts, orientation=orientation, round=round, stdeff=stdeff,
                    slack=slack, dual=dual, super=super, second=second, z=z)

  return(.dea(x, y, rts=rts, orientation=orientation, slack=slack, dual=dual, second=second, z=z,
              round=round, stdeff=stdeff, super=super, fast=fast,
              index.K=index.K, index.T=index.T, debug=debug))

}

