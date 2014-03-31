#******************************************************************************
#
# Copyright Tom Shott, ETA, PSU, 2013
# Use granted under BSD license terms
#
# R TFDEA Package
# 
# $Author: tshott $
# $Revision: 99 $
# $Date: 2013-08-10 12:53:47 -0700 (Sat, 10 Aug 2013) $
# $Id: DEA.R 99 2013-08-10 19:53:47Z tshott $
#
# NOTE: See utility.R for description of naming nomenclature, use of _l & _b suffix 
#
# DEA function
#
#
# Global Dimension Varables
#
# Global Varables
#
# Input Arrays
#
# Type          Variable    Size      iteration var     Dimensions
# Inputs        x           m         i                 i, k
# Outputs       y           s         r                 r, k
# DMU's                     n         k
#
# Output Arrays
#
# Type          Variable    Size      iteration var     Dimensions
# Effeciency    effeciency  n         k                 n
# Lambda        lambda      n                           n, n
# ux            x_weight    n                           n, m     
# vy            y_weight    n                           n, s
#
#******************************************************************************
#
# TODO
# * Add debug prints
# * check secondary_data size
# * should parm names match benchmarking?
# Also require utilities
library(lpSolveAPI) 

# Public DEA interface
# uppercase dea to not mask Benchmarking::dea
DEA <- function(x, y, rts="vrs", orientation="input", 
                round=FALSE, 
                slack=FALSE, dual=FALSE,
                second="none", z=0,
                debug=0){
  
  results <- .dea(x, y, rts, orientation, 
                  round=round, normalize=FALSE,
                  slack=slack, dual=dual,
                  super=FALSE, 
                  second=second, z=z,
                  debug=debug)
  
  return(results)
}

#
# Internal function - not public
# Function does both dea & sdea
#

.dea <- function(x, y, rts="vrs", orientation="input", 
                 round=FALSE, normalize=FALSE,
                 slack=FALSE, dual=FALSE,
                 super=FALSE, 
                 second="none", z=0,
                 debug=0){
    
  # Check for legal input values
  # Check functions also cleanup and convert to standard internal values
  #
  # Check that x & y are legal inputs & convert to standard values
  x <- .check_data(x, "I")
  y <- .check_data(y, "O")
  if (dim(x)[1] != dim(y)[1]){
    stop("Number of DMU's in inputs != number of DMU's in outputs", call. = FALSE)
  }  

  rts         <- .check_rts(rts)
  orientation <- .check_orientation(orientation)

  # ToDo convert to check functions
  if(! (is.logical(round) && length(round) == 1) )
    stop("round option must be True or False, length = 1", call. = FALSE)
  
  if(! (is.logical(normalize) && length(normalize) == 1))
    stop("normalize option must be True or False, length = 1", call. = FALSE)
  
  if(! (is.logical(super) && length(super) == 1))
    stop("super option must be True or False, length = 1", call. = FALSE)

  if(! (is.logical(dual) && length(dual) == 1))
    stop("dual option must be True or False, length = 1", call. = FALSE)

  if(! (is.logical(slack) && length(slack) == 1))
    stop("slack option must be True or False, length = 1", call. = FALSE)

  if( slack && second != "none")
    stop("Can not use both second and slack option", call. = FALSE)

  # Check secondary optimization parms
  if (second != "none"){
    if(! second %in% c("min", "max")){
      stop("second must be none, min or max", call. = FALSE)
    } else {
      # Is a secondary onjective function, check secondary data right size
      z <- .check_data(z,"Z")
      
      if(dim(x)[1] != dim(z)[1]){
        stop("secondary data size must match number of DMU's", call. = FALSE)
      }      
    }
  }
  
  #
  # Save arguments for print when done - record of what done
  #
  arguments <- list(version=tfdea_id,
                    rts=rts, orientation=orientation,
                    round=round, normalize=normalize,
                    slack=slack, dual=dual,
                    super=super,
                    second=second, z=z)
  
  #
  # Setup results names, results output arrays
  #
  m <- ncol(x)  # number of inputs
  in_names      <- colnames(x)

  n <- nrow(x)  # number of units, firms, DMUs
  dmu_names     <- rownames(x)
  
  s <- ncol(y)  # number of outputs
  out_names     <- colnames(y)

  slack_names   <- paste("S", c(in_names, out_names), sep="_") # IN's + OUT's
  

  efficiency    <- array(NA, c(n),                        # [ #DMU ]
                        list(dmu=dmu_names))
  lambda        <- array(NA, c(n, n),                     # [ #DMU x #DMU ]
                        list(dmu=dmu_names, dmu2=dmu_names))

  if (dual){      
    x_weight    <- array(NA, c(n, m),                     # [ #DMU x #IN ]
                         list(dmu=dmu_names, "in"=in_names))  # Must quote in, reserved word
    y_weight    <- array(NA, c(n, s),                     # [ #DMU x #OUT] 
                         list(dmu=dmu_names, out=out_names))  
  }
  
  if(slack){
    x_slack       <- array(NA, c(n, m),                     # [ #DMU x #IN ]
                           list(dmu=dmu_names, "in"=in_names))  # Must quote in, reserved word
    y_slack       <- array(NA, c(n, s),                     # [ #DMU x #OUT] 
                           list(dmu=dmu_names, out=out_names))  
  }
  
  # Setup model Col & Row names for later use in DMU loop
  lp_model_cname <- c(dmu_names, "Eff", slack_names)    # DMUs + Eff + Slacks
  lp_model_rname <- c(in_names, out_names, "Scale",     # INs' + OUT's + Scale + Super + Sloack
                      "Super", "Phase2")     

# Linear Model Description
#
# Layout of linear model - Input Orientation, for Kth DMU solutio
#
#           1       2             n       n+1     n+2       n+1+m   n+m+2     n+1+m+s
#           DMU1    DMU2    ...   DMU     Theta   S_X1  ..  S_Xm .. S_Y1  ..  S_Ys        RHS
                                                  
# 0 Obj     Lamda1  Lambda2 ...   Lambdan 1       0         0       0         0           
# 1         X1(1)   X1(2)                 X1(k)   -1        0       0         0     =     0
# ...       ...
# m         Xm                            Xm(k)   0         0       0         0     =     0
# m + 1     Y1                                    0         0       -1        0     =     Y1(k)
# ...       ...
# m + s     Ys                                                                      =     Ym(k)
# 
# Scale / RTS
# m + s + 1 1       1             1       0       0         0       0         0     >=    1
#
# Super Effeciency
# m + s + 2 0       0        Kth 1                                                 =     0         
#
# Slack / Secondary 
# m + s + 3 Slack
#
#
#
  
  
  #
  # Loop for each DMU
  # Solve one linear equation and extract eff for each DMU
  #
  for (k in 1:n)  {
    #
    # Create initial LP (Linear Program) model [IN's + OUT's + Scale, X #DMUS + Eff + Slacks]
    #
    # Cols (# Vars) = #DMU + 1 for EFF
    lp_model <- make.lp( m + s + 3, n + 1 + m + s)
    name.lp(lp_model, paste0(getSrcFilename(DEA), orientation, rts, sep="-", tfdea_id))
    dimnames(lp_model) <- list(lp_model_rname, lp_model_cname)
    
    # set lambdas, eff
    set.bounds (lp_model, lower = rep(0.0, n + 1 + m + s))    # >= 0
    if (k == 1 && debug >= 2) 
      lp_debug_file(lp_model, "Initial Model Created", append=FALSE) # Start new file
    
    # Setup Orientation & Theta = 1
    if (orientation == "input"){
      # Input Orientation - Phase 1:  Minimize theta, all other coefficients 0
      set.objfn(lp_model,  1, n+1)
      lp.control(lp_model, sense="min")
    } else if (orientation == "output") {
      # Output Orientation - Phase 1:  Maximize theta, all other coefficients 0
      set.objfn(lp_model,  1, n+1)
      lp.control(lp_model, sense="max")
    } else {
      stop("Only orientation=input or output supported", call. = FALSE)
    }
    if (k == 1 && debug >= 2) lp_debug_file(lp_model, "Orientation Objective Function Set")  
    
    #
    # Setup CTS & load into LP model
    #
    # ToDo: Add rts types fdh, add, maybe irs2
    # rts_options_l, rts_rhs, rts,typ global package vars defined in utility file
    # rts_options_l <- c("vrs","drs", "crs", "irs")
    # rts_rhs       <- c( 1,    1,     0,     1)
    # rts_typ       <- c( "=", "<=",  ">=",  ">=")
    
    # lookup rts in list of rts options and return numerix index number
    rts_n = match(rts, rts_options_l)
    if (is.null(rts_n))
      stop("Only rts=vrs, cts, irs, drs supported", call. = FALSE)
    
    # Scale constraint is m+s+1 row in table, row = #INs + #Outs + 1
    set.row(lp_model, m+s+1, rep(1, n), c(1:n))	        # Fill in 1's for IN & OUT
    set.rhs(lp_model,         rts_rhs[rts_n],  m+s+1)   # 1 or 0
    set.constr.type(lp_model, rts_typ[rts_n],  m+s+1)   # =, >= or <=
    if (k == 1 && debug >= 2) lp_debug_file(lp_model, "rts constraints added")  
    
    #
    # Setup constant part LP model
    # This part is the same for every K
    #
    for (i in 1:m) {                                    # Setup IN's
      set.row(lp_model, i,     -x[,i],  c(1:n))
      set.constr.type(lp_model, "=",    c(i))
    }
    for (r in 1:s) {                                    # Setup OUT's
      set.row(lp_model, m+r,    y[,r],  c(1:n))
      set.constr.type(lp_model, "=",    c(m+r))
    }

    # Setup slack varables - even if not using them - need to set so not free to be assigned
    for(i in 1:m){    # Set Input Slack to -1
      set.mat(lp_model, i, (n+1+i), -1)
    }
    for(r in 1:s){    # Set Input Slack to -1
      set.mat(lp_model, (m+r), (n+m+1+r), -1)
    }
    
    if (k == 1 && debug >= 2) lp_debug_file(lp_model, "constraint values added")  
    
    #    
    # Setup Variable part LP model that changes for each K
    #
    # Theta is stored in Col after In's & Out's - so Col is #In's + #Out's + 1
    # NOTE: Have to reset Theta since cleared by set.column; setting row 0 = Set Obj Fun
    if (orientation == "input"){
      # Input orientation
      # 1 for theta, THeta Colume Inputs set to K'th input
      set.column(lp_model, n+1, c(1,  x[k,]), c(0:m))
      # Right Hand Side - Outputs set to K'th output
      set.rhs(lp_model, y[k,], c((m+1):(m+s)))
    } else if (orientation == "output"){
      # Output Orientation
      # 1 for theta, theta column Outputs set to - K'th output
      set.column(lp_model, n+1, c(1, -y[k,]), c(0,(m+1):(m+s)))
      # Right Hand Side - inputs set to -K'th input
      set.rhs(lp_model, -x[k,], c(1:m))
    } else {
      stop("Only orientation=input or output supported", call. = FALSE)
    }

    #
    # Super Efficiency Calculation
    #
    # If super efficiency - need to add an extra constraint for K'th value = 0
    if(super){
      set.row(lp_model,                (m + s + 2), 1, c(k))   # Theta forced to 1
      set.constr.type(lp_model, "=",  c(m + s + 2) )
      set.rhs(lp_model, 0,            c(m + s + 2) )
      
      if (k == 1 && debug >= 2) lp_debug_file(lp_model, paste0("super effeciency constraint added"))
    }
    
    if(debug >= 2) 
      lp_debug_file(lp_model, paste0("Final constraints added k= ",k)) # print for every K 
    
    #
    # Solve LP model
    #
    # todo: Check for failure in solving model
    status <- solve(lp_model)  
    if(status != 0){
#       if (status != 2){
#         stop("Solver Phase 1 returned status:", status, call. = FALSE)
#       } else {
      msg <- lp_solve_error_msg(status)
      warning("Solver Phase 1 returned status: ", status, " ", msg, " - setting objective to NaN",
                call. = FALSE)
        phase1_eff    <-  NaN
        efficiency[k] <-  phase1_eff
        lambda[k,]    <-  rep(NaN, n)
#      }
    } else {
      # Normal Results from model    

      # Extract Phase 1 results for DMU K
      # Need to do before slack minimization
      phase1_eff    <-  get.objective(lp_model)           # Store theta in efficiency 
      efficiency[k] <-  phase1_eff
      lambda[k,]    <-  get.variables(lp_model)[1:n]      # Only put lambdas in, leave out theta
    }
    
    # If dual option selected, extract dual weights
    if (dual){
      if ( orientation == "output" ) sign <- -1 else sign <- 1
      x_weight[k,] <-  sign * get.dual.solution(lp_model)[ 2:(m+1)] 
      y_weight[k,] <-  sign * get.dual.solution(lp_model)[(m+2):(m + s + 1)]
    }
        
    #
    # Phase 2 - slack maximization or secondary objective
    #
    
    # Setup for slack or secondary objective is almost the same
    if( !is.nan(phase1_eff) && (slack || ( second != "none"))){
        
      # Round to effecient
       if((abs(phase1_eff-1) < epsilon))
         phase1_eff <- 1
  
      # Add constraint that effeciency must equal efficiency found in Phase 1
      set.row(lp_model,                (m + s + 3), 1, c(n+1))   # Theta forced to 1
      set.constr.type(lp_model, "=",  c(m + s + 3) )
      set.rhs(lp_model, phase1_eff,   c(m + s + 3) )

      if(slack){
        # Maxamize Slacks  - set theta to 0, slacks to -1
        if ( orientation == "output" ) sign <- 1 else sign <- -1
        set.objfn(lp_model, c(0, rep(sign, m + s)), c( (n + 1) : (n + 1 + m + s) ) )
      }

      if(second != "none"){
        if ( second == "min" ) sign <- -1 else sign <- 1
#       WARN: z is a dataframe, need to extract row of data
        set.objfn(lp_model, sign * z[,1], c(1:n))
      }       
            
       if(k == 1 && debug >= 2) lp_debug_file(lp_model, "Phase 2 - Pre slack maximizations")  

       status <- solve(lp_model)  
       if(status != 0){
         #          if (status != 2){
         #            stop("Solver Phase 2 returned status:", status, call. = FALSE)
         #          } else {
         
         msg <- lp_solve_error_msg(status)
         warning("Solver Phase 2 returned status: ", status, " ", msg, " - setting objective to NaN",
                 call. = FALSE)
         lambda[k,]    <- rep(NaN, n)
         if(slack){
           x_slack[k,] <- rep(NaN, m)
           y_slack[k,] <- rep(NaN, s)
         }
       } else {
         # Normal Solver Status   
         # Update lambdas
         lambda[k,]    <- get.variables(lp_model)[1:n]      # Only put lambdas in, leave out theta
         if(slack){
           x_slack[k,] <- get.variables(lp_model)[ (n+2)   : (n + 1 + m)]
           y_slack[k,] <- get.variables(lp_model)[ (n+m+2) : (n + 1 + m + s)]
         }
       } # End check solver      
    } # End secondary solve - phase 2
    
    
#     rm(lp_model)        # Delete lp model
#     delete.lp(lp_model)
    
  } # End of DMU Loop
  
  
  # Cleanup & return Results
  # WARN: Rounding turned off by default. To match other DEA packages, turn rounding on,
  # but does hide some numerical results.

  # Round effeciencies that are within epilson to 1 to 1.
  if(normalize && orientation == "output"){
    efficiency <- 1 / efficiency
  }
  
  if(round){
    # Round Effeciency to 0 and 1
    efficiency[abs(efficiency-1) < epsilon] <- 1
    efficiency[abs(efficiency)   < epsilon] <- 0
  }

  results_l = list(eff=efficiency, lambda=lambda)
  
  if (dual)
    results_l = c(results_l, list(vx=x_weight, uy=y_weight))
  
  if(slack)
    results_l = c(results_l, list(sx=x_slack, sy=y_slack))

  results_l = c(results_l, list(arguments=arguments))

  class(results_l) <- "DEA"
  return(results_l)
    
}
# End dea function
