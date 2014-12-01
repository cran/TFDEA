\name{SDEA}
\alias{SDEA}
\title{Super-efficiency Data Envelopment Analysis}

\description{Calculate Data Envelopment Analysis super-efficiency}
\usage{
SDEA(x, y, rts="vrs", orientation="input", slack=TRUE, dual = FALSE,
  cook=FALSE,
  second="none", z=0,
  round=FALSE, debug=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Inputs or resources used by each decision making unit}

\item{y}{Outputs or products of each decision making unit}

\item{rts}{Returns to scale for the application, production technology, or
industry studied
  \tabular{rl}{

  vrs \tab Variable returns to scale, convexity and free disposability \cr

  drs \tab Decreasing returns to scale, convexity, down-scaling and free disposability \cr

  crs \tab Constant returns to scale, convexity and free disposability \cr

  irs \tab Increasing returns to scale, (up-scaling, but not down-scaling), convexity and
  free disposability \cr
  }
}

\item{orientation}{Orientation of the DEA model - primary emphasis on input-reduction \code{input}
  or output-augmentation \code{output}}

\item{slack}{Optional: \code{slack=TRUE} indicates a secondary objective function of
  maximizing non-radial slacks}

\item{dual}{Optional: \code{dual=TRUE} reports back the dual weights (multipliers) for the inputs
and outputs}

\item{cook}{Optional: \code{cook=TRUE} enables using Cook algo. to compute super efficiency of DMUs that are infeasible under the standard model.}

\item{round}{Optional: \code{round=TRUE} rounds efficiency values to 0 and 1 if close.}

\item{second}{Optional: Enables an alternate secondary objective function based on lambda and
  the z matrix.  The default is \code{second=none}. Other options include
  \code{min} or \code{max} which will then minimize or maximize z*lambda while holding efficiency
  constant for each decision making unit. Note that this precludes slack maximization in the
  current implementation}

\item{z}{Optional: a matrix with one column and the same number of rows (decision making
  units) as x and y, it is required when \code{second=min} or \code{max}}

\item{debug}{Optional: Only for debugging. If \code{debug} is a integer
  greater then zero debug information is output.}

}

\details{
This function implements either the standard super-efficiency model or the Cook's extended model whereby each decision making unit is compared to all other decision making units but not itself. This allows for efficiency scores that are "better" than 1 for most efficient Decision Making Units.  Note that the standard super-efficiency model can cause infeasibilities especially when the Varialbe Returns to Scale is assumed. (Ex. \code{rts=vrs})
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use

\item{$status}{If the solver returned a non-zero status for each decision making unit}
\item{$eff}{Efficiency score for each decision making unit}
\item{$lambda}{Lambda values for each decision making unit}

\item{$se.eff}{Tau(Gamma) in input(output)-orientated model, only returned when cook=TRUE}
\item{$se.excess}{Beta(Delta) in input(output)-oriented model, only returned when cook=TRUE}

\item{$vx}{Input weights used by each decision making unit, only returned when dual=TRUE}
\item{$uy}{Output weights used by each decision making unit, only returned when dual=TRUE}
\item{$w}{W value for each decision making unit, only returned when dual=TRUE}

\item{$sx}{Radial input slacks, only returned when slack=TRUE}
\item{$sy}{Radial output slacks, only returned when slack=TRUE}


}
\references{
P. Andersen, N.C. Petersen, A Procedure for Ranking Efficient Units in Data Envelopment Analysis, Manage. Sci. 39 (1993) 1261~1264.

W.D. Cook, L. Liang, Y. Zha, J. Zhu, A modified super-efficiency DEA model for infeasibility, J. Oper. Res. Soc. 60 (2009) 276~281.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

\code{\link{DEA}} Data Envelopment Analysis - provides a variety of standard DEA models for
examining the efficiency of different decision making units.

\code{\link{TFDEA}} Technology Forecasting Using Data Envelopment Analysis - a method of
technology forecasting using past data to predict future capabilities

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- matrix(c(1,2,3,4),ncol=1,dimnames=list(LETTERS[1:4],"X"))
y <- matrix(c(1,3,4,3),ncol=1,dimnames=list(LETTERS[1:4],"Y"))

# Simple radial DEA super-efficiency
SDEA(x, y, rts="vrs", orientation="input")

# Simple radial DEA super-efficiency with slack maximization
SDEA(x, y, rts="vrs", orientation="input", slack=TRUE)

# Cook's super-efficiency with slack maximization
SDEA(x, y, rts="vrs", orientation="input", slack=TRUE, cook=TRUE)

# Example of secondary objective function
x <- matrix(c(8,2,4,7,10,12,2,8,4,3,2,2),ncol=2,dimnames=list(LETTERS[1:6]))
colnames(x) = c("X1", "X2")
y <- matrix(c(1,1,1,1,1,1),ncol=1,dimnames=list(LETTERS[1:6],"Y"))
z <- matrix(c(1:6),ncol=1,dimnames=list(LETTERS[1:6],"Z"))

SDEA (x, y, rts="crs", orientation="input", slack=FALSE, second="min", z=z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
