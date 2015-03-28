\name{DEA}
\alias{DEA}
%- Also NEED an '\alias' for EACH other topic documented here.
%% ToDo
%% Make RTS, common vars a subfile
%%
\title{DEA:  Data Envelopment Analysis}
\description{Estimates a Data Envelopment Analysis frontier and calculates efficiency measures}

\usage{
DEA(x, y, rts="vrs", orientation="input", slack=TRUE, dual=FALSE,
  second="none", z=0, round=FALSE, debug=1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Inputs or resources used by each decision making unit.}

\item{y}{Outputs or products of each decision making unit. Must have same number of rows
as \code{x}}

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

\item{slack}{Optional: \code{slack=TRUE} indicates a secondary objective function of maximizing
radial slacks to identify weakly efficient DMUs}

\item{dual}{Optional: \code{dual=TRUE} reports back the dual weights (multipliers) for the inputs
and outputs}

\item{round}{Optional: \code{round=TRUE} rounds efficiency values to 0 and 1 if close.}

\item{second}{Optional:  Enables an alternate secondary objective function based on lambda and
  the z argument.  The default is \code{none}.  Other options are \code{min} or \code{max}
  which will then minimize or maximize z*lambda while holding efficiency constant for each
  decision making unit. Note that this precludes slack maximization in the current implementation}

\item{z}{Optional:a matrix with one column and the same number of rows (decision making units) as x
  and y, it is only used when \code{second=min} or \code{max}}

\item{debug}{Optional: Only for debugging. If \code{debug} is a integer
  greater then zero debug information will be output.}

}

\details{
This DEA function draws inspiration from previous R packages for doing DEA including Benchmarking
and FEAR.  As such it was designed to use similar parameters and return similar results to allow
users to switch between packages.  The DEA function was developed to support a function for doing
Technology Forecasting using Data Envelopment Analysis or TFDEA.  In particular, TFDEA requires an
option to resolve multiple optima that is similar to but different from the standard slack
maximization approach in Data Envelopment Analysis.  This feature is exposed through DEA function's
second and z parameters.
}

\value{
\item{$status}{If the solver returned a non-zero status for each decision making unit}
\item{$eff}{Efficiency score for each decision making unit}
\item{$lambda}{Lambda values for each decision making unit}

\item{$vx}{Input weights used by each decision making unit, only returned when dual=TRUE}
\item{$uy}{Output weights used by each decision making unit, only returned when dual=TRUE}
\item{$w}{W value for each decision making unit, only returned when dual=TRUE}

\item{$sx}{Radial input slacks, only returned when slack=TRUE}
\item{$sy}{Radial output slacks, only returned when slack=TRUE}

}
\references{

Bogetoft and Otto; Benchmarking with DEA, SFA, and R; Springer 2011

Paul W. Wilson (2008), FEAR 1.0: A Software Package for Frontier Efficiency Analysis with R,
Socio-Economic Planning Sciences 42, 247-254

}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SDEA}} Super-efficiency - an extension to regular DEA that allows for differentiating
between efficient DMUs.

\code{\link{TFDEA}} Technology Forecasting using Data Envelopment Analysis - a method of technology
forecasting using past data to predict future capabilities
}

\examples{

x <- matrix(c(8,2,4,7,2,8,4,3),ncol=2,dimnames=list(LETTERS[1:4]))
colnames(x) = c("X1", "X2")
y <- matrix(c(1,1,1,1),ncol=1,dimnames=list(LETTERS[1:4],"Y"))

# Simple radial DEA efficiency
DEA(x, y, rts="crs", orientation="input")

# Simple radial DEA efficiency with slack maximization
DEA(x, y, rts="crs", orientation="input", slack=TRUE)

# Example of secondary objective function
x <- matrix(c(8,2,4,7,10,12,2,8,4,3,2,2),ncol=2,dimnames=list(LETTERS[1:6]))
colnames(x) = c("X1", "X2")
y <- matrix(c(1,1,1,1,1,1),ncol=1,dimnames=list(LETTERS[1:6],"Y"))
z <- matrix(c(1:6),ncol=1,dimnames=list(LETTERS[1:6],"Z"))

DEA (x,y,rts="crs", orientation="input", round=TRUE, slack=FALSE,
  second="min", z=z)

DEA (x,y,rts="crs", orientation="input", round=TRUE, slack=FALSE,
  second="max", z=z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
