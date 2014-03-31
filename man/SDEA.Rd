\name{SDEA}
\alias{SDEA}
\title{Super-efficiency Data Envelopment Analysis}

\description{Calculate Data Envelopment Analysis super-efficiency}
\usage{
SDEA(x, y, rts="vrs", orientation="input", slack=FALSE, second="none", z=0, 
  round=FALSE, debug=0) 
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

\item{round}{Optional: \code{round=TRUE} rounds values to 0 and 1 if close.}

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
This function implements super-efficiency model whereby each decision making unit is compared to 
all other decision making units but not itself. This allows for efficiency scores that are "better" 
than 1 for most efficient Decision Making Units.  Note that this can cause infeasibilities when the 
Returns to Scale is not constant.  (Ex. \code{rts=vrs})
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{$eff}{Efficiency scores for each decision making unit}
\item{$lambda}{Lambda scores for each decision making unit}
\item{$sx}{Radial input slacks, only returned when slack=TRUE}
\item{$sy}{Radial output slacks, only returned when slack=TRUE}
}
\references{
Per Andersen and Niels Christian Petersen, "A Procedure for Ranking Efficient Units in Data 
Envelopment Analysis," Management Science,  October 1993 vol. 39 no. 10 1261-1264. Per Andersen 
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
x <- matrix(c(8,2,4,7,2,8,4,3),ncol=2,dimnames=list(LETTERS[1:4]))
colnames(x) = c("X1", "X2")
y <- matrix(c(1,1,1,1),ncol=1,dimnames=list(LETTERS[1:4],"Y"))

# Simple radial DEA efficiency
SDEA(x, y, rts="crs", orientation="input") 

# Simple radial DEA efficiency with slack maximization 
SDEA(x, y, rts="crs", orientation="input", slack=TRUE) 
  
# Example of secondary objective function
x <- matrix(c(8,2,4,7,10,12,2,8,4,3,2,2),ncol=2,dimnames=list(LETTERS[1:6]))
colnames(x) = c("X1", "X2")
y <- matrix(c(1,1,1,1,1,1),ncol=1,dimnames=list(LETTERS[1:6],"Y"))
z <- matrix(c(1:6),ncol=1,dimnames=list(LETTERS[1:6],"Z"))

SDEA (x, y, rts="crs", orientation="input", second="min", z=z)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
