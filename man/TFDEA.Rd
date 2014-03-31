\name{TFDEA}
\alias{TFDEA}
\title{Technology Forecasting Using DEA}
\description{Calculate the technology rate of change (ROC) which can then be used for predicting 
future product performance, estimate new product release dates, and other purposes.}
\usage{
TFDEA(x, y, dmu_date_rel, date_forecast, rts="vrs", orientation="output",second="min", 
  mode="static", debug=0) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{Inputs or resources used by each decision making unit}

\item{y}{Outputs or products of each decision making unit}

\item{dmu_date_rel}{Date of introduction (release) for decision making unit or product}

\item{date_forecast}{Date of forecast}

\item{rts}{Returns to scale for the application, production technology, or industry studied

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

\item{second}{Optional: Enables an alternate secondary objective function based on the product 
  of lambda and dmu_date_rel.  The default is \code{min}.  Other options include \code{none} or 
  \code{max} which will then skip this step or minimize date * lambda while holding the efficiency 
  constant.}

\item{mode}{Declares if the technology forecast is done using \code{static} or \code{dynamic} 
  frontiers.  A static frontier is where the forecast is made using a fixed date of date_forecast. 
  A dynamic frontier allows for different based dates on the frontier based on the product of 
  dmu_date and lambda.  When a dynamic frontier is selected, a secondary objective function should 
  be specified to avoid problems with multiple optima (\code{second=min} or \code{max})}

\item{debug}{Optional: Only for debugging. If \code{debug} is a integer
  greater then zero debug information will be output.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{$date_soa}{List of unique release dates for conducting analyses}

\item{$dmu_eff_rel}{Efficiency at time of release}
\item{$dmu_eff}{Efficiency at time of forecast (date_forecast)}

\item{$dmu_roc}{Rate of Change for each product}

\item{$dmu_date_cur}{Date that each DMU/product uses for comparison}
\item{$dmu_eff_for}{Efficiency of each DMU/product for forecasting}
\item{$dmu_date_for}{Date that each DMU/product uses for forecasting}
\item{$roc}{Average Rate of Change}

%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DEA}} 

\code{\link{SDEA}}  

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
xt <- matrix(c(12,10,8,25),ncol=1,dimnames=list(LETTERS[1:4],"Y"))
yt <- matrix(c(16,32,32,128),ncol=1,dimnames=list(LETTERS[1:4],"Y"))
date_rel <- matrix(c(2003, 2004, 2005, 2005),ncol=1,dimnames=list(LETTERS[1:4],"t"))

tf <- TFDEA(x=xt, y=yt, dmu_date_rel=date_rel, date_forecast=2005, rts="crs")

1/tf$dmu_eff  # Displays identical results Inman, Table 13, pg. 67
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
