\name{TFDEA}
\alias{TFDEA}
\title{Technology Forecasting Using DEA}
\description{Calculate the technology rate of change (ROC) which can then be used for predicting
future product performance, estimate new product release dates, and other purposes.}
\usage{
TFDEA(x, y, dmu_date_rel, date_forecast, rts="vrs", orientation="output",
  second="min", mode="static", segroc=FALSE, debug=1)
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
or output-augmentation \code{output}. Note that unlike the DEA functions, the default
is \code{output} orientation.}

\item{second}{Optional: Enables an alternate secondary objective function based on the product
of lambda and dmu_date_rel.  The default is \code{min}.  Other options include \code{none} or
\code{max} which will then skip this step or minimize date * lambda while holding the efficiency
constant.}

\item{mode}{Declares if the technology forecast is done using \code{static} or \code{dynamic}
frontiers.  A static frontier is where the forecast is made using a fixed date of date_forecast.
A dynamic frontier allows for different based dates on the frontier based on the product of
dmu_date and lambda.  When a dynamic frontier is selected, a secondary objective function should
be specified to avoid problems with multiple optima (\code{second=min} or \code{max})}

\item{segroc}{Uses segmented rate of change if \code{TURE}.}

\item{debug}{Optional: Only for debugging. If \code{debug} is a integer greater then zero
debug information will be output.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{

TFDEA returns sets of values at three points in time, at the initial release of the product,
\code{_rel}, at the forecast data specified by \code{date_forecast}, \code{_cur}, and for
the date in the future when the product is forecast \code{for}.

The function returns a number of values per product (DMU). The standardized efficiency (all
inefficiencies are between 0 and 1, for \code{input} and \code{output} orientation) \code{eff},
and the lambda values, \code{lambda}, are returned.

A rate of technology change \code{roc} is returned for products efficient at release. At current
time a local rate of change and at forecast time an individualized rate of change is
returned - \code{sroc}. If \code{segroc = FALSE} then the \code{sroc} is the average
rate of change and is the same for every product.

Lastly a date for current and forecast is returned, \code{date}. If \code{mode = dynamic} then the
current date is the current date adjusted by what products the product is being compared to. If
\code{static} is used then the \code{date} is the \code{date_forecast} for all products. If the
product release date is in the future then a forecast for the product is returned.

Not all values are calculated for all products at all points in time. For example a a product
released in the past with that is inefficient at release would not have a \code{roc} or \code{sroc}
calculated because the product is not used to calculate the overall rate of technology change.

%%  ~Describe the value returned
%%  If it is a LIST, use
\item{$date_soa}{List of unique release dates for conducting analyses}

\item{$dmu_eff_rel}{Efficiency per product (DMU) at time of release}
\item{$dmu_lambda_rel}{Lambdas per DMU at time of release}

\item{$dmu_eff_cur}{Efficiency per product (DMU) at current time (date of forecast)}
\item{$dmu_lambda_cur}{Lambdas per DMU at current time}
\item{$dmu_date_cur}{Adjusted current date for DMU}
\item{$dmu_roc_cur}{Rate of Change for product from date of release to current date}
\item{$dmu_sroc_cur}{Local rate of change for the product}

\item{$dmu_eff_for}{Superefficiency per product (DMU)}
\item{$dmu_lambda_for}{Lambda per DMU at forecast time}
\item{$dmu_date_for}{Date forecast for product based upon supereffeciency and technology rate
of change}
\item{$dmu_sroc_for}{Individualized rate of change for product}

\item{$roc}{Average Rate of Technology Change}

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
  # Example from Inman (2004) p. 93-104, predicting flash drive introduction dates

  drive <- c("A", "B", "C", "D", "E", "F", "G")

  x           <- data.frame(c(16, 14, 8, 25, 40, 30, 40))
  rownames(x) <- drive
  colnames(x) <- c("Cost")

  y           <- data.frame(c(16, 32, 32, 128, 32, 64, 256))
  rownames(y) <- drive
  colnames(y) <- c("Capacity")


  z           <- data.frame(c(2001, 2002, 2003, 2004, 2001, 2002, 2004))
  rownames(z) <- drive
  colnames(z) <- c("Date_Intro")

  # Calc intro date for products using forecast year 2003
  results <- TFDEA(x, y, z, 2003, rts="vrs", orientation="output", mode="dynamic")

  # Examine what dates are forecast for DMU D & G
  print(results$dmu_date_for)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
