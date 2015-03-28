\name{wireless_2012}
\alias{wireless_2012}
\docType{data}
\title{
Data: Wireless Forecasting Data
}
\description{
Technical data from wireless protocol for technology forecasting.
}
\usage{data(wireless_2012)}
\format{
  A data frame containing historical wireless data for technology forecasting.
  \describe{
    \item{\code{Date}}{a numeric vector}
    \item{\code{Bandwidth}}{a numeric vector}
    \item{\code{Bitrate}}{a numeric vector}
}}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\examples{
# Reproduce the published results
data(wireless_2012)
attach(wireless_2012)
x <- data.frame(Bandwidth); row.names(x)  <- rownames(wireless_2012)
y <- data.frame(Bitrate); row.names(y) <- rownames(wireless_2012)
z <- data.frame(Date); row.names(z) <- rownames(wireless_2012)
detach(wireless_2012)
TFDEA(x, y, z, 2001, rts="crs", orientation="input", second="none",
  mode="static", debug=0)
}
\keyword{datasets}
\references{
  D.-J. Lim, T.R. Anderson, J. Kim, Forecast of wireless communication technology: A comparative study of regression and TFDEA Model, in: Technol. Manag. Emerg. Technol., PICMET, Vancouver, Canada, 2012: pp. 1247~1253.
}
