\name{fighter_jet}
\alias{fighter_jet}
\docType{data}
\title{
Data: Fighter Jet technology forecasting data
}
\description{
Technical data from fighter jets for technology forecasting.
}
\usage{data(fighter_jet)}
\format{
  A data frame containing historical jet fighter data for technology forecasting.
  \describe{
    \item{\code{CruiseSpeed}}{a numeric vector}
    \item{\code{FirstFlight}}{a numeric vector}
    \item{\code{Guns}}{a numeric vector}
    \item{\code{InstantaneousTurnRate}}{a numeric vector}
    \item{\code{MaintenanceHoursPerFlightHour}}{a numeric vector}
    \item{\code{MaximumMachNumber}}{a numeric vector}
    \item{\code{MeanFlightHoursBetweenFailure}}{a numeric vector}
    \item{\code{Name}}{a numeric vector}
    \item{\code{NumberOfBVRMissiles}}{a numeric vector}
    \item{\code{NumberOfDogfightMissiles}}{a numeric vector}
    \item{\code{NumberOfSimultaneousTargets}}{a numeric vector}
    \item{\code{Payload}}{a numeric vector}
    \item{\code{RadarRange}}{a numeric vector}
    \item{\code{Range}}{a numeric vector}
    \item{\code{RangeOfBVRMissiles}}{a numeric vector}
    \item{\code{RangeOfDogfightMissiles}}{a numeric vector}
    \item{\code{SeaLevelClimbRate}}{a numeric vector}
    \item{\code{SustainedTurnRate}}{a numeric vector}
    \item{\code{TakeoffRoll}}{a numeric vector}
}}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
O. Inman, T.R. Anderson, R. Harmon, Predicting U.S. jet fighter aircraft introductions from 1944 to 1982: A dogfight between regression and TFDEA, Technol. Forecast. Soc. Change. 73 (2006) 1178~1187.

D.-J. Lim, T.R. Anderson, O.L. Inman, Choosing effective dates from multiple optima in Technology Forecasting using Data Envelopment Analysis (TFDEA), Technol. Forecast. Soc. Change. 88 (2014) 91~97.
}
\examples{
# Reproduce the published results
data(fighter_jet)
attach(fighter_jet)
x <- data.frame(rep(1, 19)); row.names(x) <- Name
y <- data.frame(MeanFlightHoursBetweenFailure, Payload, MaximumMachNumber,
  RangeOfBVRMissiles);
row.names(y) <- Name
z <- data.frame(FirstFlight); row.names(z) <- Name
detach(fighter_jet)
TFDEA(x, y, z, 1960, rts="vrs", orientation="output", mode="dynamic",
  second="min", debug=0)
}
\keyword{datasets}
