\name{isEfficient}
\alias{isEfficient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check if DEA DMU is efficient
}
\description{
Cheks that based upon model orientation a DMU is efficient. For input orientation DMU's >= 0 are efficient, for output orientation DMU's <= 0.
}
\usage{
isEfficient(eff, orientation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eff}{
Efficiency - numeric value.
}
  \item{orientation}{
Orientation is either \code{input} or \code{output}
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Returns a logical value.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
