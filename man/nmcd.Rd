\name{nmcd}
\alias{nmcd}
\alias{plot.nmcd}
\alias{summary.nmcd}
\title{
Non-parametric Multiple Change-points Detection
}
\description{
Detect multiple change-points using  a non-parametric maximum likelihood approach.
}
\usage{
nmcd(x, kmax, cpp, ncp = length(cpp), n)
}
\arguments{
  \item{x}{
data vector
}
  \item{kmax}{
upper bound of the number of change-points
}
  \item{cpp}{
positions of candidate change-points. usually returned by  function\code{isp}.
if this parameter is not assigned, the \code{isp} function will be called to
assign it. 
}
  \item{ncp}{
the number of the candidate change-points.
}
  \item{n}{
length of the data.
}
}
\details{
NMCD use DP algorithm to select change-points, while the true number of change-points
is determined by the Bayesian information criterion(BIC). 
}
\value{
a list with class \code{nmcd} is returned with elements:
\item{npp}{ the true number of change-points }
\item{cpp}{ positions of true change-points }
\item{data}{ raw data, this is not printed on screen by default}
\item{bic}{ minimal BIC value gained.}
}
\references{
Changliang Zou, Guosheng Yin, Long Feng, Zhaojun Wang. Non-parametric Maximum Likelihood Approach to Multiple Change-points Problem
}

\note{
memory consume  may be significant with large data.
}
\examples{
t = nmcd(Blocks)
summary(t)
plot(t)

}
\keyword{ change-points }
