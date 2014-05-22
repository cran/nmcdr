\name{isp}
\alias{isp}
\alias{plot.isp}
\title{
Initial Screening Procedure for Detecting Multiple Change-points
}
\description{
 Initial screening procedure to exclude most of the inrelevant points based on
 the Cramer-Von Mises two-sample test statistics, 
 thus NMCD is implemented in a much lower dimensional space.
} 
\usage{
isp(x, n = length(x), an = as.integer(log(n)^1.5/2))
}
\arguments{
  \item{x}{
the data vector.
}
  \item{n}{
the length of the data.
}
  \item{an}{
an appropriate integer which is the length of each subsequence of the data.
}
}

\value{
a list is returned with elements:
\item{ncp }{ the number of candidate change-points}
\item{cpp }{ the position of candidate change-points }
}

\examples{
t <- isp(Blocks)
t
plot(t)
}
\keyword{ change-points }
