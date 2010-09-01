\name{tmean}
\alias{tmean}
\title{ Modified Band Depth-Based Alpha Trimmed Mean }
\description{ 
  \code{tmean} computes the mean of the deepest observations within the sample, their depths given by the Modified Band Depth, trimming out the proportion \code{alpha} of the outest observations. 
}

\usage{
tmean(x,alpha=0.2)
}

\arguments{
  \item{x}{ an \code{nxp} data matrix containing the observations (samples) by rows and the variables (genes) by columns }
  \item{alpha}{ the proportion of observations that are trimmed out when computing the mean. 0.2 by default.}
}

\value{
  \item{tm }{the alpha-trimmed mean vector of length p of matrix \code{x}}
  \item{tm.x}{the deespest points of x after removing the proportion \code{alpha} of less deep points }
}


\details{The rows of matrix \code{x}, corresponding to genes, are ordered from center outward, that is, starting with the deepest one(s) and ending with the less deep one(s), according to MBD. The alpha-trimmed mean is computed by first removing the proportion \code{alpha} of less deep points, and then computing the component-wise average of the remaining observations.}
\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}


\examples{
  set.seed(0)
  x <- matrix(rnorm(100),10,10)
  m.x<-apply(x,2,mean)
  t.x<-tmean(x)
  matplot(t(x),type="l",col=1,lty=2,lwd=.8)
  matlines(t(t.x$tm.x), type="l",col=4,lty=2)
  lines(t.x$tm,lty=1,col=2,lwd=1.2)
  legend("bottom",lty=c(2,2,1),col=c(1,4,2),legend=c("whole sample",
         "trimmed sample","trimmed mean"))

}
\keyword{ multivariate }
