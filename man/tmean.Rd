\name{tmean}
\alias{tmean}
\title{ Modified Band Depth-Based Alpha Trimmed Mean }
\description{ 
  \code{tmean} computes the mean of the deepest observations within the sample, their depths given by the Modified Band Depth, trimming out the proportion \code{alpha} of the outest observations. 
}

\usage{
tmean(x,alpha=0.2,plotting=FALSE,new=TRUE,cols=c(1,4,8),...)
}

\arguments{
  \item{x}{ an \code{nxd} data matrix containing the observations (samples) by rows and the variables (genes) by columns }
  \item{alpha}{ the proportion of observations that are trimmed out when computing the mean. 0.2 by default.}
  \item{plotting}{ a logical value. If TRUE then a plot is built. If alpha has length 1, then the trimmed mean, the samples used for its computation and the discarded ones are plotted with different colors, according to the values of cols, below. If alpha has length greater than 1, then a plot with several trimmed means is constructed. The first element in cols is used to determine a color palette, from cols[1] (for the smallest value in alpha) to 'gray' (for the greatest value in alpha).}
  \item{new}{ a logical value. If alpha has length 1 or plotting is FALSE, then it is ignored. If TRUE, a new plot is started; otherwise, the new trimmed means are added to the existing plot.}
  \item{cols}{ a vector of length 3 containing, in the following ordering, the colors for depicting the trimmed mean, the trimmed collection of samples and the samples which are not taken into account in the computation of the trimmed mean.}
  \item{...}{ graphical parameters (see 'par') and any further arguments of 'plot'.}
}

\value{
  \item{tm }{the alpha-trimmed mean vector of length p of matrix \code{x}}
  \item{tm.x}{the deespest points of x after removing the proportion \code{alpha} of less deep points }
}


\details{The rows of matrix \code{x}, corresponding to genes, are ordered from center outward, that is, starting with the deepest one(s) and ending with the less deep one(s), according to MBD. The alpha-trimmed mean is computed by first removing the proportion \code{alpha} of less deep points, and then computing the component-wise average of the remaining observations.}
\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}


\examples{
  set.seed(50)
  x <- matrix(rnorm(100),10,10)
  m.x<-apply(x,2,mean)
  t.x<-tmean(x,plotting=TRUE, lty=1)

  t.x.seq <- tmean(x,alpha=c(0,0.25,0.5,0.75),plotting=TRUE, lty=1, cols=2)
}
\keyword{ multivariate }
