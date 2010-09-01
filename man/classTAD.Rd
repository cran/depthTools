\name{classTAD}
\alias{classTAD}
\title{Weighted Trimmed Mean Distance Classification Method }
\description{Implementation of the classification technique based on assigning each observation to the group that minimizes the trimmed average distance of the given observation to the deepest points of each group in the learning set, weighted by the depth of these points in their own group. 
}

\usage{classTAD(xl,yl,xt,alpha=0)}

\arguments{
  \item{xl}{an \code{nxp} data matrix containing the observations (samples) from the learning set by rows and the variables (genes) by columns }
  \item{yl}{a vector of length \code{n} containing the class each observations in \code{xl} belongs to }
  \item{xt}{an \code{mxp} data matrix containing the observations (samples) from the test set by rows and the variables (genes) by columns }
  \item{alpha}{an optional value for the proportion of observations that are trimmed out when computing the mean. 0 by default.}
}

\value{
  \item{pred }{the vector of length \code{m} containing the predicted class of observations in matrix \code{xt}}
}


\details{This method classifies a given observation \code{x} into one of \code{g} groups, of sizes \code{n1,...,ng}, but taking into account only the \code{m=min\{n1,...,ng\}} deepest elements of each group in the learning set. Additionally, this number can be reduced in a proportion \code{alpha}. The distance of \code{x} to these \code{m} elements is averaged and weighted with the depth of each element with respect to its own group.}

\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}
\references{
Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}
\seealso{classDS}

\examples{
  ## simulated data
  set.seed(0)
  xl <- matrix(rnorm(100),10,10); xl[1:5,]<-xl[1:5,]+1
  yl <- c(rep(0,5),rep(1,5))
  xt <- matrix(rnorm(100),10,10)
  classTAD(xl,yl,xt)
  
  ## real data
  data(prostate)
  set.seed(0)
  learning <- sample(50,40,replace=FALSE) 
  yl <- prost.y[learning]   
  xl <- prost.x[learning,]
  training <- c(1:nrow(prost.x))[-learning]
  yt.real <- prost.y[training]
  xt <- prost.x[training,] 
  yt.estimated <- classTAD(xl,yl,xt)
  yt.real==yt.estimated
}
\keyword{ multivariate }
