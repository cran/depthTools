\name{classDS}
\alias{classDS}
\title{Distance to the Trimmed Mean Classification Method }
\description{Implementation of the classification technique based on assigning each observation to the group that minimizes the distance of the observation to the trimmed mean of the group.
}

\usage{classDS(xl,yl,xt,alpha=0.2)}

\arguments{
  \item{xl}{an \code{nxp} data matrix containing the observations (samples) from the learning set by rows and the variables (genes) by columns }
  \item{yl}{a vector of length \code{n} containing the class each observations in \code{xl} belongs to }
  \item{xt}{an \code{mxp} data matrix containing the observations (samples) from the test set by rows and the variables (genes) by columns }
  \item{alpha}{the proportion of observations that are trimmed out when computing the mean. 0.2 by default.}
}

\value{
  \item{pred }{the vector of length \code{m} containing the predicted class of observations in matrix \code{xt}}
}


\details{This classification method proceeds by first computing the alpha trimmed mean corresponding to each group from the learning set, then computing the distance from a new observation to each trimmed mean. The new sample will then be assigned to the group that minimizes such distance. At the moment, only the Euclidean distance is implemented.}
\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}
\references{
Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}
\seealso{classTAD}

\examples{
  ## simulated data
  set.seed(10)
  xl <- matrix(rnorm(100),10,10); xl[1:5,]<-xl[1:5,]+1
  yl <- c(rep(0,5),rep(1,5))
  xt <- matrix(rnorm(100),10,10)
  classDS(xl,yl,xt)

  ## real data
  data(prostate)
  prost.x<-prostate[,1:100]
  prost.y<-prostate[,101]
  set.seed(1)
  learning <- sample(50,40,replace=FALSE) 
  yl <- prost.y[learning]   
  xl <- prost.x[learning,]
  training <- c(1:nrow(prost.x))[-learning]
  yt.real <- prost.y[training]
  xt <- prost.x[training,] 
  yt.estimated <- classDS(xl,yl,xt)
  yt.real==yt.estimated
  
}
\keyword{ multivariate }
