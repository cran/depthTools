\name{R.test}
\alias{R.test}
\title{Rank Test Based on the Modified Band Depth }
\description{ 
  \code{R.test} performs the rank test based on the modified band depth, to decide whether two samples come from a single parent distribution. 
}

\usage{
R.test(x,y,n,m,seed=0)
}

\arguments{
  \item{x}{a data matrix containing the observations (samples) by rows and the variables (genes) by columns from the first population }
  \item{y}{a data matrix containing the observations (samples) by rows and the variables (genes) by columns from the second population }
  \item{n}{size of the first sample (less or equal than the number of rows in \code{x})}
  \item{m}{size of the second sample (less or equal than the number of rows in \code{y})}
  \item{seed}{seed to inizialize the random number generation. 0 by default}
}
\value{
  a list containing: 

  \item{p.value }{the p-value of the rank test }
  \item{statistic }{the value of the statistic W of the rank test }
}

\details{Given a population P from which a sample of \code{n} vectors is drawn, and another population P' from which a second sample of \code{m} vectors is obtained, assume there is a third reference sample (from the same population as the largest sample), whose size is also larger than \code{n} and \code{m}. \code{R.test} identifies the largest sample as the one to be split into test and reference samples and verifies if there are enough observations to run the test. Then, the rank test calculates the proportions R and R' of elements from the reference sample whose depths are less or equal than those from the other samples, relative to the reference one, respectively, and order these values from smallest to highest, giving them a rank from 1 to \code{n+m}. The statistic sum of the ranks of values R' (from the second population) has the distribution of a sum of m elements randomly drawn from 1 to \code{n+m} without replacement. }

\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}
\references{
Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}

\examples{
## Rank test for samples from the same population
  x <- matrix(rnorm(100),10,10)
  R.test(x,x,4,4)$p.value

  ## real data
  data(prostate)
  prost.x<-prostate[,1:100]
  prost.y<-prostate[,101]
  normal<-prost.x[prost.y==0,]
  R.test(normal,normal,10,10)$p.value


##  Rank test for samples from different populations
  x <- matrix(rnorm(100),10,10)
  y <- matrix(rnorm(100,5),10,10)
  R.test(x,y,4,4)$p.value

  ## real data
  tumor<-prost.x[prost.y==1,]
  R.test(normal,tumor,10,10)$p.value


}
\keyword{ multivariate }
