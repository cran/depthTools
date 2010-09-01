\name{MBD}
\alias{MBD}
\title{ Computation of the Modified Band Depth }
\description{ 
  \code{MBD} computes the modified band depth of each observation within a sample which either includes or not the given observation.   
}

\usage{
MBD(x, xRef=NULL, plotting=TRUE, type="l", lty=2, lwd=1, cex=1, 
    col=NULL, cold=NULL, ylim=NULL, colRef=NULL,...)
}

\arguments{
  \item{x}{ a data matrix containing the observations (samples) by rows and the variables (genes) by columns }
  \item{xRef}{ an optional data matrix containing the sample of observations with respect to the modified band depth is computed. If unprovided, then all 
  elements in matrix x are taken into account to compute the depth.}
  \item{plotting}{ logical value. If TRUE then the observations in the data matrix x are plotted, in parallel coordinates}
  \item{type}{ the type of plot to be used, as in 'plot'.} 
  \item{lty}{ the line type for drawing both the data set and the reference set.}
  \item{lwd}{ the line width for both the data set and the reference set. The thickness of the deepest point is increased by 0.5 with respect to the thickest line drawn.}
  \item{cex}{ the amount by which the legend text is magnified.}
  \item{col}{ the color specification for the data set, except for the deepest point. By default, the depth of each point is represented in grayscale colors, with higher intensities corresponding to largest depths.}
  \item{cold}{ the color used to plot the deepest point.}
  \item{ylim}{ numeric vector giving the y coordinates range.}
  \item{colRef}{ the color specification for the reference data set.}
  \item{...}{  graphical parameters (see 'par') and any further arguments of 'plot'}
}

\value{
  a list containing: 

  \item{ordering }{vector giving the ordering of the samples according to their corresponding depths }
  \item{MBD }{vector of the computed depths }
}

\details{The modified band depth is the average proportion of components of the considered observation that are between the corresponding components of 
  all possible pairs of elements in the sample with respect to the depth is computed. The depth is efficiently obtained using the multiplicity of each value in the data matrix ordered by columns rather than exhaustively searching for all pairs of samples.}

\author{ 
Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}
}

\references{
Lopez-Pintado, S. and Romo, J. (2009). On the concept of depth for functional data. \emph{Journal of the American Statistical Association}, 104, 486-503.

Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}
\seealso{scalecurve, R.test}
\examples{
## MBD of all elements within a sample
  
  ## simulated data
  set.seed(0)
  x <- matrix(rnorm(1000),100,10)
  x.depths.1<-MBD(x,plotting=TRUE)
 
  ## real data
  data(prostate)
  normal.depths<-MBD(prost.x[prost.y==0,],plotting=TRUE,
                     main="Normal samples")  
  tumor.depths<-MBD(prost.x[prost.y==1,],plotting=TRUE,
                    main="Tumor samples")

## MBD of a vector with respect to a set of observations

  ## simulated data
  set.seed(0)
  v <- matrix(c(2,1,0,3,-2,1,2,1,0,-2,rnorm(3)),3,5)
  xR <- matrix(rnorm(100),20,5)
  depths<-MBD(v,xR,plotting=TRUE)

  # MBD of normal prostate samples with respect to tumoral ones
  normal.depths<-MBD(prost.x[prost.y==0,],prost.x[prost.y==1,],
                     plotting=TRUE)
}
\keyword{ multivariate }
