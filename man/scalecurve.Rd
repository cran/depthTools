\name{scalecurve}
\alias{scalecurve}
\title{ Computation and Representation of the Scale Curve }
\description{ 
  \code{scalecurve} computes the scale curve of a given group, based on the modified band depth, at a given value p as the area of the band delimited by the [np] most central observations, where [np] is the largest integer smaller than np.
}

\usage{
scalecurve(x,y=NULL,xlab="p",ylab="A(p)",main="Scale curve",lwd=2,
           ...)
}

\arguments{
  \item{x}{ a data matrix containing the observations (samples) by rows and the variables (genes) by columns }
  \item{y}{ an optional vector (numeric or factor) of length equal to the number of rows in \code{x}, containing the class of each observation. If unprovided, then all the elements in \code{x} are assumed to belong to a single class }
  \item{xlab}{ label in the x axis} 
  \item{ylab}{ label in the y axis }
  \item{main}{ plot title}
  \item{lwd}{ line widths for the corresponding scale curve(s)}
  \item{...}{ graphical parameters to be passed to 'plot'}
}

\value{
  \item{r}{the value of the scale curve at equidistant values of p, determined by the number of observation within each class. If \code{y} is not provided, then \code{r} is a vector, otherwise is a list with as many components as classes described by \code{y}.}
}

\details{The scale curve measures the increase in the area of the band determined by the fraction p most central curves, where p moves from 0 to 1, thus providing a measure of the sample dispersion. If the data set is represented in parallel coordinates, then the area is computed using the trapezoid formula.}

\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 

Aurora Torrente \email{etorrent@est-econ.uc3m.es}}
\references{
Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}
\examples{
## scale curve of a single data set
  ## simulated data
  set.seed(0)  
  x <- matrix(rnorm(100),10,10)
  scalecurve(x)

  ## real data
  data(prostate)
  scalecurve(prost.x[prost.y==0,])  ## scale curve of normal samples
  scalecurve(prost.x[prost.y==1,])  ## scale curve of tumoral samples
  
## scalecurve of different groups 
  ## simulated data
  x <- matrix(rnorm(100),10,10)
  y <- c(rep("tumoral",5),rep("normal",5))
  scalecurve(x,y)

  ## real data
  labels<-prost.y 
  labels[prost.y==0]<-"normal"; labels[prost.y==1]<-"tumoral"
  scalecurve(prost.x,labels)  
}
\keyword{ multivariate }
