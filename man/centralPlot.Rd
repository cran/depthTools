\name{centralPlot}
\alias{centralPlot}
\title{ Plot of the central curves }
\description{ 
  \code{centralPlot} plots distinctly the [np] most central observations, where [np] is the largest integer smaller than np, and the remaining most external ones, according to the modified band depth. 
}

\usage{
centralPlot(x, p=0.5,col.c='red',col.e='slategray',lty=c(1,3),gradient=FALSE,
            gradient.ramp=NULL,main=NULL,cex=1,...)}

\arguments{
  \item{x}{ a data matrix containing the observations (samples) by rows and the variables (genes) by columns. }
  \item{p}{ proportion of most central samples to be displayed.}
  \item{col.c}{ the color for the central samples, either as a character string or as a number. Ignored if gradient is TRUE.} 
  \item{col.e}{ the color for the external samples.} 
  \item{lty}{ a vector of two components with the line type of the central and external curves.}
  \item{gradient}{ a logical value. If TRUE then the most central curves are plotted with colors according to the gradient.ramp parameter.}
  \item{gradient.ramp}{ an optional vector of two components containing the first and last colors of the palette used to color the most central curves.}
  \item{main}{ a character string for the plot title.}
  \item{cex}{ the magnification to be used for the legend.}
  \item{...}{ further graphical parameters to be passed to 'plot'.}
}

\details{The \code{centralPlot} allows to visualise the most central curves within the dataset.}

\author{ Sara Lopez-Pintado \email{sl2929@columbia.edu} and 
Aurora Torrente \email{etorrent@est-econ.uc3m.es}}
\references{
Lopez-Pintado, S. \emph{et al.} (2010). Robust depth-based tools for the analysis of gene expression data. \emph{Biostatistics}, 11 (2), 254-264.
}
\examples{
  ## simulated data
  set.seed(0)  
  x <- matrix(rnorm(100),10,10)
  centralPlot(x,p=0.2)

  ## real data
  data(prostate)
  prost.x<-prostate[,1:100]
  prost.y<-prostate[,101]
  centralPlot(prost.x[prost.y==0,], p=0.5)  ## 50 % most central normal samples
  centralPlot(prost.x[prost.y==1,], p=0.5, gradient=TRUE, main='Tumor samples')  
                                            ## 50 % most central tumoral samples
}
\keyword{ multivariate }
