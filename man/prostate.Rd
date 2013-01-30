\name{prostate}
\docType{data}
\alias{prostate}
\title{Gene Expression Data from Tumoral and Normal Prostate Samples and Labels}
\description{Normalized subset from Singh et al. (2002) data included in the \code{prostate} dataset. The raw data comprise the expression of 52 tumoral and 50 non-tumoral prostate samples, obtained using the Affymetrix technology. The data were preprocessed by setting thresholds at 10 and 16000 units, excluding genes whose expression varied less than 5-fold relatively or less than 500 units absolutely between the sample, applying a base 10 logarithmic transformation, and standardising each experiment to zero mean and unit variance across the genes. The 100 most variable genes were selected following the B/W criterion (Dudoit et al. (2002)) and a random selection of 25 normal samples and 25 tumour samples was performed.}
\usage{data(prostate)}
\format{
a 50x101 matrix containing in the first 100 columns the gene expression data of 25 plus 25 randomly selected tumor and normal prostate samples at the 100 most variable genes, selected by the B/W criterion; the last column contains the sample type: 0=normal, 1=tumor.
}
\source{The data are described in Singh et al. (2002).}
\references{ Singh \emph{et al}. (2002). Gene expression correlates of clinincal prostate cancer behavior, \emph{Cancer cell}, 1 (2), 203-209.
 
   Dudoit \emph{et al}. (2002). Comparison of discrimination methods for the classification of tumors using gene expression data, \emph{Journal of the American Statistical Association}, 97 (457), 77-87.
}
\keyword{datasets}
