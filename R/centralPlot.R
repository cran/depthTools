centralPlot  <- function(x, p=0.5,col=c('red','gray'),lty=c(1,3),...)
{
    p=1-p
    x<-as.matrix(x)
    if (ncol(x) == 1)  x <- t(x)
    if (length(col)<2) stop("Not enough colors defined")
    if (length(lty)<2) stop ("Not enough line types defined")
    n <- nrow(x)
    I <- MBD(x, plotting = FALSE)$ordering
    N <- n - floor(p * n)
    m1 <- x[I[1:N], ]
    if (N<n){
      m2 <- x[I[(N+1):n], ]
      matplot(t(m2),type="l", lty=lty[2],col=col[2],xlab='',ylab='',ylim=c(min(x),max(x)),...)
      matlines(t(m1),lty=lty[1],col=col[1],...)
      } 
    else matplot(t(m1),type="l", lty=3,col='gray',xlab='',ylab='',ylim=c(min(x),max(x)))
    }
