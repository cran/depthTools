tmean<-function(x,alpha=0.2)
{
n <- nrow(x);
I <- MBD(x,plotting=FALSE)$ordering;
N <- n-floor(alpha*n)
m1 <- x[I[1:N],]
tm <- apply(m1,2,mean) 
return(list(tm=tm, tm.x=m1))
}
