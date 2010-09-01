MBD<- function(x, xRef=NULL, plotting=TRUE, type="l", lty=2, lwd=1, cex=1, col=NULL, cold=NULL, ylim=NULL, colRef=NULL,...)
{
n <- nrow(x); p <- ncol(x) # n: number of observations (samples);  p: dimension of the data
x <- as.matrix(x)

if (length(xRef)==0) {  ## MBD with respect to the same sample
  if (ncol(x) == 1) {x <- t(x)}
  depth <- matrix(0,1,n)
  ordered.matrix <- x
  if (n>1) {
     for (columns in 1:p) {
        ordered.matrix[,columns] <- sort(x[,columns])
        for (element in 1:n) {
             index.1 <- length(which(ordered.matrix[,columns] < x[element,columns]))
             index.2 <- length(which(ordered.matrix[,columns] <= x[element,columns]))
             multiplicity <- index.2 - index.1
             depth[element] <- depth[element] + index.1 * (n - (index.2)) + multiplicity * (n - index.2 + index.1) + choose(multiplicity,2)
             }   ### end FOR element
        }  ### end FOR columns
     depth <- depth / (p * choose(n,2) )
     } ## end IF
  if (n==1) {deepest <- x; depth <- 0}

  ordering<-order(depth,decreasing=TRUE)
  if (plotting) {
     par(mar=c(8,5,5,5),xpd=FALSE)  	
     Gene.Expression<-t(x[ordering[n:2],])
     if (length(type)<n) type<-rep(type,length.out=n) 
     if (length(lty)<n) lty<-rep(lty,length.out=n) 
     if (length(lwd)<n) lwd<-rep(lwd,length.out=n) 
     if (is.null(col)) {color<-gray(0:(n-1)/n)}  ## end IF no color
     else {if (length(col)<n) color<-rep(col,length.out=n)[n:1]}  ## end ELSE
     if (is.null(cold)) {cold<-2}
     if (is.null(ylim)) {ylim<-range(x)}
     matplot(Gene.Expression, ylim=ylim, type=type, lty=lty, lwd=lwd, col=color[n:2],...)
     lwdd<-max(lwd)+.5
     lines(x[ordering==1,], lty=lty, lwd=lwdd, col=cold)
     par(xpd=TRUE)
     legend("bottom",inset=-0.3, legend="deepest sample",col=cold, lty=lty, lwd=lwdd,cex=cex)
     }  ## end IF plotting
  } ## end IF no reference sample
else {
  xRef <- as.matrix(xRef)
  if (ncol(xRef)!=p) {stop("Dimensions of x and xRef do not match")}
  n0 <- nrow(xRef)
  if (ncol(x) == 1) {x <- t(x)}
  depth <- matrix(0,1,n)
  ordered.matrix <- xRef
  if (n0>1) {
     for (columns in 1:p) {
        ordered.matrix[,columns] <- sort(xRef[,columns])
        for (element in 1:n) {
             index.1 <- length(which(ordered.matrix[,columns] < x[element,columns]))
             index.2 <- length(which(ordered.matrix[,columns] <= x[element,columns]))
             multiplicity <- index.2 - index.1
             depth[element] <- depth[element] + (index.1 + multiplicity ) * (n0 - index.1 - multiplicity) + multiplicity * ( index.1 + (multiplicity-1)/2)
             }   ### end FOR element
        }   ### end FOR columns
     depth <- depth / (p * choose(n0,2) )
     } ## end IF
  if (n==1) {deepest <- x; depth <- 0}

  ordering<-order(depth,decreasing=TRUE)
  if (plotting)  {
  	 if (is.null(ylim)) {ylim<-range(x,xRef)}
  	 if (length(type)<n) type<-rep(type,length.out=n) 
  	 if (length(lty)<n) lty<-rep(lty,length.out=n) 
  	 if (length(lwd)<n) lwd<-rep(lwd,length.out=n) 
  	 if (is.null(colRef)) {colRef<-4}
  	 if (is.null(cold)) {cold<-2}
  	 else {colRef<-colRef[1]}  ###only one color for the reference data set
  	 if (is.null(col)) {color<-gray(0:(n-1)/n)}  ## end IF no color
    else {if (length(col)<n) color<-rep(col,length.out=n)[n:1]}  ## end ELSE
    par(mar=c(8,5,5,5),xpd=FALSE)
    Gene.Expression <- t(xRef)
    matplot(Gene.Expression, type=type, ylim=ylim, lty=2, lwd=lwd, col=colRef,...)
    lwdd<-max(lwd)+.5
    matlines(t(x[ordering[n:2],]), lwd=lwd, col=color[n:2], lty=lty)
    lines(x[ordering==1,],lwd=lwdd, lty=1,col=cold)
    par(xpd=TRUE)
    legend("bottom",inset=-0.35, legend=c("deepest sample","reference set"),col=c(cold,colRef), lty=c(1,2),lwd=c(lwdd,lwd[1]),cex=cex)
  }  ## end IF plotting
}  ## end ELSE
return(list(ordering=ordering,MBD=depth))
}
