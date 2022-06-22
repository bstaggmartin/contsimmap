#relative likelihood of simulation is proportional to dnorm of the seed
#   MINUS the log determinant of the jacobian of any transformations
#   so the problem is that, while this is trivial to calculate for raw simulations, the relative likelihood will shift as you make new
#     traits according to a factor that is not easily calculable in R's framework!
#     could make it so that you could manually pass a derivative function, however...
library(mvnfast)
set.seed(123)
x<-matrix(rnorm(100),ncol=2)
y<-matrix(rnorm(100),ncol=2)
xmu<-rnorm(2)
ymu<-rnorm(2)
xsig<-rWishart(1,2,diag(2))[,,1]
ysig<-rWishart(1,2,diag(2))[,,1]
t.x<-t(t(chol(xsig))%*%t(x)+xmu)
t.y<-t(t(chol(ysig))%*%t(y)+ymu)
sum(dmvn(t.x,xmu,xsig,log=TRUE))-sum(dmvn(t.y,ymu,ysig,log=TRUE))
sum(dnorm(x,log=TRUE))-nrow(x)/2*determinant(xsig)[[1]]-(sum(dnorm(y,log=TRUE))-nrow(y)/2*determinant(ysig)[[1]])
