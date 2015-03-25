path = "~/Downloads/RunSimulations/package/spals"

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[Rr]$")) {
    if(trace) cat(nm,": ")           
    source(file.path(path, nm), ...)
    if(trace) cat("sourced\n")
  }
}

dyn.load(file.path(path, "alslassoNET.so"))
dyn.load(file.path(path, "cpalslassoNET.so"))
sourceDir(path)
library(Matrix)

dl2 <- function(r, tau) 2 * abs(tau - (r < 0)) * r

n=100
p=400
# n=4
# p=8
x=matrix(rnorm(n*p),n,p)
y=rnorm(n)
tau = 0.90
pf = abs(rnorm(p))
# pf = rep(1, p)
pf2 = abs(rnorm(p))
# pf2 = rep(1, p)
lambda2=1
m1 <-ernet(y=y,x=x,tau=tau,eps=1e-14,pf=pf,pf2=pf2,
          standardize=F,lambda2=lambda2)

## check KKT condition
lambda.max = max(abs(crossprod(dl2(y-m1$b0[1],tau),x))/n/pf)
m1$lambda[1]
print(lambda.max)

B <- as.matrix(m1$beta)

for (l in 1:length(m1$lambda))
{
  ri <- y - (x%*%B[,l]+m1$b0[l])
  for (j in 1:p)
  {
    L = dl2(ri,tau)
    xl <- crossprod(x[,j], L)
    if(B[j,l]!=0)
    {
      AA<- - xl/n + lambda2*B[j,l]*pf2[j]+ pf[j]*m1$lambda[l]*sign(B[j,l]) 
      if(abs(AA)>=1e-6) print(AA)
      
    }
    else
    {
      BB <- abs(-xl/n)-m1$lambda[l]*pf[j]
      if (BB > 0) print(paste("BB=",BB,"-l=",l,"-j=",j,sep=""))
    }
  }
}



m2 <-ernet(y=y,x=x,tau=tau,eps=1e-14,pf=pf,pf2=pf2,
           standardize=F,lambda2=lambda2,intercept=F)

## check KKT condition
lambda.max = max(abs(crossprod(dl2(y,tau),x))/n/pf)
m2$lambda[1]
print(lambda.max)

B <- as.matrix(m2$beta)

for (l in 1:length(m2$lambda))
{
  ri <- y - (x%*%B[,l]+m2$b0[l])
  for (j in 1:p)
  {
    L = dl2(ri,tau)
    xl <- crossprod(x[,j], L)
    if(B[j,l]!=0)
    {
      AA<- - xl/n + lambda2*B[j,l]*pf2[j]+ pf[j]*m2$lambda[l]*sign(B[j,l]) 
      if(abs(AA)>=1e-6) print(AA)
      
    }
    else
    {
      BB <- abs(-xl/n)-m2$lambda[l]*pf[j]
      if (BB > 0) print(paste("BB=",BB,"-l=",l,"-j=",j,sep=""))
    }
  }
}
