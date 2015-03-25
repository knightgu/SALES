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
tau = 0.30
pf = abs(rnorm(p))
# pf = rep(1, p)
pf2 = abs(rnorm(p))
# pf2 = rep(1, p)
w = 2.0
lambda2=1
m1 <-cpernet(y=y,x=x,w=w,tau=tau,eps=1e-14,pf.mean=pf,pf.scale=pf2,
             standardize=F,lambda2=lambda2)

## check KKT condition
lambda.beta = max(abs((t(w*dl2(y-m1$b0[1], 0.5) + 
                         dl2(y-m1$b0[1]-m1$t0[1], tau))) %*% x)/n/pf)
lambda.gamma = max(abs(t(dl2(y-m1$b0[1]-m1$t0[1], tau)) %*% x)/n/pf2)
print(lambda.max <- max(lambda.beta, lambda.gamma))
print(m1$lambda[1])

beta = as.matrix(m1$beta)
theta = as.matrix(m1$theta)

for (l in seq_along(m1$lambda)) {
  r1 = y - x %*% beta[, l] - m1$b0[l]
  r2 = y - x %*% beta[, l] - x %*% theta[, l] - 
    m1$b0[l] - m1$t0[l]
  L = dl2(r2, tau = tau)
  for (j in seq(p)) {
    xlb = crossprod(x[, j], w * r1 + L)
    xlt = crossprod(x[, j], L)
    if (beta[j, l] != 0) {
      AA = -xlb/n + lambda2 * beta[j,l] + pf[j] * m1$lambda[l] * sign(beta[j, l])
      if (abs(AA) > 1e-6) print(paste("AA=", AA,"-l=",l,"-j=",j,sep=""))
    } else {
      BB = abs(xlb/n) - pf[j] * m1$lambda[l]
      if (BB > 0) print(paste("BB=",BB,"-l=",l,"-j=",j,sep=""))
    }
    if (theta[j, l] != 0) {
      CC = -xlt/n + lambda2 * theta[j,l] + pf2[j] * m1$lambda[l] * sign(theta[j, l])
      if (abs(CC) > 1e-6)  print(paste("CC:", CC))
    } else {
      DD = abs(xlt/n) - pf2[j] * m1$lambda[l]
      if (DD > 0) print(paste("DD:", DD))
    }
  }
}




m2 <-cpernet(y=y,x=x,w=w,tau=tau,eps=1e-14,pf.mean=pf,pf.scale=pf2,
             standardize=F,lambda2=lambda2,intercept=F)

## check KKT condition
lambda.beta = max(abs(crossprod(w*dl2(y, 0.5) + dl2(y, tau),x)/n/pf))
lambda.gamma = max(abs(crossprod(dl2(y, tau),x))/n/pf2)
print(lambda.max <- max(lambda.beta, lambda.gamma))
print(m2$lambda[1])

beta = as.matrix(m2$beta)
theta = as.matrix(m2$theta)

for (l in seq_along(m2$lambda)) {
  r1 = y - x %*% beta[, l] - m2$b0[l]
  r2 = y - x %*% beta[, l] - x %*% theta[, l] - 
    m2$b0[l] - m2$t0[l]
  L = dl2(r2, tau = tau)
  for (j in seq(p)) {
    xlb = crossprod(x[, j], w * r1 + L)
    xlt = crossprod(x[, j], L)
    if (beta[j, l] != 0) {
      AA = -xlb/n + lambda2 * beta[j,l] + pf[j] * m2$lambda[l] * sign(beta[j, l])
      if (abs(AA) > 1e-6) print(paste("AA=", AA,"-l=",l,"-j=",j,sep=""))
    } else {
      BB = abs(xlb/n) - pf[j] * m2$lambda[l]
      if (BB > 0) print(paste("BB:", BB))
    }
    if (theta[j, l] != 0) {
      CC = -xlt/n + lambda2 * theta[j,l] + pf2[j] * m2$lambda[l] * sign(theta[j, l])
      if (abs(CC) > 1e-6)  print(paste("CC:", CC))
    } else {
      DD = abs(xlt/n) - pf2[j] * m2$lambda[l]
      if (DD > 0) print(paste("DD:", DD))
    }
  }
}