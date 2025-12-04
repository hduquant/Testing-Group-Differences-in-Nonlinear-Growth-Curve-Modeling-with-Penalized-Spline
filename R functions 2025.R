 
dat <- read.csv("/NA simulated.csv")  
##################################################################################
# Frequentist B-splines function
##################################################################################
require(nlme)
require(Matrix)
ctrl = lmeControl(msMaxIter=500,msMaxEval=9999)

f.bspline <- function(x=dat$x,
                    group=dat$group,
                    y=dat$y,
                    id=dat$id,
                    n.basis= K, #the total number of internal knots and two boundary knots
                    x.test,  # supply the testing points direct 
                    plot = T,  # plot the group-specific trajectories 
                    control =ctrl){
  
  range.x <- range(x) 
  n.test<-length(x.test)
  
  S <- spfa:::bspl(x = x, n_basis =n.basis, order = 2, 
                   lwr = range.x[1], upr = range.x[2])  # linear B-spline basis
  Q <- diag(n.basis)
  Q[-(1:2), ] <- spfa:::diff_mat(n.basis, 2)  # matrix Q (2nd diff; YL 05/26/24)
  Sa <- S %*% solve(Q)  # transformed basis matrix (YL 05/26/24)
  X1 <- Sa[, 1:2]  # for beta1 & beta2, no penality (YL 05/26/24)
  Z1 <- Sa[, -(1:2)]  # split into X and Z similar to Suk's code # beta2-beta1,... penlaity
  all <- rep( 1, nrow(dat) )  # define all as a column of ones
  
  #minimize the variance of beta difference (Sa)
  fit1 <- lme(fixed = y ~ X1 + group:X1 - 1, random = list( #?
    all = pdIdent(~ Z1 - 1), #penalize this variance for group 0
    all = pdIdent(~ Z1:group - 1),  #penalize this variance for (group1-group0)
    id = pdSymm(~ X1 - 1),  # individual differences in intercept (YL 05/26/24)
    id = pdIdent(~ Z1 - 1) ), # individual differences in other coefficients (YL 03/13/24)
    control=ctrl) 
  
  # (inverse) variance-covariance matrix
  sigma2.e <- fit1$sigma^2
  vc <- lapply(pdMatrix(fit1$modelStruct$reStruct), "*", sigma2.e)
  sigma2.u0 <- vc[1]$all[1, 1]
  sigma2.u1 <- vc[2]$all[1, 1]
  SIGMA <- vc[3]$id   # (YL 05/26/24)
  sigma2.v <- vc[4]$id[1, 1]
  C1 <- cbind(X1, X1 * dat$group, Z1, Z1 * dat$group)
  uid <- sort( unique(dat$id) )  # sorted unique id's
  nobsn <- length(uid)  # (YL 03/02/24)
  for ( i in seq_along(uid) )  # append design columns for persons
    C1 <- cbind(C1, Diagonal(nrow(dat), dat$id == uid[i]) %*% X1 ) # (YL 03/13/24)
  for ( i in seq_along(uid) )  # append design columns for persons
    C1 <- cbind(C1, Diagonal(nrow(dat), dat$id == uid[i]) %*% Z1 ) # (YL 03/13/24)
  C1 <- drop0(C1)  # use sparse matrix
  G1inv <- bdiag(
    Diagonal(n.basis - 2, 1 / sigma2.u0),  # u terms in g_0
    Diagonal(n.basis - 2, 1 / sigma2.u1),  # u terms in g_1
    kronecker(diag( length(uid) ), solve(SIGMA) ), # vcov for delta terms (YL 05/26/24)
    Diagonal( (n.basis - 2) * nobsn, 1 / sigma2.v )  # v terms (YL 03/13/24)
  )
  D1 <- bdiag(
    Matrix(0, 4, 4, doDiag = F),  # fixed effects (YL 05/26/24)
    G1inv  # random effects
  )
  V1inv <- (1 / sigma2.e) * crossprod(C1) + D1  # inverse cov matrix
  
  # testing
  # initialize
  beta.hat <- fixef(fit1) # group 0 beta 1 and group1-group0 beta1
  re <- ranef(fit1)
  u.hat <- unlist(re[1:2]) # basis deviation
  S.test <- spfa:::bspl(x = x.test, n_basis = n.basis, order = 2, 
                        lwr = range.x[1], upr = range.x[2])  # B-spline basis
  Sa.test <- S.test %*% solve(Q)  # transform (YL 05/26/24)
  X1.test <- Sa.test[, 1:2]  # (YL 05/26/24)
  Z1.test <- Sa.test[, -(1:2)] # (YL 05/26/24)
  # difference
  d.hat <- cbind(X1.test * 0, X1.test) %*% beta.hat + 
    cbind(Z1.test * 0, Z1.test) %*% u.hat
  
  Cd.test <- cbind(
    X1.test * 0, X1.test, Z1.test * 0, Z1.test, 
    matrix(0, n.test, n.basis * nobsn) )
  Cd.test <- drop0(Cd.test)
  se.d <- sqrt( diag( Cd.test %*% solve( V1inv, t(Cd.test) ) ) )
  
  summary.table<-cbind(d.hat,se.d,d.hat-1.96*se.d,d.hat+1.96*se.d)
  
  # plotting (YL 04/01/24)
  if (plot)
  {
    # mean curves
    plot(x, dat$y, col = "gray", type = 'n', xlab = "x", ylab = "y", 
         xlim = range(x), ylim = range(dat$y), main = "Fitted Mean (B-spline)")
    g0.hat <- cbind(X1.test, X1.test * 0) %*% beta.hat + 
      cbind(Z1.test, Z1.test * 0) %*% u.hat   # mean curve, group 0
    g1.hat <- cbind(X1.test, X1.test) %*% beta.hat + 
      cbind(Z1.test, Z1.test) %*% u.hat  # mean curve, group 1
    lines(x.test, g0.hat, type = 'l', lwd = 2, col = "red")
    lines(x.test, g1.hat, type = 'l', lwd = 2, col = "blue")
    # 95% confidence band for group 0
    C0.test <- cbind(
      X1.test, X1.test * 0, Z1.test, Z1.test * 0,
      matrix(0, n.test, n.basis * nobsn) )
    C0.test <- drop0(C0.test)
    se0 <- sqrt( diag( C0.test %*% solve( V1inv, t(C0.test) ) ) )
    lines(x.test, g0.hat + qnorm(0.975) * se0, lty = 3, col = "red")
    lines(x.test, g0.hat - qnorm(0.975) * se0, lty = 3, col = "red")
    # 95% confidence band for group 1
    C1.test <- cbind(
      X1.test, X1.test, Z1.test, Z1.test,
      matrix(0, n.test, n.basis * nobsn) )
    C1.test <- drop0(C1.test)
    se1 <- sqrt( diag( C1.test %*% solve( V1inv, t(C1.test) ) ) )
    lines(x.test, g1.hat + qnorm(0.975) * se1, lty = 3, col = "blue")
    lines(x.test, g1.hat - qnorm(0.975) * se1, lty = 3, col = "blue")
    legend("topleft", c("group = 0", "group = 1"), col = c("red", "blue"), text.col = c("red", "blue"),
           bty = 'n', lty = 1)
  }
  
  return(summary.table)
  
}

f.bspline.result=f.bspline(x=dat$x,group=dat$group,y=dat$y,
        id=dat$id,n.basis=11, x.test=c(1,3,5,7,9,11,13,15,17,19,21),
        plot = T,control=ctrl)

f.bspline(x=dat$x,group=dat$group,y=dat$y,id=dat$id,n.basis=11, 
          x.test=c(1,3,5,7,9,11,13,15,17,19,21),plot = T,control=ctrl)

##################################################################################
# Frequentist piecewise power function
##################################################################################
get.pl <- function(
    x,      # data 
    lower,  # lower bound (numeric)
    upper,  # upper bound (numeric)
    k.int   # number of interior knots 
){
  # initialize
  k.all <- k.int + 2  # number of all knots
  nx <- length(x)  # number of data points
  space <- (upper - lower) / (k.int + 1)  # space between knots
  knots <- seq(lower, upper, space)  # all knots
  
  # main loop
  s <- matrix(NA, nrow = nx, ncol = k.all)
  for ( i in seq_len(nx) )
  {
    for ( j in 1:(k.int + 1) )
    {
      if (x[i] >= knots[j] && x[i] <= knots[j + 1] )
      {
        s[i, 1:(j - 1)] <- space
        s[i, j] <- x[i] - knots[j]
        s[i, (j + 1):k.all] <- 0
      }
    }
  }
  s <- s[, -k.all]
  s
}

require(nlme)
require(Matrix)
ctrl = lmeControl(msMaxIter=500,msMaxEval=9999)

f.piecewise <- function(x=dat$x,
                      group=dat$sex2,
                      y=dat$y,
                      id=dat$id,
                      n.basis= K2,  # (YL 03/13/24)
                       x.test,  # supply the testing points direct (YL 04/01/24)
                       plot = T,  # plot the group-specific trajectories (YL 04/01/24)
                       control =ctrl){
  
  range.x <- range(x) 
  n.test<-length(x.test)
  # fitting
  B1 <- get.pl(x = x, lower = range.x[1], upper = range.x[2],
               k.int = n.basis - 2)  # piecewise linear basis
  B <- cbind(1, B1)  # add intercept
  Q <- diag(n.basis)
  Q[-(1:2), -1] <- spfa:::diff_mat(n.basis - 1, 1)  # matrix Q (YL 02/29/24)
  Ba <- B %*% solve(Q)  # transformed basis matrix (YL 05/26/24)
  X3 <- Ba[, 1:2, drop = F]
  Z3 <- Ba[, -(1:2)]  # split into X and Z similar to Suk's code
  all <- rep( 1, nrow(dat) )  # define all as a column of ones
  
  fit3 <- lme(fixed = y ~ X3 + group:X3 - 1,  random = list(
    all = pdIdent(~ Z3 - 1), # penlaize beta2-beta1, beta3-beta2 (no beta0)
    all = pdIdent(~ Z3:group - 1),
    id = pdSymm(~ X3 - 1),  # individual diff in int and 1st coef (YL 05/26/24)
    id = pdIdent(~ Z3 - 1) ) )  # individual diff in other coefficients (YL 03/13/24)
  
  # (inverse) variance-covariance matrix
  sigma2.e <- fit3$sigma^2
  vc <- lapply(pdMatrix(fit3$modelStruct$reStruct), "*", sigma2.e)
  sigma2.u0 <- vc[1]$all[1, 1]
  sigma2.u1 <- vc[2]$all[1, 1]
  SIGMA <- vc[3]$id
  sigma2.v <- vc[4]$id[1, 1]
  C3 <- cbind(X3, X3 * dat$group, Z3, Z3 * dat$group)
  uid <- sort( unique(dat$id) )  # sorted unique id's
  nobsn <- length(uid)  # (YL 03/02/24)
  for ( i in seq_along(uid) )  # append design columns for persons
    C3 <- cbind(C3, Diagonal(nrow(dat), dat$id == uid[i]) %*% X3) # (YL 03/13/24)
  for ( i in seq_along(uid) )  # append design columns for persons
    C3 <- cbind(C3, Diagonal(nrow(dat), dat$id == uid[i]) %*% Z3 )  # (YL 03/13/24)
  C3 <- drop0(C3)  # use sparse matrix
  G3inv <- bdiag(
    Diagonal(n.basis - 2, 1 / sigma2.u0),  # u terms in g_0
    Diagonal(n.basis - 2, 1 / sigma2.u1),  # u terms in g_1
    kronecker(diag( length(uid) ), solve(SIGMA) ),  # vcov for delta terms (YL 03/13/24)
    Diagonal((n.basis - 2) * nobsn, 1 / sigma2.v)  # v terms (YL 03/13/24)
  )
  D3 <- bdiag(
    Matrix(0, 4, 4, doDiag = F),  # fixed effects (YL 02/29/24)
    G3inv  # random effects
  )
  V3inv <- (1 / sigma2.e) * crossprod(C3) + D3  # inverse cov matrix
  
  # testing and plotting
  # initialize
  beta.hat <- fixef(fit3)
  re <- ranef(fit3)
  u.hat <- unlist(re[1:2])
  B1.test <- get.pl(x = x.test, lower = range.x[1], upper = range.x[2],
                    k.int = n.basis - 2)  # piecewise linear basis
  B.test <- cbind(1, B1.test)  # append intercept
  Ba.test <- B.test %*% solve(Q)  # transform
  X3.test <- Ba.test[, 1:2]
  Z3.test <- Ba.test[, -(1:2)]
  
  # difference
  d.hat <- cbind(X3.test * 0, X3.test) %*% beta.hat + 
    cbind(Z3.test * 0, Z3.test) %*% u.hat
  
  # difference
  Cd.test <- cbind(
    X3.test * 0, X3.test, Z3.test * 0, Z3.test, 
    matrix(0, n.test, n.basis * nobsn) )
  Cd.test <- drop0(Cd.test)
  se.d <- sqrt( diag( Cd.test %*% solve( V3inv, t(Cd.test) ) ) )
  
  summary.table<-cbind(d.hat,se.d,d.hat-1.96*se.d,d.hat+1.96*se.d)
  
  # plotting (YL 04/01/24)
  if (plot)
  {
    # mean curves
    plot(x, dat$y, col = "gray", type = 'n', xlab = "x", ylab = "y", 
         xlim = range(x), ylim = range(dat$y), main = "Fitted Mean (Piecewise)")
    g0.hat <- cbind(X3.test, X3.test * 0) %*% beta.hat + 
      cbind(Z3.test, Z3.test * 0) %*% u.hat   # mean curve, group 0
    g1.hat <- cbind(X3.test, X3.test) %*% beta.hat + 
      cbind(Z3.test, Z3.test) %*% u.hat  # mean curve, group 1
    lines(x.test, g0.hat, type = 'l', lwd = 2, col = "red")
    lines(x.test, g1.hat, type = 'l', lwd = 2, col = "blue")
    # 95% confidence band for group 0
    C0.test <- cbind(
      X3.test, X3.test * 0, Z3.test, Z3.test * 0,
      matrix(0, n.test, n.basis * nobsn) )
    C0.test <- drop0(C0.test)
    se0 <- sqrt( diag( C0.test %*% solve( V3inv, t(C0.test) ) ) )
    lines(x.test, g0.hat + qnorm(0.975) * se0, lty = 3, col = "red")
    lines(x.test, g0.hat - qnorm(0.975) * se0, lty = 3, col = "red")
    # 95% confidence band for group 1
    C1.test <- cbind(
      X3.test, X3.test, Z3.test, Z3.test,
      matrix(0, n.test, n.basis * nobsn) )
    C1.test <- drop0(C1.test)
    se1 <- sqrt( diag( C1.test %*% solve( V3inv, t(C1.test) ) ) )
    lines(x.test, g1.hat + qnorm(0.975) * se1, lty = 3, col = "blue")
    lines(x.test, g1.hat - qnorm(0.975) * se1, lty = 3, col = "blue")
    legend("topleft", c("group = 0", "group = 1"), col = c("red", "blue"), text.col = c("red", "blue"),
           bty = 'n', lty = 1)
  }
  
  return(summary.table)
  
}

f.piecewise.result=f.piecewise(x=dat$x,group=dat$group,y=dat$y,
          id=dat$id,n.basis=11, x.test=c(1,3,5,7,9,11,13,15,17,19,21),
          plot = T,control=ctrl)

f.piecewise(x=dat$x,group=dat$group,y=dat$y,id=dat$id,n.basis=11, 
          x.test=c(1,3,5,7,9,11,13,15,17,19,21),plot = T,control=ctrl)

##################################################################################
# Frequentist truncated power function
##################################################################################
require(nlme)
require(Matrix)
ctrl = lmeControl(msMaxIter=500,msMaxEval=9999)

f.tpbasis <- function(x=dat$x,
                    group=dat$group,
                    y=dat$y,
                    id=dat$id,
                    n.basis= K, #the total number of internal knots and two boundary knots
                    x.test,  # supply the testing points direct 
                    plot = T,  # plot the group-specific trajectories 
                    control =ctrl){
  
  range.x <- range(x) 
  n.test<-length(x.test)
  
  knots <- seq(from = range.x[1], to = range.x[2], 
               length.out = n.basis)  # (YL 02/29/24)
  Tp <- cbind(
    1, x,  # global linear part
    pmax(outer(x, knots[-c(1, n.basis)], `-`), 0)  # truncated linear part
  ) # linear TP basis
  Tp <- unname( Tp)
  X2 <- Tp[, 1:2]
  Z2 <- Tp[, -(1:2)] # beta2 -> betak penalized to 0
  all <- rep( 1, length(x))  # define all as a column of ones
  
  fit2 <- lme(fixed = y ~ X2 + group:X2 - 1, random = list(
    all = pdIdent(~ Z2 - 1), #  #penalize this variance for group 0
    all = pdIdent(~ Z2:group - 1),  #penalize this variance for group 1-group 1
    id = pdSymm(~ X2 - 1), # intercept and slope individual variances
    id = pdIdent(~ Z2 - 1) ) ) # TP slopes individual variances
  
  # (inverse) variance-covariance matrix
  sigma2.e <- fit2$sigma^2
  vc <- lapply(pdMatrix(fit2$modelStruct$reStruct), "*", sigma2.e)
  sigma2.u0 <- vc[1]$all[1, 1]
  sigma2.u1 <- vc[2]$all[1, 1]
  SIGMA <- vc[3]$id
  sigma2.v <- vc[4]$id[1, 1]
  C2 <- cbind(X2, X2 * group, Z2, Z2 * group)
  uid <- sort( unique(id) )  # sorted unique id's
  nobsn <- length(uid)  
  for ( i in seq_along(uid) )  # append design columns for persons
    C2 <- cbind(C2, Diagonal(length(x), id == uid[i]) %*% X2 )  
  for ( i in seq_along(uid) )  # append design columns for persons
    C2 <- cbind(C2, Diagonal(length(x), id == uid[i]) %*% Z2 )  
  C2 <- drop0(C2)  # use sparse matrix
  G2inv <- bdiag(
    Diagonal(n.basis - 2, 1 / sigma2.u0),  # u terms in g_0
    Diagonal(n.basis - 2, 1 / sigma2.u1),  # u terms in g_1
    kronecker(diag( length(uid) ), solve(SIGMA) ),  # vcov for delta terms
    Diagonal((n.basis - 2) * nobsn, 1 / sigma2.v)  # v terms
  )
  D2 <- bdiag(
    Matrix(0, 4, 4, doDiag = F),  # fixed effects 
    G2inv  # random effects
  )
  V2inv <- (1 / sigma2.e) * crossprod(C2) + D2  # inverse cov matrix
  
  # testing and plotting
  # initialize
  beta.hat <- fixef(fit2)
  re <- ranef(fit2)
  u.hat <- unlist(re[1:2])
  Tp.test <- cbind(
    1, x.test,  # global linear part
    pmax(outer(x.test, knots[-c(1, n.basis)], `-`), 0)  # truncated linear part
  ) # linear TP basis
  X2.test <- Tp.test[, 1:2]
  Z2.test <- Tp.test[, -(1:2)]
  
  # difference
  d.hat <- cbind(0 * X2.test, X2.test) %*% beta.hat + 
    cbind(Z2.test * 0, Z2.test) %*% u.hat
  
  # difference
  Cd.test <- cbind(
    X2.test * 0, X2.test, Z2.test * 0, Z2.test, 
    matrix(0, n.test, n.basis * nobsn) )
  Cd.test <- drop0(Cd.test)
  se.d <- sqrt( diag( Cd.test %*% solve( V2inv, t(Cd.test) ) ) )
  
  summary.table<-cbind(d.hat,se.d,d.hat-1.96*se.d,d.hat+1.96*se.d)
  
  # plotting (YL 04/01/24)
  if (plot)
  {
    # mean curves
    plot(x, y, col = "gray", type = 'n', xlab = "x", ylab = "y", 
         xlim = range(x), ylim = range(y), main = "Fitted Mean (TP)")
    g0.hat <- cbind(X2.test, X2.test * 0) %*% beta.hat + 
      cbind(Z2.test, Z2.test * 0) %*% u.hat   # mean curve, group 0
    g1.hat <- cbind(X2.test, X2.test) %*% beta.hat + 
      cbind(Z2.test, Z2.test) %*% u.hat  # mean curve, group 1
    lines(x.test, g0.hat, type = 'l', lwd = 2, col = "red")
    lines(x.test, g1.hat, type = 'l', lwd = 2, col = "blue")
    # 95% confidence band for group 0
    C0.test <- cbind(
      X2.test, X2.test * 0, Z2.test, Z2.test * 0,
      matrix(0, n.test, n.basis * nobsn) )
    C0.test <- drop0(C0.test)
    se0 <- sqrt( diag( C0.test %*% solve( V2inv, t(C0.test) ) ) )
    lines(x.test, g0.hat + qnorm(0.975) * se0, lty = 3, col = "red")
    lines(x.test, g0.hat - qnorm(0.975) * se0, lty = 3, col = "red")
    # 95% confidence band for group 1
    C1.test <- cbind(
      X2.test, X2.test, Z2.test, Z2.test,
      matrix(0, n.test, n.basis * nobsn) )
    C1.test <- drop0(C1.test)
    se1 <- sqrt( diag( C1.test %*% solve( V2inv, t(C1.test) ) ) )
    lines(x.test, g1.hat + qnorm(0.975) * se1, lty = 3, col = "blue")
    lines(x.test, g1.hat - qnorm(0.975) * se1, lty = 3, col = "blue")
    legend("topleft", c("group = 0", "group = 1"), col = c("red", "blue"), text.col = c("red", "blue"),
           bty = 'n', lty = 1)
  }
  
  colnames(summary.table)=c("estimate","se","ci.l","ci.u")
  return(summary.table)
  
}

f.tpbasis.result=f.tpbasis(x=dat$x,group=dat$group,y=dat$y,
        id=dat$id,n.basis=11, x.test=c(1,3,5,7,9,11,13,15,17,19,21),
        plot = T,control=ctrl)

f.tpbasis(x=dat$x,group=dat$group,y=dat$y,id=dat$id,n.basis=11, 
          x.test=c(1,3,5,7,9,11,13,15,17,19,21),plot = T,control=ctrl)

##################################################################################
# Bayesian B-splines function
##################################################################################
#####################
library(splines)
library(rjags)
# Posterior mode
mode.v <- function(v){
  dens.y <- density(v)
  mode <- dens.y$x[order(dens.y$y,decreasing=T)][1]
  return(mode)}

#  HPD: Monte Carlo method
emp.hpd <- function (theta, alpha = 0.95) {
  alpha <- min(alpha, 1 - alpha) 
  n <- length(theta) 
  L.start <- round(n * alpha) 
  theta <- sort(theta) 
  e <- theta[(n - L.start + 1):n] - theta[1:L.start] 
  ind <- which(e == min(e))[1] 
  return(c(theta[ind], theta[n - L.start + ind])) 
}

B.bspline <- function(x=dat$x,
                      group=dat$group,
                      y=dat$y,
                      id=dat$id,
                      n.basis= K, #the total number of internal knots and two boundary knots
                      burnInSteps = 10^4,       # Burn-in period
                      nChains = 2,
                      nIter =10^4,  
                      x.test ) { # supply the testing points direct 
  
  range.x <- range(x,na.rm = T) 
  knots <- seq(from = range.x[1], to = range.x[2], 
               length.out = n.basis)  
 J=length(unique(id))
  t=max(sort(table(id)))
  N=length(x)
  data=cbind(id,group)
  group=aggregate(group ~ id, data=data,FUN = "mean")[,2]
  
  ind=as.numeric(as.character(
    factor(id, levels = unique(id), labels = seq_along(unique(id)) )
  ))
  
  s=matrix(NA,N,n.basis)
  tlist=NULL
  j=1
  for (i in 1:N ) {
    s[i,]=bs(x[i],knots=knots[-c(1,length(knots))],degree=1,
                            Boundary.knots = c(knots[1],knots[length(knots)]),intercept = TRUE) 
  }
  
s.test<-bs(x.test,knots=knots[-c(1,length(knots))],degree=1,Boundary.knots = c(knots[1],knots[length(knots)]),intercept = TRUE)

# Step 1: Model
model = "
model {

for (i in 1:J){ # individual

beta[i,1:K]~dmnorm(mubeta[1:K]+gamma[1:K]*group[i], pre.D[1:K,1:K]) 

}

 for (j in 1:N){# time point
mu[j]<-inprod(beta[ind[j],1:K],s[j,])
y[j] ~ dnorm(mu[j],phi) 
 }

pre.D[1:K,1:K]~dwish(R[1:K,1:K],K) 
l2.Cov[1:K,1:K]<-inverse(pre.D[1:K,1:K]) 

phi ~ dgamma(0.01,0.01)
sigma2 <- 1/phi  

mubeta[1] ~ dnorm (0,0.01) # var=1/0.0001=10^3
mubeta[2] ~ dnorm (0,0.01) 
for (j in 3:K){
mubeta[j] ~ dnorm(2*mubeta[j-1]-mubeta[j-2],tau.phi)
}
tau <- 1/tau.phi
tau.phi ~ dgamma(1,0.005)

gamma[1] ~ dnorm (0,0.01) 
gamma[2] ~ dnorm (0,0.01) 
for (j in 3:K){
gamma[j] ~ dnorm(2*gamma[j-1]-gamma[j-2],tau.phi2)
}
tau2 <- 1/tau.phi2
tau.phi2 ~ dgamma(1,0.005)

#difference
for (q in 1:Q){
diff[q]<-inprod(gamma[1:K],s.test[q,])
}

}
"
# Save model
writeLines( model , con="model.txt" ) #creates a text file object
#------------------------------------------------------------------------------
# Load data 
dataList = list(    # Put the information into a list.
  N = N , 
  J = J ,
  ind = ind,
  y = y ,
  s = s ,
  K =n.basis,
  group=group,
  R=diag(n.basis),
  s.test=s.test,
  Q=length(x.test)
)
#------------------------------------------------------------------------------
# Step 2: Specifying a starting value
# You need to set a seed; otherwise, you will get a different result every time you run the code
# initsList1 = list(  phi=1, mubeta=rep(1,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=1)
# initsList2 = list(  phi=0.5, mubeta=rep(0,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=2)

#------------------------------------------------------------------------------ 
# Step 3: Adaptation and burn-in 
parameters = c( "mubeta","tau","tau2","sigma2","gamma",
                "diff","l2.Cov")      # Specify the estimated parameters 
adaptSteps =100           # Adaptive period
burnInSteps = burnInSteps      # Burn-in period
nChains = nChains
nIter = nIter    # The number of kept iterations
jagsModel = jags.model( "model.txt" , data=dataList , 
                        n.chains=nChains , n.adapt=adaptSteps )
update( jagsModel , n.iter=burnInSteps)
#------------------------------------------------------------------------------ 
# Step 4: Sample from posterior distributions
codaSamples = coda.samples( jagsModel , variable.names=parameters, 
                            n.iter=nIter , thin=1)
diag=gelman.diag(codaSamples,multivariate=FALSE)[[1]][,1]
#------------------------------------------------------------------------------ 
# Step 5: Summarize posterior distributions
mcmcChain = as.matrix( codaSamples)
mcmcChain<-mcmc(  mcmcChain)
mean <- apply(mcmcChain, 2, mean)
median <- apply(mcmcChain, 2, median)
mode <- apply(mcmcChain, 2,mode.v)
sd <- apply(mcmcChain, 2,sd)
qbp <- apply(mcmcChain, 2, function (x) quantile(x, c(0.025,0.975)))
hpd <- apply(mcmcChain, 2,function (x) emp.hpd(x, c(0.025,0.975)))

# Now we can obtain a better summary
Bayes_sum <- rbind(mean,median, mode, sd, qbp, hpd,diag)
rownames(Bayes_sum) <- c("mean", "median", "mode", "sd", "qbp.lower", "qbp.upper",
                         "hpd.lower", "hpd.upper","PSRF")

return(Bayes_sum)
}

B.bspline.result=B.bspline(x=dat$x,group=dat$group,y=dat$y,
          id=dat$id,n.basis=11, x.test=c(1,3,5,7,9,11,13,15,17,19,21),
          burnInSteps = 10^4, nChains = 2,nIter =10^4)
B.bspline.result[,1:11]

##################################################################################
# Bayesian piecewise function
##################################################################################
#####################
library(rjags)
library(dplyr)
# Posterior mode
mode.v <- function(v){
  dens.y <- density(v)
  mode <- dens.y$x[order(dens.y$y,decreasing=T)][1]
  return(mode)}

#  HPD: Monte Carlo method
emp.hpd <- function (theta, alpha = 0.95) {
  alpha <- min(alpha, 1 - alpha) 
  n <- length(theta) 
  L.start <- round(n * alpha) 
  theta <- sort(theta) 
  e <- theta[(n - L.start + 1):n] - theta[1:L.start] 
  ind <- which(e == min(e))[1] 
  return(c(theta[ind], theta[n - L.start + ind])) 
}

B.piecewise <- function(x=dat$x,
                        group=dat$group,
                        y=dat$y,
                        id=dat$id,
                        n.basis= K, #the total number of internal knots and two boundary knots
                        burnInSteps = 10^4,       # Burn-in period
                        nChains = 2,
                        nIter =10^4,  
                      x.test ) { # supply the testing points direct 
  
  range.x <- range(x,na.rm = T) 
  knots <- seq(from = range.x[1], to = range.x[2], 
               length.out = n.basis)  
  J=length(unique(id))
  t=max(sort(table(id)))
  N=length(x)
  data=cbind(id,group)
  group=aggregate(group ~ id, data=data,FUN = "mean")[,2]
  
  ind=as.numeric(as.character(
    factor(id, levels = unique(id), labels = seq_along(unique(id)) )
  ))
  
  s=matrix(NA,N,n.basis)
  space= knots[2]- knots[1]
  for (i in 1:N) {
      for (l in 1:(n.basis-1)){
        if (between(x[i],knots[l],knots[l+1])){
          s[i,1:(l-1)]=space
          s[i,l]=x[i]-knots[l]
          s[i,(l+1):n.basis]=0
        }
      }
    }

  s=s[,-n.basis]
  
  s.test=matrix(NA,length(x.test),length(knots))
  for (j in 1:length(x.test)){
    for (l in 1:(n.basis-1)){
      if (between(x.test[j],knots[l],knots[l+1])){
        s.test[j,1:(l-1)]=space
        s.test[j,l]=x.test[j]-knots[l]
        s.test[j,(l+1):n.basis]=0
      }
    }
  }
  
  s.test =  s.test[,-n.basis]
  
  # Step 1: Model
  model = "
model {

for (i in 1:J){ # individual

beta[i,1:K]~dmnorm(mubeta[1:K]+gamma[1:K]*group[i], pre.D[1:K,1:K]) 

}

 for (j in 1:N){# time point
mu[j]<-beta[ind[j],1]+inprod(beta[ind[j],2:K],s[j,])
y[j] ~ dnorm(mu[j],phi) 
 } 

pre.D[1:K,1:K]~dwish(R[1:K,1:K],K) 
l2.Cov[1:K,1:K]<-inverse(pre.D[1:K,1:K]) 

phi ~ dgamma(.01,.01)
sigma2 <- 1/phi  

mubeta[1] ~ dnorm (0,0.01) 
mubeta[2] ~ dnorm (0,0.01)
for (j in 3:K){
mubeta[j] ~ dnorm(mubeta[j-1],tau.phi) 
}
tau <- 1/tau.phi
tau.phi ~ dgamma(1,0.005)

gamma[1] ~ dnorm (0,0.01) 
gamma[2] ~ dnorm (0,0.01)
for (j in 3:K){
gamma[j] ~ dnorm(gamma[j-1],tau.phi2)
}
tau2 <- 1/tau.phi2
tau.phi2 ~ dgamma(1,0.005)

#difference
for (q in 1:Q){
diff[q]<-gamma[1]+inprod(gamma[2:K],s.test[q,])
}

}
"
# Save model
writeLines( model , con="model.txt" ) #creates a text file object
#------------------------------------------------------------------------------
# Load data 
dataList = list(    # Put the information into a list.
  N = N , 
  J = J ,
  ind = ind,
  y = y ,
  s = s ,
  K =n.basis,
  R=diag(n.basis),
  group=group,
  s.test=s.test,
  Q=length(x.test)
)
#------------------------------------------------------------------------------
# Step 2: Specifying a starting value
# You need to set a seed; otherwise, you will get a different result every time you run the code
# initsList1 = list(  phi=1, mubeta=rep(1,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=1)
# initsList2 = list(  phi=0.5, mubeta=rep(0,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=2)
#------------------------------------------------------------------------------ 
# Step 3: Adaptation and burn-in 
parameters = c( "mubeta","tau","tau2","sigma2","gamma",
                "diff","l2.Cov")      # Specify the estimated parameters    # Specify the estimated parameters 
adaptSteps =100           # Adaptive period
burnInSteps = burnInSteps       # Burn-in period
nChains = nChains
nIter = nIter     # The number of kept iterations
jagsModel = jags.model( "model.txt" , data=dataList , 
                        n.chains=nChains , n.adapt=adaptSteps )
update( jagsModel , n.iter=burnInSteps)
#------------------------------------------------------------------------------ 
# Step 4: Sample from posterior distributions
codaSamples = coda.samples( jagsModel , variable.names=parameters, 
                            n.iter=nIter , thin=1)
diag=gelman.diag(codaSamples,multivariate=FALSE)[[1]][,1]
#------------------------------------------------------------------------------ 
# Step 5: Summarize posterior distributions
mcmcChain = as.matrix( codaSamples)
mcmcChain<-mcmc(  mcmcChain)
mean <- apply(mcmcChain, 2, mean)
median <- apply(mcmcChain, 2, median)
mode <- apply(mcmcChain, 2,mode.v)
sd <- apply(mcmcChain, 2,sd)
qbp <- apply(mcmcChain, 2, function (x) quantile(x, c(0.025,0.975)))
hpd <- apply(mcmcChain, 2,function (x) emp.hpd(x, c(0.025,0.975)))

# Now we can obtain a better summary
Bayes_sum <- rbind(mean,median, mode, sd, qbp, hpd,diag)
rownames(Bayes_sum) <- c("mean", "median", "mode", "sd", "qbp.lower", "qbp.upper",
                         "hpd.lower", "hpd.upper","PSRF")

return(Bayes_sum)
}

B.piecewise.result=B.piecewise(x=dat$x,group=dat$group,y=dat$y,
            id=dat$id,n.basis=11, x.test=c(1,3,5,7,9,11,13,15,17,19,21),
            burnInSteps = 10^4, nChains = 2,nIter =10^4)
B.piecewise.result[,1:11]

##################################################################################
# Bayesian truncated power function
##################################################################################
#####################
library(rjags)
library(dplyr)
# Posterior mode
mode.v <- function(v){
  dens.y <- density(v)
  mode <- dens.y$x[order(dens.y$y,decreasing=T)][1]
  return(mode)}

#  HPD: Monte Carlo method
emp.hpd <- function (theta, alpha = 0.95) {
  alpha <- min(alpha, 1 - alpha) 
  n <- length(theta) 
  L.start <- round(n * alpha) 
  theta <- sort(theta) 
  e <- theta[(n - L.start + 1):n] - theta[1:L.start] 
  ind <- which(e == min(e))[1] 
  return(c(theta[ind], theta[n - L.start + ind])) 
}

B.tpbasis <- function(x=dat$x,
                        group=dat$group,
                        y=dat$y,
                        id=dat$id,
                        n.basis= K, #the total number of internal knots and two boundary knots
                        burnInSteps = 10^4,       # Burn-in period
                        nChains = 2,
                        nIter =10^4,  
                        x.test ) { # supply the testing points direct 
  
  range.x <- range(x,na.rm = T) 
  knots <- seq(from = range.x[1], to = range.x[2], 
               length.out = n.basis)  
  J=length(unique(id))
  t=max(sort(table(id)))
  N=length(x)
  data=cbind(id,group)
  group=aggregate(group ~ id, data=data,FUN = "mean")[,2]
  
  ind=as.numeric(as.character(
    factor(id, levels = unique(id), labels = seq_along(unique(id)) )
  ))
  
  s <- cbind(
    x,  # global linear part
    pmax(outer(x, knots[-c(1, n.basis)], `-`), 0))  # truncated linear part
  
  s.test <- cbind(
    x.test,  # global linear part
    pmax(outer(x.test, knots[-c(1, n.basis)], `-`), 0)  # truncated linear part
  ) # linear TP basis

  
  # Step 1: Model
  model = "
model {

for (i in 1:J){ # individual

beta[i,1:K]~dmnorm(mubeta[1:K]+gamma[1:K]*group[i], pre.D[1:K,1:K]) 

}

 for (j in 1:N){# time point
mu[j]<-beta[ind[j],1]+inprod(beta[ind[j],2:K],s[j,])
y[j] ~ dnorm(mu[j],phi) 
 } 

pre.D[1:K,1:K]~dwish(R[1:K,1:K],K) 
l2.Cov[1:K,1:K]<-inverse(pre.D[1:K,1:K]) 

phi ~ dgamma(.01,.01)
sigma2 <- 1/phi  

mubeta[1] ~ dnorm (0,0.01) 
mubeta[2] ~ dnorm (0,0.01)
for (j in 3:K){
mubeta[j] ~ dnorm(0,tau.phi) 
}
tau <- 1/tau.phi
tau.phi ~ dgamma(1,0.005)

gamma[1] ~ dnorm (0,0.01) 
gamma[2] ~ dnorm (0,0.01)
for (j in 3:K){
gamma[j] ~ dnorm(0,tau.phi2)
}
tau2 <- 1/tau.phi2
tau.phi2 ~ dgamma(1,0.005)

#difference
for (q in 1:Q){
diff[q]<-gamma[1]+inprod(gamma[2:K],s.test[q,])
}

}
"
# Save model
writeLines( model , con="model.txt" ) #creates a text file object
#------------------------------------------------------------------------------
# Load data 
dataList = list(    # Put the information into a list.
  N = N , 
  J = J ,
  ind = ind,
  y = y ,
  s = s ,
  K =n.basis,
  R=diag(n.basis),
  group=group,
  s.test=s.test,
  Q=length(x.test)
)
#------------------------------------------------------------------------------
# Step 2: Specifying a starting value
# You need to set a seed; otherwise, you will get a different result every time you run the code
# initsList1 = list(  phi=1, mubeta=rep(1,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=1)
# initsList2 = list(  phi=0.5, mubeta=rep(0,K2),
#                     .RNG.name="base::Super-Duper", .RNG.seed=2)
#------------------------------------------------------------------------------ 
# Step 3: Adaptation and burn-in 
parameters = c( "mubeta","tau","tau2","sigma2","gamma",
                "diff","l2.Cov")      # Specify the estimated parameters    # Specify the estimated parameters 
adaptSteps =100           # Adaptive period
burnInSteps = burnInSteps       # Burn-in period
nChains = nChains
nIter = nIter     # The number of kept iterations
jagsModel = jags.model( "model.txt" , data=dataList , 
                        n.chains=nChains , n.adapt=adaptSteps )
update( jagsModel , n.iter=burnInSteps)
#------------------------------------------------------------------------------ 
# Step 4: Sample from posterior distributions
codaSamples = coda.samples( jagsModel , variable.names=parameters, 
                            n.iter=nIter , thin=1)
diag=gelman.diag(codaSamples,multivariate=FALSE)[[1]][,1]
#------------------------------------------------------------------------------ 
# Step 5: Summarize posterior distributions
mcmcChain = as.matrix( codaSamples)
mcmcChain<-mcmc(  mcmcChain)
mean <- apply(mcmcChain, 2, mean)
median <- apply(mcmcChain, 2, median)
mode <- apply(mcmcChain, 2,mode.v)
sd <- apply(mcmcChain, 2,sd)
qbp <- apply(mcmcChain, 2, function (x) quantile(x, c(0.025,0.975)))
hpd <- apply(mcmcChain, 2,function (x) emp.hpd(x, c(0.025,0.975)))

# Now we can obtain a better summary
Bayes_sum <- rbind(mean,median, mode, sd, qbp, hpd,diag)
rownames(Bayes_sum) <- c("mean", "median", "mode", "sd", "qbp.lower", "qbp.upper",
                         "hpd.lower", "hpd.upper","PSRF")

return(Bayes_sum)
}

B.tpbasis(x=dat$x,group=dat$group,y=dat$y,id=dat$id,n.basis=11, 
          x.test=c(1,3,5,7,9,11,13,15,17,19,21), burnInSteps = 10^4, nChains = 2,nIter =10^4)



####
f.bspline.result
f.piecewise.result
f.tpbasis.result