#######################################################################################
####------------------ Classification variational approximation -------------------####
#######################################################################################
ZILPNMVA_Poisson <- function(X,n.factors,trace,maxit=1000,cv_group=NULL,level=0.95,sd.errors=TRUE) { # I delete V here
  
  n.s <- nrow(X); n.f <- ncol(X)
  
  M <- rowSums(X); sigma_beta <- matrix(0.1,n.f,n.factors); a1 <- 2; a2 <- 2 
  if (is.null(cv_group)) {
    cvsample <- matrix(0,n.s,n.f)
  }else{cvsample <- cv_group}
  
  out.list <- list()
  
  ### Initialization 1
  pzero.col <- apply(X, 2, function(x) {sum(x==0)/n.s})
  ppi <- new.pi <- t(ifelse(t(X)==0, pzero.col, 0))
  gam <- new.gam <- matrix(c(2,2), nrow=n.f, ncol=2, byrow=T)
  sigma <- new.sigma <- matrix(1,n.s,n.factors)
  lambda <- new.lambda <- matrix(0.5,n.f,n.factors)
  factor_coefs_0 <-  new.factor_coefs_0 <-  rep(1,n.f)
  
  X.rc <- scale(log(X+0.05), scale=T, center=T)
  re <- svd(X.rc, n.factors, n.factors)
  factor_coefs_j <- new.factor_coefs_j <- re$v ## this is r_j or beta_j
  if(n.factors==1){
    factor_scores <- new.factor_scores <- re$u * (re$d[1]) # this is m_i or f_i
  }else{factor_scores <- new.factor_scores <- re$u %*% diag(re$d[1:n.factors])}
  
  ### VA iteration
  cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff <- 1e5; eps <- 1e-8; max.iter <- 500
  b.cur.logfunc <- b0.cur.logfunc <- f.cur.logfunc <- -1e6 
  
  while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
    if(trace) cat("Iteration:", iter, "\n")
    
    ####-------------- Define L_ij  --------####
    L_func <- function(s,l,m,r) {
      my.func <- function(i,j) {0.5 * (-sum(log(1-s[i,]*l[j,]+1e-8)) + sum((2*m[i,]+s[i,]*r[j,])*(1/(1-s[i,]*l[j,]+1e-8))*r[j,]) + sum(m[i,]*l[j,]*(1/(1-s[i,]*l[j,]+1e-8))*m[i,]))}
      ll <- matrix(NA, nrow=n.s, ncol=n.f)
      for(i in 1:n.s){
        for(j in 1:n.f){
          ll[i,j] <- my.func(i,j) 
        }  
      }
      return(ll)
    }
    
    ## VLB or ELBO
    VLB <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL,g=NULL) {
      new.factor_coefs_0 <- x[1:n.f]
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      new.gam <- g
      
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + new.factor_scores %*% t(new.factor_coefs_j)
      lll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      y1 <- X*ll*I(cvsample==0)
      y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(lll))))*I(cvsample==0)
      y3 <- ( (a1 - new.gam[,1])*(digamma(new.gam[,1])-digamma(new.gam[,1]+new.gam[,2])) + (a2 - new.gam[,2])*(digamma(new.gam[,2])-digamma(new.gam[,1]+new.gam[,2])) )
      y4 <- sum(lbeta(new.gam[,1],new.gam[,2])) - n.f*lbeta(a1,a2)
      fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}
      fun3 <- function(j) { 0.5 * (sum(log(t(new.lambda)[,j])) - sum((1/t(sigma_beta)[,j])*t(new.lambda)[,j]) - sum((1/sigma_beta[j,])*new.factor_coefs_j[j,]^2))}
      fun4 <- function(i) {
        e1 <- exp(digamma(new.gam[,1]))/exp(digamma(new.gam[,1]+new.gam[,2]))
        e2 <- exp(digamma(new.gam[,2]))/exp(digamma(new.gam[,1]+new.gam[,2]))
        pi.mat <-  (1-(new.pi*I(cvsample==0))[i,])*log(e2/(1-(new.pi*I(cvsample==0))[i,])) + (new.pi*I(cvsample==0))[i,]*log(e1/(new.pi*I(cvsample==0))[i,])
        sum(na.omit(pi.mat))
      }
      y <- sum(y1) + sum(y2) + sum(y3) + y4 + sum(sapply(1:n.s,fun2)) + sum(sapply(1:n.f,fun3)) + sum(sapply(1:n.s,fun4))
      return(y)
    }
    
    ## VLB_rept
    VLB_rept <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL,g=NULL) {
      new.factor_coefs_0 <- x[1:n.f]
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      new.gam <- g
      
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + new.factor_scores %*% t(new.factor_coefs_j)
      lll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      y1 <- X*ll*I(cvsample==1)
      y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(lll))))*I(cvsample==1)
      y3 <- ( (a1 - new.gam[,1])*(digamma(new.gam[,1])-digamma(new.gam[,1]+new.gam[,2])) + (a2 - new.gam[,2])*(digamma(new.gam[,2])-digamma(new.gam[,1]+new.gam[,2])) )
      y4 <- sum(lbeta(new.gam[,1],new.gam[,2])) - n.f*lbeta(a1,a2)
      fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))*I(cvsample==1)}
      fun3 <- function(j) { 0.5 * (sum(log(t(new.lambda)[,j])) - sum(t(sigma_beta)[,j]^(-1)*t(new.lambda)[,j]) - sum(sigma_beta[j,]^(-1)*new.factor_coefs_j[j,]^2))*I(cvsample==1)}
      
      fun4 <- function(i) {
        e1 <- exp(digamma(new.gam[,1]))/exp(digamma(new.gam[,1]+new.gam[,2]))
        e2 <- exp(digamma(new.gam[,2]))/exp(digamma(new.gam[,1]+new.gam[,2]))
        pi.mat <-  (1-(new.pi*I(cvsample==1))[i,])*log(e2/(1-(new.pi*I(cvsample==1))[i,])) + (new.pi*I(cvsample==1))[i,]*log(e1/(new.pi*I(cvsample==1))[i,])
        sum(na.omit(pi.mat))
      }
      y <- sum(y1) + sum(y2) + sum(y3) + y4 + sum(sapply(1:n.s,fun2)) + sum(sapply(1:n.f,fun3)) + sum(sapply(1:n.s,fun4))
      return(y)
    }
    
    ####------------- Define U_ij ----------#####
    U_func <- function(s,l,m,r) {
      my.func <- function(i,j) { (1/(1-s[i,]*l[j,]+1e-8))*r[j,] + l[j,]*(1/(1-s[i,]*l[j,]+1e-8))*m[i,]}
      u <- array(NA, dim=c(n.s,n.f,n.factors))
      for(i in 1:n.s) {
        for(j in 1:n.f) {
          u[i,j,] <- my.func(i,j)
        }
      }
      return(u)
    }
    
    ####--------------- Define T_ij ------------------#####
    T_func <- function(s,l,m,r) {
      my.func <- function(i,j) { (1/(1-s[i,]*l[j,]+1e-8))*m[i,] + s[i,]*(1/(1-s[i,]*l[j,]+1e-8))*r[j,]}
      tt <- array(NA, dim=c(n.s,n.f,n.factors))
      for(i in 1:n.s) {
        for(j in 1:n.f) {
          tt[i,j,] <- my.func(i,j)
        }
      }
      return(tt)
    }
    
    ####------------- Define R_ij ----------#####
    R_func <- function(s,l,m,r) {
      my.func <- function(i,j) (l[j,]*(1-l[j,]*s[i,]+1e-8)+(r[j,]+m[i,]*l[j,])^2)/(1-l[j,]*s[i,]+1e-8)^2
      rr <- array(NA, dim=c(n.s,n.f,n.factors))
      for(i in 1:n.s) {
        for(j in 1:n.f) {
          rr[i,j,] <- my.func(i,j)
        }
      }
      return(rr)
    }
    
    ####------------- Define P_ij ----------#####
    P_func <- function(s,l,m,r) {
      my.func <- function(i,j) (s[i,]*(1-l[j,]*s[i,]+1e-8)+(m[i,]+r[j,]*s[i,])^2)/(1-l[j,]*s[i,]+1e-8)^2
      pp <- array(NA, dim=c(n.s,n.f,n.factors))
      for(i in 1:n.s) {
        for(j in 1:n.f) {
          pp[i,j,] <- my.func(i,j)
        }
      }
      return(pp)
    }
    
    #####-------------------- Update pi -----------------#####
    ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
    alp <- log(M/(rowSums((exp(ll)))))
    #alp <- log(M/(rowSums((1-new.pi)*(exp(ll)))))
    e.mat <- ll + matrix(alp,n.s,n.f)
    sum1 <- exp(-exp(e.mat))*I(cvsample==0)
    for(i in 1:n.s){
      for(j in 1:n.f) {
        new.pi[i,j] <- exp(digamma(new.gam[j,1]))/(exp(digamma(new.gam[j,1]))+(exp(digamma(new.gam[j,2])))*sum1[i,j]+1e-8)
      }
    }
    new.pi[X!=0] <- 0
    
    grad_gamma1 <- function(x,ppi=NULL,j,g=NULL){
      
      gg1 <- (a1-x)*(trigamma(x)-trigamma(x+g[j,2]))
      gg2 <- trigamma(x)*sum(I(ppi[,j]>0.5)*I(cvsample==0)[,j])
      gg3 <- n.s*trigamma(x+g[j,2])
      gg4 <- (a2-g[j,2])*trigamma(x+g[j,2])
      gg <- gg2 - gg3 + gg1 - gg4
      return(gg)
      
    }
    
    for(j in 1:n.f){
      q <- try(uniroot(grad_gamma1,interval=c(1e-4,1000),ppi=new.pi,j=j,g=new.gam,tol=1e-8), silent=T)
      if("try-error" %in% class(q)){ new.gam[j,1] <- gam[j,1]
      }else{
        new.gam[j,1] <- q$root
      }
    }
    
    grad_gamma2 <- function(x,ppi=NULL,j,g=NULL){
      
      hh1 <- (a2-x)*(trigamma(x)-trigamma(g[j,1]+x))
      hh2 <- (I((1-ppi[,j])>0.5)*(trigamma(x)))*I(cvsample==0)[,j]
      hh3 <- (a1-g[j,1])*trigamma(g[j,1]+x)
      hh4 <- n.s*trigamma(g[j,1]+x)
      hh <- sum(hh2) + hh1 - hh3 - hh4
      
      return(hh)
    }
    
    for(j in 1:n.f){
      q <- try(uniroot(grad_gamma2,interval=c(1e-4,1000),ppi=new.pi,j=j,g=new.gam,tol=1e-8), silent=T)
      if("try-error" %in% class(q)){ new.gam[j,2] <- gam[j,2]
      }else{
        new.gam[j,2] <- q$root
      }
    }
    
    new.eta <- new.gam[,1]/rowSums(new.gam)
    
    #####----------------- Update m -------------------#####
    ## Objective function
    func_m <- function(x,b=NULL,f=NULL,s=NULL,l=NULL,ppi=NULL) {
      new.factor_coefs_0 <- x[1:n.f]
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      lf <- new.factor_scores %*% t(new.factor_coefs_j)
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      y1 <- X*lf*I(cvsample==0)
      y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll))))*I(cvsample==0)
      fun2 <- function(i) { 0.5 * (- sum(new.factor_scores[i,]^2))}
      y <- sum(y1) + sum(y2) + sum(sapply(1:n.s,fun2))
      
      return(y)
    }
    ## Score function
    grad_m <- function(x,b=NULL,f=NULL,s=NULL,l=NULL,ppi=NULL) {
      new.factor_coefs_0 <- x[1:n.f]
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      uu <- U_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      sum <- I((1-new.pi)>0.5)*(exp(ll)*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*(exp(ll))*I(cvsample==0)))
      fun2 <- function(i) { -new.factor_scores[i,]+((X*I(cvsample==0))[i,])%*%new.factor_coefs_j-(M[i]*sum[i,])%*%uu[i,,]}
      f_grad <- t(sapply(1:n.s,fun2))
      return(c(f_grad))
    }
    ## Optimization 
    q <- try(optim(c(factor_scores), x=new.factor_coefs_0,b=new.factor_coefs_j,s=new.sigma,ppi=new.pi,l=new.lambda, method="BFGS", fn=func_m, gr=grad_m, control=list(trace=0, fnscale=-1, maxit=maxit)), silent=TRUE)
    if("try-error" %in% class(q)){ new.factor_scores <- factor_scores;
    }else{
      if(iter > 1 && f.cur.logfunc > q$value){ if(trace)
        cat("Optimization of m did not improve on iteration step ",iter,"\n");
        new.factor_scores <- factor_scores;
      }else{
        if(trace) cat("Variational parameters m updated","\n")
        new.factor_scores <- matrix(q$par,n.s,n.factors);
        if(q$convergence != 0) { if(trace) cat("Optimization of m did not converge on iteration step ", iter,"\n") }
      }
    }
    
    ####----------------- Update sigma -----------------####
    ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
    grad_sigma <- function(x,la=NULL,r=NULL,i,l,ppi=NULL,m=NULL){
      out1 <- function(j) (la[j,l]*(1-la[j,l]*x)+(r[j,l]+m[i,l]*la[j,l])^2)/(1-la[j,l]*x)^2
      out2 <- (I((1-ppi[i,])>0.5)*((exp(ll))*I(cvsample==0))[i,])
      out3 <- rowSums(I((1-ppi)>0.5)*((exp(ll))*I(cvsample==0)))
      out4 <- sapply(1:n.f,out1)
      out <- (1/x) - 1 - (M[i]*sum(out2*out4))/out3
      return(out)
    }
    
    for(i in 1:n.s){
      for(l in 1:n.factors) {
        q <- try(uniroot(grad_sigma,interval=c(1e-4,1),la=new.lambda,r=new.factor_coefs_j,i=i,l=l,ppi=new.pi,m=new.factor_scores,tol=1e-8), silent=T)
        if("try-error" %in% class(q)){ new.sigma[i,l] <- sigma[i,l]
        }else{
          new.sigma[i,l] <- q$root
        }
      }
    }
    
    ####---------------- Update r ----------------####
    ## Objective function
    func_r <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL) {
      new.factor_coefs_0 <- x[1:n.f];
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      lf <- new.factor_scores %*% t(new.factor_coefs_j)
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      
      y1 <- X*lf*I(cvsample==0)
      y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll))))*I(cvsample==0)
      fun3 <- function(j) { 0.5 * (- sum(diag(sigma_beta[j,])%*%new.factor_coefs_j[j,]^2))}
      y <- sum(y1) + sum(y2) + sum(sapply(1:n.f,fun3))
      
      return(y)
    }
    ## Score function
    grad_r <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL) {
      new.factor_coefs_0 <- x[1:n.f];
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      tt <- T_func(s=new.sigma,l=new.lambda,m=new.factor_scores,r=new.factor_coefs_j)
      sum <- I((1-new.pi)>0.5)*(exp(ll)*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*(exp(ll))*I(cvsample==0)))
      b3 <- NULL
      for(p in 1:n.factors){
        b31 <- sweep((X*I(cvsample==0)),1,new.factor_scores[,p],"*")
        b3 <- c(b3,(b31 - tt[,,p]*sum))
      }
      b3 <- matrix(b3,n.s,n.f*n.factors)
      b4 <- colSums(b3)
      b5=NULL
      for(j in 1:n.f) {
        b5 <- c(b5,(1/sigma_beta[j,])*new.factor_coefs_j[j,])
      }
      return(b4-b5)
    }
    ## Optimization
    q <- try(optim(c(factor_coefs_j), x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,ppi=new.pi,l=new.lambda,method="BFGS",fn=func_r, gr=grad_r, control=list(trace=0, fnscale=-1, maxit=maxit)), silent=TRUE)
    if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j;
    }else{
      if(iter > 1 && b.cur.logfunc > q$value){if(trace)
        cat("Optimization of r did not improve on iteration step ",iter,"\n");
        new.factor_coefs_j <- factor_coefs_j;
      }else{
        if(trace) cat("Variational parameters r updated","\n")
        new.factor_coefs_j <- matrix(q$par,n.f,n.factors);
        if(q$convergence != 0) { if(trace) cat("Optimization of r did not converge on iteration step ", iter,"\n") }
      }
    }
    
    ####---------------- Update lambda ------------####
    ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
    grad_lambda <- function(x,s=NULL,m=NULL,j,l,ppi=NULL,r=NULL){
      out1 <- function(i) { (s[i,l]*(1-s[i,l]*x)+(m[i,l]+r[j,l]*s[i,l])^2)/(1-s[i,l]*x)^2 }
      out2 <- (M*I((1-ppi[,j])>0.5)*((exp(ll))*I(cvsample==0))[,j])
      out3 <- rowSums(I((1-ppi)>0.5)*((exp(ll))*I(cvsample==0)))
      out4 <- sapply(1:n.s,out1)
      out <- (1/x) - (1/sigma_beta[j,l]) - sum(out2*out4/out3)
      return(out)
    }
    
    for(j in 1:n.f){
      for(l in 1:n.factors) {
        q <- try(uniroot(grad_lambda,interval=c(1e-4,1),s=new.sigma,m=new.factor_scores,j=j,l=l,ppi=new.pi,r=new.factor_coefs_j,tol=1e-8), silent=T)
        if("try-error" %in% class(q)){ new.lambda[j,l] <- lambda[j,l]
        }else{
          new.lambda[j,l] <- q$root
        }
      }
    }
    
    ####---------------- Update beta0 -------------####
    ## Objective function
    func_b0 <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL) {
      new.factor_coefs_0 <- x[1:n.f];
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      lab <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      y1 <- X*lab*I(cvsample==0)
      y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll))))*I(cvsample==0)
      y <- sum(y1) + sum(y2)
      return(y)
    }
    ## Score function
    grad_b0 <- function(x,b=NULL,f=NULL,s=NULL,ppi=NULL,l=NULL) {
      new.factor_coefs_0 <- x[1:n.f];
      new.factor_coefs_j <- matrix(b,n.f,n.factors)
      new.factor_scores <- matrix(f,n.s,n.factors)
      new.sigma <- matrix(s,n.s,n.factors)
      new.lambda <- matrix(l,n.f,n.factors)
      new.pi <- matrix(ppi,n.s,n.f)
      
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=new.sigma, l=new.lambda, m=new.factor_scores, r=new.factor_coefs_j)
      grad <- X*I(cvsample==0)-M*I((1-new.pi)>0.5)*((exp(ll))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*((exp(ll))*I(cvsample==0))))
      return(c(colSums(grad)))
    }
    ## Optimization
    q <- try(optim(c(factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, ppi=new.pi,l=new.lambda,method="BFGS", fn=func_b0, gr=grad_b0, control=list(trace=0, fnscale=-1, maxit=maxit)), silent=TRUE)
    if("try-error" %in% class(q)){ new.factor_coefs_0 <- factor_coefs_0
    }else{
      if(iter > 1 && b0.cur.logfunc > q$value){if(trace)
        cat("Optimization of beta0 did not improve on iteration step ",iter,"\n");
        new.factor_coefs_0 <- factor_coefs_0
      }else{
        if(trace) cat("Model parameters beta0 updated","\n")
        new.factor_coefs_0 <- q$par
        if(q$convergence != 0) { if(trace) cat("Optimization of beta0 did not converge on iteration step ", iter,"\n") }
      }
    }
    
    q1 <- list(value=func_r(c(new.factor_coefs_j), x=new.factor_coefs_0, f=new.factor_scores, s=new.sigma, ppi=new.pi, l=new.lambda))
    b.new.logfunc <- q1$value
    b.cur.logfunc <- b.new.logfunc
    
    q2 <- list(value=func_m(c(new.factor_scores), x=new.factor_coefs_0, b=new.factor_coefs_j, s=new.sigma, ppi=new.pi, l=new.lambda))
    new.f.cur.logfunc <- q2$value
    f.cur.logfunc <- new.f.cur.logfunc
    
    q3 <- list(value=func_b0(c(new.factor_coefs_0), f=new.factor_scores, b=new.factor_coefs_j, s=new.sigma, ppi=new.pi, l=new.lambda))
    new.b0.cur.logfunc <- q3$value
    b0.cur.logfunc <- new.b0.cur.logfunc
    
    
    ## Take values of VLB to define stopping rule
    q <- list(value=VLB(c(new.factor_coefs_0), b=new.factor_coefs_j, f=new.factor_scores, s=new.sigma, ppi=new.pi, l=new.lambda, g=new.gam))
    new.VLB <- q$value
    diff <- abs(new.VLB-cur.VLB)
    ratio <- abs(new.VLB/cur.VLB);
    if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, "Difference in VLB:", diff,"\n")
    cur.VLB <- new.VLB
    
    if(!is.null(cvsample)){
      q <- list(value=VLB_rept(c(new.factor_coefs_0), b=new.factor_coefs_j, f=new.factor_scores, s=new.sigma, ppi=new.pi, l=new.lambda, g=new.gam))
      cur.VLB_rept <- q$value
    }
    
    gam <- new.gam
    #eta <- gam[,1]/rowSums(gam)
    eta <- new.eta
    factor_coefs_0 <- new.factor_coefs_0
    factor_coefs_j <- new.factor_coefs_j
    ppi <- new.pi
    factor_scores <- new.factor_scores
    sigma <- new.sigma
    lambda <- new.lambda
    iter <- iter + 1
  }
  
  lf <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE) + factor_scores %*% t(factor_coefs_j)
  ll <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE) + L_func(s=sigma, l=lambda, m=factor_scores, r=factor_coefs_j)
  exp.mat.bis <- exp(lf)*I(cvsample==0)
  exp.mat <- exp(ll)*I(cvsample==0)
  sum <- exp.mat/(rowSums(exp.mat)) ## estimates of rho
  sum.bis <- exp.mat.bis/(rowSums(exp.mat.bis))
  sum2 <- exp(ll)*I(cvsample==0)/rowSums(exp(ll)*I(cvsample==0))
  sum2.bis <- exp(lf)*I(cvsample==0)/rowSums(exp(lf)*I(cvsample==0))
  mu_z <- (1-new.pi)*exp.mat
  mu_z.bis <- (1-new.pi)*exp.mat.bis
  mu <- exp.mat
  mu.bis <- exp.mat.bis
  
  if(iter==max.iter){
    warning("ZILPNMVA not converging!")
  }
  
  ## Print the output
  out.list$VLB <- cur.VLB
  if(!is.null(cvsample)){out.list$VLB_rept <- cur.VLB_rept}
  out.list$iter <- iter-1
  out.list$lvs$ppi <- ppi
  out.list$lvs$gam <- gam
  out.list$lvs$factor_scores <- factor_scores
  out.list$lvs$sigma <- sigma
  out.list$lvs$lambda <- lambda
  out.list$lvs$factor_coefs_j <- factor_coefs_j
  out.list$params$factor_coefs_0 <- factor_coefs_0
  out.list$params$eta <- eta
  out.list$Q <- sum
  out.list$Q.bis <- sum.bis
  out.list$Q2 <- sum2
  out.list$Q2.bis <- sum2.bis
  out.list$mu <- mu
  out.list$mu.bis <- mu.bis
  out.list$muz <- mu_z
  out.list$muz.bis <- mu_z.bis
  
  
  ## Compute H
  tt <- T_func(s=new.sigma,l=new.lambda,m=new.factor_scores,r=new.factor_coefs_j)
  uu <- U_func(s=new.sigma,l=new.lambda,m=new.factor_scores,r=new.factor_coefs_j)
  rr <- R_func(s=new.sigma,l=new.lambda,m=new.factor_scores,r=new.factor_coefs_j)
  pp <- P_func(s=new.sigma,l=new.lambda,m=new.factor_scores,r=new.factor_coefs_j)
  alpha <- log(M/(rowSums((1-new.pi)*exp(ll))))
  #alpha <- log(M/(rowSums(exp(ll))))
  lll <- matrix(alpha,n.s,n.f) + ll
  wij <- exp(lll)
  
  ss1 <- ss2 <- ss3 <- matrix(0,n.s,n.f,byrow = T)
  for(i in 1:n.s){
    ss1[i,] <- exp(digamma(new.gam[,2]))*exp(-wij[i,])*(-wij[i,])/(exp(digamma(new.gam[,1]))+exp(digamma(new.gam[,2]))*exp(-wij[i,]))^2
    ss2[i,] <- (1-wij[i,])*exp(digamma(new.gam[,1]))+exp(digamma(new.gam[,2]))*exp(-wij[i,])
    ss3[i,] <- exp(digamma(new.gam[,1])) + exp(digamma(new.gam[,2]))*exp(-wij[i,])
  }
  #browser()
  d <- ss1*ss2
  v <- ss1*ss3
  
  ## 1. H1
  wt <- ifelse(X>0,-wij,d)
  H11 <- diag(rowSums(wt))
  #browser()
  t1 <- t2 <- array(NA, dim=c(n.s,n.f,n.factors))
  for(i in 1:n.s) {
    t1[i,,] <- d[i,]*uu[i,,]
    t2[i,,] <- wij[i,]*uu[i,,]
  }
  
  wt12 <- NULL 
  for(l in 1:n.factors) {
    wt12 <- rbind(wt12,rowSums(ifelse(X>0,-t2[,,l],t1[,,l])))
  }
  
  H12 <- matrix(wt12[,1],ncol=n.factors)
  for(i in 2:n.s) {
    H12 <- Matrix::bdiag(H12,matrix(wt12[,i],ncol=n.factors))
  }
  
  t3 <- t4 <- array(NA, dim=c(n.s,n.f,n.factors))
  for(i in 1:n.s) {
    t3[i,,] <- 0.5*d[i,]*rr[i,,]
    t4[i,,] <- 0.5*wij[i,]*rr[i,,]
  }
  
  wt34 <- NULL 
  for(l in 1:n.factors) {
    wt34 <- rbind(wt34,rowSums(ifelse(X>0,-t4[,,l],t3[,,l])))
  }
  
  H13 <- matrix(wt34[,1],ncol=n.factors)
  for(i in 2:n.s) {
    H13 <- Matrix::bdiag(H13,matrix(wt34[,i],ncol=n.factors))
  }
  
  func_ar <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out1 <- (d[i,j]*tt[i,j,]*I(X[i,j]==0))+(-wij[i,j]*tt[i,j,]*I(X[i,j]>0))
      out <- c(out,out1)
    }
    return(out)
  }
  H14 <- func_ar(1)
  for(i in 2:n.s) {
    H14 <- rbind(H14,func_ar(i))
  }
  
  func_al <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out1 <- (0.5*d[i,j]*pp[i,j,]*I(X[i,j]==0))+(-0.5*wij[i,j]*pp[i,j,]*I(X[i,j]>0))
      out <- c(out,out1)
    }
    return(out)
  }
  
  H15 <- func_al(1)
  for(i in 2:n.s) {
    H15 <- rbind(H15,func_al(i))
  }
  
  func_agam1 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out0 <- -exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,1],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2
      out1 <- (out0*I(X[i,j]==0))+(0*I(X[i,j]>0))
      out <- c(out,out1)
    }
    return(out)
  }
  H16 <- func_agam1(1)
  for(i in 2:n.s) {
    H16 <- rbind(H16,func_agam1(i))
  }
  
  func_agam2 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out0 <- exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,2],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2
      out1 <- (out0*I(X[i,j]==0))+(0*I(X[i,j]>0))
      out <- c(out,out1)
    }
    return(out)
  } 
  H17 <- func_agam2(1)
  for(i in 2:n.s) {
    H17 <- rbind(H17,func_agam2(i))
  }
  
  H1 <- cbind(H11,H12,H13,H14,H15,H16,H17)
  #browser()
  ## 2. H2 with Poisson model
  #H21 <- MatrixExtra::t_shallow(H12)
  H21 <- t(as.matrix(H12))
  
  func_mm <- function(i) {
    out <- 0
    for(j in 1:n.f) {
      out1 <- (-wij[i,j]*((uu[i,j,]%*%t(uu[i,j,]))+diag(new.lambda[j,]/(1-new.sigma[i,]*new.lambda[j,]+1e-8))))
      out2 <- (d[i,j]*(uu[i,j,]%*%t(uu[i,j,]))+v[i,j]*diag(new.lambda[j,]/(1-new.sigma[i,]*new.lambda[j,]+1e-8)))
      out <- out + out1*I(X[i,j]>0) + out2*I(X[i,j]==0)
    }
    out <- out - diag(n.factors)
    return(out)
  }
  H22 <- func_mm(1)
  for(i in 2:n.s) {
    H22 <- Matrix::bdiag(H22,func_mm(i))
  }
  
  func_ms <- function(i) {
    out <- 0
    for(j in 1:n.f) {
      out1 <- (-0.5*wij[i,j]*((uu[i,j,]%*%t(rr[i,j,]))+diag(2*new.lambda[j,]*(new.factor_coefs_j[j,]+new.factor_scores[i,]*new.lambda[j,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2)))
      out2 <- (0.5*d[i,j]*(uu[i,j,]%*%t(rr[i,j,]))+0.5*v[i,j]*diag(2*new.lambda[j,]*(new.factor_coefs_j[j,]+new.factor_scores[i,]*new.lambda[j,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2))
      out <- out + out1*I(X[i,j]>0) + out2*I(X[i,j]==0)
    }
    return(out)
  }
  H23 <- func_ms(1)
  for(i in 2:n.s) {
    H23 <- Matrix::bdiag(H23,func_ms(i))
  }
  
  func_mr <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      
      q <- 1/(1-new.sigma[i,]*new.lambda[j,]+1e-8) 
      out1 <- X[i,j]-wij[i,j]*((uu[i,j,]%*%t(tt[i,j,]))+diag(q))
      out2 <- d[i,j]*(uu[i,j,]%*%t(tt[i,j,]))+v[i,j]*diag(q)
      out <- cbind(out,out1*I(X[i,j]>0)+out2*I(X[i,j]==0))
    }
    return(out)
  }
  H24 <- func_mr(1)
  for(i in 2:n.s) {
    H24 <- rbind(H24,func_mr(i))
  }
  
  func_ml <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out1 <- -0.5*wij[i,j]*((uu[i,j,]%*%t(pp[i,j,]))+diag(2*(new.factor_scores[i,]+new.factor_coefs_j[j,]*new.sigma[i,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2))
      out2 <- 0.5*d[i,j]*(uu[i,j,]%*%t(pp[i,j,]))+0.5*v[i,j]*diag(2*(new.factor_scores[i,]+new.factor_coefs_j[j,]*new.sigma[i,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2)
      out <- cbind(out,out1*I(X[i,j]>0)+out2*I(X[i,j]==0))
    }
    return(out)
  }
  
  H25 <- func_ml(1)
  for(i in 2:n.s) {
    H25 <- rbind(H25,func_ml(i))
  }
  
  func_mgam1 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out1 <- -exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,1],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*uu[i,j,] 
      out <- cbind(out,rep(0,n.factors)*I(X[i,j]>0)+out1*I(X[i,j]==0))
    }
    return(out)
  }
  
  H26 <- func_mgam1(1)
  for(i in 2:n.s) {
    H26 <- rbind(H26,func_mgam1(i))
  }
  
  func_mgam2 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      out1 <- exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,2],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*uu[i,j,]
      out <- cbind(out,rep(0,n.factors)*I(X[i,j]>0)+out1*I(X[i,j]==0))
    }
    return(out)
  } 

  H27 <- func_mgam2(1)
  for(i in 2:n.s) {
    H27 <- rbind(H27,func_mgam2(i))
  }
  
  H2 <- cbind(H21,H22,H23,H24,H25,H26,H27)
  #browser()
  ## 3. H3 with Poisson model
  #H31 <- MatrixExtra::t_shallow(H13)
  #H32 <- MatrixExtra::t_shallow(H23)
  H31 <- t(as.matrix(H13))
  H32 <- t(as.matrix(H23))
  
  func_ss <- function(i) {
    out <- 0
    for(j in 1:n.f) {
      q <- (new.lambda[j,]^2/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2)+2*new.lambda[j,]*(new.factor_coefs_j[j,]+new.factor_scores[i,]*new.lambda[j,])^2/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^3 
      q1 <- (-wij[i,j]*(0.25*rr[i,j,]%*%t(rr[i,j,])+0.5*diag(q)))
      q2 <- (0.25*d[i,j]*(rr[i,j,]%*%t(rr[i,j,]))+0.5*v[i,j]*diag(q))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    out <- out - 0.5*(1/new.sigma[i,])%*%t(1/new.sigma[i,])
    return(out)
  }
  H33 <- func_ss(1)
  for(i in 2:n.s) {
    H33 <- Matrix::bdiag(H33,func_ss(i))
  }
  
  func_sr <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      q <- (new.factor_coefs_j[j,]+new.factor_scores[i,]*new.lambda[j,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2
      q1 <- -wij[i,j]*((0.5*rr[i,j,]%*%t(tt[i,j,]))+diag(q))
      q2 <- 0.5*d[i,j]*(rr[i,j,]%*%t(tt[i,j,]))+v[i,j]*diag(q)
      out <- cbind(out,q1*I(X[i,j]>0)+q2*I(X[i,j]==0))
    }
    return(out)
  }
  H34 <- func_sr(1)
  for(i in 2:n.s) {
    H34 <- rbind(H34,func_sr(i))
  }
  
  func_sl <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      q <- (1/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2)+2*(new.factor_coefs_j[j,]+new.factor_scores[i,]*new.lambda[j,])*(new.factor_scores[i,]+new.factor_coefs_j[j,]*new.sigma[i,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^3
      q1 <- -0.5*wij[i,j]*((0.5*rr[i,j,]%*%t(pp[i,j,]))+diag(q)) 
      q2 <- 0.25*d[i,j]*(rr[i,j,]%*%t(pp[i,j,]))+0.5*v[i,j]*diag(q)
      out <- cbind(out,q1*I(X[i,j]>0)+q2*I(X[i,j]==0))
    }
    return(out)
  }
  H35 <- func_sl(1)
  for(i in 2:n.s) {
    H35 <- rbind(H35,func_sl(i))
  }
  
  func_sgam1 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      q1 <- rep(0,n.factors)
      q2 <- -0.5*exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,1],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*rr[i,j,]
      out <- cbind(out,q1*I(X[i,j]>0)+q2*I(X[i,j]==0))
    }
    return(out)
  } 
  H36 <- func_sgam1(1)
  for(i in 2:n.s) {
    H36 <- rbind(H36,func_sgam1(i))
  }
  
  func_sgam2 <- function(i) {
    out <- NULL
    for(j in 1:n.f) {
      q1 <- rep(0,n.factors)
      q2 <- 0.5*exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,2],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*rr[i,j,]
      out <- cbind(out,q1*I(X[i,j]>0)+q2*I(X[i,j]==0))
    }
    return(out)
  } 
  H37 <- func_sgam2(1)
  for(i in 2:n.s) {
    H37 <- rbind(H37,func_sgam2(i))
  }
  
  H3 <- cbind(H31,H32,H33,H34,H35,H36,H37)

  ## 4. H4 with Poisson model
  #H41 <- MatrixExtra::t_shallow(H14)
  #H42 <- MatrixExtra::t_shallow(H24)
  #H43 <- MatrixExtra::t_shallow(H34)
  H41 <- t(as.matrix(H14))
  H42 <- t(as.matrix(H24))
  H43 <- t(as.matrix(H34))
  
  func_rr <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q <- new.sigma[i,]/(1-new.sigma[i,]*new.lambda[j,]+1e-8)
      q1 <- (-wij[i,j]*((tt[i,j,]%*%t(tt[i,j,]))+diag(q)))
      q2 <- (d[i,j]*(tt[i,j,]%*%t(tt[i,j,]))+v[i,j]*diag(q))
      out <- out + q1*I(X[i,j]>0) +q2*I(X[i,j]==0)
    }
    out <- out - sigma_beta[j,]
    return(out)
  }
  H44 <- func_rr(1)
  for(j in 2:n.f) {
    H44 <- Matrix::bdiag(H44,func_rr(j))
  }
  
  func_rl <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q <- 2*new.sigma[i,]*(new.factor_scores[i,]+new.factor_coefs_j[j,]*new.sigma[i,])/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2
      q1 <- (-0.5*wij[i,j]*((tt[i,j,]%*%t(pp[i,j,]))+diag(q)))
      q2 <- 0.5*(d[i,j]*(tt[i,j,]%*%t(pp[i,j,]))+v[i,j]*diag(q))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    return(out)
  }
  H45 <- func_rl(1)
  for(j in 2:n.f) {
    H45 <- Matrix::bdiag(H45,func_rl(j))
  }
  
  func_rgam1 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- rep(0,n.factors)
      q2 <- (-exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,1],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*tt[i,j,])
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    return(out)
  }
  H46 <- func_rgam1(1)
  for(j in 2:n.f) {
    H46 <- Matrix::bdiag(H46,func_rgam1(j))
  }
  
  func_rgam2 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- rep(0,n.factors)
      q2 <- (exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,2],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*tt[i,j,])
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    return(out)
  }
  H47 <- func_rgam2(1)
  for(j in 2:n.f) {
    H47 <- Matrix::bdiag(H47,func_rgam2(j))
  }
  
  H4 <- cbind(H41,H42,H43,H44,H45,H46,H47)
  
  
  ## 5. H5 with Poisson model
  #H51 <- MatrixExtra::t_shallow(H15)
  #H52 <- MatrixExtra::t_shallow(H25)
  #H53 <- MatrixExtra::t_shallow(H35)
  #H54 <- MatrixExtra::t_shallow(H45)
  H51 <- t(as.matrix(H15))
  H52 <- t(as.matrix(H25))
  H53 <- t(as.matrix(H35))
  H54 <- t(as.matrix(H45))
  
  func_ll <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q <- (new.sigma[i,]^2/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^2)+2*new.sigma[i,]*(new.factor_scores[i,]+new.factor_coefs_j[j,]*new.sigma[i,])^2/(1-new.sigma[i,]*new.lambda[j,]+1e-8)^3 
      q1 <- (-wij[i,j]*(0.25*pp[i,j,]%*%t(pp[i,j,])+0.5*diag(q)))
      q2 <- (0.25*d[i,j]*(pp[i,j,]%*%t(pp[i,j,]))+0.5*v[i,j]*diag(q))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    out <- out - 0.5*(1/new.lambda[j,])%*%t(1/new.lambda[j,])
    return(out)
  }
  
  H55 <- func_rr(1)
  for(j in 2:n.f) {
    H55 <- Matrix::bdiag(H55,func_rr(j))
  }
  
  func_lgam1 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- rep(0,n.factors)
      q2 <- (-0.5*exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,1],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*pp[i,j,])
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    return(out)
  }
  H56 <- func_lgam1(1)
  for(j in 2:n.f) {
    H56 <- Matrix::bdiag(H56,func_lgam1(j))
  }
  
  func_lgam2 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- rep(0,n.factors)
      q2 <- (0.5*exp(digamma(new.gam[j,1])+digamma(new.gam[j,2]))*psigamma(new.gam[j,2],deriv=1)*exp(-wij[i,j])*(-wij[i,j])/(exp(digamma(new.gam[j,1]))+exp(digamma(new.gam[j,2]))*exp(-wij[i,j]))^2*pp[i,j,])
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    return(out)
  }
  H57 <- func_lgam2(1)
  for(j in 2:n.f) {
    H57 <- Matrix::bdiag(H57,func_lgam2(j))
  }
  
  H5 <- cbind(H51,H52,H53,H54,H55,H56,H57)
  
  ss4 <- ss5 <- ss6 <- matrix(0,n.s,n.f,byrow = T)
  for(i in 1:n.s){
    q1 <- (exp(digamma(new.gam[,1]))+exp(digamma(new.gam[,2]))*exp(-wij[i,]))^2
    ss4[i,] <- exp(digamma(new.gam[,1]))*(psigamma(new.gam[,1],deriv=2)*exp(digamma(new.gam[,1]))+(psigamma(new.gam[,1],deriv=2)+(trigamma(new.gam[,1]))^2)*exp(digamma(new.gam[,2]))*exp(-wij[i,]))/q1
    ss5[i,] <- (exp(digamma(new.gam[,2]))*exp(-wij[i,]))*(psigamma(new.gam[,2],deriv=2)*exp(digamma(new.gam[,2]))*exp(-wij[i,])+(psigamma(new.gam[,2],deriv=2)+(trigamma(new.gam[,2]))^2)*exp(digamma(new.gam[,1])))/q1 
    ss6[i,] <- exp(digamma(new.gam[,1])+digamma(new.gam[,2]))*psigamma(new.gam[,1],deriv=1)*psigamma(new.gam[,2],deriv=1)*exp(-wij[i,])/q1
  }
  
  ## 6. H6 with Poisson model
  #H61 <- MatrixExtra::t_shallow(H16)
  #H62 <- MatrixExtra::t_shallow(H26)
  #H63 <- MatrixExtra::t_shallow(H36)
  #H64 <- MatrixExtra::t_shallow(H46)
  #H65 <- MatrixExtra::t_shallow(H56)
  H61 <- t(as.matrix(H16))
  H62 <- t(as.matrix(H26))
  H63 <- t(as.matrix(H36))
  H64 <- t(as.matrix(H46))
  H65 <- t(as.matrix(H56))
  
  func_g1g1 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- (-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
      q2 <- (ss4[i,j]-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    out1 <- 2*(trigamma(new.gam[j,1])-trigamma(new.gam[j,1]+new.gam[j,2]))
    out2 <- (a2-new.gam[j,2])*(psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
    out3 <- (a1-new.gam[j,1])*(psigamma(new.gam[j,1],deriv=2)-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
    out <- out - out1 - out2 + out3
    return(out)
  }
  
  h66_el <- sapply(1:n.f,func_g1g1)
  H66 <- diag(h66_el)
  
  func_g1g2 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- (-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
      q2 <- (-(ss6[i,j]+psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2)))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    out1 <- trigamma(new.gam[j,1]+new.gam[j,2]) - (a1-new.gam[j,1])*psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2)
    out2 <- psigamma((new.gam[j,1]+new.gam[j,2]),deriv=1)-(a2-new.gam[j,2])*psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2)
    out <- out + out1 + out2
    return(out)
  }
  h67_el <- sapply(1:n.f,func_g1g2)
  H67 <- diag(h67_el)
  
  H6 <- cbind(H61,H62,H63,H64,H65,H66,H67)
  
  ## 7. H7 with Poisson model
  #H71 <- MatrixExtra::t_shallow(H17)
  #H72 <- MatrixExtra::t_shallow(H27)
  #H73 <- MatrixExtra::t_shallow(H37)
  #H74 <- MatrixExtra::t_shallow(H47)
  #H75 <- MatrixExtra::t_shallow(H57)
  #H76 <- MatrixExtra::t_shallow(H67)
  H71 <- t(as.matrix(H17))
  H72 <- t(as.matrix(H27))
  H73 <- t(as.matrix(H37))
  H74 <- t(as.matrix(H47))
  H75 <- t(as.matrix(H57))
  H76 <- t(as.matrix(H67))
  
  func_g2g2 <- function(j) {
    out <- 0
    for(i in 1:n.s) {
      q1 <- ((psigamma(new.gam[j,2],deriv=2)-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2)))
      q2 <- (ss5[i,j]-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
      out <- out + q1*I(X[i,j]>0) + q2*I(X[i,j]==0)
    }
    out1 <- 2*(trigamma(new.gam[j,2])-trigamma(new.gam[j,1]+new.gam[j,2]))
    out2 <- (a1-new.gam[j,1])*(psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
    out3 <- (a2-new.gam[j,2])*(psigamma(new.gam[j,2],deriv=2)-psigamma((new.gam[j,1]+new.gam[j,2]),deriv=2))
    out <- out - out1 - out2 + out3
    return(out)
  }
  h77_el <- sapply(1:n.f,func_g2g2)
  H77 <- diag(h77_el)
  
  H7 <- cbind(H71,H72,H73,H74,H75,H76,H77)
  
  H <- -rbind(H1,H2,H3,H4,H5,H6,H7)
  HH <- as.matrix(H)
  out.list$H <- HH
  #Hinv <- MASS::ginv(HH)
  Hinv <- chol2inv(HH)
  
  ## Compute g_delta
  g1 <- matrix(0,nrow=n.f*(n.factors+1),ncol=n.s) ## g for alpha
  g2 <- matrix(0,nrow=n.f*(n.factors+1),ncol=(n.s*n.factors)) ## g for m
  g3 <- matrix(0,nrow=n.f*(n.factors+1),ncol=(n.s*n.factors)) ## g for sigma
  g4 <- rbind(diag(n.f*n.factors),matrix(0,nrow=n.f,ncol=(n.f*n.factors))) ## g for r
  g5 <- matrix(0,nrow=n.f*(n.factors+1),ncol=(n.f*n.factors))
  g6 <- rbind(matrix(0,nrow=(n.f*n.factors),ncol=n.f),diag(new.gam[,2]/(rowSums(new.gam)^2)))
  g7 <- rbind(matrix(0,nrow=(n.f*n.factors),ncol=n.f),diag(-new.gam[,1]/(rowSums(new.gam)^2)))
  g_delta <- cbind(g1,g2,g3,g4,g5,g6,g7)
  out.list$g_delta <- g_delta
  
  mat_lrvb <- g_delta%*%Hinv%*%t(g_delta)
  #browser()
  ##########
  alfa <- (1 - level) / 2
  
  blow <- new.factor_coefs_j + qnorm(alfa) * (sqrt(diag(mat_lrvb))[1:(n.f*n.factors)])
  bup <- new.factor_coefs_j + qnorm(1 - alfa) * (sqrt(diag(mat_lrvb))[1:(n.f*n.factors)])
  
  elow <- new.eta + qnorm(alfa) * (sqrt(diag(mat_lrvb))[(1+n.f*n.factors):(n.f+n.f*n.factors)])
  eup <- new.eta + qnorm(1 - alfa) * (sqrt(diag(mat_lrvb))[(1+n.f*n.factors):(n.f+n.f*n.factors)])
  
  
  conf <- cbind(c(blow,elow), c(bup,eup))
  
  colnames(conf) <-  c(paste(alfa * 100, "%"), paste((1 - alfa) * 100, "%"))
  rownames(conf) <- c(rep("beta",n.f*n.factors), rep("eta",n.f))
  
  out.list$sd <- sqrt(diag(mat_lrvb))
  out.list$conf <- conf
  
  return(out.list)
  
}


############################################################################
#####------------------------ Simulate the model ----------------------#####
############################################################################
proj_simulate_model <- function(n=50, p=100, k=2, sh1=2, sh2=2, la=2, b0) {
  
  lam <- la*diag(k); rr <- rep(0,k)
  betaj <- matrix(0, nrow=p, ncol=k)
  for(j in 1:p){
    betaj[j,] <- mvrnorm(1, rr, lam)
  } 
  
  etaj <- rbeta(n=p, shape1=sh1, shape2=sh2)
  
  si <- diag(k); me <- rep(0, k)
  f <- matrix(0, nrow=n, ncol=k)
  for(i in 1:n){
    f[i,] <- mvrnorm(1, me, si)
  }
  
  l <- matrix(b0, n, p, byrow=TRUE) + f %*%  t(betaj)
  
  z <- matrix(0, n, p, byrow=TRUE)
  for(i in 1:n){
    z[i,] <- rbinom(p, size=1, prob=etaj)
  }
  
  sum <- matrix(rowSums((1-z)*exp(l)),n,p)
  Qn_z <- matrix(I(z==0)*(exp(l)/sum)+I(z==1)*0,n,p)
  
  X <- matrix(0, n, p, byrow=TRUE)
  for(i in 1:n){
    X[i,] <- rmultinom(1, size=runif(1,800,1000), prob=Qn_z[i,])
  }
  
  zerorow <- which(rowSums(X)==0)
  if(length(zerorow)>0){
    X <- X[-zerorow,]; Qn_z <- Qn_z[-zerorow,];
  }
  zerocol <- which(colSums(X)==0)
  if(length(zerocol)>0 ){
    X <- X[,-zerocol]; Qn_z <- Qn_z[,-zerocol];
  }
  return(list(dat=X,eta=etaj,beta=betaj,factors=f,rho=Qn_z))
}


