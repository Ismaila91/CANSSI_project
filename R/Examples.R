n <- 100; p <- 5; k <- 2
beta0 <- rep(2, p)
X <- proj_simulate_model(n=n, p=p, k=k, b0=beta0, la=0.1)
out <- ZILPNMVA_Poisson(X$dat,n.factors=k,trace=TRUE)

## Block LRVB for beta
g_delta_beta <- out$g_delta[1:10,501:510]
H_beta <- out$H[501:510,501:510]
H_beta_inverse <- chol2inv(H_beta)
lrvb_beta <- g_delta_beta%*%H_beta_inverse%*%t(g_delta_beta)
sd_beta <- sqrt(diag(lrvb_beta))

## Block LRVB for gam1
g_delta_gam1 <- out$g_delta[11:15,521:525]
H_gam1 <- out$H[521:525,521:525]
H_gam1_inverse <- chol2inv(H_gam1)
lrvb_gam1 <- g_delta_gam1%*%H_gam1_inverse%*%t(g_delta_gam1)
sd_gam1 <- sqrt(diag(lrvb_gam1))


## Block LRVB for gam2
g_delta_gam2 <- out$g_delta[11:15,526:530]
H_gam2 <- out$H[526:530,526:530]
H_gam2_inverse <- chol2inv(H_gam2)
lrvb_gam2 <- g_delta_gam2%*%H_gam2_inverse%*%t(g_delta_gam2)
sd_gam2 <- sqrt(diag(lrvb_gam2))

## Block LRVB for eta
sd_eta <- sd_gam1/(sd_gam1+sd_gam2)
level <- 0.95
alfa <- (1 - level) / 2
eta_low <- out$params$eta + qnorm(alfa) * sd_eta
eta_up <- out$params$eta + qnorm(1 - alfa) * sd_eta
conf_eta <- cbind(eta_low,eta_up)

out$sd
X$eta


