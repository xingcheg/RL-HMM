####################################################
####  Functions for the EM algorithm of RL-HMM  ####
####################################################

library(genlasso)
library(Matrix)
library(parallel)


############################ likelihood computation #############################
## update coefficient delta
update_delta <- function(a, s, R, delta0, beta, rho, eva_phi){
  phi0 <- eva_phi(a,s)
  dd <- rho * R - sum(delta0 * phi0)
  delta1 <- delta0 + beta * dd * phi0
  return(delta1)
}



## RL learning probability: P(At | Ut=1), (for a single subject)
PA_U1 <- function(A, S, R, delta0, beta, rho, nu, eva_phi){
  N <- length(A)
  probs <- rep(0, N)
  for (t in 1:N){
    a <- A[t]
    s <- S[t]
    r <- R[t]
    phi_a <- eva_phi(a, s)
    phi_b <- eva_phi(1-a, s)
    Q_diff <- sum((phi_a - phi_b) * delta0)
    probs[t] <- 1 / (1 + exp( -Q_diff - sign(a-0.5)*nu ))
    if (t < N){
      delta1 <- update_delta(a, s, r, delta0, beta, rho, eva_phi)
      delta0 <- delta1
    }
  }
  return(probs)
}



## RL learning probability, (for all subjects)
PA_U1_mat <- function(A_mat, S_mat, R_mat, delta0, beta, rho, nu, eva_phi, cores){
  n <- dim(A_mat)[1]
  N <- dim(A_mat)[2]
  PA_list <- mclapply(1:n, FUN = function(i){
    PA_U1(A_mat[i,], S_mat[i,], R_mat[i,], delta0, beta, rho, nu, eva_phi)}, 
    mc.cores = cores)
  PA_mat <- matrix(unlist(PA_list), n, N, byrow = TRUE)
  return(PA_mat)
}



## transition probability: p(Ut+1=k | Ut=j), (for all subjects)
P_trans_nonstat <- function(zeta0, zeta1){
  N1 <- length(zeta0)
  C_arr <- array(0, dim = c(N1,2,2))
  C01 <- 1/(1+exp(-zeta0))
  C11 <- 1/(1+exp(-zeta1))
  C_arr[,1,1] <- 1 - C01
  C_arr[,1,2] <- C01
  C_arr[,2,1] <- 1 - C11
  C_arr[,2,2] <- C11
  return(C_arr)
}



############################# EM algorithm: M-step ########################

################## update (pi_1) in HMM ############3###
## "gamma_n_array" is an n x N x 2 array.
update_pi1_0 <- function(gamma_n_array){
  pi1_new <- mean(gamma_n_array[,1,2])
  return(pi1_new)
}



############ update (zeta_0) and (zeta_1) in HMM ########
## this function returns an (n-1) x n fused lasso matrix
diff1D <- function(n){
  D1 <- cbind(-1 * diag(n-1), rep(0, n-1))
  D2 <- cbind(rep(0, n-1), diag(n-1))
  D <- as(D1+D2, "sparseMatrix")
  return(D)
}



## this function returns an (n-k) x n trend filtering matrix of order k
diffkD <- function(n, k){
  D <- diff1D(n)
  if (k==1){
    return(D)
  } else {
    for (i in 1:(k-1)){
      D1 <- diff1D(n-i)
      D <- D1 %*% D
    }
    return(D)
  }
}



## gradient of (zeta_0) and (zeta_1)
grad_zeta <- function(xi_n_array, C_array){
  ## "xi_n_array" is an n x (N-1) x 2 x 2 array.
  ## "C_array" is an (N-1) x 2 x 2 array.
  ## return length-(N-1) vectors for the gradients.
  A <- apply(xi_n_array, c(2,3,4), sum)
  out <- -A[,,2] + (A[,,1] + A[,,2]) * C_array[,,2]
  g_zeta0 <- out[,1]
  g_zeta1 <- out[,2]
  return(list(g_zeta0 = g_zeta0, g_zeta1 = g_zeta1))
}



## Hessian matrix of (zeta_0) and (zeta_1)
hessian_zeta <- function(xi_n_array, C_array){
  ## return length-(N-1) vectors for the diagonal values of the Hessian matrix.
  A <- apply((xi_n_array[,,,1] + xi_n_array[,,,2]), c(2, 3), sum)
  B <- C_array[,,2] * (1 - C_array[,,2])
  out <- A * B
  H_zeta0 <- out[,1]
  H_zeta1 <- out[,2]
  return(list(H_zeta0 = H_zeta0, H_zeta1 = H_zeta1))
}



## one step update of (zeta_0) and (zeta_1) by genlasso
update_zeta <- function(zeta0, zeta1, xi_n_array, C_array, D, lam_rank){
  ## "D" is the generalized lasso matrix
  ## "lam_rank" is the rank of penalty parameter lambda from big to small (in total N-2 lambda)
  gg <- grad_zeta(xi_n_array, C_array)
  g_zeta0 <- gg$g_zeta0
  g_zeta1 <- gg$g_zeta1
  hh <- hessian_zeta(xi_n_array, C_array)
  H_zeta0 <- hh$H_zeta0
  H_zeta1 <- hh$H_zeta1
  J_zeta0 <- sqrt(H_zeta0)
  J_zeta1 <- sqrt(H_zeta1)
  y0 <- J_zeta0 * zeta0 - (g_zeta0 / J_zeta0)
  y1 <- J_zeta1 * zeta1 - (g_zeta1 / J_zeta1)
  D0 <- t(t(D)/J_zeta0)
  D1 <- t(t(D)/J_zeta1)
  fit0 <- genlasso(y = y0, D = D0)
  zeta0_new <- fit0$beta[,lam_rank] / J_zeta0
  fit1 <- genlasso(y = y1, D = D1)
  zeta1_new <- fit1$beta[,lam_rank] / J_zeta1
  return(list(zeta0_new = zeta0_new, zeta1_new = zeta1_new))
}




############### update RL parameters #########
## objective function of RL model
obj_p1 <- function(pars, A_mat, S_mat, R_mat, gamma1_mat, eva_phi, cores){
  ## gamma1_mat = gamma_n_array[,,2] is an n x N matrix.
  beta <- pars[1]
  rho <- pars[2]
  nu <- pars[3]
  alpha <- pars[4]
  delta0 <- alpha * c(1, -1, 0, 1)
  pp <- PA_U1_mat(A_mat, S_mat, R_mat, delta0, beta, rho, nu, eva_phi, cores)
  out <- sum(-gamma1_mat * log(pp))
  return(out)
}



## m-step update using L-BFGS-B algorithm
update_RL <- function(pars, A_mat, S_mat, R_mat, 
                      gamma1_mat, eva_phi, cores, maxit=5){
  opt <- optim(par = pars, fn = obj_p1,
               lower = c(1e-6, 1e-3, -2, 0), 
               upper = c(0.5, 10, 2, 3), 
               method = "L-BFGS-B", control = list(maxit = maxit),
               A_mat = A_mat, S_mat = S_mat, R_mat = R_mat,
               gamma1_mat = gamma1_mat, eva_phi = eva_phi, cores = cores)
  pars1 <- opt$par
  return(pars1)
}



################################################################
######################## EM algorithm  #########################
################################################################
## S_mat (states): n x N matrix
## A_mat (actions): n x N matrix
## R_mat (rewards): n x N matrix
## beta0 (initial value of learning rate in RL)
## rho0 (initial value of reward sensitivity in RL)
## alpha0 (initial value of RL starting coefs)
## nu0 (initial value of RL bias)
## pi10 (initial value of pi_1 in HMM)
## zeta0_0 (initial value of zeta_0 in HMM): length-(N-1) vector
## zeta1_0 (initial value of zeta_1 in HMM): length-(N-1) vector
## order (trend filtering order in genlasso)
## lam_rank (rank of penalty parameter lambda from 1 to (N-2) from large to small) 
## eva_phi (user specified basis for Q-function)
EM_RL_HMM <- function(S_mat, A_mat, R_mat, beta0, rho0, alpha0, nu0,
                      pi10, zeta0_0, zeta1_0, order=1, lam_rank, eva_phi,
                      M_max = 1000, tol = 1e-3, cores = 8){
  problem_size <- dim(S_mat)
  n <- problem_size[1]
  N <- problem_size[2]
  D <- diffkD(N-1, order)
  pi1_all <- rep(0, M_max)
  zeta_all <- array(0, dim = c(M_max, N-1, 2))
  beta_all <- rep(0, M_max)
  rho_all <- rep(0, M_max)
  alpha_all <- rep(0, M_max)
  nu_all <- rep(0, M_max)
  
  ## EM updates
  for (kk in 1:M_max){
    cat("-----------iter = ", kk, " ------------\n")
    ## compute P(A|U) 
    
    st1 <- Sys.time()
    PA_mat <- array(1/2, dim = c(n, N, 2))
    delta00 <- alpha0 * c(1, -1, 0, 1)
    PA_mat[,,2] <- PA_U1_mat(A_mat, S_mat, R_mat, delta00, beta0, rho0, nu0, eva_phi, cores)
    st2 <- Sys.time()
    st12 <- difftime(st2, st1, units = "mins")
    cat("compute probs using:", st12, "minutes, done! \n")
    
    ## E-step
    gamma_arr <- array(0, dim = c(n, N, 2))
    xi_arr <- array(0, dim = c(n, N-1, 2, 2))
    C0_arr <- P_trans_nonstat(zeta0_0, zeta1_0)
    for (i in 1:n){
      Pa_i <- PA_mat[i,,]
      aa_out <- forward_pass(Pa_i, pi10, C0_arr)
      bb_out <- backward_pass(Pa_i, C0_arr)
      gamma_arr[i,,] <- comp_gamma(aa_out$aa_mat, bb_out$bb_mat,
                                   aa_out$digit_aa, bb_out$digit_bb)
      xi_arr[i,,,] <- comp_xi(aa_out$aa_mat, bb_out$bb_mat,
                              aa_out$digit_aa, bb_out$digit_bb, 
                              Pa_i, C0_arr)
    }
    st3 <- Sys.time()
    st23 <- difftime(st3, st2, units = "mins")
    cat("E-step computing using:", st23, "minutes, done! \n")
    
    ## M-step
    pi11 <- update_pi1_0(gamma_arr)
    pi11 <- max(1e-4, pi11)
    pi11 <- min(1-1e-4, pi11)
    zeta_out <- update_zeta(zeta0_0, zeta1_0, xi_arr, C0_arr, D, lam_rank)
    zeta0_1 <- zeta_out$zeta0_new
    zeta1_1 <- zeta_out$zeta1_new
    par(mfrow = c(1,2))
    plot(zeta0_1, type = "l", main = "zeta0")
    plot(zeta1_1, type = "l", main = "zeta1")
    par(mfrow = c(1,1))
    st4 <- Sys.time()
    st34 <- difftime(st4, st3, units = "mins")
    cat("M-step (HMM) using:", st34, "minutes, done! \n")
    
    pars0 <- c(beta0, rho0, nu0, alpha0)
    pars1 <- update_RL(pars0, A_mat, S_mat, R_mat, gamma_arr[,,2], eva_phi, cores)
    beta1 <- pars1[1]
    rho1 <- pars1[2]
    nu1 <- pars1[3]
    alpha1 <- pars1[4]
    st5 <- Sys.time()
    st45 <- difftime(st5, st4, units = "mins")
    cat("M-step (RL model) using:", st45, "minutes, done! \n")
    cat("beta = ", beta1, "; rho = ", rho1, "alpha = ", alpha1, "nu = ", nu1, "\n")
    
    pi1_all[kk] <- pi11
    zeta_all[kk,,1] <- zeta0_1
    zeta_all[kk,,2] <- zeta1_1
    beta_all[kk] <- beta1
    rho_all[kk] <- rho1
    alpha_all[kk] <- alpha1
    nu_all[kk] <- nu1
    
    ## check jump out conditions
    cond1 <-  abs(pi11 - pi10) < tol
    cond2 <- sqrt(mean((zeta0_1 - zeta0_0)^2)) < tol
    cond3 <- sqrt(mean((zeta1_1 - zeta1_0)^2)) < tol
    cond4 <- abs(nu1 - nu0) < tol
    cond5 <- abs(beta1 - beta0) < tol
    cond6 <- abs(rho1 - rho0) < tol
    cond7 <- abs(alpha1 - alpha0) < tol
    
    if ( cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 ){
      pi1_all <- pi1_all[1:kk]
      zeta_all <- zeta_all[1:kk,,]
      beta_all <- beta_all[1:kk]
      rho_all <- rho_all[1:kk]
      alpha_all <- alpha_all[1:kk]
      nu_all <- nu_all[1:kk]
      break
    }
    pi10 <- pi11
    zeta0_0 <-zeta0_1
    zeta1_0 <-zeta1_1
    beta0 <- beta1
    rho0 <- rho1
    alpha0 <- alpha1
    nu0 <- nu1
  }
  return(list(pi1_all = pi1_all, zeta_all = zeta_all,
              beta_all = beta_all,
              rho_all = rho_all, alpha_all = alpha_all,
              nu_all = nu_all,
              gamma_arr = gamma_arr))
}








########## observed likelihood for model evaluation
obs_log_lik <- function(S_mat, A_mat, R_mat, alpha, beta, rho, nu, zeta0, zeta1, pi1){
  N <- ncol(S_mat)
  n <- nrow(S_mat)
  delta0 <- alpha * c(1, -1, 0, 1)
  
  PA_mat <- array(1/2, dim = c(n, N, 2))
  PA_mat[,,2] <- PA_U1_mat(A_mat, S_mat, R_mat, delta0, beta, rho, nu, eva_phi, cores=1)
  C0_arr <- P_trans_nonstat(zeta0, zeta1)
  
  loglik_n <- 0
  for (i in 1:n){
    Pa_i <- PA_mat[i,,]
    aa_out <- forward_pass(Pa_i, pi1, C0_arr)
    loglik_n <- loglik_n + log(sum(aa_out$aa_mat[N,])) + aa_out$digit_aa[N] * log(10)
  }
  
  return(loglik_n)
}