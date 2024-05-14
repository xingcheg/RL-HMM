################################################################
####  Functions to simulate RL-HMM data in the simulations  ####
################################################################
## State space:[0,1]; Action space:{0,1}; Reward space:[0,1]. 


## mechanism of getting rewards
reward_mechanism <- function(a, s, p0=0.5, p1=1){
  ### a = 0 or 1; s in [0, 1].
  if (a==0){
    R <- rbinom(1, 1, p0*(1-s))
  }
  if (a==1){
    R <- rbinom(1, 1, p1*s)
  }
  return(R)
}



## linear basis function phi(a,s)
eva_phi <- function(a, s){
  psi <- c(1, s)
  f0 <- ifelse(a==0, 1, 0)
  f1 <- 1 - f0
  return(c(f0 * psi, f1 * psi))
}



## update coefficient delta
update_delta <- function(a, s, R, delta0, beta, rho, eva_phi){
  phi0 <- eva_phi(a,s)
  dd <- rho * R - sum(delta0 * phi0)
  delta1 <- delta0 + beta * dd * phi0
  return(delta1)
}



## probability of taking action 1 (RL model)
action_RL <- function(s, nu, delta, eva_phi){
  phi_1 <- eva_phi(1,s)
  phi_0 <- eva_phi(0,s)
  Q_diff <- sum((phi_1 - phi_0) * delta)
  prob <- 1 / (1 + exp(-Q_diff-nu))
  a <- rbinom(1, 1, prob) 
  return(list(a = a, p = prob, q_diff = Q_diff))
}



## HMM transition probability
comp_C <- function(N, zeta0, zeta1){
  C01 <- 1/(1+exp(-zeta0))
  C11 <- 1/(1+exp(-zeta1))
  return(list(C01 = C01, C11 = C11))
}



## generate HMM latent variable
gen_U <- function(u0, c01, c11){
  u1 <- ifelse(u0==1, 
               rbinom(1, 1, c11),
               rbinom(1, 1, c01))
  return(u1)
} 



## generate actions for one subject
gen_A_single <- function(states, beta, rho, nu, C01, C11, delta0, 
                         u0, eva_phi, gen_reward){
  N <- length(states)
  m <- length(delta0)
  actions <- rep(0, N)
  rewards <- rep(0, N)
  ## HMM indicator (0 for ST model, 1 for RL model)
  U <- rep(0, N)
  U[1] <- u0
  ## Q function coefficients
  delta_chain <- matrix(0, N, m)
  delta_chain[1,] <- delta0
  for (t in 1:N){
    s <- states[t]
    if (u0==0){
      a <- rbinom(1, 1, 0.5)
    }
    if (u0==1){
      a <- action_RL(states[t], nu, delta0, eva_phi)$a
    }
    r <- gen_reward(a, s)   ## user specified reward generating function:
    ## (R_t) is derived by (A_t and S_t)
    actions[t] <- a
    rewards[t] <- r
    if (t < N){
      ## update delta (RL model coefficients)
      delta1 <- update_delta(a, s, r, delta0, beta, rho, eva_phi)
      delta_chain[t+1,] <- delta1
      delta0 <- delta1
      ## update U (HMM indicator)
      u1 <- gen_U(u0, C01[t], C11[t])
      U[t+1] <- u1
      u0 <- u1
    }
  }
  out <- list(actions = actions, rewards = rewards,
              U = U, delta_chain = delta_chain)
}



## generate actions for multiple subjects
gen_A_multi <- function(n, N, S_mat, beta, rho, nu, pi1, 
                        zeta0, zeta1, alpha, eva_phi, 
                        gen_reward){
  delta0 <- alpha * c(1, -1, 0, 1)
  A_mat <- matrix(0, n, N)
  R_mat <- matrix(0, n, N)
  U_mat <- matrix(0, n, N)
  Delta_chain <- array(0, dim = c(n, N, length(delta0)))
  C_arr <- comp_C(N, zeta0, zeta1)
  C01 <- C_arr$C01
  C11 <- C_arr$C11
  U0 <- rbinom(n, 1, pi1)
  for (i in 1:n){
    states <- S_mat[i,]
    out <- gen_A_single(states, beta, rho, nu, C01, C11, delta0, 
                        U0[i], eva_phi, gen_reward)
    A_mat[i,] <- out$actions
    R_mat[i,] <- out$rewards
    U_mat[i,] <- out$U
    Delta_chain[i,,] <- out$delta_chain
  }
  Output <- list(S_mat = S_mat, A_mat = A_mat,
                 R_mat = R_mat, U_mat = U_mat,
                 Delta_chain = Delta_chain)
  return(Output)
}



