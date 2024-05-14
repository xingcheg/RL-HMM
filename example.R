###########################################################################
####  An example to show how to simulate data and estimate parameters  ####
###########################################################################
## let lam_rank = 1 for RL-HMM-fixed method.
## let lam_rank > 1 for RL-HMM method.


source("simulation_funcs.R")
source("forward_backward_passing.R")
source("RL_HMM_algo.R")


set.seed(1234)

################################################################
####################### simulate data ######################
################################################################
n <- 100
N <- 100
## RL model parameters
beta <- 0.05
rho <- 4
nu <- 0
alpha <- 2
S_mat <- matrix(rbeta(n*N, 1, 1), n, N)
## HMM parameters
pi1 <- 0.8
zeta0 <- c(rep(-1.5, N/2), rep(-3, N/2 - 1))
zeta1 <- c(rep(2, N/2), rep(4, N/2 - 1))



Simu_data <- gen_A_multi(n=n, N=N, S_mat=S_mat, beta=beta, rho=rho, 
                         nu=nu, pi1=pi1, zeta0=zeta0, zeta1=zeta1, 
                         alpha=alpha, eva_phi=eva_phi, 
                         gen_reward=reward_mechanism)

S_mat <- Simu_data$S_mat
A_mat <- Simu_data$A_mat
R_mat <- Simu_data$R_mat


################################################################################
####################### estimate parameters ######################
################################################################################


## starting values
pi10 <- 0.75
zeta0_0 <- rep(-2, N-1)
zeta1_0 <- rep(2, N-1)
rho0 <- 3
beta0 <- 0.03
alpha0 <- 1.5
nu0 <- 0
lam_rank <- 6
order <- 1



EM_out <- EM_RL_HMM(S_mat, A_mat, R_mat, beta0, rho0, alpha0, nu0,
                    pi10, zeta0_0, zeta1_0, order, lam_rank, eva_phi,
                    M_max = 100, tol = 2e-4, cores = 10)


saveRDS(EM_out, file = "example_est.rds")

