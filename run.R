###############################################################################
####  R code to produce results (Case I, n=100, T=100, RL-HMM) in Table 1  ####
###############################################################################
## Please use server to run the 200 replicates simultaneously.
## Bootstrap (by resampling) is omitted to save space .

source("simulation_funcs.R")
source("forward_backward_passing.R")
source("RL_HMM_algo.R")


slurm_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
set.seed(1234 + slurm_id)

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
                    M_max = 100, tol = 2e-4, cores = 1)


out_name <- paste0("est_out/EM_out_rep_", slurm_id, ".rds")

saveRDS(EM_out, file = out_name)
