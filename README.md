# RL-HMM


R code for the simulation study (Case I, n=100, T=100, method = 'RL-HMM') from the paper "Reinforcement Learning with Hidden Markov Models for Discovering Decision-Making Dynamics". by Xingche Guo, Donglin Zeng, and Yuanjia Wang.

<img src="RL-HMM diagram.png" alt="Model Structure" width="800"/>


## Overview


* The .R file **"simulation_funcs.R"** contains functions to simulate RL-HMM data.

* The .R file **"forward_backward_passing.R"** contains core functions for the Forward and Backward Passing algorithm.

* The .R files **"RL_HMM_algo.R"** contains the EM algorithm for parameter estimation of RL-HMM.

* The .R file **"example.R"** is an example showing how to simulate RL-HMM data and estimate parameters using our proposed EM algorithm.

* The .R file **"run.R"** contains functions to produce results (Case I, n=100, T=100, method='RL-HMM') in Table 1 of the paper. Please change 'lam_rank=1' if we want to produce results for method 'RL-HMM-fixed'. Please use cluster for parallel computing. The nonparametric bootstrap is time-consuming and hence is omiited.

* the .rds file **"RL_HMM_est_100_100.rds"** contains simulation results for (Case I, n=100, T=100, method='RL-HMM') in Table 1 of the paper.
