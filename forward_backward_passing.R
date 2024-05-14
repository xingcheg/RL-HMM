######################################################################
####  Functions for the forward-backward algorithm in the E-step  ####
######################################################################

# give a number (e.g. 1234 = 1.234 * 10^3), return (1.234, 3)
find_digit <- function(x){
  d <- floor(log10(abs(x)))
  x1 <- x * 10^(-d)
  return(c(x1, d))
}



# give two number (e.g. (1234, 123)), return (1.234, 0.123, 3)
find_digit2 <- function(x){
  x_max <- max(x)
  y <- x / x_max
  x_max_out <- find_digit(x_max)
  a <- x_max_out[1]
  d <- x_max_out[2]
  x1 <- y * a
  return(c(x1, d))
}



## forward passing
forward_pass <- function(Pa, pi1, C_array){
  ###############################################################
  # This is the function to compute (a_jt) in the forward passing,
  # for a single subject. The output is an N x 2 matrix saves the 
  # numbers of (a_jt), and a length-N vector saves the digits of 
  # (a_jt), where N is the trial length, 2 is two learning mode.
  #
  # Pa = P(At | Ut = j) is an N x 2 matrix; 
  # pi1 = P(U1 = 1);
  # C_array[t,j,k] = P(Ut+1 = k-1 | Ut = j-1) is the (N-1) x 2 x 2
  # non-stationary transition probability array.
  ###############################################################
  N <- nrow(Pa)
  pi0 <- 1 - pi1
  aa_mat <- matrix(0, N, 2)
  dd_vec <- rep(0, N)
  a10 <- Pa[1,1] * pi0
  a11 <- Pa[1,2] * pi1
  a1_x <- find_digit2(c(a10, a11))
  aa_mat[1,] <- a1_x[1:2]
  dd_vec[1] <- a1_x[3]
  for (t in 2:N){
    at0 <- 0
    at1 <- 0
    for (k in 1:2){
      at0 <- at0 + aa_mat[t-1,k] * Pa[t,1] * C_array[t-1,k,1]
      at1 <- at1 + aa_mat[t-1,k] * Pa[t,2] * C_array[t-1,k,2]
    }
    at_x <- find_digit2(c(at0, at1))
    aa_mat[t,] <- at_x[1:2]
    dd_vec[t] <- at_x[3]
  }
  digit_aa <- cumsum(dd_vec)
  return(list(aa_mat = aa_mat, digit_aa = digit_aa))
}



## backward passing
backward_pass <- function(Pa, C_array){
  ###############################################################
  # This is the function to compute (b_jt) in the backward passing,
  # for a single subject. The output is an N x 2 matrix saves the 
  # numbers of (b_jt), and a length-N vector saves the digits of 
  # (b_jt), where N is the trial length, 2 is two learning mode.
  #
  # Pa = P(At | Ut = j) is an N x 2 matrix; 
  # C_array[t,j,k] = P(Ut+1 = k-1 | Ut = j-1) is the (N-1) x 2 x 2
  # non-stationary transition probability array.
  ###############################################################
  N <- nrow(Pa)
  bb_mat <- matrix(0, N, 2)
  dd_vec <- rep(0, N)
  bb_mat[N,] <- rep(1, 2)
  dd_vec[N] <- 0
  for (t in (N-1):1){
    bt0 <- 0
    bt1 <- 0
    for (k in 1:2){
      bt0 <- bt0 + bb_mat[t+1,k] * Pa[t+1,k] * C_array[t,1,k]
      bt1 <- bt1 + bb_mat[t+1,k] * Pa[t+1,k] * C_array[t,2,k]
    }
    bt_x <- find_digit2(c(bt0, bt1))
    bb_mat[t,] <- bt_x[1:2]
    dd_vec[t] <- bt_x[3]
  }
  digit_bb <- (cumsum(dd_vec[N:1]))[N:1]
  return(list(bb_mat = bb_mat, digit_bb = digit_bb))
}



## function to compute gamma_jt = P( Ut=j | A[1:T] )
comp_gamma <- function(aa_mat, bb_mat, digit_aa, digit_bb){
  ###############################################################
  # This is the function to compute (gamma_jt) for a single subject. 
  # The output is an N x 2 matrix.
  ###############################################################
  N <- nrow(aa_mat)
  log10_D <- log10(sum(aa_mat[N,])) + digit_aa[N]
  log10_N_a <- log10(aa_mat) + digit_aa 
  log10_N_b <- log10(bb_mat) + digit_bb
  log10_N <- log10_N_a + log10_N_b
  gamma_mat <- 10^(log10_N - log10_D)
  return(gamma_mat)
}



## function to compute xi_jkt = P( Ut+1=k, Ut=j | A[1:T] )
comp_xi <- function(aa_mat, bb_mat, digit_aa, digit_bb, Pa, C_array){
  ###############################################################
  # This is the function to compute (xi_jkt) for a single subject. 
  # The output is an (N-1) x 2 x 2 array.
  ###############################################################
  N <- nrow(aa_mat)
  log10_D <- log10(sum(aa_mat[N,])) + digit_aa[N]
  xi_array <- array(0, dim = c(N-1, 2, 2))
  for (j in 1:2){
    for (k in 1:2){
      log10_N_a <- log10(aa_mat[-N,j]) + digit_aa[-N]
      log10_N_b <- log10(bb_mat[-1,k]) + digit_bb[-1]
      log10_N_c <- log10(Pa[-1,k]) 
      log10_N_d <- log10(C_array[,j,k]) 
      log10_N <- log10_N_a + log10_N_b + log10_N_c + log10_N_d
      xi_array[,j,k] <- 10^(log10_N - log10_D)
    }
  }
  return(xi_array)
}


