# this script contains function to estiamte Bayesian probit GLM for 
# individual participants and calculate Bayes factor
#
# - note that it requires stan code be in the file `probit_blink.stan`
#
# Matteo Lisi, 2019


# Note that in this version the 
bayesian_GLM_varSD <-function(d_i, prior_sd){
  
  d_i <- d_i[d_i$cond1==1,]
  
  d_stan <- list(
    Response = d_i$Response,
    JumpSize = d_i$JumpSize,
    vel1 = ifelse(d_i$Velocity=="240", 1, 0),
    bdur = d_i$bdur,
    prior_sd =  prior_sd,
    N = nrow(d_i)
  )
  
  require(rstan)
  require(mlisi)
  options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
  rstan_options(auto_write = TRUE)
  
  m_stan <- stan(file = "./probit_blink.stan", data = d_stan, iter = 2000, chains = 4, init=0)
  beta <- extract(m_stan, pars=c("beta"))$"beta"
  print(m_stan)
  BF01 <- savage.dickey.bf(beta[,5], x_0 = 0, prior.mean = 0, prior.sd = prior_sd, plot = T)
  Bi <- apply(beta,2,mean)
  
  # PSE: -beta_0/beta_1
  mu180 <- NA
  mu240 <- NA
  mu180_blink <- -Bi[1]/Bi[2]
  mu240_blink <- -(Bi[1]+Bi[3])/(Bi[2]+Bi[4])
  
  # SIGMA: 1/beta_1
  sig180 <- NA
  sig240 <- NA
  sig180_blink <- 1/(Bi[2])
  sig240_blink <- 1/(Bi[2]+Bi[4])
  
  outV <- data.frame(t(c(unique(d_i$Subject), BF01, Bi, mu180, mu240, mu180_blink, mu240_blink, sig180, sig240, sig180_blink, sig240_blink)))
  colnames(outV) <- c("id", "BF01", paste("beta_",1:5,sep=""),"mu180", "mu240", "mu180_blink", "mu240_blink", "sig180", "sig240", "sig180_blink", "sig240_blink")
  return(outV)
}
