# Binary outcomes
dnorm_new<-function(x, mean=0, sd=1, log=F, height=1){
  out<-dnorm(x=x,mean=mean,sd=sd,log=log)
  out*height
}

dnorm_2_new<-function(x, mean_1=0, mean_2=0, sd_1=1, sd_2=1, log=F, p_2=0.5){
  out_1<-dnorm_new(x, mean=mean_1, sd=sd_1, log=log, height=1-p_2)
  out_2<-dnorm_new(x, mean=mean_2, sd=sd_2, log=log, height=p_2)
  out_1+out_2
}

ccprobs.f <- function(d=0.641, prev=0.7463, n_quantile=20){
  mu_case <- d
  mu_control <- 0
  
  varPRS <- prev*(1+(d^2) - (d*prev)^2) + (1-prev)*(1 - (d*prev)^2)
  E_PRS <- d*prev
  
  by_quant<-1/n_quantile
  p_quant <- seq(by_quant, 1-by_quant, by=by_quant)
  quant_vals_PRS <- rep(0, length(p_quant))
  quant_f_solve <- function(x, prev, d, pq){prev*pnorm(x-d) + (1-prev)*pnorm(x) - pq}
  for(i in 1:length(p_quant)){
    quant_vals_PRS[i] <- unlist(uniroot(quant_f_solve, prev=prev, d=d, pq= p_quant[i], interval=c(-2.5, 2.5), extendInt = "yes", tol=6e-12)$root)
  }
  
  ul_qv_PRS <- matrix(0, ncol=2, nrow=n_quantile)
  ul_qv_PRS[1,1] <- -Inf
  ul_qv_PRS[2:n_quantile,1] <- quant_vals_PRS
  ul_qv_PRS[1:(n_quantile-1),2] <- quant_vals_PRS
  ul_qv_PRS[n_quantile,2] <- Inf
  
  ul_qv_PRS<-cbind(ul_qv_PRS, (ul_qv_PRS[,1:2]-E_PRS)/sqrt(varPRS))
  
  prob_quantile_case <- pnorm(ul_qv_PRS[,2], mean = mu_case) - pnorm(ul_qv_PRS[,1], mean = mu_case)
  prob_quantile_control <- pnorm(ul_qv_PRS[,2], mean = mu_control) - pnorm(ul_qv_PRS[,1], mean = mu_control)
  p_case_quantile <- (prob_quantile_case*prev)/by_quant
  p_cont_quantile <- (prob_quantile_control*(1-prev))/by_quant
  
  OR <- p_case_quantile/p_cont_quantile
  OR <- OR/OR[1]
  out <- cbind(ul_qv_PRS[,3:4],p_cont_quantile, p_case_quantile, OR)
  row.names(out) <- 1:n_quantile
  colnames(out) <- c("q_min", "q_max","p_control", "p_case", "OR")
  
  data.frame(out)
}

# Normally distributed outcomes
which_quant <- function(PRS_z_score=0, n_quantile=20){
  E_PRS = 0
  SD_PRS = sqrt(1)
  
  by_quant<-1/(n_quantile)
  PRS_quantile_bounds <- qnorm(p=seq(0, 1, by=by_quant))
  lower_PRS_vec <- PRS_quantile_bounds[1:n_quantile]
  upper_PRS_vec <- PRS_quantile_bounds[2:(n_quantile+1)]
  
  out<-data.frame(q=1:n_quantile,
                  q_min=lower_PRS_vec,
                  q_max=upper_PRS_vec)
  
  quant<-out$q[PRS_z_score > out$q_min & PRS_z_score <= out$q_max]
  return(quant)
}

mean_sd_quant.f <- function(PRS_R2=0.641, Outcome_mean=1, Outcome_sd=1, n_quantile=20, quant=NA){
  ### PRS quantiles with a continuous phenotype (Y)
  library(tmvtnorm)
  ###
  E_PRS = 0
  SD_PRS = sqrt(1)
  E_phenotype = Outcome_mean
  SD_phenotype = Outcome_sd 
  
  by_quant<-1/(n_quantile)
  PRS_quantile_bounds <- qnorm(p=seq(0, 1, by=by_quant), mean= E_PRS, sd= SD_PRS)
  lower_PRS_vec <- PRS_quantile_bounds[1:n_quantile]
  upper_PRS_vec <- PRS_quantile_bounds[2:(n_quantile+1)]
  
  mean_vec <- c(E_phenotype, E_PRS)
  sigma_mat <- matrix(sqrt(PRS_R2)*SD_phenotype*SD_PRS, nrow=2, ncol=2)
  sigma_mat[1,1] <- SD_phenotype^2
  sigma_mat[2,2] <- SD_PRS^2
  
  ### mean of phenotype within the truncated PRS distribution
  out_mean_Y <- rep(0, 20)
  ### SD of phenotype within the truncated PRS distribution
  out_SD_Y <- rep(0, 20)
  ### cov of Y and PRS given truncation on PRS
  out_cov_Y_PRS <- rep(0, 20)
  ### SD of PRS given truncation on PRS
  out_SD_PRS <- rep(0, 20)
  ### mean PRS given truncation on PRS
  out_mean_PRS <- rep(0, 20)
  
  if(!is.na(quant)){
    i<-quant
    
    distribution_i <- mtmvnorm(mean = mean_vec,
                               sigma = sigma_mat,
                               lower = c(-Inf, lower_PRS_vec[i]),
                               upper = c(Inf, upper_PRS_vec[i]),
                               doComputeVariance=TRUE,
                               pmvnorm.algorithm=GenzBretz())
    out_mean_Y[i] <- distribution_i$tmean[1]
    out_mean_PRS[i] <- distribution_i$tmean[2]
    out_SD_Y[i] <- sqrt(distribution_i$tvar[1,1])
    out_SD_PRS[i] <- sqrt(distribution_i$tvar[2,2])
    out_cov_Y_PRS[i] <- distribution_i$tvar[1,2]
    
    out<-data.frame(q=quant,
                    q_min=lower_PRS_vec[quant],
                    q_max=upper_PRS_vec[quant],
                    x_mean=out_mean_Y[quant],
                    x_sd=out_SD_Y[quant])
    
  } else {
    for(i in 1:n_quantile){
      distribution_i <- mtmvnorm(mean = mean_vec,
                                 sigma = sigma_mat,
                                 lower = c(-Inf, lower_PRS_vec[i]),
                                 upper = c(Inf, upper_PRS_vec[i]),
                                 doComputeVariance=TRUE,
                                 pmvnorm.algorithm=GenzBretz())
      out_mean_Y[i] <- distribution_i$tmean[1]
      out_mean_PRS[i] <- distribution_i$tmean[2]
      out_SD_Y[i] <- sqrt(distribution_i$tvar[1,1])
      out_SD_PRS[i] <- sqrt(distribution_i$tvar[2,2])
      out_cov_Y_PRS[i] <- distribution_i$tvar[1,2]
    }
    
    out<-data.frame(q=1:n_quantile,
                    q_min=lower_PRS_vec,
                    q_max=upper_PRS_vec,
                    x_mean=out_mean_Y,
                    x_sd=out_SD_Y)
  }
  
  return(out)
  
  out_mean_Y
  out_SD_Y
  
  out_mean_PRS
  out_SD_PRS
  out_cov_Y_PRS
}
