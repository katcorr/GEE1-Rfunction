#### Lipsitz et al. One-step GEE with Large Clusters

onestepGEE <- function(dataset, outcome, cluster, predictors){
  # ------------------------------------------------------------------------------
  # ---------------------------- set parameters ----------------------------------
  # ------------------------------------------------------------------------------
  
  # identify dataset and clustering variable here, so only need to change in 
  # one place if update; all further instances refer to 'dat', 'outcomevar, and 'clustervar'
  dat <- dataset %>%
    rename(clustervar := {{ cluster }}
           , outcomevar := {{ outcome }}) %>%
    filter(!is.na(outcomevar)) %>%
    arrange(clustervar)
  
  # total number of observations
  n <- nrow(dat)
  
  # number of observations per cluster
  n_i <- as.vector(unlist(dat %>% count(clustervar) %>% select(n)))
  
  # number of clusters
  c <- length(n_i)
  
  # indices for first and last row of observations in each cluster in dataset
  n_end <- cumsum(n_i)
  n_start <- lag(n_end) + 1  
  n_start[1] <- 1
  
  # model information
  model_formula <- paste0("outcomevar ~ ", predictors)
  
  # set p below (number of model parameters, including intercept)
  
  # outcome vector
  y_obs <- dat$outcomevar 
  
  # ------------------------------------------------------------------------------
  # --------------- fit naive model to estimate rho then beta_{gee1} -------------
  # ------------------------------------------------------------------------------
  
  # estimate beta under naive independence (standard logistic regression)
  model_naive <- glm(as.formula(model_formula), data = dat, family=binomial(link="logit"))
  
  # model matrix
  model_X <- model.matrix(model_naive)
  
  # number of model parameters
  p <- ncol(model_X)
  
  # get estimated error terms (Pearson residual) based on naive beta:
  # \hat{e}_{ij} = (Y_{ij} - \hat{mu}_{ij})/sqrt{\hat{v}_{ij}}.
  mu_hat_naive <- predict(model_naive, type = "response")
  #ggplot(data=as.data.frame(naive_mu), aes(x=mu_hat_naive)) + geom_histogram()
  
  # scale (dispersion) parameter is taken to be 1; var of bernoulli is p(1-p)
  var_y_naive <- mu_hat_naive*(1-mu_hat_naive)
  
  e_hat_naive <- (y_obs - mu_hat_naive) / sqrt(var_y_naive)
  #e_hat_naive_compare <- resid(model_naive, type="pearson")
  
  # ------------------------------------------------------------------------------
  # --------------------------- estimate rho based on naive beta  ----------------
  # ---------------------------         (equation 10)             ----------------
  # ------------------------------------------------------------------------------
  
  rho1 <- dat %>%
    mutate(mu_hat_naive = mu_hat_naive
           , var_y_naive = var_y_naive
           , e_hat_naive = e_hat_naive) %>%
    group_by(clustervar) %>%
    summarize(n=n()
              , mean=mean(e_hat_naive)
              , sum=sum(e_hat_naive)
              # uss (uncorrected sum of squares)
              , uss=sum(e_hat_naive^2)) %>%
    mutate(wt_ij = n*(n-1)/2
           , rho_ij = (sum^2 - uss)/2) %>%
    ungroup() 
  
  rho_hat_naive <- sum(rho1$rho_ij)/sum(rho1$wt_ij)
  
  print(rho_hat_naive)
  
  # ------------------------------------------------------------------------------
  # ---------------------   estimate W based on naive beta  ----------------------
  # ---------------------         (equation 8)              ----------------------
  # ------------------------------------------------------------------------------
  
  linear_naive <- predict(model_naive, type = "link")
  
  d_ij_multiplier <- (exp(-linear_naive)/(1+exp(-linear_naive))^2)
  
  #"# d_ij is a given column of D_i [with as many rows as there are model parameters]"
  # for our purposes:
  # D is a pxn matrix where each COLUMN is d_ij
  D_naive <- t(model_X*d_ij_multiplier)
  #dim(D_naive)
  
  W_naive_part1 <- array(dim=c(p,p,n))
  W_naive_part2a <- array(dim=c(p,n))
  for (k in 1:n){
    W_naive_part1[,,k] <- D_naive[,k] %*% t(D_naive[,k]) / var_y_naive[k]
    W_naive_part2a[,k] <- D_naive[,k] / sqrt(var_y_naive[k]) 
  }
  
  d_ij_v_sum <- array(dim=c(p,c))
  W_i_naive <- array(dim=c(p,p,c))
  for (i in 1:c){
    # this will be used again below, so save it
    d_ij_v_sum[,i] <- apply(W_naive_part2a[,n_start[i]:n_end[i],drop=F],c(1),sum) 
    W_i_naive[,,i] <- (apply(W_naive_part1[,,n_start[i]:n_end[i],drop=F],c(1,2),sum) - 
                         (rho_hat_naive/((1-rho_hat_naive) + 
                                           (n_i[i]*rho_hat_naive)))*(d_ij_v_sum[,i] %*% t(d_ij_v_sum[,i]))) 
  }
  
  # sum the matrices across clusters (W_1 + W_2 + ... + W_51)
  W_naive_tot <- apply(W_i_naive,c(1,2),sum)
  
  
  # ------------------------------------------------------------------------------
  # --------------------- estimate beta hat gee1  --------------------------------
  # ---------------------     (equation 12)       --------------------------------
  # ------------------------------------------------------------------------------
  
  beta_hat_naive <- model_naive$coeff
  
  eq12_a <- rho_hat_naive/((1-rho_hat_naive) + (n_i*rho_hat_naive))
  
  eq12_c <- array(dim=c(c))
  for (i in 1:c){
    eq12_c[i] <- sum(e_hat_naive[n_start[i]:n_end[i]])
  }
  
  eq12_abc <- array(dim=c(p,c))
  for (i in 1:c){
    eq12_abc[,i] <- eq12_a[i]*d_ij_v_sum[,i]*eq12_c[i]
  }
  
  eq12_part2 <- apply(eq12_abc,c(1),sum)
  
  beta_hat_GEE1 <- beta_hat_naive - solve(W_naive_tot)%*%eq12_part2
  
  # ------------------------------------------------------------------------------
  # ----------------- estimate asypm cov for beta hat gee ------------------------
  # ------------------           (equation 7)             ------------------------
  # ------------------------------------------------------------------------------
  
  # "The asymptotic covariance of beta hat GEE1 is consistently estimates by (7)
  #   evaluated at beta hat GEE1" : NEED TO COMPUTE ALL VALUES UNDER beta hat GEE1
  
  # ----------------- estimate rho under beta hat GEE1 ---------------------------
  
  linear_hat_gee1 <- model_X%*%beta_hat_GEE1
  mu_hat_gee1 <- exp(linear_hat_gee1)/(1+exp(linear_hat_gee1))
  #ggplot(data=as.data.frame(mu_gee1), aes(x=mu_gee1)) + geom_histogram()
  
  # scale (dispersion) parameter is taken to be 1; var of bernoulli is p(1-p)
  var_y_gee1 <- mu_hat_gee1*(1-mu_hat_gee1)
  
  e_hat_gee1 <- (y_obs - mu_hat_gee1) / sqrt(var_y_gee1)
  
  rho2 <- dat %>%
    mutate(mu_hat_gee1 = mu_hat_gee1
           , var_y_gee1 = var_y_gee1
           , e_hat_gee1 = e_hat_gee1) %>%
    group_by(clustervar) %>%
    summarize(n=n()
              , mean=mean(e_hat_gee1)
              , sum=sum(e_hat_gee1)
              # uss (uncorrected sum of squares)
              , uss=sum(e_hat_gee1^2)) %>%
    mutate(wt_ij = n*(n-1)/2
           , rho_ij = (sum^2 - uss)/2) %>%
    ungroup() 
  
  rho_hat_gee1 <- sum(rho2$rho_ij)/sum(rho2$wt_ij)
  
  print(rho_hat_gee1)
  
  # ----------------- estimate D and W under beta hat GEE1 -----------------------
  
  d_ij_multiplier_gee1 <- as.vector((exp(-linear_hat_gee1)/(1+exp(-linear_hat_gee1))^2))
  
  #"# d_ij is a given column of D_i [with as many rows as there are model parameters]"
  # for our purposes:
  # D is a pxn matrix where each COLUMN is d_ij
  D_gee1 <- t(model_X*d_ij_multiplier_gee1)
  #dim(D_gee1)
  
  W_gee1_part1 <- array(dim=c(p,p,n))
  W_gee1_part2a <- array(dim=c(p,n))
  for (k in 1:n){
    W_gee1_part1[,,k] <- D_gee1[,k] %*% t(D_gee1[,k]) / var_y_gee1[k]
    W_gee1_part2a[,k] <- D_gee1[,k] / sqrt(var_y_gee1[k]) 
  }
  
  d_ij_v_sum_gee1 <- array(dim=c(p,c))
  W_i_gee1 <- array(dim=c(p,p,c))
  for (i in 1:c){
    # this will be used again below, so save it
    d_ij_v_sum_gee1[,i] <- apply(W_gee1_part2a[,n_start[i]:n_end[i],drop=F],c(1),sum) 
    W_i_gee1[,,i] <- (apply(W_gee1_part1[,,n_start[i]:n_end[i],drop=F],c(1,2),sum) - 
                        (rho_hat_gee1/((1-rho_hat_gee1) + 
                                         (n_i[i]*rho_hat_gee1)))*(d_ij_v_sum_gee1[,i] %*% t(d_ij_v_sum_gee1[,i]))) 
  }
  
  # sum the matrices across clusters (W_1 + W_2 + ... + W_51)
  W_gee1_tot <- apply(W_i_gee1,c(1,2),sum)
  
  # ----------------- estimate scores under beta hat GEE1 -----------------------
  
  score_part1 <- array(dim=c(p,n))
  score_part2_b <- array(dim=c(p,n))
  for (k in 1:n){
    score_part1[,k] <- D_gee1[,k]*(y_obs[k] - mu_hat_gee1[k])/var_y_gee1[k]
    score_part2_b[,k] <- D_gee1[,k]/sqrt(var_y_gee1[k])
  }
  
  score_part1_sums <- array(dim=c(p,c))
  score_part2_b_sum <- array(dim=c(p,c))
  score_part2_c <- array(dim=c(c))
  for (i in 1:c){
    score_part1_sums[,i] <- apply(score_part1[,n_start[i]:n_end[i],drop=F],c(1),sum)
    score_part2_b_sum[,i] <- apply(score_part2_b[,n_start[i]:n_end[i],drop=F],c(1),sum)
    score_part2_c[i] <- sum(e_hat_gee1[n_start[i]:n_end[i]]) 
  }
  
  score_i <- array(dim=c(p,c))
  for (i in 1:c){
    score_i[,i] <- (score_part1_sums[,i] - (rho_hat_gee1/((1-rho_hat_gee1) +
                                                            n_i[i]*rho_hat_gee1))*score_part2_b_sum[,i]*score_part2_c[i])
  }
  
  sandwich_mid <- array(dim=c(p,p,c))
  for (i in 1:c){
    sandwich_mid[,,i] <- score_i[,i]%*%t(score_i[,i])
  }
  
  sandwich_mid_sum <- apply(sandwich_mid, c(1,2), sum)
  
  # ----------------- compute sandwich estimator -----------------------
  
  cov_beta_hat_gee1 <- solve(W_gee1_tot)%*%sandwich_mid_sum%*%solve(W_gee1_tot)
  
  # comparison <- data.frame(beta_hat_naive,beta_hat_GEE1=beta_hat_GEE1)
  
  # comparison[,3:4] <- data.frame(naive_se = summary(model_naive)$coeff[,2]
  #                                , gee1_se = sqrt(diag(cov_beta_hat_gee1)))
  
  # --------------------- return a list object with beta hats, covariance matrix etc. ---------
  
  namevec <- names(summary(model_naive)$coeff[,2])
  beta_hat_GEE1 <- as.vector(beta_hat_GEE1)
  names(beta_hat_GEE1) <- namevec
  se_gee1 <- sqrt(diag(cov_beta_hat_gee1))
  names(se_gee1) <- namevec
  
  rownames(cov_beta_hat_gee1) <- namevec
  colnames(cov_beta_hat_gee1) <- namevec
  
  results <- list(beta_hat_gee1 = beta_hat_GEE1
                  , se_gee1 = se_gee1
                  , covar_mat_gee1 = cov_beta_hat_gee1
                  , beta_hat_standard_logistic = beta_hat_naive
                  , se_standard_logistic = summary(model_naive)$coeff[,2]
                  , icc = rho_hat_gee1
                  , dataset = as.character(substitute(dataset))
                  , outcome = as.character(substitute(outcome))
                  , cluster = as.character(substitute(cluster))
                  , predictors = as.character(substitute(predictors))
                  )
  
  return(results)
}
