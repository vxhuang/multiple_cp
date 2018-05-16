## Changepoint detection with one or two changepoints

library(dplyr)
library(forecast)

setClass(Class = "changepoint_result",
         representation(
           plotDT = "data.frame",
           n_cp = "numeric",
           slope = "numeric", 
           intercepts = "numeric",
           c = "numeric"
         )
)


amtc <- function(dt, input_no) {
  
  # sanity check and data cleaning
  dt <- as.data.frame(dt)
  if(ncol(dt) == 1) {
    if(!is.numeric(dt[, 1])) stop("Data should be numeric. ")
    dt <- data_frame("no" = 1 : nrow(dt), "cn" = dt[, 1])
  } else if(ncol(dt) == 2) {
    names(dt) <- c("no", "cn")
    if(is.numeric(dt$no)) {
      if(!all.equal(dt$no, as.integer(dt$no))) stop("If input data have two columns, first column of input data should be integer indices. ")
    } else {
      tryCatch({
        dt$no <- as.Date(dt$no)
        first_date <- min(dt$no)
        dt$no <- as.integer(dt$no - first_date + 1)
      },error=function(e){cat("Error: If input data have two columns, first column of input data should be integer indices or dates. \n")})
    }
    if(is.unsorted(dt$no)) {
      # warning("If input data have two columns, first column of input data should be sorted. ")
      dt <- dt %>% arrange(no)
    }
    
    dt_tmp <- data.frame("no" = matrix(1:max(dt$no), nrow = max(dt$no), ncol = 1))
    dt_tmp <- dt_tmp %>% left_join(dt, by = "no")
    if(is.na(dt_tmp[1, "cn"])) dt_tmp[1, "cn"] <- 0
    for(i in 2 : (nrow(dt_tmp)-1)) {
      if(is.na(dt_tmp[i, "cn"])) {          # if there is no data for some day, we manually fill it
        if(is.na(dt_tmp[i+1, "cn"])) {      # if no data in next day, use 0
          dt_tmp[i, "cn"] <- 0
        } else {                            # if there exists data in the next day, use the data of the previous day
          dt_tmp[i, "cn"] <- dt_tmp[i-1, "cn"]
        }
      }
    }
    if(is.na(dt_tmp[nrow(dt_tmp), "cn"])) dt_tmp[nrow(dt_tmp), "cn"] <- dt_tmp[nrow(dt_tmp)-1, "cn"]
    dt <- dt_tmp
    
  } else stop("Input should not have more than 2 columns. ")
  
  
  n <- max(dt$no)
  
  
  # defind cost function: rmse
  costf <- function(c) {
    # see references for details
    beta_hat_c <- (sum(dt[1:c[1], "cn"] * (1:c[1] - (c[1]+1)/2)) + 
                     sum(dt[(c[1]+1):c[2], "cn"] * ((c[1]+1):c[2] - (c[2]*(c[2]+1)-c[1]*(c[1]+1))/(2*(c[2]-c[1])))) + 
                     sum(dt[(c[2]+1):n, "cn"] * ((c[2]+1):n - (n*(n+1)-c[2]*(c[2]+1))/(2*(n-c[2]))))) /
      (c[1]*(c[1]+1)*(c[1]-1)/12 + (c[2]-c[1])*(c[2]-c[1]+1)*(c[2]-c[1]-1)/12 + (n-c[2])*(n-c[2]+1)*(n-c[2]-1)/12)
    mu_hat_c <- sum(dt[1:c[1], "cn"] - beta_hat_c * 1:c[1]) / c[1]
    delta1_hat_c <- sum(dt[(c[1]+1):c[2], "cn"] - beta_hat_c * (c[1]+1):c[2]) / (c[2]-c[1]) - mu_hat_c
    delta2_hat_c <- sum(dt[(c[2]+1):n, "cn"] - beta_hat_c * (c[2]+1):n) / (n-c[2]) - mu_hat_c - delta1_hat_c
    
    lm.cp <- function(x) {
      mu_hat_c + beta_hat_c * x + delta1_hat_c * ifelse(x > c[1], 1, 0) + delta2_hat_c * ifelse(x > c[2], 1, 0)
    }
    dt.fit <- lm.cp(dt$no)
    cost <- accuracy(dt.fit, dt$cn)[2] # rmse
  }
  
  # GA starts here
  maxiter <- 10
  n_pop <- 200
  # step <- 5
  elite <- 0.2
  mutprob <- 0.5
  h <- .05
  # mutation operation
  mutate <- function(c, step) {
    c_prime <- c
    index <- sample(1:2, 1)
    c[index] <- c[index] + sample(-step : step, 1)
    if(c[1] >= round(n*h) && c[2] <= round(n * (1-h)) && (c[2] - c[1]) >= max(c[2] * 0.05, (n-c[1])*0.05)) {
      return(c)
    } else {
      return(c_prime)
    }
  }
  
  # crossover operation
  crossover <- function(c1, c2) {
    if((c2[2] - c1[1]) >= max(c2[2] * 0.05, (n-c1[1])*0.05)) {
      return(c(c1[1], c2[2]))
    } else if((c1[2] - c2[1]) >= max(c1[2] * 0.05, (n-c2[1])*0.05)) {
      return(c(c2[1], c1[2]))
    } else if(sample(1:2, 1) == 1) {
      return(c1)
    } else {
      return(c2)
    }
  }
  
  # build the initial population
  pop <- list()
  i_pop = 1
  while(length(pop) < n_pop) {
    pop_tmp <- sample(1:n, 2)
    if(pop_tmp[1] >= round(n*h) && pop_tmp[2] <= round(n * (1-h)) && (pop_tmp[2] - pop_tmp[1]) >= max(pop_tmp[2] * 0.05, (n-pop_tmp[1])*0.05)) {
      pop[[i_pop]] <- pop_tmp
      i_pop <- i_pop + 1
    }
  }
  
  # number of winners for each generation
  topelite <- round(n_pop * elite)
  
  # main loop
  # print(format(Sys.time(), tz = "Asia/Taipei"))
  for(i in 1 : maxiter) {
    scores <- lapply(pop, costf)
    scores <- unlist(scores)
    elite_ind <- order(scores)
    
    # start with the pure winners
    elite_pop <- pop[elite_ind[1:topelite]]
    pop <- elite_pop
    i_pop <- topelite + 1
    
    # add mutated and bred forms of the winners
    while(length(pop) < n_pop) {
      if(runif(1) < mutprob) {
        # mutation
        index <- sample(1:topelite, 1)
        pop[[i_pop]] <- mutate(elite_pop[[index]], step = maxiter - i + 1)
      } else {
        # crossover
        index <- sample(1:topelite, 2)
        pop[[i_pop]] <- crossover(elite_pop[[index[1]]], elite_pop[[index[2]]])
      }
      i_pop <- i_pop + 1
    }
    # print(pop[[1]])
  }
  # print(format(Sys.time(), tz = "Asia/Taipei"))
  
  c <- pop[[1]]
  
  # T-statistic from reference (Lund) to test if a changepoint is significant
  T_c <- function(dt, c) {
    beta_hat_c <- (sum(dt[1:c, "cn"] * (1:c - (c+1)/2)) + sum(dt[(c+1):n, "cn"] * ((c+1):n - (n*(n+1)-c*(c+1))/(2*(n-c))))) /
      (c*(c+1)*(c-1)/12 + (n-c)*(n-c+1)*(n-c-1)/12)
    mu_hat_c <- sum(dt[1:c, "cn"] - beta_hat_c * 1:c) / c
    delta_hat_c <- sum(dt[(c+1):n, "cn"] - beta_hat_c * (c+1):n) / (n-c) - mu_hat_c
    sigma_2 <- sum((dt$cn - mean(dt$cn))^2) / n
    var_delta <- sigma_2 / (c * (1-c/n) * (1-3*c*(n-c)/(n^2-1)))
    return(delta_hat_c^2 /var_delta)
  }
  test <- sapply(1:n, T_c, dt = dt)
  # print(c(test[c[1]], test[c[2]]))
  
  rss_full <- min(scores)^2 * n
  rss_red <- sum(lm(dt$cn ~ dt$no)$residual^2)
  F_c <- ((rss_red - rss_full) / 2) /
    (rss_full / (n - 4))
  
  if(F_c > qf(.99, 3, n - 4)) {    # two changepoints confirmed
    beta_hat_c <- (sum(dt[1:c[1], "cn"] * (1:c[1] - (c[1]+1)/2)) +
                     sum(dt[(c[1]+1):c[2], "cn"] * ((c[1]+1):c[2] - (c[2]*(c[2]+1)-c[1]*(c[1]+1))/(2*(c[2]-c[1])))) +
                     sum(dt[(c[2]+1):n, "cn"] * ((c[2]+1):n - (n*(n+1)-c[2]*(c[2]+1))/(2*(n-c[2]))))) /
      (c[1]*(c[1]+1)*(c[1]-1)/12 + (c[2]-c[1])*(c[2]-c[1]+1)*(c[2]-c[1]-1)/12 + (n-c[2])*(n-c[2]+1)*(n-c[2]-1)/12)
    mu_hat_c <- sum(dt[1:c[1], "cn"] - beta_hat_c * 1:c[1]) / c[1]
    delta1_hat_c <- sum(dt[(c[1]+1):c[2], "cn"] - beta_hat_c * (c[1]+1):c[2]) / (c[2]-c[1]) - mu_hat_c
    delta2_hat_c <- sum(dt[(c[2]+1):n, "cn"] - beta_hat_c * (c[2]+1):n) / (n-c[2]) - mu_hat_c - delta1_hat_c
    
    lm.cp <- function(x) {
      mu_hat_c + beta_hat_c * x + delta1_hat_c * ifelse(x > c[1], 1, 0) + delta2_hat_c * ifelse(x > c[2], 1, 0)
    }

    return(new("changepoint_result", 
               plotDT = data.frame("no" = dt$no, "fitted" = lm.cp(dt$no)),
               n_cp = 2, 
               slope = beta_hat_c,
               intercepts = c(mu_hat_c, delta1_hat_c, delta2_hat_c), 
               c = c))
    
  } else if(max(test[round(n*h) : round(n * (1-h))]) > 20.114) {          # one changepoint confirmed
    c <- which(test == max(test[round(n*h) : round(n * (1-h))]))[1]
    beta_hat_c <- (sum(dt[1:c, "cn"] * (1:c - (c+1)/2)) + sum(dt[(c+1):n, "cn"] * ((c+1):n - (n*(n+1)-c*(c+1))/(2*(n-c))))) /
      (c*(c+1)*(c-1)/12 + (n-c)*(n-c+1)*(n-c-1)/12)
    mu_hat_c <- sum(dt[1:c, "cn"] - beta_hat_c * 1:c) / c
    delta_hat_c <- sum(dt[(c+1):n, "cn"] - beta_hat_c * (c+1):n) / (n-c) - mu_hat_c
    
    lm.cp <- function(x) {
      mu_hat_c + beta_hat_c * x + delta_hat_c * ifelse(x > c, 1, 0)
    }
    
    return(new("changepoint_result", 
               plotDT = data.frame("no" = dt$no, "fitted" = lm.cp(dt$no)),
               n_cp = 1, 
               slope = beta_hat_c,
               intercepts = c(mu_hat_c, delta_hat_c), 
               c = c))
    
  } else {
    return(new("changepoint_result", 
               plotDT = data.frame("no" = dt$no, "fitted" = lm(cn ~ no, dt)$fitted),
               n_cp = 0, 
               slope = lm(cn ~ no, dt)$coefficients[2],
               intercepts = lm(cn ~ no, dt)$coefficients[1], 
               c = numeric(0)))

  }
  return(result)
}


