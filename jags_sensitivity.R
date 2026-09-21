# Conducts a sensitivity analysis on counts-effort and MCMC breakdown
# 

rm(list = ls())
library(jagsUI)
library(tidyverse)
library(posterior)
library(cmdstanr)

source("GrayWhaleCalfProduction_fcns_v2.R")

MCMC.params <- list(n.samples = 100000,
                    n.thin = 100,
                    n.burnin = 50000,
                    n.chains = 5)

n.samples <- MCMC.params$n.chains * ((MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin)

jags.params <- c("count.true",
                 "lambda", "psi", "alpha",
                 "beta1", "beta2", "eps",
                 "p.obs.corr",
                 "p.obs",
                 "Total.Calves",
                 "loglik")

year <- 2026 # This didn't work so well

out.file.name <- paste0("RData\\calf_estimates_v3",
                        "_Mv1_", year, ".rds") 

jm.out <- readRDS(file = out.file.name)

jags.data <- jm.out$jags.data

total.counts <- sum(jags.data$count.obs)
non.zero.weeks <- length(which(jags.data$weekly.max > 0))
n.weeks <- length(jags.data$weekly.max)
all.weeks <- c(1:n.weeks)

count.data <- data.frame(count = jags.data$count.obs,
                         effort = jags.data$effort,
                         week = jags.data$week) %>%
  mutate(log_offset = ifelse(effort > 0, log(effort / 3.0), 0))

weekly.sum <- count.data %>%
  group_by(week) %>%
  summarize(count = sum(count),
            effort = sum(effort),
            week = first (week))

jags.sim.data <- jags.data

weekly.counts <- week.vec.list <- list()
count.obs.list <- jm.out.sim <- list()
stan.out <- global_summary <- diag_summary <- list()
k <- 5
for (k in 1:n.weeks){
  week.vec <- vector(mode = "numeric", length = k)
  for (x in 1:k){
    week.vec[x] <- 5 + (((-1)^(x-1)) * (2*x - 1) - 1)/4 
  }
  
  assignments <- sample(as.factor(week.vec), 
                        size = total.counts,
                        replace = TRUE) %>% 
    table() %>% 
    data.frame() 
  
  colnames(assignments) <- c("week", "Counts")
  # assignments %>%
  #   mutate(week = as.numeric(week.f)) -> assignments
  # 
  weekly.counts[[k]] <- data.frame(week = as.factor(week.vec)) %>%
    left_join(assignments, by = "week") %>%
    mutate(week.num = as.numeric(as.character(week))) %>%
    arrange(by = week.num)
  
  week.vec.list[[k]] <- week.vec
  #n.per.week <- ceiling(total.counts/k)
  
  # create jags data with new weekly counts
  count.obs <- jags.data$count.obs * 0
  for (k1 in 1:nrow(weekly.counts[[k]])){
    effort.tmp <- jags.data$effort[jags.data$week == weekly.counts[[k]][k1, "week.num"]]
    effort.prob <- effort.tmp/sum(effort.tmp)
    count.tmp <- rmultinom(n = 1, size = weekly.counts[[k]][k1, "Counts"],
                           prob = effort.prob)
    
    count.obs[jags.data$week == weekly.counts[[k]][k1, "week.num"]] <- count.tmp
    
  }
  
  jags.sim.data$count.obs <- count.obs
  tmp.data <- data.frame(week = jags.sim.data$week,
                         count = jags.sim.data$count.obs)
  weekly.max <- tmp.data %>%
    group_by(week) %>%
    summarise(week = first(week),
              max = max(count))
  jags.sim.data$weekly.max <- weekly.max$max
  
  jm <- jags(jags.sim.data,
             inits = NULL,
             parameters.to.save = jags.params,
             "models\\GWCalfCount_v1.jags", 
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.thin = MCMC.params$n.thin,
             n.iter = MCMC.params$n.samples,
             DIC = T, parallel=T)
  
  # This function is in GrayWhaleCalfProduction_fcns.R
  # params.to.monitor is a string of regular expression.
  # e.g., "^BF\\.Fixed|^K\\["
  params <- "^count\\.true\\[|^lambda\\[|^p\\.obs\\.corr\\[|^p\\.obs|^Total\\.Calves\\["
       
  jm.MCMC <- MCMC.diag(jm = jm, 
                       MCMC.params = MCMC.params,
                       params.to.monitor = params)
  
  jm.out.sim[[k]] <- list(jm = jm,
                          MCMC.diag = jm.MCMC,
                          jags.data = jags.sim.data,
                          MCMC.params = MCMC.params,
                          run.date = Sys.Date())     
  

  
  stan.data <- stan_data <- list(
    n_obs      = length(jags.sim.data$count.obs),
    #n_years    = max(data.1year$year_idx),
    n_weeks    = jags.sim.data$n.weeks,
    count_obs  = jags.sim.data$count.obs,
    effort     = jags.sim.data$effort, # Passing raw effort to Stan
    log_offset = count.data$log_offset, 
    #year_idx   = data.1year$year_idx,
    week_idx   = jags.sim.data$week
    )
  
  model.file <- "models//GWCalfCount_nb_singleyear.stan"
  
  mod_ <- cmdstan_model(stan_file = model.file,
                        cpp_options = list(stan_threads = TRUE, 
                                           O = 3))
  tic <- Sys.time()
  fit_ <- mod_$sample(
    data            = stan_data,
    seed            = 12,
    chains          = 4,
    parallel_chains = 4,
    iter_warmup     = 1000,
    threads_per_chain = 2,
    iter_sampling   = 2000,
    init            = 0.1,
    adapt_delta     = 0.99,
    refresh         = 500
  )
  toc <- Sys.time() -  tic
  global_summary[[k]] <- fit_$summary(c("p_obs", "sigma_week", 
                                        "phi", "beta0", "week_eff"))
  
  
  #print(global_summary)
  # stan.out[[k]] <- stan.post.process.1year(
  #   stan.fit = fit_, 
  #   pre.stan.data = jags.data,
  #   stan.data = stan_data, 
  #   out.file.name = out.file,
  #   save.file = F)
  # 
  #PPC.out[[k]] <- PPC_counts(stan.fit = fit_,
  #                           stan.data = stan_data)
  
  diag_summary[[k]] <- fit_$diagnostic_summary()
  
}

## Summarize jags output:
posterior.summary <- list()
for (k in 1:n.weeks){
  posterior.summary[[k]] <- posterior::summarize_draws(jm.out.sim[[k]]$jm$samples) %>%
    mutate(Week = all.weeks[k])
}

posterior.summary.df <- do.call(rbind, posterior.summary) 

posterior.summary.df %>%
  group_by(Week) %>%
  summarise(Max.Rhat = max(rhat, na.rm = T),
            Min.ESS.bulk = min(ess_bulk, na.rm = T),
            Min.ESS.tail = min(ess_tail, na.rm = T)) %>%
  mutate(Model = "Poisson-Binomial") -> Rhat.ESS.summary.Jags

## Summarize Stan output
for (k in 1:n.weeks){
  global_summary[[k]] <- global_summary[[k]] %>%
    mutate(Week = all.weeks[k])
}

global_summary_df <- do.call(rbind, global_summary)

global_summary_df %>%
  group_by(Week) %>%
  summarise(Max.Rhat = max(rhat, na.rm = T),
            Min.ESS.bulk = min(ess_bulk, na.rm = T),
            Min.ESS.tail = min(ess_tail, na.rm = T)) %>%
  mutate(Model = "Negative Binomial") -> Rhat.ESS.summary.Stan

Rhat.ESS.summary <- rbind(Rhat.ESS.summary.Jags,
                          Rhat.ESS.summary.Stan)

