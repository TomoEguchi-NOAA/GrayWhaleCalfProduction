# =====================================================================
# Replicated MCMC-breakdown sensitivity analysis
# ---------------------------------------------------------------------
# Spreads a fixed total count (2026: 26 pairs) across an increasing number
# of weeks and, over many random data realizations, records how often each
# model fails to converge. Output: P(non-convergence) vs number of occupied
# weeks, for the Poisson-Binomial and Negative Binomial models.
#
# This turns the single-realization result (Poisson broke at 6-8 weeks but
# "recovered" at 9) into a dose-response curve. The apparent recovery at 9
# is almost certainly one draw's luck near a stochastic threshold; averaging
# over replicates removes it and shows the true failure probability.
#
# COMPUTE COST: n_rep * n.weeks * 2 model fits. With n_rep = 100 and 9 weeks
# that is 1800 fits — hours to overnight. START WITH A PILOT (n_rep = 20) to
# see the curve shape, then scale up. Results are checkpointed after every
# replicate, so a long run can be resumed (see the resume block below).
# =====================================================================

rm(list = ls())
library(jagsUI)
library(tidyverse)
library(posterior)
library(cmdstanr)

source("GrayWhaleCalfProduction_fcns_v2.R")

# ---- settings ------------------------------------------------------
n_rep      <- 100 #2 #20 #          # <- start with 20 as a pilot
year       <- 2026
rhat_cut   <- 1.01
ess_cut    <- 400
ckpt_file  <- "RData\\mcmc_breakdown_sweep_checkpoint.rds"
jm_data_file <- "RData//mcmc_breakdown_jags_data.rds"

# Lighter MCMC for the sweep: enough for an identified model to converge
# cleanly, while structural non-convergence (R-hat >> 1) still shows at once.
# If the curve looks noisy near the threshold, raise n.samples.
MCMC.params <- list(n.samples = 100000, n.thin = 100,
                    n.burnin = 50000, n.chains = 5)

# JAGS: monitor only what is needed to detect breakdown (lambda is where it
# shows). Fewer monitored nodes = faster to summarize, smaller objects.
jags.params <- c("lambda", "p.obs", "Total.Calves")

# ---- load the 2026 data as the template ----------------------------
jm.out       <- readRDS(paste0("RData\\calf_estimates_v3_Mv1_", year, ".rds"))
jags.data    <- jm.out$jags.data
total.counts <- sum(jags.data$count.obs)
n.weeks      <- length(jags.data$weekly.max)
all.weeks    <- 1:n.weeks

count.data <- data.frame(count  = jags.data$count.obs,
                         effort = jags.data$effort,
                         week   = jags.data$week) %>%
  mutate(log_offset = ifelse(effort > 0, log(effort / 3.0), 0))

# ---- compile the Stan model ONCE (was inside the loop before) ------
mod_ <- cmdstan_model("models//GWCalfCount_nb_singleyear.stan",
                      cpp_options = list(stan_threads = TRUE, O = 3))

# ---- helpers -------------------------------------------------------
# spread `total.counts` across k weeks, fanning out from week 5
make_week_vec <- function(k) {
  wv <- numeric(k)
  for (x in 1:k) wv[x] <- 5 + (((-1)^(x - 1)) * (2 * x - 1) - 1) / 4
  wv
}

# build one simulated count.obs vector: assign counts to weeks (multinomial
# over weeks), then within a week distribute over intervals proportional to effort
simulate_counts <- function(week.vec, jags.data, total.counts) {
  assignments <- sample(as.factor(week.vec), size = total.counts, replace = TRUE) %>%
    table() %>% as.data.frame()
  colnames(assignments) <- c("week", "Counts")
  
  wk <- data.frame(week = as.factor(week.vec)) %>%
    left_join(assignments, by = "week") %>%
    mutate(week.num = as.numeric(as.character(week)),
           Counts   = ifelse(is.na(Counts), 0, Counts)) %>%
    arrange(week.num)
  
  count.obs <- jags.data$count.obs * 0
  for (r in seq_len(nrow(wk))) {
    wn      <- wk$week.num[r]
    eff.tmp <- jags.data$effort[jags.data$week == wn]
    if (sum(eff.tmp) == 0 || wk$Counts[r] == 0) next
    count.obs[jags.data$week == wn] <-
      rmultinom(1, size = wk$Counts[r], prob = eff.tmp / sum(eff.tmp))
  }
  count.obs
}

# max R-hat / min bulk ESS over a set of draws
conv_stats <- function(summary_df) {
  list(max_rhat = max(summary_df$rhat,     na.rm = TRUE),
       min_bulk_ess  = min(summary_df$ess_bulk, na.rm = TRUE),
       min_tail_ess = min(summary_df$ess_tail, na.rm = TRUE))
}

# ---- resume from checkpoint if present -----------------------------
if (file.exists(ckpt_file)) {
  results <- readRDS(ckpt_file)
  jm.data <- readRDS(jm_data_file)
  message("Resuming: ", length(results), " cells already done.")
} else {
  results <- list()
  jm.data <- list()
}

done_key <- vapply(results, function(z) paste(z$rep[1], z$n_weeks[1]), character(1))

# ---- main sweep ----------------------------------------------------
for (rep in 1:n_rep) {
  tic <- Sys.time()
  for (k in 1:n.weeks) {
    
    key <- paste(rep, k)
    if (key %in% done_key) next                 # already computed; skip
    
    set.seed(1000 * rep + k)                     # per-cell reproducible data
    week.vec  <- make_week_vec(k)
    count.obs <- simulate_counts(week.vec, jags.data, total.counts)
    
    jags.sim.data <- jags.data
    jags.sim.data$count.obs <- count.obs
    wk.min.max <- data.frame(week = jags.sim.data$week, count = count.obs) %>%
      group_by(week) %>% summarise(max = max(count), min = min(count), .groups = "drop")
    jags.sim.data$weekly.max <- wk.min.max$max
    
    # ---------- Poisson-Binomial (JAGS) ----------
    poi <- tryCatch({
      jm <- jags(jags.sim.data, inits = NULL,
                 parameters.to.save = jags.params,
                 "models\\GWCalfCount_v1.jags",
                 n.chains = MCMC.params$n.chains, n.burnin = MCMC.params$n.burnin,
                 n.thin   = MCMC.params$n.thin,   n.iter   = MCMC.params$n.samples,
                 DIC = FALSE, parallel = TRUE, verbose = FALSE)
      conv_stats(summarize_draws(as_draws_df(jm$samples)))
    }, error = function(e) list(max_rhat = NA_real_, min_ess = NA_real_))
    
    # ---------- Negative Binomial (Stan) ----------
    nb <- tryCatch({
      stan_data <- list(
        n_obs      = length(count.obs),
        n_weeks    = jags.sim.data$n.weeks,
        count_obs  = count.obs,
        effort     = jags.sim.data$effort,
        log_offset = count.data$log_offset,
        week_idx   = jags.sim.data$week)
      fit_ <- mod_$sample(data = stan_data, seed = 12,
                          chains = 4, parallel_chains = 4,
                          iter_warmup = 1000, iter_sampling = 2000,
                          threads_per_chain = 2, init = 0.1,
                          adapt_delta = 0.99, refresh = 0)
      conv_stats(fit_$summary(c("p_obs", "sigma_week", "phi", "beta0", "week_eff")))
    }, error = function(e) list(max_rhat = NA_real_, min_ess = NA_real_))
    
    results[[length(results) + 1]] <- data.frame(
      rep = rep, n_weeks = k,
      max_weekly = max(wk.min.max$max),
      Model    = c("Poisson-Binomial", "Negative Binomial"),
      max_rhat = c(poi$max_rhat, nb$max_rhat),
      min_bulk_ess = c(poi$min_bulk_ess,  nb$min_bulk_ess),
      min_tail_ess = c(poi$min_tail_ess,  nb$min_tail_ess))
    
    jm.data[[length(jm.data) + 1]] <- jags.sim.data
            
    done_key <- c(done_key, key)
  }
  
  saveRDS(results, ckpt_file)                # checkpoint after each replicate
  saveRDS(jm.data, jm_data_file)
  cat("replicate", rep, "of", n_rep, "done  (", Sys.Date(), ")\n")
  toc <- Sys.time() - tic
  cat("Took ", toc, " ", attributes(toc)$units,  "\n")
}

# ---- assemble + convergence flag -----------------------------------
results_df <- do.call(rbind, results) %>%
  mutate(converged = (max_rhat < rhat_cut) & 
           (min_bulk_ess > ess_cut) & 
           (min_tail_ess > ess_cut))

conv_rate <- results_df %>%
  group_by(Model, n_weeks) %>%
  summarise(p_fail = mean(!converged, na.rm = TRUE),
            n_ok   = sum(!is.na(converged)), .groups = "drop")

# ---- plot: P(non-convergence) vs count spread ----------------------
p_break <- ggplot(conv_rate, aes(n_weeks, p_fail, colour = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = all.weeks) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Number of weeks the counts were spread across",
       y = "P(non-convergence)",
       title = sprintf("MCMC breakdown vs count spread (%d counts, %d replicates)",
                       total.counts, n_rep)) +
  theme_bw() +
  theme(legend.position = "top") 

print(p_break)

# optional second panel: mean max-R-hat vs spread, log scale, shows the cliff
p_rhat <- results_df %>%
  group_by(Model, n_weeks) %>%
  summarise(med_rhat = median(max_rhat, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(n_weeks, med_rhat, colour = Model)) +
  geom_hline(yintercept = rhat_cut, linetype = 2) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_x_continuous(breaks = all.weeks) +
  scale_y_log10() +
  labs(x = "Number of weeks the counts were spread across",
       y = expression("median max " * hat(R)),
       title = "Convergence diagnostic vs count spread") +
  theme_bw()+
  theme(legend.position = "top") 
