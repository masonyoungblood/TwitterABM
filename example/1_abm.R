#load packages
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(TwitterABM)

#load data and set parameters
N <- 341676
mu <- 0.4454465
load("data/overall_activity.RData")
load("data/obs_follower_counts.RData")
load("data/obs_activity_levels.RData")
load("data/obs_init_tweets.RData")

#reduce the skew of obs_follower_counts by removing anything above 100000
obs_follower_counts <- obs_follower_counts[-which(obs_follower_counts >= 100000)]

#register 50 cores for parallelization
cl <- makeForkCluster(50)
registerDoParallel(cl)

#for each of 100 iterations
foreach(i = 1:100) %dopar% {
  #set random seed
  set.seed(i)

  #set number of iterations to 2000
  iter <- 2000

  #sample from priors
  priors <- data.table(cont_bias = runif(iter, min = 0, max = 8), dem_bias = runif(iter, min = 0, max = 8), freq_bias = runif(iter, min = 0, max = 2), age_dep = runif(iter, min = 0, max = 8))

  #run simulations and save output
  sum_stats <- lapply(1:iter, function(x){twitter_ABM(N = N, mu = mu, overall_activity = overall_activity, cont_bias = priors$cont_bias[x], dem_bias = priors$dem_bias[x], freq_bias = priors$freq_bias[x], age_dep = priors$age_dep[x], obs_follower_counts = obs_follower_counts, obs_activity_levels = obs_activity_levels, obs_init_tweets = obs_init_tweets)})
  sum_stats <- data.table(do.call(rbind, sum_stats))
  colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

  #structure and save the output
  simulations <- list(priors = priors, sum_stats = sum_stats)
  save(simulations, file = paste0("data/abm_output/simulations_", i, ".RData"))
}
