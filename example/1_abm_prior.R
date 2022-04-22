#load packages
library(TwitterABM)
library(rslurm)
library(data.table)

#set seed
set.seed(12345)

#load data and set parameters
N <- 341676
load("data/obs_user_data.RData")
load("data/obs_init_tweets.RData")
load("data/overall_activity.RData")

#wrap twitter_ABM in a simpler function for slurm
twitter_ABM_slurm <- function(cont_bias, dem_bias, freq_bias, age_dep){
  #run full version of twitter_ABM
  dist <- twitter_ABM(N = N, overall_activity = overall_activity,
                      cont_bias = cont_bias, dem_bias = dem_bias, freq_bias = freq_bias, age_dep = age_dep,
                      obs_user_data = obs_user_data, obs_init_tweets = obs_init_tweets, sum_stats_TF = FALSE)$n_times
  sum_stats <- c(length(which(dist == 1))/sum(dist),
                 max(dist)/sum(dist),
                 hillR::hill_taxa(dist, q = 1),
                 hillR::hill_taxa(dist, q = 2))
  dist <- Rfast::Sort(dist, descending = TRUE)[1:100000]
  list(sum_stats, dist)
}

#number of simulations
n_sim <- 10000

#set priors
priors <- data.frame(cont_bias = runif(n_sim, min = 0, max = 8),
                     dem_bias = runif(n_sim, min = 0, max = 4),
                     freq_bias = runif(n_sim, min = 0, max = 2),
                     age_dep = runif(n_sim, min = 0, max = 8))

#save priors
save(priors, file = "data/abm_output/prior_simulations/prior_table.RData")

#run simulations
slurm <- slurm_apply(twitter_ABM_slurm, priors, jobname = "twitter",
                     nodes = 5, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))
