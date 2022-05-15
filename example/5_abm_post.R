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
load("data/predictions.RData")
load("data/abm_output/prior_table.RData")

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

#extract posteriors and trim to prior
cont_bias_post <- density(priors[, 1], weights = predictions[[1]]$weights, from = 0, to = 12)
dem_bias_post <- density(priors[, 2], weights = predictions[[2]]$weights, from = 0, to = 4)
freq_bias_post <- density(priors[, 3], weights = predictions[[3]]$weights, from = 0, to = 2)
age_dep_post <- density(priors[, 4], weights = predictions[[4]]$weights, from = 0, to = 8)

#sample priors from posteriors
post_priors <- data.frame(cont_bias = sample(cont_bias_post$x, n_sim, replace = TRUE, prob = cont_bias_post$y),
                          dem_bias = sample(dem_bias_post$x, n_sim, replace = TRUE, prob = dem_bias_post$y),
                          freq_bias = sample(freq_bias_post$x, n_sim, replace = TRUE, prob = freq_bias_post$y),
                          age_dep = sample(age_dep_post$x, n_sim, replace = TRUE, prob = age_dep_post$y))

#run simulations
slurm <- slurm_apply(twitter_ABM_slurm, post_priors, jobname = "twitter",
                     nodes = 6, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))
