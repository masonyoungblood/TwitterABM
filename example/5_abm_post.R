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
  twitter_ABM(N = N, overall_activity = overall_activity,
              cont_bias = cont_bias, dem_bias = dem_bias, freq_bias = freq_bias, age_dep = age_dep,
              obs_user_data = obs_user_data, obs_init_tweets = obs_init_tweets)
}

#wrap twitter_ABM in a simpler function for slurm, to output the first 100,000 in the distribution
twitter_ABM_dist_slurm <- function(cont_bias, dem_bias, freq_bias, age_dep){
  #run full version of twitter_ABM, sorted and subsetted to only have top 100,000
  temp <- twitter_ABM(N = N, overall_activity = overall_activity,
                      cont_bias = cont_bias, dem_bias = dem_bias, freq_bias = freq_bias, age_dep = age_dep,
                      obs_user_data = obs_user_data, obs_init_tweets = obs_init_tweets, sum_stats_TF = FALSE)
  temp <- Rfast::Sort(temp$n_times, descending = TRUE)[1:100000]
  return(temp)
}

#number of simulations
n_sim <- 1000

#extract posteriors and trim to prior
cont_bias_post <- density(priors[, 1], weights = predictions[[1]]$weights, from = 0, to = 8)
dem_bias_post <- density(priors[, 2], weights = predictions[[2]]$weights, from = 0, to = 8)
freq_bias_post <- density(priors[, 3], weights = predictions[[3]]$weights, from = 0, to = 2)
age_dep_post <- density(priors[, 4], weights = predictions[[4]]$weights, from = 0, to = 8)

#sample priors from posteriors
post_priors <- data.frame(cont_bias = sample(cont_bias_post$x, n_sim, replace = TRUE, prob = cont_bias_post$y),
                          dem_bias = sample(dem_bias_post$x, n_sim, replace = TRUE, prob = dem_bias_post$y),
                          freq_bias = sample(freq_bias_post$x, n_sim, replace = TRUE, prob = freq_bias_post$y),
                          age_dep = sample(age_dep_post$x, n_sim, replace = TRUE, prob = age_dep_post$y))

#run simulations
slurm <- slurm_apply(twitter_ABM_slurm, post_priors, jobname = "twitter",
                     nodes = 5, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))

#get output and clean files
sum_stats <- get_slurm_out(slurm)
cleanup_files(slurm)

#restructure summary statistics
sum_stats <- data.table(do.call(rbind, sum_stats))
colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

#run simulations again but get the distributions
slurm <- slurm_apply(twitter_ABM_dist_slurm, post_priors, jobname = "twitter",
                     nodes = 5, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))

#get output and clean files
dists <- get_slurm_out(slurm)
cleanup_files(slurm)

#restructure distributions
dists <- do.call(rbind, dists)

#structure and save the output
simulations <- list(priors = post_priors, sum_stats = sum_stats, dists = dists)
save(simulations, file = "data/abm_output/posterior_simulations.RData")
