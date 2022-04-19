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

#sample priors from posteriors
priors <- data.frame(cont_bias = sample(priors$cont_bias, n_sim, replace = TRUE, prob = predictions[[1]]$weights),
                     dem_bias = sample(priors$dem_bias, n_sim, replace = TRUE, prob = predictions[[2]]$weights),
                     freq_bias = sample(priors$freq_bias, n_sim, replace = TRUE, prob = predictions[[3]]$weights),
                     age_dep = sample(priors$age_dep, n_sim, replace = TRUE, prob = predictions[[4]]$weights))

#run simulations
slurm <- slurm_apply(twitter_ABM_slurm, priors, jobname = "twitter",
                     nodes = 5, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))

#get output and clean files
sum_stats <- get_slurm_out(slurm)
cleanup_files(slurm)

#restructure summary statistics
sum_stats <- data.table(do.call(rbind, sum_stats))
colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

#run simulations again but get the distributions
slurm <- slurm_apply(twitter_ABM_dist_slurm, priors, jobname = "twitter",
                     nodes = 5, cpus_per_node = 48, global_objects = objects(),
                     slurm_options = list(mem = 0))

#get output and clean files
dists <- get_slurm_out(slurm)
cleanup_files(slurm)

#restructure distributions
dists <- do.call(rbind, dists)

#structure and save the output
simulations <- list(priors = priors, sum_stats = sum_stats, dists = dists)
save(simulations, file = "data/abm_output/posterior_simulations.RData")
