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
  twitter_ABM(N = N, overall_activity = overall_activity,
              cont_bias = cont_bias, dem_bias = dem_bias, freq_bias = freq_bias, age_dep = age_dep,
              obs_user_data = obs_user_data, obs_init_tweets = obs_init_tweets)
}

#number of simulations
n_sim <- 10000

#set priors
priors <- data.frame(cont_bias = runif(n_sim, min = 0, max = 8),
                     dem_bias = runif(n_sim, min = 0, max = 8),
                     freq_bias = runif(n_sim, min = 0, max = 2),
                     age_dep = runif(n_sim, min = 0, max = 8))

#run simulations
slurm <- slurm_apply(twitter_ABM_slurm, priors, jobname = "twitter",
                     nodes = 4, cpus_per_node = 48, global_objects = objects())

#get output and clean files
sum_stats <- get_slurm_out(slurm)
cleanup_files(slurm)

#restructure summary statistics
sum_stats <- data.table(do.call(rbind, sum_stats))
colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

#structure and save the output
simulations <- list(priors = priors, sum_stats = sum_stats)
save(simulations, file = "data/abm_output/prior_simulations.RData")
