#load packages
library(data.table)
library(abc)
library(parallel)

#load files
results_0 <- readRDS("data/abm_output/results_0.RDS")
results_1 <- readRDS("data/abm_output/results_1.RDS")
results_2 <- readRDS("data/abm_output/results_2.RDS")
results_3 <- readRDS("data/abm_output/results_3.RDS")
load("data/abm_output/prior_table.RData")
load("data/tweet_distribution.RData")

#combine output
results_0 <- do.call(rbind, results_0)
results_1 <- do.call(rbind, results_1)
results_2 <- do.call(rbind, results_2)
results_3 <- do.call(rbind, results_3)
sum_stats <- do.call(rbind, list(results_0, results_1, results_2, results_3))
colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

#calculate observed summary statistics
obs_stats <- data.frame(prop_rare = length(which(tweet_distribution == 1))/sum(tweet_distribution),
                        prop_common = max(tweet_distribution)/sum(tweet_distribution),
                        hill_1 = hillR::hill_taxa(tweet_distribution, q = 1),
                        hill_2 = hillR::hill_taxa(tweet_distribution, q = 2))

#set bounds of logit transformation
logit_bounds <- matrix(data = c(0, 0, 0, 0, 12, 4, 2, 8), ncol = 2)

#tolerance levels to explore
tols <- c(0.005, 0.001, 0.0005, 0.0001, 0.00005)

#run cross validation to determine optimal tolerance level
cv_rej <- cv4abc(param = priors, sumstat = sum_stats, method = "rejection", nval = 100,
                 tols = tols, transf = c("logit", "logit", "logit", "logit"), logit.bounds = logit_bounds)
save(cv_rej, file = "data/cv_rej.RData")

#run full rejection ABC analysis at each tolerance level
predictions_rej <- list()
predictions_rej <- mclapply(1:length(tols), function(x){
  abc(target = obs_stats, param = priors, sumstat = sum_stats, method = "rejection",
      tol = tols[x], transf = c("logit", "logit", "logit", "logit"), logit.bounds = logit_bounds)
}, mc.cores = 5)
names(predictions_rej) <- tols
save(predictions_rej, file = "data/predictions_rej.RData")
