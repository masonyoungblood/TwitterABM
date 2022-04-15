#load packages
library(data.table)
library(ggplot2)
library(ggfortify)
library(abcrf)
library(EasyABC)
library(abc)
library(tuneRanger)

#load files
results_0 <- readRDS("data/abm_output/results_0.RDS")
results_1 <- readRDS("data/abm_output/results_1.RDS")
results_2 <- readRDS("data/abm_output/results_2.RDS")
results_3 <- readRDS("data/abm_output/results_3.RDS")
results_4 <- readRDS("data/abm_output/results_4.RDS")
load("data/abm_output/prior_table.RData")
load("data/tweet_distribution.RData")

#combine output
results_0 <- do.call(rbind, results_0)
results_1 <- do.call(rbind, results_1)
results_2 <- do.call(rbind, results_2)
results_3 <- do.call(rbind, results_3)
results_4 <- do.call(rbind, results_4)
sum_stats <- do.call(rbind, list(results_0, results_1, results_2, results_3, results_4))
colnames(sum_stats) <- c("prop_rare", "prop_common", "hill_1", "hill_2")

#calculate observed summary statistics
obs_stats <- data.frame(prop_rare = length(which(tweet_distribution == 1))/sum(tweet_distribution),
                        prop_common = max(tweet_distribution)/sum(tweet_distribution),
                        hill_1 = hillR::hill_taxa(tweet_distribution, q = 1),
                        hill_2 = hillR::hill_taxa(tweet_distribution, q = 2))

#detect number of cores for parallelization
ncores <- detectCores()-1

#set sample size for random forest (100% of data)
sampsize <- 0.8*nrow(sum_stats)

#construct object to hold predictions
predictions <- list()

#set value of i
i <- 1

#run the following code (between the ##########'s) in a manual loop four times

##########

#set random seed
set.seed(i)

#construct data frame for random forest abc
abcrf_data <- data.frame(param = priors[, i], sum_stats = sum_stats)
colnames(abcrf_data)[1] <- "param"
colnames(obs_stats) <- colnames(abcrf_data)[-1]

#tuning for best random forest values
task <- makeRegrTask(data = abcrf_data, target = "param")
tuning <- tuneRanger(task, num.trees = 500, parameters = list(sample.fraction = sampsize/nrow(abcrf_data)), tune.parameters = c("mtry", "min.node.size"), num.threads = ncores)

#run random forest with recommended values
reg_abcrf <- regAbcrf(formula = param ~ ., data = abcrf_data, ntree = 1000, mtry = tuning$recommended.pars$mtry, min.node.size = tuning$recommended.pars$min.node.size, sampsize = sampsize, paral = TRUE, ncores = ncores)

# #construct and save posterior distribution
# densityPlot(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores, ylab = "Density", xlab = "Parameter Value", main = "Title")
# density_plot <- recordPlot()

# #collect x and y values
# x <- density_plot[[1]][[9]][[2]][[2]]$x
# y <- density_plot[[1]][[9]][[2]][[2]]$y
#
# #restrict to prior ranges (remove plotting tails)
# if(i %in% c(1, 2, 4)){
#   y <- y[which(x >= 0 & x <= 8)]
#   x <- x[which(x >= 0 & x <= 8)]
# }
# if(i %in% c(3)){
#   y <- y[which(x >= 0 & x <= 2)]
#   x <- x[which(x >= 0 & x <= 2)]
# }

#add all output to the predictions object
predictions[[i]] <- list(OOB_MSE = reg_abcrf$model.rf$prediction.error, OOB_NMAE = reg_abcrf$model.rf$NMAE, prediction = predict(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores), var_importance = sort(reg_abcrf$model.rf$variable.importance, decreasing = TRUE), tuning = tuning$recommended.pars)
#predictions[[i]] <- list(OOB_MSE = reg_abcrf$model.rf$prediction.error, OOB_NMAE = reg_abcrf$model.rf$NMAE, prediction = predict(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores), var_importance = sort(reg_abcrf$model.rf$variable.importance, decreasing = TRUE), tuning = tuning$recommended.pars, posterior = data.frame(x, y), plot = NULL)
#predictions[[i]]$plot <- density_plot

save(predictions, file = "data/predictions.RData")

#add 1 to i
i <- i + 1

##########

#save output
save(predictions, file = "data/predictions.RData")
