#load packages
library(data.table)
library(ggplot2)
library(ggfortify)
library(abcrf)
library(EasyABC)
library(abc)
library(tuneRanger)

#get files
files <- list.files("data/abm_output")

#for each file
for(i in 1:length(files)){
  #load the file
  load(paste0("data/abm_output/", files[i]))

  #store summary statistics and priors
  if(i == 1){
    sum_stats <- simulations$sum_stats
    priors <- simulations$priors
  }
  if(i > 1){
    sum_stats <- rbindlist(list(sum_stats, simulations$sum_stats))
    priors <- rbindlist(list(priors, simulations$priors))
  }

  #remove temporary object
  rm(simulations)
}

#load the observed tweet distribution
load("data/tweet_distribution.RData")

#calculate observed summary statistics
obs_stats <- data.frame(prop_rare = length(which(tweet_distribution == 1))/sum(tweet_distribution),
                        prop_common = max(tweet_distribution)/sum(tweet_distribution),
                        hill_1 = hillR::hill_taxa(tweet_distribution, q = 1),
                        hill_2 = hillR::hill_taxa(tweet_distribution, q = 2))

#detect number of cores for parallelization
ncores <- detectCores()-1

#set sample size for random forest (100% of data)
sampsize <- 1.00*nrow(sum_stats)

#construct object to hold predictions
predictions <- list()

#set value of i
i <- 1

#run the following code (between the ##########'s) in a manual loop four times

##########

#set random seed
set.seed(i)

#construct data frame for random forest abc
abcrf_data <- data.frame(param = priors[, ..i], sum_stats = sum_stats)
colnames(abcrf_data)[1] <- "param"
colnames(obs_stats) <- colnames(abcrf_data)[-1]

#tuning for best random forest values
task <- makeRegrTask(data = abcrf_data, target = "param")
tuning <- tuneRanger(task, num.trees = 100, parameters = list(sample.fraction = sampsize/nrow(abcrf_data), mtry = 2), tune.parameters = c("min.node.size"), num.threads = ncores)

#run random forest with recommended values
reg_abcrf <- regAbcrf(formula = param ~ ., data = abcrf_data, ntree = 1000, mtry = 2, min.node.size = tuning$recommended.pars$min.node.size, sampsize = sampsize, paral = TRUE, ncores = ncores)

#construct and save posterior distribution
densityPlot(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores, ylab = "Density", xlab = "Parameter Value", main = "Title")
density_plot <- recordPlot()

#collect x and y values
x <- density_plot[[1]][[9]][[2]][[2]]$x
y <- density_plot[[1]][[9]][[2]][[2]]$y

#restrict to prior ranges (remove plotting tails)
if(i %in% c(1, 2, 4)){
  y <- y[which(x >= 0 & x <= 8)]
  x <- x[which(x >= 0 & x <= 8)]
}
if(i %in% c(3)){
  y <- y[which(x >= 0 & x <= 2)]
  x <- x[which(x >= 0 & x <= 2)]
}

#add all output to the predictions object
predictions[[i]] <- list(OOB_MSE = reg_abcrf$model.rf$prediction.error, OOB_NMAE = reg_abcrf$model.rf$NMAE, prediction = predict(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores), var_importance = sort(reg_abcrf$model.rf$variable.importance, decreasing = TRUE), tuning = tuning$recommended.pars, posterior = data.frame(x, y), plot = NULL)
predictions[[i]]$plot <- density_plot

#remove temporary objects
rm(list = c("abcrf_data", "reg_abcrf", "density_plot", "x", "y"))

#add 1 to i
i <- i + 1

##########

#save output
save(predictions, file = "data/predictions.RData")
