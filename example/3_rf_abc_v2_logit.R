#load packages
library(data.table)
library(abcrf)
library(tuneRanger)
library(rslurm)

#modify function from abcrf::densityPlot.regAbcrf for weight extraction
extract_weights <- function(object, obs, training, paral=FALSE, ncores= if(paral) max(detectCores()-1,1) else 1, ...){
  x <- obs

  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula

  mf$data <- training

  training <- mf$data

  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  resp <- model.response(mf)

  obj <- object$model.rf
  inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)

  obj[["origNodes"]] <- predict(obj, training, predict.all=TRUE, num.threads=ncores)$predictions
  obj[["origObs"]] <- model.response(mf)

  #####################

  origObs <- obj$origObs
  origNodes <- obj$origNodes

  nodes <- predict(obj, x, predict.all=TRUE, num.threads=ncores )$predictions
  if(is.null(dim(nodes))) nodes <- matrix(nodes, nrow=1)
  ntree <- obj$num.trees
  nobs <- object$model.rf$num.samples
  nnew <- nrow(x)

  weights <- abcrf:::findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
  weights.std <- weights/ntree

  return(weights.std[, 1])
}

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

#logit transform all priors (functions adopted from abc package)
logit <- function(param, logit_bounds){
  temp <- (param - logit_bounds[1])/(logit_bounds[2] - logit_bounds[1])
  return(log(temp/(1 - temp)))
}
inv_logit <- function(param, logit_bounds){
  temp <- exp(param)/(1 + exp(param))
  return((temp*(logit_bounds[2] - logit_bounds[1])) + logit_bounds[1])
}
priors[, 1] <- logit(priors[, 1], c(0, 8))
priors[, 2] <- logit(priors[, 2], c(0, 8))
priors[, 3] <- logit(priors[, 3], c(0, 2))
priors[, 4] <- logit(priors[, 4], c(0, 8))

#set sample size for random forest (80% of data)
sampsize <- 0.8*nrow(sum_stats)

#wrap random forest loop in a simpler function for rslurm
random_forest_slurm <- function(i, title){
  #set random seed
  set.seed(i)

  #detect cores
  ncores <- 20

  #construct data frame for random forest abc
  abcrf_data <- data.frame(param = priors[, i], sum_stats = sum_stats)
  colnames(abcrf_data)[1] <- "param"
  colnames(obs_stats) <- colnames(abcrf_data)[-1]

  #tuning for best random forest values
  task <- makeRegrTask(data = abcrf_data, target = "param")
  tuning <- tuneRanger(task, num.trees = 500, parameters = list(sample.fraction = sampsize/nrow(abcrf_data)), tune.parameters = c("mtry", "min.node.size"))

  #run random forest with recommended values
  reg_abcrf <- regAbcrf(formula = param ~ ., data = abcrf_data, ntree = 1000, mtry = tuning$recommended.pars$mtry, min.node.size = tuning$recommended.pars$min.node.size, sampsize = sampsize, paral = TRUE, ncores = ncores)

  #return predictions
  list(OOB_MSE = reg_abcrf$model.rf$prediction.error, OOB_NMAE = reg_abcrf$model.rf$NMAE, prediction = predict(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores),
       var_importance = sort(reg_abcrf$model.rf$variable.importance, decreasing = TRUE), tuning = tuning$recommended.pars, weights = extract_weights(object = reg_abcrf, obs = obs_stats, training = abcrf_data, paral = TRUE, ncores = ncores))
}

#set up params data frame
params <- data.frame(i = c(1:4), title = names(priors))

#run simulations without angles
slurm <- slurm_apply(random_forest_slurm, params, jobname = "abcrf",
                     nodes = 4, cpus_per_node = 1, global_objects = objects(),
                     slurm_options = list(mem = 0))

#get and save output
predictions <- get_slurm_out(slurm)
save(predictions, file = "data/predictions_logit.RData")

#cleanup files
cleanup_files(slurm)
