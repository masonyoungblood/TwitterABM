#' @title Agent-based model
#' @description An agent-based model (ABM) of cultural transmission on Twitter that incorporates content bias, frequency bias, demonstrator bias, and the level of age-dependent selection.
#'
#' @param N Overall population size.
#' @param overall_activity A vector of the total of number of tweets and retweets in each timestep. The length of this vector is used to determine the number of timesteps.
#' @param mu Probability of writing an original tweet as opposed to retweeting.
#' @param cont_bias Variation in the salience of the attractiveness of content (only positive values, where 0 is neutrality).
#' @param dem_bias Variation in the salience of the follower count (only positive values, where 0 is neutrality).
#' @param freq_bias Level of frequency bias (only positive values, where < 1 is novelty and > is conformity).
#' @param age_dep Rate of decay in age-dependent selection.
#' @param obs_follower_counts A vector of the observed follower counts, randomly drawn for the simulated users.
#' @param obs_activity_levels A vector of the observed activity levels, randomly drawn for the simulated users.
#' @param obs_init_tweets A vector of the observed retweet frequencies from the first timestep, randomly drawn to initialize the tweet data table.
#' @param sum_stats_TF Whether you want to simplify the raw data to the following summary statistics: (1) the proportion of tweets that only appear once, (2) the proportion of the most common tweet, (3) the Hill number when q = 1 (which emphasizes more rare tweets), and (4) the Hill number when q = 2 (which emphasizes more common tweets) (TRUE/FALSE).
#' @param diversity_TF Whether you want to return the Simpson's diversity index from each timepoint (TRUE/FALSE).
#'
#' @return Returns one object or a list of two objects, depending on the values for sum_stats_TF and diversity_TF.
#' @import data.table dqrng truncnorm hillR vegan
#' @export
twitter_ABM <- function(N = 1000, overall_activity, mu = 0.3, cont_bias = 0, dem_bias = 0, freq_bias = 1, age_dep = 1, obs_follower_counts, obs_activity_levels, obs_init_tweets, sum_stats_TF = TRUE, diversity_TF = FALSE){
  #optimization notes:
    #data tables used instead of data frames (set function does wonders for speed)
      #pre-allocation makes for a marginal improvement, but not much better than using rbindlist
      #pre-filling n_times, age, and content has no significant effect on runtime
      #matrix version is not faster than data table
    #wrswoR::sample_int_crank only does weighted sampling without replacement, and couldn't be used
      #even when used it does not appear to be significantly faster than base R with replacement
    #dqrng was used for unweighted sampling (does not allow probabilities yet)
    #rtruncnorm is fastest option with high sample sizes

  #remove first timestep of overall_activity, which is captured by obs_init_tweets
  overall_activity <- overall_activity[-1]

  #use the length of overall_activity to determine the number of timesteps
  t_step <- length(overall_activity)

  #generate data table of users
  users <- data.table::data.table(id = 1:N, followers = c(abs(scale(dqrng::dqsample(obs_follower_counts, N, replace = TRUE))+1))^dem_bias, activity_level = dqrng::dqsample(obs_activity_levels, N, replace = TRUE))

  #pre-allocate data table for tweets
  nrow_allocation <- sum(overall_activity)
  tweets <- data.table::data.table(user = rep(0, nrow_allocation), n_times = rep(0, nrow_allocation), age = rep(0, nrow_allocation), content = rep(0, nrow_allocation), followers = rep(0, nrow_allocation))
  data.table::set(tweets, i = 1:length(obs_init_tweets), j = "user", value = sample(users$id, length(obs_init_tweets), replace = FALSE, prob = users$activity_level))
  data.table::set(tweets, i = 1:length(obs_init_tweets), j = "n_times", value = dqrng::dqsample(obs_init_tweets, length(obs_init_tweets), replace = TRUE))
  data.table::set(tweets, i = 1:length(obs_init_tweets), j = "age", value = 1)
  data.table::set(tweets, i = 1:length(obs_init_tweets), j = "content", value = truncnorm::rtruncnorm(length(obs_init_tweets), a = 0, b = Inf, mean = 1, sd = 1)^cont_bias)
  data.table::set(tweets, i = 1:length(obs_init_tweets), j = "followers", value = users$followers[tweets$user])

  if(diversity_TF){
    #generate diversity vector
    diversity <- rep(NA, t_step)
  }

  #pre-sample active user (with replacement so that users can be active multiple times per timestep)
  active_user_list <- lapply(1:t_step, function(x){sample.int(N, overall_activity[x], replace = TRUE, prob = users$activity_level)})

  for(i in 1:t_step){
    active_users <- active_user_list[[i]]

    #set fill point for allocation
    fill_point <- which.min(tweets$user)

    #calculate retweet probabilities
    probs <- (tweets$n_times[1:(fill_point-1)]^freq_bias)*(tweets$followers[1:(fill_point-1)])*(tweets$content[1:(fill_point-1)])*(1/(tweets$age[1:(fill_point-1)]^age_dep))

    #figure out which active users will tweet new content
    tweet_events <- as.logical(stats::rbinom(length(active_users), 1, prob = mu))
    num_events <- length(which(tweet_events))

    #get tweets to be retweeted (with replacement to allow for multiple retweets)
    retweets <- sample.int(fill_point-1, length(which(!tweet_events)), replace = TRUE, prob = probs)

    #for each unique retweet, add the number of times it was retweeted to its tally (this is, bizarrely, the fastest way of doing it)
      #it's not reflected in the profiler, but the factor to integer stuff here is one of the main bottlenecks (when factors are directly passed instead of integer it's way faster, but incorrect)
      #Rfast's version of as_integer has different behavior and does not work here
    table <- data.table::data.table(table(retweets))
    data.table::set(tweets, i = as.integer(table$retweets), j = "n_times", value = tweets$n_times[as.integer(table$retweets)]+table$N)

    #add 1 to the age of each old tweet
    data.table::set(tweets, i = 1:(fill_point-1), j = "age", value = tweets$age[1:(fill_point-1)]+1)

    #add new tweets to data table
    data.table::set(tweets, i = fill_point:(fill_point+num_events-1), j = "user", value = active_users[tweet_events])
    data.table::set(tweets, i = fill_point:(fill_point+num_events-1), j = "n_times", value = 1)
    data.table::set(tweets, i = fill_point:(fill_point+num_events-1), j = "age", value = 1)
    data.table::set(tweets, i = fill_point:(fill_point+num_events-1), j = "content", value = truncnorm::rtruncnorm(num_events, a = 0, b = Inf, mean = 1, sd = 1)^cont_bias)
    data.table::set(tweets, i = fill_point:(fill_point+num_events-1), j = "followers", value = users$followers[active_users[tweet_events]])

    if(diversity_TF){
      #add diversity to vector
      diversity[i] <- vegan::diversity(tweets$n_times, "simpson")
    }

    #remove variables
    rm(list = c("active_users", "fill_point", "probs", "tweet_events", "retweets"))
  }

  #remove excess allocation
  tweets <- subset(tweets, subset = (tweets$user > 0))

  if(sum_stats_TF){
    #calculate summary statistics
    sum_stats <- c(length(which(tweets$n_times == 1))/sum(tweets$n_times), #proportion of tweets that only appear once (among all events)
                   max(tweets$n_times)/sum(tweets$n_times), #proportion of the most retweeted tweet (among all events)
                   hillR::hill_taxa(tweets$n_times, q = 1), #hill number (shannon's diversity index)
                   hillR::hill_taxa(tweets$n_times, q = 2)) #hill number (simpson's diversity index)

    if(diversity_TF){
      #return summary statistics and diversity values
      return(list(sum_stats = sum_stats, diversity = diversity))
    }

    if(!diversity_TF){
      #return summary statistics and diversity values
      return(sum_stats)
    }
  }

  if(!sum_stats_TF){
    if(diversity_TF){
      #return data table of tweets and diversity values
      return(list(tweets = tweets, diversity = diversity))
    }

    if(!diversity_TF){
      #return data table of tweets
      return(tweets)
    }
  }
}
