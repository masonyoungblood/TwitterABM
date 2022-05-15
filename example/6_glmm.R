# PREPARING DATA ----------------------------------------------------------

#load packages
library(lme4)
library(DHARMa)
library(glmmTMB)

#turn off scientific notation
options(scipen = 999)

#load user and tweet data data
load("data/promoters.RData")
load("data/promoter_tweets.RData")

#get rows of target tweets from within community, some of which are themselves quote tweets
original_tweets <- which(promoter_tweets$tweet_id %in% promoter_tweets$quote_tweet)

#trim so only includes original tweets that are not themselves quote tweets
original_tweets <- original_tweets[which(is.na(promoter_tweets$quote_tweet[original_tweets]))]

#get rows of quotes of original tweets from within community
quote_tweets <- which(promoter_tweets$quote_tweet %in% promoter_tweets$tweet_id[original_tweets])

#construct data frame for the quote glmm
quote_glmm_data <- data.frame(orig_quote = as.factor(c(rep(0, length(original_tweets)), rep(1, length(quote_tweets)))),
                              target = as.factor(c(promoter_tweets$tweet_id[original_tweets], promoter_tweets$quote_tweet[quote_tweets])),
                              user = as.factor(promoter_tweets$user_id[c(original_tweets, quote_tweets)]),
                              length = promoter_tweets$length[c(original_tweets, quote_tweets)],
                              negative = promoter_tweets$neg[c(original_tweets, quote_tweets)],
                              neutral = promoter_tweets$neu[c(original_tweets, quote_tweets)],
                              positive = promoter_tweets$pos[c(original_tweets, quote_tweets)],
                              compound = promoter_tweets$compound[c(original_tweets, quote_tweets)],
                              has_media = promoter_tweets$has_media[c(original_tweets, quote_tweets)])

#construct data frame for the main glmm
glmm_data <- data.frame(user = as.factor(promoter_tweets$user_id),
                        retweets = promoter_tweets$retweet_count,
                        positive = promoter_tweets$pos,
                        negative = promoter_tweets$neg,
                        neutral = promoter_tweets$neu,
                        compound = promoter_tweets$compound,
                        has_media = promoter_tweets$has_media,
                        length = promoter_tweets$length,
                        day = as.factor(as.Date(promoter_tweets$timestamp)),
                        hour = as.factor(lubridate::hour(promoter_tweets$timestamp)))

#binarize media
quote_glmm_data$has_media <- ifelse(quote_glmm_data$has_media == TRUE, 1, 0)
glmm_data$has_media <- ifelse(glmm_data$has_media == TRUE, 1, 0)

#add followers and verification status
glmm_data$followers <- promoters$followers_count[match(glmm_data$user, promoters$user_id)]
glmm_data$verified <- promoters$verified[match(glmm_data$user, promoters$user_id)]

#clean up verification data, including setting 1.9% of empty data to "false"
glmm_data$verified[which(glmm_data$verified == "")] <- "False"
glmm_data$verified[which(glmm_data$verified == "False")] <- 0
glmm_data$verified[which(glmm_data$verified == "True")] <- 1
glmm_data$verified <- as.numeric(glmm_data$verified)

#set 1.9% of missing followers data to 0
glmm_data$followers[which(is.na(glmm_data$followers))] <- 0

#remove unnecessary objects
rm(list = c("promoter_tweets", "promoters", "original_tweets", "quote_tweets"))

# MAIN MODEL --------------------------------------------------------------

#set random seed and get subset for model choice
set.seed(12345)
subset <- sample(nrow(glmm_data), nrow(glmm_data)*0.1)

#compare random effects models
user_model <- lme4::glmer(retweets ~ 1 + (1|user), data = glmm_data, family = poisson, subset = subset)
day_model <- lme4::glmer(retweets ~ 1 + (1|day), data = glmm_data, family = poisson, subset = subset)
hour_model <- lme4::glmer(retweets ~ 1 + (1|hour), data = glmm_data, family = poisson, subset = subset)

#save
save(user_model, file = "data/glmm_output/user_model.RData")
save(day_model, file = "data/glmm_output/day_model.RData")
save(hour_model, file = "data/glmm_output/hour_model.RData")

#compare icc
performance::icc(user_model) #0.610
performance::icc(day_model) #0.109
performance::icc(hour_model) #0.039

#compare model fit
AIC(user_model, day_model, hour_model)
lmtest::lrtest(user_model, day_model, hour_model)

#add length as control variable
length_model <- lme4::glmer(retweets ~ scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)

#save
save(length_model, file = "data/glmm_output/length_model.RData")

#compare model fit
AIC(user_model, length_model)
lmtest::lrtest(user_model, length_model)

#add media as control variable
media_model <- lme4::glmer(retweets ~ has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)

#save
save(media_model, file = "data/glmm_output/media_model.RData")

#compare model fit
AIC(length_model, media_model)
lmtest::lrtest(length_model, media_model)

#add followers and verification status
follow_model <- lme4::glmer(retweets ~ scale(followers) + has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)
verif_model <- lme4::glmer(retweets ~ verified + has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)

#save
save(follow_model, file = "data/glmm_output/follow_model.RData")
save(verif_model, file = "data/glmm_output/verif_model.RData")

#compare model fit
AIC(media_model, follow_model, verif_model)
lmtest::lrtest(media_model, follow_model)
lmtest::lrtest(media_model, verif_model)
lmtest::lrtest(follow_model, verif_model)

#add different forms of sentiment
comp_model <- lme4::glmer(retweets ~ scale(compound) + scale(followers) + has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)
neg_model <- lme4::glmer(retweets ~ scale(negative) + scale(followers) + has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)
posi_model <- lme4::glmer(retweets ~ scale(positive) + scale(followers) + has_media + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)

#save
save(comp_model, file = "data/glmm_output/comp_model.RData")
save(neg_model, file = "data/glmm_output/neg_model.RData")
save(posi_model, file = "data/glmm_output/posi_model.RData")

#compare model fit
AIC(follow_model, comp_model, neg_model, posi_model)
lmtest::lrtest(posi_model, neg_model)
lmtest::lrtest(comp_model, neg_model)
lmtest::lrtest(comp_model, posi_model)

#get AIC values for all models
AIC(user_model, day_model, hour_model, length_model, media_model, follow_model, verif_model, comp_model, neg_model, posi_model)

#best fitting has the compound score
summary(comp_model)

#re-run best model with the full dataset
full_comp_model <- lme4::glmer(retweets ~ scale(compound) + scale(followers) + has_media + scale(length) + (1|user), data = glmm_data, family = poisson)

#save
save(full_comp_model, file = "data/glmm_output/full_comp_model.RData")

#check for multicollinearity problems
car::vif(full_comp_model)

#get incidence rate ratios
exp(fixef(full_comp_model))
exp(confint(full_comp_model, method = "Wald"))

#get pseudo R squared values for the model
MuMIn::r.squaredGLMM(full_comp_model)

#check fit of best model with dharma (run on the 10%)
dharma_results <- simulateResiduals(fittedModel = comp_model, plot = F)
testDispersion(dharma_results)
testUniformity(dharma_results)
testZeroInflation(dharma_results)
testOutliers(dharma_results)

#for a partial specification curve analysis, run the best model with different combinations of predictors to ensure the results hold up
check_comp_1 <- glm(retweets ~ scale(compound), data = glmm_data, family = poisson, subset = subset)
check_comp_2 <- lme4::glmer(retweets ~ scale(compound) + (1|user), data = glmm_data, family = poisson, subset = subset)
check_comp_3 <- lme4::glmer(retweets ~ scale(compound) + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)
check_follow_1 <- glm(retweets ~ scale(followers), data = glmm_data, family = poisson, subset = subset) #no convergence
check_follow_2 <- lme4::glmer(retweets ~ scale(followers) + (1|user), data = glmm_data, family = poisson, subset = subset)
check_follow_3 <- lme4::glmer(retweets ~ scale(followers) + scale(length) + (1|user), data = glmm_data, family = poisson, subset = subset)

#save
save(check_comp_1, file = "data/glmm_output/check_comp_1.RData")
save(check_comp_2, file = "data/glmm_output/check_comp_2.RData")
save(check_comp_3, file = "data/glmm_output/check_comp_3.RData")
save(check_follow_1, file = "data/glmm_output/check_follow_1.RData")
save(check_follow_2, file = "data/glmm_output/check_follow_2.RData")
save(check_follow_3, file = "data/glmm_output/check_follow_3.RData")

#IRRs for checks (all three neg models, full, all three follow models, full)
check_effects <- c(exp(as.numeric(check_comp_1$coefficients[2])), exp(as.numeric(confint(check_comp_1, method = "Wald")[2,])),
                   exp(as.numeric(fixef(check_comp_2)[2])), exp(as.numeric(confint(check_comp_2, method = "Wald")[3,])),
                   exp(as.numeric(fixef(check_comp_3)[2])), exp(as.numeric(confint(check_comp_3, method = "Wald")[3,])),
                   exp(as.numeric(fixef(full_comp_model)[2])), exp(as.numeric(confint(full_comp_model, method = "Wald")[3,])),
                   NA, c(NA, NA), #no convergence
                   exp(as.numeric(fixef(check_follow_2)[2])), exp(as.numeric(confint(check_follow_2, method = "Wald")[3,])),
                   exp(as.numeric(fixef(check_follow_3)[2])), exp(as.numeric(confint(check_follow_3, method = "Wald")[3,])),
                   exp(as.numeric(fixef(full_comp_model)[3])), exp(as.numeric(confint(full_comp_model, method = "Wald")[4,])))
check_effects <- as.data.frame(matrix(check_effects, ncol = 3, byrow = TRUE))
colnames(check_effects) <- c("estimate", "lower", "upper")
check_effects$model <- factor(c(1:8))

# QUOTE MODEL -------------------------------------------------------------

#set random seed
set.seed(12345)

#baseline model
target_quote_model <- lme4::glmer(orig_quote ~ (1|target), data = quote_glmm_data, family = binomial)

#save
save(target_quote_model, file = "data/glmm_output/target_quote_model.RData")

#check icc of target
performance::icc(target_quote_model) #0.146

#add length as control variable
length_quote_model <- lme4::glmer(orig_quote ~ scale(length) + (1|target), data = quote_glmm_data, family = binomial)

#save
save(length_quote_model, file = "data/glmm_output/length_quote_model.RData")

#compare model fit
AIC(target_quote_model, length_quote_model)
lmtest::lrtest(target_quote_model, length_quote_model)

#add media as control variable
media_quote_model <- lme4::glmer(orig_quote ~ has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)

#save
save(media_quote_model, file = "data/glmm_output/media_quote_model.RData")

#compare model fit
AIC(length_quote_model, media_quote_model)
lmtest::lrtest(length_quote_model, media_quote_model)

#assess which sentiment measures best improve the fit of the model
comp_quote_model <- lme4::glmer(orig_quote ~ scale(compound) + has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)
abs_comp_quote_model <- lme4::glmer(orig_quote ~ scale(abs(compound)) + has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)
neg_quote_model <- lme4::glmer(orig_quote ~ scale(negative) + has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)
neu_quote_model <- lme4::glmer(orig_quote ~ scale(neutral) + has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)
posi_quote_model <- lme4::glmer(orig_quote ~ scale(positive) + has_media + scale(length) + (1|target), data = quote_glmm_data, family = binomial)

#save
save(comp_quote_model, file = "data/glmm_output/comp_quote_model.RData")
save(abs_comp_quote_model, file = "data/glmm_output/abs_comp_quote_model.RData")
save(neg_quote_model, file = "data/glmm_output/neg_quote_model.RData")
save(neu_quote_model, file = "data/glmm_output/neu_quote_model.RData")
save(posi_quote_model, file = "data/glmm_output/posi_quote_model.RData")

#compare model fit
AIC(media_quote_model, comp_quote_model, abs_comp_quote_model, neg_quote_model, neu_quote_model, posi_quote_model)
lmtest::lrtest(media_quote_model, comp_quote_model)

#get AIC values for all models
AIC(target_quote_model, length_quote_model, media_quote_model, comp_quote_model, abs_comp_quote_model, neg_quote_model, posi_quote_model, neu_quote_model)

#compound model is best fitting
summary(comp_quote_model)

#check for multicollinearity problems
car::vif(comp_quote_model)

#get odds ratios
exp(fixef(comp_quote_model))
exp(confint(comp_quote_model, method = "Wald"))

#get pseudo R squared values for the model
MuMIn::r.squaredGLMM(comp_quote_model)

#check fit of best model with dharma
dharma_quote_results <- simulateResiduals(fittedModel = comp_quote_model, plot = F)
testUniformity(dharma_quote_results)
testOutliers(dharma_quote_results)
