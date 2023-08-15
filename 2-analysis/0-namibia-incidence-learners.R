# define learners for incidence analysis 

lrn_mean <- Lrnr_mean$new()
lrn_glm <- Lrnr_glm$new()
lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
lrn_elnet <- Lrnr_glmnet$new(alpha = 0.5)
lrn_xgboost <- Lrnr_xgboost$new(eval_metric = "logloss")


interactions <- list(c("inter", "A"))
# main terms as well as the interactions above will be included
lrn_interaction <- make_learner(Lrnr_define_interactions, interactions)
# use a pipeline to combine the interaction learn with other learners 
lrn_glm_interaction <- make_learner(Pipeline, lrn_interaction, lrn_glm)
lrn_lasso_interaction <- make_learner(Pipeline, lrn_interaction, lrn_lasso)
lrn_elnet_interaction <- make_learner(Pipeline, lrn_interaction, lrn_elnet)
lrn_xgboost_interaction <- make_learner(Pipeline, lrn_interaction, lrn_xgboost)


SL_lib_Q_int <- make_learner(Stack,
                         lrn_mean,
                         lrn_glm,
                         lrn_lasso,
                         lrn_elnet,
                         lrn_xgboost,
                         lrn_glm_interaction,
                         lrn_lasso_interaction,
                         lrn_elnet_interaction,
                         lrn_xgboost_interaction)

SL_lib_g <- make_learner(Stack,
                         lrn_mean,
                         lrn_glm,
                         lrn_lasso,
                         lrn_elnet)


SL_lib_Q <- make_learner(Stack,
                         lrn_mean,
                         lrn_glm,
                         lrn_lasso,
                         lrn_elnet,
                         lrn_xgboost)


SL_lib_simple <- make_learner(Stack,
                              lrn_mean,
                              lrn_glm,
                              lrn_lasso)
