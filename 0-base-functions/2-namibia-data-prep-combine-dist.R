# combine all the dist variables to avoid colinearity
comb_distvar <- function(df_all){
  # first, make them mutually exclusive, and combine arm together
  df_all <- df_all %>% mutate(pre_int_RO_3000 = pre_int_RO_3000 - pre_int_RO_2000,
                              pre_int_RO_2000 = pre_int_RO_2000 - pre_int_RO_1000,
                              pre_int_RO_1000 = pre_int_RO_1000 - pre_int_RO_500,
                              
                              pre_int_RV_3000 = pre_int_RV_3000 - pre_int_RV_2000,
                              pre_int_RV_2000 = pre_int_RV_2000 - pre_int_RV_1000,
                              pre_int_RV_1000 = pre_int_RV_1000 - pre_int_RV_500,
                              
                              pre_int_TO_3000 = pre_int_TO_3000 - pre_int_TO_2000,
                              pre_int_TO_2000 = pre_int_TO_2000 - pre_int_TO_1000,
                              pre_int_TO_1000 = pre_int_TO_1000 - pre_int_TO_500,
                              
                              pre_int_TV_3000 = pre_int_TV_3000 - pre_int_TV_2000,
                              pre_int_TV_2000 = pre_int_TV_2000 - pre_int_TV_1000,
                              pre_int_TV_1000 = pre_int_TV_1000 - pre_int_TV_500,
                              
                              pre_int_3000 = pre_int_TV_3000 + pre_int_TO_3000 +
                                             pre_int_RV_3000 + pre_int_RO_3000,
                              pre_int_2000 = pre_int_TV_2000 + pre_int_TO_2000 +
                                             pre_int_RV_2000 + pre_int_RO_2000,
                              pre_int_1000 = pre_int_TV_1000 + pre_int_TO_1000 +
                                             pre_int_RV_1000 + pre_int_RO_1000,
                              pre_int_500 = pre_int_TV_500 + pre_int_TO_500 +
                                            pre_int_RV_500 + pre_int_RO_500,
                              
                              #
                              concur_int_RO_3000 = concur_int_RO_3000 - concur_int_RO_2000,
                              concur_int_RO_2000 = concur_int_RO_2000 - concur_int_RO_1000,
                              concur_int_RO_1000 = concur_int_RO_1000 - concur_int_RO_500,
                              
                              concur_int_RV_3000 = concur_int_RV_3000 - concur_int_RV_2000,
                              concur_int_RV_2000 = concur_int_RV_2000 - concur_int_RV_1000,
                              concur_int_RV_1000 = concur_int_RV_1000 - concur_int_RV_500,
                              
                              concur_int_TO_3000 = concur_int_TO_3000 - concur_int_TO_2000,
                              concur_int_TO_2000 = concur_int_TO_2000 - concur_int_TO_1000,
                              concur_int_TO_1000 = concur_int_TO_1000 - concur_int_TO_500,
                              
                              concur_int_TV_3000 = concur_int_TV_3000 - concur_int_TV_2000,
                              concur_int_TV_2000 = concur_int_TV_2000 - concur_int_TV_1000,
                              concur_int_TV_1000 = concur_int_TV_1000 - concur_int_TV_500,
                              
                              concur_int_3000 = concur_int_TV_3000 + concur_int_TO_3000 +
                                                concur_int_RV_3000 + concur_int_RO_3000,
                              concur_int_2000 = concur_int_TV_2000 + concur_int_TO_2000 +
                                                concur_int_RV_2000 + concur_int_RO_2000,
                              concur_int_1000 = concur_int_TV_1000 + concur_int_TO_1000 +
                                                concur_int_RV_1000 + concur_int_RO_1000,
                              concur_int_500 = concur_int_TV_500 + concur_int_TO_500 +
                                               concur_int_RV_500 + concur_int_RO_500,
                              
                              # 
                              pop_3km = pop_3km - pop_2km,
                              pop_2km = pop_2km - pop_1km,
                              pop_1km = pop_1km - pop_500m,
                              
                              # 
                              pre_index_3km = pre_index_3km - pre_index_2km,
                              pre_index_2km = pre_index_2km - pre_index_1km,
                              pre_index_1km = pre_index_1km - pre_index_500m,
                              
                              n_int_pre_recei = n_int_pre_recei_RO + 
                                                n_int_pre_recei_RV +
                                                n_int_pre_recei_TO +
                                                n_int_pre_recei_TV)
  
  # second, add indicator % of discordant arm, for each variable, each dist
  df_all_RO <- df_all %>% filter(intarm == 'RO')
  df_all_RV <- df_all %>% filter(intarm == 'RV')
  df_all_TO <- df_all %>% filter(intarm == 'TO')
  df_all_TV <- df_all %>% filter(intarm == 'TV')
  
  assert_that(nrow(df_all) == nrow(df_all_RO) + 
                nrow(df_all_RV) +
                nrow(df_all_TO) + 
                nrow(df_all_TV))
  
  df_all_RO <- df_all_RO %>% 
    mutate(pre_int_propdiffarm_500 = ifelse(pre_int_500 == 0, 
                                            0, 1-pre_int_RO_500/pre_int_500),
           
           pre_int_propdiffarm_1000 = ifelse(pre_int_1000 == 0, 
                                             0, 1-pre_int_RO_1000/pre_int_1000),
           
           pre_int_propdiffarm_2000 = ifelse(pre_int_2000 == 0, 
                                             0,1-pre_int_RO_2000/pre_int_2000),
           
           pre_int_propdiffarm_3000 = ifelse(pre_int_3000 == 0, 
                                             0, 1-pre_int_RO_3000/pre_int_3000),
           
           # 
           concur_int_propdiffarm_500 = ifelse(concur_int_500 == 0, 
                                               0, 1-concur_int_RO_500/concur_int_500),
           
           concur_int_propdiffarm_1000 = ifelse(concur_int_1000 == 0, 
                                                0, 1-concur_int_RO_1000/concur_int_1000),
           
           concur_int_propdiffarm_2000 = ifelse(concur_int_2000 == 0, 
                                                0,1-concur_int_RO_2000/concur_int_2000),
           
           concur_int_propdiffarm_3000 = ifelse(concur_int_3000 == 0, 
                                                0, 1-concur_int_RO_3000/concur_int_3000)
    )
  
  df_all_RV <- df_all_RV %>% 
    mutate(pre_int_propdiffarm_500 = ifelse(pre_int_500 == 0, 
                                            0, 1-pre_int_RV_500/pre_int_500),
           
           pre_int_propdiffarm_1000 = ifelse(pre_int_1000 == 0, 
                                             0, 1-pre_int_RV_1000/pre_int_1000),
           
           pre_int_propdiffarm_2000 = ifelse(pre_int_2000 == 0, 
                                             0,1-pre_int_RV_2000/pre_int_2000),
           
           pre_int_propdiffarm_3000 = ifelse(pre_int_3000 == 0, 
                                             0, 1-pre_int_RV_3000/pre_int_3000),
           
           # 
           concur_int_propdiffarm_500 = ifelse(concur_int_500 == 0, 
                                               0, 1-concur_int_RV_500/concur_int_500),
           
           concur_int_propdiffarm_1000 = ifelse(concur_int_1000 == 0, 
                                                0, 1-concur_int_RV_1000/concur_int_1000),
           
           concur_int_propdiffarm_2000 = ifelse(concur_int_2000 == 0, 
                                                0,1-concur_int_RV_2000/concur_int_2000),
           
           concur_int_propdiffarm_3000 = ifelse(concur_int_3000 == 0, 
                                                0, 1-concur_int_RV_3000/concur_int_3000)
    )
  
  df_all_TO <- df_all_TO %>% 
    mutate(pre_int_propdiffarm_500 = ifelse(pre_int_500 == 0, 
                                            0, 1-pre_int_TO_500/pre_int_500),
           
           pre_int_propdiffarm_1000 = ifelse(pre_int_1000 == 0, 
                                             0, 1-pre_int_TO_1000/pre_int_1000),
           
           pre_int_propdiffarm_2000 = ifelse(pre_int_2000 == 0, 
                                             0,1-pre_int_TO_2000/pre_int_2000),
           
           pre_int_propdiffarm_3000 = ifelse(pre_int_3000 == 0, 
                                             0, 1-pre_int_TO_3000/pre_int_3000),
           
           # 
           concur_int_propdiffarm_500 = ifelse(concur_int_500 == 0, 
                                               0, 1-concur_int_TO_500/concur_int_500),
           
           concur_int_propdiffarm_1000 = ifelse(concur_int_1000 == 0, 
                                                0, 1-concur_int_TO_1000/concur_int_1000),
           
           concur_int_propdiffarm_2000 = ifelse(concur_int_2000 == 0, 
                                                0,1-concur_int_TO_2000/concur_int_2000),
           
           concur_int_propdiffarm_3000 = ifelse(concur_int_3000 == 0, 
                                                0, 1-concur_int_TO_3000/concur_int_3000)
    )
  
  df_all_TV <- df_all_TV %>% 
    mutate(pre_int_propdiffarm_500 = ifelse(pre_int_500 == 0, 
                                            0, 1-pre_int_TV_500/pre_int_500),
           
           pre_int_propdiffarm_1000 = ifelse(pre_int_1000 == 0, 
                                             0, 1-pre_int_TV_1000/pre_int_1000),
           
           pre_int_propdiffarm_2000 = ifelse(pre_int_2000 == 0, 
                                             0,1-pre_int_TV_2000/pre_int_2000),
           
           pre_int_propdiffarm_3000 = ifelse(pre_int_3000 == 0, 
                                             0, 1-pre_int_TV_3000/pre_int_3000),
           
           # 
           concur_int_propdiffarm_500 = ifelse(concur_int_500 == 0, 
                                               0, 1-concur_int_TV_500/concur_int_500),
           
           concur_int_propdiffarm_1000 = ifelse(concur_int_1000 == 0, 
                                                0, 1-concur_int_TV_1000/concur_int_1000),
           
           concur_int_propdiffarm_2000 = ifelse(concur_int_2000 == 0, 
                                                0,1-concur_int_TV_2000/concur_int_2000),
           
           concur_int_propdiffarm_3000 = ifelse(concur_int_3000 == 0, 
                                                0, 1-concur_int_TV_3000/concur_int_3000)
    )
  
  df_all <- bind_rows(df_all_RO, df_all_RV, df_all_TO, df_all_TV)
  # drop the original dist variables
  dist_vars <- c("pre_int_RO_500", "pre_int_RV_500", "pre_int_TO_500", "pre_int_TV_500",
                 "pre_int_RO_1000", "pre_int_RV_1000", "pre_int_TO_1000", "pre_int_TV_1000",
                 "pre_int_RO_2000", "pre_int_RV_2000", "pre_int_TO_2000", "pre_int_TV_2000",
                 "pre_int_RO_3000", "pre_int_RV_3000", "pre_int_TO_3000", "pre_int_TV_3000",
                 "concur_int_RO_500", "concur_int_RV_500", "concur_int_TO_500", "concur_int_TV_500",
                 "concur_int_RO_1000", "concur_int_RV_1000", "concur_int_TO_1000", "concur_int_TV_1000",
                 "concur_int_RO_2000", "concur_int_RV_2000", "concur_int_TO_2000", "concur_int_TV_2000",
                 "concur_int_RO_3000", "concur_int_RV_3000", "concur_int_TO_3000", "concur_int_TV_3000",
                 "n_int_pre_recei_RO", "n_int_pre_recei_RV", "n_int_pre_recei_TO", "n_int_pre_recei_TV")
  
  df_all <- df_all[,setdiff(names(df_all), dist_vars)]
  
  return(df_all)
}


