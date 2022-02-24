adjust_ifr <- function(all_age_ifr_path, adj_file_name){
  
  df_ifr1 <- read.csv(file.path(output_dir, "ifr_prepped_input_data.csv"), as.is = TRUE) %>%
    mutate(study_id2 = paste0(location_id, "__", date, "__", bias_type))
  
  # the quantity seroprevalence / sero_sample_mean is the 
  # estimated proportion of true infections that are detected by the seroprevalence survey
  df_tmp_in <- read.csv(file.path(all_age_ifr_path, "sero_data.csv")) %>%
    mutate(adj = sero_sample_mean / seroprevalence ) %>%
    mutate(study_id2 = paste0(location_id, "__", sero_end_date, "__", bias_type))
  
  df_tmp <- df_tmp_in %>%
    mutate(study_id_tmp = paste0(location_id, "__", sero_end_date, "__", geo_accordance, "__", bias_type)) %>%
    select(study_id_tmp, location_id, is_outlier, sero_end_date, sero_sample_mean, seroprevalence, adj) %>%
    mutate(
      has_allage_seroprev = 1
    )
  
  df_tmp2 <- df_tmp %>%
    filter(duplicated(study_id_tmp))
  
  df_tmp3 <- df_tmp %>%
    filter(!duplicated(study_id_tmp)) %>%
    select(-sero_end_date, -location_id)
  
  df_ifr_adj <- df_ifr1 %>%
    mutate(sero_end_date = date) %>%
    mutate(study_id_tmp = paste0(location_id, "__", sero_end_date, "__", geo_accordance, "__", bias_type))
  
  df5 <- df_ifr_adj %>%
    left_join(df_tmp3, by = "study_id_tmp")
  
  #save file with waning immunity adjustments included
  write.csv(df5, paste0(output_dir, adj_file_name))
  
  #apply adjustments/outliering prior to modeling ifr and sero-age
  df5 <- as.data.table(df5)
  
  #filter according to all-age outlier decisions
  outliers <- df5[is.na(adj) & has_allage_seroprev==1] #n=140
  outliers2 <- df5[geo_accordance==0]
  df5 <- df5[!is.na(adj)]
  
  #remove geoaccordance == 0
  df5 <- df5[geo_accordance==1]
  
  #divide ifr and ifr se by 'adj' (seroprevalence adjustment factor)
  df5[, ifr_adj:=ifr/adj]
  df5[, ifr_se_adj:=ifr_se/adj]
  
  #multiply seroprev and seroprev se by 'adj'; se created in logit space, so first back-transform
  df5 <- as.data.frame(df5)
  df5[, c('prev_linear', 'se_linear')] <- mrbrt002::logit_to_linear(
    mean=array(df5$prev_logit),
    sd = array(df5$se_logit)
  )
  
  df5$prev_adj <- df5$prev_linear * df5$adj
  df5$se_adj <- df5$se_linear * df5$adj

  #retransform ifr and seroprev to logit space
  df5 <- as.data.frame(df5)
  df5[, c("logit_ifr_adj", "logit_ifr_se_adj")] <- mrbrt002::linear_to_logit(
    mean = array(df5$ifr_adj),
    sd = array(df5$ifr_se_adj)
  )
  df5[, c("prev_logit_adj", "se_logit_adj")] <- mrbrt002::linear_to_logit(
    mean = array(df5$prev_adj),
    sd = array(df5$se_adj)
  )
  
  #re-set names
  df5[, c('logit_ifr', 'logit_ifr_se')] <- NULL
  setnames(df5, c('logit_ifr_adj', 'logit_ifr_se_adj'), c('logit_ifr', 'logit_ifr_se'))
  df5 <- as.data.table(df5)
  
  #set up a constant se
  seroprev_logit_se_scaled_median <- median(df5$seroprev_logit_se_scaled)
  df5[, constant_median_se:=seroprev_logit_se_scaled_median]
  
  #re-save
  write.csv(df5, paste0(output_dir, 'ifr_prepped_input_data_outliered_adjusted.csv'), row.names=F)
  write.csv(outliers, paste0(output_dir, 'all_age_ifr_outliers.csv'), row.names=F)
  write.csv(outliers2, paste0(output_dir, 'geodiscordant_ifr_outliers.csv'), row.names=F)
  
  return(df5)
  
}