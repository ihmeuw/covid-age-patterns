predict_mort_age <- function(
  mod_global, cascade_fit, model_data, df_locmeta, global_mort_preds_path, mr_agepattern_preds_byloc_5yr, 
  mr_agepattern_preds_byloc_1yr, mr_agepattern_preds_global, 
  mr_agepattern_preds_sr, mr_agepattern_preds_loc
) {
  
  ## Global prediction ----------------------------------------------------------
  
  df_pred3 <- expand.grid(
    age_start = seq(0, 100, by = 1) ) %>%
    mutate(
      age_end = age_start,
      data_id = 1:nrow(.)
    )
  
  dat_pred3 <- MRData()
  
  dat_pred3$load_df(
    data = df_pred3,
    col_covs=list('age_start', 'age_end', 'data_id')
  )
  
  df_pred3$pred3 <- mod_global$predict(dat_pred3, sort_by_data_id = "data_id")
  
  write.csv(df_pred3, global_mort_preds_path, row.names = FALSE)
  
  # -- predict for super-regions
  df_pred_sr <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100,
    super_region_name = unique(model_data$super_region_name) ) %>%
    mutate(
      age_end = age_start,
      location_id = NA,
      age_mid = (age_start + age_end) / 2
    )
  
  cascade_preds_sr <- predict_spline_cascade(fit = cascade_fit, newdata = df_pred_sr)
  
  # -- predict for locations
  df_pred_loc <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100,
    location_id = unique(model_data$location_id) ) %>%
    mutate(
      age_end = age_start,
      age_mid = (age_start + age_end) / 2 ) %>%
    left_join(model_data %>% select(super_region_name, location_id, location) %>% filter(!duplicated(.)))
  
  cascade_preds_loc <- predict_spline_cascade(fit = cascade_fit, newdata = df_pred_loc)
  
  
  ####
  preds_global <- df_pred3
  preds_sr <- cascade_preds_sr
  preds_loc <- cascade_preds_loc
  
  
  df_child_ids <- data.frame(child_id = unique(df_locmeta$location_id))
  df_child_ids$location_name <- sapply(df_child_ids$child_id, function(x) {
    cid_srname <- as.character(df_locmeta[df_locmeta$location_id == x, "location_name"])
  })
  
  
  df_child_ids$mr_id <- sapply(1:nrow(df_child_ids), function(i) {
    dev <- FALSE
    if (dev) {
      i <- 43
    }
    
    cid <- df_child_ids[i, "child_id"]
    df_cid <- df_locmeta[df_locmeta$location_id == cid, ]
    cid_srname <- as.character(df_locmeta[df_locmeta$location_id == cid, "super_region_name"])
    
    parent_ids <- rev(strsplit(as.character(df_cid$path_to_top_parent), split = ",")[[1]])
    
    if (any(parent_ids %in% preds_loc$location_id)) {
      matcher <- parent_ids[min(which(parent_ids %in% preds_loc$location_id))]
      out <- paste0("location_id__", matcher)
    } else if (cid_srname %in% preds_sr$super_region_name) {
      out <- paste0("super_region_name__", cid_srname)
    } else {
      out <- "stage1__stage1"
    }
    return(out)
  })
  
  
  dat_mrpreds <- lapply(1:nrow(df_child_ids), function(i) {
    dev <- FALSE
    if (dev) {
      i <- 12
    }
    
    mr_id_tmp <- df_child_ids[i, "mr_id"]
    
    if (grepl("location_id__", mr_id_tmp)) {
      df_out <- preds_loc %>%
        filter(cascade_prediction_id == mr_id_tmp) %>%
        mutate(pred_mr = pred)
    } else if (grepl("super_region_name__", mr_id_tmp)) {
      df_out <- preds_sr %>%
        filter(cascade_prediction_id == mr_id_tmp) %>%
        mutate(pred_mr = pred)
    } else {
      df_out <- df_pred3
      df_out$pred_mr <- df_out$pred3
    }
    
    
    df_out2 <- df_out %>%
      mutate(
        age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE),
        pred_mr_ratio = inv_logit(pred_mr) / inv_logit(min(pred_mr)) )
    
    df_out3 <- df_out2 %>%
      filter(!is.na(age_group)) %>%
      group_by(age_group) %>%
      summarize(MR = mean(pred_mr_ratio)) %>%
      mutate(
        MRprob = MR / sum(MR),
        location_id = df_child_ids[i, "child_id"],
        location_name = df_child_ids[i, "location_name"],
        mr_id = mr_id_tmp,
        age_group = as.character(age_group)) %>%
      as.data.frame(.)
    
    df_out3$age_start <- as.integer(gsub("\\[", "", sapply(df_out3$age_group, function(x) strsplit(x, split = ",")[[1]][1])))
    df_out3$age_end <- as.integer(as.numeric(gsub("\\)", "", sapply(df_out3$age_group, function(x) strsplit(x, split = ",")[[1]][2]))) - 1)
    
    out <- list(
      one_year = df_out2 %>% mutate(
        location_id = df_child_ids[i, "child_id"],
        location_name = df_child_ids[i, "location_name"],
        mr_id = mr_id_tmp),
      five_year = df_out3
    )
    return(out)
  })
  
  df_mrpreds_5yr <- do.call("rbind", lapply(dat_mrpreds, function(x) x[["five_year"]]))
  df_mrpreds_1yr <- bind_rows(lapply(dat_mrpreds, function(x) x[["one_year"]]))
  
  
  write.csv(df_mrpreds_5yr, mr_agepattern_preds_byloc_5yr, row.names = FALSE)
  write.csv(df_mrpreds_1yr, mr_agepattern_preds_byloc_1yr, row.names = FALSE)
  
  saveRDS(preds_global, mr_agepattern_preds_global)
  saveRDS(preds_sr, mr_agepattern_preds_sr)
  saveRDS(preds_loc, mr_agepattern_preds_loc)
  
  return(list(preds_global, df_mrpreds_5yr, df_mrpreds_1yr, cascade_preds_sr))
}