fit_ifr_age <- function(df5, prev_best){
  
  optimize_knots <- T
  predict_ifr_from_ensemble <- F
  
  #create MRData object
  dat1 <- MRData()
  dat1$load_df(
    data = df5,  
    col_obs = "logit_ifr", 
    col_obs_se = 'constant_median_se',
    col_covs = list("age_start", "age_end", "bias"), 
    col_study_id = "study_id")
  
  #knot and spline info
  n_knots <- 4L
  spline_cov_model1 <- LinearCovModel(
    alt_cov = list("age_start", "age_end"),
    use_re = FALSE,
    use_spline = TRUE,
    spline_knots = array(seq(0, 1, length.out = n_knots + 2)),
    spline_degree = 3L,
    spline_knots_type = "domain",
    spline_r_linear = TRUE,
    spline_l_linear = FALSE
  )
  
  #Create covs
  intercept_and_bias_cov <- list(
    LinearCovModel("intercept", use_re = TRUE),
    LinearCovModel("bias", use_re = FALSE)
  )
  
  # Find optimal knots
  ifr_best_knots_file <- file.path(output_dir, "best_ifr_knots_pre_covselection.RDS")
  if (optimize_knots | predict_ifr_from_ensemble) {
    kd <- c(0.05, 0.6) # knot domain
    is <- c(0.08, 1) # interval sizes
    knot_bounds1 <- do.call("rbind", lapply(1:n_knots, function(x) kd))
    np <- import("numpy")
    np$random$seed(123L)
    
    knots_samples1 <- utils$sample_knots(
      num_intervals = n_knots + 1L, 
      knot_bounds = knot_bounds1,
      interval_sizes = do.call("rbind", lapply(1:(n_knots+1), function(x) is)),
      num_samples = 100L
    )
    
    mod1 <- MRBeRT(
      data = dat1,
      ensemble_cov_model = spline_cov_model1,
      cov_models = intercept_and_bias_cov,
      ensemble_knots = knots_samples1
    )
    
    mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
    
    py_save_object(
      object = mod1,
      filename = file.path(output_dir, "ifr_mod1.pkl"),
      pickle = "dill"
    )
    
    #save best selected knots 
    best_knots1 <- mod1$ensemble_knots[which.max(mod1$weights), ]
    saveRDS(best_knots1, ifr_best_knots_file)
  } else {
    best_knots1 <- readRDS(ifr_best_knots_file)
    mod1 <- py_load_object(
      filename = file.path(output_dir, "ifr_mod1.pkl"),
      pickle = "dill"
    )
  }
  
  if (predict_ifr_from_ensemble) {
    df_pred1$pred <- mod1$predict(data = dat_pred1)
    cat("Predicting from ensemble...\n")
  } else {
    if (optimize_knots) {
      knots_tmp <- array(best_knots1)
    } else {
      knots_tmp <- array(c(0, 0.12, 0.25, 0.55, 1))
    }
    mod2 <- MRBRT(
      data = dat1,
      cov_models = list(
        LinearCovModel("intercept", use_re = TRUE),
        LinearCovModel("bias", use_re = FALSE),
        LinearCovModel(
          alt_cov = list("age_start", "age_end"),
          use_re = FALSE,
          use_spline = TRUE,
          spline_knots = knots_tmp,
          spline_degree = 3L,
          spline_knots_type = "domain",
          spline_r_linear = TRUE,
          spline_l_linear = FALSE
        )
      )
      
    )
    mod2$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  }
  
  py_save_object(
    object = mod2,
    filename = file.path(output_dir, "ifr_mod2_global.pkl"),
    pickle = "dill"
  )
  
  #want to keep bias covariate constant through levels of the cascade; manually adjust model object here
  mod2$cov_models[[2]]$prior_beta_uniform <- matrix(c(mod2$fe_soln$bias, mod2$fe_soln$bias))
  bias_estimate <- mod2$fe_soln$bias
  
  #run cascade model from mod2 global fit, specifying high-income; high-income high nursing home; vs. other
  nursing_home_locs <- c(76, 101, 84, 529, 530, 546, 552, 556, 558, 560, 562, 570, 555)
  write.csv(data.table(location_id=nursing_home_locs), paste0(output_dir, '/nursing_home_location_ids.csv'), row.names=F)
  df5[location_id %in% nursing_home_locs, custom_hierarchy:='nursing_homes']
  df5[super_region_name=='High-income' & !location_id %in% nursing_home_locs, custom_hierarchy:='high_income']
  df5[super_region_name!='High-income' & !location_id %in% nursing_home_locs, custom_hierarchy:='all_other']
  
  #fit cascade spline
  model_label_tmp <- "cascade_ifr"
  thetas <- c(3)
  cascade_fit1 <- run_spline_cascade(
    stage1_model_object = mod2, #global model
    df = as.data.frame(df5), 
    col_obs = "logit_ifr",
    col_obs_se = 'constant_median_se',
    col_study_id = "study_id",
    stage_id_vars = c("custom_hierarchy"),
    thetas = thetas,
    gaussian_prior = TRUE,
    output_dir = output_dir,
    model_label = model_label_tmp, 
    overwrite_previous = TRUE
  )
  
  # (2) Predict -------------------------------------------------------------------------------------------------------------------
  
  ## Global prediction
  df_pred_glob <- data.frame(age_start = seq(0, 100, by = 1)) %>%
    mutate(
      age_end = age_start,
      bias = 0
    )
  dat_pred_glob <- MRData()
  dat_pred_glob$load_df(
    data = df_pred_glob,
    col_covs=list('age_start', 'age_end', 'bias')
  )
  
  df_pred_glob$pred_glob <- mod2$predict(dat_pred_glob)
  
  # -- predict for high-income vs. other
  df_pred_sr <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100,
    custom_hierarchy = unique(df5$custom_hierarchy)) %>%
    mutate(
      age_end = age_start,
      location_id = NA,
      region_name = NA,
      age_mid = (age_start + age_end) / 2,
      bias = 0
    )
  
  cascade_preds_ch <- predict_spline_cascade(fit = cascade_fit1, newdata = as.data.frame(df_pred_sr))
  
  #save
  saveRDS(df_pred_glob, paste0(output_dir, 'ifr_preds_1yr_global.rds'))
  saveRDS(cascade_preds_ch, paste0(output_dir, 'ifr_preds_1yr_3set.rds'))
  
  # (3) 5 yr predictions for GBD use ---------------------------------------------------------------------------------------------------------------
  
  prep_ifr_5yr <- function(df_pred1) {
    # Prep five year age groups
    df_ifr_new <- df_pred1 %>%
      mutate(
        age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE) ) %>%
      filter(!is.na(age_group)) %>%
      group_by(age_group, custom_hierarchy) %>%
      summarize(ifr_new = mean(pred_invlogit)) %>%
      mutate(ifr_new_logit = logit(ifr_new)) %>%
      as.data.frame(.)
    df_ifr_new$age_start <- as.integer(gsub("\\[", "", sapply(as.character(df_ifr_new$age_group), function(x) strsplit(x, split = ",")[[1]][1])))
    df_ifr_new$age_end <- as.integer(as.numeric(gsub("\\)", "", sapply(as.character(df_ifr_new$age_group), function(x) strsplit(x, split = ",")[[1]][2]))) - 1)
    df_ifr_new_out <- df_ifr_new %>%
      select(age_group_start = age_start, age_group_end = age_end, logit_ifr = ifr_new_logit, ifr = ifr_new, custom_hierarchy=custom_hierarchy)
    
    return(df_ifr_new_out)
  }
  
  df_pred1 <- readRDS(paste0(output_dir, 'ifr_preds_1yr_3set.rds'))
  df_pred1$pred_invlogit <- inv_logit(df_pred1$pred)
  preds_5yr <- prep_ifr_5yr(df_pred1)
  
  #save 5 yr predictions
  write.csv(preds_5yr, paste0(output_dir, 'ifr_preds_5yr_3set.csv'), row.names = FALSE)
  
  #GLOBAL
  prep_ifr_5yr_global <- function(df_pred1) {
    # Prep five year age groups
    df_ifr_new <- df_pred1 %>%
      mutate(
        age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE) ) %>%
      filter(!is.na(age_group)) %>%
      group_by(age_group) %>%
      summarize(ifr_new = mean(pred_invlogit)) %>%
      mutate(ifr_new_logit = logit(ifr_new)) %>%
      as.data.frame(.)
    df_ifr_new$age_start <- as.integer(gsub("\\[", "", sapply(as.character(df_ifr_new$age_group), function(x) strsplit(x, split = ",")[[1]][1])))
    df_ifr_new$age_end <- as.integer(as.numeric(gsub("\\)", "", sapply(as.character(df_ifr_new$age_group), function(x) strsplit(x, split = ",")[[1]][2]))) - 1)
    df_ifr_new_out <- df_ifr_new %>%
      select(age_group_start = age_start, age_group_end = age_end, logit_ifr = ifr_new_logit, ifr = ifr_new)
    
    return(df_ifr_new_out)
  }
  
  df_pred_glob <- readRDS(paste0(output_dir, 'ifr_preds_1yr_global.rds'))
  df_pred_glob$pred_invlogit <- inv_logit(df_pred_glob$pred_glob)
  preds_5yr_global <- prep_ifr_5yr_global(df_pred_glob)
  
  #save 5 yr predictions
  write.csv(preds_5yr_global, paste0(output_dir, 'ifr_preds_5yr_global.csv'), row.names = FALSE)
  
  # (4) CREATE PREDICTION FILE FOR EACH COVID PRODUCTION LOCATION ---------------------------------------------------------------------------------------------
  
  preds_5yr <- as.data.table(preds_5yr)
  hi <- preds_5yr[custom_hierarchy=='high_income']
  nh <- preds_5yr[custom_hierarchy=='nursing_homes']
  oth <- preds_5yr[custom_hierarchy=='all_other']
  
  #add hierarchy info to loc meta
  df_locmeta <- as.data.table(df_locmeta)
  df_locmeta[location_id %in% nursing_home_locs, custom_hierarchy:='nursing_homes']
  df_locmeta[super_region_name=='High-income' & !location_id %in% nursing_home_locs, custom_hierarchy:='high_income']
  df_locmeta[super_region_name!='High-income' & !location_id %in% nursing_home_locs, custom_hierarchy:='all_other']
  
  all_preds_5yr <- data.table()
  for (l in unique(df_locmeta[level>=3]$location_id)){
    print(l)
    ch <- df_locmeta[location_id==l]$custom_hierarchy
    preds <- preds_5yr[custom_hierarchy==ch]
    preds[, location_id:=paste0(l)]
    all_preds_5yr <- rbind(all_preds_5yr, preds)
  }
  
  preds_5yr_global$location_id <- 1
  preds_5yr_global$custom_hierarchy <- 'global'
  
  all_preds_5yr <- rbind(all_preds_5yr, preds_5yr_global)
  write.csv(all_preds_5yr, paste0(output_dir, 'ifr_preds_5yr_byloc.csv'), row.names=F)
  
  # (4) Plot ---------------------------------------------------------------------------------------------------------------------------------------
  
  model_data <- copy(df5)
  model_data <- model_data %>%
    mutate(age_mid = (age_start + age_end) / 2) %>% 
    as.data.table(model_data)
  
  df_pred_glob$result <- 'Global fit'
  cascade_preds_ch$result <- 'Super region fit'
  cascade_preds_ch <- as.data.table(cascade_preds_ch)
  cascade_preds_ch[, model:='full data']
  
  old_cascade_preds_ch <- as.data.table(readRDS(paste0('FILEPATH', prev_best, '/ifr_preds_1yr_3set.rds')))[, `:=` (result='Super region fit', model='Super region fit with all data')]
  old_global_preds <-  as.data.table(readRDS(paste0('FILEPATH', prev_best, '/ifr_preds_1yr_global.rds')))[, `:=` (result='Global fit', model='Super region fit with all data')]
  
  pdf(paste0(output_dir, 'diagnostics/ifr_age_cascading_spline_preds.pdf'), 14, 8)
  
  #then, plot by super-region
  for (hi in unique(cascade_preds_ch$custom_hierarchy)){
    gg <- ggplot()+
      geom_point(data=model_data[custom_hierarchy==hi], aes(x=age_mid, y=logit_ifr, color=as.factor(bias)), shape=1, alpha=0.5)+
      geom_line(data=cascade_preds_ch[custom_hierarchy==hi], aes(x=age_start, y=pred, linetype=result), color='blue', size=1.1)+
      geom_line(data=df_pred_glob, aes(x=age_start, y=pred_glob, linetype=result), color='blue', alpha=0.5)+
      geom_vline(xintercept=best_knots1[2]*100, linetype='dashed', color='grey')+
      geom_vline(xintercept=best_knots1[3]*100, linetype='dashed', color='grey')+
      geom_vline(xintercept=best_knots1[4]*100, linetype='dashed', color='grey')+
      geom_vline(xintercept=best_knots1[5]*100, linetype='dashed', color='grey')+
      theme_bw()+
      scale_linetype_manual(values=c('Global fit'='dashed', 'Super region fit'='solid', 'Super region fit with all data'='solid'))+
      scale_color_manual(values=c('1'='red', '0'='black'))+
      labs(title=paste0(hi, ': COVID-19 IFR by Age (Constant SE)'),
           y='IFR, logit scale', x='Age midpoint', 
           caption='Nursing Home Locs: Belgium, Canada, Ireland, Connecticut, Delaware, Minnesota, New Hampshire, North Carolina\n Ohio, Oregon, Washington, Rhode Island, New York',
           size = 'Inverse SE',
           color='Bias', 
           linetype = 'Model')
    print(gg)
  }
  
  dev.off()
  
  pdf(paste0(output_dir, 'diagnostics/ifr_age_cascading_spline_preds_compare_toprev_withglobal.pdf'), 14, 8)
  
  old_cascade_preds_ch$model <- 'previous'
  cascade_preds_ch$model <- 'new'
  old_global_preds$model <- 'previous'
  df_pred_glob$model <- 'new'
  
  #then, plot by super-region
  for (hi in unique(cascade_preds_ch$custom_hierarchy)){
    gg <- ggplot()+
      geom_point(data=model_data[custom_hierarchy==hi], aes(x=age_mid, y=logit_ifr, color=as.factor(bias)), shape=1, alpha=0.5)+
      geom_line(data=cascade_preds_ch[custom_hierarchy==hi], aes(x=age_start, y=pred, linetype=result, color=model), size=1.1)+
      geom_line(data=old_cascade_preds_ch[custom_hierarchy==hi], aes(x=age_start, y=pred, linetype=result, color=model), size=1.1)+
      geom_line(data=df_pred_glob, aes(x=age_start, y=pred_glob, linetype=result, color=model),alpha=1)+
      geom_line(data=old_global_preds, aes(x=age_start, y=pred_glob, linetype=result, color=model),alpha=1)+
      theme_bw()+
      scale_linetype_manual(values=c('Global fit'='dashed', 'Super region fit'='solid', 'Super region fit with all data'='solid'))+
      scale_color_manual(values=c('1'='red', '0'='black', 'new'='blue', 'previous'= 'orange'))+
      labs(title=paste0(hi, ': COVID-19 IFR by Age (Constant SE)'),
           y='IFR, logit scale', x='Age midpoint', 
           caption='Nursing Home Locs: Belgium, Canada, Ireland, Connecticut, Delaware, Minnesota, New Hampshire, North Carolina\n Ohio, Oregon, Washington, Rhode Island, New York',
           size = 'Inverse SE',
           color='Bias', 
           linetype = 'Model')
    print(gg)
  }
  
  dev.off()
}
