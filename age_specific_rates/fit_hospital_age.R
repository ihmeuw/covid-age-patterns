
fit_hospital_age <- function(input_dir, output_dir) {
  model_label_tmp <- "cascade_hosp_thetas3_gaussian"
  
  df_locmeta <- read.csv(file.path(
    input_dir, "locations/covariate_with_aggregates_hierarchy.csv"
  ), as.is = TRUE)
  df_hosp <- read.csv(file.path(output_dir, "df_hosp_byagedateloc_rawcumulative.csv"), stringsAsFactors = FALSE)

  df_in <- df_hosp %>%
    filter(!is.na(n_hosp_interval)) %>%
    mutate(
      hosp_rate = hospitalizations / pop
    ) %>% # NOTE! this is cumulative hospitalizations, not interval; corresponds better to seroprevalence (cumulative infections)
    filter(hosp_rate > 0) %>%
    mutate(
      age_mid = (age_start + age_end) / 2,
      hosp_rate_se = sqrt((hosp_rate / (1-hosp_rate)) / pop) ) %>%
    
    left_join(df_locmeta[, c("location_id", "super_region_name", "location_name")]) %>%
    filter(!is.nan(hosp_rate_se)) %>%
    group_by(location_id) %>%
    filter(days_from_march1 == max(days_from_march1)) %>%
    as.data.frame(.)
  
  df_in_outliers <- df_in %>%
    group_by(loctime_id) %>%
    summarize(
      hosp_allage = weighted.mean(x = hosp_rate, w = pop) ) %>%
    arrange(desc(abs(hosp_allage))) %>%
    as.data.frame(.) %>%
    group_by(.) %>%
    mutate(
      abs_deviation = abs(hosp_allage - mean(hosp_allage)),
      mad = median(abs_deviation),
      is_gt_3mads = abs_deviation > 3*mad ) %>%
    filter(is_gt_3mads)
  
  df <- df_in %>%
    filter(!loctime_id %in% unique(df_in_outliers$loctime_id)) %>%
    as.data.frame(.) %>%
    filter(!location_name == "New Jersey") %>%
    mutate(constant_se = 1)
  
  df[, c("logit_hosp_rate", "logit_hosp_rate_se")] <- mrbrt002::linear_to_logit(
    mean = array(df$hosp_rate),
    sd = array(df$hosp_rate_se)
  )
  
  dat_loc <- MRData()
  dat_loc$load_df(
    data = df,
    col_obs = "logit_hosp_rate", col_obs_se = "logit_hosp_rate_se",
    col_covs = list("age_start", "age_end"), col_study_id = "loctime_id"
  )
  
  #####
  mod3 <- MRBRT(
    data = dat_loc,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_gamma_uniform = array(c(0, 0)) ),
      LinearCovModel(
        alt_cov = list("age_start", "age_end"),
        use_re = FALSE,
        use_spline = TRUE,
        spline_knots = array(c(0, 0.15, 0.30, 0.45, 0.60, 1)),
        spline_degree = 3L,
        spline_knots_type = 'domain',
        spline_l_linear = FALSE,
        spline_r_linear = TRUE
      )
    )
    # inlier_pct = 0.95
  )
  
  mod3$fit_model(inner_print_level = 5L, inner_max_iter = 2000L)

  #####
  # make predictions
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
  
  df_pred3$pred3 <- mod3$predict(dat_pred3, sort_by_data_id = "data_id")

  system.time({
    cascade_fit1 <- run_spline_cascade(
      stage1_model_object = mod3, 
      df = df, 
      col_obs = "logit_hosp_rate",
      col_obs_se = "logit_hosp_rate_se", 
      col_study_id = "location_id", 
      stage_id_vars = c("super_region_name", "location_id"),
      thetas = c(3,3),
      gaussian_prior = TRUE,
      output_dir = output_dir,
      model_label = model_label_tmp, 
      overwrite_previous = TRUE
    )
    
  })
  
  # -- predict for super-regions
  df_pred_sr <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100, 
    super_region_name = unique(df$super_region_name) ) %>%
    mutate(
      age_end = age_start,
      location_id = NA,
      loctime_id = NA,
      age_mid = (age_start + age_end) / 2
    )
  
  
  model_path_tmp <- file.path(output_dir, model_label_tmp)
  cascade_fit1 <- list(working_dir = model_path_tmp)
  cascade_preds_sr <- predict_spline_cascade(fit = cascade_fit1, newdata = df_pred_sr)
  
  # -- predict for locations
  df_pred_loc <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100, 
    location_id = unique(df$location_id) ) %>%
    mutate(
      loctime_id = NA,
      age_end = age_start,
      age_mid = (age_start + age_end) / 2 ) %>%
    left_join(df %>% .[, c("super_region_name", "location_id", "location_name")] %>% filter(!duplicated(.)))
  
  cascade_preds_loc <- predict_spline_cascade(fit = cascade_fit1, newdata = df_pred_loc)
  
  # -- predict for location-times
  df_pred_loctime <- expand.grid(
    stringsAsFactors = FALSE,
    age_start = 0:100,
    loctime_id = unique(df$loctime_id) ) %>%
    mutate(
      age_end = age_start,
      age_mid = (age_start + age_end) / 2 ) %>%
    left_join(df %>% select(super_region_name, location_id, location_name, loctime_id) %>% filter(!duplicated(.)))
  
  cascade_preds_loctime <- predict_spline_cascade(fit = cascade_fit1, newdata = df_pred_loctime)
  
  saveRDS(cascade_fit1, file.path(output_dir, model_label_tmp, paste0(model_label_tmp, ".RDS")))
  
  
  #####
  
  
  pdf(file.path(output_dir, "/diagnostics/hospitalizations_plot_superregion_preds.pdf"))
  
  layout(1)
  mod3_knots <- mod3$cov_models[[2]]$spline_knots
  with(df, plot(
    x = age_mid,
    y = logit_hosp_rate,
    xlab = "Age",
    ylab = "Cumulative hospitalizations, logit scale",
    xlim = c(min(mod3_knots), max(mod3_knots)),
    cex = 1/logit_hosp_rate_se * 0.2,
    col = ifelse(location_name == "Netherlands", "red", ifelse(location_name == "England", "blue", adjustcolor("black", alpha.f = 0.25)))
  ))
  grid()
  
  with(df_pred3, lines(age_start, pred3, type = "l", col = "darkorange", lwd = 3))
  for (k in mod3_knots) abline(v = k)
  
  #####
  
  sr_names <- as.character(unique(df_pred_sr$super_region_name))
  colors1 <- RColorBrewer::brewer.pal(n = length(sr_names), name = "Spectral")
  df_colors <- data.frame(sr = sr_names, color1 = colors1, stringsAsFactors = FALSE)
  
  # plot global as thick black line, covering original
  with(df_pred3, lines(age_start, pred3, col = "black", lwd = 4))
  
  for (i in 1:length(sr_names)) {
    dev <- FALSE
    if (dev) {
      i <- 1
    }
    nm <- sr_names[i]
    df_tmp <- filter(cascade_preds_sr, super_region_name == nm)
    color_tmp <- df_colors[df_colors$sr == nm, "color1"]
    with(df_tmp, lines(age_start, pred, col = color_tmp, lwd = 3))
  }
  
  
  legend("bottomright", legend = c(df_colors$sr, "Global"), pch = 15, col = c(df_colors$color1, "#000000"), cex = 0.8)
  
  dev.off()
  
  outputs <- c("df", "df_pred3", "cascade_preds_sr", "cascade_preds_loc", "cascade_preds_loctime")
  outputs2 <- paste0("hosp_", outputs)
  
  for (i in 1:length(outputs)) {
    write.csv(get(outputs[i]), file.path(output_dir, model_label_tmp, paste0(outputs2[i], ".csv")), row.names = FALSE)
  }
  saveRDS(cascade_fit1, file.path(output_dir, "hosp_cascade_fit1.RDS"))
}