fit_global <- function(model_data, mort_knot_path) {
  df <- model_data %>%
    mutate(age_mid = (age_start + age_end) / 2) %>%
    mutate(
      is_indonesia = ifelse(substr(data_filename,1,3) == "IDN", TRUE, FALSE) ) %>%
    filter(!(age_start_old == 0 & age_end_old == 125)) %>%
    mutate(constant_se = 1) %>%
    filter(location != "South Sumatra")
  
  
  dat_loc <- MRData()
  dat_loc$load_df(
    data = df,
    col_obs = "logit_mort", col_obs_se = "logit_mort_se",
    col_covs = list("age_start", "age_end"), col_study_id = "location_id"
  )
  
  cov_intercept <- LinearCovModel("intercept", use_re = TRUE, prior_gamma_uniform = array(c(0, 0)) )
  
  cov_spline <- LinearCovModel(
    alt_cov = list("age_start", "age_end"),
    use_re = FALSE,
    use_spline = TRUE,
    spline_knots = array(seq(0, 1, length.out = 6)),
    spline_degree = 3L,
    spline_knots_type = 'domain',
    spline_l_linear = FALSE,
    spline_r_linear = TRUE,
    prior_spline_monotonicity = 'increasing',
    prior_spline_monotonicity_domain = tuple(0.1, 1.0)
  )
  
  kd <- c(0.0, 1.0) # knot domain
  is <- c(0.1, 0.9) # interval sizes
  
  knot_bounds1 <- do.call("rbind", list(c(0.05, 0.15), c(0.15,0.3), c(0.3,0.5), c(0.5,0.7)))
  
  np <- import("numpy")
  np$random$seed(123L)
  
  n_knots <- 4L
  knots_samples1 <- utils$sample_knots(
    num_intervals = n_knots + 1L,
    knot_bounds = knot_bounds1,
    interval_sizes = do.call("rbind", lapply(1:(n_knots+1), function(x) is)),
    num_samples = 30L
  )
  
  mod3_ensemble <- MRBeRT(
    data = dat_loc,
    ensemble_cov_model = cov_spline,
    cov_models = list(cov_intercept),
    ensemble_knots = knots_samples1
  )
  
  mod3_ensemble$fit_model(inner_print_level = 3L, inner_max_iter = 500L)
  
  best_knots1 <- mod3_ensemble$ensemble_knots[which.max(mod3_ensemble$weights), ]
  saveRDS(best_knots1, mort_knot_path)

  ##### fit global model using optimized knots
  mod3 <- MRBRT(
    data = dat_loc,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_gamma_uniform = array(c(0, 0)) ),
      LinearCovModel(
        alt_cov = list("age_start", "age_end"),
        use_re = FALSE,
        use_spline = TRUE,
        spline_knots = array(best_knots1),
        spline_degree = 3L,
        spline_knots_type = 'domain',
        spline_l_linear = FALSE,
        spline_r_linear = TRUE,
        prior_spline_monotonicity = 'increasing',
        prior_spline_monotonicity_domain = tuple(0.1, 1.0)
      )
    )
  )
  
  mod3$fit_model(inner_print_level = 5L, inner_max_iter = 5000L)
  
return(mod3)
}

fit_cascade <- function(mod_global, model_data, thetas, output_dir) {
  model_label_tmp <- "cascade_mortality"
  cascade_fit <- run_spline_cascade(
    stage1_model_object = mod_global, 
    df = model_data, 
    col_obs = "logit_mort",
    col_obs_se = "logit_mort_se", 
    col_study_id = "location_id", 
    stage_id_vars = c("super_region_name", "location_id"),
    thetas = thetas,
    gaussian_prior = TRUE,
    output_dir = output_dir,
    model_label = model_label_tmp, 
    overwrite_previous = TRUE
  )
  return(cascade_fit)
}
