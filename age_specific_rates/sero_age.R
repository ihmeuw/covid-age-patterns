fit_sero_age <- function(df5, output_dir) {
  df5_mrbrt <- df5[, c("prev_logit", "se_logit", "prev_logit_adj", "se_logit_adj", "age_start", "age_end", "age_mid", "study_id")] %>%
    mutate(se_logit = 1)
  dat5 <- MRData()
  dat5$load_df(
    data = df5_mrbrt,  
    col_obs = "prev_logit_adj", 
    col_obs_se = "se_logit_adj",
    col_covs = list("age_start", "age_end"), 
    col_study_id = "study_id" )
  
  mod5 <- MRBRT(
    data = dat5,
    cov_models = list(
      LinearCovModel("intercept", use_re = TRUE, prior_gamma_uniform = array(c(0,0))),
      LinearCovModel(
        alt_cov = list("age_start", "age_end"),
        use_spline = TRUE,
        spline_knots = array(c(0, 0.2, 0.6, 1)),
        spline_degree = 3L,
        spline_knots_type = 'domain',
        spline_r_linear = TRUE,
        spline_l_linear = TRUE
      )
    )
  )
  
  mod5$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  
  # -- save model object for use as a global cascade prior
  py_save_object(
    object = mod5,
    filename = file.path(output_dir, "global_seroprevalence_fit.pkl"),
    pickle = "dill"
  )

  df_pred5 <- data.frame(age_start = seq(0, 100, by = 1)) %>%
    mutate(age_end = age_start)
  
  dat_pred5 <- MRData()
  
  dat_pred5$load_df(
    data = df_pred5, 
    col_covs=list("age_start", "age_end")
  )
  samples5 <- mrbrt002::core$other_sampling$sample_simple_lme_beta(sample_size = 1000L, model = mod5)
  draws5 <- mod5$create_draws(
    data = dat_pred5,
    beta_samples = samples5,
    gamma_samples = matrix(rep(0, nrow(samples5)), ncol = 1),
    random_study = FALSE )
  
  df_pred5$pred5 <- mod5$predict(data = dat_pred5)
  df_pred5$pred5_lo <- apply(draws5, 1, function(x) quantile(x, 0.025))
  df_pred5$pred5_hi <- apply(draws5, 1, function(x) quantile(x, 0.975))
  
  df_pred5$pred5_invlogit <- inv_logit(df_pred5$pred5)
  df_pred5$pred5_lo_invlogit <- inv_logit(apply(draws5, 1, function(x) quantile(x, 0.025)))
  df_pred5$pred5_hi_invlogit <- inv_logit(apply(draws5, 1, function(x) quantile(x, 0.975)))
  
  return(list(df_pred5, mod5, draws5))
}

prep_sero_5yr <- function(df_pred5) {
  df_seroprev_new <- df_pred5 %>%
    mutate(
      age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE) ) %>%
    filter(!is.na(age_group)) %>%
    group_by(age_group) %>%
    summarize(seroprev_new = mean(pred5_invlogit)) %>%
    mutate(seroprev_new_logit = logit(seroprev_new)) %>%
    as.data.frame(.)
  df_seroprev_new$age_start <- as.integer(gsub("\\[", "", sapply(as.character(df_seroprev_new$age_group), function(x) strsplit(x, split = ",")[[1]][1])))
  df_seroprev_new$age_end <- as.integer(as.numeric(gsub("\\)", "", sapply(as.character(df_seroprev_new$age_group), function(x) strsplit(x, split = ",")[[1]][2]))) - 1)
  
  df_seroprev_new_out <- df_seroprev_new %>%
    select(age_group_start = age_start, age_group_end = age_end, logit_seroprev = seroprev_new_logit, seroprev = seroprev_new)
  
  return(df_seroprev_new_out)
}

plot_sero_age_fit <- function(df5, df_pred5, mod5, draws5, output_dir) {
  add_ui <- function(dat, x_var, lo_var, hi_var, color = "darkblue", opacity = 0.2) {
    polygon(
      x = c(dat[, x_var], rev(dat[, x_var])),
      y = c(dat[, lo_var], rev(dat[, hi_var])),
      col = adjustcolor(col = color, alpha.f = opacity), border = FALSE
    )
  }
  df5_mrbrt <- df5[, c("prev_logit", "se_logit", "age_start", "age_end", "age_mid", "study_id")] %>%
    mutate(se_logit = 1)
  pdf(file.path(output_dir, "diagnostics/sero_pred_age.pdf"))
  with(df5_mrbrt, plot(age_mid, prev_logit, pch = 16, cex = 0.3))
  for (stud in unique(df5$study_id)) {
    with(filter(df5, study_id == stud), lines(age_mid, prev_logit, lwd = 1/sqrt(sum(se_logit^2)) * 0.4, col = adjustcolor("black", alpha.f = 0.6)))
  }
  with(df_pred5, lines(age_start, pred5, col = "darkorange", lwd = 3))
  knots5 <- mod5$cov_models[[2]]$spline_knots
  for (k in knots5) abline(v = k)
  
  add_ui(df_pred5, "age_start", "pred5_lo", "pred5_hi", color = "darkorange")
  
  # -- normal space version
  with(df5_mrbrt, plot(age_mid, inv_logit(prev_logit), ylim = c(0, 0.1), pch = 16, cex = 0.3, ylab = "Seroprevalence", xlab = "Age"))
  for (stud in unique(df5$study_id)) {
    with(filter(df5, study_id == stud), lines(age_mid, inv_logit(prev_logit), col = adjustcolor("black", alpha.f = 0.3)))
  }
  
  with(df_pred5, lines(age_start, inv_logit(pred5), col = "darkorange", lwd = 3))
  for (k in knots5) abline(v = k)
  
  add_ui(df_pred5, "age_start", "pred5_lo_invlogit", "pred5_hi_invlogit", color = "darkorange")
  
  dev.off()
}


prep_sero_age <- function(df5, output_dir) {
  sero_fit <- fit_sero_age(df5, output_dir)
  df_pred5 <- sero_fit[[1]]
  mod5 <- sero_fit[[2]]
  draws5 <- sero_fit[[3]]
  df_seroprev_new_out <- prep_sero_5yr(df_pred5)
  plot_sero_age_fit(df5, df_pred5, mod5, draws5, output_dir)
  
  # Save outputs
  write.csv(df_pred5, file.path(output_dir, "global_sero_preds_1yr.csv"), row.names = FALSE)
  write.csv(df_pred5, file.path(output_dir, 'seroprev_preds_1yr.csv'))
  write.csv(df_seroprev_new_out, file.path(output_dir, 'seroprev_preds_5yr.csv'))
}