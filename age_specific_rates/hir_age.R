hir_age <- function(input_dir,  output_dir) {
  exclude_loc_ids <- c(60886) # western washington doesn't have population or life expectancy
  
  df_sources_in <- read.csv(file.path(input_dir, "serological_age_sources.csv"))
  df_in <- read.csv(file.path(input_dir, "raw_formatted/serological_age.csv"), as.is = T)
  df_lt3 <- read.csv(file.path(output_dir, "life_expectancy_by_1yr_age.csv"))
  df_pop_in <- read.csv(file.path(input_dir, "output_measures/population/all_populations.csv"))
  df_locmeta <- read.csv(file.path(input_dir, "locations/covariate_with_aggregates_hierarchy.csv"))
  all_age_hosps <- read.csv(all_age_data_path, as.is = TRUE)
  fit_hosp_age_pattern <- readRDS(file.path(
    output_dir,
    "cascade_hosp_thetas3_gaussian/cascade_hosp_thetas3_gaussian.RDS"
  ))
  df_hosp_stage_ids <- read.csv(file.path(fit_hosp_age_pattern$working_dir, "df_stage_ids.csv"), as.is = TRUE)
  knots_tmp <- readRDS(file.path(output_dir, "best_ifr_knots_pre_covselection.RDS"))
  
  use_goldstandard_only <- FALSE
  save_hosp_date <- FALSE
  predict_for_loctime <- FALSE
  
  if (use_goldstandard_only) {
    goldstandard_filenames <- as.character(
      filter(df_sources_in, bias == 0 & geo_accordance == 1)[, "data_filename"]
    )
    df_in <- filter(df_in, data_filename %in% goldstandard_filenames)
  }

  #####
  # 1. prep serological data
  df_sources <- df_sources_in %>%
    mutate(
      day = substr(date_epi, 1, 2),
      month = substr(date_epi, 4, 5),
      year = substr(date_epi, 7, 11),
      date = ymd(paste0(year, "-", month, "-", day)) ) %>%
    select(-day, -month, -year)
  
  logit <- function(p) log(p/(1-p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # -- calculate standard errors
  #    key variables are 'prev_logit' and 'se_logit'
  #    others are created for plotting
  # -- reclassify zero-valued observations to 1st percentile
  
  # recode to numeric
  df_in$lower <- as.numeric(df_in$lower)
  df_in$upper <- as.numeric(df_in$upper)
  df_in$mean <- as.numeric(df_in$mean)
  df_in$sample_size <- as.numeric(df_in$sample_size)
  
  lower_first_pctl <- quantile(df_in$lower[df_in$lower > 0], probs = 0.01, na.rm = TRUE)
  mean_first_pctl <- quantile(df_in$mean[df_in$mean > 0], probs = 0.01, na.rm = TRUE)
  
  vacc_dt <- fread(paste0(vacc_dir, "slow_scenario_vaccine_coverage.csv"))
  min_vacc_date <- vacc_dt[cumulative_all_effective > 0][, .(min_vacc_date = min(date)), by = .(location_id)]
  df1 <- df_in %>%
    filter(mean_percentage != 0) %>%
    mutate(
      # Adjust 0's and 1's to avoid infinite logit
      lower_is_adjusted = ifelse(lower < lower_first_pctl, TRUE, FALSE),
      mean_is_adjusted = ifelse(mean_percentage < mean_first_pctl, TRUE, FALSE),
      lower = ifelse(lower < lower_first_pctl, lower_first_pctl, lower),
      mean_percentage = ifelse(mean_percentage < mean_first_pctl, mean_first_pctl, mean_percentage),
      study_id = gsub(".csv", "", as.character(data_filename)),
      age_start = as.numeric(age_start),
      age_mid = (age_start + age_end) / 2,
      prev = mean_percentage/100,
      prev_logit = logit(prev),
      lower_logit = logit(lower/100),
      upper_logit = logit(upper/100),
      standard_error_from_ss = ifelse(is.na(lower), sqrt(prev * (1-prev)/ sample_size), NA),
      se_logit = ifelse(is.na(standard_error_from_ss), (upper_logit - lower_logit) / 3.92, NA) ) %>%
    left_join(select(df_sources, data_filename, date, location_id, location, bias, geo_accordance, bias_type), by = "data_filename") %>%
    filter(!(location_id %in% exclude_loc_ids)) %>%
    filter(!(is.na(mean_percentage))) %>%
    filter(se_logit > 0) %>%
    filter(data_filename != "JPN_Kobe_Doi2020.csv") %>% # todo: remove this, solve SE issue
    left_join(min_vacc_date) %>%
    filter(date < min_vacc_date) #filter to pre-vac rollout by location
  
  df1 <- filter(df1, age_start <= age_end) # TODO remove
  
  df1[!is.na(df1$standard_error_from_ss), "se_logit"] <- mrbrt002::linear_to_logit(
    mean = array(df1[!is.na(df1$standard_error_from_ss), "prev"]),
    sd = array(df1[!is.na(df1$standard_error_from_ss), "standard_error_from_ss"])
  )[[2]]
  
  df1$se <- mrbrt002::logit_to_linear(
    mean = array(df1$prev_logit),
    sd = array(df1$se_logit)
  )[[2]]
  
  # -- 'df1_tmp' is for plotting
  df1_tmp <- df1 %>%
    mutate(
      prev_logit_lo = prev_logit - 1.96*se_logit,
      prev_logit_hi = prev_logit + 1.96*se_logit,
      prev_lo_2 = inv_logit(prev_logit_lo),
      prev_hi_2 = inv_logit(prev_logit_hi)
    ) %>%
    filter(!is.na(mean_percentage) & !is.na(lower) & !is.na(upper))
  
  
  #####
  # MR-BRT seroprevalence model
  #####
  # 2. plot serological data
  
  # -- separately by study
  
  pdf(file.path(output_dir, "diagnostics/seroprevalence_inputdata_plots_hosp.pdf"))
  
  for (id in unique(df1_tmp$data_filename)) {
    df1_tmploc <- filter(df1_tmp, data_filename == id)
    yvals <- c(df1_tmploc$prev_lo_2, df1_tmploc$prev_hi_2)
    ylims <- c(min(yvals), max(yvals))
    if(NaN %in% ylims) ylims <- c(0,1)
    
    plot(
      x = NULL,  y = NULL,
      main = paste0(df1_tmploc[1, "location"], " - \"", gsub(".csv", "", id), "\""),
      xlab = "Age", ylab = "Seroprevalence (proportion)",
      xlim = c(0, 125),  ylim = ylims
    )
    
    for (i in 1:nrow(df1_tmploc)) {
      tmp <- df1_tmploc[i, ]
      with(tmp, lines(x = c(age_start, age_end), y = c(prev, prev)))
      with(tmp, lines(
        x = c(age_mid, age_mid),
        y = c(prev_lo_2, prev_hi_2)
        # col = ifelse(is.na(standard_error_from_ss), "black", "orange")
      ))
    }
  }
  
  dev.off()
  
  
  # -- plot all studies together
  colors1 <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(unique(df1$data_filename)))
  
  df_colors1 <- data.frame(stringsAsFactors = FALSE,
                           file = as.character(unique(df1$data_filename)),
                           color1 = colors1
  )
  
  df1_tmp_2 <- df1_tmp %>%
    left_join(df_colors1, by = c("data_filename" = "file"))
  
  yvals <- df1_tmp_2$prev
  ylims <- c(min(yvals), max(yvals))
  
  pdf(file.path(output_dir, "diagnostics/seroprevalence_inputdata_alllocs_oneplot_hosp.pdf"))
  plot(
    x = NULL,
    y = NULL,
    xlim = c(0, 125),
    ylim = ylims,
    main = "Seroprevalence by study",
    xlab = "Age",
    ylab = "Seroprevalence (proportion)"
  )
  
  for (x in unique(df1_tmp_2$data_filename)) {
    df_tmp <- filter(df1_tmp_2, data_filename == x)
    with(df_tmp, lines(age_mid, prev, col = color1))
  }
  
  
  plot.new()
  legend("topright", legend = df_colors1$file, col = df_colors1$color1, pch = 15, cex = 0.9)
  dev.off()
  
  
  
  #####
  # 3. prep population data
  df_pop <- df_pop_in %>%
    filter((age_group_years_end - age_group_years_start) == 1 | age_group_id == 235) %>%
    filter(sex_id == 3)
  
  # check for data points where age_end == age_start of the next group
  #    and subtract 1 from age_end if so;
  # -- if age_end is non-integer, round down so the whole 1-year age group is included
  df2 <- df1 %>%
    arrange(location_id, date, age_start) %>%
    group_by(location_id, date) %>%
    mutate(
      same_adjacent_age = age_end == lead(age_start, 1),
      age_end = ifelse(same_adjacent_age & row_number() != n(), age_end - 1, age_end),
      same_adjacent_age2 = age_end == lead(age_start, 1),
      age_end = floor(age_end)) %>%
    as.data.frame(.)
  
  # get population for each unique location-age group combination
  # then merge it back on to the analytic dataset
  df_age_groups <- df2 %>%
    select(location_id, age_start, age_end) %>%
    filter(!duplicated(.))
  
  age_group_dat <- lapply(1:nrow(df_age_groups), function(i) {
    # i <- 10 # dev
    # cat(i, "\n")
    ages_tmp <- df_age_groups[i, "age_start"]:df_age_groups[i, "age_end"]
    
    loc_id_tmp <- df_age_groups[i, "location_id"]
    if (loc_id_tmp == 60886) { 
      loc_id <- 570
    } else {
      loc_id <- loc_id_tmp
    }
    
    df_pop_tmp <- df_pop %>% filter(location_id == loc_id)
    
    if (any(ages_tmp >= 95)) {
      df_probs <- data.frame(prob = seq(0.01, 0.99, by = 0.01))
      age95plus_pop <- df_pop_tmp[df_pop_tmp$age_group_years_start == 95, "population"]
      lt_tmp <- filter(df_lt3, location_id == df_pop_tmp[1, "location_id"] & age_1yr == 95 )[, "mean"]
      if (length(lt_tmp) == 0) {
        lt_tmp <- filter(df_lt3, location_id == 1 & age_1yr == 95 )[, "mean"]
      }
      df_probs$age_ex <- 95 + lt_tmp
      
      get_avg_age <- function(yr_survival_prob) {
        survival_probs_tmp <- sapply(1:length(95:100), function(x) yr_survival_prob^x)
        weighted.mean(x = (95:100)+0.5, w = survival_probs_tmp)
      }
      df_probs$avg_age <- sapply(df_probs$prob, get_avg_age)
      survival_prob <- df_probs[which.min(abs(df_probs$age_ex - df_probs$avg_age)), "prob"]
      survival_probs <- sapply(1:length(95:100), function(x) survival_prob^x)
      pops <- (survival_probs / sum(survival_probs)) * age95plus_pop
      
      df_pop_tmp <- do.call("rbind", list(
        filter(df_pop_tmp, age_group_years_start %in% 0:94),
        filter(df_pop_tmp, age_group_years_start == 95)[rep(1, length(95:100)), ] %>%
          mutate(
            population = pops, # assuming annual mortality probability consistent with 95+ population and e_x
            age_group_years_start = 95:100,
            age_group_years_end = (95:100) + 1
          )))
    }
    
    df_pop_tmp2 <- df_pop_tmp %>%
      filter(age_group_years_start %in% ages_tmp)
    
    list(
      pop = sum(df_pop_tmp2$population),
      mean_age = with(df_pop_tmp2, weighted.mean(x = age_group_years_start + 0.5, w = population, na.rm = TRUE))
    )
  })
  
  df_age_groups$pop <- sapply(age_group_dat, function(x) x$pop)
  df_age_groups$mean_age <- sapply(age_group_dat, function(x) x$mean_age)
  
  df3 <- df2 %>%
    left_join(df_age_groups, by = c("location_id", "age_start", "age_end")) %>%
    group_by(data_filename) %>%
    mutate(
      age_start_orig = age_start, 
      age_end_orig = age_end ) %>%
    mutate(
      age_end = ifelse(age_end_orig == 125, age_start + 2*(mean_age - age_start), age_end),
      age_mid_bak = age_mid,
      age_mid = (age_start + age_end) / 2 ) %>%
    as.data.frame(.) %>%
    mutate(
      infections = prev * pop,
      infections_se = se * pop,
      time_minus4 = date - 4 )
  
  #####
  # 4. hospitalization data
  # -- find all-age hosp at time t-4 for each observation in the seroprevalence data
  df_hosps_in <- all_age_hosps %>%
    mutate(
      date = ymd(date),
      allage_hosps = hospitalizations2
    )
  
  
  df4 <- df3 %>%
    left_join(df_hosps_in[, c("location_id", "date", "allage_hosps")],
              by = c("location_id", "time_minus4" = "date")) %>%
    mutate(
      total_hosps = allage_hosps ) %>%
    filter(!is.na(total_hosps)) %>%
    filter(total_hosps > 0) %>%
    left_join(df_locmeta[, c("location_id", "super_region_name")]) %>%
    mutate(
      super_region_name = as.character(super_region_name),
      days_from_march1 = as.integer(ymd(date) - ymd("2020-03-01") )
    )
  if (save_hosp_date) {
    df4_tmp <- df4 %>% 
      select(data_filename, location, geo_accordance, time_plus8) %>% 
      filter(!duplicated(.)) %>%
      arrange(geo_accordance)
    
    write.csv(
      x = df4_tmp, 
      file = "FILEPATH",
      row.names = FALSE
    )
  }
  
  #####
  # 5. convert all-age hosps to age-specific based on hosp rate age pattern predictions
  #
  df_age_frame <- data.frame(age_start = 0:100) %>%
    mutate(age_end = age_start)
  
  if (predict_for_loctime) {
    
    df_hosp_stage_ids$time_mid <- sapply(1:nrow(df_hosp_stage_ids), function(i) {
      tmp <- df_hosp_stage_ids[i, "loctime_id"]
      mean(as.integer(gsub("time", "", strsplit(tmp, split = "_")[[1]])[3:4]) )
    })
    
    df4$loctime_id <- sapply(1:nrow(df4), function(i) {
      
      loc_id_tmp <- df4[i, "location_id"]
      time_tmp <- df4[i, "days_from_march1"]
      df_stage_id_tmp <- filter(df_hosp_stage_ids, location_id == loc_id_tmp)
      
      if (nrow(df_stage_id_tmp) == 0) {
        return(NA)
      } else {
        nearest_time <- with(df_stage_id_tmp, time_mid[which.min(abs(time_tmp - time_mid))] )
        df_stage_id_tmp$loctime_id[df_stage_id_tmp$time_mid == nearest_time]
      }
    })
    
    
    df_unique_locs <- df4 %>%
      select(super_region_name, location_id, loctime_id) %>%
      filter(!duplicated(.))
    
    df_pred_frame <- do.call("rbind", lapply(1:nrow(df_unique_locs), function(i) {
      if (FALSE) {
        i <- 1
      }
      df_tmp <- df_age_frame
      for (nm in names(df_unique_locs)) {
        df_tmp[, nm] <- df_unique_locs[i, nm]
      }
      return(df_tmp)
    }))
    
  } else {
    
    df_unique_locs <- df4 %>%
      select(super_region_name, location_id) %>%
      filter(!duplicated(.))
    
    df_pred_frame <- do.call("rbind", lapply(1:nrow(df_unique_locs), function(i) {
      if (FALSE) {
        i <- 1
      }
      df_tmp <- df_age_frame
      for (nm in names(df_unique_locs)) {
        df_tmp[, nm] <- df_unique_locs[i, nm]
      }
      return(df_tmp)
    }))
    
  }
  
  pred_hosp_age_pattern <- predict_spline_cascade(fit = fit_hosp_age_pattern, newdata = df_pred_frame)
  
  if (predict_for_loctime) {
    group_vars <- c("super_region_name", "location_id", "loctime_id")
  } else {
    group_vars <- c("super_region_name", "location_id")
  }
  
  df_hosp_age_pattern <- pred_hosp_age_pattern %>%
    mutate(
      location_id = as.integer(location_id) ) %>%
    left_join(df_pop, by = c("location_id", "age_start" = "age_group_years_start")) %>%
    mutate(hosp_pred = inv_logit(pred) ) %>%
    group_by(super_region_name, location_id) %>%
    group_by_at(group_vars) %>%
    mutate(
      weight = hosp_pred * population,
      weight_standardized = weight / sum(weight, na.rm = TRUE) ) %>%
    as.data.frame(.)
  
  
  # # do the age split
  df4$agespecific_hosps <- sapply(1:nrow(df4), function(i) {
    dev <- FALSE
    if (dev) {
      i <- 112
    }
    cat(i, "\n")
    
    df_i <- df4[i, ]
    df_hosp_tmp <- df_hosp_age_pattern %>%
      filter(location_id == df_i$location_id)
    df_hosp_tmp_agesubset <- df_hosp_tmp %>%
      filter(age_start %in% df4[i, "age_start"]:df4[i, "age_end"]) %>%
      mutate(
        total_hosps_tmp = df_i$total_hosps,
        age_specific_hosps_tmp = total_hosps_tmp * weight_standardized
      )
    sum(df_hosp_tmp_agesubset$age_specific_hosps_tmp, na.rm = TRUE)
  })

  #####
  # 6. calculate IHR and corresponding standard error
  # -- and include binary variable for before/after june 15
  #
  df5 <- df4 %>%
    mutate(
      agespecifc_hosps_bak = agespecific_hosps,
      agespecific_hosps = ifelse(agespecific_hosps > 0.98*infections, 0.98 * infections, agespecific_hosps) ) %>%
    mutate(
      mu_x = agespecific_hosps,
      mu_y = infections,
      se_x = 0,
      se_y = infections_se,
      hir = mu_x / mu_y,
      hir_se = sqrt((mu_x^2/mu_y^2) * (se_x^2/mu_x^2 + se_y^2/mu_y^2)) ) %>%
    mutate(
      after_june15 = ifelse(date > ymd("2020-06-15"), 1, 0) ) %>%
    filter(pop > 0) %>%
    filter(hir_se > 0) %>%
    filter(!(age_start == 0 & age_end == 0)) %>%
    mutate(group_id = gsub(".csv", "", data_filename)) %>%
    left_join(
      df_locmeta[, c("location_id", "region_id", "region_name", "super_region_id")], by = "location_id") %>%
    mutate(
      seroprev_logit_se_median = median(se_logit),
      seroprev_logit_se_scaled = se_logit * (1 / seroprev_logit_se_median),
      bias_bak = bias,
      bias = ifelse(bias_type %in% c("clinical laboratory visits", "clinical laboratory visits; adults only"), 0, bias)
    )
  
  
  df5[, c("logit_hir", "logit_hir_se")] <- mrbrt002::linear_to_logit(
    mean = array(df5$hir),
    sd = array(df5$hir_se)
  )
  
  n_u18_bins <- df5 %>%
    group_by(location_id) %>%
    filter(age_end <= 18) %>%
    select(age_start, age_end, location_id) %>%
    unique(.) %>%
    summarize(n_bins = n()) %>%
    mutate(mult_u18_bins = ifelse(n_bins > 1 | location_id > 500, 1, 0))
  
  df5 <- df5 %>%
    left_join(
      n_u18_bins %>% select(location_id, mult_u18_bins),
      by = "location_id"
    ) 
  df5 <- as.data.table(df5)
  df5[is.na(mult_u18_bins)]$mult_u18_bins <- 0
  df5 <- df5 %>% 
    filter(seroprev_logit_se_scaled != min(seroprev_logit_se_scaled))%>% 
    filter(!is.infinite(seroprev_logit_se_scaled)) %>% #removes 6 points from Scotland
    filter(mult_u18_bins == 1) %>%
    filter(upper >= mean_percentage)
  
  #adjust IHR? 
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
  
  df5_adj <- df5 %>%
    mutate(sero_end_date = date) %>%
    mutate(study_id_tmp = paste0(location_id, "__", sero_end_date, "__", geo_accordance, "__", bias_type))
  
  df6 <- df5_adj %>%
    left_join(df_tmp3, by = "study_id_tmp")
  
  #apply adjustments/outliering prior to modeling ifr and sero-age
  df6 <- as.data.table(df6)
  
  #filter according to all-age outlier decisions
  outliers <- df6[is.na(adj) & has_allage_seroprev==1] #n=140
  outliers2 <- df6[geo_accordance==0]
  df6 <- df6[!is.na(adj)]
  
  #remove geoaccordance == 0
  df6 <- df6[geo_accordance==1]
  
  #divide ihr and ihr se by 'adj' (seroprevalence adjustment factor)
  df6[, hir_adj:=hir/adj]
  df6[, hir_se_adj:=hir_se/adj]
  
  df6 <- as.data.frame(df6)
  df6[, c("logit_hir_adj", "logit_hir_se_adj")] <- mrbrt002::linear_to_logit(
    mean = array(df6$hir_adj),
    sd = array(df6$hir_se_adj)
  )
  
  #multiply seroprev and seroprev se by 'adj'; se created in logit space, so first back-transform
  df6[, c('prev_linear', 'se_linear')] <- mrbrt002::logit_to_linear(
    mean=array(df6$prev_logit),
    sd = array(df6$se_logit)
  )
  
  df6$prev_adj <- df6$prev_linear * df6$adj
  df6$se_adj <- df6$se_linear * df6$adj
  
  #retransform seroprev to logit space and create adjusted scaled se
  df6 <- as.data.frame(df6)
  df6[, c("prev_logit_adj", "se_logit_adj")] <- mrbrt002::linear_to_logit(
    mean = array(df6$prev_adj),
    sd = array(df6$se_adj)
  )
  
  df6 <- df6 %>%
    mutate(
      seroprev_logit_se_median = median(se_logit_adj),
      seroprev_logit_se_scaled = se_logit_adj * (1 / seroprev_logit_se_median)
    )

  write.csv(df6, paste0(output_dir, 'ihr_prepped_input_data_outliered_adjusted.csv'), row.names=F)
  
  #####
  # 7. MR-BRT model
  dat1 <- MRData()
  dat1$load_df(
    data = df6 ,
    col_obs = "logit_hir_adj", 
    col_obs_se = "seroprev_logit_se_scaled",
    col_covs = list("age_start", "age_end", "bias"), 
    col_study_id = "location_id"
  )
  
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
        prior_spline_monotonicity = 'increasing',
        prior_spline_monotonicity_domain = tuple(0.2181977, 1.0), #ages 20+; age range is 0-91.66; 20/91.66 = 0.2181977
        spline_knots_type = "domain",
        spline_r_linear = TRUE,
        spline_l_linear = FALSE
      )
    ) 
  )
  
  mod2$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

  # make predictions
  df_pred1 <- data.frame(age_start = seq(0, 100, by = 1)) %>%
    mutate(
      age_end = age_start,
      bias = 0
    )
  
  dat_pred1 <- MRData()
  
  dat_pred1$load_df(
    data = df_pred1,
    col_covs=list('age_start', 'age_end', 'bias')
  )
  
  df_pred1$pred <- mod2$predict(data = dat_pred1)
  df_pred1$pred_invlogit <- inv_logit(df_pred1$pred)
  
  df_hir_new <- df_pred1 %>%
    mutate(
      age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE) ) %>%
    filter(!is.na(age_group)) %>%
    group_by(age_group) %>%
    summarize(hir_new = mean(pred_invlogit)) %>%
    mutate(hir_new_logit = logit(hir_new)) %>%
    as.data.frame(.)
  df_hir_new$age_start <- as.integer(gsub("\\[", "", sapply(as.character(df_hir_new$age_group), function(x) strsplit(x, split = ",")[[1]][1])))
  df_hir_new$age_end <- as.integer(as.numeric(gsub("\\)", "", sapply(as.character(df_hir_new$age_group), function(x) strsplit(x, split = ",")[[1]][2]))) - 1)
  
  df_hir_new_out <- df_hir_new %>%
    select(age_group_start = age_start, age_group_end = age_end, logit_hir = hir_new_logit, hir = hir_new)
  
  
  # Save
  write.csv(df6, file.path(output_dir, 'hir_prepped_input_data.csv'))
  write.csv(df_pred1, file.path(output_dir, 'hir_preds_1yr.csv'))
  write.csv(df_hir_new_out, file.path(output_dir, 'hir_preds_5yr.csv'))
  
  #####
  # -- plot, log space
  pdf(file.path(output_dir, "diagnostics/log_ihr_age.pdf"))
  df5$log10_hir <- log10(inv_logit(df5$logit_hir))
  with(df5, plot(
    x = age_mid, 
    y = log10_hir, 
    col = ifelse(bias == 1, "red", adjustcolor("black", alpha.f = 0.5)),
    main = "IHR by age", 
    cex = 1 / seroprev_logit_se_scaled*0.3,
    xlim = c(0, 110),
    ylim = c(-5, 0),
    ylab = "IHR, log base 10",
    xlab = "Age"
  ))
  with(df5, arrows(y0=log10_hir, y1=log10_hir, x0=age_start, x1=age_end, code=3, angle=90, length=0.01, lwd=0.1, 
                   col = ifelse(bias  == 1, "red", adjustcolor("black", alpha.f = 0.5))))
  with(df_pred1, lines(age_start, log10(inv_logit(pred)), lty = 1, lwd = 2, col = "black"))
  df_best <- fread('FILEPATH/hir_preds_1yr.csv')
  with(df_best, lines(age_start, log10(inv_logit(pred)), lty = 1, lwd = 2, col = "blue"))
  legend("topleft",
         c("new","previous best"),
         fill=c("black","blue"))
  grid()
  
  with(df5, plot(
    x = age_mid, 
    y = log10_hir, 
    col = ifelse(bias == 1, "red", adjustcolor("black", alpha.f = 0.5)),
    main = "IHR by age", 
    cex = 1 / seroprev_logit_se_scaled*0.3,
    xlim = c(0, 110),
    ylim = c(-5, 0),
    ylab = "IHR, log base 10",
    xlab = "Age"
  ))
  with(df_pred1, lines(age_start, log10(inv_logit(pred)), lty = 1, lwd = 2, col = "black"))
  df_best <- fread('FILEPATH/hir_preds_1yr.csv')
  with(df_best, lines(age_start, log10(inv_logit(pred)), lty = 1, lwd = 2, col = "blue"))
  legend("topleft",
         c("new","previous best"),
         fill=c("black","blue"))
  grid()
  
  
  add_ifr_line <- FALSE
  if (add_ifr_line) {
    ifr_file <- "FILEPATH/ifr_preds_1yr.csv"
    df_ifr_tmp <- read.csv(ifr_file)
    
    df_pred2 <- df_pred1 %>%
      mutate(
        ifr_pred_tmp = df_ifr_tmp$pred_invlogit,
        ifr_pred_tmp_log10 = log10(ifr_pred_tmp),
        dhr = ifr_pred_tmp / pred_invlogit
      )
    
    with(df_pred2, lines(age_start, ifr_pred_tmp_log10, col = "darkorange", lwd = 3))
    legend("bottomright", legend = "IFR", lwd = 3, col = "darkorange")
    with(df_pred2, plot(age_start, dhr, type = "l", xlab = "Age", ylab = "IFR / IHR, normal space"))
    grid()
    
  }
  dev.off()
}