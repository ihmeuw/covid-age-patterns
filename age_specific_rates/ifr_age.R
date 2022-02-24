prep_ifr_age <- function(
  input_dir, mort_age_dir, output_dir, predict_ifr_from_ensemble = F,
  optimize_knots = T, exclude_loc_ids = c(532), use_goldstandard_only = F,
  all_age_data_path
) {
  ifr_covs <- c()
  df1 <- prep_sero_data(input_dir, use_goldstandard_only, exclude_loc_ids)
  plot_sero_data(df1, output_dir)
  pop_data <- add_pop_data(df1, input_dir, output_dir)
  df3 <- pop_data[[1]]
  df_pop <- pop_data[[2]]
  df4 <- add_death_data(df3, all_age_data_path)
  df4 <- add_age_specific_deaths(df4, mort_age_dir, df_pop)
  df5_tmp <- add_covs(df4, ifr_covs)
  df5 <- calc_ifr(df5_tmp, input_dir)

  # Save outputs
  write.csv(df5, file.path(output_dir, 'ifr_prepped_input_data.csv'))
  return(df5)
}

prep_sero_data <- function(input_dir, use_goldstandard_only, exclude_loc_ids) {
  df_sources_in <- read.csv(file.path(input_dir, "serological_age_sources.csv"))
  df_in <- read.csv(file.path(input_dir, "raw_formatted/serological_age.csv"), as.is = T)
  if (use_goldstandard_only) {
    goldstandard_filenames <- as.character(
      filter(df_sources_in, bias == 0 & geo_accordance == 1)[, "data_filename"]
    )
    df_in <- filter(df_in, data_filename %in% goldstandard_filenames)
  }
  ## Subset to only both sexes data
  both_sexes_filenames <- as.character(
    filter(df_sources_in, sex == "both")[, "data_filename"]
  )
  df_in <- filter(df_in, data_filename %in% both_sexes_filenames)
  
  
  df_sources <- df_sources_in %>%
    mutate(
      day = substr(date_epi, 1, 2),
      month = substr(date_epi, 4, 5),
      year = substr(date_epi, 7, 11),
      date = ymd(paste0(year, "-", month, "-", day)) ) %>%
    select(-day, -month, -year)
  
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
  upper_first_pctl <- quantile(df_in$upper[df_in$upper > 0], probs = 0.01, na.rm = TRUE)
  mean_first_pctl <- quantile(df_in$mean[df_in$mean > 0], probs = 0.01, na.rm = TRUE)
  
  df_in_checkmean0 <- df_in %>%
    filter(mean_percentage == 0)
  
  df1 <- df_in %>%
    # filter(mean_percentage != 0) %>%
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
    left_join(select(df_sources, data_filename, date, location_id, location, bias, geo_accordance, bias_type, test_name, survey_series), by = "data_filename") %>%
    filter(!(location_id %in% exclude_loc_ids)) %>%
    filter(!(is.na(mean_percentage)))
  
  table(is.na(df1$se_logit))
  
  df1_checkna <- df1 %>%
    arrange(se_logit) %>%
    select(
      age_start, age_end, mean_percentage, lower, upper, sample_size, 
      prev_logit, se_logit, lower_logit, upper_logit, standard_error_from_ss,
      data_filename, lower_is_adjusted, mean_is_adjusted
    )
  
  df1[!is.na(df1$standard_error_from_ss), "se_logit"] <- mrbrt002::linear_to_logit(
    mean = array(df1[!is.na(df1$standard_error_from_ss), "prev"]),
    sd = array(df1[!is.na(df1$standard_error_from_ss), "standard_error_from_ss"])
  )[[2]]
  
  df1 <- filter(df1, se_logit > 0)
  
  df1$se <- mrbrt002::logit_to_linear(
    mean = array(df1$prev_logit),
    sd = array(df1$se_logit)
  )[[2]]
  
  
  return(df1)
}

plot_sero_data <- function(df1, output_dir) {
  sero_data_plot_path <- file.path(output_dir, "diagnostics/seroprevalence_inputdata_plots.pdf")
  # -- 'df1_tmp' is for plotting
  df1_tmp <- df1 %>%
    mutate(
      prev_logit_lo = prev_logit - 1.96*se_logit,
      prev_logit_hi = prev_logit + 1.96*se_logit,
      prev_lo_2 = inv_logit(prev_logit_lo),
      prev_hi_2 = inv_logit(prev_logit_hi)
    ) %>%
    filter(!is.na(mean_percentage) & !is.na(lower) & !is.na(upper))
  # -- separately by study
  
  pdf(sero_data_plot_path)
  
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
  
  pdf(file.path(output_dir, "diagnostics/seroprevalence_inputdata_alllocs_oneplot.pdf"))
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
}

add_pop_data <- function(df1, input_dir, output_dir) {
  df_pop_in <- read.csv(file.path(input_dir, "output_measures/population/all_populations.csv"))
  df_lt3 <- read.csv(file.path(output_dir, "life_expectancy_by_1yr_age.csv"))
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
    ages_tmp <- df_age_groups[i, "age_start"]:df_age_groups[i, "age_end"]
    df_pop_tmp <- df_pop %>% filter(location_id == df_age_groups[i, "location_id"])
    
    if (any(ages_tmp >= 95)) {
      # print(i)
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
    mutate(age_start_orig = age_start, age_end_orig = age_end ) %>%
    mutate(age_end = ifelse(age_end_orig == 125, age_start + 2*(mean_age - age_start), age_end)) %>% 
    as.data.frame(.) %>%
    mutate(
      infections = prev * pop,
      infections_se = se * pop,
      time_plus8 = date + 8 )
  return(list(df3, df_pop))
}

add_death_data <- function(df3, all_age_data_path) {
  df_deaths_in <- read.csv(all_age_data_path, as.is = TRUE)
  # -- find all-age deaths at time t+8 for each observation in the seroprevalence data
  df_deaths_in2 <- df_deaths_in %>%
    group_by(location_id) %>%
    mutate(
      date = ymd(Date),
      allage_deaths = na.approx(Deaths, na.rm = FALSE) ) %>%
    as.data.frame(.)
  
  df4 <- df3 %>%
    left_join(df_deaths_in2[, c("location_id", "date", "allage_deaths")],
              by = c("location_id", "time_plus8" = "date")) %>%
    mutate(total_deaths =  allage_deaths)
  return(df4)
}

add_age_specific_deaths <- function(df4, mort_age_dir, df_pop) {
  df_mort_age_pattern_in <- read.csv(
    file = file.path(mort_age_dir, "mortality_agepattern_preds_byloc_1yr.csv"),
    stringsAsFactors = FALSE
  )
  df_mort_age_pattern <- df_mort_age_pattern_in %>%
    left_join(df_pop, by = c("location_id", "age_start" = "age_group_years_start")) %>%
    mutate(mort_pred = pred_mr_ratio ) %>%
    group_by(location_id) %>%
    mutate(
      mort_pred_standardized = mort_pred / sum(mort_pred),
      pop_standardized = population / sum(population, na.rm = TRUE),
      weight = mort_pred_standardized * pop_standardized,
      weight_standardized = weight / sum(weight, na.rm = TRUE),
      check1 = sum(weight_standardized, na.rm = TRUE)
    )
  # do the age split
  df4$agespecific_deaths <- sapply(1:nrow(df4), function(i) {
    df_i <- df4[i, ]
    df_mort_tmp <- df_mort_age_pattern %>%
      filter(location_id == df_i$location_id)
    df_mort_tmp_agesubset <- df_mort_tmp %>%
      filter(age_start %in% df4[i, "age_start"]:df4[i, "age_end"]) %>%
      mutate(
        total_deaths_tmp = df_i$total_deaths,
        age_specific_deaths_tmp = total_deaths_tmp * weight_standardized
      )
    sum(df_mort_tmp_agesubset$age_specific_deaths_tmp, na.rm = TRUE)
  })
  return(df4)
}

add_covs <- function(df4, ifr_covs) {
  if(length(ifr_covs) > 0) {
    source("FILEPATH/get_covariate_estimates.R")
    dat_covs <- sapply(ifr_covs, function(x) {  
      df_cov <- get_covariate_estimates(
        covariate_id = x, 
        gbd_round_id = 7, 
        decomp_step = "iterative",
        location_id = unique(df3$location_id),
        year_id = 2019
      )
      df_cov2 <- df_cov %>%
        group_by(location_id) %>%
        summarize(mean_cov_val = mean(mean_value))
      
      names(df_cov2)[names(df_cov2) == "mean_cov_val"] <- names(ifr_covs)[ifr_covs == x]
      
      return(df_cov2)
    }, simplify = FALSE)
    
    
    df5_tmp <- df4 %>%
      mutate(
        location_id_covs = ifelse(location_id == 43866, 101, location_id)
      )
    
    for (nm in names(dat_covs)) {
      df5_tmp <- left_join(df5_tmp, dat_covs[[nm]], by = c("location_id_covs" = "location_id"))
    }
  } else {
    df5_tmp <- df4
  }
  return(df5_tmp)
}

calc_ifr <- function(df5_tmp, input_dir) {
  df_locmeta <- read.csv(file.path(input_dir, "locations/covariate_with_aggregates_hierarchy.csv"))
  vacc_dt <- fread(paste0(vacc_dir, "slow_scenario_vaccine_coverage.csv"))
  min_vacc_date <- vacc_dt[cumulative_all_effective > 0][, .(min_vacc_date = min(date)), by = .(location_id)]
  df5 <- df5_tmp %>%
    mutate(
      mu_x = agespecific_deaths,
      mu_y = agespecific_deaths + infections,
      se_x = 0,
      se_y = infections_se,
      ifr = mu_x / mu_y,
      ifr_se = sqrt((mu_x^2/mu_y^2) * (se_x^2/mu_x^2 + se_y^2/mu_y^2)) ) %>%
    mutate(
      after_june15 = ifelse(date > ymd("2020-06-15"), 1, 0) ) %>%
    left_join(min_vacc_date) %>%
    filter(date < min_vacc_date) %>%
    filter(pop > 0) %>%
    filter(ifr_se > 0) %>%
    filter(!(age_start == 0 & age_end == 0)) %>%
    mutate(age_end = ifelse(age_end > 100, 100, age_end)) %>%
    mutate(group_id = gsub(".csv", "", data_filename)) %>%
    left_join(
      df_locmeta[, c("location_id", "region_id", "region_name", "super_region_id", "super_region_name")], by = "location_id") %>%
    mutate(region_name = as.character(region_name), super_region_name = as.character(super_region_name)) %>%
    mutate(is_in_latin_america = ifelse(super_region_name == "Latin America and Caribbean", 1, 0)) %>%
    mutate(
      bias_bak = bias,
      bias = ifelse(bias_type %in% c("clinical laboratory visits", "clinical laboratory visits; adults only"), 0, bias),
      constant_se = 1,
      seroprev_logit_se_median = median(se_logit),
      seroprev_logit_se_scaled = se_logit * (1 / seroprev_logit_se_median)
    )
  
  df5[, c("logit_ifr", "logit_ifr_se")] <- mrbrt002::linear_to_logit(
    mean = array(df5$ifr),
    sd = array(df5$ifr_se)
  )
  df5$age_mid <- (df5$age_start + df5$age_end) / 2
  
  return(df5)
}

prep_ifr_5yr <- function(df_pred1) {
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