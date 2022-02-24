prep_data <- function(
  df_locmeta, pop_path, nchs_dir, nchs_filename, deaths_path, deaths_add_path,
  deaths_bydate_path, deaths_byagedateloc_path, deaths_byageloc_path, 
  sero_sources_path, exclude_loc_ids
) {
  # Input data
  df_pop_in <- read.csv(pop_path)
  df_nchs_in <- read.csv(file.path(nchs_filename), as.is = TRUE)
  df_deaths_in <- read.csv(deaths_path, as.is = TRUE) %>%
    filter(!duplicated(.))
  df_deaths_add <- read.csv(deaths_add_path, as.is = TRUE) %>%
    filter(!duplicated(.))

  ## Filter pop to both-sexes and single year with ??+ as the terminal age group
  df_pop <- df_pop_in %>%
    filter((age_group_years_end - age_group_years_start) == 1 | age_group_id == 235) %>%
    filter(sex_id == 3)

  ## Clean up and filter out ages below the max age with an NA by location
  df_nchs <- df_nchs_in %>%
    rename(date_orig = date) %>%
    mutate(
      age_end = as.integer(ifelse(age_end == "MAX", 125, age_end)),
      deaths = as.integer(deaths),
      date = dmy(date_orig) ) %>%
    filter(!is.na(age_start)) %>%
    filter(sex %in% c("Male", "Female")) %>%
    group_by(location_id, state, age_start, age_end, date) %>%
    summarize(
      count = n(),
      count_nonmissing = sum(is.na(deaths)), # observations where male or female is NA become NA
      deaths2 = sum(deaths)) %>%
    arrange(location_id, age_start) %>%
    group_by(location_id) %>%
    mutate(row_with_oldest_1to9 = max(row_number() * is.na(deaths2))) %>%
    filter(row_number() > row_with_oldest_1to9) %>%
    as.data.frame(.) %>%
    mutate(sex = "both", data_filename = nchs_filename)

  # check if a location has any death data, regardless of sex
  df_deaths_in2 <- df_deaths_in %>%
    bind_rows(df_deaths_add) %>%
    filter(!duplicated(.)) %>%
    group_by(location, location_id) %>%
    mutate(
      has_death_data = !all(is.na(deaths)) ) %>%
    as.data.frame(.)
  # subset to only observations with death data, and check which type of sex data exist for each location
  df_deaths_sex <- df_deaths_in %>%
    filter(!is.na(deaths)) %>%
    group_by(location, location_id) %>%
    summarize(
      sexvals = paste(sort(unique(sex)), collapse = ",")) %>%
    mutate(
      has_both = ifelse(grepl("both", sexvals), 1, 0),
      has_only_malefemale = ifelse(sexvals == "female,male", 1, 0),
      has_only_male = ifelse(sexvals == "male", 1, 0),
      has_only_female = ifelse(sexvals == "female", 1, 0) ) %>%
    as.data.frame(.)
  # merge on sex metadata
  df_deaths_in3 <- df_deaths_in2 %>%
    left_join(df_deaths_sex %>% select(-location), by = "location_id")
  df_deaths_tmp1 <- df_deaths_in3 %>%
    filter(has_both == 1, sex == "both")
  df_deaths_tmp2 <- df_deaths_in3 %>%
    filter(has_only_malefemale == 1) %>%
    group_by(location_id, location, age_start, age_end, date) %>%
    summarize(
      cases = sum(cases, na.rm = TRUE),
      deaths = sum(deaths, na.rm = TRUE) ) %>%
    group_by(location_id) %>%
    mutate(all_deaths_are_zero = all(deaths == 0)) %>%
    as.data.frame(.)
  df_deaths <- bind_rows(df_deaths_tmp1, df_deaths_tmp2)


  # death data
  df_deaths_date <- df_deaths %>%
    group_by(location_id, data_filename) %>%
    summarize(date = first(date)) %>%
    as.data.frame(.)

  write.csv(df_deaths_date, deaths_bydate_path, row.names = FALSE)

  df_deaths$date_string <- as.character(df_deaths$date)

  # -- subset to locations with population denominators
  tmp_locs <- unique(df_deaths$location_id)
  tmp_locs2 <- tmp_locs[!tmp_locs %in% df_pop_in$location_id]
  df_tmp_locs2 <- filter(df_deaths, location_id %in% tmp_locs2)
  exclude_locs_pop <- unique(df_tmp_locs2$location)
  cat(
    "Excluding locations with no population denominator:",
    paste(exclude_locs_pop, collapse = ", "), "\n"
  )

  # Combine NCHS data
  df_nchs2 <- df_nchs %>%
    mutate(
      data_filename = nchs_filename,
      cases = NA) %>%
    select(location = state, location_id, sex, data_filename, age_start, age_end, cases, deaths = deaths2)

  df_nchs2$date_string <- strsplit(nchs_filename, split = "_")[[1]][4]
  df_nchs2$sex_specific_dataset <- FALSE

  df_deaths_in2 <- df_deaths %>%
    filter(!(location_id %in% unique(df_nchs2$location_id))) %>%
    bind_rows(., df_nchs2) %>%
    arrange(location_id)

  ## Subset to locations with non-NA deaths and population denominators
  df_deaths1 <- df_deaths_in2 %>%
    filter(!is.na(deaths)) %>%
    filter(!location %in% exclude_locs_pop) %>%
    mutate(
      date_string = ifelse(date_string == "2020723", "20200723", date_string),
      date = lubridate::ymd(date_string),
      days_from_march1 = date - ymd("20200301")
    )

  ## Adjust age_end to not overlap next age group
  df_deaths2 <- df_deaths1 %>%
    arrange(location_id, date, age_start) %>%
    group_by(location_id, date) %>%
    mutate(
      same_adjacent_age = age_end == lead(age_start, 1),
      age_end = ifelse(same_adjacent_age & row_number() != n(), age_end - 1, age_end),
      same_adjacent_age2 = age_end == lead(age_start, 1),
      age_end = floor(age_end)) %>%
    as.data.frame(.)

  ## Grab the unique age groups by locations
  df_deaths_age_groups <- df_deaths2 %>%
    select(location_id, age_start, age_end) %>%
    filter(!duplicated(.))


  # -- adjust oldest age group upper bound
  #    to be x + 2*e_x, where e_x is life expectancy at age x
  #
  # -- get life expectancy for each observation, indexed by location and age
  df_lt3 <- prep_life_expectancy(sero_sources_path, df_deaths2)


  age_group_dat <- lapply(1:nrow(df_deaths_age_groups), function(i) {
    print(i)
    ages_tmp <- df_deaths_age_groups[i, "age_start"]:df_deaths_age_groups[i, "age_end"]
    df_pop_tmp <- df_pop %>% filter(location_id == df_deaths_age_groups[i, "location_id"])
    if (any(ages_tmp >= 95)) {

      df_probs <- data.frame(prob = seq(0.01, 0.99, by = 0.01))
      age95plus_pop <- df_pop_tmp[df_pop_tmp$age_group_years_start == 95, "population"]
      lt_tmp <- filter(df_lt3, location_id == df_pop_tmp[1, "location_id"] & age_1yr == 95 )[, "mean"]
      if (nrow(lt_tmp) == 0) {
        lt_tmp <- filter(df_lt3, location_id == 1 & age_1yr == 95 )[, "mean"]
      }
      df_probs$age_ex <- 95 + lt_tmp$mean

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

  df_deaths_age_groups$pop <- sapply(age_group_dat, function(x) x$pop)
  df_deaths_age_groups$mean_age <- sapply(age_group_dat, function(x) x$mean_age)

  ## Assign population and mean age to each age group
  df_deaths3 <- df_deaths2 %>%
    left_join(df_deaths_age_groups, by = c("location_id", "age_start", "age_end"))
  
  write.csv(df_deaths3, deaths_byagedateloc_path)
  
  ##### -------------------------------------------------------------------------------------------------------------------------------------------------------------
  # analytic dataset for cumulative deaths
  ## Subset to oldest location-specific data 

  #filter by vax date
  df_deaths3 <- df_deaths3 %>%
    left_join(min_vacc_date) %>%
    filter(date < min_vacc_date)
  
  df3_cumulative <- df_deaths3 %>%
    group_by(location_id) %>%
    mutate(days_from_march1_adj = ifelse(
      data_filename == "coveragedb_CFR_database", 
      days_from_march1 - 30,
      days_from_march1)) %>% # Using coveragedb data if it is at least 30 days out from most recent IHME data
    filter(days_from_march1_adj == max(days_from_march1_adj)) %>%
    as.data.frame(.)
  
  # Check for two observations on same date and choose IHME
  n_obs <- table(unique(df3_cumulative %>% select(location_id, data_filename)) %>% .$location_id)
  if(any(n_obs > 1)) {
    doub_locs <- as.integer(names(n_obs)[which(n_obs > 1)])
    df3_cumulative_3 <- df3_cumulative %>% filter(!(location_id %in% doub_locs & data_filename == "coveragedb_CFR_database"))
  }
  
  # combine age groups so none have zero deaths
  # -- MR-BRT can't handle observations with zero deaths
  #    so we increase the upper limit of the age group until it
  #    includes at least one death
  
  df3_cumulative_2 <- do.call("rbind", lapply(unique(df3_cumulative$location_id), function(x) {
    df_x <- df3_cumulative %>% filter(location_id == x) %>%
      mutate(aggvar = NA)
    
    counter <- 1
    last_was_zero <- FALSE
    for (row in 1:nrow(df_x)) {
      # row <- 1 # dev
      df_x[row, "aggvar"] <- counter
      if (df_x[row, "deaths"] == 0) {
        last_was_zero <- TRUE
      } else {
        last_was_zero <- FALSE
        counter <- counter + 1
      }
    }
    return(df_x)
  }))
  
  df3_cumulative_3 <- df3_cumulative_2 %>%
    group_by(location_id, aggvar) %>%
    dplyr::summarize(
      location = unique(location),
      data_filename = paste(unique(data_filename), collapse = ","),
      age_start_old = paste(age_start, collapse = ","),
      age_end_old = paste(age_end, collapse = ","),
      deaths_old = paste(deaths, collapse = ","),
      pop_old = paste(pop, collapse = ","),
      mean_age_old = paste(mean_age, collapse = ","),
      age_start = min(age_start),
      age_end = max(age_end),
      pop = sum(pop),
      mean_age = first(mean_age),
      n_groups_combined = n(),
      days_from_march1 = first(days_from_march1),
      deaths = sum(deaths),
      cases = sum(cases),
      date = max(date)) %>%
    mutate(
      mean_age = ifelse(n_groups_combined > 1, NA, mean_age),
      age_mid = (age_start + age_end) / 2 ) %>%
    filter(deaths > 0 )
  
  
  # use as the reference the observation with the smallest standard error
  # -- characteristic of v5 of this script
  # -- must be the lowest standard error in log space! 
  #    otherwise SEs at the extremes will artificially appear smallest
  
  # create a temporary version of 'df3_cumulative_4' for calculating standard error,
  # then carry on as usual
  df3_cumulative_4 <- df3_cumulative_3 %>%
    filter(deaths < pop) %>%
    mutate(
      mort = deaths / pop,
      mort_se = sqrt((mort * (1-mort)) / pop))
  
  df3_cumulative_4[, c("logit_mort", "logit_mort_se")] <- mrbrt002::linear_to_logit(
    mean = array(df3_cumulative_4$mort),
    sd = array(df3_cumulative_4$mort_se)
  )
  
  df3_cumulative_5 <- df3_cumulative_4 %>%
    group_by(data_filename) %>%
    mutate(
      is_oldest = (age_start == max(age_start)) | age_end == 125,
      # is_oldest = (age_start == max(age_start)),
      age_end_bak = age_end,
      age_end = ifelse(is_oldest, age_start + 2*(mean_age - age_start), age_end) ) %>%
    as.data.frame(.) %>%
    left_join(df_locmeta[, c("location_id", "region_id", "region_name", "super_region_id", "super_region_name")]) %>%
    filter(!location_id %in% exclude_loc_ids)
  write.csv(df3_cumulative_5, deaths_byageloc_path, row.names = F)
  
  return(df3_cumulative_5)
}

prep_life_expectancy <- function(sero_sources_path, df_deaths2) {
  
  # get locations for mortality and IFR analyses
  mort_ex_locs <- unique(df_deaths2$location_id)
  ifr_sero_sources <- read.csv(sero_sources_path)
  ifr_ex_locs <- unique(ifr_sero_sources$location_id)
  ex_locs <- c(union(mort_ex_locs, ifr_ex_locs), 1) # include global
  
  # 2017 is most recent in gbd_round_id = 5; round 6 requires a decomp step to be specified
  df_lt <- get_life_table(gbd_round_id = 5, location_id = ex_locs, year_id = 2017,
                          age_group_id = "all", sex_id = 3, life_table_parameter_id = 5)
  
  df_age_meta <- get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>%
    select(age_group_id, age_start = age_group_years_start, age_end = age_group_years_end)
  
  df_lt2 <- df_lt %>%
    left_join(df_age_meta, by = "age_group_id") %>%
    mutate( 
      age_start = ifelse(age_group_id == 33, 95, age_start),
      age_end = ifelse(age_group_id == 33, 100, age_end),
      age_start = ifelse(age_group_id == 44, 100, age_start),
      age_end = ifelse(age_group_id == 44, 105, age_end) ) %>% 
    filter(age_end %in% seq(5, 105, by = 5))
  
  df_lt3 <- do.call("rbind", lapply(1:nrow(df_lt2), function(i) {
    ages_tmp <- df_lt2[i]$age_start:(df_lt2[i]$age_end-1)
    df_lt2[rep(i, length(ages_tmp)), ] %>%
      mutate(age_1yr = ages_tmp)
  }))
}
