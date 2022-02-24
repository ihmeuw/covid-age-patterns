prep_hospital_age_data <- function(input_dir, output_dir, all_age_data_path, le_out_path) {
  
  df_locmeta <- read.csv(file.path(input_dir, "locations/covariate_with_aggregates_hierarchy.csv"))
  df_pop_in <- read.csv(file.path(input_dir, "output_measures/population/all_populations.csv"))
  df_hosp <- read.csv(file.path(input_dir, "raw_formatted/hfr_age.csv"), as.is = TRUE)
  all_age_hosps <- read.csv(all_age_data_path, as.is = TRUE)
  df_lt3 <- read.csv(le_out_path, as.is = T)
  
  hosp_age_out_path <- file.path(output_dir, "df_hosp_byagedateloc_rawcumulative.csv")
  # Both sexes single year with 95+ population
  df_pop <- df_pop_in %>%
    filter((age_group_years_end - age_group_years_start) == 1 | age_group_id == 235) %>%
    filter(sex_id %in% 3) 
  
  # Merge on loc metadata and subset to both-sexes
  df_hosp_in <- df_hosp %>%
    left_join(
      df_locmeta %>% 
        select(
          location_id, location = location_name, region_name, super_region_name
        ),
      by = "location_id"
    ) %>% 
    filter(!is.na(hospitalizations)) %>%
    filter(sex_id == 3)
  
  # Check for incomplete age ranges
  age_range_df <- df_hosp_in %>% 
    group_by(data_filename) %>% 
    summarize(age_min = min(age_start), age_max = max(age_end))
  incomplete_age <- age_range_df %>% filter(age_min != 0 | age_max != 125)
  if(nrow(incomplete_age) > 0) {
    message("Incomplete age data:")
    print(incomplete_age)
  }
  
  # Look at the alignment between age groups for a location
  df_hosp_in %>% mutate(age_range = paste0(age_start, "-", age_end)) %>%
    group_by(location_id) %>%
    select(location_id, age_range) %>%
    unique(.) %>% as.data.frame(.)
  
  # Fix up the date
  df_hosp_in$date_string <- sapply(df_hosp_in$data_filename, function(x) {
    x2 <- strsplit(x, split = "_")[[1]]
    out_tmp <- gsub(".csv|-", "", x2[grepl("2020|2021", x2)])
    ifelse(length(out_tmp) == 1, out_tmp, NA)
  })
  
  # -- subset to locations with population denominators
  tmp_locs <- unique(df_hosp_in$location_id)
  exclude_locs_pop <- tmp_locs[!tmp_locs %in% df_pop_in$location_id]
  if(length(exclude_locs_pop) > 1) {
    cat(
      "Excluding locations with no population denominator:",
      paste(exclude_locs_pop, collapse = ", "), "\n"
    )
  }
  
  df_hosp1 <- df_hosp_in %>%
    filter(!is.na(hospitalizations)) %>%
    filter(!location_id %in% as.vector(exclude_locs_pop)) %>%
    mutate(
      date_tmp_ymd = ymd(date_string),
      date_tmp_dmy = dmy(date_string),
      date_tmp = ifelse(is.na(date_tmp_ymd), date_tmp_dmy, date_tmp_ymd),
      date = as_date(date_tmp),
      days_from_march1 = date - ymd("20200301")
    ) 
  
  # check for data points where age_end == age_start of the next group
  #    and subtract 1 from age_end if so;
  # -- if age_end is non-integer, round down so the whole 1-year age group is included
  df_hosp2 <- df_hosp1 %>%
    arrange(location_id, date, age_start) %>%
    group_by(location_id, date) %>%
    mutate(
      same_adjacent_age = age_end == lead(age_start, 1),
      age_end = ifelse(same_adjacent_age & row_number() != n(), age_end - 1, age_end),
      same_adjacent_age2 = age_end == lead(age_start, 1),
      age_end = floor(age_end)) %>%
    as.data.frame(.)
  
  ## Scale age-specific data to sum to total
  df_scalars <- df_hosp2 %>%
    group_by(location_id, date) %>%
    summarize(age_total = sum(hospitalizations)) %>%
    left_join(
      all_age_hosps %>%
        mutate(date = as.Date(date)) %>%
        select(location_id, date, all_age_total = hospitalizations2),
      by = c("location_id", "date")
    ) %>%
    mutate(scalar = all_age_total / age_total) %>%
    mutate(scalar = ifelse(is.na(scalar), 1, scalar))
  
  df_hosp3 <- df_hosp2 %>%
    left_join(
      df_scalars %>% select(location_id, date, scalar), 
      by = c("location_id", "date")
    ) %>%
    mutate(hospitalizations = hospitalizations * scalar)
  
  # get population for each unique location-age group combination
  # then merge it back on to the analytic dataset
  # -- also average age, in order to get reasonable upper bound on the terminal age group
  df_hosp_age_groups <- df_hosp3 %>%
    select(location_id, age_start, age_end) %>%
    filter(!duplicated(.))
  
  
  # -- adjust oldest age group upper bound 
  #    to be x + 2*e_x, where e_x is life expectancy at age x
  #

  age_group_dat <- lapply(1:nrow(df_hosp_age_groups), prep_age_group_data, df_hosp_age_groups, df_pop, df_lt3) 
  
  df_hosp_age_groups$pop <- sapply(age_group_dat, function(x) x$pop)
  df_hosp_age_groups$mean_age <- sapply(age_group_dat, function(x) x$mean_age)
  
  df_hosp4 <- df_hosp3 %>%
    left_join(df_hosp_age_groups, by = c("location_id", "age_start", "age_end")) %>%
    arrange(location_id, days_from_march1) %>%
    group_by(location_id) %>%
    mutate(
      age_end_bak = age_end,
      age_end = ifelse(
        age_end == max(age_end, na.rm = TRUE), 
        age_start + (mean_age - age_start)*2, age_end) ) %>%
    mutate(
      n_time_obs = length(unique(days_from_march1)) ) %>%
    as.data.frame(.) %>%
    filter(!duplicated(.[, c("location_id", "data_filename", "age_start", "age_end")])) %>%
    group_by(location_id, age_start) %>%
    mutate(
      prev_time_obs = lag(days_from_march1, 1), 
      prev_hosp = lag(hospitalizations, 1),
      n_days_interval = days_from_march1 - prev_time_obs,
      n_hosp_interval = hospitalizations - prev_hosp,
      loctime_id = paste0("locid_", location_id, "_time", prev_time_obs, "_", days_from_march1)
    )
  
  write.csv(df_hosp4, hosp_age_out_path)
  
}

get_avg_age <- function(yr_survival_prob) {
  survival_probs_tmp <- sapply(1:length(95:100), function(x) yr_survival_prob^x)
  weighted.mean(x = (95:100)+0.5, w = survival_probs_tmp)
}

prep_age_group_data <- function(i, df_hosp_age_groups, df_pop, df_lt3) {
  ages_tmp <- df_hosp_age_groups[i, "age_start"]:df_hosp_age_groups[i, "age_end"]
  df_pop_tmp <- df_pop %>% filter(location_id == df_hosp_age_groups[i, "location_id"])
  if (any(ages_tmp >= 95)) {
    df_probs <- data.frame(prob = seq(0.01, 0.99, by = 0.01))
    age95plus_pop <- df_pop_tmp[df_pop_tmp$age_group_years_start == 95, "population"]
    lt_tmp <- filter(df_lt3, location_id == df_pop_tmp[1, "location_id"] & age_1yr == 95 )[, "mean"]
    if (length(lt_tmp) == 0) {
      lt_tmp <- filter(df_lt3, location_id == 1 & age_1yr == 95 )[, "mean"]
    }
    df_probs$age_ex <- 95 + lt_tmp
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
  
  out_list <- list(
    pop = sum(df_pop_tmp2$population),
    mean_age = with(df_pop_tmp2, weighted.mean(x = age_group_years_start + 0.5, w = population, na.rm = TRUE))
  )
  
  return(out_list)
}