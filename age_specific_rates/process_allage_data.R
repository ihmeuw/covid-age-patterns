process_allage_data <- function(df_locmeta, cfr_in_tmp1_path, 
                                cfr_extra_hospital_path, all_age_data_path) {
  ## Read input data
  df_cfr_in_tmp1 <- read.csv(cfr_in_tmp1_path, as.is = TRUE)
  df_cfr_extra_hospital <- read.csv(cfr_extra_hospital_path)
  
  ## Substitute custom locations pulled from extra hospital
  brazil_subnats <- df_locmeta %>%
    filter(parent_id == 135) %>%
    .[, "location_id"]
  other_custom_locs <- c(
    belgium = 76
  )
  all_custom_locs <- c(brazil_subnats, 135, other_custom_locs)
  df_cfr_in_custom <-  df_cfr_extra_hospital %>%
    filter(location_id %in% all_custom_locs)
  df_cfr2 <- df_cfr_in_tmp1 %>%
    filter(!location_id %in% all_custom_locs) %>%
    bind_rows(df_cfr_in_custom) %>%
    mutate(deaths2 = Deaths, hospitalizations2 = Hospitalizations, 
           cases2 = Confirmed, date = as.Date(Date))
  
  ## Save
  write.csv(df_cfr2, all_age_data_path, row.names = FALSE)
}

#####
# some all-age seroprevalence estimates don't have corresponding death estimates for their specific location_id, 
# which we need in order to get all-age IFR (deaths/infections) corresponding to that seroprevalence study
# -- for these locations, we need to aggregate up deaths from the subnational locations
interpolate_allage_time_series <- function(dat) {
  if (FALSE) {
    dat <- df_cfr_in
  }
  
  df <- dat %>%
    mutate(date = as_date(Date) )
  
  # -- create prediction frame including all locations and dates
  df_fullframe <- expand.grid(
    location_id = unique(df$location_id),
    date = as_date(min(df$date):max(df$date)) ) %>%
    arrange(location_id, date )
  
  df2 <- df_fullframe %>%
    left_join(df, by = c("location_id", "date")) %>%
    group_by(location_id) %>%
    mutate(
      # -- interpolate
      deaths2 = na.approx(Deaths, na.rm = FALSE),
      cases2 = na.approx(Confirmed, na.rm = FALSE),
      hospitalizations2 = na.approx(Hospitalizations, na.rm = FALSE),
      is_interpolated_deaths2 = ifelse(is.na(Deaths) & !is.na(deaths2), TRUE, FALSE),
      is_interpolated_cases2 = ifelse(is.na(Confirmed) & !is.na(cases2), TRUE, FALSE),
      is_interpolated_hospitalizations2 = ifelse(is.na(Hospitalizations) & !is.na(hospitalizations2), TRUE, FALSE),
      # -- get daily counts
      diff_deaths2 = c(NA, diff(deaths2, na.rm = FALSE)),
      diff_cases2 = c(NA, diff(deaths2, na.rm = FALSE)),
      diff_hospitalizations2 = c(NA, diff(hospitalizations2, na.rm = FALSE)),
      # -- identify whether the location has any negative daily values
      has_negative_diff_deaths = any(diff_deaths2 < 0, na.rm = TRUE),
      has_negative_diff_cases = any(diff_cases2 < 0, na.rm = TRUE),
      has_negative_diff_hospitalizations = any(diff_hospitalizations2 < 0, na.rm = TRUE) ) %>%
    as.data.frame(.)
  
  return(df2)
}

#####
# aggreagte sublocations for the CFR analyses (cases/deaths)
# and sum cases and deaths within a calendar month
# -- the warnings about -Inf are okay; those are location-months that didn't have any data
# -- note that we lag the cases part of the ratio by 8 days, because case status occurs before death
aggregate_sublocs <- function(dat, locs, location_metadata) {
  
  if (FALSE) {
    dat <- df_cfr
    locs <- sort(unique(df_sero$location_id))
    location_metadata <- df_locmeta
  }
  
  missing_locs <- c(locs[!locs %in% dat$location_id], 135)
  
  df_missing_locs <- df_locmeta %>%
    filter(location_id %in% missing_locs) %>%
    select(location_id, location_name)
  
  locid_has_all_sublocs <- sapply(missing_locs, function(parent_id_tmp) {
    if (FALSE) {
      parent_id_tmp <- 86
    }
    
    location_metadata_tmp <- location_metadata %>%
      filter(parent_id == parent_id_tmp)
    
    if (!all(location_metadata_tmp$location_id %in% dat$location_id)) {
      FALSE
    } else {
      TRUE
    }
  })
  names(locid_has_all_sublocs) <- missing_locs
  
  df_locmeta2 <- df_locmeta %>%
    filter(parent_id %in% missing_locs[locid_has_all_sublocs])
  
  df_deaths_subnat <- dat %>%
    filter(location_id %in% unique(df_locmeta2$location_id) ) %>%
    left_join(df_locmeta2[, c("location_id", "parent_id")], by = "location_id") %>%
    arrange(location_id, date) %>%
    group_by(parent_id, date) %>%
    summarize(
      n_subnats = n(),
      allage_deaths1 = sum(deaths2, na.rm = TRUE),
      allage_hospitalizations1 = sum(hospitalizations2, na.rm = TRUE),
      location_id_old = paste(location_id, collapse = ",") ) %>%
    group_by(parent_id) %>%
    mutate(
      allage_deaths2 = ifelse(n_subnats == max(n_subnats), allage_deaths1, NA),
      allage_deaths_from_subnat = na.approx(allage_deaths2, na.rm = FALSE),
      subnat_is_interpolated = is.na(allage_deaths2) & !is.na(allage_deaths_from_subnat),
      deaths_tplus8 = lead(allage_deaths_from_subnat, 8) ) %>%
    mutate(
      allage_hospitalizations2 = ifelse(n_subnats == max(n_subnats), allage_hospitalizations1, NA),
      allage_hosp_from_subnat = na.approx(allage_hospitalizations2, na.rm = FALSE),
      hosp_subnat_is_interpolated = is.na(allage_hospitalizations2) & !is.na(allage_hosp_from_subnat),
      hosp_tminus4 = lag(allage_hosp_from_subnat, 4) ) %>% # hosp admission happens 4 days before seropositivity
    as.data.frame(.)
  
  return(df_deaths_subnat)
  
}
