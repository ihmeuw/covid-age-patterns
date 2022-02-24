# -- get life expectancy for each observation, indexed by location and age
prep_life_expectancy <- function(ex_locs, le_out_path) {
  # 2017 is most recent in gbd_round_id = 5; round 6 requires a decomp step to be specified
  df_lt <- get_life_table(gbd_round_id = 5, location_id = ex_locs, year_id = 2017,
                          age_group_id = "all", sex_id = 3, life_table_parameter_id = 5)
  
  df_age_meta <- get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>%
    select(age_group_id, age_start = age_group_years_start, age_end = age_group_years_end)
  
  df_lt2 <- df_lt %>%
    left_join(df_age_meta, by = "age_group_id") %>%
    mutate( # todo: this is a stopgap; find age_group_set_id with [95,100) and [100,105)
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
  
  write.csv(df_lt3, le_out_path, row.names = FALSE)
}