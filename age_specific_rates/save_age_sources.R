prep_death_sources <- function(dir) {
  deaths_byageloc_path <- file.path(dir, "df_deaths_byageloc_preppedcumulative.csv")
  dt <- fread(deaths_byageloc_path)
  data_sources_dt <- data.table( 
    unique(dt[, .(location_id, data_filename, days_from_march1)])
  )
  data_sources_dt <- data_sources_dt[!is.na(data_filename)]
  setnames(data_sources_dt, "data_filename", "file")
  data_sources_dt[, date := days_from_march1 + as.Date("2020-03-01")]
  data_sources_dt[, days_from_march1 := NULL]
  data_sources_dt[, measure := "deaths"]
  
  return(data_sources_dt[, .(measure, location_id, date, file)])
}

prep_sero_sources <- function(dir) {
  sero_path <- file.path(dir, 'ifr_prepped_input_data.csv')
  dt <- fread(sero_path)
  data_sources_dt <- data.table( 
    unique(dt[, .(location_id, data_filename, date)])
  )
  data_sources_dt <- data_sources_dt[!is.na(data_filename)]
  setnames(data_sources_dt, "data_filename", "file")
  data_sources_dt[, date := as.Date(date)]
  data_sources_dt[, measure := "seroprevalence"]
  
  return(data_sources_dt[, .(measure, location_id, date, file)])  
}

save_age_sources <- function(death_dir, sero_dir, df_locmeta) {
  mort_dt <- prep_death_sources(death_dir)
  sero_dt <- prep_sero_sources(sero_dir)
  dt <- rbind(mort_dt, sero_dt)
  dt <- merge(dt, df_locmeta %>% select(location_id, location_name), by = "location_id", all.x = T)
  dt <- dt[order(measure, location_id, date), .(measure, location_id, location_name, date, file)]
  write.csv(dt, file.path(sero_dir, "age_sources.csv"), row.names = F)
  
  outliers <- fread(paste0(output_dir, 'all_age_ifr_outliers.csv'))
  outliers2 <- fread(paste0(output_dir, 'geodiscordant_ifr_outliers.csv'))
  
  age_sources_outliers_removed <- dt[!file %in% unique(outliers$data_filename) & !file %in% unique(outliers2$data_filename)]
  write.csv(age_sources_outliers_removed, paste0(output_dir, 'age_sources_outliers_removed.csv'), row.names=F)
  
}