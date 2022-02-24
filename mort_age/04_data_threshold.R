determine_data_threshold <- function(mort_age_dir, data_rich, data_poor){
  #read in mort-age input 
  mort_input <- fread(paste0(mort_age_dir, '/df_deaths_byageloc_preppedcumulative.csv'))
  
  #read in cumulative deaths as of end of 2020
  infec_death_lag <- 24
  deaths_path <- paste0(ref_dir, "cumulative_deaths.csv")
  daily_deaths <- fread(deaths_path)
  daily_deaths[, date := as.Date(date)]
  daily_deaths <- daily_deaths[date < (as.Date("2021-01-01") + infec_death_lag) &
                                 date >= "2020-01-01"]
  daily_deaths_long <- melt.data.table(daily_deaths[date == "2020-12-31"], id.vars=c('location_id', 'date'))
  cum_deaths <- daily_deaths_long[, .(deaths=mean(value)), by='location_id']
  cum_deaths <- merge(cum_deaths, hierarchy[, .(location_id, location_name, level)], by='location_id', all.x=T)
  cum_deaths <- cum_deaths[level>=3 & location_id %in% unique(mort_input$location_id)]
  firstq <- as.integer(quantile(cum_deaths$deaths, 0.25))
  ld <- cum_deaths[deaths<=firstq]$location_id
  
  mort_input_n_bylocage <- mort_input[, .(num_data_points=.N), by=.(location_id, location, age_mid)]
  mort_input_n_bylocage[0<=age_mid & age_mid<=25, age_range:="0-25"]
  mort_input_n_bylocage[25<age_mid & age_mid<=50, age_range:="26-50"]
  mort_input_n_bylocage[50<age_mid & age_mid<=75, age_range:="51-75"]
  mort_input_n_bylocage[75<age_mid & age_mid<=125, age_range:="76-100"]
  mort_input_n_bylocage <- mort_input_n_bylocage[, .(num_data_points=sum(num_data_points)), by=.(location_id, location, age_range)]
  mort_input_n_bylocage <- dcast.data.table(mort_input_n_bylocage, location_id+location~age_range, value.var='num_data_points')
  dp <- mort_input_n_bylocage[is.na(`26-50`) | is.na(`51-75`) | is.na(`76-100`) | (`0-25`==1 & `26-50`==1 & `51-75`==1 & `76-100`==1)]$location_id
  
  #intersection
  dp_ld <- intersect(dp, ld)
  dp_or_ld <- c(dp, ld)
  dp_or_ld <- dp_or_ld[!duplicated(dp_or_ld)]
  dp_or_ld <- data.table(location_id=dp_or_ld)
  
  dp_or_ld <- merge(dp_or_ld, hierarchy[, .(location_id, location_name)], by='location_id')
  dp_or_ld[location_id %in% dp & !location_id %in% dp_ld, reason:='mort-age input data poor']
  dp_or_ld[location_id %in% ld & !location_id %in% dp_ld, reason:='low total covid deaths']
  dp_or_ld[location_id %in% dp_ld, reason:='mort-age input data poor & low total covid deaths']
  
  write.csv(dp_or_ld, paste0(mort_age_dir, '/diagnostics/data_poor_low_death_loc_ids.csv'), row.names=F)
  
  #create one file with correct mort predictions --------------------------------------------------------------------------------------------------------------------
  
  #### FIVE YEAR
  #final file format to check against
  format <- fread(paste0(mort_age_dir, '/cascade_', data_rich, '/mortality_agepattern_preds_byloc_5yr.csv'))
  dim(format)
  length(unique(format$location_id))
  
  #splice together
  data_rich_preds <- fread(paste0(mort_age_dir, '/cascade_', data_rich, '/mortality_agepattern_preds_byloc_5yr.csv'))[!location_id %in% dp_or_ld$location_id]
  data_poor_preds <- fread(paste0(mort_age_dir, '/cascade_', data_poor, '/mortality_agepattern_preds_byloc_5yr.csv'))[location_id %in% dp_or_ld$location_id]
  final_preds <- rbind(data_rich_preds, data_poor_preds)
  
  #check
  dim(final_preds)
  length(unique(final_preds$location_id))
  format[!location_id %in% unique(final_preds$location_id)]
  
  #save
  write.csv(final_preds, paste0(mort_age_dir, '/mortality_agepattern_preds_byloc_5yr.csv'), row.names=F)
  
  ### ONE YEAR
  #final file format to check against
  format <- fread(paste0(mort_age_dir, '/cascade_', data_rich, '/mortality_agepattern_preds_byloc_1yr.csv'))
  dim(format)
  length(unique(format$location_id))
  
  #splice together
  data_rich_preds <- fread(paste0(mort_age_dir, '/cascade_', data_rich, '/mortality_agepattern_preds_byloc_1yr.csv'))[!location_id %in% dp_or_ld$location_id]
  data_poor_preds <- fread(paste0(mort_age_dir, '/cascade_', data_poor, '/mortality_agepattern_preds_byloc_1yr.csv'))[location_id %in% dp_or_ld$location_id]
  final_preds <- rbind(data_rich_preds, data_poor_preds)
  
  #check
  dim(final_preds)
  length(unique(final_preds$location_id))
  format[!location_id %in% unique(final_preds$location_id)]
  
  #save
  write.csv(final_preds, paste0(mort_age_dir, '/mortality_agepattern_preds_byloc_1yr.csv'), row.names=F)
}










