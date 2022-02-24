plot_mort_preds_for_pub <- function(model_data, mort_age_dir, data_rich, data_poor){
  
  output_dir <- paste0(mort_age_dir, 'diagnostics/')
  input_dir <- mort_age_dir
  deaths_byageloc_path <- paste0(input_dir, 'df_deaths_byageloc_preppedcumulative.csv')
  
  #data rich
  df3_cumulative_5 <- model_data
  preds_global <- readRDS(paste0(input_dir, "cascade_", data_rich, "/mortality_agepattern_preds_global.RDS"))
  preds_sr <- readRDS(paste0(input_dir, "cascade_", data_rich, "/mortality_agepattern_preds_sr.RDS"))
  preds_loc <-  readRDS(paste0(input_dir, "cascade_", data_rich, "/mortality_agepattern_preds_loc.RDS"))
  
  #data poor
  preds_loc_3 <- readRDS(paste0(input_dir, "cascade_", data_poor, "/mortality_agepattern_preds_loc.RDS"))
  preds_sr_3 <- readRDS(paste0(input_dir, "cascade_", data_poor, "/mortality_agepattern_preds_sr.RDS"))
  preds_global_3 <- readRDS(paste0(input_dir, "cascade_", data_poor, "/mortality_agepattern_preds_global.RDS"))
  
  #pull in data threshold criteria to visualize
  dp_or_ld <- fread(paste0(output_dir, 'data_poor_low_death_loc_ids.csv'))
  df3_cumulative_5$data_threshold <- ifelse(df3_cumulative_5$location_id %in% dp_or_ld$location_id, 
                                            paste0('not highly impacted (', data_poor, ')'), paste0('highly impacted (', data_rich, ')'))
  
  #fxns
  logit <- function(p) log(p/(1-p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # pdf(loc_plot_path)
  pdf(paste0(output_dir, 'mortality_preds_plot_for_manuscript.pdf'))
  
  for (loc_id in unique(df3_cumulative_5$location_id)) {

      df_loc <- filter(df3_cumulative_5, location_id == loc_id) %>%
        mutate(
          age_group = paste0(age_start, "-", age_end),
          age_mid = (age_start + age_end) / 2,
          contains50 = age_start <= 50 & age_end >=50,
          age_range_yrs = age_end - age_start ) %>%
        mutate(
          logit_mort_lo = logit_mort - logit_mort_se * 1.96,
          logit_mort_hi = logit_mort + logit_mort_se * 1.96
        )
      
      loc_name <- df_loc[1, "location"]
      sr_name_tmp <- unique(df_loc$super_region_name)
      loc_id_tmp <- as.character(loc_id)
      
      df_loc2 <- df_loc
      df_loc2$plot_var <- df_loc2$logit_mort
      y_lab <- "Deaths / population, logit scale"
      main_txt <- loc_name
      
      preds_sr_tmp <- filter(preds_sr, super_region_name == sr_name_tmp)
      
      if (unique(df_loc2$data_threshold)=='highly impacted (thetas_3_3)'){
        preds_loc_tmp <- filter(preds_loc, location_id == loc_id)
      } else if (unique(df_loc2$data_threshold)=='not highly impacted (thetas_3_03)'){
        preds_loc_tmp <- filter(preds_loc_3, location_id == loc_id)
      }
      
      yvals <- do.call("c", list(
        df_loc2$logit_mort_lo, df_loc2$logit_mort_hi, 
        preds_global$pred3, preds_sr_tmp$pred, preds_loc_tmp$pred
      ))
      
      plot(
        x = c(df_loc2$age_start, df_loc$age_end),
        y = c(df_loc2$plot_var, df_loc2$plot_var),
        type = "n",
        main = main_txt,
        xlab = "Age",
        ylab = y_lab,
        xlim = c(0, 100),
        ylim = c(min(yvals), max(yvals))
      )
      title(sub=paste0('Cumulative reported deaths as of: ', unique(df_loc2$date)),adj=1,cex.sub=0.7,font=4)
      grid()
      abline(h = 0, col = "darkgray", lty = 2, lwd = 2)
      
      #global pred
      with(preds_global, lines(age_start, pred3, lwd = 3, lty = 3, col = "black"))
      
      # super-region pred
      with(preds_sr_tmp, lines(age_start, pred, lwd = 2, lty = 5, col = "black"))
      
      if (unique(df_loc2$data_threshold)=='highly impacted (thetas_3_3)'){
        # location-specific pred
        with(preds_loc_tmp, lines(age_start, pred, lwd = 2, lty = 1, col = 'blue'))
      } else if (unique(df_loc2$data_threshold)=='not highly impacted (thetas_3_03)'){
        # location-specific pred
        with(preds_loc_tmp, lines(age_start, pred, lwd = 2, lty = 1, col = 'darkorange'))
      }

      for (row in 1:nrow(df_loc2)) {
        with(df_loc2[row, ], lines(x = c(age_start, age_end), y = c(logit_mort, logit_mort)))
        with(df_loc2[row, ], lines(x = c(age_mid, age_mid), y = c(logit_mort_lo, logit_mort_hi)))
      }
      
      if (unique(df_loc2$data_threshold)=='highly impacted (thetas_3_3)'){
        # location-specific pred
        legend("bottomright",
               legend = c("Global fit", "Super-region fit", "Location fit (Prior multiplier = [3,3])"),
               lwd = c(3,2,2,2),
               lty = c(3,5,1,1),
               col = c("black", "black", "blue")
        )
      } else if (unique(df_loc2$data_threshold)=='not highly impacted (thetas_3_03)'){
        legend("bottomright",
               legend = c("Global fit", "Super-region fit", "Location fit (Prior multiplier = [3,0.3])"),
               lwd = c(3,2,2,2),
               lty = c(3,5,1,1),
               col = c("black", "black", "darkorange")
        )
      }
      

      
    }
  
  dev.off()

}