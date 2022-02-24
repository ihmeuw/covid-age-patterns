plot_hosp_loc <- function(ouput_dir, loc_plot_path) {
  
  model_dir <- file.path(output_dir, "cascade_hosp_thetas3_gaussian")
  df3_cumulative_5 <- read.csv(file.path(model_dir, "input_data.csv"))
  
  preds_global <- read.csv(file.path(model_dir, "hosp_df_pred3.csv"), as.is = TRUE)
  preds_sr <- read.csv(file.path(model_dir, "hosp_cascade_preds_sr.csv"), as.is = TRUE)
  preds_loc <-  read.csv(file.path(model_dir, "hosp_cascade_preds_loc.csv"), as.is = TRUE)
  df_pred3 <- preds_global
  
  
  df_pop_in <- read.csv(file.path(input_dir, "output_measures/population/all_populations.csv"))
  df_pop <- df_pop_in %>%
    filter((age_group_years_end - age_group_years_start) == 5 | age_group_id == 235) %>%
    filter(sex_id == 3) %>%
    mutate(
      age_start = age_group_years_start,
      age_end = ifelse(age_group_years_end == 125, 99, age_group_years_end - 1)
    )
  
  logit <- function(p) log(p/(1-p))
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  pdf(loc_plot_path)
  
  for (loc_id in unique(df3_cumulative_5$location_id)) {
    
    try({
      
      if (FALSE) {
        loc_id <- 531 # washington d.c.
        loc_id <- 527 # california
        loc_id <- 130 # mexico
      }
      
      df_loc <- filter(df3_cumulative_5, location_id == loc_id) %>%
        mutate(
          age_group = paste0(age_start, "-", age_end),
          age_mid = (age_start + age_end) / 2,
          contains50 = age_start <= 50 & age_end >=50,
          age_range_yrs = age_end - age_start ) %>%
        mutate(
          logit_hosp_rate_lo = logit_hosp_rate - logit_hosp_rate_se * 1.96,
          logit_hosp_rate_hi = logit_hosp_rate + logit_hosp_rate_se * 1.96
        )
      
      if (loc_id == 531) {
        # remove duplicate observation for washington d.c.
        df_loc <- df_loc %>%
          filter(!(age_start == 50 & age_end == 64))
      }
      
      loc_name <- df_loc[1, "location_name"]
      sr_name_tmp <- unique(df_loc$super_region_name)
      loc_id_tmp <- as.character(loc_id)
      
      df_loc2 <- df_loc
      df_loc2$plot_var <- df_loc2$logit_mort
      y_lab <- "Cumulative hospitalizations, logit scale"
      main_txt <- loc_name
      
      preds_sr_tmp <- filter(preds_sr, super_region_name == df_loc2[1, "super_region_name"])
      preds_loc_tmp <- filter(preds_loc, location_id == loc_id)
      
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
      grid()
      abline(h = 0, col = "darkgray", lty = 2, lwd = 2)
      
      # global pred
      with(preds_global, lines(age_start, pred3, lwd = 3, lty = 3, col = "blue"))
      
      # super-region pred
      with(preds_sr_tmp, lines(age_start, pred, lwd = 2, lty = 5, col = "blue"))
      
      # location-specific pred
      with(preds_loc_tmp, lines(age_start, pred, lwd = 2, lty = 1))
      
      for (row in 1:nrow(df_loc2)) {
        with(df_loc2[row, ], lines(x = c(age_start, age_end), y = c(logit_hosp_rate, logit_hosp_rate)))
        with(df_loc2[row, ], lines(x = c(age_mid, age_mid), y = c(logit_hosp_rate_lo, logit_hosp_rate_hi)))
      }
      
      
      legend("bottomright",
             legend = c("Global fit", "Super-region fit", "Location fit"),
             lwd = c(3,2,2),
             lty = c(3,5, 1),
             col = c("blue", "blue", "black")
      )
 
    })
    
  }
  
  dev.off()
}