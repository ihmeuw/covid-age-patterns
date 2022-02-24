#
# ltc_analysis.R
#

library(dplyr)
library(data.table)
library(ggplot2)
library(ggsci)
library(lubridate)
library(RColorBrewer)

if (R.Version()[["major"]] == "3") {
  library(mrbrt001, lib.loc = "FILEPATH")
} else if (R.Version()[["major"]] == "4") {
  library(mrbrt002, lib.loc = "FILEPATH")
}

path1 <- "FILEPATH"
path2 <- "FILEPATH"


list.files(path1)
list.files(path2)

df_sero <- read.csv(file.path(path1, "seroprev_preds_5yr.csv"))
df_ifr1 <- read.csv(file.path(path1, "ifr_prepped_input_data_waningadjustment.csv"), as.is = TRUE)
df_ifr_preds1 <- readRDS(file.path(path1, "ifr_preds_1yr_global.rds") )
df_ifr_3set <- readRDS(file.path(path1, "ifr_preds_1yr_3set.rds"))
df_locset_in <- read.csv(file.path(path1, "ifr_preds_5yr_byloc.csv"))

df_pop_in <- as.data.frame(fread(file.path(
  "FILEPATH", "output_measures/population/all_populations.csv"
)))

df_pop <- df_pop_in %>%
  filter((age_group_years_end - age_group_years_start) == 5 | age_group_id == 235) %>%
  filter(sex_id == 3) %>%
  mutate(
    age_start = age_group_years_start,
    age_end = ifelse(age_group_years_end == 125, 99, age_group_years_end - 1)
  )

df_locset <- df_locset_in %>%
  select(location_id, custom_hierarchy) %>%
  filter(!duplicated(.)) %>%
  arrange(custom_hierarchy)

ltc_locs <- df_locset %>%
  filter(custom_hierarchy == "nursing_homes") %>%
  .[, "location_id"] %>%
  as.vector(.)


df_ifr2 <- df_ifr1 %>%
  filter(location_id %in% 76) %>%
  select(prev_logit, se_logit, age_start, age_end, data_filename) %>%
  mutate(age_mid = (age_start + age_end) / 2) %>%
  mutate(constant_se = 1)

dat5 <- MRData()
dat5$load_df(
  data = df_ifr2,  col_obs = "prev_logit", col_obs_se = "se_logit",
  col_covs = list("age_start", "age_end"), col_study_id = "data_filename" )

mod5 <- MRBRT(
  data = dat5,
  cov_models = list(
    LinearCovModel("intercept", use_re = TRUE), 
    LinearCovModel(
      alt_cov = list("age_start", "age_end"),
      use_spline = TRUE,
      spline_knots = array(c(0, 0.3, 0.6, 0.7, 0.8, 1)),
      spline_degree = 3L,
      spline_knots_type = 'domain',
      spline_r_linear = TRUE,
      spline_l_linear = TRUE
    )
  )
)

mod5$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)


df_pred5 <- data.frame(age_start = seq(0, 100, by = 1)) %>%
  mutate(age_end = age_start)

dat_pred5 <- MRData()

dat_pred5$load_df(
  data = df_pred5, 
  col_covs=list("age_start", "age_end")
)
samples5 <- mrbrt002::core$other_sampling$sample_simple_lme_beta(sample_size = 1000L, model = mod5)
draws5 <- mod5$create_draws(
  data = dat_pred5,
  beta_samples = samples5,
  gamma_samples = matrix(rep(0, nrow(samples5)), ncol = 1),
  random_study = FALSE )

df_pred5$pred5 <- mod5$predict(data = dat_pred5)
df_pred5$pred5_lo <- apply(draws5, 1, function(x) quantile(x, 0.025))
df_pred5$pred5_hi <- apply(draws5, 1, function(x) quantile(x, 0.975))

df_pred5$pred5_invlogit <- plogis(df_pred5$pred5)
df_pred5$pred5_lo_invlogit <- plogis(apply(draws5, 1, function(x) quantile(x, 0.025)))
df_pred5$pred5_hi_invlogit <- plogis(apply(draws5, 1, function(x) quantile(x, 0.975)))

df_pred5_2 <- df_pred5 %>%
  mutate(
    age_group = cut(age_start, breaks = seq(0, 100, by = 5), right = FALSE) ) %>%
  filter(!is.na(age_group)) %>%
  group_by(age_group) %>%
  summarize(
    age_start = min(age_start),
    sero_new_logit = mean(pred5)) %>%
  mutate(
    age_mid = age_start + 2.5) %>%
  as.data.frame(.)


par(mfrow = c(1,2))
with(df_ifr2, plot(
  x = age_mid, 
  y = prev_logit, 
  xlab = "Age group mid-point",
  ylab = "Seroprevalence, logit scale",
  main = "",
  type = "n", xlim = c(0, 100)))
grid()

for (f in unique(df_ifr2$data_filename)) {
  if (FALSE) {
    f <- "belgium_serobank_1.csv"
  }
  df_ifr_f <- df_ifr2 %>%
    filter(data_filename == f)
  cat(f, "\n")
  with(filter(df_ifr_f, data_filename == f), lines(age_mid, prev_logit))
}
with(df_pred5_2, lines(age_mid, sero_new_logit, lty = 1, lwd = 3, col = "darkorange"))


x <- 72.5
df_sero_tmp <- df_sero %>%
  mutate(age_mid = age_group_start + 2.5) %>%
  filter(age_mid == x)


diff_at_x <- df_pred5_2[df_pred5_2$age_mid == x, "sero_new_logit"] - df_sero_tmp[1, "logit_seroprev"]
with(df_sero, lines(age_group_start + 2.5, logit_seroprev + diff_at_x, lwd = 5, col = "skyblue"))

legend(
  x = "bottomleft",
  lty = c(1,1,1),
  lwd = c(1, 1, 5),
  col = c("black", "darkorange", "skyblue"),
  legend = c(
    "Serologic survey in Belgium",
    "Spline fit to Belgium data",
    "Global seroprevalence age pattern, scaled to match at age [70,74)" ),
  cex = 0.8
)


df1 <- df_pop %>%
  filter(location_id %in% ltc_locs) %>%
  left_join(df_ifr1 %>% select(location_id, location) %>% filter(!duplicated(.))) %>%
  left_join(df_sero %>% select(age_group_years_start = age_group_start, logit_sero_global = logit_seroprev)) %>%  # global seroprevalence
  left_join(df_pred5_2 %>% select(age_group_years_start = age_start, logit_sero_belgium = sero_new_logit)) %>%
  mutate(
    diff_in_agegroup_65 = max(ifelse(age_group_years_start == x - 2.5, logit_sero_belgium - logit_sero_global, NA), na.rm = TRUE),
    logit_sero_global2 = logit_sero_global + diff_in_agegroup_65,
    logit_sero_counterfactual = ifelse(age_group_years_start >= x - 2.5, logit_sero_global2, logit_sero_belgium),
    sero_belgium = plogis(logit_sero_belgium),
    sero_counterfactual = plogis(logit_sero_counterfactual)
  ) %>%
  mutate(age_mid = c(age_start + age_end) / 2)

with(df_ifr2, plot(
  x = age_mid, 
  y = prev_logit, 
  xlab = "Age group mid-point",
  ylab = "Seroprevalence, logit scale",
  main = "",
  type = "n", xlim = c(0, 100)))
grid()
with(df1 %>% filter(location_id == 76), 
     lines(
       x = age_mid, y = qlogis(sero_belgium), 
       col = "darkgreen", lwd = 1
     )
)

with(df1 %>% filter(location_id == 76), 
     lines(
       x = age_mid, y = qlogis(sero_counterfactual), 
       col = "darkgreen", lwd = 3, lty = 5
     )
)

legend(
  x = "bottomleft",
  col = "darkgreen",
  lwd = c(1,3),
  lty = c(1,5),
  legend = c("Baseline: LTC epidemic", "Counterfactual: No LTC epidemic"),
  cex = 0.8
)

df2 <- df1 %>%
  select(location_id, location, age_start = age_group_years_start, population, sero_belgium, sero_counterfactual) %>%
  left_join(
    df_locset_in %>% 
      filter(custom_hierarchy == "high_income") %>% 
      select(age_start = age_group_start, ifr_highincome = ifr) %>%
      filter(!duplicated(.))) %>%
  left_join(
    df_locset_in %>% 
      filter(custom_hierarchy == "nursing_homes") %>% 
      select(age_start = age_group_start, ifr_ltc = ifr) %>%
      filter(!duplicated(.))) %>%
  mutate(
    ifr_counterfactual = ifelse(age_start >= 65, ifr_highincome, ifr_ltc),
    infecs_orig = sero_belgium * population,
    deaths_orig = ifr_ltc * infecs_orig,
    infecs_counterfactual = sero_counterfactual * population,
    deaths_counterfactual = ifr_counterfactual * infecs_counterfactual,
    infecs_change_sero_only = sero_counterfactual * population,
    deaths_change_sero_only = ifr_ltc * infecs_change_sero_only
  )

df3 <- df2 %>%
  group_by(location_id, location) %>%
  summarize(
    infecs_orig = sum(infecs_orig),
    deaths_orig = sum(deaths_orig),
    infecs_counterfactual = sum(infecs_counterfactual),
    deaths_counterfactual = sum(deaths_counterfactual),
    infecs_change_sero_only = sum(infecs_change_sero_only),
    deaths_change_sero_only = sum(deaths_change_sero_only) ) %>%
  mutate(
    ifr_orig = deaths_orig / infecs_orig,
    ifr_counterfactual = deaths_counterfactual / infecs_counterfactual,
    ifr_change_sero_only = deaths_change_sero_only / infecs_change_sero_only,
    counter_over_orig = ifr_counterfactual / ifr_orig,
    counter_over_orig2 = ifr_change_sero_only / ifr_orig
  )

df4 <- df3 %>%
  as.data.frame(.) %>%
  select(location, counter_over_orig, counter_over_orig2) %>%
  mutate(percent_reduction = round((1 - counter_over_orig) * 100, digits = 1)) %>%
  mutate(percent_reduction2 = round((1 - counter_over_orig2) * 100, digits = 1)) %>%
  select(-counter_over_orig, -counter_over_orig2)
  
names(df4) <- c("Location", "% reduction", "% reduction changing seroprevalence only")

library(kableExtra)
df4 %>%
  kbl(booktabs = T) %>%
  kable_classic(full_width = FALSE) 



