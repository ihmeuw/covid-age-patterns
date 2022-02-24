## Setup
Sys.umask("0002")
Sys.setenv(MKL_VERBOSE = 0)
suppressMessages({library(data.table); library(dplyr); library(parallel);
  library(lubridate);library(mrbrt002, lib.loc = "FILEPATH");
  library(yaml); library(RColorBrewer); library(zoo)})
setDTthreads(1)

## Install ihme.covid (always use newest version)
tmpinstall <- system("mktemp -d --tmpdir=/tmp", intern = TRUE)
.libPaths(c(tmpinstall, .libPaths()))
devtools::install_github("ihmeuw/ihme.covid", upgrade = "never", quiet = T)

## Arguments
if (interactive()) {
  user <- Sys.info()["user"]
  code_dir <- getwd()
  OUTPUT_ROOT <- "FILEPATH"
  output_dir <- ihme.covid::get_output_dir(OUTPUT_ROOT, "today")
  lsvid <- 771
} else {
  code_dir <- ihme.covid::get_script_dir()
  parser <- argparse::ArgumentParser()
  parser$add_argument("--lsvid", type = "integer", required = TRUE, help = "location set version id to use")
  parser$add_argument("--outputs-directory", required = TRUE, help = "Directory in which to write outputs")
  args <- parser$parse_args()
  lsvid <- args$lsvid
  output_dir <- args$outputs_directory
}

#settings
exclude_loc_ids <- c(143)

#vaccination roll-out information
vacc_dir <- "FILEPATH"
vacc_dt <- fread(paste0(vacc_dir, "slow_scenario_vaccine_coverage.csv"))
min_vacc_date <- vacc_dt[cumulative_all_effective > 0][, .(min_vacc_date = min(date)), by = .(location_id)]

## Paths
# In 
input_dir <- "FILEPATH" 
pop_path <- file.path(input_dir, "output_measures/population/all_populations.csv")
nchs_dir <- file.path(input_dir, "NCHS")
nchs_filename <- file.path(nchs_dir, '/NCHS_20210102_age_bins.csv') #latest pre-vac NCHS, manually selected
deaths_path <- file.path(input_dir, "raw_formatted/cfr_age_combined.csv")
deaths_add_path <- file.path(input_dir, "raw_formatted/cfr_age_additional.csv")
locmeta_path <- file.path(input_dir, "locations/covariate_with_aggregates_hierarchy.csv")
sero_sources_path <- file.path(input_dir, "serological_age_sources.csv")
ref_dir <- 'FILEPATH'

# Out
deaths_bydate_path <- file.path(output_dir, "df_deaths_date.csv")
deaths_byagedateloc_path <- file.path(output_dir, "df_deaths_byagedateloc_rawcumulative.csv")
deaths_byageloc_path <- file.path(output_dir, "df_deaths_byageloc_preppedcumulative.csv")

# Diagnostics
diagnostics_dir <- file.path(output_dir, "diagnostics")
dir.create(diagnostics_dir)
loc_plot_path <- file.path(diagnostics_dir, "mortality_age_plot_preds_by_loc.pdf")
sr_plot_path <- file.path(diagnostics_dir, "mortality_plot_superregion_preds.pdf")

## Save YAML file
metadata <- yaml::yaml.load_file(file.path(input_dir, "metadata.yaml"))
metadata_add <- list(output_path = output_dir, input_path = metadata$output_path, vacc_input_path = vacc_dir)
metadata[["mortality_age_pattern_metadata"]] <- metadata_add
yaml::write_yaml(metadata, file.path(output_dir, "metadata.yaml"))

## Functions
source(file.path("FILEPATH/get_location_metadata.R"))
source("FILEPATH/get_life_table.R")
source("FILEPATH/get_age_metadata.R")
source('FILEPATH/01_prep_data.R')
source('FILEPATH/02_model.R')
source('FILEPATH/03_predict.R')
source('FILEPATH/04_data_threshold.R')
source('FILEPATH/05_plot_mort_age.R')
source('FILEPATH/06_plot_mort_age_pub.R')

logit <- function(p) log(p/(1-p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

## Tables
df_locmeta <- fread("FILEPATH/gbd_analysis_hierarchy.csv")
hierarchy <- get_location_metadata(location_set_id = 111, location_set_version_id = lsvid, release_id=9)

## Prep and save data 
message("Prepping data...")
model_data <- prep_data(
  df_locmeta, pop_path, nchs_dir, nchs_filename, deaths_path, deaths_add_path,
  deaths_bydate_path, deaths_byagedateloc_path, deaths_byageloc_path, 
  sero_sources_path, exclude_loc_ids
)

## Fit global and cascade models with thetas 3,0.3
message("Fitting global and cascade models...")
thetas <- c(3,0.3)
version <- 'thetas_3_03'

# Model outpaths change depending on theta parameters
model_out_dir <- paste0(output_dir, '/cascade_', version)
dir.create(model_out_dir)
mort_knot_path <- file.path(model_out_dir, "mort_optimal_knot_locs.RDS")
global_mort_preds_path <- file.path(model_out_dir, "global_mort_preds_1yr.csv")
mr_agepattern_preds_byloc_5yr <- file.path(model_out_dir, "mortality_agepattern_preds_byloc_5yr.csv")
mr_agepattern_preds_byloc_1yr <- file.path(model_out_dir, "mortality_agepattern_preds_byloc_1yr.csv")
mr_agepattern_preds_global <- file.path(model_out_dir, "mortality_agepattern_preds_global.RDS")
mr_agepattern_preds_sr <- file.path(model_out_dir, "mortality_agepattern_preds_sr.RDS")
mr_agepattern_preds_loc <- file.path(model_out_dir, "mortality_agepattern_preds_loc.RDS")

#fit global
mod_global <- fit_global(model_data, mort_knot_path)


#save model object
py_save_object(
  object = mod_global,
  filename = file.path(model_out_dir, "/global_model.pkl"),
  pickle = "dill"
)

#fit cascade
cascade_fit <- fit_cascade(mod_global, model_data, thetas, model_out_dir)

## Make and save predictions
message("Making predictions...")
preds <- predict_mort_age(
  mod_global, cascade_fit, model_data, df_locmeta, global_mort_preds_path,
  mr_agepattern_preds_byloc_5yr, mr_agepattern_preds_byloc_1yr,
  mr_agepattern_preds_global, mr_agepattern_preds_sr, mr_agepattern_preds_loc
)

## Fit global and cascade models with thetas 3,3
message("Fitting global and cascade models...")
thetas <- c(3,3)
version <- 'thetas_3_3'

# Model outpaths change depending on theta parameters
model_out_dir <- paste0(output_dir, '/cascade_', version)
dir.create(model_out_dir)
mort_knot_path <- file.path(model_out_dir, "mort_optimal_knot_locs.RDS")
global_mort_preds_path <- file.path(model_out_dir, "global_mort_preds_1yr.csv")
mr_agepattern_preds_byloc_5yr <- file.path(model_out_dir, "mortality_agepattern_preds_byloc_5yr.csv")
mr_agepattern_preds_byloc_1yr <- file.path(model_out_dir, "mortality_agepattern_preds_byloc_1yr.csv")
mr_agepattern_preds_global <- file.path(model_out_dir, "mortality_agepattern_preds_global.RDS")
mr_agepattern_preds_sr <- file.path(model_out_dir, "mortality_agepattern_preds_sr.RDS")
mr_agepattern_preds_loc <- file.path(model_out_dir, "mortality_agepattern_preds_loc.RDS")

#fit global
mod_global <- fit_global(model_data, mort_knot_path)

#save model object
py_save_object(
  object = mod_global,
  filename = file.path(model_out_dir, "/global_model.pkl"),
  pickle = "dill"
)

#fit cascade
cascade_fit <- fit_cascade(mod_global, model_data, thetas, model_out_dir)

## Make and save predictions
message("Making predictions...")
preds <- predict_mort_age(
  mod_global, cascade_fit, model_data, df_locmeta, global_mort_preds_path,
  mr_agepattern_preds_byloc_5yr, mr_agepattern_preds_byloc_1yr,
  mr_agepattern_preds_global, mr_agepattern_preds_sr, mr_agepattern_preds_loc
)

## Determine threshold of which model to use and create combine final file
message('Determining data threshold and creating final preds file...')
determine_data_threshold(output_dir, data_rich='thetas_3_3', data_poor='thetas_3_03')

## Make plots
message('Plotting preds...')
output_dir <- paste0(output_dir, '/')
plot_mort_preds(model_data, output_dir, data_rich='thetas_3_3', data_poor='thetas_3_03')

## Make plots for pub
message('Plotting preds for publication...')
plot_mort_preds_for_pub(model_data, output_dir, data_rich='thetas_3_3', data_poor='thetas_3_03')

# Mark as latest
system(paste("unlink", file.path(OUTPUT_ROOT, "latest")))
system(paste0("ln -s ", output_dir, " ", file.path(OUTPUT_ROOT, "latest")))

# # 
# message(paste0(
#   "\nDone!\n",
#   "Output directory -\n",
#   output_dir, "\n",
#   "Plots -\n",
#   loc_plot_path, "\n",
#   sr_plot_path 
# ))