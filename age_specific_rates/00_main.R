## Setup --------------------------------------------------------------------------------------------------------------------------------------
Sys.umask("0002")
Sys.setenv(MKL_VERBOSE = 0)
suppressMessages({library(data.table); library(dplyr); library(parallel);
  library(lubridate);library(mrbrt002, lib.loc = "FILEPATH");
  library(yaml); library(zoo); library(RColorBrewer); library(ggsci); library(ggplot2)})
setDTthreads(1)

## Install ihme.covid (always use newest version)
tmpinstall <- system("mktemp -d --tmpdir=/tmp", intern = TRUE)
.libPaths(c(tmpinstall, .libPaths()))
devtools::install_github("ihmeuw/ihme.covid", ref = "main", upgrade = "never", quiet = T)

## Arguments
if (interactive()) {
  user <- Sys.info()["user"]
  code_dir <- getwd()
  OUTPUT_ROOT <- "FILEPATH"
  output_dir <- ihme.covid::get_output_dir(OUTPUT_ROOT, "today")
  output_dir <- paste0(output_dir, '/')
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

### Paths -----------------------------------------------------------------------------------------------------------------------------------------------

## Input directories
input_dir <- 'FILEPATH'
mort_age_dir <- 'FILEPATH'
cfr_in_tmp1_path <- 'FILEPATH'
cfr_extra_hospital_path <- 'FILEPATH'
all_age_data_path <- file.path(mort_age_dir, "prepped_cfr_file_with_aggregated_allage_seroprev_subnats.csv")
le_out_path <- file.path(output_dir, "life_expectancy_by_1yr_age.csv")
all_age_ifr_path <- "FILEPATH"
vacc_dir <- "FILEPATH"

# Diagnostics
diagnostics_dir <- file.path(output_dir, "diagnostics")
dir.create(diagnostics_dir)
loc_plot_path <- file.path(diagnostics_dir, "hospitalization_age_plot_preds_by_loc.pdf")

## Save YAML file
metadata <- yaml::yaml.load_file(file.path(mort_age_dir, "metadata.yaml"))
metadata_add <- list(
  mort_age_path = metadata$mortality_age_pattern_metadata$output_path)
metadata[["age_specific_rates_metadata"]] <- metadata_add
yaml::write_yaml(metadata, file.path(output_dir, "metadata.yaml"))

## Functions
source('FILEPATH/ifr_age.R')
source('FILEPATH/prep_life_expectancy.R')
source('FILEPATH/process_allage_data.R')
source('FILEPATH/save_age_sources.R')
source('FILEPATH/sero_age.R')
source('FILEPATH/hir_age.R')
source('FILEPATH/prep_hospital_age_data.R')
source('FILEPATH/fit_hospital_age.R')
source('FILEPATH/adjust_ifr.R')
source('FILEPATH/fit_ifr_age.R')
source('FILEPATH/plot_hosp_age.R')
source("FILEPATH/get_covariate_estimates.R")
source("FILEPATH/get_life_table.R")
source("FILEPATH/get_age_metadata.R")

logit <- function(p) log(p/(1-p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

## Tables
df_locmeta <- read.csv(file.path(input_dir, "locations/covariate_with_aggregates_hierarchy.csv"))

## Process all-age data
message("Prepping all-age data...")
process_allage_data(df_locmeta, cfr_in_tmp1_path, cfr_extra_hospital_path,
                    all_age_data_path)

## Prep life-expectancy
message("Prepping life-expectancy data...")
prep_life_expectancy(unique(df_locmeta$location_id), le_out_path)

## Prep IFR and seroprevalence by age
message("Prepping IFR by age...")
df5 <- prep_ifr_age( 
  input_dir, mort_age_dir, output_dir,
  predict_ifr_from_ensemble = F,
  optimize_knots = T, exclude_loc_ids = c(532, 143), use_goldstandard_only = F, all_age_data_path
)

## Apply IFR/seroprevalence adjustments
message('Adjusting seroprevalence data for waning immunity and all-age outliers...')
df5_adjusted <- adjust_ifr(all_age_ifr_path, paste0('ifr_prepped_input_data_v2021_', gsub('-', '_', Sys.Date()), '_waningadjustment.csv'))

# Fit IFR-age model
message('Fitting IFR by age...')
fit_ifr_age(df5_adjusted, '2021_10_12.01')

## Fit sero-age model 
message("Prepping seroprevalence by age...")
prep_sero_age(df5_adjusted, output_dir)

## Prep hospital age data
message("Prepping hospital age data...")
hospital_age_data <- prep_hospital_age_data(input_dir, output_dir, all_age_data_path, le_out_path)

## Fit hospital age model
message("Fitting hospitalization by age...")
fit_hospital_age(input_dir, output_dir)

## Prep IHR by age
message("Prepping IHR by age...")
hir_age(input_dir, output_dir)

## Plot hosp by age
message("Plotting hospitalizations by age... ")
plot_hosp_loc(output_dir, loc_plot_path)

## Save age-specific sources
save_age_sources(mort_age_dir, output_dir, df_locmeta)

# Mark as latest
OUTPUT_ROOT <- "FILEPATH"
system(paste("unlink", file.path(OUTPUT_ROOT, "latest")))
system(paste0("ln -s ", output_dir, " ", file.path(OUTPUT_ROOT, "latest")))

# Done
message(paste0(
  "\nDone!\n",
  "Output directory - ", output_dir, "\n"
))

# #mark as best
# system(paste("unlink", file.path(OUTPUT_ROOT, "best")))
# system(paste0("ln -s ", output_dir, " ", file.path(OUTPUT_ROOT, "best")))
