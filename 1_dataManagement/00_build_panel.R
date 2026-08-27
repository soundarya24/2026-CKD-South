#################################################################
##                       Soundarya S                            ##
##                       8-26                                   ##
##                       CKD-South India analysis               ##
##                       IMRG-CKD-CC                            ##
#################################################################

# ckd-cc
# Combines the 5-state 30-year temperature predictors (1994-2023) with GBD 2023
# CKD outcomes for two causes:
        #   - 589 = all-CKD aggregate
        #   - 593 = CKD due to other & unspecified causes (CKD-U, heat-relevant)
# Source GBD data: _rawData/ckd-all-measures/  (two combined pulls, Number+Rate)
#   - deaths_inc_prev.csv : Deaths(1), Prevalence(5), Incidence(6)
#   - daly_yll_yld.csv    : DALYs(2), YLDs(3), YLLs(4)
#   - all-cause DALYs from _rawData/DALY-YLL-YLD/all-cause-DALY.csv (cause 294)
#
# Author: Soundarya S
# Version: August 2026

# Packages
pacman::p_load(tidyverse, here, janitor)

#===============================================================================

# Code

south_states <- c("Andhra Pradesh", "Karnataka", "Kerala",
                  "Tamil Nadu","Telangana" )
state_map <- c(
        andhrapradesh = "Andhra Pradesh",
        karnataka     = "Karnataka",
        kerala        = "Kerala",
        tamilnadu     = "Tamil Nadu",
        telangana     = "Telangana"
)
#study_states <- unname(state_map)
yr_min <- 1994
yr_max <- 2023

# 12 temperature predictors: {mean,median,mode,max,min,range} x {2.5th,97.5th}
# 'range' (daily max-min) is reported as TEMPERATURE VARIABILITY.
pred_cols <- c(
        "mean_2_5","mean_97_5","median_2_5","median_97_5","mode_2_5","mode_97_5",
        "max_2_5","max_97_5","min_2_5","min_97_5","range_2_5","range_97_5"
)


pred_dir   <- here("_rawData", "temperature")
pred_files <- list.files(pred_dir, pattern = "_percentiles\\.csv$", full.names = TRUE)

temperature <- pred_files |>
        set_names(basename) |>
        map_dfr(read_csv, show_col_types = FALSE, .id = "file") |>
        rename(state_name = state) |>
        mutate(state_name = case_when(
                state_name == "andhrapradesh" ~ "Andhra Pradesh",
                state_name == "karnataka" ~ "Karnataka",
                state_name == "telangana" ~ "Telangana",
                state_name == "tamilnadu" ~ "Tamil Nadu",
                state_name == "kerala" ~ "Kerala"

        )) |>
        filter(year >= yr_min, year <= yr_max) |>
        select(state_name, year, all_of(pred_cols)) |>
        drop_na(state_name) |>
        arrange(state_name, year)








# ---- 2. GBD loader (combined files) -----------------------------------------


SEX_BOTH   <- 3L    # Both sexes
AGE_ALL    <- 22L   # All ages
METRIC_NUM <- 1L    # Number (counts)   -> economics
METRIC_RATE<- 3L    # Rate per 100,000  -> Aim 1 regression
CKD_CAUSES <- c(589L, 593L)            # all-CKD; CKD-other/unspecified (CKD-U)
gbd_deaths_path <- here("_rawData", "ckd-all-measures", "deaths_inc_prev.csv")
gbd_daly_path   <- here("_rawData", "ckd-all-measures", "daly_yll_yld.csv")

read_gbd_combined <- function(path, measure_ids, metric_id, cause_ids = CKD_CAUSES) {
  read_csv(path, show_col_types = FALSE) |>
    filter(
      measure_id %in% measure_ids,
      cause_id   %in% cause_ids,
      metric_id  == !!metric_id,
      sex_id     == SEX_BOTH,
      age_id     == AGE_ALL,
      location_name %in% south_states,
      year >= yr_min, year <= yr_max
    ) |>
    transmute(state_name = location_name, year, cause_id, measure_name, val)
}

# Aim 1 outcomes: RATE per 100,000 for Deaths(1), Prevalence(5), Incidence(6)
outcomes_rate <- read_gbd_combined(gbd_deaths_path, c(1L, 5L, 6L), METRIC_RATE) |>
  mutate(
    cause   = recode(as.character(cause_id), "589" = "allCKD", "593" = "CKDu"),
    measure = recode(measure_name,
                     "Deaths" = "deaths", "Incidence" = "incidence",
                     "Prevalence" = "prevalence")
  )

outcomes_wide <- outcomes_rate |>
  select(state_name, year, cause, measure, val) |>
  pivot_wider(names_from = c(measure, cause), values_from = val,
              names_glue = "{measure}_{cause}")

# ---- 3. Aim 1 panel: predictors + RATE outcomes -----------------------------

panel <- temperature |>
  left_join(outcomes_wide, by = c("state_name", "year")) |>
  arrange(state_name, year)



# ---- 5. Save + diagnostics --------------------------------------------------

out_dir <- here("2_derivedData")
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(panel,      file.path(out_dir, "ckd_temp_panel.rds"))

cat("\n================ PANEL DIAGNOSTICS ================\n")
cat("Aim 1 panel rows:", nrow(panel),
    "| states:", n_distinct(panel$state_name),
    "| years:", min(panel$year), "-", max(panel$year), "\n")
cat("Aim 1 outcomes are RATES per 100,000. Expected rows: 5 x 30 = 150\n\n")

cat("Outcome columns present:\n")
print(grep("^(deaths|incidence|prevalence)_", names(panel), value = TRUE))

cat("\nMissing values per outcome column:\n")
panel |>
  select(matches("^(deaths|incidence|prevalence)_")) |>
  summarise(across(everything(), ~ sum(is.na(.)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") |>
  print(n = Inf)

cat("\nOutcome ranges (sanity — these are per-100k rates):\n")
panel |>
  select(matches("^(deaths|incidence|prevalence)_")) |>
  summarise(across(everything(), list(min = ~min(.,na.rm=TRUE),
                                      max = ~max(.,na.rm=TRUE)))) |>
  pivot_longer(everything()) |> print(n = Inf)


cat("\nSaved: _rawData/derived/ckd_temp_panel.rds (rates) ")
cat("==================================================\n")

