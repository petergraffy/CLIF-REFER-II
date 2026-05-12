# ================================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) | PI: Peter Graffy
# Trajectory Cohort Builder (CLIF) - ARF-wide + Respiratory Support Sensitivity Cohorts
#
# Primary cohort   : Adults (>=18) with first evidence of ARF during ICU (t0 = ARF onset)
# Sensitivity      : IMV and ADVANCED respiratory support cohorts retained for comparison
#
# Outputs (written to output/run_[SITE]_[DATE]/):
#   cohort_arf72.csv                    : 1 row / hospitalization (t0 = first ARF evidence)
#   cohort_primary_imv72.csv            : sensitivity cohort (t0 = IMV start)
#   cohort_secondary_adv72.csv          : sensitivity cohort (t0 = first advanced start)
#   exclusion_*.csv                     : exclusions + first failing reason
#   flow_*.csv                          : flow counts
#
# Notes:
# - ARF-wide primary cohort avoids selecting only patients who receive IMV.
# - Evidence source flags are retained for transparent sensitivity analyses.
# - No ICU LOS >= 24h rule (allow early deaths/extubations; handle as censoring downstream).
# - Optional: require minimum respiratory_support density post-t0 (set MIN_RS_HOURS > 0).
# ================================================================================================

suppressPackageStartupMessages({
  library(fst)
  library(here)
  library(tidyverse)
  library(arrow)
  library(dplyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(glue)
  library(scales)
})

# ---------- Project config ----------
source("utils/config.R")
stopifnot(exists("config"))

repo        <- config$repo
site_name   <- config$site_name
tables_path <- config$tables_path
file_type   <- config$file_type

stopifnot(!is.null(repo), nzchar(repo))
stopifnot(!is.null(tables_path), nzchar(tables_path))

cat("Site Name:", site_name, "\n")
cat("Tables Path:", tables_path, "\n")
cat("File Type:", file_type, "\n")

tables_path <- normalizePath(tables_path, mustWork = TRUE)

# ---------- Cohort parameters ----------
START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")

ADULT_AGE_YEARS <- 18
ICU_TO_T0_MAX_H <- Inf         # include ARF that develops any time during ICU; set 24 for early-onset sensitivity
TRAJ_HOURS      <- 72          # trajectory window length [t0, t0 + 72h]
PRE_T0_BUFFER_H <- 0           # keep 0 by default; set >0 if you want a short pre-t0 baseline window
MIN_RS_HOURS    <- 0           # avoid excluding lab/vital-defined ARF before device escalation
# set to 0 to disable

FIO2_JOIN_H     <- 1
HYPERCAP_JOIN_H <- 2
ROOM_AIR_FIO2   <- 0.21
SPO2_ARF_CUTOFF <- 90
PAO2_ARF_CUTOFF <- 60
PF_ARF_CUTOFF   <- 300
PCO2_ARF_CUTOFF <- 45
PH_ARF_CUTOFF   <- 7.35

# Device categories (per your screenshot)
DEVICE_IMV      <- "IMV"
DEVICE_ADVANCED <- c("IMV", "NIPPV", "CPAP", "High Flow NC")

# ---------- File discovery / loading ----------
exts <- strsplit(file_type, "[/|,; ]+")[[1]]
exts <- exts[nzchar(exts)]
if (length(exts) == 0) exts <- c("csv","parquet","fst")
ext_pat <- paste0("\\.(", paste(unique(exts), collapse = "|"), ")$")

all_files <- list.files(
  path = tables_path,
  pattern = ext_pat,
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0) {
  stop("No files with extensions {", paste(exts, collapse = ", "), "} were found under: ", tables_path)
}

bn <- basename(all_files)
looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)
base_no_ext <- tools::file_path_sans_ext(tolower(bn))
base_norm <- ifelse(looks_clif, base_no_ext, paste0("clif_", base_no_ext))
found_map <- stats::setNames(all_files, base_norm)

required_raw <- c("patient","hospitalization","adt","respiratory_support","vitals","labs")
required_files <- paste0("clif_", required_raw)
missing <- setdiff(required_files, names(found_map))
if (length(missing) > 0) {
  cat("Detected CLIF-like files:\n"); print(sort(unique(names(found_map))))
  stop("Missing required tables: ", paste(missing, collapse = ", "))
}

clif_paths <- found_map[required_files]

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

clif_tables <- lapply(clif_paths, read_any)
names(clif_tables) <- required_files
cat("Loaded tables: ", paste(names(clif_tables), collapse = ", "), "\n")

# ---------- Helpers ----------
safe_posix <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  # robust parse; avoids adding a hard dependency on fasttime
  suppressWarnings(as.POSIXct(x, tz = "UTC"))
}

safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

get_min <- function(tbl_name, cols) {
  nm <- paste0("clif_", tbl_name)
  stopifnot(!is.null(clif_tables[[nm]]))
  out <- clif_tables[[nm]] %>% rename_with(tolower)
  cols <- tolower(cols)
  missing_cols <- setdiff(cols, names(out))
  for (missing_col in missing_cols) out[[missing_col]] <- NA
  out %>% dplyr::select(all_of(cols))
}

keep_any <- function(df, cols) dplyr::select(df, any_of(intersect(cols, names(df))))

# ---------- Minimal tables ----------
patient <- get_min("patient",
                   c("patient_id","birth_date","sex_category","race_category","ethnicity_category")) %>%
  mutate(birth_date = safe_date(birth_date))

hospitalization <- get_min("hospitalization",
                           c("patient_id","hospitalization_id","admission_dttm","discharge_dttm","age_at_admission",
                             "zipcode_nine_digit","zipcode_five_digit","census_tract","county_code")) %>%
  mutate(
    admission_dttm = safe_posix(admission_dttm),
    discharge_dttm = safe_posix(discharge_dttm)
  )

adt <- get_min("adt",
               c("hospitalization_id","in_dttm","out_dttm","location_category","location_type")) %>%
  mutate(
    in_dttm  = safe_posix(in_dttm),
    out_dttm = safe_posix(out_dttm)
  )

resp_support <- get_min("respiratory_support",
                        c("hospitalization_id","recorded_dttm","device_category","device_name",
                          "mode_category","mode_name","fio2_set","lpm_set","flow_rate_set",
                          "peep_set","peep_obs","tidal_volume_set","tidal_volume_obs",
                          "resp_rate_set","resp_rate_obs","plateau_pressure_obs",
                          "peak_inspiratory_pressure_set","peak_inspiratory_pressure_obs",
                          "mean_airway_pressure_obs","minute_vent_obs",
                          "pressure_control_set","pressure_support_set",
                          "tracheostomy","artificial_airway")) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    device_category = as.character(device_category),
    mode_category   = as.character(mode_category)
  )

vitals <- get_min("vitals",
                  c("hospitalization_id","recorded_dttm","vital_category","vital_name",
                    "vital_value","meas_site_name")) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    vital_category = as.character(vital_category),
    vital_name = as.character(vital_name),
    vital_value = suppressWarnings(as.numeric(vital_value))
  )

labs <- get_min("labs",
                c("hospitalization_id","lab_result_dttm","lab_collect_dttm",
                  "lab_category","lab_name","lab_value","lab_value_numeric",
                  "reference_unit","lab_specimen_category")) %>%
  mutate(
    lab_result_dttm = safe_posix(lab_result_dttm),
    lab_collect_dttm = safe_posix(lab_collect_dttm),
    lab_category = as.character(lab_category),
    lab_name = as.character(lab_name),
    lab_value_num = suppressWarnings(as.numeric(coalesce(lab_value_numeric, lab_value)))
  )

# ---------- ICU bounds (first ICU in) ----------
icu_segments <- adt %>%
  mutate(is_icu = str_detect(tolower(coalesce(location_category, "")), "icu")) %>%
  filter(is_icu)

icu_bounds <- icu_segments %>%
  group_by(hospitalization_id) %>%
  summarize(
    first_icu_in = suppressWarnings(min(in_dttm, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(out_dttm, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
    last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours"))
  )

# ---------- Base candidate hospitalizations ----------
base <- hospitalization %>%
  inner_join(icu_bounds, by = "hospitalization_id") %>%
  filter(!is.na(first_icu_in),
         first_icu_in >= START_DATE,
         first_icu_in <= END_DATE) %>%
  left_join(patient %>% select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
            by = "patient_id") %>%
  mutate(
    age_years = coalesce(
      suppressWarnings(as.numeric(age_at_admission)),
      ifelse(!is.na(birth_date),
             as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
    )
  ) %>%
  mutate(
    has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
    adult    = !is.na(age_years) & age_years >= ADULT_AGE_YEARS,
    has_geo  = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code)
  )

# ---------- t0 identification ----------
normalize_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  case_when(
    is.na(x) ~ NA_real_,
    x > 1.5 ~ x / 100,
    TRUE ~ x
  )
}

in_icu_window <- function(events, time_col) {
  events %>%
    inner_join(base %>% select(hospitalization_id, first_icu_in, last_icu_out),
               by = "hospitalization_id") %>%
    filter(.data[[time_col]] >= first_icu_in,
           .data[[time_col]] <= last_icu_out)
}

rs_icu <- resp_support %>%
  mutate(
    device_category = as.character(device_category),
    mode_category = as.character(mode_category),
    fio2 = normalize_fio2(fio2_set)
  ) %>%
  in_icu_window("recorded_dttm")

# First IMV time retained as sensitivity cohort.
t0_imv <- rs_icu %>%
  filter(device_category == DEVICE_IMV) %>%
  group_by(hospitalization_id) %>%
  summarize(
    t0 = min(recorded_dttm, na.rm = TRUE),
    arf_evidence = "imv",
    .groups = "drop"
  ) %>%
  mutate(cohort = "primary_imv")

# First advanced support time retained as sensitivity cohort.
t0_adv <- rs_icu %>%
  filter(device_category %in% DEVICE_ADVANCED) %>%
  group_by(hospitalization_id) %>%
  summarize(
    t0 = min(recorded_dttm, na.rm = TRUE),
    arf_evidence = "advanced_support",
    .groups = "drop"
  ) %>%
  mutate(cohort = "secondary_adv")

advanced_events <- rs_icu %>%
  filter(device_category %in% DEVICE_ADVANCED) %>%
  transmute(hospitalization_id, t0 = recorded_dttm, arf_evidence = "advanced_support")

spo2_events <- vitals %>%
  filter(vital_category == "spo2", !is.na(vital_value), vital_value < SPO2_ARF_CUTOFF) %>%
  in_icu_window("recorded_dttm") %>%
  transmute(hospitalization_id, t0 = recorded_dttm, arf_evidence = "spo2_lt_90")

lab_events <- labs %>%
  filter(!is.na(lab_result_dttm), !is.na(lab_value_num)) %>%
  in_icu_window("lab_result_dttm")

fio2_ref <- rs_icu %>%
  filter(!is.na(fio2)) %>%
  transmute(hospitalization_id, fio2_time = recorded_dttm, fio2)

pair_nearest <- function(left_df, left_time, right_df, right_time, max_gap_h) {
  if (nrow(left_df) == 0 || nrow(right_df) == 0) {
    out <- bind_cols(left_df[0, , drop = FALSE], right_df[0, setdiff(names(right_df), "hospitalization_id"), drop = FALSE])
    out$left_time_keep <- as.POSIXct(character())
    out$gap_h <- numeric()
    return(out)
  }
  left_dt <- as.data.table(left_df)
  right_dt <- as.data.table(right_df)
  left_dt[, left_time_keep := get(left_time)]
  setkeyv(left_dt, c("hospitalization_id", left_time))
  setkeyv(right_dt, c("hospitalization_id", right_time))
  out <- right_dt[left_dt, roll = "nearest", nomatch = 0L]
  out[, gap_h := abs(as.numeric(difftime(left_time_keep, get(right_time), units = "hours")))]
  as_tibble(out[gap_h <= max_gap_h])
}

pao2_labs <- lab_events %>%
  filter(lab_category == "po2_arterial") %>%
  transmute(hospitalization_id, lab_time = lab_result_dttm, pao2 = lab_value_num)

pao2_fio2 <- pair_nearest(pao2_labs, "lab_time", fio2_ref, "fio2_time", FIO2_JOIN_H)

pf_events <- pao2_fio2 %>%
  filter(!is.na(pao2), !is.na(fio2), fio2 > 0, pao2 / fio2 <= PF_ARF_CUTOFF) %>%
  transmute(hospitalization_id, t0 = left_time_keep, arf_evidence = "pf_ratio_le_300")

room_air_pao2_events <- pao2_fio2 %>%
  filter(!is.na(pao2), !is.na(fio2), fio2 <= ROOM_AIR_FIO2 + 1e-6,
         pao2 <= PAO2_ARF_CUTOFF) %>%
  transmute(hospitalization_id, t0 = left_time_keep, arf_evidence = "room_air_pao2_le_60")

pco2_labs <- lab_events %>%
  filter(lab_category == "pco2_arterial") %>%
  transmute(hospitalization_id, pco2_time = lab_result_dttm, pco2 = lab_value_num)

ph_labs <- lab_events %>%
  filter(lab_category == "ph_arterial") %>%
  transmute(hospitalization_id, ph_time = lab_result_dttm, ph = lab_value_num)

hyper_pairs <- pair_nearest(pco2_labs, "pco2_time", ph_labs, "ph_time", HYPERCAP_JOIN_H)

hyper_events <- hyper_pairs %>%
  filter(!is.na(pco2), !is.na(ph), pco2 >= PCO2_ARF_CUTOFF, ph < PH_ARF_CUTOFF) %>%
  transmute(hospitalization_id, t0 = left_time_keep, arf_evidence = "hypercapnic_acidosis")

t0_arf_events <- bind_rows(
  advanced_events,
  spo2_events,
  pf_events,
  room_air_pao2_events,
  hyper_events
) %>%
  filter(!is.na(t0))

t0_arf <- t0_arf_events %>%
  group_by(hospitalization_id) %>%
  arrange(t0, .by_group = TRUE) %>%
  summarize(
    t0 = first(t0),
    arf_evidence = first(arf_evidence),
    arf_evidence_all = paste(sort(unique(arf_evidence)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(cohort = "arf72")

# ---------- Shared cohort builder ----------
build_traj_cohort <- function(base_df, t0_df, cohort_tag) {
  
  # join t0 and compute timing rule relative to ICU admit
  df <- base_df %>%
    left_join(t0_df %>% select(hospitalization_id, t0, any_of(c("arf_evidence", "arf_evidence_all"))),
              by = "hospitalization_id") %>%
    mutate(
      t0_within_icu24h = as.numeric(difftime(t0, first_icu_in, units = "hours")) <= ICU_TO_T0_MAX_H &
        as.numeric(difftime(t0, first_icu_in, units = "hours")) >= -6  # allow small negative due to timestamp quirks
    )
  
  # define trajectory window
  win <- df %>%
    transmute(
      hospitalization_id,
      t0,
      win_start = t0 - dhours(PRE_T0_BUFFER_H),
      win_end   = t0 + dhours(TRAJ_HOURS)
    )
  
  # respiratory_support records inside window (for density QC + downstream features)
  rs_win <- resp_support %>%
    inner_join(win, by = "hospitalization_id") %>%
    filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
    select(hospitalization_id, recorded_dttm, device_category, mode_category,
           fio2_set, peep_set, peep_obs, tidal_volume_set, tidal_volume_obs,
           lpm_set, flow_rate_set, pressure_control_set, pressure_support_set,
           plateau_pressure_obs, peak_inspiratory_pressure_set,
           peak_inspiratory_pressure_obs, mean_airway_pressure_obs,
           minute_vent_obs, resp_rate_set, resp_rate_obs,
           tracheostomy, artificial_airway)
  
  # density: number of distinct hours with any RS record in [t0, t0+TRAJ_HOURS]
  rs_density <- rs_win %>%
    mutate(hour_bin = floor_date(recorded_dttm, "hour")) %>%
    group_by(hospitalization_id) %>%
    summarize(rs_hours_observed = n_distinct(hour_bin), .groups = "drop") %>%
    mutate(meets_rs_density = ifelse(MIN_RS_HOURS <= 0, TRUE, rs_hours_observed >= MIN_RS_HOURS))
  
  # inclusion flags (trajectory-focused)
  flags <- df %>%
    left_join(rs_density, by = "hospitalization_id") %>%
    mutate(
      rs_hours_observed = coalesce(rs_hours_observed, 0L),
      meets_rs_density  = coalesce(meets_rs_density, FALSE),
      include =
        coalesce(adult, FALSE) &
        coalesce(has_demo, FALSE) &
        coalesce(has_geo, FALSE) &
        coalesce(t0_within_icu24h, FALSE) &
        (if (MIN_RS_HOURS <= 0) TRUE else coalesce(meets_rs_density, FALSE))
    )
  
  # final cohort (1 row per hospitalization)
  cohort <- flags %>%
    filter(include) %>%
    transmute(
      cohort = cohort_tag,
      patient_id, hospitalization_id,
      admission_dttm, discharge_dttm,
      first_icu_in, last_icu_out, icu_los_hours,
      age_years, sex_category, race_category, ethnicity_category,
      census_tract, county_code, zipcode_five_digit, zipcode_nine_digit,
      arf_evidence = coalesce(arf_evidence, cohort_tag),
      arf_evidence_all = coalesce(arf_evidence_all, arf_evidence),
      t0,
      data_window_start = t0 - dhours(PRE_T0_BUFFER_H),
      data_window_end   = t0 + dhours(TRAJ_HOURS),
      rs_hours_observed
    )
  
  # exclusions
  exclusions <- flags %>%
    filter(!include) %>%
    mutate(reason = case_when(
      !coalesce(adult, FALSE) ~ "Under 18 or missing age",
      !coalesce(has_demo, FALSE) ~ "Missing demographics",
      !coalesce(has_geo, FALSE)  ~ "Missing geo code",
      is.na(t0) ~ "No qualifying ARF/support t0 found",
      !coalesce(t0_within_icu24h, FALSE) ~ "t0 outside configured ICU onset window",
      (MIN_RS_HOURS > 0) & !coalesce(meets_rs_density, FALSE) ~ glue("Insufficient respiratory_support density (<{MIN_RS_HOURS} hourly bins)"),
      TRUE ~ "Other"
    )) %>%
    select(patient_id, hospitalization_id, reason)
  
  # flow table
  cand_n <- nrow(df)
  
  step1 <- df %>% filter(coalesce(adult, FALSE))
  step2 <- step1 %>% filter(coalesce(has_demo, FALSE))
  step3 <- step2 %>% filter(coalesce(has_geo, FALSE))
  step4 <- step3 %>% filter(coalesce(t0_within_icu24h, FALSE))
  step5 <- if (MIN_RS_HOURS <= 0) step4 else step4 %>%
    left_join(rs_density, by = "hospitalization_id") %>%
    filter(coalesce(meets_rs_density, FALSE))
  
  flow <- tibble(
    step = c(
      "ICU candidates (date range)",
      glue(">= {ADULT_AGE_YEARS} years"),
      "Demographics present",
      "Geography present",
      if (is.infinite(ICU_TO_T0_MAX_H)) "t0 during ICU stay"
      else glue("t0 within +{ICU_TO_T0_MAX_H}h of first ICU admit"),
      if (MIN_RS_HOURS <= 0) "Respiratory_support density rule (disabled)"
      else glue(">= {MIN_RS_HOURS} distinct hourly bins in trajectory window")
    ),
    remaining = c(cand_n, nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5))
  ) %>%
    mutate(excluded_at_step = lag(remaining, default = remaining[1]) - remaining)
  
  list(
    cohort = cohort,
    exclusions = exclusions,
    flow = flow,
    rs_win = rs_win
  )
}

# ---------- Build ARF-wide primary cohort + support sensitivity cohorts ----------
res_arf       <- build_traj_cohort(base, t0_arf, "arf72")
res_primary   <- build_traj_cohort(base, t0_imv, "primary_imv72")
res_secondary <- build_traj_cohort(base, t0_adv, "secondary_adv72")

cohort_arf       <- res_arf$cohort
cohort_primary   <- res_primary$cohort
cohort_secondary <- res_secondary$cohort

excluded_arf       <- res_arf$exclusions
excluded_primary   <- res_primary$exclusions
excluded_secondary <- res_secondary$exclusions

flow_arf       <- res_arf$flow
flow_primary   <- res_primary$flow
flow_secondary <- res_secondary$flow

# ---------- Quick counts ----------
cat("\nCohort selection summary:\n")
cat("  ICU candidates:          ", nrow(base), "\n", sep = "")
cat("  ARF-wide included:       ", nrow(cohort_arf), "\n", sep = "")
cat("  IMV sensitivity included:", nrow(cohort_primary), "\n", sep = "")
cat("  Secondary ADV included:  ", nrow(cohort_secondary), "\n", sep = "")
cat("  ARF-wide excluded:       ", nrow(excluded_arf), "\n", sep = "")
cat("  Primary excluded:        ", nrow(excluded_primary), "\n", sep = "")
cat("  Secondary excluded:      ", nrow(excluded_secondary), "\n", sep = "")

# ---------- Save outputs ----------
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "SITE" else x
}
SITE_NAME   <- sanitize_tag(site_name)
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")

make_name <- function(result_name, ext = "csv") {
  paste0(result_name, "_", SITE_NAME, "_", SYSTEM_DATE, ".", ext)
}

out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# cohorts
write_csv(cohort_arf,       file.path(out_dir, make_name("cohort_arf72")))
write_csv(cohort_primary,   file.path(out_dir, make_name("cohort_primary_imv72")))
write_csv(cohort_secondary, file.path(out_dir, make_name("cohort_secondary_adv72")))

# exclusions
write_csv(excluded_arf,       file.path(out_dir, make_name("exclusion_arf72")))
write_csv(excluded_primary,   file.path(out_dir, make_name("exclusion_primary_imv72")))
write_csv(excluded_secondary, file.path(out_dir, make_name("exclusion_secondary_adv72")))

# flow
write_csv(flow_arf,       file.path(out_dir, make_name("flow_arf72")))
write_csv(flow_primary,   file.path(out_dir, make_name("flow_primary_imv72")))
write_csv(flow_secondary, file.path(out_dir, make_name("flow_secondary_adv72")))

cat("\nSaved outputs to: ", out_dir, "\n", sep = "")
cat("  - ", make_name("cohort_arf72"), "\n", sep = "")
cat("  - ", make_name("cohort_primary_imv72"), "\n", sep = "")
cat("  - ", make_name("cohort_secondary_adv72"), "\n", sep = "")
cat("  - ", make_name("exclusion_arf72"), "\n", sep = "")
cat("  - ", make_name("exclusion_primary_imv72"), "\n", sep = "")
cat("  - ", make_name("exclusion_secondary_adv72"), "\n", sep = "")
cat("  - ", make_name("flow_arf72"), "\n", sep = "")
cat("  - ", make_name("flow_primary_imv72"), "\n", sep = "")
cat("  - ", make_name("flow_secondary_adv72"), "\n", sep = "")

# ---------- Optional: keep minimal objects ----------
keep_vars <- c("clif_tables", "cohort_arf", "cohort_primary", "cohort_secondary", "repo", "out_dir")
rm(list = setdiff(ls(envir = .GlobalEnv), keep_vars), envir = .GlobalEnv)

message("\nCohort identification complete.")

