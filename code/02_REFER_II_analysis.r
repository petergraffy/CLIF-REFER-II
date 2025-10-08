# ===============================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) Index
# PI: Peter Graffy (graffy@uchicago.edu)
# Run this script after running: 01_REFER_cohort_identification.R
# Purpose: Build trajectories per patient, 
# ===============================================================================================


# Core
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)
library(purrr)
library(data.table)
library(dtwclust)   # DTW-based time-series clustering
library(nnet)       # multinomial regression for class ~ exposures
library(ggplot2)
library(fasttime)

exposome_dir <- "exposome"

pm25_yr <- read_csv(file.path(exposome_dir, "pm25_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, pm25 = pm25_mean)

no2_yr  <- read_csv(file.path(exposome_dir, "no2_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, no2 = no2_mean)

svi_yr  <- read_csv(file.path(exposome_dir, "svi_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, svi = svi_overall)

# Daymet annuals (optional, if you’d like climate covariates in the membership models)
tmax_yr <- read_csv(file.path(exposome_dir, "daymet_tmax_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, tmax = tmax_mean)

tmin_yr <- read_csv(file.path(exposome_dir, "daymet_tmin_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, tmin = tmin_mean)

vp_yr   <- read_csv(file.path(exposome_dir, "daymet_vp_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, vp = vp_mean)

prcp_yr <- read_csv(file.path(exposome_dir, "daymet_prcp_county_year.csv")) %>%
  rename(county_fips = GEOID, year = year, prcp = prcp_mean)

no2_mo <- read_csv(file.path(exposome_dir, "no2_county_month.csv")) %>%
  rename(county_fips = county_fips, year = year, no2 = no2)

cohort <- cohort_min %>%
  mutate(
    admit_dt = as_date(admission_dttm),
    admit_year = year(admit_dt),
    county_fips = county_code
  )

# ---- Annual 1–3y rolling exposure (example) ----
# Long stack annual exposome
expo_annual <- pm25_yr %>%
  full_join(no2_yr,  by = c("county_fips","year")) %>%
  full_join(svi_yr,  by = c("county_fips","year")) %>%
  full_join(tmax_yr, by = c("county_fips","year")) %>%
  full_join(tmin_yr, by = c("county_fips","year")) %>%
  full_join(vp_yr,   by = c("county_fips","year")) %>%
  full_join(prcp_yr, by = c("county_fips","year"))

# Helper: attach past 3y (y-1, y-2, y-3)
attach_annual_lags <- function(df, expo_tbl) {
  out <- df %>%
    # Build rows for y-1, y-2, y-3
    tidyr::expand_grid(lag = 1:3) %>%
    mutate(expo_year = admit_year - lag) %>%
    left_join(expo_tbl, by = c("county_fips" = "county_fips", "expo_year" = "year")) %>%
    group_by(patient_id, hospitalization_id) %>%
    summarise(
      pm25_mean_3y = mean(pm25, na.rm = TRUE),
      no2_mean_3y  = mean(no2,  na.rm = TRUE),
      svi_last     = dplyr::last(na.omit(svi)),        # take most recent available
      tmax_mean_3y = mean(tmax, na.rm = TRUE),
      tmin_mean_3y = mean(tmin, na.rm = TRUE),
      vp_mean_3y   = mean(vp,   na.rm = TRUE),
      prcp_sum_3y  = mean(prcp, na.rm = TRUE),
      .groups = "drop"
    )
  out
}

cohort_expo_annual <- attach_annual_lags(cohort, expo_annual)

# ---- Monthly NO2: 12-month pre-admission mean (if present) ----
# Path to the new monthly NO2 panel you created
no2_path <- file.path(exposome_dir, "no2_county_month.csv")
has_no2_monthly <- file.exists(no2_path)

if (has_no2_monthly) {
  no2win <- no2_mo %>%
    mutate(ym = make_date(year, month, 1)) %>%
    dplyr::select(county_fips, ym, no2 = no2)
  
  cohort_no2_12mo <- cohort %>%
    transmute(patient_id, hospitalization_id, county_fips, start = floor_date(admit_dt %m-% months(12), "month"),
              end = floor_date(admit_dt %m-% months(1), "month")) %>%
    rowwise() %>%
    mutate(ym_seq = list(seq.Date(start, end, by = "month"))) %>%
    unnest(ym_seq) %>%
    left_join(no2win, by = c("county_fips","ym_seq" = "ym")) %>%
    group_by(patient_id, hospitalization_id) %>%
    summarise(no2_mean_12mo = mean(no2, na.rm = TRUE), .groups = "drop")
  
  cohort_expo <- cohort %>%
    left_join(cohort_expo_annual, by = c("patient_id","hospitalization_id")) %>%
    left_join(cohort_no2_12mo,    by = c("patient_id","hospitalization_id"))
} else {
  cohort_expo <- cohort %>%
    left_join(cohort_expo_annual, by = c("patient_id","hospitalization_id"))
}


# Convert to POSIXct just in case
cohort_stays <- cohort %>%
  mutate(
    first_icu_in  = as_datetime(first_icu_in),
    last_icu_out  = as_datetime(last_icu_out)
  ) %>%
  filter(!is.na(first_icu_in), !is.na(last_icu_out), last_icu_out > first_icu_in)

# 3A) Day-level grid per hospitalization
# Safer day-grid builder
make_day_grid <- function(df, cap_days = NA_integer_) {
  df %>%
    transmute(
      hospitalization_id,
      start = as_date(floor_date(first_icu_in,  unit = "day")),
      end   = as_date(floor_date(last_icu_out,  unit = "day"))
    ) %>%
    mutate(end = if_else(end < start, start, end)) %>%  # guard weird clocks
    rowwise() %>%
    mutate(day = list(seq.Date(start, end, by = "day"))) %>%
    unnest(day) %>%
    group_by(hospitalization_id) %>%
    mutate(icu_day = as.integer(day - min(day)) + 1L) %>%
    { if (!is.na(cap_days)) filter(., icu_day <= cap_days) else . } %>%
    ungroup()
}

icu_days <- make_day_grid(cohort_stays)

# --- Helpers ---------------------------------------------------------------

# Find a table inside `clif_tables` by name pattern(s)
find_tbl <- function(tbl_list, patterns) {
  stopifnot(is.list(tbl_list))
  nms <- names(tbl_list)
  idx <- which(Reduce(`|`, lapply(patterns, function(p) str_detect(tolower(nms), tolower(p)))))
  if (length(idx) == 0) stop("Could not find a table matching patterns: ", paste(patterns, collapse = " | "))
  tbl_list[[idx[1]]]
}

# Pick the first existing column name from candidates
pick_name <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm) == 0) stop("None of the candidate columns found: ", paste(candidates, collapse = ", "))
  nm[1]
}

# Coerce to POSIXct safely (accepts character/POSIXct)
as_posix <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(x)
  suppressWarnings(lubridate::ymd_hms(x, tz = tz))
}

# --- Pull & standardize needed CLIF tables --------------------------------

# Respiratory support table (names may vary: "respiratory_support", "respiratory", "ventilation", etc.)
rs_raw <- find_tbl(clif_tables, c("respiratory_support", "respiratory", "vent"))

# Medication continuous infusion table (vasopressors)
med_cont_raw <- find_tbl(clif_tables, c("medication_admin_continuous", "med_admin_continuous", "infusion"))

# Vitals/labs table(s) for SpO2; if you have a separate vitals table use that term too
vitals_raw <- find_tbl(clif_tables, c("vitals", "vital", "flowsheet"))

# --- Build RESP summary per day -------------------------------------------

# helper: map various encodings -> 0/1 (integer), NA if unknown
as01 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- dplyr::case_when(
    x %in% c("1","true","t","yes","y") ~ 1L,
    x %in% c("0","false","f","no","n","") ~ 0L,
    suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) > 0) ~ 1L,
    suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) == 0) ~ 0L,
    TRUE ~ NA_integer_
  )
  out
}

# prefer recorded_dttm; fallback to recorded_time
rs_time_col <- if ("recorded_dttm" %in% names(rs_raw)) "recorded_dttm" else "recorded_time"

resp <- rs_raw %>%
  transmute(
    hospitalization_id,
    time  = suppressWarnings(lubridate::ymd_hms(.data[[rs_time_col]], tz = "UTC")),
    day   = as_date(time),
    device_name_lc      = tolower(as.character(device_name)),
    device_category_lc  = tolower(as.character(device_category)),
    mode_name_lc        = tolower(as.character(mode_name)),
    mode_category_lc    = tolower(as.character(mode_category)),
    artificial_airway   = as01(artificial_airway),
    tracheostomy        = as01(tracheostomy),
    fio2                = suppressWarnings(as.numeric(fio2_set)),
    peep_set            = suppressWarnings(as.numeric(peep_set)),
    peep_obs            = suppressWarnings(as.numeric(peep_obs))
  ) %>%
  # normalize FiO2 to % if needed (0–1 -> *100)
  mutate(
    fio2 = dplyr::case_when(
      is.na(fio2) ~ NA_real_,
      fio2 <= 1   ~ fio2 * 100,     # fraction -> percent
      TRUE        ~ fio2            # already percent
    )
  ) %>%
  # robust invasive signal
  mutate(
    invasive_flag = (artificial_airway == 1L) | (tracheostomy == 1L) |
      str_detect(device_category_lc, "invasive|mechanical|vent") |
      str_detect(mode_category_lc,   "invasive|ac|vc|pcv|simv")  |
      str_detect(device_name_lc,     "vent|invasive")            |
      str_detect(mode_name_lc,       "invasive|ac|vc|pcv|simv")
  ) %>%
  group_by(hospitalization_id, day) %>%
  summarise(
    any_imv    = as.integer(any(invasive_flag, na.rm = TRUE)),
    max_fio2   = suppressWarnings(max(fio2, na.rm = TRUE)),
    median_peep = suppressWarnings(median(coalesce(peep_obs, peep_set), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    max_fio2   = ifelse(is.infinite(max_fio2), NA_real_, max_fio2),
    median_peep = ifelse(is.infinite(median_peep), NA_real_, median_peep)
  )

# --- Build VASOPRESSOR summary per day ------------------------------------

med_id_col    <- pick_name(med_cont_raw, c("hospitalization_id", "encounter_id", "stay_id"))
med_start_col <- pick_name(med_cont_raw, c("admin_dttm", "start_time", "start_datetime", "start_dt"))
med_name_col  <- pick_name(med_cont_raw, c("medication_name", "drug_name", "med_name", "infusion_name"))

vaso_names <- c("norepinephrine","epinephrine","vasopressin","phenylephrine","dopamine")
vaso <- med_cont_raw %>%
  transmute(
    hospitalization_id = .data[[med_id_col]],
    start_dt = as_posix(.data[[med_start_col]]),
    med_name = tolower(as.character(.data[[med_name_col]]))
  ) %>%
  mutate(
    day = as_date(start_dt),
    vaso_flag = as.integer(str_detect(med_name, str_c(vaso_names, collapse = "|")))
  ) %>%
  group_by(hospitalization_id, day) %>%
  summarise(any_vaso = as.integer(any(vaso_flag == 1L, na.rm = TRUE)), .groups = "drop")

# --- Build SpO2 daily summary ---------------------------------------------

vit_id_col   <- pick_name(vitals_raw, c("hospitalization_id", "encounter_id", "stay_id"))
# helper: tolerant numeric parser (handles "95", "95%", " 95 % ", etc.)
num_safely <- function(x) {
  x_chr <- as.character(x)
  # keep digits, dot, and minus; drop everything else (e.g., %)
  x_num <- suppressWarnings(as.numeric(str_replace_all(x_chr, "[^0-9\\.-]", "")))
  x_num
}

# time column present in your schema
vit_time_col <- "recorded_dttm"

# patterns that typically indicate SpO2 in name/category
spo2_pattern <- regex("(^|\\b)(spo2|sp[ _-]?o2|o2[ _-]?sat|oxygen[ _-]?saturation|pulse[ _-]?ox(?!imeter))($|\\b)",
                      ignore_case = TRUE)

spo2_daily <- vitals_raw %>%
  transmute(
    hospitalization_id,
    time  = suppressWarnings(lubridate::ymd_hms(.data[[vit_time_col]], tz = "UTC")),
    day   = as_date(time),
    vname = tolower(trimws(vital_name)),
    vcat  = tolower(trimws(vital_category)),
    val   = num_safely(vital_value)
  ) %>%
  filter(
    # keep rows that look like oxygen saturation by name OR category
    str_detect(vname, spo2_pattern) | str_detect(vcat, spo2_pattern)
  ) %>%
  filter(!is.na(val)) %>%
  group_by(hospitalization_id, day) %>%
  summarise(
    median_spo2 = median(val, na.rm = TRUE),
    .groups = "drop"
  )

# ----- 1) Build a clean, fixed-length ICU panel (e.g., first 7 days) -----
horizon <- 7L

icu_days <- cohort_stays %>%
  transmute(
    hospitalization_id,
    start = as.Date(lubridate::floor_date(first_icu_in, "day")),
    end   = as.Date(lubridate::floor_date(last_icu_out, "day"))
  ) %>%
  mutate(end = if_else(end < start, start, end)) %>%
  rowwise() %>%
  mutate(day = list(seq.Date(start, end, by = "day"))) %>%
  unnest(day) %>%
  group_by(hospitalization_id) %>%
  mutate(icu_day = as.integer(day - min(day)) + 1L) %>%
  ungroup() %>%
  filter(icu_day <= horizon)

# Join your daily summaries (built earlier)
icu_panel <- icu_days %>%
  left_join(resp,       by = c("hospitalization_id","day")) %>%
  left_join(vaso,       by = c("hospitalization_id","day")) %>%
  left_join(spo2_daily, by = c("hospitalization_id","day")) %>%
  group_by(hospitalization_id, icu_day) %>%
  summarise(
    any_imv    = as.integer(any(any_imv == 1L, na.rm = TRUE)),
    any_vaso   = as.integer(any(any_vaso == 1L, na.rm = TRUE)),
    hypox_burden = 1 - pmin(median_spo2, 100)/100,  # 0..1; NA if absent
    .groups = "drop"
  ) %>%
  # ensure complete day grid per encounter
  group_by(hospitalization_id) %>%
  tidyr::complete(
    icu_day = 1:horizon,
    fill = list(any_imv = 0L, any_vaso = 0L, hypox_burden = 0)
  ) %>%
  ungroup()

# ----- 2) Build per-encounter vectors of equal length -----
series_df <- icu_panel %>%
  arrange(hospitalization_id, icu_day) %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_vec   = list(as.numeric(any_imv)),
    vaso_vec  = list(as.numeric(any_vaso)),
    hypox_vec = list(as.numeric(replace_na(hypox_burden, 0))),
    .groups = "drop"
  ) %>%
  # drop degenerate series (all zeros on all three)
  mutate(nonzero = map_lgl(imv_vec, ~any(.x != 0)) |
           map_lgl(vaso_vec, ~any(.x != 0)) |
           map_lgl(hypox_vec, ~any(.x != 0))) %>%
  filter(nonzero) %>%
  dplyr::select(-nonzero)

# quick guardrails
if (nrow(series_df) == 0) {
  stop("All series are empty/degenerate after filling. Check joins into icu_panel (resp/vaso/spo2) and horizon.")
}

# ----- 3) Compose a single trajectory vector per encounter -----
stack_series <- function(imv, vaso, hypox) {
  z <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) return(rep(0, length(x)))
    s <- stats::sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - m) / s
  }
  rowMeans(cbind(z(unlist(imv)), z(unlist(vaso)), z(unlist(hypox))), na.rm = TRUE)
}

series_df$traj <- purrr::pmap(
  list(series_df$imv_vec, series_df$vaso_vec, series_df$hypox_vec),
  stack_series
)

# ----- 4) Run DTW clustering safely -----
ts_list <- series_df$traj
# ensure order alignment between ts_list and series_df rows
# (pmap already preserves row order; no need to set names)
n_series <- length(ts_list)
k <- min(4L, max(2L, n_series))  # at least 2, at most 4, not exceeding series count

set.seed(123)
clust <- tsclust(
  ts_list,
  type = "partitional",
  k = k,
  distance = "dtw_basic",
  centroid = "dba",
  seed = 123,
  trace = TRUE,
  args = tsclust_args(dist = list(window.size = 2))
)

# attach cluster labels by position (not by names)
series_df$traj_cluster <- factor(clust@cluster, levels = 1:k, labels = paste0("C", 1:k))

# sanity: how many per cluster?
print(series_df %>% count(traj_cluster))

# ----- 5) Join clusters back to cohort -----
cohort_traj <- cohort %>%
  dplyr::select(patient_id, hospitalization_id, hypoxemic_arf, hypercapnic_arf, mixed_arf,
         age_years, sex_category, race_category, ethnicity_category, county_fips, admission_dttm) %>%
  left_join(series_df %>% dplyr::select(hospitalization_id, traj_cluster), by = "hospitalization_id")

# Now this should be > 0 if clustering succeeded:
analytic <- cohort_traj %>%
  left_join(cohort_expo, by = c("patient_id","hospitalization_id")) %>%
  filter(!is.na(traj_cluster))

print(nrow(analytic))

analytic_clean <- analytic %>%
  # drop duplicate .y versions we don’t need
  dplyr::select(
    -ends_with(".y"),
    -admission_dttm.y,
    -county_fips.y
  ) %>%
  # rename .x columns to remove suffix
  rename(
    hypoxemic_arf   = hypoxemic_arf.x,
    hypercapnic_arf = hypercapnic_arf.x,
    mixed_arf       = mixed_arf.x,
    age_years       = age_years.x,
    sex_category    = sex_category.x,
    race_category   = race_category.x,
    ethnicity_category = ethnicity_category.x,
    county_fips     = county_fips.x,
    admission_dttm  = admission_dttm.x
  )

analytic_clean <- analytic_clean %>%
  mutate(
    race_name = race_category,   # keep naming consistent
    re_low = str_to_lower(paste(ethnicity_category, race_name)),
    
    # flags
    is_nonhisp = str_detect(re_low, "\\bnon[- ]?hispanic\\b"),
    is_hisp    = str_detect(re_low, "\\bhispanic\\b") & !is_nonhisp,
    is_white   = str_detect(re_low, "white"),
    is_black   = str_detect(re_low, "black"),
    is_asian_any = str_detect(
      re_low,
      "asian|mideast|filipino|chinese|korean|vietnamese|pacific islander|samoan"
    ),
    
    race_ethnicity_simple = case_when(
      is_white & is_hisp    ~ "Hispanic White",
      is_white & is_nonhisp ~ "Non-Hispanic White",
      is_black & is_hisp    ~ "Hispanic Black",
      is_black & is_nonhisp ~ "Non-Hispanic Black",
      is_asian_any          ~ "Asian",
      TRUE                  ~ "Other"
    )
  ) %>%
  dplyr::select(-re_low, -is_nonhisp, -is_hisp, -is_white, -is_black, -is_asian_any) %>%
  mutate(
    sex_category = factor(sex_category),
    race_ethnicity_simple = factor(
      race_ethnicity_simple,
      levels = c("Non-Hispanic White", "Hispanic White",
                 "Non-Hispanic Black", "Hispanic Black",
                 "Asian", "Other")
    )
  )


# Example fixed set of predictors; swap in your preferred set
expo_covars <- c("no2_mean_3y", "pm25_mean_3y", "svi_last",
                 "tmax_mean_3y", "vp_mean_3y")
if (has_no2_monthly) expo_covars <- unique(c(expo_covars, "no2_mean_12mo"))

demo_covars <- c("age_years", "sex_category", "race_ethnicity_simple",
                 "hypoxemic_arf", "hypercapnic_arf", "mixed_arf")

formula_mn <- as.formula(
  paste("traj_cluster ~", paste(c(expo_covars, demo_covars), collapse = " + "))
)

set.seed(123)
fit_mn <- nnet::multinom(formula_mn, data = analytic_clean, trace = FALSE)

# Tidy effect estimates (log-odds). To get RRRs (relative risk ratios):
coef_mn <- broom::tidy(fit_mn, exponentiate = TRUE, conf.int = TRUE)
coef_mn

# Cluster size
analytic_clean %>%
  count(traj_cluster) %>%
  ggplot(aes(traj_cluster, n)) + geom_col() + labs(x = "Trajectory class", y = "N")

# Exposure by cluster (e.g., monthly NO2 if present)
if (has_no2_monthly) {
  analytic_clean %>%
    ggplot(aes(traj_cluster, no2_mean_12mo)) +
    geom_boxplot() +
    labs(x = "Trajectory class", y = "12-mo mean NO₂ (ppb)")
}


# keep only non-intercept terms
plot_df <- coef_mn %>%
  filter(term != "(Intercept)") %>%
  # reorder for nicer plotting
  mutate(
    term = fct_reorder(term, estimate),
    cluster = y.level
  )

ggplot(plot_df, aes(x = estimate, y = term, color = cluster)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 position = position_dodge(width = 0.6), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  labs(
    x = "Relative Risk Ratio (log scale)",
    y = NULL,
    title = "Environmental and Demographic Predictors \nof ARF Trajectory Cluster",
    subtitle = "Multinomial regression; reference cluster = C1",
    color = "Trajectory cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(face = "bold")
  )

traj_profiles <- icu_panel %>%
  left_join(series_df %>% dplyr::select(hospitalization_id, traj_cluster), 
            by = "hospitalization_id") %>%
  group_by(traj_cluster, icu_day) %>%
  summarise(
    prop_imv   = mean(any_imv, na.rm = TRUE),
    prop_vaso  = mean(any_vaso, na.rm = TRUE),
    mean_fio2  = mean(resp$max_fio2, na.rm = TRUE),   # adjust if stored elsewhere
    mean_spo2  = mean(100 * (1 - hypox_burden), na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

######################## HOURLY Traj


h_cap <- 168L  # 7 days of hours; change to 240L for 10 days, etc.

cohort_stays_hr <- cohort %>%
  mutate(
    first_icu_in = as_datetime(first_icu_in, tz = "UTC"),
    last_icu_out = as_datetime(last_icu_out, tz = "UTC")
  ) %>%
  filter(!is.na(first_icu_in), !is.na(last_icu_out), last_icu_out > first_icu_in) %>%
  transmute(hospitalization_id,
            start_hr = floor_date(first_icu_in, "hour"),
            end_hr   = floor_date(last_icu_out, "hour"))

hour_grid <- cohort_stays_hr %>%
  rowwise() %>%
  mutate(hr_seq = list(seq.POSIXt(start_hr, end_hr, by = "hour"))) %>%
  unnest(hr_seq) %>%
  group_by(hospitalization_id) %>%
  mutate(hour_idx = as.integer(difftime(hr_seq, min(hr_seq), units = "hours")) + 1L) %>%
  ungroup() %>%
  filter(hour_idx <= h_cap) %>%
  rename(hour_ts = hr_seq)

# helper: 0/1 coercer
as01 <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("1","true","t","yes","y") ~ 1L,
    x %in% c("0","false","f","no","n","") ~ 0L,
    suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) > 0) ~ 1L,
    suppressWarnings(!is.na(as.numeric(x)) & as.numeric(x) == 0) ~ 0L,
    TRUE ~ NA_integer_
  )
}

# normalize FiO2 to percent
norm_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  case_when(is.na(x) ~ NA_real_, x <= 1 ~ x*100, TRUE ~ x)
}

# Pre-filter rs_raw to ICU range to cut rows early
rs_small <- rs_raw %>%
  semi_join(
    cohort %>%
      mutate(icu_start = as_datetime(first_icu_in, tz="UTC"),
             icu_end   = as_datetime(last_icu_out, tz="UTC")) %>%
      select(hospitalization_id, icu_start, icu_end),
    by = "hospitalization_id"
  ) %>%
  # cheap filter by recorded_dttm before parse (if available)
  filter(!is.na(recorded_dttm))

# Precompute regex
re_ac_vc <- regex("assist.?control.*volume|ac.?vc", ignore_case = TRUE)
re_ac_pc <- regex("assist.?control.*pressure|ac.?pc", ignore_case = TRUE)
re_simv  <- fixed("simv", ignore_case = TRUE)
re_psv   <- regex("psv|pressure support", ignore_case = TRUE)
re_inv   <- regex("invasive|mechanical|vent", ignore_case = TRUE)

# robust parser that tries several formats + epoch
parse_ts <- function(x, tz = "UTC") {
  x_chr <- as.character(x)
  
  # if looks like pure digits of epoch seconds or milliseconds
  is_epoch <- str_detect(x_chr, "^\\d{10}(?:\\d{3})?$")
  out <- rep(NA_real_, length(x_chr))
  
  # epoch seconds/millis
  if (any(is_epoch)) {
    xe <- x_chr[is_epoch]
    # 13-digit -> ms; 10-digit -> s
    out[is_epoch] <- as.numeric(ifelse(nchar(xe) >= 13, as.numeric(xe)/1000, as.numeric(xe)))
  }
  
  # try multiple text formats for the rest
  need <- is.na(out) & !is.na(x_chr)
  if (any(need)) {
    out[need] <- suppressWarnings(lubridate::parse_date_time(
      x_chr[need],
      orders = c(          # generous and fast enough
        "ymd HMS","ymd HM","ymd H",
        "Ymd HMS","Ymd HM","Ymd H",
        "mdy HMS","mdy HM","mdy H",
        "dmy HMS","dmy HM","dmy H"
      ),
      tz = tz,
      truncated = 3
    )) |> as.numeric()
  }
  
  # return POSIXct
  as.POSIXct(out, origin = "1970-01-01", tz = tz)
}

# ---- start here instead of ymd_hms() ----
rs_small_parsed <- rs_small %>%
  mutate(
    time_raw = coalesce(as.character(recorded_dttm), as.character(recorded_time)),
    time     = parse_ts(time_raw, tz = "UTC")
  )

# quick diagnostic (optional)
n_failed <- sum(is.na(rs_small_parsed$time) & !is.na(rs_small_parsed$time_raw))
message("Unparseable timestamps: ", n_failed)

rs_step1 <- rs_small_parsed %>%
  filter(!is.na(time)) %>%
  select(
    hospitalization_id, time,
    recorded_dttm, recorded_time,
    device_category, mode_category, mode_name,
    artificial_airway, tracheostomy,
    fio2_set, tidal_volume_set, resp_rate_set,
    pressure_control_set, pressure_support_set, peep_set,
    peak_inspiratory_pressure_obs, plateau_pressure_obs, peep_obs, minute_vent_obs
  ) %>%
  mutate(hour_ts = floor_date(time, "hour"))

# pre-cast numerics once; suppressWarnings in one place
num_cols <- c("fio2_set","tidal_volume_set","resp_rate_set","pressure_control_set",
              "pressure_support_set","peep_set","peak_inspiratory_pressure_obs",
              "plateau_pressure_obs","peep_obs","minute_vent_obs")

rs_step2 <- rs_step1 %>%
  mutate(
    device_category_lc = tolower(device_category),
    mode_category_lc   = tolower(mode_category),
    mode_name_lc       = tolower(mode_name),
    across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))),
    artificial_airway  = as01(artificial_airway),
    tracheostomy       = as01(tracheostomy)
  ) %>%
  # normalize FiO2 to percent
  mutate(fio2_set = if_else(is.na(fio2_set), NA_real_, if_else(fio2_set <= 1, fio2_set*100, fio2_set)))

rs_step3 <- rs_step2 %>%
  mutate(
    mode_ac_vc = str_detect(mode_category_lc, re_ac_vc),
    mode_ac_pc = str_detect(mode_category_lc, re_ac_pc),
    mode_simv  = str_detect(mode_category_lc, re_simv),
    mode_psv   = str_detect(mode_category_lc, re_psv),
    mode_invasive = (artificial_airway == 1L | tracheostomy == 1L |
                       str_detect(device_category_lc, re_inv))
  )


rs_step4 <- rs_step3 %>%
  transmute(
    hospitalization_id, hour_ts, mode_invasive,
    set_rr   = if_else(mode_ac_vc | mode_ac_pc | mode_simv, resp_rate_set, NA_real_),
    set_vt   = if_else(mode_ac_vc | mode_simv, tidal_volume_set, NA_real_),
    set_pc   = if_else(mode_ac_pc, pressure_control_set, NA_real_),
    set_ps   = if_else(mode_psv | mode_simv, pressure_support_set, NA_real_),
    set_peep = coalesce(peep_obs, peep_set),
    set_fio2 = fio2_set,
    pip_obs  = peak_inspiratory_pressure_obs,
    pplat_obs= plateau_pressure_obs,
    mv_obs   = minute_vent_obs
  )

# Helpers: no warnings, return NA if a group has no finite values
safe_max <- function(x) {
  v <- x[is.finite(x)]
  if (length(v) == 0) NA_real_ else max(v)
}
safe_med <- function(x) {
  v <- x[is.finite(x)]
  if (length(v) == 0) NA_real_ else stats::median(v)
}

rs_hr <- rs_step4 %>%
  summarise(
    any_imv   = as.integer(any(mode_invasive == 1L, na.rm = TRUE)),
    fio2_max  = safe_max(set_fio2),
    peep_med  = safe_med(set_peep),
    rr_med    = safe_med(set_rr),
    vt_med    = safe_med(set_vt),
    pc_med    = safe_med(set_pc),
    ps_med    = safe_med(set_ps),
    pip_med   = safe_med(pip_obs),
    pplat_med = safe_med(pplat_obs),
    mv_med    = safe_med(mv_obs),
    .by = c(hospitalization_id, hour_ts)
  )


### clean up environment

# 2) Decide what to KEEP for the next steps (edit if needed)
keep <- c(
  "cohort",             # cohort-level metadata
  "rs_raw",             # respiratory support raw
  "vitals_raw",         # vitals raw (long form)
  "clif_labs",          # labs (for ABGs)
  "exposome_dir",       # path to exposome files
  "no2_mo", "pm25_mo",  # monthly exposures if already loaded (optional)
  "cohort_expo",        # exposure-joined cohort (optional)
  "rs_hr",              # already computed hourly resp support
  "hour_grid",          # ICU hourly grid
  "parse_ts", "as01",   # helper functions you've defined
  "re_ac_vc","re_ac_pc","re_simv","re_psv","re_inv"  # compiled regex
)

# 3) Drop everything else and GC
rm(list = setdiff(ls(), keep))
gc(full = TRUE)

library(DBI)
library(duckdb)

# Connect (reuse your existing connection if you already have one)
con <- dbConnect(duckdb::duckdb(), dbdir = tempfile())

# If you haven't already copied vitals_raw into DuckDB:
copy_to(con, vitals_raw, "vitals", temporary = FALSE, overwrite = TRUE)

# Build the query with dbplyr; use sql() where we need explicit DuckDB functions
spo2_tbl <- tbl(con, "vitals") %>%
  transmute(
    hospitalization_id,
    recorded_dttm,
    vital_name,
    vital_category,
    vital_value
  ) %>%
  mutate(
    vname_lc = sql("lower(vital_name)"),
    vcat_lc  = sql("lower(vital_category)"),
    looks_spo2 = sql(
      "regexp_matches(vname_lc, 'spo2|sp[ _-]?o2|o2[ _-]?sat|oxygen[ _-]?saturation')
       OR regexp_matches(vcat_lc,  'spo2|o2|oxygen')"
    )
  ) %>%
  filter(looks_spo2) %>%
  mutate(
    ts = sql("try_cast(recorded_dttm AS TIMESTAMP)"),
    hour_ts = sql("date_trunc('hour', ts)"),
    val = sql("try_cast(regexp_replace(CAST(vital_value AS VARCHAR), '[^0-9\\.-]', '', 'g') AS DOUBLE)")
  ) %>%
  filter(!is.na(hour_ts), !is.na(val)) %>%
  group_by(hospitalization_id, hour_ts) %>%
  summarise(
    spo2_med = sql("median(val)"),  # NA-safe in SQL, avoids dbplyr NA warning
    .groups = "drop"
  )

# Materialize to R (already aggregated)
spo2_hr <- collect(spo2_tbl)

# Optional: close DB when done
dbDisconnect(con, shutdown = TRUE)


# Optional: drop objects used above
rm(vitals_min, ids, id_chunks); gc()

clif_labs <- clif_tables[[4]]

# ABGs (pH, PaCO2, PaO2)
labs_hr <- clif_labs %>%
  transmute(
    hospitalization_id,
    time = ymd_hms(result_dttm, tz = "UTC"),
    hour_ts = floor_date(time, "hour"),
    lname = tolower(lab_name),
    val   = num_safely(result_value)
  ) %>%
  mutate(
    ph   = if_else(str_detect(lname, "\\bph\\b|arterial ph|blood gas ph"), val, NA_real_),
    paco2= if_else(str_detect(lname, "paco2|pa co2|arterial co2"), val, NA_real_),
    pao2 = if_else(str_detect(lname, "pao2|pa o2|arterial o2"), val, NA_real_)
  ) %>%
  group_by(hospitalization_id, hour_ts) %>%
  summarise(
    ph_med    = median(ph,    na.rm = TRUE),
    paco2_med = median(paco2, na.rm = TRUE),
    pao2_med  = median(pao2,  na.rm = TRUE),
    .groups = "drop"
  )



icu_hourly <- hour_grid %>%
  left_join(rs_hr,   by = c("hospitalization_id","hour_ts")) %>%
  left_join(spo2_hr, by = c("hospitalization_id","hour_ts")) %>%
  left_join(labs_hr, by = c("hospitalization_id","hour_ts")) %>%
  arrange(hospitalization_id, hour_idx) %>%
  group_by(hospitalization_id) %>%
  # carry last obs forward within the capped window
  mutate(across(c(any_imv, fio2_max, peep_med, rr_med, vt_med, pc_med, ps_med,
                  pip_med, pplat_med, mv_med, spo2_med, ph_med, paco2_med, pao2_med),
                ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
  ungroup() %>%
  # replace remaining NAs with clinically neutral values
  mutate(
    any_imv   = replace_na(any_imv, 0L),
    spo2_med  = replace_na(spo2_med, 97),   # benign saturation
    fio2_max  = replace_na(fio2_max, 40),   # room air-ish to low O2
    peep_med  = replace_na(peep_med, 5),
    rr_med    = replace_na(rr_med, 14),
    vt_med    = replace_na(vt_med, 450),
    pc_med    = replace_na(pc_med, 0),
    ps_med    = replace_na(ps_med, 0),
    pip_med   = replace_na(pip_med, NA_real_),
    pplat_med = replace_na(pplat_med, NA_real_),
    mv_med    = replace_na(mv_med, NA_real_),
    ph_med    = replace_na(ph_med, 7.40),
    paco2_med = replace_na(paco2_med, 40),
    pao2_med  = replace_na(pao2_med, 80)
  )



no2_mo <- readr::read_csv(file.path(exposome_dir, "no2_county_month.csv"), show_col_types=FALSE) %>%
  transmute(county_fips = as.character(county_fips),
            ym = make_date(year, month, 1),
            no2 = as.double(no2))

pm25_mo <- readr::read_csv(file.path(exposome_dir, "pm25_county_month.csv"), show_col_types=FALSE) %>%
  transmute(county_fips = as.character(county_fips),
            ym = make_date(year, month, 1),
            pm25 = as.double(pm25))

# build 12-month monthly sequence ending the month before admission, then interpolate to hourly
build_exposure_embedding <- function(fips, admit_ts) {
  end_m   <- floor_date(admit_ts %m-% months(1), "month")
  start_m <- end_m %m-% months(11)
  mo_grid <- tibble(ym = seq.Date(start_m, end_m, by = "month"))
  # join county series
  series <- mo_grid %>%
    left_join(filter(no2_mo, county_fips == fips), by="ym") %>%
    left_join(filter(pm25_mo, county_fips == fips), by="ym")
  # fallback fill with county means if gaps (rare)
  if (anyNA(series$no2))  series$no2  <- zoo::na.locf(series$no2, na.rm = FALSE)
  if (anyNA(series$pm25)) series$pm25 <- zoo::na.locf(series$pm25, na.rm = FALSE)
  
  # hourly time grid for past h_cap hours before admission
  hr_grid <- tibble(hour_ts = seq.POSIXt(admit_ts - hours(h_cap-1), admit_ts, by = "hour"))
  
  # linear interpolation month→hour using approx
  # map monthly points to their POSIXct timestamps (start of month)
  x_no2  <- as.numeric(as.POSIXct(series$ym, tz = "UTC"))
  y_no2  <- series$no2
  x_pm25 <- as.numeric(as.POSIXct(series$ym, tz = "UTC"))
  y_pm25 <- series$pm25
  xout   <- as.numeric(as.POSIXct(hr_grid$hour_ts, tz = "UTC"))
  
  no2_hr  <- approx(x = x_no2,  y = y_no2,  xout = xout, method = "linear", rule = 2)$y
  pm25_hr <- approx(x = x_pm25, y = y_pm25, xout = xout, method = "linear", rule = 2)$y
  
  tibble(hour_ts = hr_grid$hour_ts, no2_hist = no2_hr, pm25_hist = pm25_hr)
}

# build once per encounter
admit_tbl <- cohort %>%
  transmute(hospitalization_id, county_fips = as.character(county_fips),
            admit_ts = as_datetime(admission_dttm, tz = "UTC"))

embeddings <- admit_tbl %>%
  mutate(embed = pmap(list(county_fips, admit_ts), build_exposure_embedding)) %>%
  select(hospitalization_id, embed) %>%
  unnest(embed)



traj_hr <- icu_hourly %>%
  left_join(embeddings, by = c("hospitalization_id","hour_ts")) %>%
  mutate(
    no2_hist  = replace_na(no2_hist, median(no2_hist, na.rm = TRUE)),
    pm25_hist = replace_na(pm25_hist, median(pm25_hist, na.rm = TRUE))
  ) %>%
  # build a composite signal for DTW (z-score each feature, then average)
  group_by(hospitalization_id) %>%
  arrange(hour_idx, .by_group = TRUE) %>%
  mutate(
    z = function(x){ s <- sd(x, na.rm=TRUE); m <- mean(x, na.rm=TRUE); ifelse(is.na(s) || s==0, 0, (x-m)/s) },
    sig = rowMeans(cbind(
      z(as.numeric(any_imv)),
      z(fio2_max), z(peep_med), z(rr_med), z(vt_med), z(pc_med), z(ps_med),
      z(spo2_med * -1),   # lower SpO2 worse → invert
      z(ph_med * -1),     # acidemia worse → invert
      z(paco2_med),       # hypercapnia worse
      z(pao2_med * -1),   # lower PaO2 worse
      z(no2_hist), z(pm25_hist)   # EMBEDDED exposures
    ), na.rm = TRUE)
  ) %>%
  ungroup()



series_df_hr <- traj_hr %>%
  select(hospitalization_id, hour_idx, sig) %>%
  arrange(hospitalization_id, hour_idx) %>%
  group_by(hospitalization_id) %>%
  summarise(traj = list(as.numeric(sig)), .groups = "drop") %>%
  # drop degenerate series
  filter(map_lgl(traj, ~ any(is.finite(.x) & .x != 0)))

ts_list <- series_df_hr$traj
set.seed(123)
k <- 4  # tune as needed

clust_hr <- tsclust(
  ts_list, type = "partitional", k = k,
  distance = "dtw_basic", centroid = "dba",
  seed = 123, trace = TRUE,
  args = tsclust_args(dist = list(window.size = 6))  # small Sakoe-Chiba band for hourly misalignment
)

series_df_hr$traj_cluster_hr <- factor(clust_hr@cluster, levels = 1:k, labels = paste0("C",1:k))

# Attach back to cohort for modeling
cohort_traj_hr <- cohort %>%
  left_join(series_df_hr %>% select(hospitalization_id, traj_cluster_hr), by = "hospitalization_id")






