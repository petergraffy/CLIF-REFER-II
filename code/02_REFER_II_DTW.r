# ===============================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) Index — Full Pipeline
# PI: Peter Graffy (graffy@uchicago.edu)
# Purpose: Build hourly -> daily clinical trajectories, daily NO2/PM2.5 exposures, 
#          multivariate DTW clustering, prototypes, and multinomial membership model.
# Run after you have: cohort_min, clif_tables
# ===============================================================================================

# ---- Libraries -------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(stringr)
  library(purrr)
  library(data.table)
  library(dtwclust)
  library(nnet)
  library(ggplot2)
  library(DBI)
  library(duckdb)
  library(future)
  library(forcats)
  library(zoo)
  library(arrow)
  library(tibble)
  library(cluster)
})

# ---- Helpers ----------------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

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

num_safely <- function(x) as.numeric(str_replace_all(as.character(x), "[^0-9\\.-]", ""))

as_posix <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) return(x)
  suppressWarnings(lubridate::ymd_hms(x, tz = tz))
}

pick_name <- function(df, candidates) {
  nm <- intersect(candidates, names(df))
  if (length(nm) == 0) stop("None of the candidate columns found: ", paste(candidates, collapse = ", "))
  nm[1]
}

find_tbl <- function(tbl_list, patterns) {
  stopifnot(is.list(tbl_list))
  nms <- names(tbl_list)
  idx <- which(Reduce(`|`, lapply(patterns, function(p) str_detect(tolower(nms), tolower(p)))))
  if (length(idx) == 0) stop("Could not find a table matching patterns: ", paste(patterns, collapse = " | "))
  tbl_list[[idx[1]]]
}

# Robust timestamp parser that handles multiple formats + epoch
parse_ts <- function(x, tz = "UTC") {
  x_chr <- as.character(x)
  is_epoch <- str_detect(x_chr, "^\\d{10}(?:\\d{3})?$")
  out <- rep(NA_real_, length(x_chr))
  if (any(is_epoch)) {
    xe <- x_chr[is_epoch]
    out[is_epoch] <- as.numeric(ifelse(nchar(xe) >= 13, as.numeric(xe)/1000, as.numeric(xe)))
  }
  need <- is.na(out) & !is.na(x_chr)
  if (any(need)) {
    out[need] <- suppressWarnings(lubridate::parse_date_time(
      x_chr[need],
      orders = c("ymd HMS","ymd HM","ymd H","Ymd HMS","Ymd HM","Ymd H","mdy HMS","mdy HM","mdy H","dmy HMS","dmy HM","dmy H"),
      tz = tz, truncated = 3
    )) |> as.numeric()
  }
  as.POSIXct(out, origin = "1970-01-01", tz = tz)
}

safe_max <- function(x) { v <- x[is.finite(x)]; if (length(v) == 0) NA_real_ else max(v) }
safe_med <- function(x) { v <- x[is.finite(x)]; if (length(v) == 0) NA_real_ else stats::median(v) }
z_within <- function(x) { m <- mean(x, na.rm=TRUE); s <- sd(x, na.rm=TRUE); if (!is.finite(s) || s==0) rep(0, length(x)) else (x-m)/s }

# Interpolate monthly series to daily within a date window
interp_monthly_to_daily <- function(df_mo, value_col, start_date, end_date) {
  mo_key <- df_mo %>% arrange(ym)
  if (nrow(mo_key) == 0) return(tibble(date = seq.Date(start_date, end_date, "day"), val = NA_real_))
  max_m <- max(mo_key$ym)
  anchors <- bind_rows(
    mo_key %>% transmute(date = ym, val = .data[[value_col]]),
    tibble(date = seq.Date(max_m, by = "month", length.out = 2)[2], val = tail(mo_key[[value_col]], 1))
  ) %>% arrange(date)
  daily <- tibble(date = seq.Date(start_date, end_date, by = "day"))
  daily$val <- approx(x = as.numeric(anchors$date), y = anchors$val,
                      xout = as.numeric(daily$date), method = "linear",
                      rule = 2, ties = "ordered")$y
  daily
}

# DTW helpers (daily)
safe_dtw_daily <- function(x, y, band_primary = 1L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = asymmetricP2, window.type = "sakoechiba",
             window.size = band_primary, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = asymmetricP2, distance.only = TRUE)$distance,
           error = function(e) Inf)
}

get_medoid_idx_mv <- function(series_mats, dist_fun) {
  m <- length(series_mats)
  if (m <= 1) return(1L)
  dsum <- numeric(m)
  for (i in seq_len(m)) {
    xi1 <- series_mats[[i]][1,]; xi2 <- series_mats[[i]][2,]
    s <- 0
    for (j in seq_len(m)) if (i != j) {
      yj1 <- series_mats[[j]][1,]; yj2 <- series_mats[[j]][2,]
      s <- s + dist_fun(xi1, yj1) + dist_fun(xi2, yj2)
    }
    dsum[i] <- s
  }
  which.min(dsum)
}

# ---- Exposome inputs (Annual + Monthly + optional Daily) --------------------------------------
exposome_dir <- "exposome"

pm25_yr <- read_csv(file.path(exposome_dir, "pm25_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, pm25 = pm25_mean)

no2_yr  <- read_csv(file.path(exposome_dir, "no2_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, no2 = no2_mean)

svi_yr  <- read_csv(file.path(exposome_dir, "svi_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, svi = svi_overall)

tmax_yr <- read_csv(file.path(exposome_dir, "daymet_tmax_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, tmax = tmax_mean)

tmin_yr <- read_csv(file.path(exposome_dir, "daymet_tmin_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, tmin = tmin_mean)

vp_yr   <- read_csv(file.path(exposome_dir, "daymet_vp_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, vp = vp_mean)

prcp_yr <- read_csv(file.path(exposome_dir, "daymet_prcp_county_year.csv"), show_col_types = FALSE) %>%
  rename(county_fips = GEOID, year = year, prcp = prcp_mean)

# Monthly
no2_mo <- readr::read_csv(file.path(exposome_dir, "no2_county_month.csv"), show_col_types = FALSE) %>%
  transmute(county_fips = as.character(county_fips), ym = make_date(year, month, 1), no2 = as.double(no2))
pm25_mo <- readr::read_csv(file.path(exposome_dir, "pm25_county_month.csv"), show_col_types = FALSE) %>%
  transmute(county_fips = as.character(county_id), ym = make_date(year, month, 1), pm25 = as.double(pm25))

# Optional daily files (if present)
pm25_daily_path <- file.path(exposome_dir, "pm25_county_day.csv")
no2_daily_path  <- file.path(exposome_dir, "no2_county_day.csv")
pm25_day <- if (file.exists(pm25_daily_path)) readr::read_csv(pm25_daily_path, show_col_types = FALSE) %>%
  transmute(county_fips = as.character(county_fips), date = as.Date(date), pm25 = as.double(pm25)) else NULL
no2_day  <- if (file.exists(no2_daily_path)) readr::read_csv(no2_daily_path, show_col_types = FALSE) %>%
  transmute(county_fips = as.character(county_fips), date = as.Date(date), no2  = as.double(no2)) else NULL

# National month means + overall medians (fallback)
no2_nat_by_month  <- no2_mo  %>% group_by(ym) %>% summarise(no2_nat  = mean(no2,  na.rm = TRUE), .groups = "drop")
pm25_nat_by_month <- pm25_mo %>% group_by(ym) %>% summarise(pm25_nat = mean(pm25, na.rm = TRUE), .groups = "drop")
no2_overall_med   <- median(no2_mo$no2,   na.rm = TRUE)
pm25_overall_med  <- median(pm25_mo$pm25, na.rm = TRUE)

# ---- Cohort & annual exposures (3y lags) ------------------------------------------------------
cohort <- cohort_min %>%
  mutate(
    admit_dt   = as_date(admission_dttm),
    admit_year = year(admit_dt),
    county_fips = county_code
  )

expo_annual <- pm25_yr %>%
  full_join(no2_yr,  by = c("county_fips","year")) %>%
  full_join(svi_yr,  by = c("county_fips","year")) %>%
  full_join(tmax_yr, by = c("county_fips","year")) %>%
  full_join(tmin_yr, by = c("county_fips","year")) %>%
  full_join(vp_yr,   by = c("county_fips","year")) %>%
  full_join(prcp_yr, by = c("county_fips","year"))

attach_annual_lags <- function(df, expo_tbl) {
  df %>%
    tidyr::expand_grid(lag = 1:3) %>%
    mutate(expo_year = admit_year - lag) %>%
    left_join(expo_tbl, by = c("county_fips","expo_year" = "year")) %>%
    group_by(patient_id, hospitalization_id) %>%
    summarise(
      pm25_mean_3y = mean(pm25, na.rm = TRUE),
      no2_mean_3y  = mean(no2,  na.rm = TRUE),
      svi_last     = dplyr::last(na.omit(svi)),
      tmax_mean_3y = mean(tmax, na.rm = TRUE),
      tmin_mean_3y = mean(tmin, na.rm = TRUE),
      vp_mean_3y   = mean(vp,   na.rm = TRUE),
      prcp_mean_3y = mean(prcp, na.rm = TRUE),
      .groups = "drop"
    )
}

cohort_expo_annual <- attach_annual_lags(cohort, expo_annual)

# Monthly NO2: 12-mo pre-admission mean
no2win <- no2_mo %>% select(county_fips, ym, no2)
cohort_no2_12mo <- cohort %>%
  transmute(patient_id, hospitalization_id, county_fips, 
            start = floor_date(admit_dt %m-% months(12), "month"),
            end   = floor_date(admit_dt %m-% months(1), "month")) %>%
  rowwise() %>% mutate(ym_seq = list(seq.Date(start, end, by = "month"))) %>%
  unnest(ym_seq) %>%
  left_join(no2win, by = c("county_fips","ym_seq" = "ym")) %>%
  group_by(patient_id, hospitalization_id) %>%
  summarise(no2_mean_12mo = mean(no2, na.rm = TRUE), .groups = "drop")

cohort_expo <- cohort %>%
  left_join(cohort_expo_annual, by = c("patient_id","hospitalization_id")) %>%
  left_join(cohort_no2_12mo,    by = c("patient_id","hospitalization_id"))

# ---- ICU stay windows (hours & days) ----------------------------------------------------------
cohort_stays <- cohort %>%
  mutate(
    first_icu_in  = as_datetime(first_icu_in, tz = "UTC"),
    last_icu_out  = as_datetime(last_icu_out, tz = "UTC")
  ) %>%
  filter(!is.na(first_icu_in), !is.na(last_icu_out), last_icu_out > first_icu_in)

# Hours grid (cap to h_cap hours)
h_cap <- 24L * 7L  # 7 days
cohort_stays_hr <- cohort_stays %>%
  transmute(hospitalization_id,
            start_hr = floor_date(first_icu_in, "hour"),
            end_hr   = floor_date(last_icu_out, "hour"))
hour_grid <- cohort_stays_hr %>%
  rowwise() %>% mutate(hr_seq = list(seq.POSIXt(start_hr, end_hr, by = "hour"))) %>%
  unnest(hr_seq) %>%
  group_by(hospitalization_id) %>%
  mutate(hour_idx = as.integer(difftime(hr_seq, min(hr_seq), units = "hours")) + 1L) %>%
  ungroup() %>%
  filter(hour_idx <= h_cap) %>%
  rename(hour_ts = hr_seq)

# Days grid (cap to D days)
make_day_grid <- function(df, cap_days = NA_integer_) {
  df %>%
    transmute(
      hospitalization_id,
      start = as_date(floor_date(first_icu_in,  unit = "day")),
      end   = as_date(floor_date(last_icu_out,  unit = "day"))
    ) %>%
    mutate(end = if_else(end < start, start, end)) %>%
    rowwise() %>% mutate(day = list(seq.Date(start, end, by = "day"))) %>%
    unnest(day) %>%
    group_by(hospitalization_id) %>%
    mutate(icu_day = as.integer(day - min(day)) + 1L) %>%
    { if (!is.na(cap_days)) filter(., icu_day <= cap_days) else . } %>%
    ungroup()
}
D <- 7L
icu_days <- make_day_grid(cohort_stays, cap_days = D)

# ---- Pull & standardize needed CLIF tables ----------------------------------------------------
rs_raw       <- find_tbl(clif_tables, c("respiratory_support","respiratory","vent"))
med_cont_raw <- find_tbl(clif_tables, c("medication_admin_continuous","med_admin_continuous","infusion"))
vitals_raw   <- find_tbl(clif_tables, c("vitals","vital","flowsheet"))
clif_labs    <- find_tbl(clif_tables, c("lab","labs","results"))

# ---- Respiratory support -> hourly ------------------------------------------------------------
rs_small <- rs_raw %>%
  filter(!is.na(recorded_dttm)) %>%
  mutate(
    time_raw = coalesce(as.character(recorded_dttm), as.character(recorded_time)),
    time     = parse_ts(time_raw, tz = "UTC")
  )

n_failed <- sum(is.na(rs_small$time) & !is.na(rs_small$time_raw))
message("Unparseable RS timestamps: ", n_failed)

re_ac_vc <- regex("assist.?control.*volume|ac.?vc", ignore_case = TRUE)
re_ac_pc <- regex("assist.?control.*pressure|ac.?pc", ignore_case = TRUE)
re_simv  <- fixed("simv", ignore_case = TRUE)
re_psv   <- regex("psv|pressure support", ignore_case = TRUE)
re_inv   <- regex("invasive|mechanical|vent", ignore_case = TRUE)

rs_step <- rs_small %>%
  filter(!is.na(time)) %>%
  transmute(
    hospitalization_id, time,
    device_category, mode_category, mode_name,
    artificial_airway, tracheostomy,
    fio2_set, tidal_volume_set, resp_rate_set,
    pressure_control_set, pressure_support_set, peep_set,
    peak_inspiratory_pressure_obs, plateau_pressure_obs, peep_obs, minute_vent_obs
  ) %>%
  mutate(
    hour_ts = floor_date(time, "hour"),
    device_category_lc = tolower(device_category),
    mode_category_lc   = tolower(mode_category),
    mode_name_lc       = tolower(mode_name),
    across(c(fio2_set, tidal_volume_set, resp_rate_set, pressure_control_set,
             pressure_support_set, peep_set, peak_inspiratory_pressure_obs,
             plateau_pressure_obs, peep_obs, minute_vent_obs),
           ~ suppressWarnings(as.numeric(.x))),
    artificial_airway  = as01(artificial_airway),
    tracheostomy       = as01(tracheostomy),
    fio2_set = if_else(is.na(fio2_set), NA_real_, if_else(fio2_set <= 1, fio2_set*100, fio2_set)),
    mode_ac_vc = str_detect(mode_category_lc, re_ac_vc),
    mode_ac_pc = str_detect(mode_category_lc, re_ac_pc),
    mode_simv  = str_detect(mode_category_lc, re_simv),
    mode_psv   = str_detect(mode_category_lc, re_psv),
    mode_invasive = (artificial_airway == 1L | tracheostomy == 1L | str_detect(device_category_lc, re_inv))
  ) %>%
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

rs_hr <- rs_step %>%
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

# ---- SpO2 -> hourly via DuckDB ----------------------------------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = tempfile())
copy_to(con, vitals_raw, "vitals", temporary = FALSE, overwrite = TRUE)

spo2_tbl <- tbl(con, "vitals") %>%
  transmute(
    hospitalization_id,
    recorded_dttm,
    vital_name,
    vital_category,
    vital_value
  ) %>%
  mutate(
    vname_lc   = sql("lower(vital_name)"),
    vcat_lc    = sql("lower(vital_category)"),
    looks_spo2 = sql("(regexp_matches(lower(vital_name), '(^|\\b)(spo2|sp[ _-]?o2|o2[ _-]?sat|oxygen[ _-]?saturation)($|\\b)') OR regexp_matches(lower(vital_category), 'spo2|o2|oxygen'))"),
    ts         = sql("try_cast(recorded_dttm AS TIMESTAMP)"),
    hour_ts    = sql("date_trunc('hour', ts)"),
    val        = sql("try_cast(regexp_replace(CAST(vital_value AS VARCHAR), '[^0-9\\.-]', '', 'g') AS DOUBLE)")
  ) %>%
  filter(looks_spo2, !is.na(hour_ts), !is.na(val)) %>%
  group_by(hospitalization_id, hour_ts) %>%
  summarise(spo2_med = sql('median(val)')) %>%
  ungroup()

spo2_hr <- collect(spo2_tbl)
dbDisconnect(con, shutdown = TRUE)

# ---- ABGs -> hourly ----------------------------------------------------------------------------
# Objects we definitely still need
keep_objs <- c(
  "cohort", "clif_tables",       # cohort metadata + CLIF raw tables
  "rs_hr", "spo2_hr",            # hourly respiratory + SpO₂ summaries
  "hour_grid",                   # ICU hour grid
  "con",                         # active DuckDB connection
  "spo2_tbl", "spo2_hr",         # SPO2 already processed
  "parse_ts", "as01",            # helper functions
  "re_ac_vc","re_ac_pc","re_simv","re_psv","re_inv" # regexes
)

# Drop everything else, garbage collect
rm(list = setdiff(ls(), keep_objs))
gc(full = TRUE)

clif_labs <- clif_tables[[4]]

# Keep only columns we actually need
labs_min <- as.data.table(clif_labs)[
  , .(
    hospitalization_id,
    lab_result_dttm = if ("lab_result_dttm" %in% names(clif_labs)) lab_result_dttm else NA,
    result_dttm     = if ("result_dttm"     %in% names(clif_labs)) result_dttm     else NA,
    lab_name,
    lab_value
  )
]

# Early filter to ABG-like names — this slashes the rows dramatically
labs_min[, lname_lc := tolower(as.character(lab_name))]
labs_abg <- labs_min[
  grepl("(^|\\b)ph(\\b|$)|arterial ph|blood gas ph|paco2|pa co2|arterial co2|pao2|pa o2|arterial o2",
        lname_lc, perl = TRUE)
][, lname_lc := NULL]  # drop helper column

# Free the big ones
rm(clif_labs, labs_min); gc()

# 2) Write the small ABG subset to a temp Parquet (disk-backed, zero-copy in DuckDB)
parquet_path <- file.path(tempdir(), "labs_abg_min.parquet")
arrow::write_parquet(labs_abg, parquet_path, compression = "zstd")
rm(labs_abg); gc()

# 3) Use a fresh DuckDB connection just for labs (avoid "connection busy" errors)
# 0) Turn off parallel workers so nothing touches DuckDB from other processes
if ("future" %in% loadedNamespaces()) future::plan(sequential)

# 1) Clear lazy dbplyr objects that might keep the connection busy
suppressWarnings(rm(list = ls(pattern = "(?:_tbl$)|^(spo2_tbl)$"), envir = .GlobalEnv))

# 2) Close ANY existing DuckDB connections
try({
  for (nm in ls(envir = .GlobalEnv)) {
    obj <- get(nm, envir = .GlobalEnv)
    if (inherits(obj, "duckdb_connection")) {
      try(DBI::dbDisconnect(obj, shutdown = TRUE), silent = TRUE)
    }
  }
}, silent = TRUE)
gc()

# 3) Open one fresh connection for the labs step
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tempfile())
# DBI::dbExecute(con, paste0("PRAGMA threads=", max(1L, parallel::detectCores()-1L)))
# DBI::dbExecute(con, "PRAGMA memory_limit = '2GB'")

# ---- Scan Parquet in DuckDB and compute hourly medians ----
labs_tbl <- dplyr::tbl(
  con,
  dplyr::sql(sprintf("
    SELECT
      hospitalization_id,
      date_trunc('hour',
        try_cast(
          COALESCE(CAST(lab_result_dttm AS VARCHAR), CAST(result_dttm AS VARCHAR))
          AS TIMESTAMP
        )
      ) AS hour_ts,
      lower(CAST(lab_name AS VARCHAR)) AS lname_lc,
      try_cast(regexp_replace(CAST(lab_value AS VARCHAR), '[^0-9\\.-]', '', 'g') AS DOUBLE) AS val
    FROM read_parquet('%s')
    WHERE (lab_result_dttm IS NOT NULL OR result_dttm IS NOT NULL)
  ", gsub("'", "''", parquet_path)))
) %>%
  dplyr::mutate(
    ph    = dplyr::sql("CASE WHEN regexp_matches(lname_lc, '(^|\\b)ph(\\b|$)|arterial ph|blood gas ph') THEN val END"),
    paco2 = dplyr::sql("CASE WHEN regexp_matches(lname_lc, 'paco2|pa co2|arterial co2') THEN val END"),
    pao2  = dplyr::sql("CASE WHEN regexp_matches(lname_lc, 'pao2|pa o2|arterial o2') THEN val END")
  ) %>%
  dplyr::filter(!is.na(hour_ts), !is.na(val)) %>%
  dplyr::group_by(hospitalization_id, hour_ts) %>%
  dplyr::summarise(
    ph_med    = dplyr::sql("median(ph)"),
    paco2_med = dplyr::sql("median(paco2)"),
    pao2_med  = dplyr::sql("median(pao2)")
  ) %>%
  dplyr::ungroup()

labs_hr <- dplyr::collect(labs_tbl)

# Clean up connection & temp file
DBI::dbDisconnect(con, shutdown = TRUE)
unlink(parquet_path)
gc()

# ---- Join to hourly grid & compute hourly signature --------------------------------------------
# NA-safe, within-encounter z-score helper
z_within <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    # if all NA or constant, return zeros (so it won't blow up rowMeans)
    rep(0, length(x))
  } else {
    (x - m) / s
  }
}

icu_hourly <- hour_grid %>%
  left_join(rs_hr,   by = c("hospitalization_id","hour_ts")) %>%
  left_join(spo2_hr, by = c("hospitalization_id","hour_ts")) %>%
  left_join(labs_hr, by = c("hospitalization_id","hour_ts")) %>%
  arrange(hospitalization_id, hour_ts, hour_idx) %>%
  mutate(
    fio2_max  = ifelse(!is.na(fio2_max),  pmin(pmax(fio2_max, 21), 100), NA_real_),
    peep_med  = ifelse(!is.na(peep_med),  pmin(pmax(peep_med, 0), 25),   NA_real_),
    rr_med    = ifelse(!is.na(rr_med),    pmin(pmax(rr_med, 4), 40),     NA_real_),
    vt_med    = ifelse(!is.na(vt_med),    pmin(pmax(vt_med, 200), 900),  NA_real_),
    pc_med    = ifelse(!is.na(pc_med),    pmin(pmax(pc_med, 0), 40),     NA_real_),
    ps_med    = ifelse(!is.na(ps_med),    pmin(pmax(ps_med, 0), 40),     NA_real_),
    pip_med   = ifelse(!is.na(pip_med),   pmin(pmax(pip_med, 0), 80),    NA_real_),
    pplat_med = ifelse(!is.na(pplat_med), pmin(pmax(pplat_med, 0), 60),  NA_real_),
    mv_med    = ifelse(!is.na(mv_med),    pmin(pmax(mv_med, 0), 30),     NA_real_),
    spo2_med  = ifelse(!is.na(spo2_med),  pmin(pmax(spo2_med, 50), 100), NA_real_),
    ph_med    = ifelse(!is.na(ph_med),    pmin(pmax(ph_med, 6.8), 7.8),  NA_real_),
    paco2_med = ifelse(!is.na(paco2_med), pmin(pmax(paco2_med, 20), 120),NA_real_),
    pao2_med  = ifelse(!is.na(pao2_med),  pmin(pmax(pao2_med, 30), 500), NA_real_)
  ) %>%
  mutate(
    spo2_inv  = 100 - spo2_med,                 
    fio2_frac = fio2_max / 100,
    pf_ratio  = ifelse(fio2_frac > 0, pao2_med / fio2_frac, NA_real_),
    hypox_burden = pmax(0, 92 - spo2_med)
  ) %>%
  group_by(hospitalization_id) %>%
  arrange(hour_ts, .by_group = TRUE) %>%
  mutate(
    z_imv   = z_within(as.numeric(any_imv)),
    z_fio2  = z_within(fio2_max),
    z_peep  = z_within(peep_med),
    z_rr    = z_within(rr_med),
    z_vt    = z_within(vt_med),
    z_pc    = z_within(pc_med),
    z_ps    = z_within(ps_med),
    z_spo2  = z_within(100 - spo2_med),
    z_ph    = z_within(7.40 - ph_med),
    z_paco2 = z_within(paco2_med),
    z_pao2  = z_within(80 - pao2_med)
  ) %>%
  mutate(
    sig = rowMeans(cbind(
      z_imv, z_fio2, z_peep, z_rr, z_vt, z_pc, z_ps, z_spo2, z_ph, z_paco2, z_pao2
    ), na.rm = TRUE)
  ) %>%
  filter(is.finite(sig)) %>%
  ungroup()

# ---- Aggregate to DAILY clinical signature -----------------------------------------------------
# Define a scalar horizon (and keep compatibility if other code still uses D)
horizon_days <- 7L
D <- horizon_days

icu_daily <- icu_hourly %>%
  mutate(day = as.Date(hour_ts)) %>%
  group_by(hospitalization_id, day) %>%
  summarise(
    sig_daily = if (all(is.na(sig))) NA_real_ else mean(sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(hospitalization_id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(day_idx = row_number()) %>%
  ungroup() %>%
  filter(day_idx <= horizon_days)

# ---- Build aligned DAILY exposure series per stay ----------------------------------------------

# ================================
# Monthly exposome (12m pre-admit)
# ================================

# 0) Inputs (robust renaming across common schemas)
exposome_dir <- get0("exposome_dir", ifnotfound = "exposome")

no2_mo <- readr::read_csv(file.path(exposome_dir, "no2_county_month.csv"), show_col_types = FALSE) %>%
  rename_with(~tolower(.x)) %>%
  transmute(
    county_fips = as.character(.data$county_fips),
    ym          = make_date(as.integer(year), as.integer(month), 1),
    no2         = as.double(.data$no2)
  ) %>%
  filter(!is.na(county_fips), !is.na(ym))

pm25_mo <- readr::read_csv(file.path(exposome_dir, "pm25_county_month.csv"), show_col_types = FALSE) %>%
  rename_with(~tolower(.x)) %>%
  transmute(
    county_fips = as.character(.data$county_id),
    ym          = make_date(as.integer(year), as.integer(month), 1),
    pm25        = as.double(.data$pm25)
  ) %>%
  filter(!is.na(county_fips), !is.na(ym))

# National month means + global medians for fallback
no2_nat_by_month  <- no2_mo  %>% group_by(ym) %>% summarise(no2_nat  = mean(no2,  na.rm = TRUE), .groups = "drop")
pm25_nat_by_month <- pm25_mo %>% group_by(ym) %>% summarise(pm25_nat = mean(pm25, na.rm = TRUE), .groups = "drop")
no2_overall_med   <- median(no2_mo$no2,   na.rm = TRUE)
pm25_overall_med  <- median(pm25_mo$pm25, na.rm = TRUE)

# 1) Helper: build 12-month window (months admit-12 .. admit-1), keep MONTHLY, z-score within-window
build_exposure_12m_monthly <- function(fips, admit_date) {
  end_m   <- floor_date(as.Date(admit_date) %m-% months(1), "month")
  start_m <- end_m %m-% months(11)
  mo_grid <- tibble(ym = seq.Date(start_m, end_m, by = "month"))
  
  # Join county series (may be sparse) + national means for those months
  no2_12 <- mo_grid %>%
    left_join(filter(no2_mo, county_fips == fips) %>% select(ym, no2), by = "ym") %>%
    left_join(no2_nat_by_month, by = "ym") %>%
    mutate(no2 = coalesce(no2, no2_nat, no2_overall_med)) %>%
    pull(no2)
  
  pm25_12 <- mo_grid %>%
    left_join(filter(pm25_mo, county_fips == fips) %>% select(ym, pm25), by = "ym") %>%
    left_join(pm25_nat_by_month, by = "ym") %>%
    mutate(pm25 = coalesce(pm25, pm25_nat, pm25_overall_med)) %>%
    pull(pm25)
  
  # Standardize within-window so shape dominates level
  z <- function(x) { m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - m) / s }
  no2_z  <- z(no2_12)
  pm25_z <- z(pm25_12)
  
  list(
    ym_seq    = mo_grid$ym,
    no2_vec   = no2_12,
    pm25_vec  = pm25_12,
    no2_z     = no2_z,
    pm25_z    = pm25_z,
    combo_vec = rowMeans(cbind(no2_z, pm25_z), na.rm = TRUE)  # single monthly exposure trajectory
  )
}

# 2) Encounters we’ll build exposures for (restrict to those with clinical series if you like)
# Assumes you already have `series_df_hr2` with hourly clinical trajectories (traj) ready.
# If not, swap to `distinct(cohort, hospitalization_id, county_fips, admission_dttm)`.
admit_tbl <- cohort %>%
  transmute(
    hospitalization_id,
    county_fips = as.character(county_fips),
    admit_date  = as.Date(admission_dttm)
  )

# If you want to restrict to those with usable hourly clinical series:
if (exists("series_df_hr2")) {
  admit_tbl <- admit_tbl %>% semi_join(series_df_hr2 %>% select(hospitalization_id), by = "hospitalization_id")
}

# 3) Build monthly 12m trajectories per encounter
expo_list <- admit_tbl %>%
  mutate(expo = map2(county_fips, admit_date, build_exposure_12m_monthly))

expo_series <- expo_list %>%
  transmute(
    hospitalization_id,
    expo_no2   = map(expo, "no2_z"),    # length-12 z-scored NO2
    expo_pm25  = map(expo, "pm25_z"),   # length-12 z-scored PM2.5
    expo_combo = map(expo, "combo_vec") # length-12 combined exposure trajectory
  )

stopifnot(all(c("hospitalization_id","hour_ts","sig") %in% names(icu_hourly)))

min_hours <- 12L

series_df_hr2 <- icu_hourly %>%
  arrange(hospitalization_id, hour_ts) %>%
  group_by(hospitalization_id) %>%
  summarise(
    hour_ts = list(hour_ts),
    traj    = list(as.numeric(sig)),
    n_hours = dplyr::n(),
    .groups = "drop"
  ) %>%
  filter(n_hours >= min_hours) %>%
  mutate(ok = map_lgl(traj, ~ all(is.finite(.x)))) %>%
  filter(ok) %>%
  select(-ok)

# Encounters we’ll build exposures for (restrict to those with clinical series)
admit_tbl <- cohort %>%
  transmute(
    hospitalization_id,
    county_fips = as.character(county_fips),
    admit_date  = as.Date(admission_dttm)
  ) %>%
  semi_join(series_df_hr2 %>% select(hospitalization_id), by = "hospitalization_id")

expo_series <- admit_tbl %>%
  mutate(expo = purrr::map2(county_fips, admit_date, build_exposure_12m_monthly)) %>%
  transmute(
    hospitalization_id,
    expo_vec = purrr::map(expo, "combo_vec")  # length-12 monthly exposure trajectory
  )


dat_for_dist <- series_df_hr2 %>%
  mutate(traj = purrr::map(traj, as.numeric)) %>%
  inner_join(expo_series, by = "hospitalization_id")

n <- nrow(dat_for_dist)
stopifnot(n >= 2)

ts_clin <- dat_for_dist$traj       # hourly (variable length)
ts_expo <- dat_for_dist$expo_vec   # monthly (length 12)

safe_dtw_hourly <- function(x, y, band_primary = 6L, band_fallback = 12L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = dtw::asymmetricP2, window.type = "sakoechiba",
             window.size = band_primary, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = dtw::asymmetricP2, window.type = "sakoechiba",
             window.size = band_fallback, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = dtw::asymmetricP2, distance.only = TRUE)$distance,
           error = function(e) Inf)
}

safe_dtw_monthly <- function(x, y, band_primary = 1L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = dtw::asymmetricP2, window.type = "sakoechiba",
             window.size = band_primary, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = dtw::asymmetricP2, distance.only = TRUE)$distance,
           error = function(e) Inf)
}

# Distance matrices
D_clin <- matrix(0, n, n)
D_expo <- matrix(0, n, n)
for (i in seq_len(n)) {
  for (j in i:n) {
    if (i == j) { D_clin[i,j] <- 0; D_expo[i,j] <- 0 } else {
      D_clin[i,j] <- D_clin[j,i] <- safe_dtw_hourly(ts_clin[[i]], ts_clin[[j]], 6L, 12L)
      D_expo[i,j] <- D_expo[j,i] <- safe_dtw_monthly(ts_expo[[i]], ts_expo[[j]], 1L)
    }
  }
}

# Scale & combine
scale_dist <- function(M) {
  upp <- M[upper.tri(M, diag = FALSE)]
  s <- stats::median(upp[is.finite(upp)], na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
  M / s
}
D_comb <- 0.7 * scale_dist(D_clin) + 0.3 * scale_dist(D_expo)

# PAM clustering on combined distance
k <- max(2L, min(6L, floor(sqrt(n))))
set.seed(123)
pam_fit <- cluster::pam(as.dist(D_comb), k = k, diss = TRUE)
clusters <- factor(pam_fit$clustering, levels = seq_len(k), labels = paste0("C", seq_len(k)))

# Attach labels + keep exposure vector
series_df_hr2 <- dat_for_dist %>%
  transmute(hospitalization_id, traj, expo_vec, traj_cluster_mv = clusters)




















# Hourly clinical list (assumes series_df_hr2$traj exists and is numeric)
# If you have `series_df_hr` not `series_df_hr2`, adapt the name.
stopifnot(exists("series_df_hr2"), "traj" %in% names(series_df_hr2))

dat_for_dist <- series_df_hr2 %>%
  mutate(traj = map(traj, as.numeric)) %>%
  inner_join(expo_series %>% transmute(hospitalization_id, expo_vec = expo_combo), by = "hospitalization_id")

# Safe DTW wrappers (hourly vs monthly)
safe_dtw_hourly <- function(x, y, band_primary = 6L, band_fallback = 12L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = asymmetricP2, window.type = "sakoechiba",
             window.size = band_primary, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = asymmetricP2, window.type = "sakoechiba",
             window.size = band_fallback, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = asymmetricP2, distance.only = TRUE)$distance,
           error = function(e) Inf)
}

safe_dtw_monthly <- function(x, y, band_primary = 2L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = asymmetricP2, window.type = "sakoechiba",
             window.size = band_primary, distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = asymmetricP2, distance.only = TRUE)$distance,
           error = function(e) Inf)
}

# Build distance matrices
n <- nrow(dat_for_dist)
ts_clin <- dat_for_dist$traj
ts_expo <- dat_for_dist$expo_vec

D_clin <- matrix(0, n, n)
D_expo <- matrix(0, n, n)

for (i in seq_len(n)) {
  for (j in i:n) {
    if (i == j) {
      D_clin[i, j] <- 0; D_expo[i, j] <- 0
    } else {
      d1 <- safe_dtw_hourly(ts_clin[[i]], ts_clin[[j]], band_primary = 6L, band_fallback = 12L)
      d2 <- safe_dtw_monthly(ts_expo[[i]], ts_expo[[j]], band_primary = 1L)
      D_clin[i, j] <- D_clin[j, i] <- d1
      D_expo[i, j] <- D_expo[j, i] <- d2
    }
  }
}

# Robust scale each matrix so they’re comparable, then combine
scale_dist <- function(M) {
  upp <- M[upper.tri(M, diag = FALSE)]
  s <- stats::median(upp[is.finite(upp)], na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- 1
  M / s
}
D_clin_s <- scale_dist(D_clin)
D_expo_s <- scale_dist(D_expo)

w_clin <- 0.7; w_expo <- 0.3  # weights; tweak as needed
D_comb <- w_clin * D_clin_s + w_expo * D_expo_s

# Cluster with PAM on the combined distance
k <- max(2L, min(6L, floor(sqrt(n))))  # simple heuristic
pam_fit <- cluster::pam(as.dist(D_comb), k = k, diss = TRUE)
clusters <- factor(pam_fit$clustering, levels = seq_len(k), labels = paste0("C", seq_len(k)))

# Attach labels back to encounters
series_df_hr2 <- dat_for_dist %>% mutate(traj_cluster_mv = clusters)

# (Optional) Keep exposure vectors alongside for later summaries
series_df_hr2 <- series_df_hr2 %>% select(hospitalization_id, traj, expo_vec, traj_cluster_mv)

# Now you can proceed to join labels to cohort, fit membership models, etc.



































clin_daily_series <- icu_daily %>%
  group_by(hospitalization_id) %>%
  summarise(
    clin_vec = list(z_within(sig_daily)),
    n = length(clin_vec[[1]]),
    .groups = "drop"
  ) %>%
  filter(n >= 2) %>% select(-n)

mv_series <- clin_daily_series %>%
  inner_join(expo_daily_series, by = "hospitalization_id") %>%
  mutate(
    len_ok = map_int(clin_vec, length) >= D & map_int(expo_combo, length) >= D
  ) %>%
  filter(len_ok) %>%
  mutate(
    clin_vec = map(clin_vec, ~ .x[seq_len(D)]),
    expo_vec = map(expo_combo, ~ .x[seq_len(D)]),
    mv_mat   = map2(clin_vec, expo_vec, ~ rbind(.x, .y))
  ) %>%
  select(hospitalization_id, mv_mat)

stopifnot(nrow(mv_series) > 1)

# ---- Multivariate DTW clustering (bounded warp, parallel) -------------------------------------
future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))

ts_mv <- mv_series$mv_mat
k <- max(2L, min(6L, floor(sqrt(length(ts_mv)))))  # heuristic

set.seed(123)
cl_mv <- dtwclust::tsclust(
  tslist   = ts_mv,
  type     = "partitional",
  k        = k,
  distance = "dtw_basic",
  centroid = "dba",
  seed     = 123,
  trace    = TRUE,
  args     = dtwclust::tsclust_args(dist = list(window.size = 1)),  # ±1 day warp
  control  = dtwclust::partitional_control(iter.max = 100L, nrep = 3L)
)

mv_series$traj_cluster_mv <- factor(dtwclust::partition(cl_mv),
                                    levels = seq_len(k), labels = paste0("C", seq_len(k)))
series_labels_mv <- mv_series %>% select(hospitalization_id, traj_cluster_mv)

cat("Cluster sizes (MV-DTW):\n")
print(series_labels_mv %>% count(traj_cluster_mv))

# ---- Prototypes (medoids) per cluster ---------------------------------------------------------
prototypes_mv <- mv_series %>%
  left_join(series_labels_mv, by = "hospitalization_id") %>%
  group_by(traj_cluster_mv) %>%
  summarise(
    midx = get_medoid_idx_mv(mv_mat, safe_dtw_daily),
    clin_proto = list(mv_mat[[midx]][1, ]),
    expo_proto = list(mv_mat[[midx]][2, ]),
    .groups = "drop"
  )

# ---- Membership model (multinomial) -----------------------------------------------------------
# Simple demo + exposure covariates
expo_covars <- c("no2_mean_3y", "pm25_mean_3y", "svi_last", "tmax_mean_3y", "vp_mean_3y", "no2_mean_12mo")
demo_covars <- c("age_years", "sex_category", "race_ethnicity_simple",
                 "hypoxemic_arf", "hypercapnic_arf", "mixed_arf")

# Recode race/ethnicity from cohort
analytic_mv <- cohort %>%
  select(patient_id, hospitalization_id, hypoxemic_arf, hypercapnic_arf, mixed_arf,
         age_years, sex_category, race_category, ethnicity_category, county_fips, admission_dttm) %>%
  left_join(series_labels_mv, by = "hospitalization_id") %>%
  left_join(cohort_expo, by = c("patient_id","hospitalization_id")) %>%
  filter(!is.na(traj_cluster_mv)) %>%
  mutate(
    race_name = race_category,
    re_low = str_to_lower(paste(ethnicity_category, race_name)),
    is_nonhisp = str_detect(re_low, "\\bnon[- ]?hispanic\\b"),
    is_hisp    = str_detect(re_low, "\\bhispanic\\b") & !is_nonhisp,
    is_white   = str_detect(re_low, "white"),
    is_black   = str_detect(re_low, "black"),
    is_asian_any = str_detect(re_low, "asian|mideast|filipino|chinese|korean|vietnamese|pacific islander|samoan"),
    race_ethnicity_simple = case_when(
      is_white & is_hisp    ~ "Hispanic White",
      is_white & is_nonhisp ~ "Non-Hispanic White",
      is_black & is_hisp    ~ "Hispanic Black",
      is_black & is_nonhisp ~ "Non-Hispanic Black",
      is_asian_any          ~ "Asian",
      TRUE                  ~ "Other"
    ),
    sex_category = factor(sex_category),
    race_ethnicity_simple = factor(
      race_ethnicity_simple,
      levels = c("Non-Hispanic White", "Hispanic White",
                 "Non-Hispanic Black", "Hispanic Black",
                 "Asian", "Other")
    )
  ) %>%
  select(-re_low, -is_nonhisp, -is_hisp, -is_white, -is_black, -is_asian_any)

# Build formula and fit
has_cov <- intersect(c(expo_covars, demo_covars), names(analytic_mv))
formula_mn <- as.formula(paste("traj_cluster_mv ~", paste(has_cov, collapse = " + ")))

set.seed(123)
fit_mn <- nnet::multinom(formula_mn, data = analytic_mv, trace = FALSE)
coef_mn <- broom::tidy(fit_mn, exponentiate = TRUE, conf.int = TRUE)

cat("\nMultinomial membership model (RRRs):\n")
print(coef_mn)

# ---- Quick summaries/plots --------------------------------------------------------------------
# Cluster counts
p_counts <- analytic_mv %>%
  count(traj_cluster_mv) %>%
  ggplot(aes(traj_cluster_mv, n)) +
  geom_col() +
  labs(x = "Trajectory class (MV-DTW)", y = "N", title = "Cluster sizes") +
  theme_minimal(base_size = 13)

print(p_counts)

# Exposure by cluster (12-mo NO2 if present)
if ("no2_mean_12mo" %in% names(analytic_mv)) {
  p_no2 <- analytic_mv %>%
    ggplot(aes(traj_cluster_mv, no2_mean_12mo)) +
    geom_boxplot() +
    labs(x = "Trajectory class (MV-DTW)", y = "12-mo mean NO\u2082 (ppb)") +
    theme_minimal(base_size = 13)
  print(p_no2)
}

# Clean daily profiles per cluster for inspection (clinical channel only)
daily_profiles <- icu_daily %>%
  left_join(series_labels_mv, by = "hospitalization_id") %>%
  group_by(traj_cluster_mv, day_idx) %>%
  summarise(clin_mean = mean(sig_daily, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

p_prof <- daily_profiles %>%
  ggplot(aes(day_idx, clin_mean, color = traj_cluster_mv)) +
  geom_line(size = 1) +
  labs(x = "ICU Day", y = "Mean daily clinical severity (z)",
       title = "Clinical trajectory profiles by cluster") +
  theme_minimal(base_size = 13)
print(p_prof)

# ---- Return key objects in environment --------------------------------------------------------
# mv_series: per-encounter 2xD matrices
# series_labels_mv: cluster labels per hospitalization_id
# prototypes_mv: medoid-based prototypes per cluster
# analytic_mv: covariates joined with cluster labels
invisible(TRUE)
