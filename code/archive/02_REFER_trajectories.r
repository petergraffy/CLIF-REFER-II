# ===============================================================================================
# ICU Respiratory Failure Environmental Risk (REFER-II) — Full Pipeline
# PI: Peter Graffy (graffy@uchicago.edu)
# Purpose: Build daily clinical trajectories, daily NO2/PM2.5 exposures via monthly→daily
#          interpolation, multivariate DTW clustering (exposure + clinical), and labels.
# Inputs:  cohort_min, clif_tables; exposome_dir with pm25_county_month.csv, no2_county_month.csv
# ===============================================================================================

# ---- Libraries -------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate); library(readr); library(stringr)
  library(purrr); library(data.table); library(DBI); library(duckdb); library(arrow)
  library(zoo); library(tibble); library(cluster); library(dtw); library(forcats)
  library(ragg); library(svglite)
})

# Load configuration utility
source("utils/config.R")
repo <- config$repo
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type

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
      orders = c("ymd HMS","ymd HM","ymd H","Ymd HMS","Ymd HM","Ymd H",
                 "mdy HMS","mdy HM","mdy H","dmy HMS","dmy HM","dmy H"),
      tz = tz, truncated = 3
    )) |> as.numeric()
  }
  as.POSIXct(out, origin = "1970-01-01", tz = tz)
}
safe_max <- function(x) { v <- x[is.finite(x)]; if (length(v) == 0) NA_real_ else max(v) }
safe_med <- function(x) { v <- x[is.finite(x)]; if (length(v) == 0) NA_real_ else stats::median(v) }
z_win    <- function(x) { m <- mean(x, na.rm=TRUE); s <- sd(x, na.rm=TRUE); if (!is.finite(s) || s==0) rep(0,length(x)) else (x-m)/s }
clamp    <- function(x, lo, hi) ifelse(is.na(x), NA_real_, pmin(pmax(x, lo), hi))

# Monthly -> Daily interpolation over [start_date, end_date]
interp_monthly_to_daily <- function(df_mo, value_col, start_date, end_date) {
  mo_key <- df_mo %>% arrange(ym)
  daily <- tibble(date = seq.Date(start_date, end_date, by = "day"))
  if (nrow(mo_key) == 0) { daily$val <- NA_real_; return(daily) }
  anchors <- bind_rows(
    mo_key %>% transmute(date = ym, val = .data[[value_col]]),
    tibble(date = seq.Date(max(mo_key$ym), by = "month", length.out = 2)[2],
           val  = tail(mo_key[[value_col]], 1))
  ) %>% arrange(date)
  daily$val <- approx(
    x = as.numeric(anchors$date), y = anchors$val,
    xout = as.numeric(daily$date),
    method = "linear", rule = 2, ties = "ordered"
  )$y
  daily
}

# DTW helpers (daily)
safe_dtw_daily <- function(x, y, band_primary = 2L, band_fallback = 5L) {
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = dtw::asymmetricP2,
             window.type = "sakoechiba", window.size = band_primary,
             distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  out <- tryCatch(
    dtw::dtw(x, y, step.pattern = dtw::asymmetricP2,
             window.type = "sakoechiba", window.size = band_fallback,
             distance.only = TRUE)$distance,
    error = function(e) NA_real_
  )
  if (is.finite(out)) return(out)
  tryCatch(dtw::dtw(x, y, step.pattern = dtw::asymmetricP2,
                    distance.only = TRUE)$distance,
           error = function(e) Inf)
}
scale_dist <- function(M) {
  upp <- M[upper.tri(M, diag = FALSE)]
  s <- stats::median(upp[is.finite(upp)], na.rm = TRUE); if (!is.finite(s) || s==0) s <- 1
  M / s
}

# ---- Cohort & ICU windows ---------------------------------------------------------------------
cohort <- cohort_min %>%
  mutate(
    admit_dt     = as_date(admission_dttm),
    admit_year   = year(admit_dt),
    county_fips  = county_code,
    first_icu_in = as_datetime(first_icu_in, tz = "UTC"),
    last_icu_out = as_datetime(last_icu_out, tz = "UTC")
  ) %>%
  filter(!is.na(first_icu_in), !is.na(last_icu_out), last_icu_out > first_icu_in)

# Hourly ICU grid (cap if desired; here unlimited)
cohort_stays_hr <- cohort %>%
  transmute(hospitalization_id,
            start_hr = floor_date(first_icu_in, "hour"),
            end_hr   = floor_date(last_icu_out, "hour"))
hour_grid <- cohort_stays_hr %>%
  rowwise() %>% mutate(hr_seq = list(seq.POSIXt(start_hr, end_hr, by = "hour"))) %>%
  unnest(hr_seq) %>%
  group_by(hospitalization_id) %>%
  mutate(hour_idx = as.integer(difftime(hr_seq, min(hr_seq), units = "hours")) + 1L) %>%
  ungroup() %>%
  rename(hour_ts = hr_seq)

# ---- Pull & standardize needed CLIF tables ----------------------------------------------------
rs_raw       <- { # respiratory support
  idx <- which(Reduce(`|`, lapply(c("respiratory_support","respiratory","vent"), function(p)
    str_detect(tolower(names(clif_tables)), tolower(p)))))
  clif_tables[[idx[1]]]
}
vitals_raw   <- { # vitals/flowsheets
  idx <- which(Reduce(`|`, lapply(c("vitals","vital","flowsheet"), function(p)
    str_detect(tolower(names(clif_tables)), tolower(p)))))
  clif_tables[[idx[1]]]
}
clif_labs    <- { # labs/results
  idx <- which(Reduce(`|`, lapply(c("lab","labs","results"), function(p)
    str_detect(tolower(names(clif_tables)), tolower(p)))))
  clif_tables[[idx[1]]]
}

# ---- Respiratory support -> hourly summaries ---------------------------------------------------
rs_small <- rs_raw %>%
  filter(!is.na(recorded_dttm)) %>%
  mutate(
    time_raw = as.character(recorded_dttm),
    time     = parse_ts(time_raw, tz = "UTC")
  ) %>% filter(!is.na(time))

re_ac_vc <- regex("assist.?control.*volume|ac.?vc", ignore_case = TRUE)
re_ac_pc <- regex("assist.?control.*pressure|ac.?pc", ignore_case = TRUE)
re_simv  <- fixed("simv", ignore_case = TRUE)
re_psv   <- regex("psv|pressure support", ignore_case = TRUE)
re_inv   <- regex("IMV", ignore_case = TRUE)

rs_step <- rs_small %>%
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
# con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tempfile())
# DBI::dbWriteTable(con, "vitals", vitals_raw, overwrite = TRUE)
# 
# spo2_tbl <- dplyr::tbl(con, "vitals") %>%
#   transmute(
#     hospitalization_id,
#     recorded_dttm,
#     vital_name,
#     vital_category,
#     vital_value
#   ) %>%
#   mutate(
#     vname_lc   = sql("lower(vital_name)"),
#     vcat_lc    = sql("lower(vital_category)"),
#     looks_spo2 = sql("regexp_matches(lower(vital_name), 'sp[o0]2|oxygen saturation|oximetry')"),
#     ts         = sql("try_cast(recorded_dttm AS TIMESTAMP)"),
#     hour_ts    = sql("date_trunc('hour', ts)"),
#     val        = sql("try_cast(regexp_replace(CAST(vital_value AS VARCHAR), '[^0-9\\.-]', '', 'g') AS DOUBLE)")
#   ) %>%
#   filter(looks_spo2, !is.na(hour_ts), !is.na(val)) %>%
#   group_by(hospitalization_id, hour_ts) %>%
#   summarise(spo2_med = sql('median(val)')) %>%
#   ungroup()
# 
# spo2_hr <- collect(spo2_tbl)
# DBI::dbDisconnect(con, shutdown = TRUE)

# ================================================
# Read ABGs (and optional SpO2) directly from Parquet
# ================================================

stopifnot(!is.null(config$tables_path), dir.exists(config$tables_path))

# -- Tidy up any old DuckDB state to avoid "connection working on another query"
# if ("future" %in% loadedNamespaces()) future::plan(sequential)
# try({
#   objs <- ls(.GlobalEnv, all.names = TRUE)
#   for (nm in objs) {
#     obj <- get(nm, envir = .GlobalEnv)
#     if (inherits(obj, "tbl_sql")) rm(list = nm, envir = .GlobalEnv)
#     if (inherits(obj, "duckdb_connection")) DBI::dbDisconnect(obj, shutdown = TRUE)
#   }
# }, silent = TRUE); gc()

# --------------------------
# Connect & set pragmas
# --------------------------
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tempfile())
DBI::dbExecute(con, sprintf("PRAGMA threads=%d", max(1L, parallel::detectCores()-1L)))
# tune as needed; lower if RAM is tight
DBI::dbExecute(con, "PRAGMA memory_limit='6GB'")
# optional: spill temp files here
DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", tempdir()))

# ==========================================================
# A) ABGs -> hourly medians (fully in DuckDB from Parquet)
#    Adjust 'labs_glob' to your actual Parquet names.
# ==========================================================
labs_glob <- file.path(config$tables_path, "clif_labs*.parquet")  # e.g., "labs.parquet" or "labs_*.parquet"

labs_hr <- DBI::dbGetQuery(con, sprintf("
WITH base AS (
  SELECT
    hospitalization_id,
    date_trunc('hour',
      try_cast(
        CAST(lab_result_dttm AS VARCHAR)
        AS TIMESTAMP
      )
    ) AS hour_ts,
    lower(CAST(lab_name AS VARCHAR)) AS lname_lc,
    try_cast(regexp_replace(CAST(lab_value AS VARCHAR), '[^0-9\\.-]', '', 'g') AS DOUBLE) AS val
  FROM read_parquet('%s')
  WHERE lab_result_dttm IS NOT NULL
),
picked AS (
  SELECT hospitalization_id, hour_ts,
    CASE WHEN regexp_matches(lname_lc, '(^|\\b)ph(\\b|$)|arterial ph|blood gas ph') THEN val END AS ph,
    CASE WHEN regexp_matches(lname_lc, 'paco2|pa co2|arterial co2') THEN val END AS paco2,
    CASE WHEN regexp_matches(lname_lc, 'pao2|pa o2|arterial o2') THEN val END AS pao2
  FROM base
  WHERE hour_ts IS NOT NULL AND val IS NOT NULL
)
SELECT hospitalization_id, hour_ts,
       median(ph)    AS ph_med,
       median(paco2) AS paco2_med,
       median(pao2)  AS pao2_med
FROM picked
GROUP BY 1,2
ORDER BY 1,2
", gsub("'", "''", labs_glob)))

# ==========================================================
# B) (Optional) SpO2 directly from Parquet as well
#    Skip if you already built spo2_hr earlier.
# ==========================================================
# Change the glob to your vitals parquet files (flowsheets, vitals, etc.).
vitals_glob <- file.path(config$tables_path, "clif_vitals*.parquet")  # e.g., "flowsheets_*.parquet"

if (length(Sys.glob(vitals_glob)) > 0) {
  spo2_hr <- DBI::dbGetQuery(con, sprintf("
  WITH base AS (
    SELECT
      hospitalization_id,
      date_trunc('hour', try_cast(recorded_dttm AS TIMESTAMP)) AS hour_ts,
      lower(CAST(vital_name AS VARCHAR))     AS vname_lc,
      lower(CAST(vital_category AS VARCHAR)) AS vcat_lc,
      try_cast(
        regexp_replace(CAST(vital_value AS VARCHAR), '[^0-9\\.-]', '', 'g'
      ) AS DOUBLE) AS val
    FROM read_parquet('%s')
  ),
  picked AS (
    SELECT hospitalization_id, hour_ts, val
    FROM base
    WHERE hour_ts IS NOT NULL
      AND val IS NOT NULL
      AND (
        regexp_matches(vname_lc, 'sp[o0]2|oxygen saturation|oximetry')
        OR regexp_matches(vcat_lc, 'oximetry|pulse ox')
      )
  )
  SELECT hospitalization_id, hour_ts, median(val) AS spo2_med
  FROM picked
  GROUP BY 1,2
  ORDER BY 1,2
  ", gsub("'", "''", vitals_glob)))
}

# All done with DuckDB for these steps
DBI::dbDisconnect(con, shutdown = TRUE); gc()

# ---- Join hourly + compute clinical signature -------------------------------------------------
icu_hourly <- hour_grid %>%
  left_join(rs_hr,   by = c("hospitalization_id","hour_ts")) %>%
  left_join(spo2_hr, by = c("hospitalization_id","hour_ts")) %>%
  left_join(labs_hr, by = c("hospitalization_id","hour_ts")) %>%
  arrange(hospitalization_id, hour_ts, hour_idx) %>%
  mutate(
    fio2_max  = clamp(fio2_max, 21, 100),
    peep_med  = clamp(peep_med, 0, 25),
    rr_med    = clamp(rr_med, 4, 40),
    vt_med    = clamp(vt_med, 200, 900),
    pc_med    = clamp(pc_med, 0, 40),
    ps_med    = clamp(ps_med, 0, 40),
    pip_med   = clamp(pip_med, 0, 80),
    pplat_med = clamp(pplat_med, 0, 60),
    mv_med    = clamp(mv_med, 0, 30),
    spo2_med  = clamp(spo2_med, 50, 100),
    ph_med    = clamp(ph_med, 6.8, 7.8),
    paco2_med = clamp(paco2_med, 20, 120),
    pao2_med  = clamp(pao2_med, 30, 500)
  ) %>%
  mutate(
    spo2_inv  = 100 - spo2_med,
    fio2_frac = fio2_max / 100,
    pf_ratio  = ifelse(fio2_frac > 0, pao2_med / fio2_frac, NA_real_)
  ) %>%
  group_by(hospitalization_id) %>%
  arrange(hour_ts, .by_group = TRUE) %>%
  mutate(
    z_imv   = z_win(as.numeric(any_imv)),
    z_fio2  = z_win(fio2_max),
    z_peep  = z_win(peep_med),
    z_rr    = z_win(rr_med),
    z_vt    = z_win(vt_med),
    z_pc    = z_win(pc_med),
    z_ps    = z_win(ps_med),
    z_spo2  = z_win(100 - spo2_med),
    z_ph    = z_win(7.40 - ph_med),
    z_paco2 = z_win(paco2_med),
    z_pao2  = z_win(80 - pao2_med),
    sig     = rowMeans(cbind(z_imv, z_fio2, z_peep, z_rr, z_vt, z_pc, z_ps, z_spo2, z_ph, z_paco2, z_pao2), na.rm = TRUE)
  ) %>%
  filter(is.finite(sig)) %>%
  ungroup()

# ---- Aggregate to DAILY clinical signature -----------------------------------------------------
icu_daily_full <- icu_hourly %>%
  mutate(day = as.Date(hour_ts)) %>%
  group_by(hospitalization_id, day) %>%
  summarise(sig_daily = if (all(is.na(sig))) NA_real_ else mean(sig, na.rm = TRUE), .groups = "drop") %>%
  group_by(hospitalization_id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(day_idx = row_number()) %>%
  ungroup()

build_clinical_matrix <- function(df_daily, horizons = c(7L, 14L, 30L)) {
  mats <- list()
  for (H in horizons) {
    mat <- df_daily %>%
      filter(day_idx <= H) %>%
      group_by(hospitalization_id) %>%
      summarise(vec = list(if (n() < H) c(sig_daily, rep(NA_real_, H - n())) else sig_daily[1:H]), .groups = "drop") %>%
      mutate(vec = map(vec, ~ z_win(.x)))
    V <- do.call(rbind, mat$vec)
    rownames(V) <- mat$hospitalization_id
    mats[[paste0("clin_", H, "d_z")]] <- V
  }
  mats
}
clin_mats <- build_clinical_matrix(icu_daily_full, horizons = c(7L, 14L, 30L))
X_clin_7d_z  <- clin_mats$clin_7d_z
X_clin_14d_z <- clin_mats$clin_14d_z
X_clin_30d_z <- clin_mats$clin_30d_z

# ---- Exposome inputs (monthly) ----------------------------------------------------------------
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

no2_nat_by_month  <- no2_mo  %>% group_by(ym) %>% summarise(no2_nat  = mean(no2,  na.rm = TRUE), .groups = "drop")
pm25_nat_by_month <- pm25_mo %>% group_by(ym) %>% summarise(pm25_nat = mean(pm25, na.rm = TRUE), .groups = "drop")
no2_overall_med   <- median(no2_mo$no2,   na.rm = TRUE)
pm25_overall_med  <- median(pm25_mo$pm25, na.rm = TRUE)

# ---- Build daily exposure windows (90d & 365d pre-admit) --------------------------------------
build_exposure_daily_window <- function(fips, admit_date, window_days = 90) {
  end_m   = floor_date(as.Date(admit_date) %m-% months(1), "month")
  start_m = floor_date((as.Date(admit_date) - days(window_days)), "month")
  mo_grid = tibble(ym = seq.Date(start_m, end_m, by = "month"))
  
  no2_series <- mo_grid %>%
    left_join(filter(no2_mo, county_fips == fips) %>% select(ym, no2), by="ym") %>%
    left_join(no2_nat_by_month, by="ym") %>%
    mutate(no2 = coalesce(no2, no2_nat, no2_overall_med)) %>%
    select(ym, no2)
  
  pm25_series <- mo_grid %>%
    left_join(filter(pm25_mo, county_fips == fips) %>% select(ym, pm25), by="ym") %>%
    left_join(pm25_nat_by_month, by="ym") %>%
    mutate(pm25 = coalesce(pm25, pm25_nat, pm25_overall_med)) %>%
    select(ym, pm25)
  
  daily_grid <- tibble(date = seq.Date(as.Date(admit_date) - days(window_days),
                                       as.Date(admit_date) - days(1), by = "day"))
  
  no2_daily <- interp_monthly_to_daily(
    df_mo = no2_series %>% rename(val = no2),
    value_col = "val",
    start_date = min(daily_grid$date),
    end_date   = max(daily_grid$date)
  ) %>% pull(val)
  
  pm25_daily <- interp_monthly_to_daily(
    df_mo = pm25_series %>% rename(val = pm25),
    value_col = "val",
    start_date = min(daily_grid$date),
    end_date   = max(daily_grid$date)
  ) %>% pull(val)
  
  no2_z   <- z_win(no2_daily)
  pm25_z  <- z_win(pm25_daily)
  combo_z <- rowMeans(cbind(no2_z, pm25_z), na.rm = TRUE)
  
  list(
    dates    = daily_grid$date,
    no2_z    = no2_z,
    pm25_z   = pm25_z,
    combo_z  = combo_z,
    no2_raw  = no2_daily,
    pm25_raw = pm25_daily
  )
}

admit_tbl <- cohort %>%
  transmute(
    hospitalization_id,
    county_fips = as.character(county_fips),
    admit_date  = as.Date(admission_dttm)
  )

expo_daily_90  <- admit_tbl %>% mutate(expo = map2(county_fips, admit_date, \(f,d) build_exposure_daily_window(f,d,90)))
expo_daily_365 <- admit_tbl %>% mutate(expo = map2(county_fips, admit_date, \(f,d) build_exposure_daily_window(f,d,365)))

list_to_mat <- function(lst, field) {
  vecs <- map(lst$expo, field)
  maxlen <- max(lengths(vecs))
  stopifnot(all(lengths(vecs) == maxlen))
  M <- do.call(rbind, vecs)
  rownames(M) <- lst$hospitalization_id
  M
}

X_no2_90_z    <- list_to_mat(expo_daily_90,  "no2_z")
X_pm25_90_z   <- list_to_mat(expo_daily_90,  "pm25_z")
X_combo_90_z  <- list_to_mat(expo_daily_90,  "combo_z")

X_no2_365_z   <- list_to_mat(expo_daily_365, "no2_z")
X_pm25_365_z  <- list_to_mat(expo_daily_365, "pm25_z")
X_combo_365_z <- list_to_mat(expo_daily_365, "combo_z")

# ---- DTW: Combine clinical & exposure windows flexibly ----------------------------------------
# 1) Optional: downsample exposures with PAA for speed
clean_vec <- function(x) {
  x <- as.numeric(x)
  # interpolate interior NAs
  if (all(is.na(x))) return(rep(0, length(x)))
  xi <- suppressWarnings(na.approx(x, na.rm = FALSE))
  # carry edges
  xi <- na.locf(xi, na.rm = FALSE)
  xi <- na.locf(xi, fromLast = TRUE, na.rm = FALSE)
  # any stubborn NA (e.g., all NA originally) -> 0 (neutral after z-scoring)
  xi[!is.finite(xi)] <- 0
  xi[is.na(xi)] <- 0
  xi
}

paa <- function(x, m) {
  x <- clean_vec(x)
  n <- length(x)
  if (m >= n) return(x)
  # map indices to m bins
  idx <- floor(((seq_len(n) - 0.5) * m) / n) + 1L
  as.numeric(tapply(x, idx, function(v) mean(v, na.rm = TRUE)))
}

downsample_rows <- function(M, m) {
  ids <- rownames(M)
  V <- lapply(seq_len(nrow(M)), function(i) paa(M[i, ], m))
  V <- do.call(rbind, V)
  rownames(V) <- ids
  V
}

# ---- Very robust DTW distance ----
dtw_dist_robust <- function(x, y, band_hint = 2L) {
  x <- clean_vec(x); y <- clean_vec(y)
  nx <- length(x); ny <- length(y)
  dlen <- abs(nx - ny)
  # band must be >= length difference for Sakoe–Chiba to allow any path
  band1 <- max(band_hint, dlen)
  try_seq <- list(
    list(win = "sakoechiba", band = band1, step = dtw::asymmetricP2, open = c(FALSE,FALSE)),
    list(win = "sakoechiba", band = max(band1, floor(0.15 * max(nx,ny))), step = dtw::asymmetricP2, open = c(FALSE,FALSE)),
    list(win = NULL,         band = NA, step = dtw::asymmetricP2, open = c(FALSE,FALSE)),
    list(win = NULL,         band = NA, step = dtw::asymmetricP2, open = c(TRUE, TRUE)),
    list(win = NULL,         band = NA, step = dtw::symmetric2,   open = c(TRUE, TRUE))
  )
  for (cfg in try_seq) {
    res <- try({
      dtw::dtw(
        x, y,
        step.pattern = cfg$step,
        window.type  = cfg$win,
        window.size  = if (is.null(cfg$win)) NULL else cfg$band,
        open.end     = cfg$open[2],
        open.begin   = cfg$open[1],
        distance.only = TRUE
      )$distance
    }, silent = TRUE)
    if (is.numeric(res) && is.finite(res)) return(res)
  }
  Inf
}

# Choose matrices (same as before)
A_mat <- X_clin_14d_z                      # clinical (equal length)
B_mat <- X_combo_90_z                      # exposure (longer)
B_mat_ds <- downsample_rows(B_mat, m = 30) # 90d -> 30 bins (use 60 for 365d)

build_medoids_bank <- function(A_mat, B_mat_ds, wA = 0.7, wB = 0.3, k = NULL, subsample = 500,
                               seed = 123L, bandA = 2L, bandB = 2L) {
  set.seed(seed)
  ids_all <- intersect(rownames(A_mat), rownames(B_mat_ds))
  s <- min(subsample, length(ids_all))
  ids_sub <- sample(ids_all, s)
  A_sub <- A_mat[ids_sub, , drop = FALSE]
  B_sub <- B_mat_ds[ids_sub, , drop = FALSE]
  
  sN <- length(ids_sub)
  D_A <- matrix(0, sN, sN)
  D_B <- matrix(0, sN, sN)
  
  for (i in seq_len(sN)) {
    xiA <- clean_vec(A_sub[i, ])
    xiB <- clean_vec(B_sub[i, ])
    for (j in i:sN) {
      if (i == j) { D_A[i,j] <- 0; D_B[i,j] <- 0 } else {
        d1 <- dtw_dist_robust(xiA, clean_vec(A_sub[j, ]), band_hint = bandA)
        d2 <- dtw_dist_robust(xiB, clean_vec(B_sub[j, ]), band_hint = bandB)
        D_A[i,j] <- D_A[j,i] <- d1
        D_B[i,j] <- D_B[j,i] <- d2
      }
    }
  }
  
  scale_dist <- function(M) {
    upp <- M[upper.tri(M, diag = FALSE)]
    s <- stats::median(upp[is.finite(upp)], na.rm = TRUE); if (!is.finite(s) || s==0) s <- 1
    M / s
  }
  D_comb <- wA * scale_dist(D_A) + wB * scale_dist(D_B)
  
  if (is.null(k)) k <- max(2L, min(6L, floor(sqrt(sN))))
  pam_fit <- cluster::pam(stats::as.dist(D_comb), k = k, diss = TRUE)
  list(ids_sub = ids_sub, medoid_ids = ids_sub[pam_fit$id.med], D_sub = D_comb)
}

assign_to_medoids <- function(A_mat, B_mat_ds, medoid_ids, wA = 0.7, wB = 0.3,
                              bandA = 2L, bandB = 2L) {
  ids_all <- intersect(rownames(A_mat), rownames(B_mat_ds))
  A <- A_mat[ids_all, , drop = FALSE]
  B <- B_mat_ds[ids_all, , drop = FALSE]
  
  A_med <- lapply(medoid_ids, function(id) clean_vec(A_mat[id, ]))
  B_med <- lapply(medoid_ids, function(id) clean_vec(B_mat_ds[id, ]))
  
  DA_proj <- matrix(NA_real_, nrow = length(ids_all), ncol = length(medoid_ids),
                    dimnames = list(ids_all, medoid_ids))
  DB_proj <- DA_proj
  
  for (i in seq_along(ids_all)) {
    xiA <- clean_vec(A[i, ])
    xiB <- clean_vec(B[i, ])
    for (m in seq_along(medoid_ids)) {
      DA_proj[i, m] <- dtw_dist_robust(xiA, A_med[[m]], band_hint = bandA)
      DB_proj[i, m] <- dtw_dist_robust(xiB, B_med[[m]], band_hint = bandB)
    }
  }
  
  # Row-wise robust scaling to balance rows with differing absolute DTWs
  rob_scale <- function(M) {
    s <- apply(M, 1, function(r) median(r[is.finite(r)], na.rm = TRUE))
    s[!is.finite(s) | s <= 0] <- 1
    M / s
  }
  DA_s <- rob_scale(DA_proj); DB_s <- rob_scale(DB_proj)
  D_comb_proj <- wA * DA_s + wB * DB_s
  
  nearest_idx <- apply(D_comb_proj, 1, which.min)
  tibble::tibble(
    hospitalization_id = ids_all,
    traj_cluster_mv = factor(nearest_idx, labels = paste0("C", seq_len(length(medoid_ids)))),
    nearest_medoid = medoid_ids[nearest_idx]
  )
}

# ---- Run it ----
set.seed(123)
bank <- build_medoids_bank(A_mat, B_mat_ds, wA = 0.7, wB = 0.3, k = NULL, subsample = 500,
                           seed = 123, bandA = 2L, bandB = 2L)
labels_fast <- assign_to_medoids(A_mat, B_mat_ds, bank$medoid_ids, wA = 0.7, wB = 0.3,
                                 bandA = 2L, bandB = 2L)


# ----- Join cluster labels to cohort -----
stopifnot(all(c("hospitalization_id","traj_cluster_mv") %in% names(labels_fast)))
prof_df <- cohort %>%
  select(hospitalization_id, patient_id, admit_dt, county_fips,
         tidyselect::any_of(c(
           # put your canonical fields here if you have them
           "age","sex","sex_category","race","race_ethnicity_simple","svi_overall",
           "acs_median_income","acs_pct_lt_hs","acs_pct_insured","acs_unemp_rate_pct",
           # outcomes if available
           "in_hosp_death","death_30d","icu_los_days","mech_vent_days"
         ))) %>%
  right_join(labels_fast, by = "hospitalization_id") %>%
  mutate(traj_cluster_mv = forcats::fct_inorder(traj_cluster_mv))

# ----- Helpers for robust summaries -----
summ_num <- function(x) {
  c(n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    p25 = quantile(x, 0.25, na.rm = TRUE),
    p50 = quantile(x, 0.50, na.rm = TRUE),
    p75 = quantile(x, 0.75, na.rm = TRUE))
}
summ_bin <- function(x) c(n = sum(!is.na(x)), pct = mean(as.numeric(x) == 1, na.rm = TRUE) * 100)

# Detect common columns present
num_cols <- intersect(c("age","svi_overall","acs_median_income","acs_pct_lt_hs",
                        "acs_pct_insured","acs_unemp_rate_pct","icu_los_days","mech_vent_days"),
                      names(prof_df))
bin_cols <- intersect(c("in_hosp_death","death_30d"), names(prof_df))
cat_cols <- intersect(c("sex","sex_category","race","race_ethnicity_simple"), names(prof_df))

# ----- Size + numerics -----
prof_num <- prof_df %>%
  group_by(traj_cluster_mv) %>%
  summarise(
    cluster_n = n(),
    across(all_of(num_cols), ~list(summ_num(.x))),
    .groups = "drop"
  )

# ----- Binaries -----
prof_bin <- NULL
if (length(bin_cols) > 0) {
  prof_bin <- prof_df %>%
    group_by(traj_cluster_mv) %>%
    summarise(across(all_of(bin_cols), ~list(summ_bin(.x))), .groups = "drop")
}

# ----- Categoricals (top levels) -----
top_k_table <- function(x, k = 5) {
  tb <- sort(table(x), decreasing = TRUE)
  head(round(100 * tb / sum(tb), 1), k)
}
prof_cat <- NULL
if (length(cat_cols) > 0) {
  prof_cat <- prof_df %>%
    group_by(traj_cluster_mv) %>%
    summarise(across(all_of(cat_cols), ~list(top_k_table(.x))), .groups = "drop")
}

# Combine readable profile: you can print these or explode them to a nice table later.
cluster_sizes <- prof_df %>% count(traj_cluster_mv, name = "n")

message("Cluster sizes:")
print(cluster_sizes)

message("\nNumeric summaries (mean/sd/p25/p50/p75) by cluster:")
print(prof_num)

if (!is.null(prof_bin)) {
  message("\nBinary outcome rates (%%) by cluster:")
  print(prof_bin)
}
if (!is.null(prof_cat)) {
  message("\nTop categorical levels (%%) by cluster:")
  print(prof_cat)
}

# ---------- Utilities ----------
get_medoid_id <- function(series_mat, ids_in_cluster, dist_fun, band_hint = 2L,
                          pb_title = NULL, show_progress = TRUE) {
  idx <- ids_in_cluster[ids_in_cluster %in% rownames(series_mat)]
  if (length(idx) <= 1) return(idx)
  
  k <- length(idx)
  D <- matrix(0, k, k)
  
  # total unique pairs i<j
  total_pairs <- k * (k - 1) / 2
  if (show_progress) {
    message(sprintf("%s  |  n=%d  |  pairs=%d", pb_title %||% "Medoid search", k, total_pairs))
    pb <- utils::txtProgressBar(min = 0, max = total_pairs, style = 3)
  }
  cdone <- 0L
  
  for (i in seq_len(k)) {
    xi <- as.numeric(series_mat[idx[i], ])
    for (j in i:k) {
      if (i == j) {
        D[i, j] <- 0
      } else {
        d <- dist_fun(xi, as.numeric(series_mat[idx[j], ]), band_hint)
        D[i, j] <- D[j, i] <- d
        cdone <- cdone + 1L
        if (show_progress) utils::setTxtProgressBar(pb, cdone)
      }
    }
  }
  if (show_progress) close(pb)
  
  idx[which.min(rowSums(D, na.rm = TRUE))]
}
to_long_df <- function(M, label_key, varname = "value", name = "t") {
  stopifnot(!is.null(rownames(M)))
  df <- as.data.frame(M)
  df[[ "hospitalization_id" ]] <- rownames(M)
  df_long <- tidyr::pivot_longer(df, -hospitalization_id, names_to = name, values_to = varname)
  df_long <- df_long %>%
    mutate(!!name := as.integer(factor(.data[[name]], levels = unique(.data[[name]])))) %>%
    left_join(label_key, by = "hospitalization_id")
  df_long
}

# Choose your working matrices (match what you clustered on)
A_mat <- X_clin_14d_z       # clinical signatures (z)
B_mat <- X_combo_90_z       # 90d combined exposure (z)
# (Optionally also have X_no2_90_z and X_pm25_90_z for split plots)
B_no2 <- if (exists("X_no2_90_z")) X_no2_90_z else NULL
B_pm25 <- if (exists("X_pm25_90_z")) X_pm25_90_z else NULL

label_key <- labels_fast %>%
  select(hospitalization_id, traj_cluster_mv) %>%
  mutate(traj_cluster_mv = forcats::fct_inorder(traj_cluster_mv))

# ---------- Compute medoids per cluster ----------
clusters <- levels(label_key$traj_cluster_mv)
medoids <- purrr::map_dfr(clusters, function(cl) {
  ids_cl <- label_key %>% dplyr::filter(traj_cluster_mv == cl) %>% dplyr::pull(hospitalization_id)
  med_clin <- get_medoid_id(A_mat, ids_cl,
                            function(x, y, band) dtw_dist_robust(x, y, band_hint = band),
                            band_hint = 2L,
                            pb_title = sprintf("Cluster %s — clinical", cl),
                            show_progress = TRUE)
  med_expo <- get_medoid_id(B_mat, ids_cl,
                            function(x, y, band) dtw_dist_robust(x, y, band_hint = band),
                            band_hint = 2L,
                            pb_title = sprintf("Cluster %s — exposure", cl),
                            show_progress = TRUE)
  tibble::tibble(traj_cluster_mv = cl, medoid_clin = med_clin, medoid_expo = med_expo)
})
print(medoids)

# ---------- Build long data for plotting ----------
set.seed(123)
# Thin spaghetti: up to 40 random stays / cluster
thin_ids <- label_key %>%
  group_by(traj_cluster_mv) %>%
  reframe(
    hospitalization_id = sample(hospitalization_id, size = min(40, dplyr::n()), replace = FALSE)
  ) %>%
  ungroup()

# Clinical long
clin_long <- to_long_df(A_mat, label_key, varname = "clin_z", name = "day")
clin_long_thin <- clin_long %>%
  semi_join(thin_ids, by = c("hospitalization_id","traj_cluster_mv"))

# Exposure long (combined)
expo_long <- to_long_df(B_mat, label_key, varname = "expo_z", name = "day")
expo_long_thin <- expo_long %>%
  semi_join(thin_ids, by = c("hospitalization_id","traj_cluster_mv"))

# Optional NO2/PM2.5 split
expo_no2_long <- if (!is.null(B_no2)) to_long_df(B_no2, label_key, varname = "no2_z", name = "day") else NULL
expo_pm25_long <- if (!is.null(B_pm25)) to_long_df(B_pm25, label_key, varname = "pm25_z", name = "day") else NULL

# ---------- Summaries (mean & IQR) ----------
summ_iqr <- function(x) c(mean = mean(x, na.rm = TRUE),
                          p25 = quantile(x, 0.25, na.rm = TRUE),
                          p75 = quantile(x, 0.75, na.rm = TRUE))

clin_summ <- clin_long %>%
  group_by(traj_cluster_mv, day) %>%
  summarise(across(clin_z, ~ list(summ_iqr(.x))), .groups = "drop") %>%
  tidyr::unnest_wider(clin_z)

expo_summ <- expo_long %>%
  group_by(traj_cluster_mv, day) %>%
  summarise(across(expo_z, ~ list(summ_iqr(.x))), .groups = "drop") %>%
  tidyr::unnest_wider(expo_z)

if (!is.null(expo_no2_long)) {
  expo_no2_summ <- expo_no2_long %>%
    group_by(traj_cluster_mv, day) %>%
    summarise(across(no2_z, ~ list(summ_iqr(.x))), .groups = "drop") %>%
    tidyr::unnest_wider(no2_z)
}
if (!is.null(expo_pm25_long)) {
  expo_pm25_summ <- expo_pm25_long %>%
    group_by(traj_cluster_mv, day) %>%
    summarise(across(pm25_z, ~ list(summ_iqr(.x))), .groups = "drop") %>%
    tidyr::unnest_wider(pm25_z)
}

# ---------- Pull medoid series for overlay ----------
get_row_vec <- function(M, id) as.numeric(M[id, ])
medoid_overlay_clin <- medoids %>%
  mutate(vec = purrr::map(medoid_clin, ~ get_row_vec(A_mat, .x))) %>%
  transmute(traj_cluster_mv, day = map(vec, ~ seq_along(.x)),
            value = vec) %>%
  tidyr::unnest(c(day, value)) %>%
  rename(clin_medoid = value)

medoid_overlay_expo <- medoids %>%
  mutate(vec = purrr::map(medoid_expo, ~ get_row_vec(B_mat, .x))) %>%
  transmute(traj_cluster_mv, day = map(vec, ~ seq_along(.x)),
            value = vec) %>%
  tidyr::unnest(c(day, value)) %>%
  rename(expo_medoid = value)

# ---------- Plots ----------
p_clin <- ggplot() +
  geom_ribbon(data = clin_summ,
              aes(x = day, ymin = `p25.25%`, ymax = `p75.75%`, fill = traj_cluster_mv),
              alpha = 0.18, show.legend = FALSE) +
  geom_line(data = clin_summ,
            aes(x = day, y = mean, color = traj_cluster_mv),
            linewidth = 1, show.legend = FALSE) +
  geom_line(data = medoid_overlay_clin,
            aes(x = day, y = clin_medoid),
            linewidth = 1.2, linetype = 2, color = "black") +
  geom_line(data = clin_long_thin,
            aes(x = day, y = clin_z, group = hospitalization_id),
            alpha = 0.15, linewidth = 0.3, show.legend = FALSE) +
  facet_wrap(~ traj_cluster_mv, scales = "free_y") +
  labs(title = "Clinical 14-day signature by cluster",
       x = "ICU day", y = "Within-horizon z-score") +
  theme_minimal(base_size = 12)

p_expo <- ggplot() +
  geom_ribbon(data = expo_summ,
              aes(x = day, ymin = `p25.25%`, ymax = `p75.75%`, fill = traj_cluster_mv),
              alpha = 0.18, show.legend = FALSE) +
  geom_line(data = expo_summ,
            aes(x = day, y = mean, color = traj_cluster_mv),
            linewidth = 1, show.legend = FALSE) +
  geom_line(data = medoid_overlay_expo,
            aes(x = day, y = expo_medoid),
            linewidth = 1.2, linetype = 2, color = "black") +
  geom_line(data = expo_long_thin,
            aes(x = day, y = expo_z, group = hospitalization_id),
            alpha = 0.12, linewidth = 0.25, show.legend = FALSE) +
  facet_wrap(~ traj_cluster_mv, scales = "free_y") +
  labs(title = "Exposure (90-day combined) by cluster",
       x = "Day before admission", y = "Within-window z-score") +
  theme_minimal(base_size = 12)

# Optional: split NO2 and PM2.5 panels
if (!is.null(expo_no2_long) && !is.null(expo_pm25_long)) {
  p_no2 <- ggplot() +
    geom_ribbon(data = expo_no2_summ,
                aes(x = day, ymin = `p25.25%`, ymax = `p75.75%`, fill = traj_cluster_mv),
                alpha = 0.18, show.legend = FALSE) +
    geom_line(data = expo_no2_summ,
              aes(x = day, y = mean, color = traj_cluster_mv),
              linewidth = 1, show.legend = FALSE) +
    facet_wrap(~ traj_cluster_mv, scales = "free_y") +
    labs(title = "NO\u2082 (90-day) by cluster", x = "Day before admission", y = "z") +
    theme_minimal(base_size = 12)
  
  p_pm25 <- ggplot() +
    geom_ribbon(data = expo_pm25_summ,
                aes(x = day, ymin = `p25.25%`, ymax = `p75.75%`, fill = traj_cluster_mv),
                alpha = 0.18, show.legend = FALSE) +
    geom_line(data = expo_pm25_summ,
              aes(x = day, y = mean, color = traj_cluster_mv),
              linewidth = 1, show.legend = FALSE) +
    facet_wrap(~ traj_cluster_mv, scales = "free_y") +
    labs(title = "PM\u2082.\u2085 (90-day) by cluster", x = "Day before admission", y = "z") +
    theme_minimal(base_size = 12)
}

print(p_clin)
print(p_expo)
if (exists("p_no2")) print(p_no2)
if (exists("p_pm25")) print(p_pm25)

stopifnot(all(c("hospitalization_id","traj_cluster_mv") %in% names(labels_fast)))
# ---- 1) Build a robust label_key (cluster labels per hospitalization) -------------------------

# ---------- 1) Build a robust label_key ----------
# Find a labels object
labels_obj <- NULL
for (nm in c("labels_fast","labels_14_90","labels_30_365","series_df_hr2")) {
  if (exists(nm, inherits = FALSE)) { labels_obj <- get(nm); break }
}
stopifnot(!is.null(labels_obj))

# Detect cluster column
cluster_col <- intersect(names(labels_obj),
                         c("traj_cluster_mv","cluster","traj_cluster","cluster_label","cluster_id"))
if (length(cluster_col) == 0) stop("Couldn't find a cluster column in your labels object.")
cluster_col <- cluster_col[1]

# Normalize id types and build label_key
label_key <- labels_obj %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  select(hospitalization_id, !!cluster_col) %>%
  rename(traj_cluster_mv = !!cluster_col) %>%
  distinct(hospitalization_id, .keep_all = TRUE) %>%
  mutate(traj_cluster_mv = fct_inorder(traj_cluster_mv))

# ---------- 2) Prep hourly table & choose real columns ----------
icu_hourly_idx <- icu_hourly %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  arrange(hospitalization_id, hour_ts) %>%
  group_by(hospitalization_id) %>%
  mutate(hour_idx = row_number()) %>%
  ungroup()

candidate_vars <- c(
  "any_imv","fio2_max","peep_med","rr_med","vt_med","pc_med","ps_med",
  "pip_med","pplat_med","mv_med","spo2_med","ph_med","paco2_med","pao2_med",
  "pf_ratio","sig"
)
key_vars <- intersect(candidate_vars, names(icu_hourly_idx))
if (!length(key_vars)) stop("No expected clinical columns found in icu_hourly.")

icu_hourly_idx <- icu_hourly_idx %>%
  select(hospitalization_id, hour_idx, all_of(key_vars))

# ---------- 3) Helpers ----------
robust_slope <- function(y, dt = 1) {
  x <- seq_along(y); y <- as.numeric(y)
  ok <- is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  coef(stats::lm(y[ok] ~ x[ok]))[[2]] / dt
}
q_na <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, p, names = FALSE, na.rm = TRUE))
}

# ---------- 4) Horizon summariser (labels joined inside) ----------
summarise_horizon <- function(H, icu_hourly_idx, label_key, key_vars) {
  dfH <- icu_hourly_idx %>%
    filter(hour_idx <= H) %>%
    inner_join(label_key, by = "hospitalization_id")   # ensures traj_cluster_mv exists
  
  # Per-stay medians + slopes within horizon
  per_stay <- dfH %>%
    group_by(hospitalization_id, traj_cluster_mv) %>%
    summarise(
      across(all_of(key_vars), ~median(.x, na.rm = TRUE), .names = "{.col}_med"),
      slope_spo2   = if ("spo2_med" %in% names(dfH))  robust_slope(spo2_med)  else NA_real_,
      slope_pf     = if ("pf_ratio" %in% names(dfH))  robust_slope(pf_ratio)  else NA_real_,
      slope_fio2   = if ("fio2_max" %in% names(dfH))  robust_slope(fio2_max)  else NA_real_,
      slope_peep   = if ("peep_med" %in% names(dfH))  robust_slope(peep_med)  else NA_real_,
      slope_paco2  = if ("paco2_med"%in% names(dfH))  robust_slope(paco2_med) else NA_real_,
      .groups = "drop"
    )
  
  # Determine which columns to summarise (everything except id/cluster)
  measure_cols <- setdiff(names(per_stay), c("hospitalization_id","traj_cluster_mv"))
  
  by_cluster <- per_stay %>%
    group_by(traj_cluster_mv) %>%
    summarise(
      n = n(),
      across(all_of(measure_cols), list(
        m   = ~mean(.x, na.rm = TRUE),
        p50 = ~median(.x, na.rm = TRUE),
        p25 = ~q_na(.x, .25),
        p75 = ~q_na(.x, .75)
      )),
      .groups = "drop"
    ) %>%
    mutate(horizon_hr = H)
  
  list(per_stay = per_stay, by_cluster = by_cluster)
}

# ---------- 5) Run 24/48/72 ----------
hz_list <- lapply(c(24L, 48L, 72L),
                  summarise_horizon,
                  icu_hourly_idx = icu_hourly_idx,
                  label_key = label_key,
                  key_vars = key_vars)

phenos_by_cluster <- dplyr::bind_rows(lapply(hz_list, `[[`, "by_cluster"))

# Optional: heatmap-ready
heatmap_tbl <- phenos_by_cluster %>%
  select(traj_cluster_mv, horizon_hr, n,
         tidyselect::matches("_med_(m|p50|p25|p75)$"),
         tidyselect::starts_with("slope_") & tidyselect::matches("_(m|p50|p25|p75)$")) %>%
  pivot_longer(-c(traj_cluster_mv, horizon_hr, n), names_to = "metric", values_to = "value") %>%
  arrange(horizon_hr, traj_cluster_mv)


heatmap_tbl2 <- heatmap_tbl %>%
  rename(value = any_of(c("mean_value","value"))) %>%
  mutate(
    traj_cluster_mv = as.factor(traj_cluster_mv),
    horizon_hr = as.integer(horizon_hr)
  )

# ---- Pretty metric names (edit to taste) ----
pretty_metric <- function(x) {
  x <- str_remove(x, "_med$")
  x <- str_replace_all(x, c(
    "fio2_max"   = "FiO\u2082 (%)",
    "peep_med"   = "PEEP (cmH\u2082O)",
    "rr_med"     = "Resp Rate (bpm)",
    "vt_med"     = "Vt (mL)",
    "pc_med"     = "Pcontrol (cmH\u2082O)",
    "ps_med"     = "PSupport (cmH\u2082O)",
    "pip_med"    = "PIP (cmH\u2082O)",
    "pplat_med"  = "Plateau (cmH\u2082O)",
    "mv_med"     = "Minute Vent (L/min)",
    "spo2_med"   = "SpO\u2082 (%)",
    "ph_med"     = "pH",
    "paco2_med"  = "PaCO\u2082 (mmHg)",
    "pao2_med"   = "PaO\u2082 (mmHg)",
    "pf_ratio"   = "P/F ratio",
    "sig"        = "Composite sig",
    "^slope_spo2$"  = "\u0394 SpO\u2082 /h",
    "^slope_pf$"    = "\u0394 P/F /h",
    "^slope_fio2$"  = "\u0394 FiO\u2082 /h",
    "^slope_peep$"  = "\u0394 PEEP /h",
    "^slope_paco2$" = "\u0394 PaCO\u2082 /h"
  ))
  x
}

# ---- Z-score within metric × horizon (to compare clusters fairly) ----
heatmap_scaled <- heatmap_tbl2 %>%
  group_by(metric, horizon_hr) %>%
  mutate(
    z = {
      m <- mean(value, na.rm = TRUE)
      s <- sd(value, na.rm = TRUE)
      if (!is.finite(s) || s == 0) 0 else (value - m) / s
    }
  ) %>%
  ungroup() %>%
  mutate(metric_label = pretty_metric(metric))

# ---- Order metrics by variability (within horizon) and clusters by size (optional) ----
metric_order <- heatmap_scaled %>%
  group_by(metric_label) %>%
  summarise(varz = var(z, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(varz)) %>% pull(metric_label)

cluster_order <- heatmap_scaled %>%
  count(traj_cluster_mv, name = "n") %>%
  arrange(desc(n)) %>% pull(traj_cluster_mv)

heatmap_scaled <- heatmap_scaled %>%
  mutate(
    metric_label = factor(metric_label, levels = metric_order),
    traj_cluster_mv = factor(traj_cluster_mv, levels = cluster_order)
  )

# ---- Symmetric color scale based on robust range ----
lim <- quantile(abs(heatmap_scaled$z), 0.98, na.rm = TRUE) %>% as.numeric()
lim <- if (!is.finite(lim) || lim == 0) 1 else lim

p_heat <- ggplot(heatmap_scaled,
                 aes(x = traj_cluster_mv, y = metric_label, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ horizon_hr, nrow = 1, labeller = label_both) +
  scale_fill_gradient2(
    name = "z (within metric × horizon)",
    low = "#3B4CC0", mid = "white", high = "#B40426",
    midpoint = 0, limits = c(-lim, lim), oob = scales::squish
  ) +
  labs(x = "Cluster", y = "Clinical metric",
       title = "Cluster phenotypes by clinical metrics",
       subtitle = "Values z-scored within each metric and horizon") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

print(p_heat)

is_slope <- grepl("^slope_", heatmap_tbl2$metric)

heat_levels <- heatmap_scaled %>% filter(!grepl("^slope_", metric))
heat_slopes <- heatmap_scaled %>% filter(grepl("^slope_", metric))

# (Reuse lims separately for each figure)
limL <- quantile(abs(heat_levels$z), 0.98, na.rm = TRUE) %>% as.numeric(); if (!is.finite(limL) || limL==0) limL <- 1
limS <- quantile(abs(heat_slopes$z), 0.98, na.rm = TRUE) %>% as.numeric(); if (!is.finite(limS) || limS==0) limS <- 1

p_levels <- ggplot(heat_levels, aes(traj_cluster_mv, metric_label, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ horizon_hr, nrow = 1, labeller = label_both) +
  scale_fill_gradient2("z (levels)", low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-limL, limL), oob = scales::squish) +
  labs(x = "Cluster", y = "Metric (level)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

p_slopes <- ggplot(heat_slopes, aes(traj_cluster_mv, metric_label, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ horizon_hr, nrow = 1, labeller = label_both) +
  scale_fill_gradient2("z (slopes)", low = "#3B4CC0", mid = "white", high = "#B40426",
                       midpoint = 0, limits = c(-limS, limS), oob = scales::squish) +
  labs(x = "Cluster", y = "Metric (slope)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

 print(p_levels); print(p_slopes)

 mat_feat <- heatmap_scaled %>%
   select(metric_label, traj_cluster_mv, horizon_hr, z) %>%
   tidyr::unite(col = "col", traj_cluster_mv, horizon_hr, sep = "_") %>%
   tidyr::pivot_wider(names_from = col, values_from = z) %>%
   tibble::column_to_rownames("metric_label") %>% as.matrix()
 
 # hclust on features
 row_ord <- hclust(dist(mat_feat[apply(mat_feat,1,function(v) any(is.finite(v))), , drop=FALSE]))$order
 metric_ord2 <- rownames(mat_feat)[row_ord]
 
 heatmap_clustered <- heatmap_scaled %>%
   mutate(metric_label = factor(metric_label, levels = metric_ord2),
          traj_cluster_mv = factor(traj_cluster_mv))  # keep current cluster order or do similar for columns
 
 p_heat_clustered <- ggplot(heatmap_clustered,
                            aes(x = traj_cluster_mv, y = metric_label, fill = z)) +
   geom_tile(color = "white", linewidth = 0.2) +
   facet_wrap(~ horizon_hr, nrow = 1, labeller = label_both) +
   scale_fill_gradient2(
     name = "z",
     low = "#3B4CC0", mid = "white", high = "#B40426",
     midpoint = 0, limits = c(-lim, lim), oob = scales::squish
   ) +
   labs(x = "Cluster", y = "Clinical metric",
        title = "Clustered heatmap of clinical phenotypes") +
   theme_minimal(base_size = 12) +
   theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))
 
 print(p_heat_clustered)

 # ---- High-res save helpers ----------------------------------------------------
 dir.create("figs", showWarnings = FALSE)
 
 save_plot <- function(p, filename, width = 12, height = 8, dpi = 600) {
   stopifnot(inherits(p, "ggplot"))
   # High-quality PNG (raster)
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".png")),
     plot = p, width = width, height = height, dpi = dpi, units = "in",
     device = ragg::agg_png, bg = "white", limitsize = FALSE
   )
   # Vector copies (resolution-independent for journals)
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".pdf")),
     plot = p, width = width, height = height, units = "in",
     device = cairo_pdf, bg = "white", limitsize = FALSE
   )
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".svg")),
     plot = p, width = width, height = height, units = "in",
     device = svglite::svglite, bg = "white", limitsize = FALSE
   )
   message(sprintf("Saved: figs/%s.[png|pdf|svg]", filename))
 }
 
 # ---- Suggested sizes (edit to taste) ------------------------------------------
 # Facet exposure figure (your 6-panel 90-day plot)
 if (exists("p_expo"))   save_plot(p_expo,   "exposure_90d_by_cluster", width = 12, height = 8, dpi = 600)
 
 # Clinical 14-day signature plot
 if (exists("p_clin"))   save_plot(p_clin,   "clinical_14d_signature_by_cluster", width = 12, height = 8, dpi = 600)
 
 # Early clinical signatures (24/48/72h)
 if (exists("p_sig24"))  save_plot(p_sig24,  "clinical_signature_24h", width = 12, height = 7, dpi = 600)
 if (exists("p_sig48"))  save_plot(p_sig48,  "clinical_signature_48h", width = 12, height = 7, dpi = 600)
 if (exists("p_sig72"))  save_plot(p_sig72,  "clinical_signature_72h", width = 12, height = 7, dpi = 600)
 
 # Heatmaps (phenotypes)
 if (exists("p_heat"))         save_plot(p_heat,         "phenotype_heatmap_all", width = 14, height = 7.5, dpi = 600)
 if (exists("p_levels"))       save_plot(p_levels,       "phenotype_heatmap_levels", width = 14, height = 7.5, dpi = 600)
 if (exists("p_slopes"))       save_plot(p_slopes,       "phenotype_heatmap_slopes", width = 14, height = 7.5, dpi = 600)
 if (exists("p_heat_clustered")) save_plot(p_heat_clustered, "phenotype_heatmap_clustered", width = 14, height = 7.5, dpi = 600)
 
 
 # ---- Checkpoint helpers -------------------------------------------------------
 cache_root <- file.path(getwd(), "cache")
 dir.create(cache_root, showWarnings = FALSE)
 
 stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
 cache_dir <- file.path(cache_root, paste0("referii_", stamp))
 dir.create(cache_dir, showWarnings = FALSE)
 
 save_if_exists <- function(obj_name, dir = cache_dir) {
   if (exists(obj_name, inherits = FALSE)) {
     f <- file.path(dir, paste0(obj_name, ".rds"))
     saveRDS(get(obj_name, inherits = FALSE), f)
     message("Saved: ", basename(f))
     return(invisible(f))
   } else {
     message("Skip (not found): ", obj_name)
     return(invisible(NULL))
   }
 }
 
 obj_size_MB <- function(x) as.numeric(format(object.size(x), units = "Mb"))
 
 # ---- Core clustering / DTW artifacts -----------------------------------------
 to_save <- c(
   # Medoid bank + labels / assignments
   "bank",               # list with medoid_ids, ids_sub, etc. (from fast pipeline)
   "labels_fast",        # hospitalization_id, traj_cluster_mv, nearest_medoid
   "medoids",            # per-cluster medoid ids used for plotting
   
   # Time-series matrices (expensive to rebuild)
   "X_clin_14d_z",       # clinical 14-day matrix (z)
   "X_combo_90_z",       # exposure 90-day combined (z)
   "X_no2_90_z", "X_pm25_90_z",
   "X_combo_365_z",      # if you built it
   "B_mat_ds",           # downsampled exposure matrix used for fast DTW
   
   # Hourly & daily derived tables
   "icu_hourly",         # hourly clinical with sig
   "icu_daily",          # daily signature (if used)
   "hour_grid",          # ICU hour grid
   
   # Membership at horizons
   "Xsig24", "Xsig48", "Xsig72",   # hourly signature matrices (24/48/72)
   "P24", "P48", "P72",            # soft membership frames
   
   # Features & phenotyping
   "feat_24", "feat_48", "feat_72", "features_all",
   "phenos_by_cluster", "heatmap_tbl",
   
   # Risk modeling objects
   "dat24", "dat48", "dat72",
   "mods"                  # list of glmnet models & xcols
 )
 
 invisible(lapply(to_save, save_if_exists))
 
 # ---- Plots (optional; you can re-create from saved data, but cheap to keep) ---
 plot_objs <- c("p_clin","p_expo","p_no2","p_pm25","p_sig24","p_sig48","p_sig72",
                "p_heat","p_levels","p_slopes","p_heat_clustered",
                "p_risk24","p_risk48","p_risk72","p_cal_48")
 invisible(lapply(plot_objs, save_if_exists))
 
 # Collect sizes for what we just saved
 saved_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
 manifest <- lapply(saved_files, function(fp) {
   nm <- sub("\\.rds$", "", basename(fp))
   obj <- try(readRDS(fp), silent = TRUE)
   sz  <- if (inherits(obj, "try-error")) NA_real_ else obj_size_MB(obj)
   list(object = nm, file = basename(fp), size_MB = sz)
 })
 manifest <- do.call(rbind, lapply(manifest, as.data.frame))
 
 # Key parameters / versions you might want later
 meta <- list(
   timestamp = stamp,
   seed = if (exists(".Random.seed", inherits = FALSE)) TRUE else FALSE,
   dtw_version = as.character(utils::packageVersion("dtw")),
   cluster_version = as.character(utils::packageVersion("cluster")),
   dplyr_version = as.character(utils::packageVersion("dplyr")),
   fast_pipeline = list(
     subsample = if (exists("bank")) length(bank$ids_sub) else NA_integer_,
     k = if (exists("bank")) length(bank$medoid_ids) else NA_integer_,
     weights = c(wA = 0.7, wB = 0.3),
     bands = c(bandA = 2L, bandB = 2L),
     paa_bins_90d = if (exists("B_mat_ds")) ncol(B_mat_ds) else NA_integer_
   )
 )
 
 # Write manifest + meta
 saveRDS(manifest, file.path(cache_dir, "_manifest.rds"))
 saveRDS(meta,     file.path(cache_dir, "_meta.rds"))
 
 # Also a quick human-readable text manifest
 writeLines(
   c(
     paste0("REFER-II cache @ ", stamp),
     paste0("Path: ", cache_dir),
     "",
     "Saved objects:",
     paste(sprintf("- %s (%.2f MB)", manifest$object, manifest$size_MB), collapse = "\n"),
     "",
     "Params:",
     paste(capture.output(str(meta)), collapse = "\n")
   ),
   con = file.path(cache_dir, "_manifest.txt")
 )
 
 message("Checkpoint complete in: ", cache_dir)
 
 # ---- Restore everything in a cache folder ------------------------------------
 # restore_cache <- function(dir) {
 #   stopifnot(dir.exists(dir))
 #   rds <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
 #   # Skip the meta/manifest if you don't want them as objects
 #   rds <- rds[!basename(rds) %in% c("_manifest.rds","_meta.rds")]
 #   for (fp in rds) {
 #     nm <- sub("\\.rds$", "", basename(fp))
 #     obj <- readRDS(fp)
 #     assign(nm, obj, envir = .GlobalEnv)
 #     message("Loaded: ", basename(fp))
 #   }
 #   invisible(TRUE)
 # }
 # 
 # # Example: pick the latest cache automatically
 # latest_cache <- tail(sort(list.dirs(cache_root, recursive = FALSE, full.names = TRUE)), 1)
 # restore_cache(latest_cache)
 # 
 # # Inspect what was saved
 # readRDS(file.path(latest_cache, "_manifest.rds"))
 # readRDS(file.path(latest_cache, "_meta.rds"))
 
 
 
 
 # --- Label key (cluster labels per stay) ---
 stopifnot(all(c("hospitalization_id","traj_cluster_mv") %in% names(labels_fast)))
 label_key <- labels_fast %>%
   distinct(hospitalization_id, .keep_all = TRUE) %>%
   mutate(hospitalization_id = as.character(hospitalization_id),
          traj_cluster_mv = fct_inorder(traj_cluster_mv))
 
 # --- Clinical variables to summarize ---
 clin_vars <- intersect(c(
   "any_imv","fio2_max","peep_med","rr_med","vt_med","pc_med","ps_med",
   "pip_med","pplat_med","mv_med","spo2_med","ph_med","paco2_med","pao2_med",
   "pf_ratio", "sig"
 ), names(icu_hourly))
 
 # --- Helper: robust slope in hours 1..H ---
 robust_slope <- function(y) {
   x <- seq_along(y); y <- as.numeric(y); ok <- is.finite(y)
   if (sum(ok) < 3) return(NA_real_)
   coef(stats::lm(y[ok] ~ x[ok]))[[2]]
 }
 
 # --- NA-safe stats (return NA if no finite values) ---------------------------------------------
 safe_min   <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else min(x) }
 safe_max   <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }
 safe_median<- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else stats::median(x) }
 
 # (keep your robust_slope as-is; it already returns NA when <3 finite points)
 
 # --- Rebuild the feature constructor using the safe stats --------------------------------------
 build_features_H <- function(H) {
   icu_hourly %>%
     arrange(hospitalization_id, hour_ts) %>%
     group_by(hospitalization_id) %>%
     mutate(hour_idx = dplyr::row_number()) %>%
     ungroup() %>%
     filter(hour_idx <= H) %>%
     group_by(hospitalization_id) %>%
     summarise(
       across(all_of(clin_vars), list(
         med = ~safe_median(.x),
         min = ~safe_min(.x),
         max = ~safe_max(.x)
       ), .names = "{.col}_{.fn}"),
       slope_spo2   = if ("spo2_med" %in% names(cur_data())) robust_slope(spo2_med) else NA_real_,
       slope_pf     = if ("pf_ratio" %in% names(cur_data()))  robust_slope(pf_ratio) else NA_real_,
       slope_fio2   = if ("fio2_max" %in% names(cur_data()))  robust_slope(fio2_max) else NA_real_,
       slope_peep   = if ("peep_med" %in% names(cur_data()))  robust_slope(peep_med) else NA_real_,
       slope_paco2  = if ("paco2_med"%in% names(cur_data()))  robust_slope(paco2_med) else NA_real_,
       .groups = "drop"
     ) %>%
     mutate(horizon_hr = H) %>%
     left_join(label_key, by = "hospitalization_id")
 }
 
 # Recompute features (warnings should disappear)
 feat_24 <- build_features_H(24L)
 feat_48 <- build_features_H(48L)
 feat_72 <- build_features_H(72L)
 features_all <- dplyr::bind_rows(feat_24, feat_48, feat_72)

 # ------------------------------------ 1) Table Access & Helpers ---------------------------------
 get_tbl <- function(nm) {
   ct <- get0("clif_tables", inherits = TRUE)
   if (is.null(ct)) stop("Couldn't find 'clif_tables' in your environment.")
   key <- if (nm %in% names(ct)) nm else {
     ci <- names(ct)[tolower(names(ct)) == tolower(nm)]
     if (length(ci) == 1) ci else nm
   }
   if (!key %in% names(ct)) stop(sprintf("Table '%s' not in clif_tables. Available: %s", nm, paste(names(ct), collapse = ", ")))
   janitor::clean_names(ct[[key]])
 }
 
 pick_col <- function(df, candidates, required = TRUE) {
   for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
   if (required) stop(sprintf("None of these columns found: %s", paste(candidates, collapse = ", ")))
   rep(NA, nrow(df))
 }
 
 coalesce_any <- function(data, candidates) {
   cols <- dplyr::select(data, dplyr::any_of(candidates))
   if (ncol(cols) == 0) return(rep(NA_character_, nrow(data)))
   dplyr::coalesce(!!!cols)
 }
 
 safe_ts <- function(x, tz = time_zone) {
   if (inherits(x, "POSIXt")) return(x)
   if (is.numeric(x)) {
     x2 <- ifelse(x > 1e12, x/1000, x)
     return(lubridate::as_datetime(x2, tz = tz))
   }
   suppressWarnings(lubridate::parse_date_time(
     x,
     orders = c("ymd_HMS","ymd_HM","ymd","ymdTz","ymdT","mdy_HMS","mdy_HM","mdy","dmy_HMS","dmy_HM","dmy","HMS"),
     tz = tz, quiet = TRUE
   ))
 }
 
 add_index_fields <- function(df) {
   if (!all(c("admission_dttm","discharge_dttm") %in% names(df))) {
     df <- df |>
       dplyr::left_join(
         hospitalization |>
           dplyr::select(patient_id, hospitalization_id,
                         admission_dttm, discharge_dttm,
                         admitting_service, discharge_service, zip_code),
         by = c("patient_id","hospitalization_id")
       )
   }
   df |>
     dplyr::mutate(
       index_admit     = safe_ts(admission_dttm),
       index_discharge = safe_ts(discharge_dttm),
       index_year      = lubridate::year(index_admit),
       index_date      = as.Date(index_admit)
     )
 }
 
 patient             <- get_tbl("clif_patient")
 hospitalization     <- get_tbl("clif_hospitalization")
 diagnosis           <- get_tbl("clif_hospital_diagnosis")
 support             <- get_tbl("clif_respiratory_support")
 med_admin           <- get_tbl("clif_medication_admin_continuous")
 icu_stay            <- get_tbl("clif_adt")
 vitals              <- get_tbl("clif_vitals")
 labs_df             <- get_tbl("clif_labs")   
 
 adt_tmp <- icu_stay |>
   dplyr::left_join(
     hospitalization |> dplyr::select(hospitalization_id, patient_id),
     by = join_by(hospitalization_id)
   )
 
 icu_segs <- adt_tmp |>
   mutate(
     in_raw  = in_dttm,
     out_raw = out_dttm,
     in_ts   = safe_ts(in_raw),
     out_ts  = safe_ts(out_raw),
     loccat  = tolower(location_category)
   ) |>
   filter(loccat == "icu") |>
   filter(!is.na(patient_id), !is.na(in_ts), !is.na(out_ts), out_ts > in_ts)
 

 safe_ts <- function(x) suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE))
 
 # ---- ICU LOS ---------------------------------------------------------------
 icu_los <- icu_segs %>%
   semi_join(cohort, by = "hospitalization_id") %>%
   mutate(seg_days = as.numeric(difftime(out_ts, in_ts, units = "days"))) %>%
   group_by(hospitalization_id) %>%
   summarise(icu_los_days = sum(seg_days, na.rm = TRUE), .groups = "drop")
 
 # ---- Hospital LOS (via vitals) ---------------------------------------------
 vitals_dttm <- vitals %>%
   filter(hospitalization_id %in% cohort$hospitalization_id) %>%
   mutate(vital_recorded_ts = safe_ts(recorded_dttm)) %>%
   filter(!is.na(vital_recorded_ts)) %>%
   group_by(hospitalization_id) %>%
   summarise(
     first_vital_dttm = min(vital_recorded_ts),
     last_vital_dttm  = max(vital_recorded_ts),
     .groups = "drop"
   )
 
 hosp_los <- cohort %>%
   select(hospitalization_id) %>%
   left_join(vitals_dttm, by = "hospitalization_id") %>%
   mutate(
     hosp_los_days = as.numeric(difftime(last_vital_dttm, first_vital_dttm, units = "days")),
     hosp_los_days = ifelse(is.finite(hosp_los_days), pmax(hosp_los_days, 0), NA_real_)
   ) %>%
   left_join(
     hospitalization %>% select(hospitalization_id, discharge_category, county_code),
     by = "hospitalization_id"
   )
 
 # ---- Final outcome times (for mortality linkage) ---------------------------
 final_outcome_times <- hospitalization %>%
   select(patient_id, hospitalization_id, discharge_category, discharge_dttm) %>%
   filter(hospitalization_id %in% cohort$hospitalization_id) %>%
   mutate(
     discharge_cat_low = tolower(discharge_category),
     discharge_time    = safe_ts(discharge_dttm)
   ) %>%
   left_join(patient %>% select(patient_id, death_dttm), by = "patient_id") %>%
   left_join(vitals_dttm, by = "hospitalization_id") %>%
   mutate(
     death_dttm_final = case_when(
       discharge_cat_low %in% c("expired", "hospice") & is.na(death_dttm) ~ last_vital_dttm,
       TRUE ~ death_dttm
     )
   )
 
 # ---- Mortality -------------------------------------------------------------
 mortality_instay <- cohort %>%
   left_join(final_outcome_times %>% select(hospitalization_id, death_dttm_final), by = "hospitalization_id") %>%
   mutate(
     death_ts      = safe_ts(death_dttm_final),
     in_hosp_death = as.integer(!is.na(death_ts) & death_ts >= first_icu_in & death_ts <= last_icu_out),
     death_30d     = as.integer(!is.na(death_ts) & death_ts <= (first_icu_in + days(30)))
   ) %>%
   select(hospitalization_id, in_hosp_death, death_30d)
 
 # ---- Ventilation flags & durations -----------------------------------------
 vent_flag <- support %>%
   mutate(dev_low = tolower(device_category)) %>%
   filter(str_detect(dev_low, "imv")) %>%
   semi_join(cohort, by = "hospitalization_id") %>%
   distinct(hospitalization_id) %>%
   mutate(vent_proc_flag = 1L)
 
 support_tmp <- support %>%
   left_join(hospitalization %>% select(hospitalization_id, patient_id), by = "hospitalization_id") %>%
   mutate(rec_time = safe_ts(recorded_dttm),
          dev_low  = tolower(device_category)) %>%
   filter(!is.na(rec_time)) %>%
   semi_join(cohort, by = "hospitalization_id")
 
 support_class <- support_tmp %>%
   mutate(
     is_niv         = str_detect(dev_low, "nippv|cpap|high flow nc"),
     has_vent_token = str_detect(dev_low, "imv"),
     is_invasive_vent = has_vent_token & !is_niv
   )
 
 gap_hours <- 6
 vent_durations <- support_class %>%
   arrange(hospitalization_id, rec_time) %>%
   group_by(hospitalization_id) %>%
   mutate(
     next_time   = lead(rec_time),
     next_invas  = lead(is_invasive_vent),
     gap_hr      = as.numeric(difftime(next_time, rec_time, units = "hours")),
     add_hours   = if_else(is_invasive_vent & next_invas & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours, gap_hr, 0),
     next_niv    = lead(is_niv),
     add_niv_hrs = if_else(is_niv & next_niv & !is.na(gap_hr) & gap_hr > 0 & gap_hr <= gap_hours, gap_hr, 0)
   ) %>%
   summarise(
     vent_hours = sum(add_hours, na.rm = TRUE),
     niv_hours  = sum(add_niv_hrs, na.rm = TRUE),
     .groups = "drop"
   ) %>%
   mutate(vent_proc_flag = as.integer(vent_hours > 0))
 
 # ---- AKI flag --------------------------------------------------------------
 aki_flag <- labs_df %>%
   mutate(name_low = tolower(lab_category)) %>%
   filter(str_detect(name_low, "creatinine")) %>%
   semi_join(cohort, by = "hospitalization_id") %>%
   group_by(hospitalization_id) %>%
   summarise(
     aki_flag = as.integer(
       (max(lab_value_numeric, na.rm = TRUE) - min(lab_value_numeric, na.rm = TRUE)) >= 0.3
     ),
     .groups = "drop"
   )
 
 # ---- Vasoactive medication flag -------------------------------------------
 vaso_flag <- med_admin %>%
   mutate(med_low = tolower(med_category)) %>%
   filter(str_detect(med_low, "norepinephrine|epinephrine|phenylephrine|vasopressin|dopamine")) %>%
   semi_join(cohort, by = "hospitalization_id") %>%
   distinct(hospitalization_id) %>%
   mutate(vaso_flag = 1L)
 
 # ---- Combine all outcomes --------------------------------------------------
 cohort_outcomes <- cohort %>%
   select(patient_id, hospitalization_id, first_icu_in, last_icu_out, admit_dt,
          census_tract, county_code) %>%
   left_join(icu_los,  by = "hospitalization_id") %>%
   left_join(hosp_los, by = "hospitalization_id") %>%
   left_join(mortality_instay, by = "hospitalization_id") %>%
   left_join(vaso_flag, by = "hospitalization_id") %>%
   left_join(vent_flag, by = "hospitalization_id") %>%
   left_join(vent_durations, by = "hospitalization_id") %>%
   left_join(aki_flag,  by = "hospitalization_id") %>%
   mutate(
     across(c(aki_flag, in_hosp_death, death_30d), ~ replace_na(.x, 0L)),
     vent_hours = coalesce(vent_hours, 0),
     niv_hours  = coalesce(niv_hours, 0)
   )
 
 # Outcomes from cohort (adapt names as needed)
 cohort_out <- cohort_outcomes %>%
   transmute(
     hospitalization_id = as.character(hospitalization_id),
     death_inhosp = as.integer(in_hosp_death %||% NA),
     death_30d    = as.integer(death_30d %||% NA),
     icu_los_days = as.numeric(hosp_los_days %||% NA),
     mech_vent_hours = as.numeric(vent_hours %||% NA),
     icu_los_gt5 = ifelse(is.finite(icu_los_days), as.integer(icu_los_days >= 5), NA_integer_),
     vent_gt3d   = ifelse(is.finite(mech_vent_hours), as.integer(mech_vent_hours >= 72), NA_integer_)
   )
 
 dat <- features_all %>%
   mutate(hospitalization_id = as.character(hospitalization_id),
          traj_cluster_mv = forcats::fct_inorder(traj_cluster_mv)) %>%
   left_join(cohort_out, by = "hospitalization_id")
 
 # Keep a compact set of robust features (feel free to expand)
 feat_keep <- names(dat) %>% grep(
   pattern = paste(c(
     # medians (levels)
     "sig_med","pf_ratio_med","fio2_max_med","peep_med","spo2_med","paco2_med",
     # simple dynamics
     "^slope_"
   ), collapse="|"), value = TRUE
 )
 
 # Quick missingness snapshot
 miss_tbl <- sapply(dat[feat_keep], function(x) mean(!is.finite(x))) %>% sort(decreasing = TRUE)
 message("Feature finite-rates (top missing first):")
 print(round(1 - miss_tbl, 3)[1:10])
 
 # Cluster sizes per horizon
 cluster_sizes <- dat %>%
   count(horizon_hr, traj_cluster_mv, name = "n") %>%
   arrange(horizon_hr, desc(n))
 print(cluster_sizes)
 
 q_na <- function(x, p) {
   x <- x[is.finite(x)]
   if (!length(x)) return(NA_real_)
   as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
 }
 
 # Identify columns by role
 cols_med <- grep("_med$", names(dat), value = TRUE)          # per-stay medians already
 cols_min <- grep("_min$", names(dat), value = TRUE)
 cols_max <- grep("_max$", names(dat), value = TRUE)
 cols_slp <- grep("^slope_", names(dat), value = TRUE)
 
 # Optionally, restrict to a curated subset
 feat_keep <- unique(c(cols_med, cols_min, cols_max, cols_slp))
 
 # Summarize across stays: mean + IQR (no p50 for *_med to avoid "median of medians")
 pheno_tbl <- dat %>%
   group_by(horizon_hr, traj_cluster_mv) %>%
   summarise(
     n = dplyr::n(),
     across(all_of(feat_keep),
            list(
              mean = ~mean(.x, na.rm = TRUE),
              p25  = ~q_na(.x, 0.25),
              p75  = ~q_na(.x, 0.75)
            ),
            .names = "{.col}_{.fn}"),
     .groups = "drop"
   )
 
 # If you *do* want medians across stays for certain features (e.g., slopes), add this:
 # pheno_tbl <- pheno_tbl %>%
 #   left_join(
 #     dat %>%
 #       group_by(horizon_hr, traj_cluster_mv) %>%
 #       summarise(across(all_of(cols_slp), ~median(.x, na.rm = TRUE), .names = "{.col}_p50"),
 #                 .groups = "drop"),
 #     by = c("horizon_hr","traj_cluster_mv")
 #   )
 
 # Example: a compact table with key anchors
 pheno_anchor <- pheno_tbl %>%
   select(horizon_hr, traj_cluster_mv, n,
          matches("^pf_ratio_med_(mean|p25|p75)$"),
          matches("^fio2_max_med_(mean|p25|p75)$"),
          matches("^peep_med_med_(mean|p25|p75)$"),
          matches("^spo2_med_med_(mean|p25|p75)$"),
          matches("^paco2_med_med_(mean|p25|p75)$"),
          matches("^slope_(pf|fio2|peep|paco2)_(mean|p25|p75)$")) %>%
   arrange(horizon_hr, traj_cluster_mv)
 
 print(pheno_anchor)

 # Observed risks by cluster (per horizon)
 risk_by_cluster <- function(outcome, H) {
   dat %>%
     filter(horizon_hr == H, is.finite(.data[[outcome]])) %>%
     group_by(traj_cluster_mv) %>%
     summarise(n = n(),
               events = sum(.data[[outcome]] == 1, na.rm = TRUE),
               risk_pct = 100 * events / n, .groups = "drop") %>%
     arrange(desc(risk_pct))
 }
 
 print(risk_by_cluster("death_inhosp", 24L))
 print(risk_by_cluster("death_inhosp", 48L))
 print(risk_by_cluster("death_inhosp", 72L))
 
 # Barplot helper
 plot_risk_bars <- function(outcome, H) {
   df <- risk_by_cluster(outcome, H)
   ggplot(df, aes(traj_cluster_mv, risk_pct)) +
     geom_col() +
     geom_text(aes(label = sprintf("%.1f%%", risk_pct)), vjust = -0.25, size = 3) +
     labs(title = paste0("Observed ", outcome, " by cluster (", H, "h)"),
          x = "Cluster", y = "Risk (%)") +
     theme_minimal(base_size = 12)
 }
 p_risk_48_death <- plot_risk_bars("death_inhosp", 48L)
 print(p_risk_48_death)
 

 # Minimal, interpretable model per horizon (logistic GLM)
 fit_glm_bin <- function(H, outcome = "death_inhosp") {
   df <- dat %>% filter(horizon_hr == H, is.finite(.data[[outcome]]))
   # Impute missing predictors by horizon-wise median for simplicity
   impute_med <- function(x) { ifelse(is.finite(x), x, median(x, na.rm = TRUE)) }
   xcols <- intersect(c("traj_cluster_mv",
                        "sig_med","pf_ratio_med","fio2_max_med","peep_med_med","spo2_med_med","paco2_med_med",
                        "slope_pf","slope_fio2","slope_peep","slope_paco2"),
                      names(df))
   df_small <- df %>%
     select(all_of(c("hospitalization_id", outcome, xcols))) %>%
     mutate(across(where(is.numeric), impute_med),
            traj_cluster_mv = fct_drop(traj_cluster_mv))
   # Fit
   form <- as.formula(paste(outcome, "~ traj_cluster_mv + pf_ratio_med + fio2_max_med + peep_med_med + spo2_med_med + paco2_med_med + slope_pf + slope_fio2 + slope_peep"))
   fit <- glm(form, family = binomial(), data = df_small)
   list(fit = fit, data = df_small)
 }
 
 m24 <- fit_glm_bin(24L, "death_inhosp")
 m48 <- fit_glm_bin(48L, "death_inhosp")
 m72 <- fit_glm_bin(72L, "death_inhosp")
 
 # Quick performance (resubstitution AUROC; replace with CV if you like)
 auroc_glm <- function(modobj, outcome) {
   phat <- predict(modobj$fit, type = "response")
   y <- modobj$data[[outcome]]
   as.numeric(pROC::auc(y, phat))
 }
 # install.packages("pROC") if needed
 library(pROC)
 message(sprintf("AUROC 24h: %.3f  48h: %.3f  72h: %.3f",
                 auroc_glm(m24,"death_inhosp"),
                 auroc_glm(m48,"death_inhosp"),
                 auroc_glm(m72,"death_inhosp")))
 
 or_table <- function(modobj) {
   est <- coef(summary(modobj$fit))
   tibble::tibble(
     term = rownames(est),
     OR = exp(est[, "Estimate"]),
     LCL = exp(est[, "Estimate"] - 1.96*est[, "Std. Error"]),
     UCL = exp(est[, "Estimate"] + 1.96*est[, "Std. Error"]),
     p = est[, "Pr(>|z|)"]
   )
 }
 print(or_table(m48))
 
 predict_patient_card <- function(hosp_id, H = 48L, outcome = "death_inhosp") {
   modobj <- list(`24`=m24, `48`=m48, `72`=m72)[[as.character(H)]]
   stopifnot(!is.null(modobj))
   row <- modobj$data %>% filter(hospitalization_id == as.character(hosp_id))
   if (nrow(row) == 0) stop("No row for this patient/horizon in features_all.")
   
   phat <- as.numeric(predict(modobj$fit, newdata = row, type = "response"))
   cl   <- as.character(row$traj_cluster_mv[[1]])
   
   list(
     hospitalization_id = hosp_id,
     horizon_hr = H,
     phenotype_cluster = cl,
     predicted_risk = phat[1],
     drivers = names(sort(abs(coef(modobj$fit))[-1], decreasing = TRUE))[1:5]
   )
 }
 
# Example:
 card48 <- predict_patient_card("632384141", 48L, "death_inhosp")
 str(card48)
 
 dir.create("analysis_out", showWarnings = FALSE)
 readr::write_csv(pheno_anchor, "analysis_out/phenotype_profiles_anchor.csv")
 readr::write_csv(cluster_sizes, "analysis_out/cluster_sizes_by_horizon.csv")
 if (exists("p_risk_48_death")) {
   ggsave("analysis_out/risk_by_cluster_48h.png", p_risk_48_death,
          width = 8, height = 5, dpi = 600, device = ragg::agg_png)
 }
 
 
 
 # ---- Exposure summaries per stay (90d) ----
 expo_sum_90 <- expo_daily_90 %>%
   transmute(
     hospitalization_id = as.character(hospitalization_id),
     no2_mean_90  = map_dbl(expo, ~ mean(.x$no2_raw,  na.rm = TRUE)),
     pm25_mean_90 = map_dbl(expo, ~ mean(.x$pm25_raw, na.rm = TRUE)),
     no2_med_90   = map_dbl(expo, ~ median(.x$no2_raw,  na.rm = TRUE)),
     pm25_med_90  = map_dbl(expo, ~ median(.x$pm25_raw, na.rm = TRUE))
   )
 
 # ---- Choose cutpoints (quartiles shown; change to tertiles if you prefer) ----
 mk_qtile <- function(x, k = 4) {
   q <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
   q <- unique(q)
   list(q = q, k = length(q) - 1)
 }
 
 no2_cuts  <- mk_qtile(expo_sum_90$no2_mean_90,  k = 4)
 pm25_cuts <- mk_qtile(expo_sum_90$pm25_mean_90, k = 4)
 
 expo_sum_90 <- expo_sum_90 %>%
   mutate(
     no2_q  = cut(no2_mean_90,
                  breaks = no2_cuts$q,
                  include.lowest = TRUE,
                  labels = paste0("Q", seq_len(no2_cuts$k))),
     pm25_q = cut(pm25_mean_90,
                  breaks = pm25_cuts$q,
                  include.lowest = TRUE,
                  labels = paste0("Q", seq_len(pm25_cuts$k))),
     joint_hi = case_when(
       no2_q == paste0("Q", no2_cuts$k) & pm25_q == paste0("Q", pm25_cuts$k) ~ "High NO2 & High PM2.5",
       no2_q == "Q1" & pm25_q == "Q1" ~ "Low NO2 & Low PM2.5",
       TRUE ~ "Mixed"
     )
   )
 
 
 
 # ---- Build plotting panel for first 72h ----
 plot_vars <- intersect(c("sig","fio2_max","pf_ratio","peep_med","paco2_med","spo2_med"), names(icu_hourly))
 
 icu_panel_72 <- icu_hourly %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   arrange(hospitalization_id, hour_ts) %>%
   group_by(hospitalization_id) %>%
   mutate(hour_idx = row_number()) %>%
   ungroup() %>%
   filter(hour_idx <= 72) %>%
   select(hospitalization_id, hour_idx, any_of(plot_vars)) %>%
   left_join(expo_sum_90 %>% select(hospitalization_id, no2_q, pm25_q, joint_hi,
                                    no2_mean_90, pm25_mean_90),
             by = "hospitalization_id")
 
 library(ggplot2)
 
 summ_iqr_std <- function(x) {
   x <- x[is.finite(x)]
   if (!length(x)) return(c(mean = NA_real_, p25 = NA_real_, p75 = NA_real_, n = 0))
   c(
     mean = mean(x, na.rm = TRUE),
     p25  = unname(quantile(x, 0.25, na.rm = TRUE)),
     p75  = unname(quantile(x, 0.75, na.rm = TRUE)),
     n    = length(x)
   )
 }
 
 make_traj_summ <- function(df, y, grp) {
   df %>%
     dplyr::filter(is.finite(.data[[y]]), !is.na(.data[[grp]])) %>%
     dplyr::group_by(group = .data[[grp]], hour_idx) %>%
     dplyr::summarise(stats = list(summ_iqr_std(.data[[y]])), .groups = "drop") %>%
     tidyr::unnest_wider(stats)
 }
 
 traj_sig_no2  <- make_traj_summ(icu_panel_72, "sig", "no2_q")
 traj_sig_pm25 <- make_traj_summ(icu_panel_72, "sig", "pm25_q")
 traj_pf_joint <- make_traj_summ(icu_panel_72, "pf_ratio", "joint_hi")
 
 p_sig_no2 <- ggplot(traj_sig_no2, aes(hour_idx, mean, color = group)) +
   geom_ribbon(aes(ymin = p25, ymax = p75, fill = group), alpha = 0.18, color = NA) +
   geom_line(linewidth = 1.1) +
   labs(x = "ICU hour", y = "Composite severity signature (sig)",
        title = "Hour-to-hour ICU trajectories stratified by chronic NO\u2082 (90-day mean)",
        subtitle = "Lines = mean; ribbon = IQR (p25–p75)") +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank())
 
 p_sig_pm25 <- ggplot(traj_sig_pm25, aes(hour_idx, mean, color = group)) +
   geom_ribbon(aes(ymin = p25, ymax = p75, fill = group), alpha = 0.18, color = NA) +
   geom_line(linewidth = 1.1) +
   labs(x = "ICU hour", y = "Composite severity signature (sig)",
        title = "Hour-to-hour ICU trajectories stratified by chronic PM\u2082.\u2085 (90-day mean)",
        subtitle = "Lines = mean; ribbon = IQR (p25–p75)") +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank())
 
 p_pf_joint <- ggplot(traj_pf_joint, aes(hour_idx, mean, color = group)) +
   geom_ribbon(aes(ymin = p25, ymax = p75, fill = group), alpha = 0.18, color = NA) +
   geom_line(linewidth = 1.1) +
   labs(x = "ICU hour", y = "P/F ratio",
        title = "P/F trajectories by joint exposure group (90-day mean)",
        subtitle = "Low-Low vs Mixed vs High-High") +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank())
 
 print(p_sig_no2)
 print(p_sig_pm25)
 print(p_pf_joint)
 
 # ============================================================
 # Smoothed (GAM) mean trajectories over first 72 ICU hours
 # - Smooths ONLY the mean curve (not the IQR ribbon)
 # - Produces consistent, publication-ready plots
 # ============================================================
 
 suppressPackageStartupMessages({
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   library(mgcv)
   library(ragg)
 })
 
 # ---------- 0) Standard summariser + builder (mean + IQR) ----------
 summ_iqr_std <- function(x) {
   x <- x[is.finite(x)]
   if (!length(x)) return(c(mean = NA_real_, p25 = NA_real_, p75 = NA_real_, n = 0))
   c(
     mean = mean(x, na.rm = TRUE),
     p25  = unname(quantile(x, 0.25, na.rm = TRUE)),
     p75  = unname(quantile(x, 0.75, na.rm = TRUE)),
     n    = length(x)
   )
 }
 
 make_traj_summ <- function(df, y, grp, hour = "hour_idx") {
   stopifnot(all(c(y, grp, hour) %in% names(df)))
   df %>%
     filter(is.finite(.data[[y]]), !is.na(.data[[grp]]), is.finite(.data[[hour]])) %>%
     group_by(group = as.factor(.data[[grp]]), hour_idx = .data[[hour]]) %>%
     summarise(stats = list(summ_iqr_std(.data[[y]])), .groups = "drop") %>%
     tidyr::unnest_wider(stats)
 }
 
 # ---------- 1) GAM smoother for the mean trajectory ----------
 smooth_traj_gam <- function(traj_df, k = 10) {
   # traj_df must have: group, hour_idx, mean
   traj_df %>%
     group_by(group) %>%
     group_modify(~{
       d <- .x %>% arrange(hour_idx)
       # If too sparse, leave unsmoothed
       if (sum(is.finite(d$mean)) < 6) {
         d$mean_smooth <- d$mean
         return(d)
       }
       fit <- mgcv::gam(mean ~ s(hour_idx, k = k), data = d, method = "REML")
       d$mean_smooth <- as.numeric(predict(fit, newdata = d))
       d
     }) %>%
     ungroup()
 }
 
 # ---------- 2) Plot helper ----------
 plot_traj_smoothed <- function(traj_df_s, title, ylab,
                                subtitle = "Smoothed mean (GAM); ribbon = empirical IQR",
                                min_n_per_hour = 30,
                                alpha_ribbon = 0.18,
                                lw = 1.3) {
   
   d <- traj_df_s %>%
     filter(is.finite(mean_smooth), is.finite(p25), is.finite(p75)) %>%
     # optional: drop hours where a group has tiny N (end-of-stay attrition)
     filter(n >= min_n_per_hour)
   
   ggplot(d, aes(x = hour_idx, y = mean_smooth, color = group)) +
     geom_ribbon(aes(ymin = p25, ymax = p75, fill = group),
                 alpha = alpha_ribbon, color = NA) +
     geom_line(linewidth = lw) +
     labs(x = "ICU hour", y = ylab, title = title, subtitle = subtitle) +
     theme_minimal(base_size = 12) +
     theme(legend.title = element_blank())
 }
 
 # ---------- 3) Build & plot: SIG by NO2 quartile ----------
 traj_sig_no2  <- make_traj_summ(icu_panel_72, y = "sig", grp = "no2_q")
 traj_sig_no2_s <- smooth_traj_gam(traj_sig_no2, k = 10)
 
 p_sig_no2_smooth <- plot_traj_smoothed(
   traj_sig_no2_s,
   title = "Smoothed 72-hour ICU trajectories by chronic NO\u2082 exposure (90-day mean)",
   ylab  = "Composite severity signature (sig)"
 )
 
 # ---------- 4) Build & plot: SIG by PM2.5 quartile ----------
 traj_sig_pm25  <- make_traj_summ(icu_panel_72, y = "sig", grp = "pm25_q")
 traj_sig_pm25_s <- smooth_traj_gam(traj_sig_pm25, k = 10)
 
 p_sig_pm25_smooth <- plot_traj_smoothed(
   traj_sig_pm25_s,
   title = "Smoothed 72-hour ICU trajectories by chronic PM\u2082.\u2085 exposure (90-day mean)",
   ylab  = "Composite severity signature (sig)"
 )
 
 # ---------- 5) Build & plot: P/F ratio by JOINT exposure group ----------
 # (Low-Low vs Mixed vs High-High)
 traj_pf_joint  <- make_traj_summ(icu_panel_72, y = "pf_ratio", grp = "joint_hi")
 traj_pf_joint_s <- smooth_traj_gam(traj_pf_joint, k = 10)
 
 p_pf_joint_smooth <- plot_traj_smoothed(
   traj_pf_joint_s,
   title = "Smoothed P/F trajectories by joint exposure group (90-day mean)",
   ylab  = "P/F ratio"
 )
 
 # ---------- 6) OPTIONAL: Additional physiologic plots (same style) ----------
 # Add/remove vars here depending on availability in icu_panel_72
 extra_vars <- intersect(
   c("fio2_max", "peep_med", "paco2_med", "spo2_med"),
   names(icu_panel_72)
 )
 
 # Build plots for each extra variable by NO2 quartile and PM2.5 quartile
 plots_extra <- list()
 
 for (v in extra_vars) {
   # by NO2 quartile
   t_no2  <- make_traj_summ(icu_panel_72, y = v, grp = "no2_q")
   t_no2s <- smooth_traj_gam(t_no2, k = 10)
   
   plots_extra[[paste0(v, "_no2")]] <- plot_traj_smoothed(
     t_no2s,
     title = paste0("Smoothed ", v, " over 72h by chronic NO\u2082 quartile"),
     ylab  = v
   )
   
   # by PM2.5 quartile
   t_pm   <- make_traj_summ(icu_panel_72, y = v, grp = "pm25_q")
   t_pms  <- smooth_traj_gam(t_pm, k = 10)
   
   plots_extra[[paste0(v, "_pm25")]] <- plot_traj_smoothed(
     t_pms,
     title = paste0("Smoothed ", v, " over 72h by chronic PM\u2082.\u2085 quartile"),
     ylab  = v
   )
 }
 
 # ---------- 7) Print everything ----------
 print(p_sig_no2_smooth)
 print(p_sig_pm25_smooth)
 print(p_pf_joint_smooth)
 
 # Optional extras:
 if (length(plots_extra)) {
   for (nm in names(plots_extra)) {
     print(plots_extra[[nm]])
   }
 }
 
 # ---------- 8) Save all plots (high-res PNG + vector PDF) ----------
 dir.create("figs", showWarnings = FALSE)
 
 save_plot <- function(p, filename, width = 8.5, height = 5.5, dpi = 600) {
   stopifnot(inherits(p, "ggplot"))
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".png")),
     plot = p, width = width, height = height, dpi = dpi, units = "in",
     device = ragg::agg_png, bg = "white", limitsize = FALSE
   )
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".pdf")),
     plot = p, width = width, height = height, units = "in",
     device = cairo_pdf, bg = "white", limitsize = FALSE
   )
   message("Saved: figs/", filename, ".[png|pdf]")
 }
 
 save_plot(p_sig_no2_smooth,  "traj72_sig_by_no2_quartile_smoothed")
 save_plot(p_sig_pm25_smooth, "traj72_sig_by_pm25_quartile_smoothed")
 save_plot(p_pf_joint_smooth, "traj72_pf_by_joint_exposure_smoothed")
 
 if (length(plots_extra)) {
   for (nm in names(plots_extra)) {
     save_plot(plots_extra[[nm]], paste0("traj72_", nm, "_smoothed"))
   }
 }

 
 # --- Select clinically interpretable variables ---
 clin_vars <- c("pf_ratio", "fio2_max", "peep_med")
 
 medoid_ids <- medoids %>%
   select(traj_cluster_mv, medoid_clin)
 
 medoid_icu <- icu_hourly %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(medoid_ids, by = c("hospitalization_id" = "medoid_clin")) %>%
   arrange(traj_cluster_mv, hour_ts) %>%
   group_by(traj_cluster_mv) %>%
   mutate(hour_idx = row_number()) %>%
   ungroup() %>%
   filter(hour_idx <= 72) %>%
   select(traj_cluster_mv, hour_idx, all_of(clin_vars)) %>%
   pivot_longer(-c(traj_cluster_mv, hour_idx),
                names_to = "variable", values_to = "value")
 

 
 p_dtw_clin <- ggplot(medoid_icu,
                      aes(hour_idx, value, color = traj_cluster_mv)) +
   geom_line(linewidth = 1.4) +
   facet_wrap(~ variable, scales = "free_y",
              labeller = as_labeller(c(
                pf_ratio = "P/F ratio",
                fio2_max = "FiO\u2082 (%)",
                peep_med = "PEEP (cmH\u2082O)"
              ))) +
   labs(
     x = "ICU hour",
     y = NULL,
     title = "DTW-derived ICU trajectory phenotypes",
     subtitle = "Clinical medoid trajectories over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     legend.title = element_blank(),
     strip.text = element_text(face = "bold")
   )
 
 print(p_dtw_clin)
 
 clin_vars <- c("pf_ratio", "fio2_max", "peep_med")
 
 # Helper: clamp numeric to [lo, hi]
 clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
 
 icu_hourly_clean <- icu_hourly %>%
   mutate(
     # Ensure numeric (protect against character import)
     pf_ratio = suppressWarnings(as.numeric(pf_ratio)),
     fio2_max = suppressWarnings(as.numeric(fio2_max)),
     peep_med = suppressWarnings(as.numeric(peep_med)),
     
     # ---- Harmonize FiO2 units ----
     # If FiO2 looks like a fraction (<= 1.5), convert to percent
     fio2_pct = case_when(
       is.na(fio2_max) ~ NA_real_,
       fio2_max <= 1.5 ~ fio2_max * 100,
       TRUE            ~ fio2_max
     ),
     
     # ---- Clamp to clinically plausible ranges ----
     fio2_pct  = clamp(fio2_pct, 21, 100),
     peep_med  = clamp(peep_med, 0, 30),
     pf_ratio  = clamp(pf_ratio, 0, 600)
   )
 
 
 medoid_ids <- medoids %>% select(traj_cluster_mv, medoid_clin)
 
 medoid_icu <- icu_hourly_clean %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(medoid_ids, by = c("hospitalization_id" = "medoid_clin")) %>%
   group_by(traj_cluster_mv) %>%
   arrange(hour_ts, .by_group = TRUE) %>%
   mutate(hour_idx = as.integer(difftime(hour_ts, min(hour_ts), units = "hours")) + 1L) %>%
   ungroup() %>%
   filter(hour_idx >= 1, hour_idx <= 72) %>%
   select(traj_cluster_mv, hour_idx,
          pf_ratio, fio2_pct, peep_med) %>%
   pivot_longer(
     cols = c(pf_ratio, fio2_pct, peep_med),
     names_to = "variable", values_to = "value"
   )
 
 
 # Join only needed columns and filter early
 traj_wide <- icu_hourly_clean %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(
     feat_72 %>%
       mutate(hospitalization_id = as.character(hospitalization_id)) %>%
       select(hospitalization_id, traj_cluster_mv),   # keep only what you need
     by = "hospitalization_id"
   ) %>%
   select(hospitalization_id, traj_cluster_mv, hour_ts, pf_ratio, fio2_pct, peep_med) %>%
   group_by(traj_cluster_mv, hospitalization_id) %>%
   arrange(hour_ts, .by_group = TRUE) %>%
   mutate(hour_idx = as.integer(difftime(hour_ts, min(hour_ts), units = "hours")) + 1L) %>%
   ungroup() %>%
   filter(hour_idx >= 1, hour_idx <= 72) %>%
   # optional: drop rows where all three are NA
   filter(!(is.na(pf_ratio) & is.na(fio2_pct) & is.na(peep_med)))
 
 traj_sum_wide <- traj_wide %>%
   group_by(traj_cluster_mv, hour_idx) %>%
   summarise(
     pf_med  = median(pf_ratio, na.rm = TRUE),
     pf_q25  = quantile(pf_ratio, 0.25, na.rm = TRUE),
     pf_q75  = quantile(pf_ratio, 0.75, na.rm = TRUE),
     
     fio2_med = median(fio2_pct, na.rm = TRUE),
     fio2_q25 = quantile(fio2_pct, 0.25, na.rm = TRUE),
     fio2_q75 = quantile(fio2_pct, 0.75, na.rm = TRUE),
     
     peep_med2 = median(peep_med, na.rm = TRUE),
     peep_q25  = quantile(peep_med, 0.25, na.rm = TRUE),
     peep_q75  = quantile(peep_med, 0.75, na.rm = TRUE),
     .groups = "drop"
   )
 
 # Now pivot the much smaller summary table
 traj_sum <- bind_rows(
   traj_sum_wide %>% transmute(traj_cluster_mv, hour_idx, variable="pf_ratio",  med=pf_med,  q25=pf_q25,  q75=pf_q75),
   traj_sum_wide %>% transmute(traj_cluster_mv, hour_idx, variable="fio2_pct",  med=fio2_med, q25=fio2_q25, q75=fio2_q75),
   traj_sum_wide %>% transmute(traj_cluster_mv, hour_idx, variable="peep_med",  med=peep_med2,q25=peep_q25, q75=peep_q75)
 )
 
 var_labs <- c(
   pf_ratio = "P/F ratio",
   fio2_pct = "FiO\u2082 (%)",
   peep_med = "PEEP (cmH\u2082O)"
 )
 
 p_dtw_clin <- ggplot() +
   # IQR ribbon for each cluster
   geom_ribbon(
     data = traj_sum,
     aes(x = hour_idx, ymin = q25, ymax = q75, fill = traj_cluster_mv),
     alpha = 0.15,
     colour = NA
   ) +
   # Median trajectory per cluster
   geom_line(
     data = traj_sum,
     aes(x = hour_idx, y = med, color = traj_cluster_mv),
     linewidth = 1.1
   ) +
   # Overlay medoid trajectory
   geom_line(
     data = medoid_icu,
     aes(x = hour_idx, y = value, color = traj_cluster_mv),
     linewidth = 1.8
   ) +
   facet_wrap(
     ~ variable,
     scales = "free_y",
     labeller = as_labeller(var_labs)
   ) +
   labs(
     x = "ICU hour",
     y = NULL,
     title = "DTW-derived ICU trajectory phenotypes",
     subtitle = "Cluster median (line) and IQR (ribbon) with clinical medoid overlay over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     legend.title = element_blank(),
     strip.text = element_text(face = "bold")
   )
 
 p_dtw_clin
 
 ggsave(
   filename = "figures/aim1_dtw_trajectory_phenotypes.png",
   plot     = p_dtw_clin,
   width    = 11,
   height   = 4.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 
 
 # -------------------- Settings --------------------
 cluster_levels <- paste0("C", 1:6)
 
 var_labs <- c(
   pf_ratio = "P/F ratio",
   fio2_pct = "FiO\u2082 (%)",
   peep_med = "PEEP (cmH\u2082O)"
 )
 
 # -------------------- Ensure clean ordering (C1-C6) --------------------
 traj_sum <- traj_sum %>%
   mutate(traj_cluster_mv = factor(traj_cluster_mv, levels = cluster_levels))
 
 medoid_icu <- medoid_icu %>%
   mutate(
     traj_cluster_mv = factor(traj_cluster_mv, levels = cluster_levels),
     # ensure variable names match facet labels
     variable = case_when(
       variable %in% c("pf_ratio", "fio2_pct", "peep_med") ~ variable,
       variable == "fio2_max" ~ "fio2_pct",   # in case your medoid table still uses fio2_max
       TRUE ~ variable
     )
   )
 
 # -------------------- Smooth cluster median lines (rolling median) --------------------
 # k=5 ~ smooth but preserves step-changes; adjust to 3 or 7 if desired
 traj_sum <- traj_sum %>%
   group_by(traj_cluster_mv, variable) %>%
   arrange(hour_idx, .by_group = TRUE) %>%
   mutate(med_smooth = zoo::rollmedian(med, k = 5, fill = NA, align = "center")) %>%
   ungroup()
 
 # -------------------- Plot: ribbon + smoothed median + medoid overlay --------------------
 p_dtw_clin <- ggplot() +
   geom_ribbon(
     data = traj_sum,
     aes(
       x = hour_idx,
       ymin = q25,
       ymax = q75,
       fill = traj_cluster_mv,
       group = interaction(traj_cluster_mv, variable)
     ),
     alpha = 0.15,
     colour = NA
   ) +
   geom_line(
     data = traj_sum,
     aes(
       x = hour_idx,
       y = med_smooth,
       color = traj_cluster_mv,
       group = interaction(traj_cluster_mv, variable)
     ),
     linewidth = 1.2
   ) +
   geom_line(
     data = medoid_icu,
     aes(
       x = hour_idx,
       y = value,
       color = traj_cluster_mv,
       group = interaction(traj_cluster_mv, variable)
     ),
     linewidth = 1.8,
     alpha = 0.8
   ) +
   facet_wrap(
     ~ variable,
     scales = "free_y",
     labeller = as_labeller(var_labs)
   ) +
   # enforce legend order C1-C6 for both color and fill
   scale_color_discrete(limits = cluster_levels) +
   scale_fill_discrete(limits = cluster_levels) +
   labs(
     x = "ICU hour",
     y = NULL,
     title = "DTW-derived ICU trajectory phenotypes",
     subtitle = "Cluster median (line) and IQR (ribbon) with clinical medoid overlay over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     legend.title = element_blank(),
     legend.position = "right",
     strip.text = element_text(face = "bold"),
     panel.grid.minor = element_blank()
   )
 
 print(p_dtw_clin)
 
 ggsave(
   filename = "figures/aim1_dtw_trajectory_phenotypes2.png",
   plot     = p_dtw_clin,
   width    = 11,
   height   = 4.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 

 p_dtw_pf <- ggplot() +
   geom_ribbon(
     data = traj_sum %>% filter(variable == "pf_ratio"),
     aes(
       x = hour_idx,
       ymin = q25,
       ymax = q75,
       fill = traj_cluster_mv,
       group = traj_cluster_mv
     ),
     alpha = 0.15,
     colour = NA
   ) +
   geom_line(
     data = traj_sum %>% filter(variable == "pf_ratio"),
     aes(
       x = hour_idx,
       y = med_smooth,
       color = traj_cluster_mv,
       group = traj_cluster_mv
     ),
     linewidth = 1.3
   ) +
   geom_line(
     data = medoid_icu %>% filter(variable == "pf_ratio"),
     aes(
       x = hour_idx,
       y = value,
       color = traj_cluster_mv,
       group = traj_cluster_mv
     ),
     linewidth = 2.0,
     alpha = 0.8
   ) +
   scale_color_discrete(limits = cluster_levels) +
   scale_fill_discrete(limits = cluster_levels) +
   labs(
     x = "ICU hour",
     y = "P/F ratio",
     title = "DTW-derived ICU trajectory phenotypes",
     subtitle = "P/F ratio trajectories over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     legend.title = element_blank(),
     legend.position = "right",
     panel.grid.minor = element_blank()
   )
 
 
 ggsave(
   filename = "figures/aim1_dtw_pf_ratio_trajectory.png",
   plot     = p_dtw_pf,
   width    = 6.5,
   height   = 4.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 
 
 
 
 
 
 
 
 
 medoid_ids <- medoids %>%
   select(traj_cluster_mv, medoid_clin, medoid_expo) %>%
   mutate(
     medoid_clin = as.character(medoid_clin),
     medoid_expo = as.character(medoid_expo)
   )
 
 # Extract exposure medoids (90-day)
 expo_list_to_long <- function(traj_cluster, expo_obj) {
   n <- length(expo_obj$no2_raw)
   tibble::tibble(
     traj_cluster_mv = traj_cluster,
     day = seq_len(n) - n,  # -89..-1 for 90d
     no2_raw  = as.numeric(expo_obj$no2_raw),
     pm25_raw = as.numeric(expo_obj$pm25_raw)
   )
 }
 
 medoid_expo <- expo_daily_90 %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(medoid_ids, by = c("hospitalization_id" = "medoid_expo")) %>%
   select(traj_cluster_mv, expo) %>%
   mutate(expo_long = purrr::map2(traj_cluster_mv, expo, expo_list_to_long)) %>%
   select(expo_long) %>%
   tidyr::unnest(expo_long) %>%
   pivot_longer(c(no2_raw, pm25_raw),
                names_to = "pollutant", values_to = "value")
 
 
 p_dtw_expo <- ggplot(medoid_expo,
                      aes(day, value, color = traj_cluster_mv)) +
   geom_line(linewidth = 1.3) +
   facet_wrap(~ pollutant, scales = "free_y",
              labeller = as_labeller(c(
                no2_raw  = "NO\u2082 (ppb)",
                pm25_raw = "PM\u2082.\u2085 (\u03bcg/m\u00b3)"
              ))) +
   labs(x = "Days before ICU admission", y = NULL,
        title = "DTW exposure medoids by clinical phenotype") +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank())
 
 print(p_dtw_expo)
 
 ggsave(
   filename = "figures/aim1_p_dtw_expo.png",
   plot     = p_dtw_expo,
   width    = 6.5,
   height   = 4.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 
 
 # Ensure hospitalization_id types match
 prof_df <- prof_df %>%
   mutate(hospitalization_id = as.character(hospitalization_id))
 
 if (exists("cohort_out")) {
   cohort_out2 <- cohort_out %>%
     mutate(hospitalization_id = as.character(hospitalization_id))
   
   prof_df <- prof_df %>%
     left_join(cohort_out2, by = "hospitalization_id")
 } else if (exists("cohort_outcomes")) {
   # Fallback: derive the same fields from cohort_outcomes
   cohort_out2 <- cohort_outcomes %>%
     transmute(
       hospitalization_id = as.character(hospitalization_id),
       death_inhosp       = as.integer(in_hosp_death),
       death_30d          = as.integer(death_30d),
       icu_los_days       = as.numeric(icu_los_days),
       hosp_los_days      = as.numeric(hosp_los_days),
       mech_vent_hours    = as.numeric(vent_hours %||% 0),
       vent_days          = mech_vent_hours / 24
     )
   
   prof_df <- prof_df %>%
     left_join(cohort_out2, by = "hospitalization_id")
 } else {
   stop("Neither 'cohort_out' nor 'cohort_outcomes' exists yet. Run the outcomes construction block first.")
 }
 
 # Quick check
 dplyr::glimpse(prof_df)
 
 
 
 
 
 # helper: first existing name
 pick1 <- function(df, candidates) {
   hit <- intersect(candidates, names(df))
   if (length(hit) == 0) NA_character_ else hit[1]
 }
 
 death_col <- pick1(prof_df, c("in_hosp_death", "death_inhosp", "death_in_hosp", "mortality_in_hosp"))
 vent_col  <- pick1(prof_df, c("mech_vent_days", "mech_vent_hours", "vent_days", "vent_hours"))
 los_col   <- pick1(prof_df, c("icu_los_days", "hosp_los_days", "icu_los"))
 
 message("Using columns: ",
         "\n  death_col = ", death_col,
         "\n  vent_col  = ", vent_col,
         "\n  los_col   = ", los_col)
 
 stopifnot(!is.na(death_col))  # at least mortality should exist for the outcomes panel
 
 p_dtw_outcomes <- prof_df %>%
   mutate(
     death = as.numeric(.data[[death_col]]),
     vent  = if (!is.na(vent_col)) as.numeric(.data[[vent_col]]) else NA_real_,
     los   = if (!is.na(los_col))  as.numeric(.data[[los_col]])  else NA_real_,
     vent_days = if (!is.na(vent_col) && grepl("hour", vent_col, ignore.case = TRUE)) vent / 24 else vent
   ) %>%
   group_by(traj_cluster_mv) %>%
   summarise(
     Mortality_pct = 100 * mean(death == 1, na.rm = TRUE),
     Vent_days     = mean(vent_days, na.rm = TRUE),
     ICU_LOS_days  = mean(los, na.rm = TRUE),
     .groups = "drop"
   ) %>%
   pivot_longer(-traj_cluster_mv, names_to = "outcome", values_to = "value") %>%
   # drop outcomes that are entirely NA
   group_by(outcome) %>%
   filter(any(is.finite(value))) %>%
   ungroup() %>%
   ggplot(aes(traj_cluster_mv, value, fill = traj_cluster_mv)) +
   geom_col(show.legend = FALSE) +
   facet_wrap(~ outcome, scales = "free_y") +
   labs(
     x = "DTW phenotype",
     y = NULL,
     title = "Clinical outcomes by DTW-derived phenotype"
   ) +
   theme_minimal(base_size = 12)
 
 print(p_dtw_outcomes)

 
 # -------------------- Make outcomes panel grant-ready --------------------
 
 # 1) Enforce cluster order and ensure colors match your prior plots
 cluster_levels <- paste0("C", 1:6)
 
 # Use the SAME palette mapping you used previously (if you already have it).
 # If you do NOT already have a named palette, the safest approach is to
 # extract the discrete palette ggplot assigned in the previous plot and reuse it.
 # Here: define a named palette explicitly (edit if you already have one).
 # Example: dtw_pal <- setNames(scales::hue_pal()(6), cluster_levels)
 # Better: if you already used scale_color_manual/scale_fill_manual earlier,
 # paste that SAME named vector here.
 dtw_pal <- setNames(scales::hue_pal()(6), cluster_levels)
 
 
 # 3) Build grant-ready outcome labels (no variable names) + plot
 outcome_labels <- c(
   Mortality_pct = "In-hospital mortality (%)",
   Vent_days     = "Ventilation duration (days)",
   ICU_LOS_days  = "ICU length of stay (days)"
 )
 
 p_dtw_outcomes <- prof_df %>%
   mutate(
     traj_cluster_mv = factor(traj_cluster_mv, levels = cluster_levels),
     death = as.numeric(.data[[death_col]]),
     vent  = if (!is.na(vent_col)) as.numeric(.data[[vent_col]]) else NA_real_,
     los   = if (!is.na(los_col))  as.numeric(.data[[los_col]])  else NA_real_,
     vent_days = if (!is.na(vent_col) && grepl("hour", vent_col, ignore.case = TRUE)) vent / 24 else vent
   ) %>%
   group_by(traj_cluster_mv) %>%
   summarise(
     Mortality_pct = 100 * mean(death == 1, na.rm = TRUE),
     Vent_days     = mean(vent_days, na.rm = TRUE),
     ICU_LOS_days  = mean(los, na.rm = TRUE),
     n             = dplyr::n(),
     .groups = "drop"
   ) %>%
   pivot_longer(
     cols = c(Mortality_pct, Vent_days, ICU_LOS_days),
     names_to = "outcome",
     values_to = "value"
   ) %>%
   # drop outcomes that are entirely NA
   group_by(outcome) %>%
   filter(any(is.finite(value))) %>%
   ungroup() %>%
   mutate(
     outcome = factor(outcome, levels = names(outcome_labels),
                      labels = unname(outcome_labels))
   ) %>%
   ggplot(aes(x = traj_cluster_mv, y = value, fill = traj_cluster_mv)) +
   geom_col(width = 0.75, show.legend = FALSE) +
   facet_wrap(~ outcome, scales = "free_y") +
   scale_fill_manual(values = dtw_pal, limits = cluster_levels) +
   labs(
     x = "DTW phenotype",
     y = NULL,
     title = "Clinical outcomes by DTW-derived phenotype"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     strip.text = element_text(face = "bold"),
     panel.grid.minor = element_blank()
   )
 
 print(p_dtw_outcomes)
 
 
 ggsave(
   filename = "figures/aim1_p_dtw_outcomes.png",
   plot     = p_dtw_outcomes,
   width    = 11,
   height   = 4.5,
   units    = "in",
   dpi      = 600,
   bg       = "white"
 )
 
 
 
 
 
 
 
 
 
 
 q25 <- function(x) as.numeric(quantile(x, 0.25, na.rm = TRUE))
 q75 <- function(x) as.numeric(quantile(x, 0.75, na.rm = TRUE))
 
 p_dtw_outcomes_med <- prof_df %>%
   mutate(
     death = as.numeric(.data[[death_col]]),
     vent  = if (!is.na(vent_col)) as.numeric(.data[[vent_col]]) else NA_real_,
     los   = if (!is.na(los_col))  as.numeric(.data[[los_col]])  else NA_real_,
     vent_days = if (!is.na(vent_col) && grepl("hour", vent_col, ignore.case = TRUE)) vent / 24 else vent
   ) %>%
   group_by(traj_cluster_mv) %>%
   summarise(
     Mortality_pct = 100 * mean(death == 1, na.rm = TRUE),
     Vent_days_p50 = median(vent_days, na.rm = TRUE),
     Vent_days_p25 = q25(vent_days),
     Vent_days_p75 = q75(vent_days),
     ICU_LOS_p50   = median(los, na.rm = TRUE),
     ICU_LOS_p25   = q25(los),
     ICU_LOS_p75   = q75(los),
     .groups = "drop"
   )
 

 dir.create("figs", showWarnings = FALSE)
 
 save_plot <- function(p, fname,
                       width = 8.5, height = 5.5, dpi = 600) {
   stopifnot(inherits(p, "ggplot"))
   
   # High-res raster (journals often want this)
   ggsave(
     filename = file.path("figs", paste0(fname, ".png")),
     plot = p,
     width = width, height = height, units = "in",
     dpi = dpi,
     device = ragg::agg_png,
     bg = "white",
     limitsize = FALSE
   )
   
   # Vector (for Illustrator / Inkscape / LaTeX)
   ggsave(
     filename = file.path("figs", paste0(fname, ".pdf")),
     plot = p,
     width = width, height = height, units = "in",
     device = cairo_pdf,
     bg = "white",
     limitsize = FALSE
   )
   
   message("Saved: figs/", fname, ".[png|pdf]")
 }
 
 # ------------------------------------------------------------
 # Smoothed exposure–trajectory plots (72h)
 # ------------------------------------------------------------
 if (exists("p_sig_no2_smooth"))
   save_plot(p_sig_no2_smooth,
             "traj72_sig_by_no2_quartile_smoothed")
 
 if (exists("p_sig_pm25_smooth"))
   save_plot(p_sig_pm25_smooth,
             "traj72_sig_by_pm25_quartile_smoothed")
 
 if (exists("p_pf_joint_smooth"))
   save_plot(p_pf_joint_smooth,
             "traj72_pf_by_joint_exposure_smoothed")
 
 # ------------------------------------------------------------
 # DTW clinical medoid trajectories
 # ------------------------------------------------------------
 if (exists("p_dtw_clin"))
   save_plot(p_dtw_clin,
             "dtw_clinical_medoid_trajectories_72h",
             width = 10, height = 6)
 
 # ------------------------------------------------------------
 # DTW exposure medoid trajectories (pre-admission)
 # ------------------------------------------------------------
 if (exists("p_dtw_expo"))
   save_plot(p_dtw_expo,
             "dtw_exposure_medoid_trajectories_90d",
             width = 10, height = 6)
 
 # ------------------------------------------------------------
 # DTW outcomes panel
 # ------------------------------------------------------------
 if (exists("p_dtw_outcomes"))
   save_plot(p_dtw_outcomes,
             "dtw_outcomes_by_phenotype",
             width = 9, height = 5)
 
 # ------------------------------------------------------------
 # Optional: save any extra smoothed physiologic plots
 # (from the earlier loop, if you ran it)
 # ------------------------------------------------------------
 if (exists("plots_extra") && length(plots_extra)) {
   for (nm in names(plots_extra)) {
     save_plot(plots_extra[[nm]],
               paste0("traj72_", nm, "_smoothed"))
   }
 }
 
 message("All available plots saved to ./figs/") 
 
 
 # ============================================================
 # FORCE k = 3 DTW CLUSTERS  →  LABELS  →  MEDOIDS  →  PLOTS
 # Assumes you already built:
 #   - A_mat      (e.g., X_clin_14d_z)  rows = hospitalization_id
 #   - B_mat_ds   (downsampled exposure matrix, e.g., 90d -> 30 bins) rows = hospitalization_id
 #   - icu_hourly (hour-level clinical table with hour_ts + pf_ratio/fio2_max/peep_med/etc.)
 #   - expo_daily_90 (tibble with hospitalization_id + list-col expo containing no2_raw/pm25_raw)
 #   - cohort_out OR cohort_outcomes (outcomes per hospitalization_id)
 #   - dtw_dist_robust(), clean_vec(), build_medoids_bank(), assign_to_medoids()
 # ============================================================
 
 suppressPackageStartupMessages({
   library(dplyr)
   library(tidyr)
   library(purrr)
   library(forcats)
   library(ggplot2)
   library(mgcv)
   library(ragg)
 })
 
 stopifnot(exists("A_mat"), exists("B_mat_ds"),
           exists("build_medoids_bank"), exists("assign_to_medoids"),
           exists("icu_hourly"), exists("expo_daily_90"))
 
 # ------------------------------------------------------------
 # 1) Run DTW clustering with k = 3 (forced)
 # ------------------------------------------------------------
 set.seed(123)
 
 bank3 <- build_medoids_bank(
   A_mat, B_mat_ds,
   wA = 0.7, wB = 0.3,
   k = 3,                 # <-- force 3
   subsample = 500,
   seed = 123,
   bandA = 2L, bandB = 2L
 )
 
 labels3 <- assign_to_medoids(
   A_mat, B_mat_ds,
   bank3$medoid_ids,
   wA = 0.7, wB = 0.3,
   bandA = 2L, bandB = 2L
 )
 
 labels3 <- labels3 %>%
   mutate(
     hospitalization_id = as.character(hospitalization_id),
     traj_cluster_mv = fct_inorder(traj_cluster_mv)
   )
 
 message("Cluster sizes (k=3):")
 print(labels3 %>% count(traj_cluster_mv, name = "n"))
 
 label_key3 <- labels3 %>% select(hospitalization_id, traj_cluster_mv)
 
 # ------------------------------------------------------------
 # 2) Recompute medoids per cluster (clinical + exposure)
 #    Uses your existing get_medoid_id() helper if present;
 #    otherwise defines a minimal version here.
 # ------------------------------------------------------------
 if (!exists("get_medoid_id")) {
   get_medoid_id <- function(series_mat, ids_in_cluster, dist_fun, band_hint = 2L) {
     idx <- ids_in_cluster[ids_in_cluster %in% rownames(series_mat)]
     if (length(idx) <= 1) return(idx)
     k <- length(idx)
     D <- matrix(0, k, k)
     for (i in seq_len(k)) {
       xi <- as.numeric(series_mat[idx[i], ])
       for (j in i:k) {
         if (i == j) {
           D[i, j] <- 0
         } else {
           d <- dist_fun(xi, as.numeric(series_mat[idx[j], ]), band_hint)
           D[i, j] <- D[j, i] <- d
         }
       }
     }
     idx[which.min(rowSums(D, na.rm = TRUE))]
   }
 }
 
 clusters3 <- levels(label_key3$traj_cluster_mv)
 
 medoids3 <- purrr::map_dfr(clusters3, function(cl) {
   ids_cl <- label_key3 %>% filter(traj_cluster_mv == cl) %>% pull(hospitalization_id)
   med_clin <- get_medoid_id(
     A_mat, ids_cl,
     function(x, y, band) dtw_dist_robust(x, y, band_hint = band),
     band_hint = 2L
   )
   med_expo <- get_medoid_id(
     B_mat_ds, ids_cl,
     function(x, y, band) dtw_dist_robust(x, y, band_hint = band),
     band_hint = 2L
   )
   tibble(traj_cluster_mv = cl, medoid_clin = as.character(med_clin), medoid_expo = as.character(med_expo))
 })
 
 print(medoids3)
 
 # ------------------------------------------------------------
 # 3) DTW clinical medoid trajectories (first 72 ICU hours)
 # ------------------------------------------------------------
 clin_vars <- intersect(c("pf_ratio", "fio2_max", "peep_med"), names(icu_hourly))
 stopifnot(length(clin_vars) >= 2)  # require at least 2 for a meaningful panel
 
 medoid_clin_ids <- medoids3 %>% select(traj_cluster_mv, medoid_clin)
 
 medoid_icu <- icu_hourly %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(medoid_clin_ids, by = c("hospitalization_id" = "medoid_clin")) %>%
   arrange(traj_cluster_mv, hour_ts) %>%
   group_by(traj_cluster_mv) %>%
   mutate(hour_idx = row_number()) %>%
   ungroup() %>%
   filter(hour_idx <= 72) %>%
   select(traj_cluster_mv, hour_idx, all_of(clin_vars)) %>%
   pivot_longer(-c(traj_cluster_mv, hour_idx), names_to = "variable", values_to = "value")
 
 var_labs <- c(
   pf_ratio = "P/F ratio",
   fio2_max = "FiO\u2082 (%)",
   peep_med = "PEEP (cmH\u2082O)"
 )
 
 p_dtw_clin3 <- ggplot(medoid_icu, aes(hour_idx, value, color = traj_cluster_mv)) +
   geom_line(linewidth = 1.35) +
   facet_wrap(~ variable, scales = "free_y",
              labeller = as_labeller(var_labs[names(var_labs) %in% unique(medoid_icu$variable)])) +
   labs(
     x = "ICU hour",
     y = NULL,
     title = "DTW-derived ICU trajectory phenotypes (k=3)",
     subtitle = "Clinical medoid trajectories over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank(), strip.text = element_text(face = "bold"))
 
 print(p_dtw_clin3)
 
 # ------------------------------------------------------------
 # 4) DTW exposure medoid trajectories (90 days pre-admission)
 #    (Uses expo_daily_90 list-column 'expo' with $no2_raw and $pm25_raw)
 # ------------------------------------------------------------
 expo_list_to_long <- function(traj_cluster, expo_obj) {
   n <- length(expo_obj$no2_raw)
   tibble(
     traj_cluster_mv = traj_cluster,
     day = seq_len(n) - n,  # e.g., -89..-1
     no2_raw  = as.numeric(expo_obj$no2_raw),
     pm25_raw = as.numeric(expo_obj$pm25_raw)
   )
 }
 
 medoid_expo_ids <- medoids3 %>% select(traj_cluster_mv, medoid_expo)
 
 medoid_expo_long <- expo_daily_90 %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   inner_join(medoid_expo_ids, by = c("hospitalization_id" = "medoid_expo")) %>%
   select(traj_cluster_mv, expo) %>%
   mutate(expo_long = map2(traj_cluster_mv, expo, expo_list_to_long)) %>%
   select(expo_long) %>%
   unnest(expo_long) %>%
   pivot_longer(c(no2_raw, pm25_raw), names_to = "pollutant", values_to = "value")
 
 poll_labs <- c(
   no2_raw  = "NO\u2082 (ppb)",
   pm25_raw = "PM\u2082.\u2085 (\u03bcg/m\u00b3)"
 )
 
 p_dtw_expo3 <- ggplot(medoid_expo_long, aes(day, value, color = traj_cluster_mv)) +
   geom_line(linewidth = 1.25) +
   facet_wrap(~ pollutant, scales = "free_y",
              labeller = as_labeller(poll_labs)) +
   labs(
     x = "Days before ICU admission",
     y = NULL,
     title = "DTW exposure medoids by clinical phenotype (k=3)",
     subtitle = "90-day pre-admission exposure trajectories for medoid stays"
   ) +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank(), strip.text = element_text(face = "bold"))
 
 print(p_dtw_expo3)
 
 # ------------------------------------------------------------
 # 5) Outcomes by DTW phenotype (attach outcomes then plot)
 # ------------------------------------------------------------
 `%||%` <- function(a, b) if (!is.null(a)) a else b
 
 prof_df3 <- labels3 %>%
   select(hospitalization_id, traj_cluster_mv, nearest_medoid) %>%
   mutate(hospitalization_id = as.character(hospitalization_id))
 
 if (exists("cohort_out")) {
   cohort_out2 <- cohort_out %>%
     mutate(hospitalization_id = as.character(hospitalization_id))
   prof_df3 <- prof_df3 %>% left_join(cohort_out2, by = "hospitalization_id")
 } else if (exists("cohort_outcomes")) {
   cohort_out2 <- cohort_outcomes %>%
     transmute(
       hospitalization_id = as.character(hospitalization_id),
       death_inhosp       = as.integer(in_hosp_death),
       death_30d          = as.integer(death_30d),
       icu_los_days       = as.numeric(icu_los_days),
       hosp_los_days      = as.numeric(hosp_los_days),
       mech_vent_hours    = as.numeric(vent_hours %||% 0),
       vent_days          = mech_vent_hours / 24
     )
   prof_df3 <- prof_df3 %>% left_join(cohort_out2, by = "hospitalization_id")
 } else {
   warning("No cohort outcomes object found. Skipping outcomes plot.")
 }
 
 pick1 <- function(df, candidates) {
   hit <- intersect(candidates, names(df))
   if (length(hit) == 0) NA_character_ else hit[1]
 }
 
 if (any(c("death_inhosp","in_hosp_death","death_30d") %in% names(prof_df3))) {
   
   death_col <- pick1(prof_df3, c("death_inhosp","in_hosp_death"))
   vent_col  <- pick1(prof_df3, c("vent_days","mech_vent_days","mech_vent_hours","vent_hours"))
   los_col   <- pick1(prof_df3, c("icu_los_days","icu_los"))
   
   prof_df3b <- prof_df3 %>%
     mutate(
       death = as.numeric(.data[[death_col]]),
       vent  = if (!is.na(vent_col)) as.numeric(.data[[vent_col]]) else NA_real_,
       los   = if (!is.na(los_col))  as.numeric(.data[[los_col]])  else NA_real_,
       vent_days = case_when(
         is.na(vent_col) ~ NA_real_,
         grepl("hour", vent_col, ignore.case = TRUE) ~ vent / 24,
         TRUE ~ vent
       )
     )
   
   p_dtw_outcomes3 <- prof_df3b %>%
     group_by(traj_cluster_mv) %>%
     summarise(
       Mortality_pct = 100 * mean(death == 1, na.rm = TRUE),
       Vent_days     = mean(vent_days, na.rm = TRUE),
       ICU_LOS_days  = mean(los, na.rm = TRUE),
       n             = n(),
       .groups = "drop"
     ) %>%
     pivot_longer(cols = c(Mortality_pct, Vent_days, ICU_LOS_days),
                  names_to = "outcome", values_to = "value") %>%
     group_by(outcome) %>%
     filter(any(is.finite(value))) %>%
     ungroup() %>%
     ggplot(aes(traj_cluster_mv, value, fill = traj_cluster_mv)) +
     geom_col(show.legend = FALSE) +
     facet_wrap(~ outcome, scales = "free_y",
                labeller = as_labeller(c(
                  Mortality_pct = "In-hospital mortality (%)",
                  Vent_days     = "Ventilation duration (days)",
                  ICU_LOS_days  = "ICU length of stay (days)"
                ))) +
     labs(
       x = "DTW phenotype",
       y = NULL,
       title = "Clinical outcomes by DTW-derived phenotype (k=3)"
     ) +
     theme_minimal(base_size = 12)
   
   print(p_dtw_outcomes3)
 }
 
 # ------------------------------------------------------------
 # 6) OPTIONAL: Smooth the clinical medoid trajectories (GAM)
 #    (Useful if you want a cleaner clinician-facing panel.)
 # ------------------------------------------------------------
 smooth_medoid_panel <- function(df_long, k = 10) {
   df_long %>%
     group_by(traj_cluster_mv, variable) %>%
     group_modify(~{
       d <- .x %>% arrange(hour_idx)
       
       # Fit only on finite pairs
       d_fit <- d %>% filter(is.finite(hour_idx), is.finite(value))
       
       # If too sparse, no smoothing
       n_ux <- length(unique(d_fit$hour_idx))
       if (nrow(d_fit) < 6 || n_ux < 6) {
         d$value_smooth <- d$value
         return(d)
       }
       
       # Choose k adaptively (must be < #unique x)
       k_use <- min(k, max(5, n_ux - 1))
       
       # Try GAM; if it fails, fall back to LOESS; if that fails, moving average
       yhat <- tryCatch({
         fit <- mgcv::gam(
           value ~ s(hour_idx, k = k_use, bs = "tp"),
           data = d_fit,
           method = "REML"
         )
         as.numeric(predict(fit, newdata = d))
       }, error = function(e1) {
         # LOESS fallback
         tryCatch({
           fit2 <- stats::loess(value ~ hour_idx, data = d_fit, span = 0.25,
                                control = loess.control(surface = "direct"))
           as.numeric(predict(fit2, newdata = d))
         }, error = function(e2) {
           # moving average fallback (centered)
           # NOTE: this requires zoo; you already load zoo earlier in your pipeline
           v <- d$value
           zoo::rollmean(v, k = 7, fill = NA, align = "center")
         })
       })
       
       d$value_smooth <- yhat
       d
     }) %>%
     ungroup()
 }
 
 
 
 medoid_icu_s <- smooth_medoid_panel(medoid_icu, k = 10)
 
 p_dtw_clin3_smooth <- ggplot(medoid_icu_s,
                              aes(hour_idx, value_smooth, color = traj_cluster_mv)) +
   geom_line(linewidth = 1.35) +
   facet_wrap(~ variable, scales = "free_y",
              labeller = as_labeller(var_labs[names(var_labs) %in% unique(medoid_icu_s$variable)])) +
   labs(
     x = "ICU hour",
     y = NULL,
     title = "DTW-derived ICU trajectory phenotypes (k=3)",
     subtitle = "Smoothed clinical medoid trajectories over the first 72 ICU hours"
   ) +
   theme_minimal(base_size = 12) +
   theme(legend.title = element_blank(), strip.text = element_text(face = "bold"))
 
 print(p_dtw_clin3_smooth)
 
 # ------------------------------------------------------------
 # 7) Save all k=3 plots
 # ------------------------------------------------------------
 dir.create("figs", showWarnings = FALSE)
 
 save_plot <- function(p, fname, width = 10, height = 6, dpi = 600) {
   stopifnot(inherits(p, "ggplot"))
   ggsave(
     filename = file.path("figs", paste0(fname, ".png")),
     plot = p, width = width, height = height, units = "in",
     dpi = dpi, device = ragg::agg_png, bg = "white", limitsize = FALSE
   )
   ggsave(
     filename = file.path("figs", paste0(fname, ".pdf")),
     plot = p, width = width, height = height, units = "in",
     device = cairo_pdf, bg = "white", limitsize = FALSE
   )
   message("Saved: figs/", fname, ".[png|pdf]")
 }
 
 save_plot(p_dtw_clin3,        "dtw_k3_clinical_medoids_72h", width = 10, height = 6)
 save_plot(p_dtw_clin3_smooth, "dtw_k3_clinical_medoids_72h_smoothed", width = 10, height = 6)
 save_plot(p_dtw_expo3,        "dtw_k3_exposure_medoids_90d", width = 10, height = 6)
 
 if (exists("p_dtw_outcomes3")) {
   save_plot(p_dtw_outcomes3,  "dtw_k3_outcomes_by_phenotype", width = 9, height = 5)
 }
 
 message("Done (k=3).")
 
 
 
 suppressPackageStartupMessages({
   library(dplyr); library(tidyr); library(purrr); library(ggplot2)
   library(mgcv); library(ragg)
 })
 
 # # ---------- Exposure summary per stay (90d) ----------
 # if (!exists("expo_sum_90")) {
 #   stopifnot(exists("expo_daily_365"))
 #   expo_sum_90 <- expo_daily_365 %>%
 #     transmute(
 #       hospitalization_id = as.character(hospitalization_id),
 #       no2_mean_90  = map_dbl(expo, ~ mean(.x$no2_raw,  na.rm = TRUE)),
 #       pm25_mean_90 = map_dbl(expo, ~ mean(.x$pm25_raw, na.rm = TRUE)),
 #       no2_med_90   = map_dbl(expo, ~ median(.x$no2_raw,  na.rm = TRUE)),
 #       pm25_med_90  = map_dbl(expo, ~ median(.x$pm25_raw, na.rm = TRUE))
 #     )
 # }
 # 
 # mk_qtile <- function(x, k = 4) {
 #   q <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
 #   q <- unique(q)
 #   list(q = q, k = length(q) - 1)
 # }
 # 
 # no2_cuts  <- mk_qtile(expo_sum_90$no2_mean_90,  k = 4)
 # pm25_cuts <- mk_qtile(expo_sum_90$pm25_mean_90, k = 4)
 # 
 # expo_sum_90 <- expo_sum_90 %>%
 #   mutate(
 #     no2_q  = cut(no2_mean_90,
 #                  breaks = no2_cuts$q,
 #                  include.lowest = TRUE,
 #                  labels = paste0("Q", seq_len(no2_cuts$k))),
 #     pm25_q = cut(pm25_mean_90,
 #                  breaks = pm25_cuts$q,
 #                  include.lowest = TRUE,
 #                  labels = paste0("Q", seq_len(pm25_cuts$k))),
 #     joint_hi = case_when(
 #       no2_q == paste0("Q", no2_cuts$k) & pm25_q == paste0("Q", pm25_cuts$k) ~ "High NO2 & High PM2.5",
 #       no2_q == "Q1" & pm25_q == "Q1" ~ "Low NO2 & Low PM2.5",
 #       TRUE ~ "Mixed"
 #     )
 #   )
 # 
 # # ---------- Build 72h panel ----------
 # stopifnot(exists("icu_hourly"))
 # plot_vars <- intersect(c("pf_ratio", "fio2_max", "peep_med", "spo2_med", "sig", "paco2_med"), names(icu_hourly))
 # stopifnot("pf_ratio" %in% plot_vars)
 # 
 # icu_panel_72 <- icu_hourly %>%
 #   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
 #   arrange(hospitalization_id, hour_ts) %>%
 #   group_by(hospitalization_id) %>%
 #   mutate(hour_idx = row_number()) %>%
 #   ungroup() %>%
 #   filter(hour_idx <= 72) %>%
 #   select(hospitalization_id, hour_idx, any_of(plot_vars)) %>%
 #   left_join(
 #     expo_sum_90 %>% select(hospitalization_id, no2_q, pm25_q, joint_hi),
 #     by = "hospitalization_id"
 #   )
 # 
 
 # ---------- Exposure summary per stay (365d) ----------
 if (!exists("expo_sum_365")) {
   stopifnot(exists("expo_daily_365"))
   expo_sum_365 <- expo_daily_365 %>%
     transmute(
       hospitalization_id = as.character(hospitalization_id),
       no2_mean_365  = map_dbl(expo, ~ mean(.x$no2_raw,  na.rm = TRUE)),
       pm25_mean_365 = map_dbl(expo, ~ mean(.x$pm25_raw, na.rm = TRUE)),
       no2_med_365   = map_dbl(expo, ~ median(.x$no2_raw,  na.rm = TRUE)),
       pm25_med_365  = map_dbl(expo, ~ median(.x$pm25_raw, na.rm = TRUE))
     )
 }
 
 mk_qtile <- function(x, k = 4) {
   q <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
   q <- unique(q)
   list(q = q, k = length(q) - 1)
 }
 
 no2_cuts  <- mk_qtile(expo_sum_365$no2_mean_365,  k = 3)
 pm25_cuts <- mk_qtile(expo_sum_365$pm25_mean_365, k = 3)
 
 expo_sum_365 <- expo_sum_365 %>%
   mutate(
     no2_q  = cut(no2_mean_365,
                  breaks = no2_cuts$q,
                  include.lowest = TRUE,
                  labels = paste0("Q", seq_len(no2_cuts$k))),
     pm25_q = cut(pm25_mean_365,
                  breaks = pm25_cuts$q,
                  include.lowest = TRUE,
                  labels = paste0("Q", seq_len(pm25_cuts$k))),
     joint_hi = case_when(
       no2_q == paste0("Q", no2_cuts$k) & pm25_q == paste0("Q", pm25_cuts$k) ~ "High NO2 & High PM2.5",
       no2_q == "Q1" & pm25_q == "Q1" ~ "Low NO2 & Low PM2.5",
       TRUE ~ "Mixed"
     )
   )
 
 # ---------- Build 72h panel ----------
 stopifnot(exists("icu_hourly"))
 plot_vars <- intersect(c("pf_ratio", "fio2_max", "peep_med", "spo2_med", "sig", "paco2_med"), names(icu_hourly))
 stopifnot("pf_ratio" %in% plot_vars)
 
 icu_panel_72 <- icu_hourly %>%
   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
   arrange(hospitalization_id, hour_ts) %>%
   group_by(hospitalization_id) %>%
   mutate(hour_idx = row_number()) %>%
   ungroup() %>%
   filter(hour_idx <= 72) %>%
   select(hospitalization_id, hour_idx, any_of(plot_vars)) %>%
   left_join(
     expo_sum_365 %>% select(hospitalization_id, no2_q, pm25_q, joint_hi),
     by = "hospitalization_id"
   )
 
 
 
 # ---------- Fixed, evidence-based cutpoints ----------
 # # NO2 in ppb (example bins): 0–2, 2–5, 5–8, 8+
 # no2_breaks <- c(-Inf, 2, 5, 8, Inf)
 # no2_labels <- c("0–2", "2–5", "5–8", "8+")
 # 
 # # PM2.5 in µg/m^3 (example bins): 0–2, 2–5, 5–8, 8+
 # pm25_breaks <- c(-Inf, 2, 5, 8, Inf)
 # pm25_labels <- c("0–2", "2–5", "5–8", "8+")
 # 
 # expo_sum_365 <- expo_sum_365 %>%
 #   mutate(
 #     no2_cat = cut(
 #       no2_mean_365,
 #       breaks = no2_breaks,
 #       labels = no2_labels,
 #       include.lowest = TRUE,
 #       right = FALSE  # intervals are [a, b)
 #     ),
 #     pm25_cat = cut(
 #       pm25_mean_365,
 #       breaks = pm25_breaks,
 #       labels = pm25_labels,
 #       include.lowest = TRUE,
 #       right = FALSE
 #     ),
 #     joint_hi = dplyr::case_when(
 #       no2_cat == "8+" & pm25_cat == "8+" ~ "High NO2 & High PM2.5",
 #       no2_cat == "0–2" & pm25_cat == "0–2" ~ "Low NO2 & Low PM2.5",
 #       TRUE ~ "Mixed"
 #     )
 #   )
 # 
 # # Optional: sanity check distribution
 # # table(expo_sum_365$no2_cat, useNA = "ifany")
 # # table(expo_sum_365$pm25_cat, useNA = "ifany")
 # # table(expo_sum_365$joint_hi, useNA = "ifany")
 # 
 # # ---------- Build 72h panel ----------
 # stopifnot(exists("icu_hourly"))
 # plot_vars <- intersect(c("pf_ratio", "fio2_max", "peep_med", "spo2_med", "sig", "paco2_med"), names(icu_hourly))
 # stopifnot("pf_ratio" %in% plot_vars)
 # 
 # icu_panel_72 <- icu_hourly %>%
 #   mutate(hospitalization_id = as.character(hospitalization_id)) %>%
 #   arrange(hospitalization_id, hour_ts) %>%
 #   group_by(hospitalization_id) %>%
 #   mutate(hour_idx = dplyr::row_number()) %>%
 #   ungroup() %>%
 #   filter(hour_idx <= 72) %>%
 #   select(hospitalization_id, hour_idx, any_of(plot_vars)) %>%
 #   left_join(
 #     expo_sum_365 %>% select(hospitalization_id, no2_cat, pm25_cat, joint_hi),
 #     by = "hospitalization_id"
 #   )
 # 
 # 
 
 summ_iqr_std <- function(x) {
   x <- x[is.finite(x)]
   if (!length(x)) return(c(mean = NA_real_, p25 = NA_real_, p75 = NA_real_, n = 0))
   c(
     mean = mean(x, na.rm = TRUE),
     p25  = unname(quantile(x, 0.25, na.rm = TRUE)),
     p75  = unname(quantile(x, 0.75, na.rm = TRUE)),
     n    = length(x)
   )
 }
 
 make_traj_summ <- function(df, y, grp, hour = "hour_idx") {
   stopifnot(all(c(y, grp, hour) %in% names(df)))
   df %>%
     filter(is.finite(.data[[y]]), !is.na(.data[[grp]]), is.finite(.data[[hour]])) %>%
     group_by(group = as.factor(.data[[grp]]), hour_idx = .data[[hour]]) %>%
     summarise(stats = list(summ_iqr_std(.data[[y]])), .groups = "drop") %>%
     unnest_wider(stats)
 }
 
 smooth_traj_gam <- function(traj_df, k = 10) {
   traj_df %>%
     group_by(group) %>%
     group_modify(~{
       d <- .x %>% arrange(hour_idx)
       if (sum(is.finite(d$mean)) < 6) {
         d$mean_smooth <- d$mean
         return(d)
       }
       fit <- mgcv::gam(mean ~ s(hour_idx, k = k), data = d, method = "REML")
       d$mean_smooth <- as.numeric(predict(fit, newdata = d))
       d
     }) %>%
     ungroup()
 }
 
 plot_traj_smoothed <- function(traj_df_s, title, ylab,
                                subtitle = "Smoothed mean (GAM); ribbon = empirical IQR",
                                min_n_per_hour = 30,
                                alpha_ribbon = 0.18,
                                lw = 1.3) {
   
   d <- traj_df_s %>%
     filter(is.finite(mean_smooth), is.finite(p25), is.finite(p75)) %>%
     filter(n >= min_n_per_hour)
   
   ggplot(d, aes(x = hour_idx, y = mean_smooth, color = group)) +
     geom_ribbon(aes(ymin = p25, ymax = p75, fill = group),
                 alpha = alpha_ribbon, color = NA) +
     geom_line(linewidth = lw) +
     labs(x = "ICU hour", y = ylab, title = title, subtitle = subtitle) +
     theme_minimal(base_size = 12) +
     theme(legend.title = element_blank())
 }
 
 # --- P/F by NO2 quartile ---
 traj_pf_no2   <- make_traj_summ(icu_panel_72, y = "pf_ratio", grp = "no2_q")
 traj_pf_no2_s <- smooth_traj_gam(traj_pf_no2, k = 10)
 
 p_pf_no2_smooth <- plot_traj_smoothed(
   traj_pf_no2_s,
   title = "Smoothed P/F ratio over first 72 ICU hours by chronic NO\u2082 quartile (365-day mean)",
   ylab  = "P/F ratio"
 )
 
 # --- P/F by PM2.5 quartile ---
 traj_pf_pm25   <- make_traj_summ(icu_panel_72, y = "pf_ratio", grp = "pm25_q")
 traj_pf_pm25_s <- smooth_traj_gam(traj_pf_pm25, k = 10)
 
 p_pf_pm25_smooth <- plot_traj_smoothed(
   traj_pf_pm25_s,
   title = "Smoothed P/F ratio over first 72 ICU hours by chronic PM\u2082.\u2085 quartile (365-day mean)",
   ylab  = "P/F ratio"
 )
 
 # --- P/F by JOINT exposure group (Low-Low / Mixed / High-High) ---
 traj_pf_joint   <- make_traj_summ(icu_panel_72, y = "pf_ratio", grp = "joint_hi")
 traj_pf_joint_s <- smooth_traj_gam(traj_pf_joint, k = 10)
 
 p_pf_joint_smooth <- plot_traj_smoothed(
   traj_pf_joint_s,
   title = "Smoothed P/F ratio over first 72 ICU hours by joint exposure group (365-day mean)",
   ylab  = "P/F ratio",
   min_n_per_hour = 30
 )
 
 print(p_pf_no2_smooth)
 print(p_pf_pm25_smooth)
 print(p_pf_joint_smooth)
 
 
 dir.create("figs", showWarnings = FALSE)
 
 save_plot <- function(p, filename, width = 8.5, height = 5.5, dpi = 600) {
   stopifnot(inherits(p, "ggplot"))
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".png")),
     plot = p, width = width, height = height, dpi = dpi, units = "in",
     device = ragg::agg_png, bg = "white", limitsize = FALSE
   )
   ggplot2::ggsave(
     filename = file.path("figs", paste0(filename, ".pdf")),
     plot = p, width = width, height = height, units = "in",
     device = cairo_pdf, bg = "white", limitsize = FALSE
   )
   message(sprintf("Saved: figs/%s.[png|pdf]", filename))
 }
 
 save_plot(p_pf_no2_smooth,   "traj72_pf_by_no2_quartile_smoothed")
 save_plot(p_pf_pm25_smooth,  "traj72_pf_by_pm25_quartile_smoothed")
 save_plot(p_pf_joint_smooth, "traj72_pf_by_joint_exposure_smoothed")
 
 
 
 
  
 


