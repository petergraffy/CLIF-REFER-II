# ================================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) | PI: Peter Graffy
# Trajectory Cohort Builder (CLIF) — Primary + Secondary Respiratory Support Cohorts
#
# Primary cohort   : Adults (>=18) with IMV start within +24h of first ICU admit (t0 = IMV start)
# Secondary cohort : Adults (>=18) with ADVANCED support start within +24h of first ICU admit
#                   (ADVANCED = IMV, NIPPV, CPAP, High Flow NC; t0 = first advanced start)
#
# Outputs (written to output/run_[SITE]_[DATE]/):
#   cohort_primary_imv72.csv            : 1 row / hospitalization (t0 = IMV start; traj window 0–72h)
#   cohort_secondary_adv72.csv          : 1 row / hospitalization (t0 = first advanced start; 0–72h)
#   exclusion_primary_imv72.csv         : exclusions + first failing reason (primary)
#   exclusion_secondary_adv72.csv       : exclusions + first failing reason (secondary)
#   flow_primary_imv72.csv              : flow counts (primary)
#   flow_secondary_adv72.csv            : flow counts (secondary)
#
# Notes:
# - No ABG/SpO2-density inclusion rule (avoid measurement-intensity selection bias for trajectories)
# - No ICU LOS >= 24h rule (allow early deaths/extubations; handle as censoring downstream)
# - Optional: require minimum respiratory_support density post-t0 (set MIN_RS_HOURS > 0)
# ================================================================================================

pkgs <- c("tidyverse", "lubridate", "data.table", "dtwclust", "RcppRoll", "zoo", "ggplot2")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(data.table)
  library(dtwclust)
  library(RcppRoll)
  library(zoo)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})

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


clif_paths <- found_map

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

clif_tables <- lapply(clif_paths, read_any)
cat("Loaded tables: ", paste(names(clif_tables), collapse = ", "), "\n")

cohort_use <- cohort_primary   # or cohort_secondary
H <- 72
bin_unit <- "hour"

feat_spec <- list(
  fio2 = list(col = "fio2_set", prefer = c("fio2_set")),
  peep = list(col = "peep", prefer = c("peep_set", "peep_obs"))
)

# Optional additional features once you see coverage:
feat_spec <- list(
  fio2 = list(col = "fio2", prefer = c("fio2_set")),
  peep = list(col = "peep", prefer = c("peep_set","peep_obs")),
  plat = list(col = "plat", prefer = c("plateau_pressure_obs")),
  mv   = list(col = "mv",   prefer = c("minute_vent_obs")),
  rr   = list(col = "rr",   prefer = c("resp_rate_obs","resp_rate_set"))
)


build_rs_hourly <- function(cohort_df, clif_tables, H = 72, bin_unit = "hour",
                            feat_spec, max_locf_gap_hours = 6) {
  
  rs <- clif_tables[["clif_respiratory_support"]] %>%
    rename_with(tolower)
  
  stopifnot(all(c("hospitalization_id","recorded_dttm") %in% names(rs)))
  stopifnot(all(c("hospitalization_id","t0") %in% names(cohort_df)))
  
  prefer_cols <- unique(unlist(lapply(feat_spec, `[[`, "prefer")))
  keep_cols <- unique(c("hospitalization_id","recorded_dttm", prefer_cols))
  keep_cols <- intersect(keep_cols, names(rs))
  
  rs <- rs %>% select(all_of(keep_cols)) %>%
    mutate(recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC"))
  
  df <- rs %>%
    inner_join(cohort_df %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
    mutate(t0 = as.POSIXct(t0, tz = "UTC"),
           t_end = t0 + dhours(H)) %>%
    filter(recorded_dttm >= t0, recorded_dttm <= t_end) %>%
    mutate(hour = floor_date(recorded_dttm, unit = bin_unit),
           t = as.integer(difftime(hour, t0, units = "hours"))) %>%
    filter(t >= 0, t <= H)
  
  # ---- robust fallback feature builder ----
  for (nm in names(feat_spec)) {
    prefs <- feat_spec[[nm]]$prefer
    prefs <- prefs[prefs %in% names(df)]
    if (length(prefs) == 0) {
      df[[nm]] <- NA_real_
    } else if (length(prefs) == 1) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[prefs]]))
    } else {
      cols <- lapply(prefs, function(p) df[[p]])
      df[[nm]] <- suppressWarnings(as.numeric(Reduce(dplyr::coalesce, cols)))
    }
  }
  
  hourly <- df %>%
    group_by(hospitalization_id, t0, hour, t) %>%
    summarize(across(all_of(names(feat_spec)), ~ median(.x, na.rm = TRUE)),
              .groups = "drop")
  
  grid <- cohort_df %>%
    select(hospitalization_id, t0) %>%
    distinct() %>%
    mutate(t0 = as.POSIXct(t0, tz = "UTC")) %>%
    tidyr::expand_grid(t = 0:H) %>%
    mutate(hour = t0 + dhours(t))
  
  hourly_full <- grid %>%
    left_join(hourly, by = c("hospitalization_id","t0","t","hour")) %>%
    arrange(hospitalization_id, t)
  
  fill_one <- function(x) {
    x_f <- zoo::na.locf(x, na.rm = FALSE)
    x_f <- zoo::na.locf(x_f, fromLast = TRUE, na.rm = FALSE)
    
    is_na <- is.na(x)
    r <- rle(is_na)
    idx <- inverse.rle(list(values = seq_along(r$lengths), lengths = r$lengths))
    long_na_pos <- is_na & (r$lengths[idx] > max_locf_gap_hours)
    x_f[long_na_pos] <- NA_real_
    x_f
  }
  
  hourly_filled <- hourly_full %>%
    group_by(hospitalization_id) %>%
    mutate(across(all_of(names(feat_spec)), fill_one)) %>%
    ungroup()
  
  hourly_filled
}


rs_hourly <- build_rs_hourly(
  cohort_df = cohort_use,
  clif_tables = clif_tables,
  H = H,
  bin_unit = bin_unit,
  feat_spec = feat_spec,
  max_locf_gap_hours = 6
)

rs_hourly %>% glimpse()

qc <- rs_hourly %>%
  group_by(hospitalization_id) %>%
  summarize(
    n_t = n(),
    fio2_obs = sum(!is.na(fio2)),
    peep_obs = sum(!is.na(peep)),
    plat_obs = sum(!is.na(plat)),
    mv_obs = sum(!is.na(mv)),
    rr_obs = sum(!is.na(rr)),
    frac_complete = mean(!is.na(fio2) & !is.na(peep) & !is.na(plat) & !is.na(mv) & !is.na(rr)),
    .groups = "drop"
  )

summary(qc$frac_complete)

# Recommended filter: >= 0.6 of hours have both (fio2, peep)
keep_ids <- qc %>% filter(frac_complete >= 0.60) %>% pull(hospitalization_id)

rs_hourly_keep <- rs_hourly %>% filter(hospitalization_id %in% keep_ids)
length(unique(rs_hourly_keep$hospitalization_id))

make_series_list <- function(rs_hourly_long, feats = c("fio2","peep"), H = 72) {
  
  # global scaling parameters (across all patients/timepoints)
  mu <- rs_hourly_long %>% summarize(across(all_of(feats), ~ mean(.x, na.rm = TRUE)))
  sdv <- rs_hourly_long %>% summarize(across(all_of(feats), ~ sd(.x, na.rm = TRUE)))
  
  rs_scaled <- rs_hourly_long %>%
    mutate(across(all_of(feats), ~ (.x - as.numeric(mu[[cur_column()]])) /
                    as.numeric(sdv[[cur_column()]]) ))
  
  # split into list of matrices, each (H+1) x p
  ids <- sort(unique(rs_scaled$hospitalization_id))
  
  series <- lapply(ids, function(id) {
    m <- rs_scaled %>%
      filter(hospitalization_id == id) %>%
      arrange(t) %>%
      select(all_of(feats)) %>%
      as.matrix()
    # enforce expected length
    if (nrow(m) != (H + 1)) {
      # pad/truncate safely
      m2 <- matrix(NA_real_, nrow = H + 1, ncol = length(feats))
      colnames(m2) <- feats
      rr <- min(nrow(m), H + 1)
      m2[1:rr, ] <- m[1:rr, , drop = FALSE]
      m <- m2
    }
    m
  })
  names(series) <- ids
  list(series = series, mu = mu, sd = sdv)
}

# =========================
# DTW-safe filter + impute
# =========================

# 1) choose features (must match your build_rs_hourly feat_spec output names)
feats <- c("fio2", "peep", "plat", "mv", "rr")

# 2) QC missingness by hospitalization
na_qc <- rs_hourly_keep %>%
  group_by(hospitalization_id) %>%
  summarize(
    # any NA across the feature set at a given hour
    frac_any_na = mean(if_any(all_of(feats), is.na)),
    # all NA across the feature set at a given hour (completely missing row)
    frac_all_na = mean(if_all(all_of(feats), is.na)),
    .groups = "drop"
  )

# 3) filter rule (tune as needed)
#    keep patients where <=40% of hourly rows have any NA across feats
MAX_FRAC_ANY_NA <- 0.40

keep_ids2 <- na_qc %>%
  filter(frac_any_na <= MAX_FRAC_ANY_NA) %>%
  pull(hospitalization_id)

rs_hourly_dtw <- rs_hourly_keep %>%
  filter(hospitalization_id %in% keep_ids2)

cat("DTW input N (after missingness filter): ",
    length(unique(rs_hourly_dtw$hospitalization_id)), "\n", sep = "")

# 4) DTW-safe imputation:
#    fill remaining NA with within-patient median (fallback to global median)
impute_for_dtw <- function(df, feats) {
  # global medians as fallback
  gmed <- df %>% summarize(across(all_of(feats), ~ median(.x, na.rm = TRUE)))
  
  df %>%
    group_by(hospitalization_id) %>%
    mutate(
      across(all_of(feats), \(x) {
        pmed <- suppressWarnings(median(x, na.rm = TRUE))
        if (is.na(pmed)) pmed <- as.numeric(gmed[[cur_column()]])
        replace(x, is.na(x), pmed)
      })
    ) %>%
    ungroup()
}

rs_hourly_dtw_imp <- impute_for_dtw(rs_hourly_dtw, feats)

# 5) sanity checks: confirm there are no NA left in features
na_counts <- rs_hourly_dtw_imp %>%
  summarize(across(all_of(feats), ~ sum(is.na(.x))))
print(na_counts)

stopifnot(all(as.numeric(na_counts[1,]) == 0))

# =========================
# Rebuild series list & run DTW clustering
# =========================

# IMPORTANT: rebuild from the *imputed* data
prep2 <- make_series_list(rs_hourly_dtw_imp, feats = feats, H = H)
series_list2 <- prep2$series

# Final check: any NA inside any series matrix?
anyNA_series <- function(mats) any(vapply(mats, anyNA, logical(1)))
stopifnot(!anyNA_series(series_list2))

# 6) DTW clustering
set.seed(1)
k <- 4

cl <- tsclust(
  series_list2,
  type     = "partitional",
  k        = k,
  distance = "dtw_basic",
  centroid = "dba",
  trace    = TRUE,
  args = tsclust_args(
    dist = list(window.size = 12),  # Sakoe-Chiba band (hours); tune 6–18
    cent = list(max.iter = 20)
  )
)

print(table(cl@cluster))

# Map cluster to ids
cluster_df <- tibble(
  hospitalization_id = as.character(names(series_list2)),
  cluster = as.integer(cl@cluster)
)

# Join to original (unscaled) rs_hourly_keep for plotting
plot_df <- rs_hourly_keep %>%
  inner_join(cluster_df, by = "hospitalization_id")

# Plot median trajectory per cluster
plot_proto <- function(df, y, ylab, ylim_min, ylim_max) {
  ggplot(df, aes(x = t, y = .data[[y]], group = hospitalization_id)) +
    geom_line(alpha = 0.03, color = "grey50") +
    stat_summary(aes(group = 1),
                 fun = median,
                 geom = "line",
                 linewidth = 1.3,
                 color = "black") +
    facet_wrap(~cluster, ncol = 2) +
    coord_cartesian(ylim = c(ylim_min, ylim_max)) +
    labs(x = "Hours from t0",
         y = ylab,
         title = paste0("Trajectories: ", ylab)) +
    theme_minimal(base_size = 14)
}

p1<-plot_proto(plot_df, "fio2", "FiO2", 0.21, 1.0)
p2<-plot_proto(plot_df, "peep", "PEEP (cmH2O)", 0, 20)
p3<-plot_proto(plot_df, "plat", "Plateau Pressure (cmH2O)", 0, 40)
p4<-plot_proto(plot_df, "mv", "Minute Ventilation (L/min)", 4, 18)
p5<-plot_proto(plot_df, "rr", "Respiratory Rate (breaths/min)", 10, 40)

top4 <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv")
fig5 <- cowplot::plot_grid(top4, p5, ncol = 1, rel_heights = c(2, 1))

# optional title
fig5_titled <- cowplot::ggdraw(fig5) +
  cowplot::draw_label("ARF Severity Trajectories by Cluster (0–72h from t0)",
                      x = 0.5, y = 0.99, hjust = 0.5, vjust = -3,
                      fontface = "bold", size = 16)

print(fig5_titled)

# ---- 4) Save ----
ggsave("traj_5panel_clusters.png", fig5_titled, width = 14, height = 15, dpi = 300)

cohort_primary_clusters <- cohort_use %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(cluster_df, by = "hospitalization_id")

cohort_primary_clusters %>% count(cluster)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# ---- choose features you clustered/plot ----
feats_main <- c("fio2","peep","plat","mv","rr")

# helper: compute slope of y ~ t
slope_lm <- function(t, y) {
  ok <- is.finite(t) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  coef(lm(y[ok] ~ t[ok]))[2]
}

traj_summary <- rs_hourly_dtw_imp %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  inner_join(cluster_df %>% mutate(hospitalization_id = as.character(hospitalization_id)),
             by = "hospitalization_id") %>%
  group_by(hospitalization_id, cluster) %>%
  summarize(
    # baseline at t0 (t==0)
    across(all_of(feats_main), ~ median(.x[t == 0], na.rm = TRUE), .names = "t0_{.col}"),
    
    # early window 0-12h
    across(all_of(feats_main), ~ median(.x[t >= 0 & t <= 12], na.rm = TRUE), .names = "h0_12_med_{.col}"),
    across(all_of(feats_main), ~ max(.x[t >= 0 & t <= 12], na.rm = TRUE),    .names = "h0_12_max_{.col}"),
    
    # full window 0-72h
    across(all_of(feats_main), ~ median(.x, na.rm = TRUE), .names = "h0_72_med_{.col}"),
    across(all_of(feats_main), ~ max(.x, na.rm = TRUE),    .names = "h0_72_max_{.col}"),
    
    # slopes (overall + early)
    across(all_of(feats_main), ~ slope_lm(t, .x),                 .names = "slope_0_72_{.col}"),
    across(all_of(feats_main), ~ slope_lm(t[t <= 12], .x[t <= 12]), .names = "slope_0_12_{.col}"),
    
    # data density (how much is “real” vs imputed is hard to know now, but we can track variability)
    n_hours = n(),
    .groups = "drop"
  )

cluster_profile <- traj_summary %>%
  group_by(cluster) %>%
  summarize(
    n = n(),
    across(starts_with("t0_"),      ~ median(.x, na.rm = TRUE), .names = "{.col}_med"),
    across(starts_with("h0_72_med"),~ median(.x, na.rm = TRUE), .names = "{.col}_med"),
    across(starts_with("h0_72_max"),~ median(.x, na.rm = TRUE), .names = "{.col}_med"),
    across(starts_with("slope_0_72"),~ median(.x, na.rm = TRUE), .names = "{.col}_med"),
    .groups = "drop"
  ) %>%
  arrange(cluster)

print(cluster_profile)

hosp <- clif_tables[["clif_hospitalization"]] %>% rename_with(tolower)

# try common CLIF fields; adjust if yours differ
outcome_cols <- intersect(
  c("hospitalization_id","discharge_dttm","discharge_category"),
  names(hosp)
)

outcomes <- hosp %>%
  select(any_of(outcome_cols)) %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  mutate(
    died = case_when(
      "discharge_category" %in% names(.) ~ tolower(discharge_category) %in% c("expired", "hospice"),
      TRUE ~ NA
    )
  )

cluster_outcomes <- traj_summary %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(outcomes, by = "hospitalization_id") %>%
  group_by(cluster) %>%
  summarize(
    n = n(),
    died_n = sum(died %in% TRUE, na.rm = TRUE),
    died_pct = round(100 * mean(died %in% TRUE, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(cluster_outcomes)

vitals <- clif_tables[["clif_vitals"]]
stopifnot(all(c("hospitalization_id","recorded_dttm","vital_category","vital_value") %in% names(vitals)))

# Pull SpO2 only
spo2_hourly <- vitals %>%
  filter(tolower(vital_category) == "spo2") %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC"),
    spo2 = suppressWarnings(as.numeric(vital_value))
  ) %>%
  inner_join(cohort_use %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                      t0 = as.POSIXct(t0, tz="UTC")),
             by = "hospitalization_id") %>%
  mutate(
    hour = floor_date(recorded_dttm, "hour"),
    t = as.integer(difftime(hour, t0, units = "hours"))
  ) %>%
  filter(t >= 0, t <= H) %>%
  group_by(hospitalization_id, t) %>%
  summarize(spo2 = median(spo2, na.rm = TRUE), .groups = "drop")

# Merge spo2 + compute S/F
rs_plus_spo2 <- rs_hourly_dtw_imp %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(spo2_hourly, by = c("hospitalization_id","t")) %>%
  mutate(
    # SpO2 fraction (0-1) if it was 90-100; keep NA otherwise
    spo2_frac = ifelse(!is.na(spo2) & spo2 > 1.5, spo2/100, spo2),
    sf_ratio = ifelse(!is.na(spo2_frac) & !is.na(fio2) & fio2 > 0, spo2_frac / fio2, NA_real_)
  )
































