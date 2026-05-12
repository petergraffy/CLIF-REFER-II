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
  library(ggalluvial)
})

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

repo_root <- if (exists("repo")) repo else getwd()
trajectory_out_dir <- if (exists("out_dir")) out_dir else file.path(repo_root, "output", "trajectory_run")
trajectory_fig_dir <- file.path(trajectory_out_dir, "figures")
dir.create(trajectory_fig_dir, recursive = TRUE, showWarnings = FALSE)

cohort_use <- if (exists("cohort_arf")) cohort_arf else cohort_secondary
H <- 72
bin_unit <- "hour"

feat_spec <- list(
  fio2 = list(col = "fio2", prefer = c("fio2_set")),
  peep = list(col = "peep", prefer = c("peep_set","peep_obs")),
  pplat = list(col = "pplat", prefer = c("plateau_pressure_obs")),
  pip = list(col = "pip", prefer = c("peak_inspiratory_pressure_obs","peak_inspiratory_pressure_set")),
  mapaw = list(col = "mapaw", prefer = c("mean_airway_pressure_obs")),
  mv   = list(col = "mv",   prefer = c("minute_vent_obs")),
  rr   = list(col = "rr",   prefer = c("resp_rate_obs","resp_rate_set")),
  vt   = list(col = "vt",   prefer = c("tidal_volume_obs","tidal_volume_set")),
  pc   = list(col = "pc",   prefer = c("pressure_control_set")),
  ps   = list(col = "ps",   prefer = c("pressure_support_set")),
  lpm  = list(col = "lpm",  prefer = c("lpm_set")),
  flow_rate = list(col = "flow_rate", prefer = c("flow_rate_set"))
)


build_rs_hourly <- function(cohort_df, clif_tables, H = 72, bin_unit = "hour",
                            feat_spec, max_locf_gap_hours = 6) {
  
  rs <- clif_tables[["clif_respiratory_support"]] %>%
    rename_with(tolower)
  
  stopifnot(all(c("hospitalization_id","recorded_dttm") %in% names(rs)))
  stopifnot(all(c("hospitalization_id","t0") %in% names(cohort_df)))
  
  prefer_cols <- unique(unlist(lapply(feat_spec, `[[`, "prefer")))
  keep_cols <- unique(c("hospitalization_id","recorded_dttm","device_category","mode_category", prefer_cols))
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

  df <- df %>%
    mutate(
      support_level = case_when(
        device_category == "Room Air" ~ 0,
        device_category == "Nasal Cannula" ~ 1,
        device_category %in% c("Face Mask", "Trach Collar", "Other") ~ 2,
        device_category %in% c("High Flow NC", "CPAP", "NIPPV") ~ 3,
        device_category == "IMV" ~ 4,
        TRUE ~ NA_real_
      ),
      any_imv = as.integer(device_category == "IMV"),
      any_advanced_support = as.integer(device_category %in% c("IMV", "NIPPV", "CPAP", "High Flow NC")),
      positive_pressure = as.integer(device_category %in% c("IMV", "NIPPV", "CPAP") |
                                       mode_category %in% c("Assist Control-Volume Control",
                                                            "Pressure Control",
                                                            "Pressure-Regulated Volume Control",
                                                            "SIMV",
                                                            "Pressure Support/CPAP"))
    )
  
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
    if (nm == "fio2") df[[nm]] <- if_else(df[[nm]] > 1.5, df[[nm]] / 100, df[[nm]])
  }
  
  hourly <- df %>%
    group_by(hospitalization_id, t0, hour, t) %>%
    summarize(
      support_level = max(support_level, na.rm = TRUE),
      any_imv = as.integer(any(any_imv == 1L, na.rm = TRUE)),
      any_advanced_support = as.integer(any(any_advanced_support == 1L, na.rm = TRUE)),
      positive_pressure = as.integer(any(positive_pressure == 1L, na.rm = TRUE)),
      across(all_of(names(feat_spec)), ~ median(.x, na.rm = TRUE)),
              .groups = "drop")
  hourly <- hourly %>%
    mutate(support_level = if_else(is.infinite(support_level), NA_real_, support_level))
  
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
    mutate(across(c("support_level", all_of(names(feat_spec))), fill_one),
           across(c(any_imv, any_advanced_support, positive_pressure), ~ replace_na(.x, 0L))) %>%
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

build_vitals_hourly <- function(cohort_df, clif_tables, H = 72) {
  vitals <- clif_tables[["clif_vitals"]] %>%
    rename_with(tolower)
  
  needed <- c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value")
  if (!all(needed %in% names(vitals))) return(tibble())
  
  vitals %>%
    filter(vital_category %in% c("spo2", "respiratory_rate", "heart_rate", "map")) %>%
    transmute(
      hospitalization_id = as.character(hospitalization_id),
      recorded_dttm = as.POSIXct(recorded_dttm, tz = "UTC"),
      vital_category,
      vital_value = suppressWarnings(as.numeric(vital_value))
    ) %>%
    inner_join(cohort_df %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                       t0 = as.POSIXct(t0, tz = "UTC")),
               by = "hospitalization_id") %>%
    mutate(t = as.integer(difftime(floor_date(recorded_dttm, "hour"), t0, units = "hours"))) %>%
    filter(t >= 0, t <= H) %>%
    group_by(hospitalization_id, t, vital_category) %>%
    summarize(value = median(vital_value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = vital_category, values_from = value, names_prefix = "vital_")
}

build_labs_hourly <- function(cohort_df, clif_tables, H = 72) {
  labs <- clif_tables[["clif_labs"]] %>%
    rename_with(tolower)
  
  needed <- c("hospitalization_id", "lab_result_dttm", "lab_category")
  if (!all(needed %in% names(labs))) return(tibble())
  
  lab_val_col <- intersect(c("lab_value_numeric", "lab_value"), names(labs))[1]
  if (is.na(lab_val_col)) return(tibble())
  
  lung_labs <- c(
    "po2_arterial", "pco2_arterial", "ph_arterial", "so2_arterial",
    "pco2_venous", "ph_venous", "so2_mixed_venous", "so2_central_venous",
    "bicarbonate", "lactate", "hemoglobin", "wbc", "platelet_count",
    "crp", "procalcitonin", "ferritin", "ldh"
  )
  
  labs %>%
    filter(lab_category %in% lung_labs) %>%
    transmute(
      hospitalization_id = as.character(hospitalization_id),
      lab_result_dttm = as.POSIXct(lab_result_dttm, tz = "UTC"),
      lab_category,
      lab_value = suppressWarnings(as.numeric(.data[[lab_val_col]]))
    ) %>%
    inner_join(cohort_df %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                       t0 = as.POSIXct(t0, tz = "UTC")),
               by = "hospitalization_id") %>%
    mutate(t = as.integer(difftime(floor_date(lab_result_dttm, "hour"), t0, units = "hours"))) %>%
    filter(t >= 0, t <= H) %>%
    group_by(hospitalization_id, t, lab_category) %>%
    summarize(value = median(lab_value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = lab_category, values_from = value, names_prefix = "lab_")
}

vitals_hourly <- build_vitals_hourly(cohort_use, clif_tables, H)
labs_hourly <- build_labs_hourly(cohort_use, clif_tables, H)

add_missing_cols <- function(df, cols) {
  for (col in setdiff(cols, names(df))) df[[col]] <- NA_real_
  df
}

rs_hourly <- rs_hourly %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(vitals_hourly, by = c("hospitalization_id", "t")) %>%
  left_join(labs_hourly, by = c("hospitalization_id", "t")) %>%
  add_missing_cols(c("vital_spo2", "lab_po2_arterial", "lab_pco2_arterial", "lab_ph_arterial")) %>%
  mutate(
    spo2_frac = if_else(!is.na(vital_spo2) & vital_spo2 > 1.5, vital_spo2 / 100, vital_spo2),
    sf_ratio = if_else(!is.na(spo2_frac) & !is.na(fio2) & fio2 > 0, spo2_frac / fio2, NA_real_),
    pf_ratio = if_else(!is.na(lab_po2_arterial) & !is.na(fio2) & fio2 > 0, lab_po2_arterial / fio2, NA_real_),
    ventilatory_acidosis = as.numeric(!is.na(lab_pco2_arterial) & !is.na(lab_ph_arterial) &
                                        lab_pco2_arterial >= 45 & lab_ph_arterial < 7.35)
  )

rs_hourly %>% glimpse()

candidate_feats <- c(
  "support_level", "any_imv", "any_advanced_support", "positive_pressure",
  "fio2", "peep", "pplat", "pip", "mapaw", "mv", "rr", "vt", "pc", "ps",
  "lpm", "flow_rate", "vital_spo2", "vital_respiratory_rate", "vital_heart_rate",
  "vital_map", "sf_ratio", "pf_ratio", "lab_po2_arterial", "lab_pco2_arterial",
  "lab_ph_arterial", "lab_so2_arterial", "lab_bicarbonate", "lab_lactate",
  "lab_hemoglobin", "lab_wbc", "lab_platelet_count", "lab_crp",
  "lab_procalcitonin", "lab_ferritin", "lab_ldh", "ventilatory_acidosis"
)
candidate_feats <- intersect(candidate_feats, names(rs_hourly))

feature_coverage <- rs_hourly %>%
  summarize(across(all_of(candidate_feats), ~ mean(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "coverage") %>%
  arrange(desc(coverage), feature)

print(feature_coverage)
write_csv(feature_coverage, file.path(trajectory_out_dir, "trajectory_feature_coverage.csv"))

required_feats <- c("support_level", "fio2", "vital_spo2", "sf_ratio")
coverage_selected <- feature_coverage %>%
  filter(coverage >= 0.20) %>%
  pull(feature)
feats <- unique(c(required_feats[required_feats %in% candidate_feats], coverage_selected))

if (length(feats) < 3) {
  stop("Fewer than 3 trajectory features passed coverage checks. Review CLIF mappings and feature coverage output.")
}

cat("DTW feature set: ", paste(feats, collapse = ", "), "\n", sep = "")

qc <- rs_hourly %>%
  group_by(hospitalization_id) %>%
  summarize(
    n_t = n(),
    frac_any_na = mean(if_any(all_of(feats), is.na)),
    frac_all_na = mean(if_all(all_of(feats), is.na)),
    .groups = "drop"
  )

summary(qc$frac_any_na)

MAX_FRAC_ANY_NA <- 0.60
keep_ids <- qc %>% filter(frac_any_na <= MAX_FRAC_ANY_NA) %>% pull(hospitalization_id)

rs_hourly_keep <- rs_hourly %>% filter(hospitalization_id %in% keep_ids)
length(unique(rs_hourly_keep$hospitalization_id))

make_series_list <- function(rs_hourly_long, feats = c("fio2","peep"), H = 72) {
  
  # global scaling parameters (across all patients/timepoints)
  mu <- rs_hourly_long %>% summarize(across(all_of(feats), ~ mean(.x, na.rm = TRUE)))
  sdv <- rs_hourly_long %>% summarize(across(all_of(feats), ~ sd(.x, na.rm = TRUE)))
  sdv <- sdv %>% mutate(across(everything(), ~ if_else(is.na(.x) | .x == 0, 1, .x)))
  
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
p3<-plot_proto(plot_df, "pplat", "Plateau Pressure (cmH2O)", 0, 40)
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
ggsave(file.path(trajectory_fig_dir, "traj_5panel_clusters.png"),
       fig5_titled, width = 14, height = 15, dpi = 300)

cohort_clusters <- cohort_use %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(cluster_df, by = "hospitalization_id")

cohort_clusters %>% count(cluster)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# ---- choose features you clustered/plot ----
feats_main <- feats

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

plot_df2 <- rs_plus_spo2 %>%
  inner_join(cluster_df %>% mutate(hospitalization_id = as.character(hospitalization_id)),
             by = "hospitalization_id")

# Examples:
p6<-plot_proto(plot_df2, "spo2_frac", "SpO2 (fraction)", 0.70, 1.00)
p7<-plot_proto(plot_df2, "sf_ratio", "S/F (SpO2/FiO2)", 1.0, 4.5)

# Find likely category names at your site first (one-time):
sort(unique(tolower(vitals$vital_category)))[1:50]

# Adjust these if your site uses different strings
ht_cats <- c("height_cm","height")
wt_cats <- c("weight_kg","weight")

hw_baseline <- vitals %>%
  mutate(vital_category = tolower(vital_category)) %>%
  filter(vital_category %in% c(ht_cats, wt_cats)) %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = as.POSIXct(recorded_dttm, tz="UTC"),
    vital_category,
    value = suppressWarnings(as.numeric(vital_value))
  ) %>%
  inner_join(cohort_use %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                      t0 = as.POSIXct(t0, tz="UTC")),
             by = "hospitalization_id") %>%
  mutate(dt_h = abs(as.numeric(difftime(recorded_dttm, t0, units="hours")))) %>%
  filter(dt_h <= 24) %>%
  group_by(hospitalization_id, vital_category) %>%
  slice_min(order_by = dt_h, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(vital_category = case_when(
    vital_category %in% ht_cats ~ "height_cm",
    vital_category %in% wt_cats ~ "weight_kg",
    TRUE ~ vital_category
  )) %>%
  select(hospitalization_id, vital_category, value) %>%
  pivot_wider(names_from = vital_category, values_from = value)

# join onto traj_summary or cohort
traj_summary2 <- traj_summary %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(hw_baseline, by = "hospitalization_id")

labs <- clif_tables[["clif_labs"]] %>% rename_with(tolower)
stopifnot(all(c("hospitalization_id","lab_result_dttm","lab_category") %in% names(labs)))

# pick numeric value column name (varies by export)
lab_val_col <- intersect(c("lab_value_numeric","lab_value"), names(labs))[1]
stopifnot(!is.na(lab_val_col))

pao2_hourly <- labs %>%
  filter(tolower(lab_category) %in% c("po2_arterial","pao2","po2")) %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm = as.POSIXct(lab_result_dttm, tz="UTC"),
    pao2 = suppressWarnings(as.numeric(.data[[lab_val_col]]))
  ) %>%
  inner_join(cohort_use %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                      t0 = as.POSIXct(t0, tz="UTC")),
             by = "hospitalization_id") %>%
  mutate(
    hour = floor_date(lab_result_dttm, "hour"),
    t = as.integer(difftime(hour, t0, units="hours"))
  ) %>%
  filter(t >= 0, t <= H) %>%
  group_by(hospitalization_id, t) %>%
  summarize(pao2 = median(pao2, na.rm = TRUE), .groups = "drop")

rs_plus_abg <- rs_plus_spo2 %>%
  left_join(pao2_hourly, by = c("hospitalization_id","t")) %>%
  mutate(
    pf_ratio = ifelse(!is.na(pao2) & !is.na(fio2) & fio2 > 0, pao2 / fio2, NA_real_)
  )

# plot PF only if it looks sane
plot_df3 <- rs_plus_abg %>%
  inner_join(cluster_df %>% mutate(hospitalization_id = as.character(hospitalization_id)),
             by = "hospitalization_id")

p8<-plot_proto(plot_df3, "pao2", "PaO2 (mmHg)", 40, 200)
p9<-plot_proto(plot_df3, "pf_ratio", "P/F ratio", 50, 400)




cut_points   <- c(0, 12, 24, 48, 72)
block_labels <- c("0-12","12-24","24-48","48-72")

collapse_device <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x) | !nzchar(trimws(x)), NA_character_, x)
  case_when(
    x %in% c("Nasal Cannula", "Face Mask") ~ "Low Flow O2",
    x %in% c("CPAP", "NIPPV")             ~ "NIV (CPAP/NIPPV)",
    x %in% c("High Flow NC")              ~ "HFNC",
    x %in% c("IMV")                       ~ "IMV",
    x %in% c("Trach Collar")              ~ "Trach Collar",
    x %in% c("Room Air")                  ~ "Room Air",
    TRUE                                  ~ x
  )
}

# -----------------------------
# 1) Build death time using discharge_category + discharge_dttm
# -----------------------------
hosp <- clif_tables[["clif_hospitalization"]] %>% rename_with(tolower)

stopifnot(all(c("hospitalization_id", "discharge_dttm", "discharge_category") %in% names(hosp)))

death_map <- hosp %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    discharge_dttm = as.POSIXct(discharge_dttm, tz = "UTC"),
    discharge_category = tolower(trimws(as.character(discharge_category)))
  ) %>%
  mutate(
    died_or_hospice = discharge_category %in% c("expired", "hospice"),
    death_dttm = ifelse(died_or_hospice, discharge_dttm, as.POSIXct(NA, tz = "UTC"))
  ) %>%
  mutate(death_dttm = as.POSIXct(death_dttm, origin = "1970-01-01", tz = "UTC")) %>%
  select(hospitalization_id, death_dttm, died_or_hospice)

# -----------------------------
# 2) Add death timing (hours from t0)
# -----------------------------
cohort_key <- cohort_use %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    t0 = as.POSIXct(t0, tz = "UTC")
  ) %>%
  left_join(death_map, by = "hospitalization_id") %>%
  mutate(
    death_h = ifelse(!is.na(death_dttm),
                     as.numeric(difftime(death_dttm, t0, units = "hours")),
                     NA_real_)
  )


grid_blocks <- cohort_key %>%
  tidyr::expand_grid(block = factor(block_labels, levels = block_labels)) %>%
  mutate(
    block_start = c(0, 12, 24, 48)[match(block, block_labels)],
    block_end   = c(12, 24, 48, 72)[match(block, block_labels)],
    # death occurs within this block
    death_in_block = !is.na(death_h) & death_h >= block_start & death_h < block_end,
    # already dead at the start of this block
    dead_before_block = !is.na(death_h) & death_h < block_start
  )

# -----------------------------
# 3) Dominant respiratory support state per block (if any rows exist)
# -----------------------------
rs_raw <- clif_tables[["clif_respiratory_support"]] %>% rename_with(tolower)

rs_block_dom <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = as.POSIXct(recorded_dttm, tz="UTC"),
    device_category = collapse_device(device_category)
  ) %>%
  inner_join(cohort_key %>% select(hospitalization_id, t0), by = "hospitalization_id") %>%
  mutate(t = as.numeric(difftime(recorded_dttm, t0, units="hours"))) %>%
  filter(t >= 0, t < max(cut_points)) %>%
  mutate(
    block = cut(
      t,
      breaks = cut_points,
      include.lowest = TRUE,
      right = FALSE,
      labels = block_labels
    ) |> factor(levels = block_labels),
    device_category = ifelse(is.na(device_category), "Unknown device", device_category)
  ) %>%
  filter(!is.na(block)) %>%
  count(hospitalization_id, block, device_category, name = "n") %>%
  group_by(hospitalization_id, block) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

# -----------------------------
# 4) Combine: death overrides device states
# -----------------------------
rs_blocks <- grid_blocks %>%
  left_join(rs_block_dom, by = c("hospitalization_id","block")) %>%
  mutate(
    device_state = case_when(
      dead_before_block ~ "Dead",               # absorbing after death
      death_in_block    ~ "Death",              # event occurs during block
      is.na(device_category) ~ "No RS record",  # no RS rows in block
      TRUE ~ device_category
    )
  ) %>%
  select(hospitalization_id, block, device_state) %>%
  pivot_wider(names_from = block, values_from = device_state)

# Join clusters
rs_blocks_cl <- rs_blocks %>%
  left_join(cluster_df %>% mutate(hospitalization_id = as.character(hospitalization_id)),
            by = "hospitalization_id") %>%
  filter(!is.na(cluster))

alluv <- rs_blocks_cl %>%
  count(cluster, `0-12`, `12-24`, `24-48`, `48-72`, name = "n")

# -----------------------------
# 5) Plot
# -----------------------------
ggplot(alluv,
       aes(axis1 = `0-12`, axis2 = `12-24`, axis3 = `24-48`, axis4 = `48-72`, y = n)) +
  geom_alluvium(aes(fill = `0-12`), alpha = 0.7, width = 1/12) +
  geom_stratum(width = 1/12, color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = block_labels, expand = c(.05, .05)) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    x = "Hours from t0",
    y = "Patients",
    title = "Respiratory support transitions by DTW cluster (death integrated)",
    subtitle = "Death is treated as a competing event (Death during block; Dead thereafter)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


# Store plot object first (important)
p_sankey <- ggplot(alluv,
                   aes(axis1 = `0-12`, axis2 = `12-24`, axis3 = `24-48`, axis4 = `48-72`, y = n)) +
  geom_alluvium(aes(fill = `0-12`), alpha = 0.7, width = 1/12) +
  geom_stratum(width = 1/12, color = "grey30") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = block_labels, expand = c(.05, .05)) +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    x = "Hours from t0",
    y = "Patients",
    title = "Respiratory Support Transitions by DTW Cluster",
    subtitle = "Death treated as competing event (Expired or Hospice)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )

# ---- Save large, high resolution ----
ggsave(
  filename = file.path(trajectory_fig_dir, "sankey_respiratory_transitions_clusters.png"),
  plot     = p_sankey,
  width    = 18,     # inches
  height   = 14,     # inches
  dpi      = 600,    # journal-ready
  device   = "png",
  bg       = "white"
)



# ---- Optional: enforce consistent theme tweaks across panels ----
# (Only if you want to standardize axis/title sizes across p1-p9)
tweak_panel <- function(p) {
  p + theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    plot.margin = margin(4, 4, 4, 4)
  )
}

p_list <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9) %>% lapply(tweak_panel)

# ---- 3x3 grid ----
grid_3x3 <- cowplot::plot_grid(
  plotlist = p_list,
  ncol = 3,
  align = "hv"
)

# ---- Title (optional) ----
fig9 <- cowplot::ggdraw(grid_3x3) +
  cowplot::draw_label(
    "ARF Severity Trajectories by DTW Cluster (0–72h from t0)",
    x = 0.5, y = 0.995, hjust = 0.5, vjust = 1,
    fontface = "bold", size = 18
  )

print(fig9)

# ---- Save large + high-res ----
ggsave(file.path(trajectory_fig_dir, "traj_9panel_clusters.png"),
       plot = fig9,
       width = 20, height = 16, dpi = 600, bg = "white")

# Optional vector (best for manuscripts if all elements are ggplot)
ggsave(file.path(trajectory_fig_dir, "traj_9panel_clusters.pdf"),
       plot = fig9,
       width = 20, height = 16, device = cairo_pdf)




# ---- Load exposome data ----
no2_exp  <- read_csv(file.path(repo_root, "exposome", "no2_county_year.csv"), show_col_types = FALSE)
pm25_exp <- read_csv(file.path(repo_root, "exposome", "pm25_county_year.csv"), show_col_types = FALSE)

# Standardize names
no2_exp  <- no2_exp  %>% rename_with(tolower)
pm25_exp <- pm25_exp %>% rename_with(tolower)

# Expecting: county_code, year, no2 (or similar)
# Check column names:
# names(no2_exp); names(pm25_exp)

# Harmonize variable names.
no2_exp <- no2_exp %>%
  transmute(
    county_code = as.character(geoid),
    year = as.integer(year),
    no2 = as.numeric(no2_mean)
  )

pm25_exp <- pm25_exp %>%
  transmute(
    county_code = as.character(geoid),
    year = as.integer(year),
    pm25 = as.numeric(pm25_mean)
  )

# ---- Extract admission year from cohort ----
cohort_exp <- cohort_use %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admit_year = lubridate::year(admission_dttm),
    county_code = as.character(county_code)
  ) %>%
  left_join(cluster_df %>% mutate(hospitalization_id = as.character(hospitalization_id)),
            by = "hospitalization_id") %>%
  left_join(no2_exp,  by = c("county_code","admit_year" = "year")) %>%
  left_join(pm25_exp, by = c("county_code","admit_year" = "year"))

summary(cohort_exp$no2)
summary(cohort_exp$pm25)

iqr_no2  <- IQR(cohort_exp$no2,  na.rm = TRUE)
iqr_pm25 <- IQR(cohort_exp$pm25, na.rm = TRUE)

cohort_exp <- cohort_exp %>%
  mutate(
    no2_iqr  = no2  / iqr_no2,
    pm25_iqr = pm25 / iqr_pm25
  )

library(nnet)

cohort_exp$cluster <- factor(cohort_exp$cluster)

# Choose reference cluster (example: 1)
cohort_exp$cluster <- relevel(cohort_exp$cluster, ref = "1")

model_multinom <- multinom(
  cluster ~ no2_iqr + pm25_iqr +
    age_years + sex_category + race_category +
    factor(admit_year),
  data = cohort_exp
)

summary(model_multinom)

coefs <- summary(model_multinom)$coefficients
ses   <- summary(model_multinom)$standard.errors

zvals <- coefs / ses
pvals <- 2 * (1 - pnorm(abs(zvals)))

or_table <- exp(coefs)
ci_low   <- exp(coefs - 1.96*ses)
ci_high  <- exp(coefs + 1.96*ses)

results <- tibble(
  cluster = rep(rownames(coefs), each = ncol(coefs)),
  variable = rep(colnames(coefs), times = nrow(coefs)),
  OR = as.vector(or_table),
  CI_low = as.vector(ci_low),
  CI_high = as.vector(ci_high),
  p = as.vector(pvals)
)

results %>% filter(str_detect(variable,"no2|pm25"))

library(ggeffects)

pred_no2 <- ggpredict(model_multinom, terms = "pm25_iqr [all]")
plot(pred_no2)

pred_df <- as.data.frame(pred_no2)

ggplot(pred_df,
       aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  labs(
    x = "NO2 (IQR-scaled)",
    y = "Predicted Probability",
    color = "Cluster",
    title = "Predicted ARF Trajectory Cluster by NO2 Exposure"
  ) +
  theme_minimal(base_size = 14)










