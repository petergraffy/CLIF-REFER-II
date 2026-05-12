# ================================================================================================
# ICU Respiratory Failure Environmental Risk (REFER) | PI: Peter Graffy
# Patient Profile Builder (CLIF)
#
# Purpose:
#   Add site-local patient characterization for the ARF cohort, including Charlson comorbidity
#   components/scores and first-24-hour SOFA severity from both ICU admission and ARF onset.
#
# Outputs (written to output/run_[SITE]_[DATE]/ when available):
#   patient_profiles_arf72_[SITE]_[DATE].csv
#   patient_profile_numeric_summary_[SITE]_[DATE].csv
#   patient_profile_categorical_summary_[SITE]_[DATE].csv
#   patient_profile_cluster_numeric_summary_[SITE]_[DATE].csv       when clusters are available
#   patient_profile_cluster_categorical_summary_[SITE]_[DATE].csv  when clusters are available
#
# Note:
#   The patient-level profile file is for site-local analysis and validation. Federated exports should
#   use the aggregate summary files unless patient-level sharing is explicitly approved.
# ================================================================================================

pkgs <- c("tidyverse", "lubridate", "arrow", "fst", "data.table", "comorbidity")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(arrow)
  library(fst)
  library(data.table)
  library(comorbidity)
})

source("utils/config.R")
source("utils/sofa_calculator.R")
stopifnot(exists("config"))

repo        <- config$repo
site_name   <- config$site_name
tables_path <- config$tables_path
file_type   <- config$file_type

stopifnot(!is.null(repo), nzchar(repo))
stopifnot(!is.null(tables_path), nzchar(tables_path))

SITE_NAME <- if (!is.null(site_name) && nzchar(site_name)) site_name else "SITE"
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")

safe_ts <- function(x) {
  if (inherits(x, "POSIXt")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  suppressWarnings(as.POSIXct(x, tz = "UTC"))
}

safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

make_name <- function(stem) paste0(stem, "_", SITE_NAME, "_", SYSTEM_DATE, ".csv")

find_latest_run_dir <- function(repo_root) {
  run_dirs <- list.dirs(file.path(repo_root, "output"), recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("run_", basename(run_dirs), fixed = TRUE)]
  if (!length(run_dirs)) return(file.path(repo_root, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE)))
  run_dirs[which.max(file.info(run_dirs)$mtime)]
}

find_latest_output <- function(repo_root, pattern) {
  out_files <- list.files(file.path(repo_root, "output"), pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(out_files)) return(NA_character_)
  out_files[which.max(file.info(out_files)$mtime)]
}

add_missing_cols <- function(df, cols) {
  missing_cols <- setdiff(cols, names(df))
  for (missing_col in missing_cols) df[[missing_col]] <- NA
  df
}

select_cols <- function(df, cols) {
  df <- df %>% rename_with(tolower)
  df <- add_missing_cols(df, tolower(cols))
  df %>% dplyr::select(all_of(tolower(cols)))
}

load_clif_tables <- function(tables_path, file_type) {
  tables_path <- normalizePath(tables_path, mustWork = TRUE)
  exts <- strsplit(file_type, "[/|,; ]+")[[1]]
  exts <- exts[nzchar(exts)]
  if (length(exts) == 0) exts <- c("csv", "parquet", "fst")
  ext_pat <- paste0("\\.(", paste(unique(exts), collapse = "|"), ")$")

  all_files <- list.files(
    path = tables_path,
    pattern = ext_pat,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (!length(all_files)) stop("No CLIF files found under: ", tables_path)

  bn <- basename(all_files)
  looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)
  base_no_ext <- tools::file_path_sans_ext(tolower(bn))
  base_norm <- ifelse(looks_clif, base_no_ext, paste0("clif_", base_no_ext))
  found_map <- stats::setNames(all_files, base_norm)

  desired <- paste0("clif_", c(
    "patient", "hospitalization", "adt", "vitals", "labs", "respiratory_support",
    "hospital_diagnosis", "medication_admin_continuous", "patient_assessments"
  ))
  available <- intersect(desired, names(found_map))
  out <- lapply(found_map[available], read_any)
  names(out) <- available
  out
}

load_cohort <- function(repo_root) {
  if (exists("cohort_arf", envir = .GlobalEnv)) {
    return(get("cohort_arf", envir = .GlobalEnv))
  }

  cohort_path <- find_latest_output(repo_root, "^cohort_arf72_.*\\.csv$")
  if (is.na(cohort_path)) {
    stop("Could not find cohort_arf in memory or a saved cohort_arf72 CSV under output/. Run code/01_cohort_identification.R first.")
  }
  readr::read_csv(cohort_path, show_col_types = FALSE)
}

load_clusters <- function(repo_root) {
  if (exists("cluster_df", envir = .GlobalEnv)) {
    return(get("cluster_df", envir = .GlobalEnv) %>% mutate(hospitalization_id = as.character(hospitalization_id)))
  }
  cluster_path <- find_latest_output(repo_root, "^trajectory_cluster_assignments.*\\.csv$")
  if (is.na(cluster_path)) return(NULL)
  readr::read_csv(cluster_path, show_col_types = FALSE) %>%
    mutate(hospitalization_id = as.character(hospitalization_id))
}

compute_charlson <- function(hospital_dx, cohort_ids) {
  if (is.null(hospital_dx) || !nrow(hospital_dx)) {
    return(tibble(hospitalization_id = cohort_ids, charlson_score = NA_real_, charlson_dx_available = FALSE))
  }

  dx <- hospital_dx %>%
    rename_with(tolower) %>%
    add_missing_cols(c("hospitalization_id", "diagnosis_code", "diagnosis_code_format",
                       "diagnosis_category", "diagnosis_type", "poa_present")) %>%
    mutate(
      hospitalization_id = as.character(hospitalization_id),
      diagnosis_code = as.character(diagnosis_code),
      diagnosis_code_format = toupper(coalesce(as.character(diagnosis_code_format), "")),
      diagnosis_category = toupper(coalesce(as.character(diagnosis_category), "")),
      diagnosis_type = toupper(coalesce(as.character(diagnosis_type), "")),
      poa_present = suppressWarnings(as.integer(poa_present)),
      code_clean = gsub("[^A-Za-z0-9]", "", diagnosis_code)
    ) %>%
    filter(hospitalization_id %in% cohort_ids, !is.na(code_clean), nzchar(code_clean))

  if (!nrow(dx)) {
    return(tibble(hospitalization_id = cohort_ids, charlson_score = NA_real_, charlson_dx_available = FALSE))
  }

  poa_is_informative <- any(!is.na(dx$poa_present))
  if (poa_is_informative) dx <- dx %>% filter(poa_present == 1)

  dx <- dx %>%
    mutate(
      icd_version = case_when(
        str_detect(diagnosis_code_format, "10") | str_detect(diagnosis_category, "ICD.?10") ~ "icd10",
        str_detect(diagnosis_code_format, "9") | str_detect(diagnosis_category, "ICD.?9") ~ "icd9",
        str_detect(code_clean, "^[A-Za-z]") ~ "icd10",
        str_detect(code_clean, "^[0-9]") ~ "icd9",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(icd_version))

  score_one_map <- function(dx_subset, map_name) {
    if (!nrow(dx_subset)) return(tibble())
    tryCatch(
      comorbidity::comorbidity(
        x = dx_subset %>% dplyr::select(id = hospitalization_id, code = code_clean),
        id = "id",
        code = "code",
        map = map_name,
        assign0 = FALSE
      ),
      error = function(e) {
        warning("Charlson map failed for ", map_name, ": ", conditionMessage(e))
        tibble()
      }
    )
  }

  charlson_components <- bind_rows(
    score_one_map(dx %>% filter(icd_version == "icd10"), "charlson_icd10_quan"),
    score_one_map(dx %>% filter(icd_version == "icd9"), "charlson_icd9_quan")
  )

  if (!nrow(charlson_components)) {
    return(tibble(hospitalization_id = cohort_ids, charlson_score = NA_real_, charlson_dx_available = nrow(dx) > 0))
  }

  component_cols <- setdiff(names(charlson_components), "id")
  charlson_all <- charlson_components %>%
    mutate(id = as.character(id)) %>%
    group_by(id) %>%
    summarise(across(all_of(component_cols), ~ max(.x, na.rm = TRUE)), .groups = "drop")

  charlson_score_input <- charlson_all %>% dplyr::select(all_of(component_cols))
  class(charlson_score_input) <- c("comorbidity", class(charlson_score_input))
  attr(charlson_score_input, "map") <- "charlson_icd10_quan"

  charlson_all <- charlson_all %>%
    mutate(
      charlson_score = comorbidity::score(charlson_score_input, weights = "charlson", assign0 = FALSE),
      charlson_dx_available = TRUE,
      charlson_poa_filter = if_else(poa_is_informative, "present_on_admission_only", "all_diagnoses_poa_unavailable")
    ) %>%
    rename(hospitalization_id = id)

  tibble(hospitalization_id = cohort_ids) %>%
    left_join(charlson_all, by = "hospitalization_id") %>%
    mutate(
      charlson_score = coalesce(charlson_score, 0),
      charlson_dx_available = coalesce(charlson_dx_available, FALSE),
      charlson_poa_filter = coalesce(charlson_poa_filter, if_else(poa_is_informative, "present_on_admission_only", "all_diagnoses_poa_unavailable"))
    )
}

prep_sofa_inputs <- function(clif_tables) {
  vitals_df <- clif_tables[["clif_vitals"]] %>%
    select_cols(c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value"))

  labs_df <- clif_tables[["clif_labs"]] %>%
    select_cols(c("hospitalization_id", "lab_result_dttm", "lab_category", "lab_value", "lab_value_numeric")) %>%
    mutate(lab_value_numeric = coalesce(suppressWarnings(as.numeric(lab_value_numeric)),
                                        suppressWarnings(as.numeric(lab_value))))

  support_df <- clif_tables[["clif_respiratory_support"]] %>%
    select_cols(c("hospitalization_id", "recorded_dttm", "device_category", "fio2_set"))

  med_admin_df <- clif_tables[["clif_medication_admin_continuous"]]
  if (is.null(med_admin_df)) med_admin_df <- tibble()
  med_admin_df <- med_admin_df %>%
    select_cols(c("hospitalization_id", "admin_dttm", "med_name", "med_category", "med_dose", "med_dose_unit"))

  scores_df <- clif_tables[["clif_patient_assessments"]]
  if (is.null(scores_df)) scores_df <- tibble()
  scores_df <- scores_df %>%
    select_cols(c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value"))

  list(
    vitals_df = vitals_df,
    labs_df = labs_df,
    support_df = support_df,
    med_admin_df = med_admin_df,
    scores_df = scores_df
  )
}

calculate_window_sofa <- function(cohort_data, sofa_inputs, window_start_col, prefix) {
  sofa_cohort <- cohort_data %>%
    mutate(
      hospitalization_id = as.character(hospitalization_id),
      icu_admit_time = safe_ts(.data[[window_start_col]])
    ) %>%
    filter(!is.na(icu_admit_time)) %>%
    dplyr::select(hospitalization_id, icu_admit_time)

  if (!nrow(sofa_cohort)) {
    return(tibble(hospitalization_id = character()))
  }

  out <- calculate_sofa(
    cohort_data = sofa_cohort,
    vitals_df = sofa_inputs$vitals_df,
    labs_df = sofa_inputs$labs_df,
    support_df = sofa_inputs$support_df,
    med_admin_df = sofa_inputs$med_admin_df,
    scores_df = sofa_inputs$scores_df,
    window_hours = 24,
    safe_ts = safe_ts
  )

  rename_with(out, ~ paste0(prefix, "_", .x), -hospitalization_id)
}

numeric_summary <- function(df, vars, strata = NULL) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(tibble())
  base <- if (is.null(strata)) df %>% mutate(.stratum = "overall") else df %>% mutate(.stratum = as.character(.data[[strata]]))

  base %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
    group_by(.stratum, variable) %>%
    summarise(
      n = n(),
      n_observed = sum(!is.na(value)),
      missing_pct = mean(is.na(value)) * 100,
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    mutate(across(c(mean, sd, median, q25, q75), ~ ifelse(is.nan(.x), NA_real_, .x))) %>%
    rename(stratum = .stratum)
}

categorical_summary <- function(df, vars, strata = NULL) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(tibble())
  base <- if (is.null(strata)) df %>% mutate(.stratum = "overall") else df %>% mutate(.stratum = as.character(.data[[strata]]))

  base %>%
    pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "level") %>%
    mutate(level = coalesce(as.character(level), "MISSING")) %>%
    group_by(.stratum, variable, level) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    rename(stratum = .stratum)
}

cat("Loading CLIF tables for patient profiling...\n")
if (exists("clif_tables", envir = .GlobalEnv)) {
  clif_tables <- get("clif_tables", envir = .GlobalEnv)
  disk_tables <- load_clif_tables(tables_path, file_type)
  for (nm in setdiff(names(disk_tables), names(clif_tables))) {
    clif_tables[[nm]] <- disk_tables[[nm]]
  }
} else {
  clif_tables <- load_clif_tables(tables_path, file_type)
}

required_tables <- paste0("clif_", c("patient", "hospitalization", "adt", "vitals", "labs", "respiratory_support"))
missing_required <- setdiff(required_tables, names(clif_tables))
if (length(missing_required)) stop("Missing required profiling tables: ", paste(missing_required, collapse = ", "))

out_dir <- if (exists("out_dir", envir = .GlobalEnv)) get("out_dir", envir = .GlobalEnv) else find_latest_run_dir(repo)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort <- load_cohort(repo) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    patient_id = as.character(patient_id),
    t0 = safe_ts(t0),
    first_icu_in = safe_ts(first_icu_in)
  )

cohort_ids <- unique(cohort$hospitalization_id)

patient <- clif_tables[["clif_patient"]] %>%
  select_cols(c("patient_id", "birth_date", "sex_category", "race_category", "ethnicity_category", "preferred_language", "death_dttm")) %>%
  mutate(patient_id = as.character(patient_id), birth_date = safe_date(birth_date), death_dttm = safe_ts(death_dttm))

hospitalization <- clif_tables[["clif_hospitalization"]] %>%
  select_cols(c("patient_id", "hospitalization_id", "admission_dttm", "discharge_dttm", "age_at_admission",
                "discharge_name", "discharge_category")) %>%
  mutate(
    patient_id = as.character(patient_id),
    hospitalization_id = as.character(hospitalization_id),
    admission_dttm = safe_ts(admission_dttm),
    discharge_dttm = safe_ts(discharge_dttm)
  )

charlson <- compute_charlson(clif_tables[["clif_hospital_diagnosis"]], cohort_ids)
sofa_inputs <- prep_sofa_inputs(clif_tables)

cat("Calculating ICU-admission SOFA...\n")
sofa_icu24 <- calculate_window_sofa(cohort, sofa_inputs, "first_icu_in", "sofa_icu24")

cat("Calculating ARF-onset SOFA...\n")
sofa_arf24 <- calculate_window_sofa(cohort, sofa_inputs, "t0", "sofa_arf24")

clusters <- load_clusters(repo)

patient_profiles <- cohort %>%
  left_join(patient, by = "patient_id") %>%
  left_join(hospitalization, by = c("hospitalization_id", "patient_id"), suffix = c("", "_hosp")) %>%
  left_join(charlson, by = "hospitalization_id") %>%
  left_join(sofa_icu24, by = "hospitalization_id") %>%
  left_join(sofa_arf24, by = "hospitalization_id")

if (!is.null(clusters)) {
  patient_profiles <- patient_profiles %>%
    left_join(clusters %>% dplyr::select(hospitalization_id, cluster), by = "hospitalization_id")
}

numeric_vars <- c(
  "age_at_admission", "charlson_score",
  "sofa_icu24_sofa_total", "sofa_icu24_sofa_cv", "sofa_icu24_sofa_coag", "sofa_icu24_sofa_liver",
  "sofa_icu24_sofa_renal", "sofa_icu24_sofa_resp", "sofa_icu24_sofa_cns",
  "sofa_arf24_sofa_total", "sofa_arf24_sofa_cv", "sofa_arf24_sofa_coag", "sofa_arf24_sofa_liver",
  "sofa_arf24_sofa_renal", "sofa_arf24_sofa_resp", "sofa_arf24_sofa_cns"
)

charlson_components <- setdiff(
  names(patient_profiles)[str_detect(names(patient_profiles), "^(mi|chf|pvd|cevd|dementia|cpd|rheumd|pud|mld|diab|diabwc|hp|rend|canc|msld|metacanc|aids)$")],
  character()
)

categorical_vars <- c(
  "sex_category", "race_category", "ethnicity_category", "preferred_language",
  "discharge_category", "discharge_name", "charlson_dx_available", "charlson_poa_filter",
  charlson_components
)

profile_numeric_summary <- numeric_summary(patient_profiles, numeric_vars)
profile_categorical_summary <- categorical_summary(patient_profiles, categorical_vars)

write_csv(patient_profiles, file.path(out_dir, make_name("patient_profiles_arf72")))
write_csv(profile_numeric_summary, file.path(out_dir, make_name("patient_profile_numeric_summary")))
write_csv(profile_categorical_summary, file.path(out_dir, make_name("patient_profile_categorical_summary")))

if ("cluster" %in% names(patient_profiles)) {
  cluster_numeric_summary <- patient_profiles %>%
    filter(!is.na(cluster)) %>%
    numeric_summary(numeric_vars, strata = "cluster")
  cluster_categorical_summary <- patient_profiles %>%
    filter(!is.na(cluster)) %>%
    categorical_summary(categorical_vars, strata = "cluster")

  write_csv(cluster_numeric_summary, file.path(out_dir, make_name("patient_profile_cluster_numeric_summary")))
  write_csv(cluster_categorical_summary, file.path(out_dir, make_name("patient_profile_cluster_categorical_summary")))
}

cat("\nSaved patient profile outputs to: ", out_dir, "\n", sep = "")
cat("  - ", make_name("patient_profiles_arf72"), "\n", sep = "")
cat("  - ", make_name("patient_profile_numeric_summary"), "\n", sep = "")
cat("  - ", make_name("patient_profile_categorical_summary"), "\n", sep = "")
if ("cluster" %in% names(patient_profiles)) {
  cat("  - ", make_name("patient_profile_cluster_numeric_summary"), "\n", sep = "")
  cat("  - ", make_name("patient_profile_cluster_categorical_summary"), "\n", sep = "")
}
