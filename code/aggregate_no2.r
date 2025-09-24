
library(ncdf4)
library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(purrr)


# Point to your folder
nc_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/Data"

# Grab all monthly NO2 files 2019–2024
nc_files <- list.files(nc_dir,
                       pattern = "monthly_mean_tropomi_lur_conus_surface_no2_.*\\.nc$",
                       full.names = TRUE)

# Read in your county shapefile (FIPS codes included)
# Replace with your actual shapefile/crosswalk path
county_shp <- st_read("/Users/saborpete/Desktop/Peter/Postdoc/Data/cb_2018_us_county_500k.shp") %>%
  st_transform(4326)  # match raster CRS

process_no2_file <- function(nc_file, counties) {
  r <- raster(nc_file)  # read NO2 raster
  r <- projectRaster(r, crs = st_crs(counties)$proj4string) # align CRS
  
  # Extract mean per county
  vals <- exactextractr::exact_extract(r, counties, 'mean')
  
  # Parse year-month from filename
  fname <- basename(nc_file)
  ym <- str_match(fname, "_([0-9]{2})_([A-Za-z]+)_([0-9]{4})")[, -1]
  month_num <- match(ym[2], month.name) # convert month name to number
  year_num <- as.integer(ym[3])
  
  tibble(
    county_fips = counties$GEOID,
    year = year_num,
    month = month_num,
    no2 = vals
  )
}

no2_monthly <- map_dfr(nc_files, process_no2_file, counties = county_shp)

# 1) Clean NaNs and impossible values
no2_monthly_clean <- no2_monthly %>%
  mutate(
    no2  = ifelse(is.nan(no2), NA_real_, no2),
    no2  = ifelse(is.infinite(no2), NA_real_, no2),
    no2  = ifelse(no2 < 0, NA_real_, no2),  # clamp negatives
    date = make_date(year, month, 1)
  ) %>%
  arrange(county_fips, date)

seasonal_means_county <- no2_monthly_clean %>%
  group_by(county_fips, month) %>%
  summarise(no2_c_mo_mean = mean(no2, na.rm = TRUE), .groups = "drop")

overall_mean_county <- no2_monthly_clean %>%
  group_by(county_fips) %>%
  summarise(no2_c_overall = mean(no2, na.rm = TRUE), .groups = "drop")

seasonal_means_national <- no2_monthly_clean %>%
  group_by(month) %>%
  summarise(no2_nat_mo_mean = mean(no2, na.rm = TRUE), .groups = "drop")

overall_mean_national <- tibble(no2_nat_overall =
                                  mean(no2_monthly_clean$no2, na.rm = TRUE))

no2_monthly_filled <- no2_monthly_clean %>%
  left_join(seasonal_means_county,   by = c("county_fips","month")) %>%
  left_join(overall_mean_county,     by = "county_fips") %>%
  left_join(seasonal_means_national, by = "month") %>%
  mutate(
    no2_imputed = case_when(
      !is.na(no2)                 ~ no2,
      is.na(no2) & !is.na(no2_c_mo_mean)   ~ no2_c_mo_mean,
      is.na(no2) &  is.na(no2_c_mo_mean) & !is.na(no2_c_overall) ~ no2_c_overall,
      is.na(no2) &  is.na(no2_c_mo_mean) &  is.na(no2_c_overall) & !is.na(no2_nat_mo_mean) ~ no2_nat_mo_mean,
      TRUE                        ~ overall_mean_national$no2_nat_overall
    )
  ) %>%
  dplyr::select(county_fips, year, month, date, no2 = no2_imputed) %>%
  arrange(county_fips, date)

# 3) Fit per-county linear trend (use only observed/imputed 2019–2024)
#    Guardrails: require >= 18 months of non-missing to estimate a slope.
no2_for_slope <- no2_monthly_filled %>%
  filter(date >= ymd("2019-01-01") & date <= ymd("2024-12-01")) %>%
  group_by(county_fips) %>%
  mutate(time_idx = as.numeric(difftime(date, min(date), units = "days"))) %>%
  ungroup()

slopes <- no2_for_slope %>%
  group_by(county_fips) %>%
  group_modify(~{
    dat <- .x %>% filter(!is.na(no2))
    if (nrow(dat) >= 18) {
      fit <- lm(no2 ~ time_idx, data = dat)
      tibble(
        have_slope = TRUE,
        intercept  = unname(coef(fit)[1]),
        slope      = unname(coef(fit)[2]),
        t0         = min(dat$date)  # origin for time_idx
      )
    } else {
      tibble(have_slope = FALSE, intercept = NA_real_, slope = NA_real_, t0 = min(.x$date))
    }
  }) %>%
  ungroup()

# 4) Build 2018 months and back-fill using slope if available; otherwise seasonal mean
dates_2018 <- tibble(date = seq(ymd("2018-01-01"), ymd("2018-12-01"), by = "1 month")) %>%
  mutate(year = year(date), month = month(date))

# Precompute seasonal means from (filled) 2019–2024 for fallback
seasonal_2019_2024 <- no2_monthly_filled %>%
  filter(year >= 2019 & year <= 2024) %>%
  group_by(county_fips, month) %>%
  summarise(no2_c_mo_mean_19_24 = mean(no2, na.rm = TRUE), .groups = "drop")

# Counties present in your 2019–2024 panel
all_counties <- no2_monthly_filled %>%
  distinct(county_fips)

# Build county × 2018-month grid, then attach slopes and seasonal fallbacks
no2_2018 <- crossing(all_counties, dates_2018) %>%   # cartesian product, no warnings
  rename(date = date, month = month, year = year) %>%
  # bring slope params (may be NA if county had <18 months)
  left_join(slopes %>% dplyr::select(county_fips, have_slope, intercept, slope, t0),
            by = "county_fips") %>%
  # seasonal fallback from 2019–2024 (county × month mean)
  left_join(seasonal_2019_2024, by = c("county_fips","month")) %>%
  mutate(
    # time since slope origin (per county)
    time_idx_2018 = as.numeric(difftime(date, t0, units = "days")),
    no2_trend     = intercept + slope * time_idx_2018,
    # choose value: trend if slope exists & finite, else seasonal mean, else national mean
    no2 = dplyr::case_when(
      isTRUE(have_slope) & is.finite(no2_trend)               ~ no2_trend,
      !is.na(no2_c_mo_mean_19_24)                             ~ no2_c_mo_mean_19_24,
      TRUE                                                    ~ overall_mean_national$no2_nat_overall
    ),
    # final clamps
    no2 = ifelse(is.na(no2), overall_mean_national$no2_nat_overall, no2),
    no2 = pmax(no2, 0)
  ) %>%
  dplyr::select(county_fips, year, month, date, no2) %>%
  arrange(county_fips, date)

# 5) Bind 2018 + cleaned/imputed 2019–2024, final sanity clamps
no2_full <- bind_rows(no2_monthly_filled, no2_2018) %>%
  arrange(county_fips, date) %>%
  mutate(
    no2 = ifelse(no2 < 0, 0, no2),           # non-negative
    no2 = ifelse(is.infinite(no2), NA_real_, no2)
  )

# Optional: report counties that still have NAs (should be none)
leftovers <- no2_full %>% filter(is.na(no2)) %>% distinct(county_fips)
if (nrow(leftovers) > 0) {
  message("Counties with remaining NA after imputation: ", paste(leftovers$county_fips, collapse = ", "))
}


write.csv(no2_full, "county_monthly_no2_2018_2024_clean.csv", row.names = FALSE)








