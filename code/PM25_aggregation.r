
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(lubridate)
library(stringr)
library(glue)
library(fs)

# paths
base_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/Data/pm25"

# load counties shapefile (adjust path to your file)
counties_sf <- st_read("/Users/saborpete/Desktop/Peter/Postdoc/Data/conus county shapefiles/cb_2018_us_county_500k.shp")

# pick a county ID column
county_id_col <- intersect(c("GEOID","geoid","GEOIDFIPS","FIPS","GEOID10","COUNTYFP"), names(counties_sf))[1]
stopifnot(!is.na(county_id_col))

# --------------------------
# Helpers
# --------------------------

# Parse year & month from names like ...201801-201801.nc
parse_ym_from_fname <- function(x) {
  b  <- basename(x)
  m  <- str_match(b, "\\.(\\d{6})-(\\d{6})\\.nc$")[, 2]   # first YYYYMM
  if (is.na(m)) stop("Could not parse YYYYMM from filename: ", b)
  year  <- as.integer(substr(m, 1, 4))
  month <- as.integer(substr(m, 5, 6))
  list(year = year, month = month)
}

# Extract county means from a single NetCDF (uses first band by default)
extract_pm25_file <- function(nc_path, counties_sf) {
  r <- rast(nc_path)
  r <- r[[1]]  # first layer (rename here if you know var name)
  csf <- st_transform(counties_sf, crs(r))
  vals <- terra::extract(r, vect(csf), fun = mean, na.rm = TRUE)
  # vals: first col is ID, second col is raster values
  tibble(
    county_id = csf[[county_id_col]],
    pm25      = vals[[2]],
    file      = basename(nc_path)
  )
}

# --------------------------
# Index all .nc files 2018–2023
# --------------------------
year_dirs <- file.path(base_dir, as.character(2018:2023))
missing <- year_dirs[!dir_exists(year_dirs)]
if (length(missing)) message("These year folders are missing (will be skipped): ", paste(basename(missing), collapse = ", "))

nc_index <- map_dfr(year_dirs[dir_exists(year_dirs)], function(yd) {
  tibble(nc = dir_ls(yd, recurse = TRUE, type = "file", glob = "*.nc"))
}) %>%
  mutate(
    ym   = map(nc, parse_ym_from_fname),
    year = vapply(ym, `[[`, integer(1), "year"),
    month= vapply(ym, `[[`, integer(1), "month")
  ) %>%
  dplyr::select(-ym) %>%
  filter(year %in% 2018:2023)

stopifnot(nrow(nc_index) > 0)

# --------------------------
# Extract county means for each file
# --------------------------
pm25_df <- map_dfr(seq_len(nrow(nc_index)), function(i) {
  row <- nc_index[i, ]
  out <- extract_pm25_file(row$nc, counties_sf)
  out$year  <- row$year
  out$month <- row$month
  out
})

pm25_long <- pm25_df %>%
  dplyr::select(county_id, year, month, pm25)

# --------------------------
# Trend to 2024 (per county x month)
# --------------------------
trend_2024 <- pm25_long %>%
  filter(year >= 2018, year <= 2023) %>%
  group_by(county_id, month) %>%
  group_modify(~{
    dat <- .x
    if (sum(!is.na(dat$pm25)) >= 3) {
      fit <- lm(pm25 ~ year, data = dat)
      tibble(year = 2024, pm25 = predict(fit, newdata = tibble(year = 2024)))
    } else {
      tibble(year = 2024, pm25 = NA_real_)
    }
  }) %>%
  ungroup()

pm25_all <- bind_rows(pm25_long, trend_2024) %>%
  arrange(county_id, year, month)

# --------------------------
# Save
# --------------------------
out_dir <- file.path(base_dir, "outputs")
dir_create(out_dir)

write_csv(pm25_all,  file.path(out_dir, "county_pm25_monthly_2018_2024_long.csv"))

pm25_wide <- pm25_all %>%
  pivot_wider(id_cols = c(county_id, year),
              names_from = month, values_from = pm25,
              names_prefix = "m") %>%
  arrange(county_id, year)

write_csv(pm25_wide, file.path(out_dir, "county_pm25_monthly_2018_2024_wide.csv"))

message("Done. Wrote: ",
        file.path(out_dir, "county_pm25_monthly_2018_2024_long.csv"), " and _wide.csv")


library(tigris)
library(ggplot2)
options(tigris_use_cache = TRUE)

# -----------------------------
# Paths
# -----------------------------
data_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/Data/conus county level"
out_dir  <- file.path(data_dir, "outputs"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pm25_csv <- file.path(data_dir, "pm25_county_month.csv")
no2_csv  <- file.path(data_dir, "no2_county_month.csv")

# -----------------------------
# Read data (rename if needed)
# -----------------------------
pm25 <- read_csv(pm25_csv, show_col_types = FALSE) %>%
  rename(pm25 = !!(names(.)[tolower(names(.)) %in% c("pm25","pm2_5","pm_25")][1])) %>%
  mutate(county_id = sprintf("%05d", as.integer(county_id)))

no2 <- read_csv(no2_csv, show_col_types = FALSE) %>%
  rename(no2 = !!(names(.)[tolower(names(.)) %in% c("no2","nitrogen_dioxide")][1])) %>%
  mutate(county_id = sprintf("%05d", as.integer(county_fips)))

# -----------------------------
# Annual means per county
# -----------------------------
pm25_ann <- pm25 %>%
  group_by(county_id, year) %>%
  summarise(pm25 = mean(pm25, na.rm = TRUE), .groups = "drop")

no2_ann <- no2 %>%
  group_by(county_id, year) %>%
  summarise(no2 = mean(no2, na.rm = TRUE), .groups = "drop")

# -----------------------------
# Trend (slope) and change (Δ) helpers
# -----------------------------
trend_change <- function(df, value_col) {
  v <- rlang::ensym(value_col)
  df %>%
    arrange(year) %>%
    group_by(county_id) %>%
    group_modify(~{
      dat <- .x %>% filter(!is.na(!!v)) %>% arrange(year)
      if (nrow(dat) == 0) return(tibble(slope = NA_real_, delta = NA_real_, pct_delta = NA_real_, start_year = NA_integer_, end_year = NA_integer_))
      slope <- if (nrow(dat) >= 3) coef(lm(rlang::eval_tidy(v) ~ year, data = dat))[2] else NA_real_
      first_val <- dat[[rlang::as_string(v)]][1]
      last_val  <- dat[[rlang::as_string(v)]][nrow(dat)]
      delta <- last_val - first_val
      pct_delta <- if (!is.na(first_val) && first_val != 0) 100 * delta / first_val else NA_real_
      tibble(
        slope = slope,
        delta = delta,
        pct_delta = pct_delta,
        start_year = dat$year[1],
        end_year = dat$year[nrow(dat)]
      )
    }) %>%
    ungroup()
}

pm25_tr <- trend_change(pm25_ann, pm25) %>% mutate(pollutant = "PM2.5")
no2_tr  <- trend_change(no2_ann,  no2 ) %>% mutate(pollutant = "NO₂")

trends <- bind_rows(pm25_tr, no2_tr)

# -----------------------------
# Geometry (CONUS counties)
# -----------------------------
counties <- counties(cb = TRUE, year = 2020, class = "sf") %>%
  filter(!STATEFP %in% c("02","15","72","78","69","66","60")) %>%  # drop AK, HI, PR & territories
  st_transform(5070) %>% # NAD83 / Conus Albers
  dplyr::select(GEOID, NAME, STATEFP)

# Join trends to geometry
geo_trends <- counties %>%
  left_join(trends, by = c("GEOID" = "county_id"))

# -----------------------------
# Plot helpers
# -----------------------------
theme_map <- theme_void(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# 1) Slope maps (β per year)
p_slope <- ggplot(geo_trends) +
  geom_sf(aes(fill = slope), color = NA) +
  scale_fill_gradient2(name = "Slope (per year)", na.value = "grey90") +
  facet_wrap(~ pollutant, ncol = 2) +
  labs(title = "County PM2.5 and NO₂ Trend Slopes (β/year)") +
  theme_map
ggsave(file.path(out_dir, "trend_slopes_pm25_no2_conus.png"), p_slope, width = 12, height = 7, dpi = 300)

# 2) Absolute change (last – first)
p_delta <- ggplot(geo_trends) +
  geom_sf(aes(fill = delta), color = NA) +
  scale_fill_gradient2(name = "Δ (last – first)", na.value = "grey90") +
  facet_wrap(~ pollutant, ncol = 2) +
  labs(title = "County PM2.5 and NO₂ Absolute Change (last – first year)") +
  theme_map
ggsave(file.path(out_dir, "absolute_change_pm25_no2_conus.png"), p_delta, width = 12, height = 7, dpi = 300)

# 3) Percent change
p_pct <- ggplot(geo_trends) +
  geom_sf(aes(fill = pct_delta), color = NA) +
  scale_fill_gradient2(name = "% Change", na.value = "grey90", labels = scales::label_percent(accuracy = 1, scale = 1)) +
  facet_wrap(~ pollutant, ncol = 2) +
  labs(title = "County PM2.5 and NO₂ Percent Change (last – first year)") +
  theme_map
ggsave(file.path(out_dir, "percent_change_pm25_no2_conus.png"), p_pct, width = 12, height = 7, dpi = 300)

# Also save as a GeoPackage for GIS use
st_write(geo_trends, file.path(out_dir, "county_trends_pm25_no2.gpkg"), delete_dsn = TRUE, quiet = TRUE)


library(scales)


# Helper: symmetric limits around 0 using quantile trimming
sym_lims <- function(x, q = 0.99) {
  if (all(is.na(x))) return(c(NA, NA))
  rng <- quantile(x, probs = c(1 - q, q), na.rm = TRUE)
  m <- max(abs(rng))
  c(-m, m)
}

# Aesthetic settings: blue (neg) → white → red (pos), long legend
grad2 <- function(name, limits, label_fun = waiver()) {
  scale_fill_gradient2(
    name   = name,
    low    = "#2B6CB0",  # blue = decrease
    mid    = "white",
    high   = "#C53030",  # red  = increase
    limits = limits,
    oob    = squish,
    labels = label_fun
  )
}

theme_map <- theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

# Split by pollutant
geo_pm25 <- geo_trends %>% filter(pollutant == "PM2.5")
geo_no2  <- geo_trends %>% filter(pollutant == "NO₂")

long_bar <- guides(
  fill = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    barwidth       = unit(20, "cm"),  # ← longer legend bar
    barheight      = unit(0.6, "cm"),
    ticks          = TRUE
  )
)
# Units
pm25_unit <- "µg/m³"
no2_unit  <- "ppb"
# ---------- SLOPE (unit per year) ----------
lims_pm25_slope <- sym_lims(geo_pm25$slope, q = 0.995)
lims_no2_slope  <- sym_lims(geo_no2$slope,  q = 0.995)

p_pm25_slope <- ggplot(geo_pm25) +
  geom_sf(aes(fill = slope), color = NA) +
  grad2(paste0("Slope (", pm25_unit, " per year)"), limits = lims_pm25_slope) +
  labs(title = "PM2.5 Trend Slope (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "PM25_trend_slope_CONUS.png"),
       p_pm25_slope, width = 12, height = 7, dpi = 300)

p_no2_slope <- ggplot(geo_no2) +
  geom_sf(aes(fill = slope), color = NA) +
  grad2(paste0("Slope (", no2_unit, " per year)"), limits = lims_no2_slope) +
  labs(title = "NO₂ Trend Slope (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "NO2_trend_slope_CONUS.png"),
       p_no2_slope, width = 12, height = 7, dpi = 300)

# ---------- ABSOLUTE CHANGE (last – first) ----------
lims_pm25_delta <- sym_lims(geo_pm25$delta, q = 0.995)
lims_no2_delta  <- sym_lims(geo_no2$delta,  q = 0.995)

p_pm25_delta <- ggplot(geo_pm25) +
  geom_sf(aes(fill = delta), color = NA) +
  grad2(paste0("Δ (", pm25_unit, ")"), limits = lims_pm25_delta) +
  labs(title = "PM2.5 Absolute Change (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "PM25_absolute_change_CONUS.png"),
       p_pm25_delta, width = 12, height = 7, dpi = 300)

p_no2_delta <- ggplot(geo_no2) +
  geom_sf(aes(fill = delta), color = NA) +
  grad2(paste0("Δ (", no2_unit, ")"), limits = lims_no2_delta) +
  labs(title = "NO₂ Absolute Change (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "NO2_absolute_change_CONUS.png"),
       p_no2_delta, width = 12, height = 7, dpi = 300)

# ---------- PERCENT CHANGE ----------
lims_pm25_pct <- sym_lims(geo_pm25$pct_delta, q = 0.99)
lims_no2_pct  <- sym_lims(geo_no2$pct_delta,  q = 0.99)

p_pm25_pct <- ggplot(geo_pm25) +
  geom_sf(aes(fill = pct_delta), color = NA) +
  grad2("% Change", limits = lims_pm25_pct, label_fun = scales::label_percent(accuracy = 1, scale = 1)) +
  labs(title = "PM2.5 Percent Change (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "PM25_percent_change_CONUS.png"),
       p_pm25_pct, width = 12, height = 7, dpi = 300)

p_no2_pct <- ggplot(geo_no2) +
  geom_sf(aes(fill = pct_delta), color = NA) +
  grad2("% Change", limits = lims_no2_pct, label_fun = scales::label_percent(accuracy = 1, scale = 1)) +
  labs(title = "NO₂ Percent Change (2018-2024)") +
  theme_map + long_bar
ggsave(file.path(out_dir, "NO2_percent_change_CONUS.png"),
       p_no2_pct, width = 12, height = 7, dpi = 300)











