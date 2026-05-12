## Code Directory

This directory now contains the canonical site-facing workflow for the CLIF REFER-II acute respiratory failure trajectory analysis.

### Canonical Workflow

0. Restore the R environment.

```r
source("code/00_renv_restore.R")
```

1. Build the trajectory cohorts.

```r
source("code/01_cohort_identification.R")
```

This script identifies adult ICU hospitalizations from 2018-2024 and creates an ARF-wide primary cohort plus respiratory-support sensitivity cohorts:

- `cohort_arf72`: first evidence of ARF during the ICU stay; `t0` is earliest qualifying ARF evidence.
- `cohort_primary_imv72`: invasive mechanical ventilation sensitivity cohort; `t0` is IMV start.
- `cohort_secondary_adv72`: advanced respiratory support sensitivity cohort; `t0` is first IMV, NIPPV, CPAP, or high-flow nasal cannula.

ARF evidence sources include advanced respiratory support, SpO2 <90%, PaO2 <=60 on room air, P/F ratio <=300, or hypercapnic acidosis with PaCO2 >=45 and arterial pH <7.35.

It writes cohort, exclusion, and flow files under `output/run_[SITE]_[DATE]/`.

2. Build and cluster ARF trajectories.

```r
source("code/02_trajectories.R")
```

This script currently expects the cohort objects from step 1 to remain in the R session. It defaults to `cohort_arf72` when available. It builds 0-72 hour respiratory trajectories from `t0`, uses CLIF/mCIDE-aligned respiratory support, vital sign, ABG, and biomarker features, encodes active care, missing-active-care, discharge, and death as explicit trajectory state channels, applies missingness checks and state-aware imputation for DTW compatibility, clusters trajectories with dynamic time warping, and creates trajectory and respiratory support transition figures.

### Exposome Utilities

- `PM25_aggregation.r`: aggregates PM2.5 inputs to county-year and county-month panels.
- `aggregate_no2.r`: aggregates NO2 inputs to county-year and county-month panels.

These scripts support maintenance of the included `exposome/` files and are not required at every participating CLIF site if the shared exposome files are already distributed with the project.

### Archived Exploratory Code

Older exploratory scripts were moved to `code/archive/` so the main workflow is easier to follow. They are preserved for provenance, but they should not be treated as the current analysis entry points.

### Near-Term Refactor Target

For federated execution, `02_trajectories.R` should be split into:

- `02_build_trajectory_features.R`: builds and saves site-level hourly trajectory features.
- `03_cluster_and_profile_trajectories.R`: fits or applies DTW phenotypes, summarizes clusters, and creates exportable site-level results.
- `04_federated_models.R`: runs site-level association and prognostic models with pre-specified output tables.
