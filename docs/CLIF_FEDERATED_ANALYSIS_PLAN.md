# CLIF Federated Analysis Plan

## Project Frame

REFER-II will phenotype acute respiratory failure trajectories among adult ICU patients using CLIF respiratory support, vital sign, laboratory, and hospitalization tables. The primary trajectory window is 0-72 hours after respiratory support initiation. Dynamic time warping will be used to identify trajectory phenotypes that are clinically interpretable, prognostically useful, and suitable for federated evaluation across CLIF sites.

The environmental aim will evaluate whether chronic or recent ambient exposures are associated with trajectory phenotype membership and with outcomes conditional on phenotype.

## Federated Principle

Patient-level clinical data stay at each site. Sites run the same versioned scripts locally and export only approved aggregate products. Central pooling uses site-level counts, model coefficients, variance estimates, quality control summaries, and non-identifying trajectory prototypes.

## Cohorts

### Primary Cohort

Adult ICU hospitalizations from 2018-2024 with invasive mechanical ventilation starting within 24 hours of first ICU admission.

- Index time: first IMV record, `t0`.
- Trajectory window: `t0` through `t0 + 72h`.
- Rationale: high specificity for severe ARF and dense respiratory support measurements.

### Secondary Cohort

Adult ICU hospitalizations from 2018-2024 with advanced respiratory support starting within 24 hours of first ICU admission.

- Advanced support: IMV, NIPPV, CPAP, or high-flow nasal cannula.
- Index time: first advanced support record, `t0`.
- Rationale: broader ARF phenotype and sensitivity analysis for support escalation before intubation.

### Sensitivity Cohorts

- ARF criteria-based cohort using hypoxemia, hypercapnia, and room-air or oxygen-adjusted definitions where data support it.
- First ICU stay per hospitalization only versus all qualifying ICU stays.
- Exclusion or censoring approaches for death or discharge before 72 hours.

## Required CLIF Tables

Minimum for primary and secondary trajectory phenotyping:

- `patient`
- `hospitalization`
- `adt`
- `respiratory_support`
- `vitals`
- `labs`

Additional tables for outcomes and adjustment:

- `hospital_diagnosis`
- `medication_admin_continuous`
- optional site-specific severity or comorbidity inputs if harmonized later

## Trajectory Features

Core DTW features for the 0-72 hour window:

- FiO2
- PEEP
- plateau pressure
- minute ventilation
- respiratory rate

Secondary descriptive features:

- SpO2
- SpO2/FiO2 ratio
- PaO2
- PaO2/FiO2 ratio
- respiratory support state transitions
- death or discharge before 72 hours as competing events

Feature processing:

- Bin measurements hourly from `t0`.
- Aggregate multiple values per hour using median for continuous variables and dominant state for categorical support.
- Normalize FiO2 to a consistent fraction or percent scale before modeling.
- Track observedness before imputation.
- Use limited carry-forward/carry-backward only within pre-specified gaps.
- Impute remaining DTW-required values after preserving missingness QC metrics.

## DTW Phenotyping Strategy

### Development Phase

Use one or more high-volume development sites to fit candidate phenotypes.

1. Build hourly feature matrices.
2. Evaluate `k = 3` through `k = 8`.
3. Compare cluster size, medoid interpretability, silhouette or other internal validity metrics, stability under bootstrap/subsampling, and outcome separation.
4. Select a small number of clinically interpretable phenotypes for the federated primary analysis.

### Federated Application Phase

Two approaches should be compared:

- **Fixed medoid assignment:** central medoids from the development phase are distributed to sites; each site assigns patients to nearest medoid using the pre-specified DTW distance.
- **Local clustering sensitivity:** each site performs the same local DTW clustering, then maps local clusters to central phenotypes using medoid similarity and clinical profiles.

The primary federated analysis should use fixed medoid assignment if it performs acceptably, because it produces phenotype labels with consistent meaning across sites.

## Environmental Exposure Linkage

Exposure linkage should occur locally using distributed county-level exposure files and site-held geography.

Primary exposures:

- PM2.5
- NO2

Candidate exposure windows:

- 12-month pre-admission mean.
- 3-year pre-admission annual mean.
- 90-day pre-admission mean as sensitivity.

Additional contextual covariates:

- SVI
- Daymet temperature and humidity summaries
- calendar year and season

Exported exposure summaries should not include patient-level geography. Site outputs can include exposure distribution summaries and model coefficients.

## Outcomes

Primary prognostic outcomes:

- in-hospital death or hospice discharge
- ventilator-free days through day 28, if derivable
- duration of mechanical ventilation
- ICU length of stay
- hospital length of stay

Secondary outcomes:

- transition to lower support by 72 hours
- extubation by 72 hours
- re-escalation after initial improvement
- discharge disposition categories if harmonized

## Federated Modeling

### Phenotype Membership Models

Outcome: DTW phenotype membership.

Candidate model:

```text
multinomial phenotype ~ PM2.5 + NO2 + SVI + temperature + humidity
                       + age + sex + race + ethnicity + calendar year + season
```

Site export:

- coefficient estimates
- robust standard errors or model variance-covariance matrices
- model N and cluster counts
- exposure scaling parameters used at the site

Central pooling:

- fixed-effect inverse-variance meta-analysis for primary estimates
- random-effects meta-analysis as sensitivity
- heterogeneity summaries by site

### Prognostic Models

Models should test whether phenotype membership improves early prognostication beyond demographics and early severity features.

Candidate outcomes and models:

- death or hospice: logistic regression
- ICU or hospital length of stay: negative binomial, Poisson with robust SEs, or time-to-discharge model
- ventilator duration: competing-risk-aware or censored time model where feasible

Compare:

- baseline demographic model
- early clinical model
- early clinical model plus DTW phenotype
- exposure plus phenotype models

Site export:

- model coefficients
- variance estimates
- discrimination and calibration summaries
- confusion or risk-stratum summaries if used

## Site Output Package

Each site should upload a single run folder containing:

- cohort flow tables
- missingness and data density summaries
- phenotype counts
- medoid/profile summaries
- de-identified aggregate trajectory summaries by phenotype and hour
- respiratory support transition count tables by phenotype
- exposure distribution summaries
- site-level model result tables
- run metadata including git commit, CLIF version, config fields, and package versions

No patient-level rows should leave the site unless an IRB-approved sharing model is added later.

## Immediate Engineering Tasks

1. Make `02_trajectories.R` standalone or split it into feature-building and clustering scripts.
2. Save cluster assignments, feature missingness summaries, and trajectory profiles to `output/run_[SITE]_[DATE]/`.
3. Remove unclustered patients from cluster-specific figures or show them only in separate QC panels.
4. Add a central medoid export/import format for fixed phenotype assignment.
5. Add a site export manifest that defines every output allowed to leave a CLIF site.
6. Add script-level parameters for cohort choice, `k`, DTW window, feature set, and exposure window.
7. Add QC checks for FiO2 units, implausible ventilator values, duplicate timestamps, early death/discharge, and geography linkage completeness.

## Open Analytic Decisions

- Whether the primary phenotype model should use IMV-only or advanced support as the main cohort.
- Whether DTW should cluster clinical features only, with exposures modeled afterward, or include exposure histories as part of the distance metric.
- Whether early death should be represented as an absorbing trajectory state or handled outside the DTW feature matrix.
- Whether the final phenotype count should prioritize interpretability over internal cluster metrics.
- Whether site-level medoid assignment is stable enough for a primary federated analysis.
