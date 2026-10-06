# Early-onset cancer trends in Europe

Interactive dashboard of early-onset (<50) and late-onset (50+) cancer incidence and mortality trends
for European countries: truncated age-standardised rates, annual percent change, Lexis diagrams,
age-period-cohort estimable functions and comparative age-period-cohort analysis.

**Dashboard:** https://filhoalm.github.io/eu40-cancer-trends/ - choose a country, cancer site and sex;
each country's data (`data/<M49 code>.js`) loads when it is selected.

## Data sources

- Incidence: European Cancer Information System (ECIS), Joint Research Centre, European Commission -
  historical registry data, restricted per country to a common registry panel.
- Mortality: WHO Mortality Database, World Health Organization - national deaths and population.

Rates and model estimates on this site are derived from these sources; see Methods in the dashboard.

## Code

`code/` holds the analysis pipeline (R, data.table, Rcan). The Lexis, APC and comparative APC methods come
from [matisseR](https://github.com/filhoalm/matisseR), an R port of the MATISSE toolbox (Rosenberg;
Miranda Filho & Rosenberg, Front Oncol 2024): `remotes::install_github("filhoalm/matisseR")`.

| Script | Output |
|---|---|
| `functions_asr.R` | ECIS incidence and common registry panel, WHO mortality, truncated ASR (Rcan), APC |
| `10_country_eligibility.R` | `country_eligibility.csv`: panel >= 15 years and >= 15 consecutive WHO mortality years |
| `08_dashboard.R` + `dashboard_template.html` | dashboard, per-country data and CSV outputs with a data dictionary |
| `07_asr_onset.Rmd`, `09_lexis.Rmd` | single-country reports (ASR; MATISSE-style Lexis diagrams) |
| `11_site.R` | this site |
| `test_functions_*.R` | checks of the APC and CAPRICORN models on Germany incidence |

`data/dict/site_map.csv` maps ECIS cancer entities to WHO mortality codes (ICD-10) for the 42 sites.
Raw data are not included: the ECIS historical file (`data/ECIS/ecis_historical_data.csv`, from ECIS),
WHO deaths and population for all years (`data/who_mortality/*_full_history.rds`, prepared from the WHO
Mortality Database) and a country dictionary with UN M49 codes (`data/dict/g2022_country_info.csv`).
`renv.lock` records the R package versions.

Generated 2026-10-06.
