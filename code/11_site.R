## Static site for GitHub Pages: the trends dashboard (report/) as site/index.html + site/data/*.js, the
## analysis code (code/, current pipeline; no raw data), the site crosswalk, renv.lock, README and .nojekyll.
## Run after 08_dashboard.R. Usage (project root): Rscript code/11_site.R

library(here)

dir_site <- here("site")
code_files <- c("functions_asr.R", "07_asr_onset.Rmd", "07_render_asr_onset.R", "08_dashboard.R", "dashboard_template.html",
                "09_lexis.Rmd", "09_render_lexis.R", "10_country_eligibility.R", "11_site.R",
                "test_functions_apc.R", "test_functions_capricorn.R")
for (d in c("data", "data/dict", "code")) dir.create(file.path(dir_site, d), recursive = TRUE, showWarnings = FALSE)
unlink(c(list.files(file.path(dir_site, "data"), "\\.js$", full.names = TRUE), list.files(file.path(dir_site, "code"), full.names = TRUE)))
stopifnot(file.copy(here("report", "eu40_trends_dashboard.html"), file.path(dir_site, "index.html"), overwrite = TRUE),
          all(file.copy(list.files(here("report", "data"), "\\.js$", full.names = TRUE), file.path(dir_site, "data"))),
          all(file.copy(here("code", code_files), file.path(dir_site, "code"))),
          file.copy(here("data", "dict", "site_map.csv"), file.path(dir_site, "data", "dict"), overwrite = TRUE),
          file.copy(here("renv.lock"), dir_site, overwrite = TRUE))
file.create(file.path(dir_site, ".nojekyll"))                       # serve files as they are (no Jekyll build)

writeLines(c(
  "# Early-onset cancer trends in Europe",
  "",
  "Interactive dashboard of early-onset (<50) and late-onset (50+) cancer incidence and mortality trends",
  "for European countries: truncated age-standardised rates, annual percent change, Lexis diagrams,",
  "age-period-cohort estimable functions and comparative age-period-cohort analysis.",
  "",
  "**Dashboard:** https://filhoalm.github.io/eu40-cancer-trends/ - choose a country, cancer site and sex;",
  "each country's data (`data/<M49 code>.js`) loads when it is selected.",
  "",
  "## Data sources",
  "",
  "- Incidence: European Cancer Information System (ECIS), Joint Research Centre, European Commission -",
  "  historical registry data, restricted per country to a common registry panel.",
  "- Mortality: WHO Mortality Database, World Health Organization - national deaths and population.",
  "",
  "Rates and model estimates on this site are derived from these sources; see Methods in the dashboard.",
  "",
  "## Code",
  "",
  "`code/` holds the analysis pipeline (R, data.table, Rcan). The Lexis, APC and comparative APC methods come",
  "from [matisseR](https://github.com/filhoalm/matisseR), an R port of the MATISSE toolbox (Rosenberg;",
  "Miranda Filho & Rosenberg, Front Oncol 2024): `remotes::install_github(\"filhoalm/matisseR\")`.",
  "",
  "| Script | Output |",
  "|---|---|",
  "| `functions_asr.R` | ECIS incidence and common registry panel, WHO mortality, truncated ASR (Rcan), APC |",
  "| `10_country_eligibility.R` | `country_eligibility.csv`: panel >= 15 years and >= 15 consecutive WHO mortality years |",
  "| `08_dashboard.R` + `dashboard_template.html` | dashboard, per-country data and CSV outputs with a data dictionary |",
  "| `07_asr_onset.Rmd`, `09_lexis.Rmd` | single-country reports (ASR; MATISSE-style Lexis diagrams) |",
  "| `11_site.R` | this site |",
  "| `test_functions_*.R` | checks of the APC and CAPRICORN models on Germany incidence |",
  "",
  "`data/dict/site_map.csv` maps ECIS cancer entities to WHO mortality codes (ICD-10) for the 42 sites.",
  "Raw data are not included: the ECIS historical file (`data/ECIS/ecis_historical_data.csv`, from ECIS),",
  "WHO deaths and population for all years (`data/who_mortality/*_full_history.rds`, prepared from the WHO",
  "Mortality Database) and a country dictionary with UN M49 codes (`data/dict/g2022_country_info.csv`).",
  "`renv.lock` records the R package versions.",
  "",
  sprintf("Generated %s.", format(Sys.Date()))),
  file.path(dir_site, "README.md"))

cat(sprintf("site/: index.html + %d country files + %d code files, %.1f MB\n", length(list.files(file.path(dir_site, "data"), "\\.js$")),
            length(list.files(file.path(dir_site, "code"))),
            sum(file.size(list.files(dir_site, recursive = TRUE, full.names = TRUE, all.files = TRUE)[
              !grepl("/\\.git/", list.files(dir_site, recursive = TRUE, full.names = TRUE, all.files = TRUE))])) / 1e6))
