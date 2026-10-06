## Countries eligible for the trends dashboard (08_dashboard.R): EU-40 countries with
## ECIS incidence from a common registry panel of >= 15 years (20 preferred) and
## >= 15 consecutive years of WHO mortality (>= 3 five-year periods for the APC models).
## Usage (project root): Rscript code/10_country_eligibility.R
## Writes output/dashboard/country_eligibility.csv

library(data.table); library(here)
source(here("code", "functions_asr.R"))

min_years <- 15; pref_years <- 20
panel_entity <- "All cancer entities excluding keratinocytic skin cancers"
file_inc  <- here("data", "ECIS", "ecis_historical_data.csv")
file_mort <- here("data", "who_mortality", "data_full_history.rds")
file_ctry <- here("data", "dict", "g2022_country_info.csv")

dt_inc  <- ecis_read(file_inc, panel_entity)
dt_mort <- unique(readRDS(file_mort)[, .(id_code, year)])
eu40    <- fread(file_ctry, encoding = "Latin-1")[continent == 5, country_code]

dt_elig <- rbindlist(lapply(sort(unique(dt_inc$country)), function(cn) {
  code <- tryCatch(country_m49(file_ctry, cn), error = function(e) NA_integer_)
  d    <- ecis_incidence(dt_inc, cn, panel_entity)
  pnl  <- NULL
  for (k in c(pref_years, min_years)) if (is.null(pnl)) pnl <- tryCatch(c(registry_panel(d, k)$sel, min_years = k), error = function(e) NULL)
  m <- year_run(dt_mort[id_code == code, year])
  data.table(country = cn, country_code = code, eu40 = code %in% eu40, n_registry_ecis = uniqueN(d$registry),
             inc_y0 = pnl$y0, inc_y1 = pnl$y1, inc_years = pnl$n_years, inc_registries = pnl$n_registry,
             panel_min_years = pnl$min_years, mort_y0 = m[1], mort_y1 = m[2])
}), fill = TRUE)
dt_elig[, `:=`(mort_years = mort_y1 - mort_y0 + 1L, common_years = pmax(0L, pmin(inc_y1, mort_y1) - pmax(inc_y0, mort_y0) + 1L))]
dt_elig[, reason := fcase(is.na(country_code) | !eu40, "not in the EU-40 list",
                          is.na(inc_years), sprintf("no common registry panel of >= %d years", min_years),
                          is.na(mort_years) | mort_years < min_years, sprintf("< %d consecutive years of WHO mortality", min_years),
                          default = "")][, eligible := reason == ""]
setcolorder(dt_elig, c("country", "country_code", "eligible", "reason"))
fwrite(dt_elig, here("output", "dashboard", "country_eligibility.csv"))
print(dt_elig[, .(country, eligible, inc = sprintf("%s-%s (%s reg)", inc_y0, inc_y1, inc_registries), panel_min_years,
                  mort = sprintf("%s-%s", mort_y0, mort_y1), common_years, reason)], nrows = 50)
