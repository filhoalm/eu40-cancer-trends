## Trends dashboard: early/late-onset truncated ASR and APC by cancer site.
## Incidence: ECIS historical data, common registry panel. Mortality: WHO, national.
## Usage (project root): Rscript code/08_dashboard.R [--refresh] [Latvia Slovenia ...]
##   includes all eligible countries (output/dashboard/country_eligibility.csv, 10_country_eligibility.R) plus any named;
##   results are cached per country in output/dashboard/cache/: named countries are recomputed, --refresh recomputes all
##   (needed after changing the analysis functions; template-only changes just reassemble from the cache)
## Writes output/dashboard/*.csv, report/eu40_trends_dashboard.html and report/data/<M49>.js (one per country;
## keep the data folder beside the HTML)

library(data.table); library(Rcan); library(here); library(jsonlite)
source(here("code", "functions_asr.R"))
library(matisseR)   # MATISSE port: Lexis/rates, APC, CAPRICORN (remotes::install_github("filhoalm/matisseR"))

args      <- commandArgs(trailingOnly = TRUE)
refresh   <- "--refresh" %in% args; named <- setdiff(args, "--refresh")
file_elig <- here("output", "dashboard", "country_eligibility.csv")
countries <- union(if (file.exists(file_elig)) fread(file_elig)[eligible == TRUE, country], named)
if (!length(countries)) stop("Run code/10_country_eligibility.R first, or name the countries")
age_cut <- 50; pop_base <- "SEGI"; panel_years <- c(20, 15)   # common registry panel: preferred, minimum length
lexis_ages <- 5:17   # ages 20-84: 5 x 5 MATISSE cells, and 1 x 1 cells interpolated with i5

file_inc  <- here("data", "ECIS", "ecis_historical_data.csv")
file_mort <- here("data", "who_mortality", "data_full_history.rds")
file_pop  <- here("data", "who_mortality", "population_full_history.rds")
file_ctry <- here("data", "dict", "g2022_country_info.csv")
file_site <- here("data", "dict", "site_map.csv")
dir_out   <- here("output", "dashboard")
dir_cache <- file.path(dir_out, "cache")
dir.create(dir_cache, showWarnings = FALSE, recursive = TRUE)

dt_site  <- fread(file_site, colClasses = "character", na.strings = NULL)
dt_ent   <- dt_site[, .(entity = strsplit(ecis_entity, ";")[[1]]), by = site]               # colorectum = Colon + Rectum
dt_codes <- dt_site[who_codes != "", .(cancer_code = as.integer(strsplit(who_codes, ";")[[1]])), by = site]
dt_sex   <- dt_site[, .(sex = as.integer(strsplit(sexes, ";")[[1]])), by = site]
site_both <- dt_sex[, .N, by = site][N == 2, site]
dt_ecis  <- NULL                                       # ECIS read once, only if a country is (re)computed
skip <- function(what, expr) tryCatch(expr, error = function(e) { message("  skipped ", what, ": ", conditionMessage(e)); NULL })

country_res <- function(country_select) skip(country_select, {
  message(format(Sys.time(), "%H:%M:%S "), country_select)
  if (is.null(dt_ecis)) dt_ecis <<- ecis_read(file_inc, unique(c(dt_ent$entity, "All cancer entities excluding keratinocytic skin cancers")))
  code <- country_m49(file_ctry, country_select)
  inc  <- NULL
  for (k in panel_years) if (is.null(inc)) inc <- tryCatch(ecis_country_incidence(dt_ecis, country_select, dt_ent, k), error = function(e) NULL)
  if (is.null(inc)) stop(sprintf("no common registry panel of >= %d years", min(panel_years)))
  pnl  <- inc$panel
  mort <- who_mortality(file_mort, file_pop, code, dt_codes)
  mort <- mort[year %between% year_run(year)]                                       # longest run of consecutive years

  dt <- rbind(inc$data[, .(indicator = "Incidence", site, year, sex, age, cases, py)],
              mort[, .(indicator = "Mortality", site, year, sex, age, cases, py)], use.names = TRUE)
  dt <- merge(dt, dt_sex, by = c("site", "sex"))
  dt <- rbind(dt, dt[site %in% site_both, .(sex = 0L, cases = sum(cases), py = sum(py)), by = .(site, indicator, year, age)],
              use.names = TRUE)

  dt_asr <- smooth_onset(asr_onset(dt, c("site", "indicator", "sex", "year"), age_cut, pop_base),
                         c("site", "indicator", "sex", "onset"))

  ## 5 x 5 Lexis diagrams (MATISSE rates: chunk to 5-year periods, fill zero cells); events kept unfilled
  dt_lex <- dt[age %in% lexis_ages, skip(paste("Lexis", paste(.BY, collapse = " ")), {
    R0 <- rates_chunk(rates_make(.SD, lexis_ages), per_block = 5)
    R  <- rates_fill(R0)
    rates_long(R)[, .(age = age_lo, period = per_lo, cohort = coh_mid, events = rates_long(R0)$events,
                      py = offset, rate, fill = R$fill_value)]
  }), by = .(site, indicator, sex)]

  ## APC model (MATISSE @apc, WLS) on the 5 x 5 cells: the four basic estimable functions and local drifts
  dt_ef <- dt[age %in% lexis_ages, skip(paste("APC", paste(.BY, collapse = " ")), {
    R <- rates_fill(rates_chunk(rates_make(.SD, lexis_ages), per_block = 5))
    if (ncol(R$Events) >= 3) { M <- apc_fit(R); apc_efs(M)[, `:=`(net_drift = M$net_drift, s2 = M$s2)] }
  }), by = .(site, indicator, sex)]

  ## comparative APC (MATISSE @capricorn) on the 5 x 5 cells: males vs females by site and indicator, and
  ## mortality vs incidence by site and sex over their common years; tests as in MATISSE (Poisson) and
  ## overdispersion-adjusted; estimable functions under NPH; the last stratum is the reference
  rt <- function(d) rates_fill(rates_chunk(rates_make(d, lexis_ages), per_block = 5))
  cap_tabs <- function(d, by, strata, id) skip(paste(c("CAPRICORN", unlist(id)), collapse = " "), {
    Rs <- if (all(strata %in% d[[by]])) lapply(strata, function(s) rt(d[get(by) == s]))
    if (!is.null(Rs) && ncol(Rs[[1]]$Events) >= 3) {
      id$years <- sprintf("%d-%d", Rs[[1]]$periods[1], tail(Rs[[1]]$periods, 1) - 1)
      S <- lapply(c(FALSE, TRUE), capricorn, Rs = Rs)
      tst <- rbindlist(lapply(1:2, function(i) rbind(rbindlist(S[[i]]$tests, idcol = "type", fill = TRUE),
                                                     S[[i]]$aic[, .(type = "aic", test = model, stat = AICc, df, delta, rank)], fill = TRUE)[
                                                       , `:=`(overdispersion = i == 2, phi = S[[i]]$phi)]))
      list(tests = tst[, (names(id)) := id], ef = rbindlist(lapply(c("lac", "cac", "ftt", "fcp", "ld"), cap_ef, S = S[[1]]))[, (names(id)) := id])
    }
  })
  dc  <- dt[age %in% lexis_ages]
  yrs <- dc[, .(y0 = min(year), y1 = max(year)), by = .(site, sex, indicator)][
    , .(y0 = max(y0), y1 = min(y1), .N), by = .(site, sex)][N == 2 & y1 - y0 >= 14]          # >= 3 common 5-year periods
  cap <- Filter(Negate(is.null), c(
    unlist(lapply(site_both, function(s) lapply(c("Incidence", "Mortality"), function(ind)
      cap_tabs(dc[site == s & indicator == ind], "sex", c(Male = 1L, Female = 2L),
               list(comparison = "sex", site = s, indicator = ind, sex = NA_integer_)))), recursive = FALSE),
    lapply(seq_len(nrow(yrs)), function(i)
      cap_tabs(dc[yrs[i], on = .(site, sex)][year %between% c(yrs$y0[i], yrs$y1[i])], "indicator",
               c(Mortality = "Mortality", Incidence = "Incidence"),
               list(comparison = "mortality_incidence", site = yrs$site[i], indicator = NA_character_, sex = yrs$sex[i])))))

  ## 1 x 1 Lexis diagrams: single-year ages interpolated from 5-year groups (MATISSE i5, cubic), annual periods
  dt_lex1 <- dt[age %in% lexis_ages, skip(paste("Lexis 1x1", paste(.BY, collapse = " ")), {
    R <- rates_fill(rates_i5(rates_make(.SD, lexis_ages), "cubic"))
    rates_long(R)[, .(age = age_lo, period = per_lo, cohort = coh_mid, events, py = offset, rate)]
  }), by = .(site, indicator, sex)]

  list(asr   = dt_asr[, `:=`(country = country_select, country_code = code)],
       apc   = apc_onset(dt_asr, c("site", "indicator", "sex", "onset"))[, `:=`(country = country_select, country_code = code)],
       lexis = dt_lex[, `:=`(country = country_select, country_code = code)],
       lexis1 = dt_lex1[, `:=`(country = country_select, country_code = code)],
       ef    = dt_ef[, `:=`(country = country_select, country_code = code)],
       cap_tests = rbindlist(lapply(cap, `[[`, "tests"))[, `:=`(country = country_select, country_code = code)],
       cap_ef    = rbindlist(lapply(cap, `[[`, "ef"))[, `:=`(country = country_select, country_code = code)],
       panel = pnl$sel[, .(country = country_select, country_code = code, y0, y1, n_registry,
                           registries = sapply(registries, paste, collapse = "; "))])
})

## cached per country; named countries and --refresh recompute
ls_res <- Filter(Negate(is.null), lapply(countries, function(cn) {
  f <- file.path(dir_cache, paste0(gsub("[^A-Za-z]+", "_", cn), ".rds"))
  if (!refresh && !cn %in% named && file.exists(f)) return(readRDS(f))
  r <- country_res(cn)
  if (!is.null(r)) saveRDS(r, f)
  r
}))

dt_asr   <- rbindlist(lapply(ls_res, `[[`, "asr"))[, .(country, country_code, site, indicator, sex, onset, year,
                                                       cases, py, asr, se, asr_lo, asr_hi, asr_smooth)]
dt_apc   <- rbindlist(lapply(ls_res, `[[`, "apc"))[, .(country, country_code, site, indicator, sex, onset, y0, y1,
                                                       apc, apc_lo, apc_hi, fit0, fit1)]
dt_panel <- rbindlist(lapply(ls_res, `[[`, "panel"))
dt_ef    <- rbindlist(lapply(ls_res, `[[`, "ef"))[, .(country, country_code, site, indicator, sex, ef, ref, x, value, lo, hi,
                                                    net_drift, s2)]
dt_cap_t <- rbindlist(lapply(ls_res, `[[`, "cap_tests"))[, .(country, country_code, comparison, site, indicator, sex, years,
                                                            type, test, stat, df, p, delta, rank, overdispersion, phi)]
dt_cap_e <- rbindlist(lapply(ls_res, `[[`, "cap_ef"))[, .(country, country_code, comparison, site, indicator, sex, years,
                                                         model, ef, stratum, type, ref, x, value, lo, hi)]
dt_lex1  <- rbindlist(lapply(ls_res, `[[`, "lexis1"))[, .(country, country_code, site, indicator, sex, age, period, cohort,
                                                        events, py, rate)]
dt_lex   <- rbindlist(lapply(ls_res, `[[`, "lexis"))[, .(country, country_code, site, indicator, sex, age, period, cohort,
                                                       events, py, rate, fill)]
setorder(dt_asr, country, site, indicator, sex, onset, year)

missing_site <- dt_site[!site %in% dt_asr[indicator == "Incidence", site], site]
if (length(missing_site)) warning("No ECIS incidence for site(s): ", paste(missing_site, collapse = ", "))

fwrite(dt_asr, file.path(dir_out, "asr_onset_sites.csv"))
fwrite(dt_apc, file.path(dir_out, "apc_onset_sites.csv"))
fwrite(dt_panel, file.path(dir_out, "incidence_registry_panel.csv"))
fwrite(dt_lex, file.path(dir_out, "lexis_cells.csv"))
fwrite(dt_lex1, file.path(dir_out, "lexis_cells_1x1.csv.gz"))
fwrite(dt_ef, file.path(dir_out, "apc_ef_sites.csv"))
fwrite(dt_cap_t, file.path(dir_out, "capricorn_tests.csv"))
fwrite(dt_cap_e, file.path(dir_out, "capricorn_ef.csv"))
fwrite(data.table(
  file = c(rep("asr_onset_sites.csv", 14), rep("apc_onset_sites.csv", 13), rep("lexis_cells.csv", 12),
           rep("lexis_cells_1x1.csv.gz", 11), rep("apc_ef_sites.csv", 13), rep("capricorn_tests.csv", 16), rep("capricorn_ef.csv", 16)),
  variable = c(names(dt_asr), names(dt_apc), names(dt_lex), names(dt_lex1), names(dt_ef), names(dt_cap_t), names(dt_cap_e)),
  description = c("ECIS country name", "UN M49 country code", "Site id (data/dict/site_map.csv)",
                  "Incidence (ECIS registry panel) or Mortality (WHO, national)", "0 = both, 1 = male, 2 = female",
                  sprintf("early = ages 0-%d, late = ages %d+", age_cut - 1, age_cut), "Calendar year",
                  "Cases or deaths in the age range", "Person-years in the age range",
                  "Truncated ASR per 100,000, Segi world standard", "Standard error of asr",
                  "Lower 95% CI (floored at 0)", "Upper 95% CI",
                  "Smoothed ASR: loess (degree 1, span 0.75) on log ASR; NA if < 5 non-zero years",
                  "ECIS country name", "UN M49 country code", "Site id", "Incidence or Mortality", "0 = both, 1 = male, 2 = female",
                  "Onset group", "First year of the series", "Last year of the series",
                  "Annual percent change: weighted log-linear fit of ASR on year", "Lower 95% CI of apc", "Upper 95% CI of apc",
                  "Fitted ASR in the first year (y0)", "Fitted ASR in the last year (y1)",
                  "ECIS country name", "UN M49 country code", "Site id", "Incidence or Mortality", "0 = both, 1 = male, 2 = female",
                  "Lower bound of the 5-year age group", "First year of the 5-year period",
                  "Cohort midpoint birth year (period midpoint - age midpoint); a 5 x 5 cell spans 10 birth years",
                  "Cases or deaths (not filled)", "Person-years",
                  "Rate per 100,000; zero cells use the fill value", "MATISSE fill value for zero cells of this Lexis diagram (>= 0.5)",
                  "ECIS country name", "UN M49 country code", "Site id", "Incidence or Mortality", "0 = both, 1 = male, 2 = female",
                  "Single year of age (interpolated from 5-year groups, MATISSE i5 cubic)", "Calendar year",
                  "Cohort midpoint birth year (year - age); a 1 x 1 cell spans 2 birth years", "Interpolated cases or deaths (zero cells filled)",
                  "Interpolated person-years", "Rate per 100,000",
                  "ECIS country name", "UN M49 country code", "Site id", "Incidence or Mortality", "0 = both, 1 = male, 2 = female",
                  "Estimable function of the APC model (MATISSE): lac = longitudinal age curve, cac = cross-sectional age curve, ftt = fitted temporal trends, fcp = fitted cohort pattern, ld = local drifts",
                  "Reference value: cohort (lac), period (cac) or age (ftt, fcp) midpoint; MATISSE default = central group; NA for ld",
                  "Age (lac, cac, ld), period (ftt) or cohort (fcp) midpoint",
                  "Rate per 100,000 (lac, cac, ftt, fcp; 5 x 5 cells, ages 20-84) or annual % change (ld)", "Lower 95% CI (delta method)", "Upper 95% CI",
                  "Net drift, % per year", "Overdispersion: max(1, deviance / df)",
                  "ECIS country name", "UN M49 country code",
                  "Comparison (MATISSE CAPRICORN): sex = males vs females; mortality_incidence = mortality vs incidence", "Site id",
                  "Incidence or Mortality (sex comparison); NA for mortality_incidence",
                  "Sex for mortality_incidence (0 = both, 1 = male, 2 = female); NA for the sex comparison",
                  "Calendar years of the 5-year periods (mortality_incidence: years common to both, most recent kept)",
                  "global = PH model vs NPH; homogeneity = one APC block equal vs NPH; composite = Bonferroni combinations (MATISSE); aic = model selection",
                  "Model or block: PHL, PHT, PHX, PHA, NPH; LAT, NetDrift, Quad/HiOrd Age/Per/Coh, CAT; Par-LAC ... Par-FCP",
                  "Likelihood-ratio statistic (deviance difference, divided by phi) or AICc", "Degrees of freedom (aic: model df)",
                  "p-value", "AICc minus the smallest AICc (aic rows)", "AICc rank (aic rows)",
                  "FALSE = Poisson, as in MATISSE; TRUE = deviances divided by phi (quasi-Poisson, QAICc)",
                  "NPH dispersion max(1, deviance / df) used when overdispersion = TRUE",
                  "ECIS country name", "UN M49 country code",
                  "Comparison (MATISSE CAPRICORN): sex = males vs females; mortality_incidence = mortality vs incidence", "Site id",
                  "Incidence or Mortality (sex comparison); NA for mortality_incidence",
                  "Sex for mortality_incidence (0 = both, 1 = male, 2 = female); NA for the sex comparison",
                  "Calendar years of the 5-year periods (mortality_incidence: years common to both, most recent kept)",
                  "Joint model the estimates come from (NPH)",
                  "Estimable function: lac, cac, ftt, fcp, ld (see apc_ef_sites.csv)",
                  "Stratum (Male, Female, Mortality, Incidence) or contrast 'A vs B' (reference B: Female, Incidence)",
                  "ef = estimable function of the stratum; ratio = rate ratio A / B (lac, cac, ftt, fcp) or ratio of annual changes in % per year (ld)",
                  "Reference value (see apc_ef_sites.csv)", "Age, period or cohort midpoint",
                  "Rate per 100,000, annual % change, rate ratio or % per year (see type)", "Lower 95% CI (delta method, Poisson)", "Upper 95% CI"))
  , file.path(dir_out, "dictionary.csv"))

## dashboard: one script per country, EU40_LOAD(<compact JSON>), loaded on demand with a <script> tag so it also
## works from file:// (series as arrays, 1 x 1 Lexis rates as
## differences of round(100 * log(rate)), i.e. 1% precision); sites, panels and countries embedded in the HTML
dir_data <- here("report", "data")
dir.create(dir_data, showWarnings = FALSE)
unlink(list.files(dir_data, "\\.(js|json|json\\.gz)$", full.names = TRUE))
dt_cty <- unique(dt_asr[, .(country, code = country_code)])[, file := sprintf("data/%d.js", code)][order(country)]
for (i in seq_len(nrow(dt_cty))) {
  cn <- dt_cty$country[i]
  d  <- list(
    country = cn,
    asr     = dt_asr[country == cn][order(year), { stopifnot(all(diff(year) == 1))
                .(y0 = year[1], a = list(round(asr, 3)), s = list(signif(asr_smooth, 4)), n = list(round(cases))) },
                by = .(site, indicator, sex, onset)],
    apc     = dt_apc[country == cn, .(site, indicator, sex, onset, y0, y1, apc = round(apc, 2), lo = round(apc_lo, 2), hi = round(apc_hi, 2),
                                      fit0 = signif(fit0, 6), fit1 = signif(fit1, 6))],
    lexis   = dt_lex[country == cn, .(p0 = min(period), P = uniqueN(period), f = signif(fill[1], 4),
                                      e = list(round(events)), py = list(round(py))), by = .(site, indicator, sex)],
    apcef   = dt_ef[country == cn, .(x = list(x), r = list(signif(value, 4))), by = .(site, indicator, sex, ef, ref, nd = round(net_drift, 2))],
    capef   = dt_cap_e[country == cn & type == "ratio" & ef %in% c("fcp", "ld"), .(x = list(x), r = list(signif(value, 4))),
                       by = .(comparison, site, indicator = fcoalesce(indicator, "-"), sex = fcoalesce(sex, -1L), ef, ref, years)],
    captest = merge(dt_cap_t[country == cn & type == "aic", .(comparison, site, indicator = fcoalesce(indicator, "-"), sex = fcoalesce(sex, -1L),
                                                              years, od = overdispersion, model = test, d = round(delta, 1))],
                    dt_cap_t[country == cn & type == "global", .(comparison, site, indicator = fcoalesce(indicator, "-"), sex = fcoalesce(sex, -1L),
                                                                 years, od = overdispersion, model = test, p = signif(p, 3))], all.x = TRUE),
    lexis1  = dt_lex1[country == cn, .(p0 = min(period), P = uniqueN(period), r = list(diff(c(0L, as.integer(round(100 * log(rate))))))),
                      by = .(site, indicator, sex)])
  writeLines(paste0("EU40_LOAD(", toJSON(d, dataframe = "columns", auto_unbox = TRUE, na = "null"), ");"),
             file.path(dir_data, sprintf("%d.js", dt_cty$code[i])), useBytes = TRUE)
}
payload <- list(
  generated = format(Sys.Date()), age_cut = age_cut, pop_base = "World (Segi 1960)", lexis_a0 = 5 * (min(lexis_ages) - 1),
  sites  = dt_site[, .(id = site, label = site_label, group, eo = eo_rising == "1", sexes, comparability, note, icd10,
                       ecis = gsub(";", " + ", ecis_entity))],
  panels = dt_panel, countries = dt_cty[, .(country, code, file)])
json <- gsub("</", "<\\/", toJSON(payload, dataframe = "columns", auto_unbox = TRUE, na = "null"), fixed = TRUE)
html <- strsplit(paste(readLines(here("code", "dashboard_template.html"), encoding = "UTF-8"), collapse = "\n"),
                 "__DASHBOARD_DATA__", fixed = TRUE)[[1]]
stopifnot(length(html) == 2)
writeLines(paste0(html[1], json, html[2]), here("report", "eu40_trends_dashboard.html"), useBytes = TRUE)

cat(sprintf("%d countries, %d sites | %s rows ASR, %d APC series\n", nrow(dt_cty), uniqueN(dt_asr$site),
            format(nrow(dt_asr), big.mark = ","), nrow(dt_apc)))
print(dt_panel[, .(country, y0, y1, n_registry)])
