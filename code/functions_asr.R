## Shared functions: ECIS incidence, common registry panel, WHO mortality,
## early/late-onset truncated ASR (Rcan) and log-linear APC.

## UN M49 code for an ECIS country name (3 ECIS names differ from the dictionary)
country_m49 <- function(file_dict, country_select) {
  fix <- c("France" = "France (metropolitan)", "Netherlands" = "The Netherlands",
           "Bosnia and Herzegovina" = "Bosnia Herzegovina")
  lab <- if (country_select %in% names(fix)) fix[[country_select]] else country_select
  code <- fread(file_dict, encoding = "Latin-1")[country_label == lab, country_code]
  if (length(code) != 1) stop(sprintf("No UN M49 code for '%s'", country_select))
  code
}

## longest run of consecutive years: c(first, last)
year_run <- function(y) {
  if (!length(y)) return(c(NA_integer_, NA_integer_))
  y <- sort(unique(y)); g <- cumsum(c(1, diff(y) != 1)); r <- y[g == which.max(tabulate(g))]
  c(min(r), max(r))
}

## ECIS historical incidence rows for some entities, all countries (read the 1.5 GB file once)
ecis_read <- function(file, entities) {
  fread(file, select = c("indicator", "registry", "country", "sex", "cancer_entity_label",
                         "year", "age_min", "age_max", "cases", "population"))[
    indicator == "Incidence" & cancer_entity_label %in% entities]
}

## ECIS historical incidence, cleaned: one row per registry, entity, year, sex, age (1-18);
## src = file path or the output of ecis_read
ecis_incidence <- function(src, country_select, entities) {
  dt <- (if (is.data.table(src)) src else ecis_read(src, entities))[country == country_select & cancer_entity_label %in% entities]
  if (!nrow(dt)) stop(sprintf("No ECIS incidence for '%s'", country_select))
  dt <- dt[age_max < 85 | (age_min == 85 & age_max == 998)]   # 85-89/90-94/95+ rows repeat the 85+ total
  dt <- dt[!dt[is.na(cases)], on = .(registry, cancer_entity_label, year, sex)]   # not reported (e.g. Latvia, 2 entities from 2014)
  dt[, `:=`(age = age_min %/% 5L + 1L, sex = fifelse(sex == "Male", 1L, 2L))]
  dt[, n_age := uniqueN(age), by = .(registry, cancer_entity_label, year, sex)]
  n_drop <- uniqueN(dt[n_age < 18], by = c("registry", "cancer_entity_label", "year", "sex"))
  dt <- dt[n_age == 18, .(registry, entity = cancer_entity_label, year, sex, age, cases, py = population)]
  stopifnot(!anyNA(dt$py))
  setattr(dt, "n_drop", n_drop)
}

## Common registry panel: registry units reporting every year of [y0, y1];
## selected = most person-years among periods of >= min_years. A registry whose
## covered population changes by > max_jump in a year (e.g. North Rhine-Westphalia
## going from one region to the whole state in 2010) becomes a new unit from that year.
registry_panel <- function(dt, min_years, max_jump = 0.2) {
  dt_ry <- dt[, .(py = sum(py)), by = .(registry, year)][order(registry, year)]
  dt_ry[, seg := cumsum(c(TRUE, abs(py[-1] / py[-.N] - 1) > max_jump | diff(year) > 1)), by = registry]
  dt_ry[, n_seg := max(seg), by = registry]
  dt_ry[, unit := if (n_seg[1] == 1) registry else sprintf("%s (%d-%d)", registry, min(year), max(year)), by = .(registry, seg)]
  yrs   <- sort(unique(dt_ry$year))
  dt_cand <- CJ(y0 = yrs, y1 = yrs)[y1 >= y0][, {
    r <- dt_ry[year >= y0 & year <= y1, .(n = .N, py = sum(py)), by = unit][n == y1 - y0 + 1L]
    .(n_registry = nrow(r), py = sum(r$py), registries = list(sort(r$unit)))
  }, by = .(y0, y1)][n_registry > 0][, n_years := y1 - y0 + 1L]
  dt_sel <- dt_cand[n_years >= min_years][order(-py)][1]
  if (is.na(dt_sel$y0)) stop(sprintf("No common registry panel of >= %d years", min_years))
  dt_ry[, used := unit %in% dt_sel$registries[[1]] & year %between% c(dt_sel$y0, dt_sel$y1)]
  list(cand = dt_cand, sel = dt_sel, coverage = dt_ry[, c("seg", "n_seg") := NULL][])
}

## ECIS incidence for a country: sites (dt_ent: site, entity) summed over the common registry panel.
## Sites made of several entities add cases but count each registry's population once.
ecis_country_incidence <- function(file, country_select, dt_ent, panel_min_years,
                                   panel_entity = "All cancer entities excluding keratinocytic skin cancers") {
  dt  <- ecis_incidence(file, country_select, unique(c(dt_ent$entity, panel_entity)))
  pnl <- registry_panel(dt[entity == panel_entity], panel_min_years)
  dt  <- merge(dt, dt_ent, by = "entity", allow.cartesian = TRUE)
  stopifnot(dt[, uniqueN(py), by = .(site, registry, year, sex, age)][, all(V1 == 1)])
  dt  <- dt[, .(cases = sum(cases), py = py[1]), by = .(site, registry, year, sex, age)]
  dt  <- merge(dt, pnl$coverage[used == TRUE, .(registry, year)], by = c("registry", "year"))
  list(data = dt[, .(cases = sum(cases), py = sum(py)), by = .(site, year, sex, age)], panel = pnl)
}

## WHO deaths (GLOBOCAN codes summed per site) with WHO population; zero-filled age grid
who_mortality <- function(file_mort, file_pop, code, dt_codes) {
  dt_d <- merge(readRDS(file_mort)[id_code == code], dt_codes, by = "cancer_code", allow.cartesian = TRUE)[
    , .(cases = sum(cases)), by = .(site, year, sex, age)]
  if (!nrow(dt_d)) stop(sprintf("No WHO mortality for country_code %d", code))
  dt <- CJ(site = unique(dt_codes$site), year = unique(dt_d$year), sex = 1:2, age = 1:18)
  dt <- merge(dt, dt_d, by = c("site", "year", "sex", "age"), all.x = TRUE)[is.na(cases), cases := 0]
  merge(dt, readRDS(file_pop)[id_code == code, .(year, sex, age, py)], by = c("year", "sex", "age"))
}

## Truncated ASR with Rcan::csu_asr: early (< age_cut) and late (>= age_cut) onset, with 95% CI
asr_onset <- function(dt, group_by, age_cut = 50, pop_base = "SEGI") {
  grp_cut <- age_cut %/% 5
  stopifnot(age_cut %% 5 == 0, grp_cut %in% 1:17)
  dt <- dt[, c(group_by, "age", "cases", "py"), with = FALSE]
  rbind(Rcan::csu_asr(dt, group_by = group_by, first_age = 1, last_age = grp_cut, pop_base = pop_base, var_st_err = "se")[, onset := "early"],
        Rcan::csu_asr(dt, group_by = group_by, first_age = grp_cut + 1, last_age = 18, pop_base = pop_base, var_st_err = "se")[, onset := "late"])[
    , `:=`(asr_lo = pmax(asr - 1.96 * se, 0), asr_hi = asr + 1.96 * se)][]
}

## Smoothed ASR: LOWESS-type local linear fit (loess, degree 1) on log ASR; zero rates excluded
smooth_onset <- function(dt_asr, group_by, span = 0.75) {
  dt_asr[, asr_smooth := NA_real_]
  dt_asr[asr > 0, asr_smooth := if (.N >= 5) exp(predict(loess(log(asr) ~ year, span = span, degree = 1))) else NA_real_,
         by = group_by][]
}

## Annual percent change: weighted log-linear fit of ASR on year (weights = inverse variance of log ASR);
## fit0/fit1 = fitted ASR in the first/last year (the fitted line is straight on a log scale)
apc_onset <- function(dt_asr, group_by) {
  dt_asr[asr > 0 & se > 0, {
    if (.N < 5) list(apc = NA_real_, apc_lo = NA_real_, apc_hi = NA_real_, y0 = min(year), y1 = max(year),
                     fit0 = NA_real_, fit1 = NA_real_)
    else {
      f <- lm(log(asr) ~ year, weights = (asr / se)^2)
      b <- coef(summary(f))[2, 1:2]
      list(apc = 100 * expm1(b[[1]]), apc_lo = 100 * expm1(b[[1]] - 1.96 * b[[2]]),
           apc_hi = 100 * expm1(b[[1]] + 1.96 * b[[2]]), y0 = min(year), y1 = max(year),
           fit0 = exp(sum(coef(f) * c(1, min(year)))), fit1 = exp(sum(coef(f) * c(1, max(year)))))
    }
  }, by = group_by]
}
