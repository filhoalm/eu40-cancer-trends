## Checks of matisseR CAPRICORN (MATISSE @capricorn port) on Germany incidence 5 x 5:
## colon vs rectum (both sexes) and males vs females (all cancers).
## Usage (project root): Rscript code/test_functions_capricorn.R

library(data.table); library(here)
source(here("code", "functions_asr.R")); library(matisseR)

ent <- data.table(site = c("colon", "rectum", "all"), entity = c("Colon", "Rectum", "All cancer entities excluding keratinocytic skin cancers"))
inc <- ecis_country_incidence(here("data", "ECIS", "ecis_historical_data.csv"), "Germany", ent, 20)$data
rt  <- function(d) rates_fill(rates_chunk(rates_make(d[, .(year, age, cases, py)], 5:17), per_block = 5))
tol <- 1e-6

for (cmp in list(list(colon = rt(inc[site == "colon"]), rectum = rt(inc[site == "rectum"])),
                 list(male = rt(inc[site == "all" & sex == 1]), female = rt(inc[site == "all" & sex == 2])))) {
  S <- capricorn(cmp); G <- S$G; M <- ncol(G$X)

  ## 1. NPH = separate APC models (Poisson)
  sep <- lapply(cmp, apc_fit, method = "poisson")
  stopifnot(max(abs(S$fits$NPH$B - unlist(lapply(sep, `[[`, "B")))) < tol,
            abs(S$fits$NPH$DEV - sum(sapply(sep, `[[`, "DEV"))) < tol)

  ## 2. PH models = factor models with shared / stratum-specific age, period, cohort effects
  L <- rbindlist(lapply(names(cmp), function(s) rates_long(cmp[[s]])[order(per_lo, age_lo)][, s := s]))
  L[, `:=`(fa = factor(age_mid), fp = factor(per_mid), fc = factor(coh_mid), s = factor(s, names(cmp)))]
  dev <- function(fm) suppressWarnings(glm(fm, poisson, data = L, offset = log(offset)))$deviance
  ref <- c(PHL = dev(events ~ fa + fp + s:fc), PHT = dev(events ~ s:fa + fp + fc), PHX = dev(events ~ fa + s:fp + fc),
           PHA = dev(events ~ s + fa + fp + fc), NPH = dev(events ~ s:(fa + fp + fc)))
  stopifnot(max(abs(S$aic[match(names(ref), model), DEV] - ref)) < 1e-4)

  ## 3. stratum EFs under NPH = single-model EFs; ratios = EF ratios
  for (cp in c("lac", "cac", "ftt", "fcp", "ld")) {
    e <- cap_ef(S, cp)
    for (g in seq_along(cmp)) stopifnot(max(abs(e[type == "ef" & stratum == names(cmp)[g], value] / apc_ef(sep[[g]], cp)$value - 1)) < tol)
    r <- e[type == "ratio", value]; a <- e[stratum == names(cmp)[1], value]; b <- e[stratum == names(cmp)[2], value]
    stopifnot(max(abs(r - if (cp == "ld") 100 * ((1 + a / 100) / (1 + b / 100) - 1) else a / b)) < tol)
  }
  cat(sprintf("%s: ok | best (AICc) %s | global PH p: %s\n", paste(names(cmp), collapse = " vs "), S$best,
              paste(sprintf("%s %.2g", S$tests$global$test, S$tests$global$p), collapse = ", ")))
}
