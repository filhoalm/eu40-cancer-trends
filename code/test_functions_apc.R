## Checks of the matisseR APC model (MATISSE @apc port) on Germany all-cancer incidence 5 x 5
## (the package's own tests use the MATISSE example data).
## Usage (project root): Rscript code/test_functions_apc.R

library(data.table); library(here)
source(here("code", "functions_asr.R")); library(matisseR)

ent <- data.table(site = "all", entity = "All cancer entities excluding keratinocytic skin cancers")
inc <- ecis_country_incidence(here("data", "ECIS", "ecis_historical_data.csv"), "Germany", ent, 20)$data
R   <- rates_fill(rates_chunk(rates_make(inc[, .(year, age, cases, py)], 5:17), per_block = 5))
L   <- rates_long(R)[order(per_lo, age_lo)]                      # column-major, as c(R$Events)
tol <- 1e-8

for (m in c("wls", "poisson")) {
  M <- apc_fit(R, m); G <- M$G
  u <- drop(G$X %*% M$B)
  ## 1. the design spans the full APC model: same fit as age + period + cohort factors
  f <- if (m == "wls") fitted(lm(log(events / offset) ~ factor(age_mid) + factor(per_mid) + factor(coh_mid), weights = events, data = L))
       else predict(glm(events ~ factor(age_mid) + factor(per_mid) + factor(coh_mid) + offset(log(offset)), poisson, data = L)) - log(L$offset)
  stopifnot(qr(G$X)$rank == G$A + G$P + G$C - 3, max(abs(u - f)) < tol)

  ## 2. each basic EF equals the fitted log rates on its reference line, other deviations removed
  E  <- apc_efs(M); r <- apc_ref(M); lu <- log(1e5) + u
  pd <- drop(G$X[, c(5, G$pp)] %*% M$B[c(5, G$pp)]); cd <- drop(G$X[, c(6, G$pc)] %*% M$B[c(6, G$pc)])
  chk <- function(e, on, xv, dev) { v <- E[ef == e]; max(abs(log(v$value[match(xv[on], v$x)]) - (lu - dev)[on])) }
  stopifnot(chk("lac", L$coh_mid == G$coh[r["c"]], L$age_mid, pd) < tol,
            chk("cac", L$per_mid == G$per[r["p"]], L$age_mid, cd) < tol,
            chk("ftt", L$age_mid == G$age[r["a"]], L$per_mid, cd) < tol,
            chk("fcp", L$age_mid == G$age[r["a"]], L$coh_mid, pd) < tol)

  ## 3. local drifts = least-squares slope over periods of the fitted log rates at each age
  sl <- sapply(seq_len(G$A), function(j) coef(lm(u[seq(j, G$A * G$P, G$A)] ~ G$per))[2])
  stopifnot(max(abs(E[ef == "ld", value] - 100 * (exp(sl) - 1))) < tol)
  cat(sprintf("%s: ok | net drift %.2f%%/yr | overdispersion %.1f\n", m, M$net_drift, M$s2))
}
