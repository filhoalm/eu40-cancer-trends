# Early-onset cancer trends in Europe

Interactive dashboard of early-onset (<50) and late-onset (50+) cancer incidence and mortality trends
for European countries: truncated age-standardised rates, annual percent change, Lexis diagrams,
age-period-cohort estimable functions and comparative age-period-cohort analysis (MATISSE methods, ported to R).

Open `index.html` (GitHub Pages) and choose a country, cancer site and sex. Each country's data
(`data/<M49 code>.js`) loads when it is selected.

## Data sources

- Incidence: European Cancer Information System (ECIS), Joint Research Centre, European Commission -
  historical registry data, restricted per country to a common registry panel.
- Mortality: WHO Mortality Database, World Health Organization - national deaths and population.

Rates and model estimates on this site are derived from these sources; see the Methods section of the dashboard.

## Methods and code

Built with R (data.table, Rcan) by the Eu40_APC project: `code/10_country_eligibility.R`, `code/08_dashboard.R`
and `code/11_site.R`; Lexis, APC and comparative APC functions port the MATISSE toolbox (Rosenberg;
Miranda Filho & Rosenberg, Front Oncol 2024).

Generated 2026-10-06.
