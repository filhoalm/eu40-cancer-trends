## Render code/07_asr_onset.Rmd for one or more ECIS countries.
## Usage (from the project root): Rscript code/07_render_asr_onset.R Germany "United Kingdom"

library(here)

countries <- commandArgs(trailingOnly = TRUE)
if (!length(countries)) countries <- "Germany"

for (country in countries)
  rmarkdown::render(here("code", "07_asr_onset.Rmd"), params = list(country = country),
                    output_file = sprintf("asr_onset_%s.html", tolower(gsub("\\W+", "_", country))),
                    output_dir = here("report"), intermediates_dir = tempdir(), envir = new.env())
