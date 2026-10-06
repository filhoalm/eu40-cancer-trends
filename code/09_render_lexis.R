## Render code/09_lexis.Rmd (MATISSE-style Lexis diagrams) for one or more ECIS countries.
## Usage (from the project root): Rscript code/09_render_lexis.R Germany [site] [sex]
##   site: id from data/dict/site_map.csv (default "all"); sex: 0 both, 1 male, 2 female (default 0)

library(here)

args    <- commandArgs(trailingOnly = TRUE)
country <- if (length(args) >= 1) args[1] else "Germany"
site    <- if (length(args) >= 2) args[2] else "all"
sex     <- if (length(args) >= 3) as.integer(args[3]) else 0L

rmarkdown::render(here("code", "09_lexis.Rmd"), params = list(country = country, site = site, sex = sex),
                  output_file = sprintf("lexis_%s_%s_%s.html", tolower(gsub("\\W+", "_", country)), site, c("MF", "M", "F")[sex + 1]),
                  output_dir = here("report"), intermediates_dir = tempdir(), envir = new.env())
