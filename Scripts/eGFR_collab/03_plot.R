#!/usr/bin/Rscript
# Plot SumStatTune + IndivTune R values per (method, source, target).
#
# Reads the *.sumstat_tune.tsv and *.indiv_tune.tsv outputs from
# 02_evaluate_pgs.R for one or more target populations and renders a
# facet_grid plot (rows = target, cols = source GWAS) like the figure
# example.
#
# Usage:
#   Rscript 03_plot.R --inputs "EUR=results_EUR,AFR=results_AFR" \
#                     --output figure.png
#
# --inputs   Comma-separated list of <label>=<prefix> pairs. Each prefix
#            is the --output value passed to 02_evaluate_pgs.R for that
#            target population (script reads <prefix>.sumstat_tune.tsv
#            and <prefix>.indiv_tune.tsv).
# --output   Output figure path (.png or .pdf).
# --title    Figure title.
# --width    Inches (default 10).
# --height   Inches (default 6).

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt_list <- list(
  make_option("--inputs", type = "character"),
  make_option("--output", type = "character"),
  make_option("--title",  type = "character", default = ""),
  make_option("--width",  type = "double", default = 10),
  make_option("--height", type = "double", default = 6)
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$inputs), !is.null(opt$output))

pairs <- strsplit(opt$inputs, ",")[[1]]
all_rows <- list()
for (p in pairs) {
  kv <- strsplit(p, "=")[[1]]
  if (length(kv) != 2) stop("Bad --inputs entry: ", p)
  label <- kv[1]; prefix <- kv[2]
  for (suffix in c(".sumstat_tune.tsv", ".indiv_tune.tsv")) {
    f <- paste0(prefix, suffix)
    if (!file.exists(f)) { message("Missing (skipping): ", f); next }
    d <- fread(f)
    if (!"target" %in% names(d)) d[, target := label]
    all_rows[[length(all_rows) + 1L]] <- d
  }
}
dt <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
if (nrow(dt) == 0) stop("No data loaded from --inputs.")

# Order methods on x-axis: single-source-likely-on-left + prscsx last
method_order <- c("ptclump", "lassosum", "ldpred2", "megaprs", "prscs",
                  "quickprs", "sbayesrc", "prscsx")
dt[, method := factor(method, levels = intersect(method_order, unique(method)))]
dt[, source := factor(source, levels = c("AFR", "EUR", "EUR+AFR"))]
dt[, target := factor(target)]
dt[, tune   := factor(tune,   levels = c("SumStatTune", "IndivTune"))]

p <- ggplot(dt, aes(x = method, y = R, colour = tune, shape = tune)) +
  geom_errorbar(aes(ymin = R - R_se, ymax = R + R_se), width = 0.25,
                position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_grid(rows = vars(target), cols = vars(source),
             labeller = labeller(target = function(x) paste0(x, " target"),
                                 source = function(x) paste0(x, " GWAS"))) +
  scale_colour_manual(values = c(SumStatTune = "#1FB7C1", IndivTune = "#D9706C")) +
  scale_shape_manual(values = c(SumStatTune = 16, IndivTune = 18)) +
  labs(x = NULL, y = "R (predicted vs observed phenotype)",
       colour = NULL, shape = NULL, title = opt$title) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave(opt$output, p, width = opt$width, height = opt$height, dpi = 150)
message("Wrote: ", opt$output)
