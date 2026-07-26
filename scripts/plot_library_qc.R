#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_library_qc.R RESULT_DIR OUTPUT_DIR")
}

result_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_count_file <- file.path(
  result_dir,
  "combined",
  "ALL_ASSIGNED_variant_counts.csv"
)
tsv_count_file <- file.path(
  result_dir,
  "combined",
  "ALL_ASSIGNED_variant_counts.tsv"
)
if (file.exists(csv_count_file)) {
  count_file <- csv_count_file
  counts <- read.csv(
    count_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else if (file.exists(tsv_count_file)) {
  count_file <- tsv_count_file
  counts <- read.delim(
    count_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  stop(paste(
    "Count table not found. Expected:",
    csv_count_file,
    "or",
    tsv_count_file
  ))
}

required <- c(
  "variant_ids",
  "length",
  "gc_percent",
  "exact_count",
  "near_count",
  "total_count"
)
missing <- setdiff(required, colnames(counts))
if (length(missing) > 0) {
  stop(paste("Missing required column(s):", paste(missing, collapse = ", ")))
}
if (nrow(counts) == 0) {
  stop("The count table contains no variants.")
}

numeric_columns <- c(
  "length",
  "gc_percent",
  "exact_count",
  "near_count",
  "total_count"
)
for (column in numeric_columns) {
  counts[[column]] <- suppressWarnings(as.numeric(counts[[column]]))
  if (any(!is.finite(counts[[column]]))) {
    stop(paste("Non-numeric or missing value in column:", column))
  }
}
if (any(counts$total_count < 0)) {
  stop("total_count contains a negative value.")
}

format_integer <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

coverage_levels <- c("Dropout", "1-9 reads", "10-99 reads", "\u2265100 reads")
coverage_class <- cut(
  counts$total_count,
  breaks = c(-Inf, 0, 9, 99, Inf),
  labels = coverage_levels,
  right = TRUE
)
coverage_table <- table(factor(coverage_class, levels = coverage_levels))
coverage_summary <- data.frame(
  coverage_class = coverage_levels,
  variants = as.integer(coverage_table),
  percent = 100 * as.integer(coverage_table) / nrow(counts),
  stringsAsFactors = FALSE
)
write.csv(
  coverage_summary,
  file.path(output_dir, "coverage_class_summary.csv"),
  quote = FALSE,
  row.names = FALSE
)

detected <- sum(counts$total_count > 0)
at_least_10 <- sum(counts$total_count >= 10)
at_least_100 <- sum(counts$total_count >= 100)
total_assigned <- sum(counts$total_count)
median_count <- median(counts$total_count)
mean_count <- mean(counts$total_count)
count_cv <- if (mean_count > 0) sd(counts$total_count) / mean_count else NA_real_

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

length_rho <- safe_spearman(counts$length, counts$total_count)
gc_rho <- safe_spearman(counts$gc_percent, counts$total_count)

plot_statistics <- data.frame(
  metric = c(
    "reference_variants",
    "detected_variants",
    "detected_percent",
    "variants_at_least_10_reads",
    "variants_at_least_10_reads_percent",
    "variants_at_least_100_reads",
    "variants_at_least_100_reads_percent",
    "total_assigned_reads",
    "median_reads_per_variant",
    "mean_reads_per_variant",
    "count_cv",
    "length_count_spearman_rho",
    "gc_count_spearman_rho"
  ),
  value = c(
    nrow(counts),
    detected,
    100 * detected / nrow(counts),
    at_least_10,
    100 * at_least_10 / nrow(counts),
    at_least_100,
    100 * at_least_100 / nrow(counts),
    total_assigned,
    median_count,
    mean_count,
    count_cv,
    length_rho,
    gc_rho
  ),
  stringsAsFactors = FALSE
)
write.csv(
  plot_statistics,
  file.path(output_dir, "plot_statistics.csv"),
  quote = FALSE,
  row.names = FALSE
)

palette <- c(
  navy = "#173F5F",
  blue = "#20639B",
  teal = "#3CAEA3",
  gold = "#F6D55C",
  orange = "#ED8B3A",
  red = "#C94C4C",
  purple = "#7C5AA6",
  ink = "#263238",
  grid = "#DCE3E8",
  pale = "#F4F7F9"
)

set_plot_theme <- function(margins = c(5.0, 5.2, 4.8, 1.8)) {
  par(
    mar = margins,
    mgp = c(3.0, 0.8, 0),
    tcl = -0.25,
    las = 1,
    bty = "n",
    family = "sans",
    fg = palette[["ink"]],
    col.axis = palette[["ink"]],
    col.lab = palette[["ink"]],
    col.main = palette[["ink"]]
  )
}

add_title <- function(main_title, subtitle = "") {
  title(main = main_title, adj = 0, font.main = 2, cex.main = 1.15, line = 2.2)
  if (nzchar(subtitle)) {
    mtext(subtitle, side = 3, adj = 0, line = 0.75, cex = 0.82, col = "#5F6B73")
  }
}

plot_coverage <- function() {
  set_plot_theme(c(4.8, 7.2, 5.2, 2.2))
  values <- rev(coverage_summary$variants)
  labels <- rev(coverage_summary$coverage_class)
  percentages <- rev(coverage_summary$percent)
  bar_colors <- rev(c(
    palette[["red"]],
    palette[["orange"]],
    palette[["gold"]],
    palette[["blue"]]
  ))
  x_max <- max(1, max(values))
  label_room <- max(1, x_max * 0.30)
  positions <- barplot(
    values,
    names.arg = labels,
    horiz = TRUE,
    col = bar_colors,
    border = NA,
    xlim = c(0, x_max + label_room),
    xlab = "Number of reference variants",
    axes = FALSE,
    cex.names = 0.95
  )
  abline(v = pretty(c(0, x_max)), col = palette[["grid"]], lwd = 0.8)
  axis(1, at = pretty(c(0, x_max)), labels = format_integer(pretty(c(0, x_max))))
  text(
    x = values + max(1, x_max * 0.02),
    y = positions,
    labels = sprintf("%s  (%.1f%%)", format_integer(values), percentages),
    adj = 0,
    cex = 0.92,
    font = 2,
    xpd = FALSE,
    col = palette[["ink"]]
  )
  add_title(
    "5'UTR library coverage",
    sprintf(
      "%s variants; %s (%.1f%%) have at least 10 assigned reads",
      format_integer(nrow(counts)),
      format_integer(at_least_10),
      100 * at_least_10 / nrow(counts)
    )
  )
}

plot_rank_abundance <- function() {
  set_plot_theme()
  ranked <- sort(counts$total_count, decreasing = TRUE)
  x <- seq_along(ranked)
  y <- ranked + 1
  plot(
    x,
    y,
    type = "l",
    log = "y",
    lwd = 2.2,
    col = palette[["blue"]],
    xlab = "Variant rank",
    ylab = "Assigned reads + 1 (log scale)",
    axes = FALSE
  )
  x_ticks <- pretty(range(x))
  y_ticks <- 10^(0:ceiling(log10(max(y))))
  abline(h = y_ticks, col = palette[["grid"]], lwd = 0.8)
  axis(1, at = x_ticks, labels = format_integer(x_ticks))
  axis(2, at = y_ticks, labels = format_integer(y_ticks))
  abline(h = c(11, 101), lty = 3, lwd = 1.0, col = c(palette[["orange"]], palette[["red"]]))
  legend(
    "topright",
    legend = c("10 reads", "100 reads"),
    lty = 3,
    lwd = 1.2,
    col = c(palette[["orange"]], palette[["red"]]),
    bty = "n",
    cex = 0.82
  )
  add_title(
    "Rank-abundance curve",
    sprintf(
      "Median %s reads/variant; total %s assigned reads",
      format_integer(median_count),
      format_integer(total_assigned)
    )
  )
}

plot_distribution <- function() {
  set_plot_theme()
  transformed <- log10(counts$total_count + 1)
  breaks <- unique(pretty(range(transformed), n = 35))
  if (length(breaks) < 2) {
    breaks <- c(transformed[[1]] - 0.5, transformed[[1]] + 0.5)
  }
  histogram <- hist(transformed, breaks = breaks, plot = FALSE)
  plot(
    histogram,
    col = palette[["teal"]],
    border = "white",
    freq = TRUE,
    main = "",
    xlab = expression(log[10] * "(assigned reads + 1)"),
    ylab = "Number of variants",
    axes = FALSE
  )
  abline(h = pretty(c(0, max(histogram$counts))), col = palette[["grid"]], lwd = 0.8)
  axis(1)
  axis(2, at = pretty(c(0, max(histogram$counts))), labels = format_integer(pretty(c(0, max(histogram$counts)))))
  add_title(
    "Variant read-count distribution",
    sprintf(
      "Detected %s of %s variants (%.1f%%)",
      format_integer(detected),
      format_integer(nrow(counts)),
      100 * detected / nrow(counts)
    )
  )
}

plot_relationship <- function(x, x_label, color, rho, panel_title) {
  set_plot_theme()
  y <- log10(counts$total_count + 1)
  plot(
    x,
    y,
    pch = 16,
    cex = 0.72,
    col = grDevices::adjustcolor(color, alpha.f = 0.50),
    xlab = x_label,
    ylab = expression(log[10] * "(assigned reads + 1)"),
    axes = FALSE
  )
  abline(
    h = pretty(range(y)),
    v = pretty(range(x)),
    col = palette[["grid"]],
    lwd = 0.7
  )
  axis(1)
  axis(2)
  rho_label <- if (is.finite(rho)) sprintf("Spearman rho = %.3f", rho) else "Spearman rho = NA"
  add_title(panel_title, rho_label)
}

write_png <- function(filename, width, height, plot_function) {
  png(
    file.path(output_dir, filename),
    width = width,
    height = height,
    res = 180,
    bg = "white"
  )
  on.exit(dev.off())
  plot_function()
}

write_png("01_coverage_classes.png", 1800, 1100, plot_coverage)
write_png("02_rank_abundance.png", 1800, 1100, plot_rank_abundance)
write_png("03_count_distribution.png", 1800, 1100, plot_distribution)

png(
  file.path(output_dir, "04_length_gc_vs_count.png"),
  width = 2400,
  height = 1100,
  res = 180,
  bg = "white"
)
par(mfrow = c(1, 2))
plot_relationship(
  counts$length,
  "5'UTR length (nt)",
  palette[["blue"]],
  length_rho,
  "Length vs abundance"
)
plot_relationship(
  counts$gc_percent,
  "GC content (%)",
  palette[["purple"]],
  gc_rho,
  "GC content vs abundance"
)
dev.off()

pdf(
  file.path(output_dir, "library_qc_figures.pdf"),
  width = 10,
  height = 6.5,
  onefile = TRUE,
  family = "Helvetica",
  paper = "special"
)
plot_coverage()
plot_rank_abundance()
plot_distribution()
par(mfrow = c(1, 2))
plot_relationship(
  counts$length,
  "5'UTR length (nt)",
  palette[["blue"]],
  length_rho,
  "Length vs abundance"
)
plot_relationship(
  counts$gc_percent,
  "GC content (%)",
  palette[["purple"]],
  gc_rho,
  "GC content vs abundance"
)
dev.off()

message("R figures written to: ", normalizePath(output_dir, mustWork = TRUE))
