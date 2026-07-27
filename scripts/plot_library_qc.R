#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_library_qc.R RESULT_DIR OUTPUT_DIR")
}

required_packages <- c(
  "ggplot2",
  "scales",
  "patchwork",
  "viridisLite",
  "ragg"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the plotting environment before running this script."
  )
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
  stop(
    "Count table not found. Expected: ",
    csv_count_file,
    " or ",
    tsv_count_file
  )
}

required_columns <- c(
  "variant_ids",
  "length",
  "gc_percent",
  "exact_count",
  "near_count",
  "total_count"
)
missing_columns <- setdiff(required_columns, colnames(counts))
if (length(missing_columns) > 0) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
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
    stop("Non-numeric or missing value in column: ", column)
  }
}
if (any(counts$total_count < 0)) {
  stop("total_count contains a negative value.")
}

format_integer <- function(x) {
  trimws(formatC(
    as.numeric(x),
    format = "f",
    digits = 0,
    big.mark = ","
  ))
}

result_name <- basename(result_dir)
is_one_percent <- grepl(
  "1pct|1percent|1_percent",
  result_name,
  ignore.case = TRUE
)
is_full_dataset <- grepl("full", result_name, ignore.case = TRUE)

dataset_context <- if (is_one_percent) {
  "1% paired-end subsample"
} else if (is_full_dataset) {
  "Full FASTQ dataset"
} else {
  paste("Result:", result_name)
}

if (is_one_percent) {
  coverage_levels <- c(
    "Dropout",
    "1-9 reads",
    "10-99 reads",
    "Reads >= 100"
  )
  coverage_label_math <- c(
    "plain('Dropout')",
    "plain('1-9 reads')",
    "plain('10-99 reads')",
    "plain(Reads) >= 100"
  )
  coverage_breaks <- c(-Inf, 0, 9, 99, Inf)
  coverage_thresholds <- c(10, 100)
} else {
  coverage_levels <- c(
    "0 reads",
    "1 <= Reads < 10^3",
    "10^3 <= Reads < 10^4",
    "10^4 <= Reads < 10^5",
    "Reads >= 10^5"
  )
  coverage_label_math <- c(
    "plain('0 reads')",
    "1 <= plain(Reads) * ' < ' * 10^3",
    "10^3 <= plain(Reads) * ' < ' * 10^4",
    "10^4 <= plain(Reads) * ' < ' * 10^5",
    "plain(Reads) >= 10^5"
  )
  coverage_breaks <- c(-Inf, 0, 999, 9999, 99999, Inf)
  coverage_thresholds <- c(1000, 10000, 100000)
}

coverage_axis_labeler <- function(x) {
  parse(text = coverage_label_math[match(x, coverage_levels)])
}

coverage_class <- cut(
  counts$total_count,
  breaks = coverage_breaks,
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
coverage_summary$coverage_class <- factor(
  coverage_summary$coverage_class,
  levels = rev(coverage_levels)
)
coverage_summary$label <- sprintf(
  "%s  (%.1f%%)",
  format_integer(coverage_summary$variants),
  coverage_summary$percent
)

write.csv(
  transform(
    coverage_summary,
    coverage_class = as.character(coverage_class)
  )[, c("coverage_class", "variants", "percent")],
  file.path(output_dir, "coverage_class_summary.csv"),
  quote = FALSE,
  row.names = FALSE
)

detected <- sum(counts$total_count > 0)
coverage_threshold_counts <- vapply(
  coverage_thresholds,
  function(threshold) sum(counts$total_count >= threshold),
  integer(1)
)
total_assigned <- sum(counts$total_count)
median_count <- median(counts$total_count)
mean_count <- mean(counts$total_count)
count_cv <- if (mean_count > 0) {
  stats::sd(counts$total_count) / mean_count
} else {
  NA_real_
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (
    sum(keep) < 3 ||
      length(unique(x[keep])) < 2 ||
      length(unique(y[keep])) < 2
  ) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = "spearman"))
}

length_rho <- safe_spearman(counts$length, counts$total_count)
gc_rho <- safe_spearman(counts$gc_percent, counts$total_count)

threshold_statistics <- do.call(
  rbind,
  lapply(seq_along(coverage_thresholds), function(index) {
    threshold <- coverage_thresholds[[index]]
    threshold_count <- coverage_threshold_counts[[index]]
    data.frame(
      metric = c(
        paste0("variants_at_least_", threshold, "_reads"),
        paste0("variants_at_least_", threshold, "_reads_percent")
      ),
      value = c(
        threshold_count,
        100 * threshold_count / nrow(counts)
      ),
      stringsAsFactors = FALSE
    )
  })
)

plot_statistics <- rbind(
  data.frame(
    metric = c(
      "reference_variants",
      "detected_variants",
      "detected_percent"
    ),
    value = c(
      nrow(counts),
      detected,
      100 * detected / nrow(counts)
    ),
    stringsAsFactors = FALSE
  ),
  threshold_statistics,
  data.frame(
    metric = c(
      "total_assigned_reads",
      "median_reads_per_variant",
      "mean_reads_per_variant",
      "count_cv",
      "length_count_spearman_rho",
      "gc_count_spearman_rho",
      "log_scale_pseudocount"
    ),
    value = c(
      total_assigned,
      median_count,
      mean_count,
      count_cv,
      length_rho,
      gc_rho,
      1
    ),
    stringsAsFactors = FALSE
  )
)
write.csv(
  plot_statistics,
  file.path(output_dir, "plot_statistics.csv"),
  quote = FALSE,
  row.names = FALSE
)

colors <- c(
  dropout = "#D55E00",
  low = "#E69F00",
  medium = "#56B4E9",
  high = "#0072B2",
  teal = "#009E73",
  purple = "#7B61A8",
  ink = "#24323D",
  muted = "#667580",
  grid = "#DDE5EA"
)

coverage_palette <- if (is_one_percent) {
  setNames(
    c(
      colors[["dropout"]],
      colors[["low"]],
      colors[["medium"]],
      colors[["high"]]
    ),
    coverage_levels
  )
} else {
  setNames(
    c(
      colors[["dropout"]],
      colors[["low"]],
      colors[["medium"]],
      colors[["high"]],
      colors[["teal"]]
    ),
    coverage_levels
  )
}
rank_line_colors <- if (is_one_percent) {
  c(colors[["low"]], colors[["dropout"]])
} else {
  c(colors[["low"]], colors[["medium"]], colors[["dropout"]])
}
rank_line_labels <- if (is_one_percent) {
  c("10~plain(reads)", "100~plain(reads)")
} else {
  c(
    "10^3~plain(reads)",
    "10^4~plain(reads)",
    "10^5~plain(reads)"
  )
}

theme_library_qc <- function() {
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = colors[["ink"]],
        face = "bold",
        size = 16,
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        color = colors[["muted"]],
        size = 10.5,
        lineheight = 1.15,
        margin = ggplot2::margin(b = 12)
      ),
      plot.caption = ggplot2::element_text(
        color = colors[["muted"]],
        size = 9,
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      axis.title = ggplot2::element_text(
        color = colors[["ink"]],
        face = "bold",
        size = 11
      ),
      axis.text = ggplot2::element_text(
        color = colors[["ink"]],
        size = 10
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        color = colors[["grid"]],
        linewidth = 0.35
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(18, 30, 16, 18)
    )
}

count_axis_values <- function(max_count) {
  max_count <- max(0, max_count)
  if (max_count == 0) {
    return(c(0, 1))
  }
  max_power <- max(1, ceiling(log10(max_count)))
  values <- sort(unique(c(0, 1, 10^(seq_len(max_power)))))
  values[values <= max_count | values %in% c(0, 1)]
}

count_breaks <- count_axis_values(max(counts$total_count))
count_break_positions <- count_breaks + 1
count_break_labels <- format_integer(count_breaks)

coverage_max <- max(coverage_summary$variants)
coverage_padding <- max(1, coverage_max * 0.025)
coverage_summary$label_x <- coverage_summary$variants + coverage_padding
coverage_limit <- max(
  1,
  max(coverage_summary$label_x) + max(12, coverage_max * 0.28)
)

p_coverage <- ggplot2::ggplot(
  coverage_summary,
  ggplot2::aes(
    x = variants,
    y = coverage_class,
    fill = coverage_class
  )
) +
  ggplot2::geom_col(width = 0.66, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(x = label_x, label = label),
    hjust = 0,
    color = colors[["ink"]],
    fontface = "bold",
    size = 4.0
  ) +
  ggplot2::scale_fill_manual(
    values = coverage_palette
  ) +
  ggplot2::scale_y_discrete(
    labels = coverage_axis_labeler
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, coverage_limit),
    labels = format_integer,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "5'UTR library coverage",
    subtitle = sprintf(
      "%s | %s reference variants | %s (%.1f%%) have at least %s assigned reads",
      dataset_context,
      format_integer(nrow(counts)),
      format_integer(coverage_threshold_counts[[1]]),
      100 * coverage_threshold_counts[[1]] / nrow(counts),
      format_integer(coverage_thresholds[[1]])
    ),
    x = "Number of reference variants",
    y = NULL
  ) +
  theme_library_qc()

ranked <- counts[order(counts$total_count, decreasing = TRUE), , drop = FALSE]
ranked$rank <- seq_len(nrow(ranked))
ranked$plot_count <- ranked$total_count + 1

p_rank <- ggplot2::ggplot(
  ranked,
  ggplot2::aes(x = rank, y = plot_count)
) +
  ggplot2::geom_line(
    color = colors[["high"]],
    linewidth = 0.85
  ) +
  ggplot2::geom_hline(
    yintercept = coverage_thresholds + 1,
    color = rank_line_colors,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  ggplot2::annotate(
    "text",
    x = nrow(ranked) * 0.98,
    y = coverage_thresholds + 1,
    label = rank_line_labels,
    parse = TRUE,
    hjust = 1,
    vjust = -0.45,
    color = rank_line_colors,
    size = 3.3
  ) +
  ggplot2::scale_y_log10(
    breaks = count_break_positions,
    labels = count_break_labels,
    expand = ggplot2::expansion(mult = c(0.03, 0.12))
  ) +
  ggplot2::scale_x_continuous(
    labels = format_integer,
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::labs(
    title = "Rank-abundance curve",
    subtitle = sprintf(
      "%s | Median %s reads/variant | Total %s assigned reads",
      dataset_context,
      format_integer(median_count),
      format_integer(total_assigned)
    ),
    caption = "Tick labels show actual read counts; spacing is logarithmic. Zero-count variants are plotted using total_count + 1.",
    x = "Variant rank",
    y = "Assigned reads per variant"
  ) +
  theme_library_qc()

counts$log10_count_plus_one <- log10(counts$total_count + 1)
distribution_break_values <- count_axis_values(max(counts$total_count))

p_distribution <- ggplot2::ggplot(
  counts,
  ggplot2::aes(x = log10_count_plus_one)
) +
  ggplot2::geom_histogram(
    bins = 35,
    boundary = 0,
    closed = "left",
    fill = colors[["teal"]],
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_x_continuous(
    breaks = log10(distribution_break_values + 1),
    labels = format_integer(distribution_break_values),
    expand = ggplot2::expansion(mult = c(0.01, 0.03))
  ) +
  ggplot2::scale_y_continuous(
    labels = format_integer,
    expand = ggplot2::expansion(mult = c(0, 0.08))
  ) +
  ggplot2::labs(
    title = "Variant read-count distribution",
    subtitle = sprintf(
      "%s | Detected %s of %s variants (%.1f%%)",
      dataset_context,
      format_integer(detected),
      format_integer(nrow(counts)),
      100 * detected / nrow(counts)
    ),
    caption = "The x-axis uses log10(total_count + 1), so dropout variants remain visible at 0 reads.",
    x = "Assigned reads",
    y = "Number of variants"
  ) +
  theme_library_qc()

counts$plot_count <- counts$total_count + 1
relationship_plot <- function(
  x_column,
  x_label,
  point_color,
  rho,
  panel_title
) {
  rho_label <- if (is.finite(rho)) {
    sprintf("Spearman rho = %.3f", rho)
  } else {
    "Spearman rho = NA"
  }

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data[[x_column]], y = plot_count)
  ) +
    ggplot2::geom_point(
      color = point_color,
      alpha = 0.46,
      size = 1.55
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = colors[["ink"]],
      fill = point_color,
      alpha = 0.12,
      linewidth = 0.7
    ) +
    ggplot2::scale_y_log10(
      breaks = count_break_positions,
      labels = count_break_labels,
      expand = ggplot2::expansion(mult = c(0.03, 0.16))
    ) +
    ggplot2::labs(
      title = panel_title,
      subtitle = rho_label,
      x = x_label,
      y = "Assigned reads per variant"
    ) +
    theme_library_qc() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14),
      plot.margin = ggplot2::margin(16, 18, 14, 16)
    )
}

p_length <- relationship_plot(
  "length",
  "5'UTR length (nt)",
  colors[["high"]],
  length_rho,
  "Length vs abundance"
)
p_gc <- relationship_plot(
  "gc_percent",
  "GC content (%)",
  colors[["purple"]],
  gc_rho,
  "GC content vs abundance"
)

p_relationship <- (
  p_length + p_gc +
    patchwork::plot_layout(ncol = 2)
) +
  patchwork::plot_annotation(
    title = "Sequence properties versus abundance",
    subtitle = dataset_context,
    caption = "Tick labels show actual read counts; spacing is logarithmic. Zero-count variants are retained using total_count + 1.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = colors[["ink"]],
        face = "bold",
        size = 16
      ),
      plot.subtitle = ggplot2::element_text(
        color = colors[["muted"]],
        size = 10.5
      ),
      plot.caption = ggplot2::element_text(
        color = colors[["muted"]],
        size = 9,
        hjust = 0
      ),
      plot.margin = ggplot2::margin(14, 20, 12, 16)
    )
  )

save_png <- function(filename, plot, width, height) {
  ragg::agg_png(
    filename = file.path(output_dir, filename),
    width = width,
    height = height,
    units = "in",
    res = 320,
    background = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

save_png("01_coverage_classes.png", p_coverage, 10.5, 6.2)
save_png("02_rank_abundance.png", p_rank, 10.5, 6.2)
save_png("03_count_distribution.png", p_distribution, 10.5, 6.2)
save_png("04_length_gc_vs_count.png", p_relationship, 13.0, 6.3)

grDevices::pdf(
  file.path(output_dir, "library_qc_figures.pdf"),
  width = 10.5,
  height = 6.5,
  onefile = TRUE,
  family = "Helvetica",
  paper = "special"
)
print(p_coverage)
print(p_rank)
print(p_distribution)
print(p_relationship)
grDevices::dev.off()

message(
  "R figures written to: ",
  normalizePath(output_dir, mustWork = TRUE)
)
