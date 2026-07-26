#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_library_qc.R RESULT_DIR OUTPUT_DIR")
}

result_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(
  result_dir,
  "combined",
  "ALL_ASSIGNED_variant_counts.tsv"
)
if (!file.exists(count_file)) {
  stop(paste("Count table not found:", count_file))
}

counts <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required <- c("variant_ids", "length", "gc_percent", "exact_count", "near_count", "total_count")
missing <- setdiff(required, colnames(counts))
if (length(missing) > 0) {
  stop(paste("Missing required column(s):", paste(missing, collapse = ", ")))
}

coverage_class <- cut(
  counts$total_count,
  breaks = c(-Inf, 0, 9, 99, Inf),
  labels = c("Dropout", "1-9", "10-99", ">=100"),
  right = TRUE
)
coverage_levels <- c("Dropout", "1-9", "10-99", ">=100")
coverage_table <- table(factor(coverage_class, levels = coverage_levels))
coverage_summary <- data.frame(
  coverage_class = coverage_levels,
  variants = as.integer(coverage_table),
  percent = 100 * as.integer(coverage_table) / nrow(counts)
)
write.table(
  coverage_summary,
  file.path(output_dir, "coverage_class_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

png(
  file.path(output_dir, "01_coverage_classes.png"),
  width = 1600,
  height = 1100,
  res = 180
)
colors <- c("#C0392B", "#E67E22", "#F1C40F", "#2E86C1")
bars <- barplot(
  coverage_summary$variants,
  names.arg = coverage_summary$coverage_class,
  col = colors,
  border = NA,
  ylab = "Number of variants",
  xlab = "Assigned read count",
  main = "5'UTR library coverage classes"
)
text(
  bars,
  coverage_summary$variants,
  labels = sprintf(
    "%s\n(%.1f%%)",
    coverage_summary$variants,
    coverage_summary$percent
  ),
  pos = 3,
  cex = 0.9
)
dev.off()

ranked <- sort(counts$total_count, decreasing = TRUE)
png(
  file.path(output_dir, "02_rank_abundance.png"),
  width = 1600,
  height = 1100,
  res = 180
)
plot(
  seq_along(ranked),
  ranked + 1,
  type = "l",
  log = "y",
  lwd = 2,
  col = "#21618C",
  xlab = "Variant rank",
  ylab = "Assigned reads + 1 (log scale)",
  main = "Rank-abundance curve"
)
grid(col = "grey85")
dev.off()

png(
  file.path(output_dir, "03_count_distribution.png"),
  width = 1600,
  height = 1100,
  res = 180
)
hist(
  log10(counts$total_count + 1),
  breaks = 40,
  col = "#5DADE2",
  border = "white",
  xlab = "log10(assigned reads + 1)",
  ylab = "Number of variants",
  main = "Variant count distribution"
)
dev.off()

png(
  file.path(output_dir, "04_length_gc_vs_count.png"),
  width = 1800,
  height = 900,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
plot(
  counts$length,
  log10(counts$total_count + 1),
  pch = 16,
  col = rgb(33 / 255, 97 / 255, 140 / 255, 0.55),
  xlab = "5'UTR length (nt)",
  ylab = "log10(assigned reads + 1)",
  main = "Length vs abundance"
)
grid(col = "grey88")
plot(
  counts$gc_percent,
  log10(counts$total_count + 1),
  pch = 16,
  col = rgb(125 / 255, 60 / 255, 152 / 255, 0.55),
  xlab = "GC (%)",
  ylab = "log10(assigned reads + 1)",
  main = "GC vs abundance"
)
grid(col = "grey88")
dev.off()

message("R figures written to: ", output_dir)
