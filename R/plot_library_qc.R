#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript R/plot_library_qc.R <reference_counts.csv> <output_prefix>")
}

input_csv <- args[[1]]
out_prefix <- args[[2]]

df <- read_csv(input_csv, show_col_types = FALSE) %>%
  mutate(
    status = if_else(raw_count > 0, "Detected", "Dropout"),
    rank = rank(-raw_count, ties.method = "first"),
    log10_count = log10(raw_count + 1)
  )

p1 <- ggplot(df, aes(x = log10_count)) +
  geom_histogram(bins = 30, boundary = 0) +
  labs(
    title = "Reference count distribution",
    x = "log10(count + 1)",
    y = "Number of reference sequences"
  ) +
  theme_classic(base_size = 12)

p2 <- ggplot(df, aes(x = rank, y = raw_count, shape = status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_y_continuous(trans = "log1p", labels = label_number()) +
  labs(
    title = "Ranked library abundance",
    x = "Abundance rank",
    y = "Raw count",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

label_df <- df %>%
  arrange(desc(raw_count)) %>%
  slice_head(n = 3) %>%
  bind_rows(df %>% filter(raw_count == 0)) %>%
  distinct()

p3 <- ggplot(df, aes(x = rank, y = frequency, shape = status)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(
    data = label_df,
    aes(label = .data[[names(df)[1]]]),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 0.01)) +
  labs(
    title = "Variant frequency and dropout",
    x = "Abundance rank",
    y = "Library frequency",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

combined <- p1 / (p2 | p3) +
  plot_annotation(title = "NGS Library QC")

ggsave(paste0(out_prefix, ".ggplot_qc.png"), combined, width = 12, height = 9, dpi = 300)
ggsave(paste0(out_prefix, ".ggplot_qc.pdf"), combined, width = 12, height = 9)

message("Created: ", out_prefix, ".ggplot_qc.png")
message("Created: ", out_prefix, ".ggplot_qc.pdf")
