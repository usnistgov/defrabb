#!/usr/bin/env Rscript
################################################################################
####
#### Plot DeFrABB evaluation comparisons.
####
#### Consumes the tidy TSV written by scripts/compare_evaluations.py and renders
#### comparison plots with ggplot2. Plotting lives here (R) rather than in the
#### Python utility so it shares the project's tidyverse/quarto reporting stack.
####
#### Usage:
####   Rscript scripts/plot_evaluations.R <comparison.tsv> <out_dir>
####
#### Outputs (PNG) in <out_dir>:
####   metrics_by_analysis.png   precision/recall/F1 by analysis & variant type
####   delta_<metric>.png        deltas vs baseline (only if delta_* present)
####
#### Developed with assistance from Claude (Anthropic); reviewed by the
#### primary author.
####
################################################################################

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_evaluations.R <comparison.tsv> <out_dir>", call. = FALSE)
}
tsv_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_tsv(tsv_path, show_col_types = FALSE)
if (nrow(df) == 0) {
  stop("No rows in comparison table: ", tsv_path, call. = FALSE)
}

## --- Plot 1: precision / recall / F1 by analysis and variant type ------------
metrics_long <- df %>%
  select(analysis_id, variant_type, precision, recall, f1) %>%
  pivot_longer(
    c(precision, recall, f1),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c("precision", "recall", "f1")))

p_metrics <- ggplot(
  metrics_long,
  aes(x = value, y = analysis_id, fill = metric)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~variant_type, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    title = "Evaluation metrics by analysis",
    x = "value", y = NULL, fill = "metric"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(out_dir, "metrics_by_analysis.png"),
  p_metrics,
  width = 10,
  height = max(3, 0.4 * length(unique(df$analysis_id)) + 2),
  dpi = 150,
  limitsize = FALSE
)

## --- Plot 2: deltas vs baseline (only when compare ran with --baseline) ------
delta_cols <- grep("^delta_", names(df), value = TRUE)
if (length(delta_cols) > 0 && "status" %in% names(df)) {
  delta_long <- df %>%
    filter(status != "baseline") %>%
    select(analysis_id, variant_type, status, all_of(delta_cols)) %>%
    pivot_longer(
      all_of(delta_cols),
      names_to = "metric", values_to = "delta",
      names_prefix = "delta_"
    ) %>%
    filter(!is.na(delta))

  if (nrow(delta_long) > 0) {
    p_delta <- ggplot(
      delta_long,
      aes(x = delta, y = analysis_id, fill = status)
    ) +
      geom_col(width = 0.7) +
      facet_grid(variant_type ~ metric, scales = "free_y") +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      scale_fill_manual(
        values = c(
          improved = "#2c7fb8", regressed = "#d95f0e",
          same = "grey70", `n/a` = "grey90"
        )
      ) +
      labs(
        title = "Change vs baseline",
        x = "delta (analysis - baseline)", y = NULL, fill = "status"
      ) +
      theme_minimal(base_size = 11)

    ggsave(
      file.path(out_dir, "delta_metrics.png"),
      p_delta,
      width = 11,
      height = max(3, 0.4 * length(unique(delta_long$analysis_id)) + 2),
      dpi = 150,
      limitsize = FALSE
    )
  }
}

message("Wrote plots to ", out_dir)
