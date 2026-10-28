#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Trace plots of the joint log-likelihood for the beta and kappa sensitivity fits
#
# Plots logProb_Z1 + logProb_Z2 after burn-in (iterations 501 onwards), every
# tenth iteration, one panel per grid value, in the style of the trace plots in
# scripts/application/03_fit_models.R. The panel titles do not show R-hat; it is
# reported in the supplement table.
#
# Run from the project root after the "summarise" step of
# scripts/sensitivity/sensitivity_beta.R and sensitivity_kappa.R.
# Reads : results/sensitivity/{beta,kappa}/summary/tables/logProb_trace.rds
# Writes: results/figures/sens_trace_beta.png, results/figures/sens_trace_kappa.png
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
suppressMessages({library(tidyverse); library(posterior)})
burn <- 500
thin <- 10

default_theme <- theme(legend.position = "bottom",
                       legend.key.size = unit(0.5, "cm"),
                       axis.text = element_text(size = 8),
                       axis.title = element_text(size = 10),
                       legend.title = element_text(size = 10),
                       legend.spacing.x = unit(0.02, "cm"),
                       legend.margin = margin(0, 0, 0, 0),
                       plot.margin = margin(0, 10, 0, 0),
                       strip.background = element_blank(),
                       strip.text = element_text(size = 8, margin = margin(0, 0, 2, 0)),
                       legend.text = element_text(size = 9))

trace_fig <- function(study, symbol) {
  tr <- readRDS(file.path("results/sensitivity", study, "summary/tables/logProb_trace.rds")) %>%
    filter(iteration > burn) %>%
    mutate(loglik = logProb_Z1 + logProb_Z2)
  # R-hat on all post-burn-in iterations, printed for checking against the table
  rh <- tr %>% group_by(kappa_or_beta = label) %>%
    summarise(rhat = rhat(sapply(split(loglik, chain), identity)), .groups = "drop")
  print(as.data.frame(rh), digits = 3)
  tr %>%
    filter((iteration - burn) %% thin == 0) %>%
    mutate(value = as.numeric(label),
           value_txt = if_else(value == round(value, 1), sprintf("%.1f", value), sprintf("%.2f", value)),
           value_txt = if_else(value == 1.238, "1.238", value_txt),
           panel = fct_reorder(sprintf('%s == "%s"', symbol, value_txt), value)) %>%
    ggplot() +
    geom_line(aes(x = iteration, y = loglik, col = factor(chain)), alpha = 0.6, linewidth = 0.2) +
    facet_wrap(vars(panel), scales = "free_y", ncol = 3, labeller = label_parsed) +
    scale_x_continuous(breaks = \(lim) seq(500, lim[2], by = if (lim[2] > 3000) 1500 else 500),
                       labels = \(x) format(x, big.mark = ",")) +
    scale_y_continuous(labels = \(x) format(x, big.mark = ",")) +
    labs(x = "Iteration", y = "Joint log-likelihood", col = "Chain") +
    theme_bw() + default_theme
}

dir.create("results/figures", showWarnings = FALSE)
g_beta <- trace_fig("beta", "beta")
ggsave("results/figures/sens_trace_beta.png", g_beta, width = 6.5, height = 5.5, dpi = 300)
g_kappa <- trace_fig("kappa", "kappa")
ggsave("results/figures/sens_trace_kappa.png", g_kappa, width = 6.5, height = 4, dpi = 300)
message("Saved results/figures/sens_trace_beta.png and sens_trace_kappa.png")
