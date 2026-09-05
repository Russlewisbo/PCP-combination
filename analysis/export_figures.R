# Export high-resolution (300 dpi) PNG copies of every figure on the site.
# Each block mirrors the plotting code from the corresponding .qmd page, so the PNGs
# match what is rendered online. (The two PRISMA diagrams have their own scripts:
# analysis/prisma_pcp.R and analysis/prisma2020_flow.R.)
# Usage:  source("analysis/export_figures.R")
source("analysis/meta_helpers.R")
library(coda)

ma <- load_ma()
dir.create("figures", showWarnings = FALSE)
save_fig <- function(plot, name, width = 8, height = 5)
  ggsave(file.path("figures", name), plot, width = width, height = height,
         dpi = 300, bg = "white")

## ---- Frequentist (frequentist.qmd) ----
re <- rma(yi, vi, data = ma, method = "REML")
mh <- rma.mh(ai = e_comb, n1i = n_comb, ci = e_ref, n2i = n_ref, data = ma, measure = "OR")

est  <- ma |> transmute(label = study, kind = "study",
                        OR = exp(yi), lo = exp(yi - 1.96 * sei), hi = exp(yi + 1.96 * sei))
pool <- tibble(label = c("Random-effects (REML)", "Mantel-Haenszel (fixed)"), kind = "pooled",
               OR = c(exp(re$b[1]), exp(mh$b[1])),
               lo = c(exp(re$ci.lb), exp(mh$ci.lb)),
               hi = c(exp(re$ci.ub), exp(mh$ci.ub)))
p_freq_forest <- bind_rows(est, pool) |>
  mutate(label = factor(label, levels = rev(c(ma$study,
           "Random-effects (REML)", "Mantel-Haenszel (fixed)")))) |>
  ggplot(aes(OR, label, color = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = kind == "pooled"), linewidth = 0.6) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4)) +
  scale_color_manual(values = c(study = jama_blue_grey[["medium"]],
                                pooled = jama_blue_grey[["accent"]]), guide = "none") +
  scale_shape_manual(values = c(16, 18), guide = "none") +
  labs(x = "Odds ratio (log scale) \u2014 <1 favours combination", y = NULL,
       title = "Frequentist forest plot",
       subtitle = "Per-study estimates with random-effects and fixed-effect pooled odds ratios") +
  theme_jama()
save_fig(p_freq_forest, "freq_forest.png", 8, 5)

p_funnel <- ggplot(ma, aes(x = OR, y = sei)) +
  geom_vline(xintercept = exp(re$b[1]), linetype = "dashed", color = "grey50") +
  geom_point(size = 2.4, color = jama_blue_grey[["steel"]]) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  scale_y_reverse() +
  labs(x = "Odds ratio (log scale)", y = "Standard error", title = "Funnel plot",
       subtitle = "Dashed line = random-effects pooled estimate") +
  theme_jama()
save_fig(p_funnel, "freq_funnel.png", 7, 5)

## ---- Bayesian (bayesian.qmd) ----
bm_sens <- fit_bm(ma)

set.seed(7183)
ps <- bm_sens$rposterior(n = 6000)
theta_draws <- purrr::map_dfr(seq_len(nrow(ma)), function(i) {
  pd <- 1 / ma$sei[i]^2; pp <- 1 / ps[, "tau"]^2; v <- 1 / (pd + pp)
  m  <- v * (ma$yi[i] * pd + ps[, "mu"] * pp)
  tibble(label = ma$study[i], draw = exp(rnorm(nrow(ps), m, sqrt(v))))
})
pooled_draws <- tibble(label = "Pooled (\u03bc)", draw = exp(ps[, "mu"]))
alldraws <- bind_rows(theta_draws, pooled_draws) |>
  mutate(label = factor(label, levels = rev(c(ma$study, "Pooled (\u03bc)"))),
         is_pooled = label == "Pooled (\u03bc)")
obs <- ma |> transmute(label = factor(study, levels = levels(alldraws$label)), OR)
p_post_forest <- ggplot(alldraws, aes(draw, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_violin(aes(fill = is_pooled), orientation = "y", color = NA,
              scale = "width", alpha = 0.85) +
  geom_point(data = obs, aes(OR, label), shape = 21, fill = "white",
             color = jama_blue_grey[["dark"]], size = 2) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4)) +
  scale_fill_manual(values = c(`FALSE` = jama_blue_grey[["medium"]],
                               `TRUE`  = jama_blue_grey[["accent"]]), guide = "none") +
  labs(x = "Odds ratio (log scale) \u2014 <1 favours combination", y = NULL,
       title = "Posterior forest plot",
       subtitle = "Posterior densities of shrunken study effects and the pooled effect; points = observed ORs") +
  theme_jama()
save_fig(p_post_forest, "bayes_posterior_forest.png", 8, 5)

mu_grid <- seq(-2.2, 1.6, length.out = 500)
p_post_mu <- tibble(OR = exp(mu_grid), density = bm_sens$dposterior(mu = mu_grid)) |>
  ggplot(aes(OR, density)) +
  geom_area(fill = jama_blue_grey[["light"]], color = jama_blue_grey[["steel"]]) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = or_ci(bm_sens)[1], color = jama_blue_grey[["accent"]]) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4)) +
  labs(x = "Pooled odds ratio (log scale)", y = "Posterior density",
       title = "Posterior of the pooled effect",
       subtitle = "Solid line = posterior median; dashed = no effect") +
  theme_jama()
save_fig(p_post_mu, "bayes_posterior_mu.png", 7, 4.5)

tau_grid <- seq(0, 2, length.out = 500)
p_post_tau <- tibble(tau = tau_grid, density = bm_sens$dposterior(tau = tau_grid)) |>
  ggplot(aes(tau, density)) +
  geom_area(fill = jama_blue_grey[["light"]], color = jama_blue_grey[["steel"]]) +
  labs(x = "Between-study SD (\u03c4, log-OR scale)", y = "Posterior density",
       title = "Posterior of between-study heterogeneity (\u03c4)") +
  theme_jama()
save_fig(p_post_tau, "bayes_posterior_tau.png", 7, 4.5)

# Clinical response / cure
resp <- tribble(
  ~study,    ~e_comb, ~n_comb, ~e_ref, ~n_ref, ~confounded,
  "Jin2019",      24,      35,     54,     91, FALSE,
  "Qi2023",       33,      43,     11,     21, FALSE,
  "Li2024a",      20,      20,     12,     18, TRUE)
respE <- escalc(measure = "OR", ai = e_comb, n1i = n_comb, ci = e_ref, n2i = n_ref,
                data = resp, slab = study)
respE$sei <- sqrt(respE$vi); respE$OR <- exp(respE$yi)
bm_resp_clean <- fit_bm(filter(respE, !confounded))
bm_resp_all   <- fit_bm(respE)
est  <- respE |> transmute(label = study, kind = "study",
                           OR = exp(yi), lo = exp(yi - 1.96 * sei), hi = exp(yi + 1.96 * sei))
pool <- tibble(label = c("Pooled (k=2, clean)", "Pooled (k=3, incl. Li2024a)"), kind = "pooled",
               OR = c(or_ci(bm_resp_clean)[1], or_ci(bm_resp_all)[1]),
               lo = c(or_ci(bm_resp_clean)[2], or_ci(bm_resp_all)[2]),
               hi = c(or_ci(bm_resp_clean)[3], or_ci(bm_resp_all)[3]))
p_resp_forest <- bind_rows(est, pool) |>
  mutate(label = factor(label, levels = rev(c("Jin2019", "Qi2023", "Li2024a",
           "Pooled (k=2, clean)", "Pooled (k=3, incl. Li2024a)")))) |>
  ggplot(aes(OR, label, color = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi, shape = kind == "pooled"), linewidth = 0.6) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  scale_color_manual(values = c(study = jama_blue_grey[["medium"]],
                                pooled = jama_blue_grey[["accent"]]), guide = "none") +
  scale_shape_manual(values = c(16, 18), guide = "none") +
  labs(x = "Odds ratio (log scale) \u2014 >1 favours combination", y = NULL,
       title = "Clinical response / cure", subtitle = "Event = clinical response") +
  theme_jama()
save_fig(p_resp_forest, "bayes_response_forest.png", 8, 4.5)

## ---- Sensitivity (sensitivity.qmd) ----
loo <- purrr::map_dfr(ma$study, function(s_) {
  bm <- fit_bm(dplyr::filter(ma, study != s_))
  tibble(removed = paste0("\u2212 ", s_),
         OR = or_ci(bm)[1], lo = or_ci(bm)[2], hi = or_ci(bm)[3])
})
loo <- bind_rows(
  tibble(removed = "Full model (k=6)",
         OR = or_ci(bm_sens)[1], lo = or_ci(bm_sens)[2], hi = or_ci(bm_sens)[3]),
  loo)
p_loo <- loo |>
  mutate(kind = ifelse(removed == "Full model (k=6)", "full", "loo"),
         removed = factor(removed, levels = rev(removed))) |>
  ggplot(aes(OR, removed, color = kind)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), shape = 18, linewidth = 0.6) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2)) +
  scale_color_manual(values = c(full = jama_blue_grey[["accent"]],
                                loo = jama_blue_grey[["medium"]]), guide = "none") +
  labs(x = "Pooled odds ratio (log scale) \u2014 <1 favours combination", y = NULL,
       title = "Leave-one-out sensitivity", subtitle = "Top row = full model") +
  theme_jama()
save_fig(p_loo, "sensitivity_loo.png", 8, 4.5)

## ---- Risk of bias (risk-of-bias.qmd) ----
rob <- tribble(
  ~study,   ~D1,        ~D2,        ~D3,        ~D4,        ~D5,        ~D6,     ~D7,        ~Overall,
  "Qi2025",  "Serious",  "Moderate", "Low",      "Serious",  "Moderate", "Low",   "Serious",  "Serious",
  "Xu2025",  "Serious",  "Moderate", "Low",      "Serious",  "Moderate", "Low",   "Moderate", "Serious",
  "Lu2017",  "Serious",  "Serious",  "Low",      "Moderate", "Moderate", "Low",   "Moderate", "Serious",
  "Jin2019", "Serious",  "Moderate", "Low",      "Moderate", "Low",      "Low",   "Moderate", "Serious",
  "Qi2023",  "Serious",  "Moderate", "Low",      "Moderate", "Low",      "Low",   "Moderate", "Serious",
  "Li2024a", "Serious",  "Moderate", "Moderate", "Critical", "Low",      "Low",   "Moderate", "Critical")
domain_labels <- c(
  D1 = "D1 Confounding", D2 = "D2 Selection of participants",
  D3 = "D3 Classification of interventions", D4 = "D4 Deviations (co-interventions/timing)",
  D5 = "D5 Missing data", D6 = "D6 Measurement of outcome",
  D7 = "D7 Selection of reported result", Overall = "Overall")
levs <- c("Low", "Moderate", "Serious", "Critical")
rob_long <- rob |>
  pivot_longer(-study, names_to = "domain", values_to = "judgment") |>
  mutate(judgment = factor(judgment, levels = levs),
         domain = factor(domain, levels = names(domain_labels), labels = domain_labels),
         study = factor(study, levels = rev(c("Qi2025","Xu2025","Lu2017","Jin2019","Qi2023","Li2024a"))))
pal <- c(Low = "#4CAF50", Moderate = "#FBC02D", Serious = "#E53935", Critical = "#B71C1C")
txt <- c(Low = "white", Moderate = "#333333", Serious = "white", Critical = "white")
p_rob <- ggplot(rob_long, aes(domain, study, fill = judgment)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = substr(judgment, 1, 1), color = judgment),
            size = 3.4, fontface = "bold") +
  scale_fill_manual(values = pal, name = "Concern") +
  scale_color_manual(values = txt, guide = "none") +
  labs(x = NULL, y = NULL, title = "Preliminary ROBINS-I traffic-light plot",
       subtitle = "Green = low, yellow = moderate, red = serious, dark red = critical (L/M/S/C)") +
  theme_jama() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_rob, "rob_trafficlight.png", 9, 4.5)

## ---- Model fit / MCMC (model-fit.qmd) ----
s <- fit_mcmc(ma)
draws <- purrr::imap_dfr(s, function(chain, i)
  as_tibble(as.matrix(chain)[, c("mu", "tau")]) |>
    mutate(chain = factor(i), iter = row_number()))
draws_long <- pivot_longer(draws, c(mu, tau), names_to = "parameter", values_to = "value")
p_trace <- ggplot(draws_long, aes(iter, value, color = chain)) +
  geom_line(alpha = 0.6, linewidth = 0.25) +
  facet_wrap(~ parameter, ncol = 1, scales = "free_y") +
  scale_color_manual(values = jama_chains) +
  labs(x = "Iteration", y = NULL, title = "Trace plots (4 chains)") +
  theme_jama()
save_fig(p_trace, "mcmc_trace.png", 8, 5)
p_density <- ggplot(draws_long, aes(value, color = chain)) +
  geom_density() +
  facet_wrap(~ parameter, ncol = 2, scales = "free") +
  scale_color_manual(values = jama_chains) +
  labs(x = NULL, y = "Density", title = "Posterior density by chain") +
  theme_jama()
save_fig(p_density, "mcmc_density.png", 8, 4)

## ---- Goodness of fit (goodness-of-fit.qmd) ----
# Baujat influence plot
tmpf <- tempfile(fileext = ".pdf"); grDevices::pdf(tmpf)
bj <- metafor::baujat(re); grDevices::dev.off(); unlink(tmpf)
bj$study <- ma$study
p_baujat <- ggplot(bj, aes(x, y)) +
  geom_point(size = 2.8, color = jama_blue_grey[["steel"]]) +
  geom_text(aes(label = study), size = 3, hjust = -0.15) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = "Contribution to overall heterogeneity (Q)",
       y = "Influence on pooled estimate",
       title = "Baujat plot",
       subtitle = "Upper-right = studies driving both heterogeneity and the pooled OR") +
  theme_jama()
save_fig(p_baujat, "gof_baujat.png", 7.5, 5)

# Posterior predictive check
set.seed(5127)
D <- bm_sens$rposterior(n = 10000)
ppc <- purrr::map_dfr(seq_len(nrow(ma)), function(i) {
  theta <- rnorm(nrow(D), D[, "mu"], D[, "tau"])
  yrep  <- rnorm(nrow(D), theta, ma$sei[i])
  tibble(study = ma$study[i], obs_OR = exp(ma$yi[i]),
         pp_lo = exp(unname(quantile(yrep, .025))),
         pp_med = exp(unname(quantile(yrep, .5))),
         pp_hi = exp(unname(quantile(yrep, .975))))
})
p_ppc <- ppc |>
  mutate(study = factor(study, levels = rev(ma$study))) |>
  ggplot(aes(y = study)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = pp_lo, xmax = pp_hi),
                 color = jama_blue_grey[["light"]], linewidth = 3, alpha = 0.8) +
  geom_point(aes(x = pp_med), shape = 124, size = 4, color = jama_blue_grey[["steel"]]) +
  geom_point(aes(x = obs_OR), shape = 21, fill = "white",
             color = jama_blue_grey[["accent"]], size = 2.8) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8)) +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "Posterior predictive check",
       subtitle = "Bars = 95% posterior predictive interval; white points = observed ORs") +
  theme_jama()
save_fig(p_ppc, "gof_ppc.png", 8, 5)

message("Exported: freq_forest, freq_funnel, bayes_posterior_forest, bayes_posterior_mu, ",
        "bayes_posterior_tau, bayes_response_forest, sensitivity_loo, rob_trafficlight, ",
        "mcmc_trace, mcmc_density, gof_baujat, gof_ppc (300 dpi) to figures/")
