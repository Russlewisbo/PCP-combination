# Shared setup for the PCP echinocandin meta-analysis website.
# Sourced at the top of each page; loads packages, data, and fitting helpers.

suppressPackageStartupMessages({
  library(tidyverse)
  library(metafor)
  library(bayesmeta)
})
source("analysis/jama_style.R")

# Load the verified analysis dataset and compute per-study log odds ratios.
load_ma <- function(path = "analysis/pcp_primary_meta_dataset.csv") {
  d <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::select(study, tier, timing, timepoint, e_comb, n_comb, e_ref, n_ref)
  d <- metafor::escalc(measure = "OR", ai = e_comb, n1i = n_comb,
                       ci = e_ref, n2i = n_ref, data = d, slab = study)
  d$sei <- sqrt(d$vi)
  d$OR  <- exp(d$yi)
  d
}

# Bayesian random-effects (normal-normal) with weakly informative priors:
# mu ~ N(0, 1.5) on the log-OR; tau ~ half-normal(0.5).
fit_bm <- function(d) {
  bayesmeta::bayesmeta(y = d$yi, sigma = d$sei, labels = d$study,
                       mu.prior  = c(mean = 0, sd = 1.5),
                       tau.prior = function(t) bayesmeta::dhalfnormal(t, scale = 0.5))
}

or_ci <- function(bm) c(exp(bm$summary["median", "mu"]),
                        exp(bm$summary["95% lower", "mu"]),
                        exp(bm$summary["95% upper", "mu"]))

fmt <- function(bm) { v <- or_ci(bm); sprintf("%.2f (95%% CrI %.2f\u2013%.2f)", v[1], v[2], v[3]) }

# --- JAGS: the same normal-normal model, fitted by MCMC for convergence checks ---
jags_model_text <- "model{
  for (i in 1:k) {
    y[i]     ~ dnorm(theta[i], prec[i])   # prec[i] = 1 / sei^2 (within-study)
    theta[i] ~ dnorm(mu, inv_tau2)         # study effect
  }
  mu  ~ dnorm(0, 0.4444444)                # N(0, 1.5) on log-OR  (precision 1/1.5^2)
  tau ~ dnorm(0, 4) T(0, )                 # half-normal(0.5) on tau
  inv_tau2 <- 1 / (tau * tau)
  OR <- exp(mu)
}"

fit_jags <- function(d, monitor = c("mu", "tau", "OR", "theta"),
                     n.chains = 4, n.adapt = 1000, n.burn = 2000,
                     n.iter = 10000, seed = 4242) {
  library(rjags); library(coda)
  set.seed(seed)
  jm <- jags.model(textConnection(jags_model_text),
                   data = list(y = d$yi, prec = 1 / d$sei^2, k = nrow(d)),
                   n.chains = n.chains, n.adapt = n.adapt, quiet = TRUE)
  update(jm, n.burn, progress.bar = "none")
  coda.samples(jm, monitor, n.iter = n.iter, progress.bar = "none")
}

# Four-colour blue-grey ramp for MCMC chains.
jama_chains <- c("#A9BCCF", "#6E8CA8", "#3C5F80", "#243B53")
