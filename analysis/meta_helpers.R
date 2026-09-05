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

# --- MCMC: the same normal-normal model, refit by a self-contained random-walk
# Metropolis sampler in base R (no JAGS/Stan dependency) for convergence checks.
# Samples (mu, log-tau) from the marginal posterior with study effects integrated
# out, using a Haario-style adaptive proposal covariance during burn-in.
fit_mcmc <- function(d, n.chains = 4, n.burn = 3000, n.iter = 10000, seed = 7314) {
  y <- d$yi; s2 <- d$sei^2
  lpost <- function(mu, eta) {
    tau <- exp(eta); v <- s2 + tau^2
    sum(dnorm(y, mu, sqrt(v), log = TRUE)) +
      dnorm(mu, 0, 1.5, log = TRUE) +                    # N(0, 1.5) on log-OR
      (log(2) + dnorm(tau, 0, 0.5, log = TRUE)) + eta    # half-normal(0.5) + Jacobian
  }
  set.seed(seed)
  inits <- list(c(0.6, log(0.6)), c(-0.6, log(0.1)),
                c(0.3, log(0.3)), c(-0.3, log(0.2)))
  run_chain <- function(init) {
    x <- init; lp <- lpost(x[1], x[2])
    S <- diag(c(0.1, 0.1)); mean_x <- x; cov_x <- diag(0, 2); n_ad <- 0
    keep <- matrix(NA_real_, n.iter, 2)
    for (t in seq_len(n.burn + n.iter)) {
      prop <- x + as.numeric(chol(S) %*% rnorm(2))
      lp_p <- lpost(prop[1], prop[2])
      if (log(runif(1)) < lp_p - lp) { x <- prop; lp <- lp_p }
      if (t <= n.burn) {                                 # adapt proposal during burn-in
        n_ad <- n_ad + 1; d_x <- x - mean_x; mean_x <- mean_x + d_x / n_ad
        cov_x <- cov_x + (d_x %o% (x - mean_x))
        if (n_ad > 50) S <- (2.38^2 / 2) * (cov_x / (n_ad - 1)) + diag(1e-6, 2)
      } else keep[t - n.burn, ] <- x
    }
    keep
  }
  coda::mcmc.list(lapply(inits[seq_len(n.chains)], function(init) {
    draws <- run_chain(init)
    coda::mcmc(cbind(mu = draws[, 1], tau = exp(draws[, 2]), OR = exp(draws[, 1])))
  }))
}

# Four-colour blue-grey ramp for MCMC chains.
jama_chains <- c("#A9BCCF", "#6E8CA8", "#3C5F80", "#243B53")
