# PCP combination therapy — Bayesian random-effects meta-analysis
# PROSPERO question: In HIV-negative/immunocompromised adults with PCP, comparative
# effectiveness of combination/adjunctive regimens vs monotherapy.
#
# PRIMARY CONTRAST (this script): TMP-SMX + adjunctive echinocandin/caspofungin
#   vs TMP-SMX monotherapy. Primary outcome: all-cause mortality (in-hospital / ~30-day).
#
# Data source: PCP_combination_extraction_database_v0.2_machine-draft.xlsx
#   NB: v0.2 is a MACHINE FIRST-PASS extraction — every value needs dual human verification.
#
# Scope decisions (locked with maintainer 2026-08-06):
#   - Primary = adjunctive echinocandin vs TMP-SMX monotherapy.
#   - Exclude confounded studies (corticosteroid dose/timing comparisons) and the
#     Taniguchi studies with back-calculated events from the primary pool.
#   - Effect measure = odds ratio; timepoint = 30-day preferred, else in-hospital.
#   - Duplicate Qi2025 == Yanmeng2026 (same PUMCH cohort) -> keep Qi2025 only.

library(tidyverse)
library(metafor)
library(bayesmeta)

# ---- Source verification (counts checked against full-text extractions) -----
# All arm counts below were confirmed against the source PDFs' plain text.
# Two studies were EXCLUDED at verification (not eligibility/analysis errors in
# the counts themselves, but in the machine draft's classification):
#   - Tian2020: published in HIV Medicine; cohort is HIV-INFECTED patients
#       ("we only included HIV-infected patients") -> violates PROSPERO
#       HIV-negative/non-HIV inclusion criterion. EXCLUDED.
#   - Li2024a: "synergistic" arm = TMP-SMX + caspofungin + glucocorticoid
#       (+ lower TMP dose) vs TMP-SMX alone -> confounded; cannot isolate a
#       caspofungin effect. RE-INCLUDED per user decision (2026-08-06), kept in
#       the sensitivity tier and flagged confounded. NB: its reference arm is the
#       only steroid-FREE comparator in the set.
# Additional caveats (do not affect inclusion but affect interpretation):
#   - Qi2025 == Yanmeng2026 (byte-identical file; same PUMCH cohort) -> keep one.
#     In the OVERALL cohort the caspofungin arm had HIGHER mortality (51.4% vs
#     39.2%); the "benefit" appears only in the ventilated subgroup. Unadjusted,
#     exploratory subgroup from a risk-factor paper.
#   - Xu2025: observational; monotherapy arm was sicker (more solid tumour,
#     more aspergillosis co-infection) -> confounding by indication. Raw counts
#     are unadjusted; authors' Cox (HR 2.37) and IPTW (HR 2.97) agree in direction.
#   - Lu2017: n=5 vs 6 case series; mortality timepoint not explicitly defined.

# ---- Analysis dataset (verified, full-cohort rows) --------------------------
# Arm counts are events/N for the combination (echinocandin) and reference (mono) arms.
# `timing` (initial / salvage / mixed) classified per Yang et al. (2026) Table 1.
# No pure salvage-vs-monotherapy study exists in our verified set.
ma <- tribble(
  ~study,     ~timepoint,          ~e_comb, ~n_comb, ~e_ref, ~n_ref, ~tier,          ~timing,
  "Qi2025",   "in-hospital",            18,      35,     31,     79, "primary",      "mixed",
  "Xu2025",   "28-day",                 33,      89,     19,     31, "primary",      "mixed",
  "Lu2017",   "episode/in-hosp",         1,       5,      2,      6, "primary",      "mixed",
  # Pre-planned sensitivity: longer mortality timepoints
  "Jin2019",  "3-month",                 8,      35,     33,     91, "sensitivity",  "initial",
  "Qi2023",   "90-day",                 17,      43,     10,     21, "sensitivity",  "initial",
  "Li2024a",  "episode (unspec)",        5,      20,      3,     18, "sensitivity",  "initial"  # confounded (user-included)
)

# ---- Per-study log odds ratios ----------------------------------------------
ma <- escalc(measure = "OR",
             ai = e_comb, n1i = n_comb,   # combination (echinocandin) arm
             ci = e_ref,  n2i = n_ref,    # reference (monotherapy) arm
             data = ma, slab = study) |>
  mutate(sei = sqrt(vi), OR = exp(yi))

# ---- Bayesian random-effects meta-analysis (normal-normal) ------------------
# Weakly informative priors: mu ~ N(0, 1.5) on log-OR; tau ~ half-normal(0.5).
fit_bm <- function(d) {
  bayesmeta(y = d$yi, sigma = d$sei, labels = d$study,
            mu.prior = c(mean = 0, sd = 1.5),
            tau.prior = function(t) dhalfnormal(t, scale = 0.5))
}

bm_primary <- fit_bm(filter(ma, tier == "primary"))
bm_all     <- fit_bm(ma)

report <- function(bm, label) {
  mu <- bm$summary[, "mu"]; tau <- bm$summary[, "tau"]
  cat("\n==== ", label, " (k = ", bm$k, ") ====\n", sep = "")
  cat(sprintf("Pooled OR (median): %.2f  95%% CrI [%.2f, %.2f]\n",
              exp(mu["median"]), exp(mu["95% lower"]), exp(mu["95% upper"])))
  cat(sprintf("tau (median): %.2f  95%% CrI [%.2f, %.2f]\n",
              tau["median"], tau["95% lower"], tau["95% upper"]))
  cat(sprintf("P(OR < 1 | data): %.3f\n", bm$pposterior(mu = 0)))
}
report(bm_primary, "PRIMARY (in-hospital / ~30-day mortality)")
report(bm_all,     "SENSITIVITY (+ longer-timepoint studies)")

# ---- Timing subgroups (initial vs mixed), per Yang et al. (2026) --------------
bm_initial <- fit_bm(filter(ma, timing == "initial"))  # Jin2019, Qi2023, Li2024a
bm_mixed   <- fit_bm(filter(ma, timing == "mixed"))     # Qi2025, Xu2025, Lu2017
report(bm_initial, "INITIAL-strategy subgroup")
report(bm_mixed,   "MIXED-strategy subgroup")
