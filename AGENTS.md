# PCP Combination — Project Memory & Session Log

> **How to use this file:** Read the **Current status** and **Next steps** sections first
> to re-orient. At the **end of every working session**, append a row to the
> **Session log** and refresh **Current status** / **Next steps**. Anything provisional
> should be flagged so it can be reconciled later.

Bayesian meta-analysis of combination/adjunctive therapy for *Pneumocystis jirovecii*
pneumonia (PCP/PJP) in **HIV-negative / immunocompromised adults** (PROSPERO-registered).
Primary question: comparative effectiveness/safety of combination/adjunctive regimens
vs monotherapy (TMP-SMX). Primary outcome: all-cause mortality (30-day and/or in-hospital),
OR with 95% credible intervals. Plan: Bayesian random-effects meta-analysis (NMA if data permit).

**Live site:** https://Russlewisbo.github.io/PCP-combination/ (published from the `gh-pages`
branch via `quarto publish gh-pages --no-prompt`).

## Current status (as of 2026-08-06)
- Verified primary analysis complete: **adjunctive echinocandin/caspofungin + TMP-SMX vs
  TMP-SMX monotherapy**, all-cause mortality.
- Headline: primary (k=3) **OR 0.78 (95% CrI 0.31–1.93)**; sensitivity (k=6)
  **OR 0.76 (0.43–1.37)**. Suggestive of benefit but inconclusive (CrI includes 1).
- Timing subgroups added (per Yang et al. 2026 Table 1): initial (k=3) **OR 0.74 (0.34–1.62)**;
  mixed (k=3) **OR 0.78 (0.31–1.93)**. No salvage-vs-monotherapy study available on our side.
- Benchmarked against the Yang et al. 2026 echinocandin meta-analysis (see below).

## Figure style (ALWAYS APPLY)
All figures must be **JAMA style with a blue-grey palette**. `source("analysis/jama_style.R")`,
then add `theme_jama()` and use `jama_blue_grey` / `jama_tier` for colours.
No panel gridlines; classic axes; sans font.

## Environment / gotchas
- **Stan is broken on this machine** (rstan TBB linker error; no cmdstanr). Fit Bayesian
  models with **`bayesmeta`** (semi-analytic, no MCMC). `rjags` is available as a fallback.
- Priors used throughout: `mu ~ N(0, 1.5)` on log-OR; `tau ~ half-normal(0.5)`.
- Helper `fit_bm()` wraps `bayesmeta()` with these priors (defined in the analysis script).

## Key files
- `analysis/pcp_meta_analysis.R` — reproducible meta-analysis script (data, log-ORs, models, subgroups).
- `analysis/pcp_primary_meta_dataset.csv` — verified analysis dataset (has `timing` column).
- `analysis/jama_style.R` — figure theme + palette.
- `PCP_combination_extraction_database_v0.2_machine-draft.xlsx` — extraction DB
  (**MACHINE FIRST-PASS, unverified**; needs dual human verification).
- Source PDFs: `~/main/PCP_combo_thesis/Review_PDFs/`; plain text in
  `~/main/PCP_combo_thesis/scratch/pdf_text/` (used to verify counts against source).
- `STATUS.md` — older running log (superseded by this file for analysis work).

## Reading the Excel DB
Sheets have a title row above the header → `read_excel(path, sheet, skip = 1)`.
First data row (`Smith2019`) is an EXAMPLE row — drop it. Key sheets:
`1_Study_characteristics`, `2_Arms_covariates`, `3_Outcomes`, `4_Beta_D_glucan`,
`5_Corticosteroids`, `8_Study_index` (eligibility).

## Locked analysis decisions
- Primary contrast: adjunctive echinocandin/caspofungin + TMP-SMX vs TMP-SMX monotherapy.
- Effect measure: odds ratio. Timepoint: 30-day preferred, else in-hospital.
- Duplicate `Qi2025 == Yanmeng2026` (byte-identical, same PUMCH cohort) → keep Qi2025 only.

## Study inclusion (after SOURCE verification)
Verified included set — primary = **Qi2025, Xu2025, Lu2017**; sensitivity adds
**Jin2019, Qi2023, Li2024a**. All extracted counts confirmed accurate against source text.
- **Excluded Tian2020** — HIV-POSITIVE cohort (published in *HIV Medicine*); machine draft
  misclassified it. (Yang et al. independently excluded it too.)
- **Li2024a** — confounded (caspofungin bundled with corticosteroid + lower TMP dose); kept
  in the SENSITIVITY tier per user decision, flagged confounded. Its reference arm is the
  only steroid-free comparator in the set.
- **Qi2025 extraction confirmed**: full-cohort Table 4 = 18/35 (caspo) vs 31/79 (mono).
  Yang et al.'s 18/27 & 62/82 come from the *ventilated subgroup* (Table 5) + an incorrect
  pooling of monotherapy+clindamycin deaths — i.e. Yang analysed a subgroup as the whole study.

## Sub-analyses NOT estimable
- Corticosteroid presence/absence: near-universal steroid use across studies; only Li2024a
  has a steroid-free reference arm (single confounded data point).
- β-D-glucan yes/no: only one small study (Lu2017) omitted BDG → no informative contrast.

## Comparison with Yang et al. 2026 (SSRN preprint, echinocandin + TMP-SMX vs TMP-SMX)
- We share **9 of their 16** studies; missing 7 (mostly Chinese-language papers / master's
  theses: Li2015, Liu2015, Wang2021, Wu2023, Xiang2015, Yu2017, Wang-ZG-2019 SLE) because
  their search included CNKI/Wanfang.
- Yang overall OR 0.93 (0.71–1.21) is null; their **initial-strategy OR 0.50 (0.32–0.77)**.
- **Timing explains the gap.** Restricting our set to initial-strategy studies (dropping the
  mixed Qi2025) gives OR 0.59 (0.32–1.09), consistent with Yang's initial signal.
- Caveats on Yang: frequentist FIXED-effect M-H (narrow CIs); visible extraction errors
  (Fig 4 ORs 21.32 & 19.00; the Qi2025 subgroup mix-up above). Useful comparator, not ground truth.
- Nuance: Yang plotted Xu2025 as *initial* (Fig 3A) but Table 1 calls it initial+salvage; we
  classified it **mixed** (SOURCE-CONFIRMED: Xu2025 never restricts caspofungin to first-line
  use, reports no time-to-caspofungin, and flags "selection bias in caspofungin use"). This
  makes our strict initial subgroup (0.74) less protective than Yang's (0.50).

## Retrieval status of the 7 missing studies (checked 2026-08-06)
None are in our local corpus; none could be added to the analysis this session.
- **Wang-ZG-2019 (SLE, Medicine Baltimore; PMID 31169741, PMC6571266, open access)** —
  OBTAINED (Europe PMC). **INELIGIBLE**: single-arm case series (9 patients, all combination,
  0 deaths); no monotherapy comparator. Yang's "0/9 vs 2/6" invents a comparator absent from
  the source — another Yang extraction error. Do NOT add.
- **Yu-2017 (Int J Clin Exp Med 2017;10(1):1234–1242)** — low retrievability; that journal
  ceased publication (~2018), no reliable DOI/PMID. Try ijcem.com archive / ResearchGate.
- **Chinese-language, CNKI/Wanfang only (need institutional access; 2 are theses):**
  Li-T-2015 (thesis, Capital Medical Univ), Liu-2015 (Chin J Mod Drug Appl 9(17):90–91),
  Wang-Z-2021 (thesis, Hebei Medical Univ; DOI 10.27111/d.cnki.ghyku.2021.000212),
  Wu-2023 (J Internal Intensive Med 29(5):393–396), Xiang-2015 (Chin J Clin Rational Drug Use
  8(22):69–70). Route: team's Chinese-reading members via cnki.net / wanfangdata.com.cn.

## Next steps / open questions
1. Obtain the 5 CNKI/Wanfang Chinese studies via institutional access (team's Chinese readers);
   Wang-ZG-2019 is ineligible (single-arm) and Yu-2017 is low-accessibility.
2. Complete dual human verification of the extraction DB; have a reviewer sign off the AI
   disclaimer at the top of `pcp_echinocandin_meta_writeup.qmd`.
3. Consider a binomial-likelihood hierarchical model (needs a working Stan toolchain).
4. (Optional) Render the write-up to Word/PDF for the thesis team.

### Done (2026-08-06)
- Xu2025 timing resolved → **mixed** (source-confirmed).
- Leave-one-out + prior-sensitivity checks complete (see status below).
- Reproducible JAMA-styled write-up drafted: `pcp_echinocandin_meta_writeup.qmd` (renders clean).

## Session log
| Date | Who | What happened / decisions | Where we left off |
|---|---|---|---|
| 2026-08-06 | Russ + Assistant | Set up Bayesian meta-analysis from the v0.2 machine draft; verified all counts against source PDFs; excluded Tian2020 (HIV+), de-duplicated Qi2025/Yanmeng2026; re-included confounded Li2024a per user; established JAMA blue-grey figure style; benchmarked vs Yang et al. (timing explains the difference); confirmed our Qi2025 extraction is correct and Yang's is a ventilated-subgroup error; added `timing` variable + initial/mixed subgroups; built side-by-side forest plot. | Primary OR 0.78 (0.31–1.93); analysis reproducible in `analysis/pcp_meta_analysis.R`. Pick up at "Next steps" above. |
| 2026-08-06 | Russ + Assistant | Attempted retrieval of the 7 missing studies (none obtainable; Wang-ZG-2019 obtained but ineligible single-arm — another Yang extraction error). Ran robustness on the 6-study sensitivity model: **leave-one-out** pooled OR 0.59 (drop Qi2025) → 0.94 (drop Xu2025), all CrIs cross 1, Xu2025 most influential; **prior sensitivity** OR stable 0.75–0.77 across 7 priors (only CrI width responds to the tau prior). Drafted reproducible JAMA-styled Quarto write-up `pcp_echinocandin_meta_writeup.qmd` (loads from CSV + jama_style.R; renders clean; missing Chinese studies flagged as a limitation). | Headline unchanged: sensitivity OR 0.76 (0.43–1.37), inconclusive & robust. Write-up awaits human review/sign-off. Pick up at "Next steps". |
| 2026-08-07 | Russ + Assistant | Built a Quarto **website** (`_quarto.yml`, `index.qmd`, `styles.css`, `assets/translate-*.html`) with a persistent cookie-based Google Translate **EN/IT toggle** (hardened to https+async; note: toggle needs an external browser — Positron's Viewer blocks it). Pushed initial site, then **expanded to 10 pages**: Overview, Included studies, PRISMA (Mermaid), Frequentist (metafor REML OR 0.75, I²=42%, + Mantel-Haenszel, forest+funnel), Bayesian (posterior forest + marginal posteriors), **Model fit & MCMC** (genuine JAGS refit: R̂≈1.00, ESS in thousands, trace/density, JAGS≈bayesmeta cross-check), Sensitivity, preliminary **ROBINS-I** risk-of-bias (traffic-light), Conclusions, + retained full write-up. Shared helpers in `analysis/meta_helpers.R`. Enabled site-wide **code-folding** (`code-fold: true`, "Show code"). Committed + pushed (75cd8e9, e741bec). | Site renders clean to `_site/` (gitignored). RoB page is PRELIMINARY (needs human ROBINS-I); MCMC page uses JAGS because headline bayesmeta fit is semi-analytic. Preview via `quarto preview` in an external browser. |
