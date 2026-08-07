# PCP Combination Therapy — Bayesian meta-analysis

**Live site:** https://Russlewisbo.github.io/PCP-combination/

Interim systematic review and meta-analysis of **adjunctive echinocandin/caspofungin +
trimethoprim-sulfamethoxazole (TMP-SMX) versus TMP-SMX monotherapy** for all-cause
mortality in HIV-negative / immunocompromised adults with *Pneumocystis jirovecii*
pneumonia (PCP/PJP). PROSPERO-registered; primary outcome is all-cause mortality
(30-day and/or in-hospital), reported as an odds ratio with 95% credible intervals.

## Headline result

Pooling six source-verified studies, adjunctive echinocandin shows a **possible but
inconclusive** mortality benefit (Bayesian sensitivity OR 0.76, 95% CrI 0.43–1.37;
frequentist REML OR 0.75, 95% CI 0.42–1.36) — the interval includes no effect under both
frameworks.

## Website

The published site (built with Quarto) presents the full analysis across pages for the
included studies, PRISMA-style selection flow, frequentist and Bayesian results, MCMC
convergence diagnostics (JAGS), sensitivity analyses, a preliminary ROBINS-I risk-of-bias
assessment, and conclusions. It includes an **EN / IT** Google Translate toggle and
collapsible code.

- Live: https://Russlewisbo.github.io/PCP-combination/
- Build locally: `quarto preview` (open in an external browser for the translate toggle)
- Republish: `quarto publish gh-pages --no-prompt`

## Repository layout

| Path | Contents |
|---|---|
| `index.qmd`, `studies.qmd`, `prisma.qmd`, `frequentist.qmd`, `bayesian.qmd`, `model-fit.qmd`, `sensitivity.qmd`, `risk-of-bias.qmd`, `conclusions.qmd` | Website pages |
| `pcp_echinocandin_meta_writeup.qmd` | Single-page full write-up |
| `analysis/meta_helpers.R` | Shared data-loading and model-fitting helpers |
| `analysis/jama_style.R` | JAMA blue-grey figure theme + palette |
| `analysis/pcp_meta_analysis.R` | Reproducible analysis script |
| `analysis/pcp_primary_meta_dataset.csv` | Verified analysis dataset |
| `_quarto.yml` | Quarto website configuration |
| `AGENTS.md` | Project memory + session log |

## Status

This is an **interim, machine-assisted** analysis pending human verification: dual
extraction verification and a formal ROBINS-I assessment are still to be completed, and
several Chinese-language studies remain to be retrieved. See `AGENTS.md` for details.
