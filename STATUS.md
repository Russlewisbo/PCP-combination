# PCP Combination — Project STATUS

**Project:** Bayesian meta-analysis of combination / adjunctive therapy for *Pneumocystis jirovecii* pneumonia (PCP/PJP)
**Maintainer:** R.E.L. (Russ)
**Project home:** `~/main/PCP_combo_thesis` (connected Cowork folder — this is the persistent home; save all outputs here).
**This file:** running log of decisions, what was done, and where to pick up. **Update at the end of every working session.**

---

## 0. How to use this file
- Read this first at the start of each session to re-orient.
- At the end of each session, update **Section 5 (Session log)** and **Section 4 (Next steps)**.
- Anything marked **[PROVISIONAL]** is an assumption Claude made because the PROSPERO record / protocol was not available in-session — confirm and lock these against the real protocol.

---

## 1. Current status (as of 2026-08-06)
- **Stage:** Setting up the data-collection instrument for extraction.
- **Done today:** Built v0.1 of the structured extraction database (`PCP_combination_extraction_database_v0.1.xlsx`), modelled on the MSSA cefazolin-vs-ASP NMA database, and created this STATUS.md.
- **Blocker:** The project's `STATUS.md` and `PROSPERO_field_answers` / protocol did **not** sync into the Cowork session (only project metadata came through). The extraction form was therefore built from standard PCP-combination-therapy assumptions and must be reconciled with the actual PROSPERO record before extraction begins.

## 2. Scope — PICO [PROVISIONAL, confirm vs PROSPERO/protocol]
- **Population:** Adults with proven/probable PCP. HIV status recorded as a primary effect modifier (HIV-positive vs non-HIV immunocompromised). *Confirm inclusion of paediatric / mixed populations.*
- **Intervention (nodes):** Combination / adjunctive anti-PCP regimens — the default network scaffold covers:
  - TMP-SMX (reference)
  - TMP-SMX **+ adjunctive echinocandin** (caspofungin / micafungin / anidulafungin)
  - TMP-SMX + caspofungin only (if echinocandins not pooled)
  - Clindamycin–primaquine
  - Atovaquone
  - IV pentamidine
  - Dapsone ± TMP
  - "Combination (other, specify)" and "Monotherapy (pooled)" for sensitivity analyses
  - *Confirm the exact node list and whether this is a full NMA or pairwise combination-vs-monotherapy.*
- **Comparator:** TMP-SMX monotherapy (and/or standard of care).
- **Primary outcome [PROVISIONAL]:** All-cause mortality (in-hospital / 30-day / 90-day). *Confirm.*
- **Secondary outcomes:** Treatment failure (study-defined), respiratory failure / need for mechanical ventilation, ICU admission, change of therapy, length of stay, adverse-event discontinuation, β-D-glucan response/kinetics.
- **Designs / RoB [PROVISIONAL]:** RCTs (RoB 2) + non-randomised studies (ROBINS-I). *Confirm whether observational studies are eligible.*
- **Effect modifiers / meta-regression covariates:** HIV status, adjunctive corticosteroid use, baseline severity (PaO₂ / A–a gradient / respiratory support at baseline), baseline β-D-glucan, time from admission to treatment.

## 3. Requested additions (user, this session) — INCLUDED in v0.1
- **β-D-glucan (1,3-β-D-glucan / BDG):** captured on a dedicated sheet — baseline level, peak, follow-up/kinetics, assay, cut-off, units, and its relationship to outcome. Per-arm baseline BDG also on the Arms sheet as a covariate.
- **Concomitant corticosteroids:** captured per arm on the Arms sheet (% receiving adjunctive corticosteroids, regimen, indication) and available as an effect-modifier column.

## 4. Next steps
1. **Reconcile scope with the real PROSPERO record + protocol** — lock the node list, primary outcome, eligible designs, and effect modifiers; clear all **[PROVISIONAL]** flags.
2. Rename the workbook to v1.0 once scope is confirmed and freeze the Codebook.
3. Populate `7_Study_index` with the included study list (from screening) and assign reviewers.
4. Run the dual-extraction calibration (first ~5–10 studies), then reconcile.
5. Decide β-D-glucan handling in the model (covariate vs outcome; continuous vs threshold).
6. Confirm corticosteroid variable definition (any use vs guideline-indicated use; dose threshold).

## 5. Session log
| Date | Who | What happened / decisions | Where we left off |
|---|---|---|---|
| 2026-08-06 | Claude + R.E.L. | Project docs (STATUS/PROSPERO) did not sync. Built extraction DB v0.1 from PCP assumptions incl. β-D-glucan + corticosteroid fields. Created this STATUS.md. Transferred deliverables + scratch into project home `~/main/PCP_combo_thesis`. | Awaiting PROSPERO/protocol to confirm PICO and lock scope. |
| 2026-08-06 | Claude + R.E.L. | Review_PDFs initially held the wrong (S. aureus) papers; R.E.L. swapped in the 48 PCP full texts. Ran machine first-pass extraction of all 48 via parallel agents → `PCP_combination_extraction_database_v0.2_machine-draft.xlsx` (48 studies, 85 arms, 278 outcome rows, 63 BDG, 66 steroid rows). Recalc clean, vocab conforms. | v0.2 is a MACHINE DRAFT — needs dual human verification. Next: dedup Qi2025/Yanmeng2026, verify the 3 Taniguchi % → event conversions, confirm which contrasts enter the network, then build R analysis dataset. |

## 5a. Extraction summary (v0.2 machine draft) — KEY CAVEATS
- **This is a machine first-pass (Reviewer = "Claude (machine draft)"), NOT verified.** Two human reviewers must check every value; quotes + pages are provided on 3_Outcomes/4/5 for that purpose.
- **Eligibility split:** Include-network **9**, Include-pairwise **8**, Include-sensitivity **5**, Background/context-only **26**. Most of the corpus is single-arm / prognostic / steroid-timing / dose studies without an eligible antifungal-regimen contrast.
- **Genuine antifungal combination-vs-monotherapy contrasts with per-arm counts:** Gu2022, Jin2019, Liu2022, Qi2023, Qi2025(=Yanmeng2026 dup), Tian2020, Xu2025, Lu2017, Wang2019, Zhang2018, Taniguchi2024b, Li2024a.
- **Duplicate:** `Qi2025` and `Yanmeng2026` appear to be the SAME PUMCH cohort/paper — dedupe before analysis.
- **Computed (not reported) events:** the 3 Taniguchi studies (2024a low vs conventional dose; 2025 pentamidine vs TMP-SMX; 2024b + echinocandins) report only %/risk-differences; events were back-calculated and flagged `"uncertain - verify"`.
- **Confounding flags:** several "combination" studies are really corticosteroid-dose or treatment-timing comparisons (Lemiale2025, LiX2024, WangL2026, Reizine2026, Mundo2020, Wang2022, Zhang2018) — node labels use "Other/(specify)" and notes flag this for reviewers.

## 6. Key files (in project home `~/main/PCP_combo_thesis`)
- `PCP_combination_extraction_database_v0.2_machine-draft.xlsx` — **populated** extraction DB (48 studies, machine first-pass).
- `PCP_combination_extraction_database_v0.1.xlsx` — empty template (blank instrument).
- `STATUS.md` — this file.
- `Review_PDFs/` — 48 PCP full-text PDFs (the extraction corpus).
- `review_PCP.pdf` — Li et al. 2024 (a key included study).
- `Paulina_177_priority_for_Covidence.ris` — 177 references for Covidence screening.
- `scratch/pdf_text/` — plain-text of each PDF; `scratch/extraction/` — per-batch JSON + extraction spec; `scratch/load_extraction.py`, `scratch/build_pcp_db.py` — scripts.
- *(to add)* PROSPERO record / protocol; screening decisions; R analysis dataset.

## 6. Key files (in project home `~/main/PCP_combo_thesis`)
- `PCP_combination_extraction_database_v0.1.xlsx` — extraction instrument (this session).
- `STATUS.md` — this file.
- `review_PCP.pdf` — existing review document (pre-existing, 24 Jul 2026).
- `Paulina_177_priority_for_Covidence.ris` — 177 priority references for Covidence screening (pre-existing).
- `Review_PDFs/` — folder of ~56 full-text PDFs for extraction (pre-existing).
- `scratch/build_pcp_db.py` — script that generated the extraction workbook.
- *(to add)* PROSPERO record / protocol; screening/eligibility decisions.
