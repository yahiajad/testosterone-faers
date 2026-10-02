# Testosterone × sleep apnoea / polycythaemia — FAERS disproportionality analysis

Code for the study *"Testosterone and reports of sleep apnoea and polycythaemia in men: an age-stratified pharmacovigilance study of the FDA Adverse Event Reporting System, 2012–2026."*

The pipeline downloads the public FDA Adverse Event Reporting System (FAERS) quarterly ASCII files, deduplicates them, and computes four disproportionality measures (ROR, PRR, IC, EBGM) for testosterone with sleep apnoea and polycythaemia in men. It includes age strata, sensitivity analyses, time to onset and figures. No patient-level data are stored in this repository.

## Requirements
- Python ≥ 3.11 and `pip install -r requirements.txt`
- `bash`, `curl`, `unzip`
- About 3 GB of disk for the raw files, plus about 6 GB for the working database

## Run
```bash
export FAERS_DIR=$PWD/data          # optional; this is the default
bash scripts/01_download.sh          # FAERS 2012Q4 → latest (8 parallel downloads, re-runnable)
python scripts/02_build.py           # zip → parquet per table
python scripts/03_analysis.py --rebuild   # dedupe + ROR/PRR/IC, strata, sensitivity, time to onset
python scripts/04_ebgm.py            # MGPS / EBGM
python scripts/05_descriptives.py    # case characteristics, age-share comparisons
python scripts/06_figures.py         # Figures 1–4
```
Outputs go to `$FAERS_DIR/results/` and `$FAERS_DIR/figures/`. The study used the quarterly extracts from 2012Q4 to 2026Q2. Later quarters will change the counts slightly.

## Definitions (as in the paper)
- **Exposure:** testosterone (esters and brand names; methyltestosterone and oestrogen/progesterone combinations excluded) as primary suspect drug.
- **Sleep apnoea:** MedDRA PTs sleep apnoea syndrome, obstructive sleep apnoea syndrome, central sleep apnoea syndrome.
- **Polycythaemia:** MedDRA PTs polycythaemia, secondary polycythaemia, haematocrit increased, haemoglobin increased, red blood cell count increased.
- **Population:** male reports. Age strata 18–39, 40–64 and 65–119 years.
- **Signal:** ≥ 3 cases, lower ROR 95% CI > 1 and IC025 > 0. EB05 ≥ 2 is reported for MGPS.

## Licence
MIT. See `LICENSE`. FAERS data are public-domain US government data.

## AI assistance
Analysis code was written with the assistance of Claude (Anthropic) and reviewed by the author.
