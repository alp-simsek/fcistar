# FCI-star Live Website — Claude Code Instructions

## What This Project Is

We are building a live, auto-updating public website for **FCI\*** — the neutral level of a
financial conditions index consistent with output at potential. FCI\* is the financial-conditions
analog of r\*: unlike r\*, it is insulated from financial market fluctuations, and FCI gaps
(FCI minus FCI\*) provide real-time guidance on monetary policy stance.

**Reference design:** The NY Fed r-star page (https://www.newyorkfed.org/research/policy/rstar).
Our goal is to replicate and improve on it — with interactive, auto-updating charts rather than
static images, and eventually multiple series (FCI\*, FCI-T, FCI gaps).

**Paper:** `backend/refs/FCIstarPublic.pdf` — read this to understand the model and estimation.

**Reference code:** `backend/refs/matlab/` — the most recent stable Matlab codebase (dated `26_01_10`) that produces the figures in the paper. This code uses static data, not live feeds. It is a read-only reference — do not modify it.

---

## Current Status (as of 2026-04-26)

**Site is live and fully auto-updating.** https://fcistar.org shows three interactive Plotly charts (FCI & FCI\*, FCI gap, output gap), a global date-range selector (1Y / 5Y / 10Y / 20Y / All) that rescales all charts simultaneously, and a header showing "Estimates through YYYY QN · Last updated Month YYYY".

**Pipeline runs automatically on the 1st of each month** (06:00 UTC). It pulls fresh FRED, FCI-G, and HLW data, runs Kentaro's Kalman estimation, commits any output changes, and auto-redeploys the site via the `workflow_run` chain in `deploy.yml`.

**Current sample end:** 2025Q4. Advances dynamically — the `sample_window` function in `estimate_CCS_struct.py` now reads the latest available quarter from `data_hlm.csv` when `recent_data=3` (commit `aa3401a`).

**Vintages archive live and link shipped (2026-04-26).** Kentaro's `write_outputs` (commits `faa2aa7` / `055fd2b`) now also writes `fcistar_YYYY-MM-DD.csv` and `metadata_YYYY-MM-DD.json` into `backend/data/output/vintages/`, keyed by pipeline run date. Existing vintage files are never overwritten — the function raises `FileExistsError` if a file for the run date already exists. Folder was first populated by a manual local run today (commit `ac57496`). The header now includes a "Past vintages" link pointing to that folder on GitHub (commit `2acd548`). Email sent to Kentaro acknowledging; status update sent to Ricardo and Tomas.

### Daily FCI nowcast + one-quarter-ahead FCI\* forecast + sensitivity (live, 2026-06)

A large body of work shipped in June 2026, all auto-updating on a **weekday cron**
(`.github/workflows/update-nowcast.yml`, 22:00 UTC Mon–Fri). The daily pipeline (in
`backend/pipeline/`, run with `working-directory: backend/pipeline`) is:

```
get_fcig.py → daily_backcast.py → nowcast_fcig.py   (daily FCI-G nowcast)
get_spf.py                                            (SPF medians + 25/75 percentiles)
forecast_fcistar.py                                  (one-quarter FCI* forecast + 3×3 sensitivity)
build_fci_nowcast.py                                 (frontend contract + metadata)
```

What's live on https://fcistar.org:

1. **Daily FCI nowcast.** FCI (the FCI-G index) is extended past the last official monthly
   release to the latest market day, as a dotted light-blue line on Figure 1 (and the gap on
   Figure 2). Detailed in the "Daily FCI-G nowcast" section below.
2. **One-quarter-ahead FCI\* forecast.** FCI\* is extended one quarter past the last estimated
   quarter using SPF median GDP-growth / core-PCE forecasts fed through the Kalman filter one
   step with fixed parameters (no re-estimation). Shown as a dark-blue dotted interpolation line
   to a filled dot on Figure 1. Detailed in "One-quarter-ahead FCI\* forecast" below.
3. **Summary boxes** above the charts ("Our \<quarter\> estimate"): FCI, FCI\*, FCI gap at the
   **quarter-end** of the forecast quarter. FCI = the nowcast carried **flat (random walk)** to
   quarter-end; FCI\* = the forecast; gap = their difference. Both lines end in a filled dot at
   quarter-end, so the boxes are exactly the right edge of the charts. See the "DESIGN DECISION"
   section for why we report the quarter-end forecast (not the interpolated daily value or the
   quarter-average), and the gap-interpolation note.
4. **Sensitivity page** (`sensitivity.html`, linked from the summary cards): a 3×3 heatmap of
   FCI\* over the SPF 25/50/75 percentiles of GDP growth × core PCE. See the "Sensitivity grid"
   paragraph below.

**Canonical research copies (FCIstar Dropbox, NOT in this repo):** the FCI-G nowcast scripts are
ported from `code/empirics/fcig/`; `get_spf.py` from `code/empirics/SPF/`. Those are canonical —
change them there first, then port (the ports differ only in I/O paths). `forecast_fcistar.py` and
`build_fci_nowcast.py` are website-native (they depend on the estimation code).

**Monthly cron interaction:** the monthly estimation (`estimate_CCS_struct.py`) commits a
`forecast_inputs/` artifact (`data_hlm.csv`, `covid_dummies.csv`, `filter_settings.json` = exact
`xi_0`/`P_0`) that the daily `forecast_fcistar.py` reads, so the daily forecast runs the filter with
the same fixed parameters **without re-assembling FRED or re-estimating**. Both crons MERGE into
`metadata.json` (never overwrite) so neither clobbers the other's fields — see the 2026-06 fix that
added this after a monthly run wiped the nowcast fields.

**Robustness fixes from this work:** added FRED 429 retry/backoff + inter-series sleep to
`assemble_data.py` (monthly cron had been failing); bumped `actions/checkout@v5` /
`setup-python@v6` (Node 24) across workflows; `metadata.json` merge in both crons.

**Repo secret `FRED_API_KEY`** is set. Used by both the monthly-update and daily-nowcast workflows.

---

## Outstanding Items

- **Done (2026-06):** daily FCI nowcast; one-quarter-ahead FCI\* forecast; quarter-end summary
  boxes; FCI\* sensitivity over SPF percentiles. See the dated status block above.
- **(Natural next extension)** **Multi-quarter FCI\* forecast.** The current forecast is one quarter
  only. Going to two+ quarters needs: (a) SPF forecasts at further horizons (already in
  `spf_forecasts.csv`, horizons 0–4), and (b) the FCI nowcast in the **lagged** regressor — because
  FCI enters the measurement lagged, the step for quarter Q+2 needs `FCI_{Q+1}`, which is not yet
  realized, so the FCI random-walk nowcast finally feeds the filter itself (not just the gap). It's
  a clean layer on top of `one_step()` in `forecast_fcistar.py` (iterate, carrying the appended
  quarter forward). Also think about how to present the widening uncertainty.
- **(Roadmap, not yet scheduled)** Confidence intervals for FCI\* and the FCI gap; mixed-frequency
  display with monthly FCI alongside quarterly FCI\*; eventually FCI-T.
- **(Not for Claude — PI discussion)** Tomas raised a research question about adapting the model to Goldman vs. FCI-G and the levels-vs-differences issue. That's for Alp, Ricardo, and Tomas to work through.

---

## Team and Ownership

| Person | Role |
|---|---|
| Alp Simsek (Yale SOM) | PI; owns the frontend and overall architecture |
| Ricardo Caballero (MIT) | PI |
| Tomas Caravello (MIT) | Estimation lead; wrote the original Matlab code |
| Kentaro Sakata | Yale Econ PhD RA; owns the backend pipeline |

**Related repo (static estimates):** https://github.com/tcaravello/fcistar — Tomas's repo
with static estimates and figures from the paper. The live site replicates and extends this
with auto-updating data and interactive charts.

---

## Folder Structure

```
website/
├── CLAUDE.md               ← you are here
├── backend/
│   ├── refs/
│   │   ├── FCIstarPublic.pdf     ← paper
│   │   └── matlab/               ← most recent stable Matlab code (26_01_10); READ-ONLY
│   │       ├── _aux_fun/         ← auxiliary functions
│   │       ├── assemble_data/    ← data assembly (static, not live)
│   │       ├── construct_fcistar_ygap/
│   │       ├── estimate_CCS/
│   │       ├── estimate_HLW/
│   │       └── quick_play/
│   ├── data/
│   │   ├── raw/                  ← raw downloaded data (gitignored; regenerated by pipeline)
│   │   └── output/               ← estimation outputs consumed by frontend (tracked in git)
│   ├── pipeline/                 ← live Python pipeline
│   │   ├── assemble_data.py      ← pulls FRED + FCI-G + HLW → data/raw/assembled/
│   │   ├── estimate_CCS_struct.py ← Kalman estimation (multistart); writes data/output/
│   │   ├── _aux_fun/             ← Kalman filter helpers
│   │   ├── requirements.txt
│   │   └── README.md
│   └── README.md
├── frontend/
│   ├── refs/
│   │   └── FCIstarPublic.pdf     ← paper (for context when building visualizations)
│   ├── src/                      ← website source code
│   └── README.md
```

---

## Data Flow

```
External data sources
        ↓
backend/pipeline/   [pulls data → runs estimation → writes structured outputs]
        ↓
backend/data/output/   [CSV or JSON files]
        ↓
frontend/src/   [reads output files → renders interactive charts]
```

The backend and frontend are deliberately decoupled. The contract between them is the
file format in `backend/data/output/`. Before building either side, agree on and document
that format in `backend/README.md`.

---

## Backend Status

**The pipeline is live.** Kentaro's Python implementation reproduces the paper's Matlab results with a multistart routine added for robustness. Runs end-to-end in ~2 minutes.

**Entry points:**
- `backend/pipeline/assemble_data.py` — pulls FRED series + FCI-G CSVs from the Fed Board + HLW COVID index from the bundled xlsx, aggregates to a quarterly panel, writes intermediates to `backend/data/raw/assembled/` (`data_hlm.csv`, `covid_dummies.csv`, `data_fred_raw.csv`, `assemble_data_metadata.json`).
- `backend/pipeline/estimate_CCS_struct.py` — loads the intermediates, runs the Kalman filter + smoother via `_aux_fun/`, writes `fcistar.csv`, `metadata.json`, and `theta_opt3.json` to `backend/data/output/`, and archives dated copies of the CSV + metadata into `backend/data/output/vintages/`.

**Key behavior:**
- `--recent-data` defaults to `3`, the live-site spec. In this mode, the sample end is read from the latest available quarter in `data_hlm.csv` (not hardcoded). Other `recent_data` values (0, 1, 2, 4) use the hardcoded windows in `sample_window()` for paper replication.
- `--sample-end-decimal` overrides the auto-detected end if needed (useful for reproducing a historical vintage).
- `--compute-se` (0/1, default 0) toggles bootstrap SE computation.
- The COVID index is set to zero from 2023Q1 onwards, per page 12 of the paper.
- `write_outputs` writes the live frontend files (`fcistar.csv`, `metadata.json`, `theta_opt3.json`) and also archives `vintages/fcistar_YYYY-MM-DD.csv` + `vintages/metadata_YYYY-MM-DD.json`. Vintages are keyed by run date (not sample end), and the function raises `FileExistsError` if the run-date files already exist — so a same-day rerun fails loudly rather than silently overwriting history. To force a same-day rerun, delete the offending vintage files first.

**Running locally:** see "Local dev" section below. Requires `FRED_API_KEY` env var (from `.env` in the repo root, gitignored).

**Do not modify** `backend/refs/matlab/` — it's the read-only Matlab spec used to validate the Python port.

### Daily FCI-G nowcast (extends FCI on Figure 1 to the most recent day)

A second, **daily/weekday** pipeline produces a high-frequency nowcast of FCI-G that extends the FCI
line (and FCI gap) past the last official monthly release. Scripts in `backend/pipeline/`
(`get_fcig.py`, `reproduce_factors.py`, `daily_backcast.py`, `nowcast_fcig.py`) are **PORTED FROM THE
RESEARCH COPY** at `code/empirics/fcig/` in the FCIstar Dropbox project (NOT in this repo). That
research copy is canonical: make nowcast changes there first, then port here. Keep the two in sync —
they should differ only in I/O paths and the output written for the frontend.

Idea: keep the Fed's published FCI-G weights, feed daily public proxies, evaluate with trailing
3-month windows; anchor to each official release. Output (interface contract): `fci_nowcast.csv`
(date, fci, fci_gap, kind ∈ {official, nowcast}) for dates after the last estimated quarter, plus
nowcast fields in `metadata.json`. `fcistar.csv` stays quarterly. A weekday cron runs the daily
pipeline, commits if changed, redeploys.

**Equity source — `^DWCF` primary, `^W5000` pre-1995 backfill (fixed 2026-07-02).** The stock
proxy (`wilshire5000.csv`, from Yahoo's keyless chart API) sets the daily nowcast grid in
`daily_backcast.py` (`days = _load_daily("wilshire5000").index`), so whatever the equity series'
last date is caps the whole nowcast. Yahoo **froze the `^W5000` (FT Wilshire 5000) ticker on
2026-06-26**, which silently stalled the nowcast at 06-26 — the weekday cron kept *succeeding* but
produced no new dates, so "commit if changed" made no commit and the site froze without any error.
Fix: `fetch_yahoo_equity` now uses **`^DWCF` (Dow Jones U.S. Total Stock Market — the exact index
the FCI-G uses) as the primary/live source**, and splices `^W5000` on only for the **pre-1995
history** `^DWCF` lacks (anchored by ratio at `^DWCF`'s 1995 start; only 3-month log-changes enter
the FCI, so the ~1% level offset cancels and the join is seamless). A dead `^W5000` no longer
affects the live edge; falls back to `^DWCF`-only if `^W5000` disappears. Same-index check: daily
log-return corr 0.998, level ratio ~1.00; month-end validation vs official FCI-G corr ~0.998.
Watch-out for the future: if `^DWCF` ever also stalls, the whole nowcast freezes again silently —
consider a cron alert if `nowcast_through` doesn't advance for N weekdays.

### One-quarter-ahead FCI\* forecast (extends FCI\* one quarter past the last estimate)

The same daily cron also extends **FCI\*** one quarter past the last estimated quarter using Survey
of Professional Forecasters (SPF) medians. `get_spf.py` (ported from the canonical research copy
`code/empirics/SPF/`) fetches the Phil Fed median real-GDP-growth and core-PCE forecasts and selects
the most recent survey covering the target quarter (falling back one quarter when the current survey
isn't out yet, so the current quarter always has a forecast). `forecast_fcistar.py` (website-native)
appends that quarter to the panel (`gdp_tot += drgdp/400`, `pce_core` as-is; FCI enters the
measurement *lagged*, so no contemporaneous FCI is needed) and runs the Kalman filter **one step with
the fixed parameters** (`theta_opt3.json` + the committed `forecast_inputs/` = `data_hlm.csv`,
`covid_dummies.csv`, `filter_settings.json` with the exact `xi_0`/`P_0`). No re-estimation — that
stays the monthly job. A reproduction gate asserts the in-sample filtered states match `fcistar.csv`
to ~1e-12. Output: `fcistar_forecast.csv` (target_quarter, date, fcistar, y_gap, drgdp, corepce,
survey_quarter, horizon) + metadata fields `fcistar_forecast`, `fcistar_forecast_through`,
`spf_target_quarter`, `spf_survey_quarter` (null when no valid forecast → gap falls back to
`fcistar_last`). One quarter only; multi-quarter would need the FCI nowcast in the lagged regressor.

**Sensitivity grid.** `get_spf.py` also pulls the SPF cross-sectional dispersion files
(`Dispersion_RGDP.xlsx` sheet D2 = Q/Q-growth 25th/75th percentiles; `Dispersion_COREPCE.xlsx`
sheet D1 = level 25th/75th), so `spf_forecasts.csv` carries `drgdp/corepce` at the 25/50/75
percentiles. `forecast_fcistar.py` then runs the one-step filter over the 3×3 grid (rows = core-PCE
percentile, cols = GDP-growth percentile) and writes `fcistar_sensitivity.json`
(`{target_quarter, survey_quarter, gdp{p25,p50,p75}, corepce{p25,p50,p75}, fcistar[infl][gdp]}`);
the p50/p50 cell equals the headline forecast. `sensitivity.html` + `sensitivity.js` render it as a
heatmap table, linked from the summary cards on the main page. The percentiles are forecaster
disagreement, not an outcome distribution.

### DESIGN DECISION — the nowcast/forecast FCI gap (Figure 2)

The plotted gap in the nowcast/forecast region is **`FCI_nowcast(t) − FCI*_interp(t)`**, where
`FCI*_interp` is a **linear interpolation** between the last quarterly estimate (`fcistar_last`, at
`last_quarter`) and the one-quarter-ahead forecast (`fcistar_forecast`, at `fcistar_forecast_through`),
held flat beyond the forecast quarter-end. This makes Figure 2 exactly the vertical distance between
the FCI and FCI\* lines in Figure 1, and matches how the historical charts already draw FCI\* (quarterly
points connected by straight lines).

This interpolation is done **in the frontend only** — the backend does NOT track an interpolated FCI\*.
`build_fci_nowcast.py` still writes a `fci_gap` column in `fci_nowcast.csv` using a constant anchor
(`fcistar_forecast` held flat), but the **chart recomputes the gap from `fci` + interpolated FCI\***,
so the plotted gap is not the CSV column. The metadata endpoints (`fcistar_last`, `fcistar_forecast`
and their dates) are what the frontend interpolates between.

Known, intentional discrepancy (recorded so we don't rediscover it): this plotted nowcast gap —
**daily FCI nowcast minus interpolated FCI\*** — is a *different object* from the historical quarterly
gap, which is **(quarterly-average FCI) minus (quarterly FCI\*)**. They agree at the quarter-end knots
but differ within a quarter, both because of frequency (daily vs quarter-average FCI) and because FCI\*
is interpolated rather than the single quarterly value. This is deliberate and invisible to users; it
is **not** explained on the methodology page.

---

## Frontend Status

**Stack:** Plain HTML/CSS/JS + **Plotly.js** (CDN, no build step).

**Files in `frontend/src/`:** `index.html`, `main.js`, `style.css`, `CNAME`, plus two extra pages:
`nowcast.html` ("How the FCI nowcast is computed") and `sensitivity.html` + `sensitivity.js`
(the FCI\* sensitivity heatmap). All four HTML/JS files and `style.css` are copied + path-rewritten
by `deploy.yml`.

**Charts (Plotly), all in `index.html` / `main.js`:**
- `chart-fci-and-star` — FCI and FCI\* (Figure 1). Two views toggled by clicking the line or the
  inline pill: line view (FCI + FCI\*, each solid → dotted forecast extension) and decomposition
  view (7 signed stacked bars + black total). FCI extends as a daily nowcast then flat RW to
  quarter-end; FCI\* extends as a dotted interpolation to a forecast dot. Hover is `x unified`;
  the forecast traces are dense (one point per date) with the dots' own hover skipped so each
  value shows once.
- `chart-fci-gap` — FCI gap (Figure 2). Gap = FCI − interpolated FCI\* (recomputed in the
  frontend, `gapOf()`), so it equals the vertical distance between the two Figure 1 lines.
- `chart-ygap` — Output gap (Figure 3), unchanged.

**Summary boxes** (`#summary-cards`, shown when a forecast exists) sit above the charts; a
"Sensitivity analysis" pill under them links to `sensitivity.html`.

**Header:** title, description, authors, link row, then
`Estimates through YYYY QN · Last updated Month YYYY · FCI nowcast through <date> · FCI* forecast for <quarter>`.

**Interactivity:**
- Global date-range selector above the charts: `1Y / 5Y / 10Y / 20Y / All`, rescales all three charts simultaneously via `Plotly.relayout` on `CHART_IDS`.
- Each chart's y-axis recomputes from data within the visible x-window (Plotly's built-in `yaxis.autorange` considers all data, which is why we compute y-range manually in `visibleYRange()`).
- X-axis tick labels adapt to zoom via `tickformatstops`: year-only at long zoom, "Mon YYYY" at short zoom. Hover always shows full quarter-end date via `hoverformat: '%B %-d, %Y'`.

**Data loading:** `main.js` fetches `../../backend/data/output/fcistar.csv` + `metadata.json` relative to the page. The deploy workflow rewrites these to `data/` for production. **Never hardcode estimates.**

**Planned extensions (to land when backend produces the data):**
- Confidence intervals for FCI\* and the FCI gap
- 1–4 quarter forecasts
- Mixed-frequency display: monthly FCI alongside quarterly FCI\*, gap, y\_gap
- Eventually: FCI-T and related series
- "Past vintages" link in the header (vintages archive shipped 2026-04-25; folder fills on the next monthly run)

---

## Hosting and Deployment

**Live site:** https://fcistar.org (custom domain, purchased March 2026, registered for 5 years)

**GitHub repo:** https://github.com/alp-simsek/fcistar

**Hosting:** GitHub Pages from the repo above

**Three workflows, chained:**

1. `.github/workflows/update-data.yml` ("Monthly data update") — monthly cron (`0 18 1 * *`, 18:00 UTC on the 1st — timed after the Fed republishes FCI-G ~14:00 UTC) plus `workflow_dispatch`. Installs deps, runs `assemble_data.py` then `estimate_CCS_struct.py` (with `working-directory: backend/pipeline` so its `from _aux_fun.…` imports resolve), stages `backend/data/output/` (which now includes `forecast_inputs/`), commits as `github-actions[bot]` only if there's a diff, then pushes.

2. `.github/workflows/update-nowcast.yml` ("Daily FCI nowcast") — weekday cron (`0 22 * * 1-5`, 22:00 UTC after the US market close) plus `workflow_dispatch`. Runs the daily pipeline (get_fcig → daily_backcast → nowcast_fcig → get_spf → forecast_fcistar → build_fci_nowcast) in `backend/pipeline`, needs `FRED_API_KEY`, and commits `fci_nowcast.csv, fci_components.csv, fcistar_forecast.csv, fcistar_sensitivity.json, metadata.json` if changed.

3. `.github/workflows/deploy.yml` — assembles `_site/` (copies the frontend HTML/JS + `style.css` + the `backend/data/output/` files into a `data/` subdir, rewrites `../../backend/data/output/` → `data/` via sed on `main.js` / `sensitivity.js` / `index.html`), deploys to GitHub Pages. Triggers on `push: main` **and** on `workflow_run` completion of "Monthly data update" **and** "Daily FCI nowcast" — the `workflow_run` paths are required because pushes made with `GITHUB_TOKEN` don't fire the push trigger (GitHub anti-loop rule). A guard skips deploy if the upstream run failed.

**If the folder structure changes, update the "Assemble site" step in `deploy.yml` accordingly.**

**DNS:** Namecheap → GitHub Pages, configured. HTTPS via Let's Encrypt (GitHub handles automatically).

**Repo secrets:** `FRED_API_KEY` is set at https://github.com/alp-simsek/fcistar/settings/secrets/actions. Rotatable from fredaccount.stlouisfed.org if needed.

---

## Backend-Frontend Interface Contract

The backend delivers files to `backend/data/output/`. The frontend reads them at page load.
**Any change to file formats must be coordinated between Kentaro (backend) and Alp (frontend),
and the frontend JavaScript (main.js) and page descriptions (index.html) updated accordingly.**

**Current files (all in `backend/data/output/`, all tracked in git):**

**`fcistar.csv`** — one row per quarter, columns:
```
date, fci, fcistar, fci_gap, y_gap
```
(dates in YYYY-MM-DD, end-of-quarter; `fci_gap` = FCI − FCI\*)

**`metadata.json`** — run information from the monthly estimation AND the daily nowcast/forecast
(both crons MERGE into this file, never overwrite):
```json
{
  "last_updated": "YYYY-MM-DD", "sample_start": "YYYY-MM-DD", "sample_end": "YYYY-MM-DD",
  "last_quarter": "YYYY-MM-DD", "last_official_date": "YYYY-MM-DD", "nowcast_through": "YYYY-MM-DD",
  "fcistar_last": -0.78,
  "fcistar_forecast": -0.77, "fcistar_forecast_through": "YYYY-MM-DD",
  "spf_target_quarter": "YYYYQn", "spf_survey_quarter": "YYYYQn"
}
```
`sample_end` → "Estimates through YYYY QN"; `last_updated` → "Last updated Month YYYY";
`nowcast_through` → "FCI nowcast through <date>"; `spf_target_quarter` → "FCI\* forecast for <quarter>".
The four `fcistar_forecast*` / `spf_*` keys are **null** when there is no valid forecast (then the gap
falls back to `fcistar_last`).

**`fci_nowcast.csv`** — `date, fci, fci_gap, kind` (kind ∈ {official, nowcast}) for dates after the
last estimated quarter (recent monthly official FCI-G, then daily nowcast). Written by
`build_fci_nowcast.py`. The frontend recomputes the displayed gap via interpolated FCI\*, so the
`fci_gap` column here is informational (constant-anchor).

**`fci_components.csv`** — `date, kind, ffr, t10yr, mortgage, bbb, stock, house, dollar` (kind ∈
{quarter, official, nowcast}) — the 7-factor FCI decomposition for Figure 1's decomposition view.

**`fcistar_forecast.csv`** — one row: `target_quarter, date, fcistar, y_gap, drgdp, corepce,
survey_quarter, horizon` — the one-quarter-ahead FCI\* forecast (median SPF). Written by
`forecast_fcistar.py`.

**`fcistar_sensitivity.json`** — the 3×3 FCI\* sensitivity grid:
`{target_quarter, survey_quarter, gdp{p25,p50,p75}, corepce{p25,p50,p75}, fcistar[inflPct][gdpPct]}`.
Written by `forecast_fcistar.py`; rendered by `sensitivity.html`/`sensitivity.js`.

**`forecast_inputs/`** — `data_hlm.csv`, `covid_dummies.csv`, `filter_settings.json` (exact
`xi_0`/`P_0` + spec). Committed by the MONTHLY estimation; read by the daily `forecast_fcistar.py`
so it can run the one-step filter with fixed parameters without re-assembling/re-estimating.

**`theta_opt3.json`** — fitted parameter vector from the Kalman estimation. Read by
`forecast_fcistar.py` (the fixed θ); kept in git for reproducibility/debugging.

**`fci_star_results.xlsx`** — Excel bundle written by the estimation. Gitignored (pattern `backend/data/output/*.xlsx`), regenerated each run.

**`vintages/fcistar_YYYY-MM-DD.csv`** and **`vintages/metadata_YYYY-MM-DD.json`** — dated archive of every pipeline run. Same schemas as the live files above. Keyed by the pipeline **run date**, not the sample end (so two runs with the same sample end produce two separate vintages, capturing data revisions). Never overwritten — the pipeline raises `FileExistsError` on a same-day rerun.

**Planned extensions to the output files:**

- Confidence intervals: add `fcistar_upper`, `fcistar_lower`, `fci_gap_upper`,
  `fci_gap_lower` columns to `fcistar.csv`
- Forecasts: add columns `fci_gap_f1` through `fci_gap_f4` (1–4 quarter ahead forecasts)
- Mixed frequency: consider splitting into two files —
  `fci_monthly.csv` (date, fci at monthly frequency) and
  `fcistar_quarterly.csv` (date, fcistar, fci_gap, y_gap at quarterly frequency).
  The frontend will plot them on the same chart with frequency clearly labeled.

These extensions require coordinated changes to both the backend pipeline and the
frontend charts. Do not add columns silently — update this contract and main.js together.

---

## Key Conventions

- `backend/data/raw/` is never hand-edited — always regenerated by the pipeline (gitignored)
- `backend/data/output/` is the single source of truth for the frontend
- Document all external data sources (series names, FRED codes, vintage) in `backend/README.md`
- The output file format is the interface contract — any change must be coordinated between backend and frontend
- Commit output files to the repo so the frontend can be developed and tested without running the full pipeline

---

## Local Dev / Picking Up Work

**Venv and deps** (first time, or fresh clone):
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r backend/pipeline/requirements.txt
```

**FRED API key:** stored in `.env` at repo root (gitignored). Load before running the pipeline:
```bash
set -a && source .env && set +a
```

**Run the pipeline locally:**
```bash
source venv/bin/activate
set -a && source .env && set +a
python backend/pipeline/assemble_data.py
cd backend/pipeline && python estimate_CCS_struct.py
```
`cd` is needed for the estimation step because its `from _aux_fun.…` imports require the script's directory on `sys.path`. Total runtime ~2 minutes on a laptop.

**Preview frontend locally:** `file://` URLs won't work because browsers block `fetch()` against local files. Serve over HTTP:
```bash
python3 -m http.server 8000  # run from repo root
```
Then open `http://localhost:8000/frontend/src/`. The relative data paths (`../../backend/data/output/…`) resolve against the HTTP root.

**After deploys — browser cache gotcha:** HTML and JS cache independently. If a site change is deployed and you see the new markup but old JS behavior, the browser is serving a cached `main.js`. Test in an incognito window or DevTools → "Empty Cache and Hard Reload" to confirm.

**Working from a new machine:** clone with an auth token in the URL, then set identity:
```bash
git clone https://alp-simsek:YOUR_TOKEN@github.com/alp-simsek/fcistar.git
git config --global user.name  "Alp Simsek"
git config --global user.email "alp.simsek@yale.edu"
```

**Manual pipeline trigger:** https://github.com/alp-simsek/fcistar/actions/workflows/update-data.yml → "Run workflow". Useful for forcing a fresh estimation without waiting for the monthly cron.
