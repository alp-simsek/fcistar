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

## Team and Ownership

| Person | Role |
|---|---|
| Alp Simsek (Yale SOM) | PI; owns the frontend and overall architecture |
| Ricardo Caballero (MIT) | PI |
| Tomas Caravello (MIT) | Estimation lead; wrote the original Matlab code |
| Kentaro [last name] | Yale Econ PhD RA; owns the backend pipeline |

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
│   │   ├── raw/                  ← raw downloaded data (never hand-edit)
│   │   └── output/               ← estimation outputs consumed by frontend (CSV or JSON)
│   ├── pipeline/                 ← scripts that pull data, run estimation, write output/
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

## Backend Instructions (Kentaro)

**Goal:** A fully automated pipeline that requires no manual intervention after setup.

**Steps:**
1. Pull the most recent data from external sources (document all sources and series IDs in `backend/README.md`)
2. Run the estimation (Kalman filter / state-space model — see paper and `refs/estimation.m`)
3. Write structured output to `backend/data/output/` in an agreed format (CSV or JSON)

**Language:** Python preferred for the pipeline (automation, scheduling, FRED API access).
The reference Matlab code in `refs/matlab/` (version `26_01_10`) is the stable implementation
that reproduces the paper's figures using static data. Kentaro's job is to build a new
Python pipeline that replicates this logic but pulls live data automatically. Treat the
Matlab code as the specification — understand what each module does, then reimplement:

- `assemble_data/` → data ingestion (replace with live API pulls)
- `estimate_CCS/` and `estimate_HLW/` → Kalman filter estimation (translate to Python)
- `construct_fcistar_ygap/` → FCI* and output gap construction
- `_aux_fun/` → helper functions needed across modules
- `quick_play/` → scratch/exploratory; lower priority

**Key reference:** Laubach-Williams at the NY Fed does something similar. Their replication
code and data release schedule are at the URL above — useful for understanding the pipeline
pattern we want to match.

---

## Frontend Instructions (Alp)

**Goal:** An interactive, clean public-facing website that auto-updates when backend outputs change.

**Design reference:** NY Fed r-star page — but with interactive charts (not static images),
a cleaner layout, and room to add series over time.

**Stack:** Plain HTML/CSS/JS + **Plotly.js** for charts. No build step, no framework.
Deployable as a static site via GitHub Pages.

**Charts should show at minimum:**
- FCI\* over time
- FCI over time
- FCI gap (FCI minus FCI\*) — the policy stance indicator
- Eventually: FCI-T and other related series

**Page structure (modeled on NY Fed r-star):**
- Header: title, brief description, last-updated date
- Chart panels (interactive, zoomable, with hover tooltips)
- "Download data" link pointing to the CSV in `backend/data/output/`
- "About" section linking to the paper

**Data source:** Read from `backend/data/output/` at page load — do not hardcode any estimates.

---

## Hosting and Deployment

**Live site:** https://fcistar.org (custom domain, purchased March 2026, registered for 5 years)

**GitHub repo:** https://github.com/alp-simsek/fcistar

**Hosting:** GitHub Pages from the repo above

**Auto-update mechanism:** GitHub Actions cron job runs Kentaro's Python pipeline on a
schedule (e.g., weekly or quarterly), commits updated output files to the repo, and GitHub
Pages redeploys automatically.

**DNS:** Namecheap → GitHub Pages (to be configured once the site is ready to go live).
GitHub handles HTTPS automatically via Let's Encrypt.

**Working from a new machine:** Clone with:
`git clone https://alp-simsek:YOUR_TOKEN@github.com/alp-simsek/fcistar.git`
Then set identity: `git config --global user.name "Alp Simsek"` and `git config --global user.email "alp.simsek@yale.edu"`.
Token is stored in the remote URL in `.git/config` — valid from any machine.

---

## Backend-Frontend Interface Contract

The backend delivers exactly two files to `backend/data/output/`:

**`fcistar.csv`** — one row per time period, columns:
```
date, fci, fcistar, fci_gap, fcistar_upper, fcistar_lower
```
(dates in YYYY-MM-DD format; confidence bands for FCI\*)

**`metadata.json`** — pipeline run information:
```json
{
  "last_updated": "YYYY-MM-DD",
  "sample_start": "YYYY-MM-DD",
  "sample_end": "YYYY-MM-DD"
}
```

The `last_updated` field is displayed prominently on the site — if stale, it signals a
pipeline failure. Any change to this format must be coordinated between Kentaro and Alp.

---

## Key Conventions

- `backend/data/raw/` is never hand-edited — always regenerated by the pipeline
- `backend/data/output/` is the single source of truth for the frontend
- Document all external data sources (series names, FRED codes, vintage) in `backend/README.md`
- The output file format is the interface contract — any change must be coordinated between backend and frontend
- Commit output files to the repo so the frontend can be developed and tested without running the full pipeline
