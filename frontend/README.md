# Frontend

## Overview
Static public-facing website for FCI*. Reads estimation outputs from `backend/data/output/`
and renders interactive charts. Live at https://fcistar.org.
Also accessible at https://alp-simsek.github.io/fcistar/

## Stack
- Plain HTML/CSS/JS — no build step, no framework
- Plotly.js (loaded from CDN) for interactive charts
- Deployed via GitHub Actions → GitHub Pages (see `.github/workflows/deploy.yml`)

## Folder structure
- `refs/` — paper PDF for reference
- `src/` — website source files:
  - `index.html` — page structure and text
  - `style.css` — styling (navy/gold header, Georgia serif, chart layout)
  - `main.js` — data loading, CSV parsing, Plotly chart rendering
  - `CNAME` — custom domain declaration (fcistar.org)

## Data source
Reads `backend/data/output/fcistar.csv` and `metadata.json` at page load.
Do not hardcode any estimates. Data paths are rewritten by the deploy workflow
from `../../backend/data/output/` (local dev) to `data/` (production).

## Charts
1. **FCI and FCI\*** — the main chart. FCI-G baseline index (3-year lookback, Ajello et al. 2023)
   and our estimate of its neutral level. Units: percentage points of next-year GDP growth.
2. **Financial Conditions Gap (FCI − FCI\*)** — policy stance indicator. Positive = tighter than
   neutral (restraining output); negative = looser than neutral (pushing output above potential).
3. **Output gap** — real-time estimate of GDP relative to potential, in percent.

## Page structure
- Header: title, two-paragraph description, author names, links (paper, data download),
  sample end date ("Estimates through YYYY QN")
- Three interactive chart panels (zoomable, hover tooltips, downloadable as PNG)
- Minimal footer with copyright

## Local development
Run a local server from the repo root:
```
python3 -m http.server 8000
```
Then open `http://localhost:8000/frontend/src/`
