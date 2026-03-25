# Frontend

## Overview
Static public-facing website for FCI*. Reads estimation outputs from `backend/data/output/` and renders interactive charts. Deployed via GitHub Pages at https://fcistar.org.

## Stack
- Plain HTML/CSS/JS — no build step, no framework
- Plotly.js for interactive charts
- Deployable as a static site

## Folder structure
- `refs/` — paper PDF for reference when building visualizations
- `src/` — website source code

## Data source
Reads `../backend/data/output/fcistar.csv` and `metadata.json` at page load. Do not hardcode any estimates.

## Charts (minimum)
- FCI* over time (with confidence band)
- FCI over time
- FCI gap (FCI minus FCI*) — the policy stance indicator

## Page structure
- Header: title, brief description, last-updated date
- Chart panels (interactive, zoomable, hover tooltips)
- Download data link → `backend/data/output/fcistar.csv`
- About section linking to the paper
