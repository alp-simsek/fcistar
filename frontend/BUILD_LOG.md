# Frontend Build Log

## Step 1 — HTML skeleton (2026-03-24)

Created `src/index.html` with the basic page structure:
- Header section: title, subtitle, last-updated date
- Three chart containers (div elements) for Plotly to draw into
- About/download section at the bottom
- Linked Plotly.js from CDN (no install needed)
- Linked `style.css` and `main.js` (to be created next)

**Key concept:** The HTML file is just a skeleton. The chart `<div>` elements are empty
boxes — JavaScript will fill them with charts once the page loads.

## Step 2 — CSS styling (2026-03-24)

Created `src/style.css` with sections:
1. Reset: zero out browser defaults for consistency
2. Layout: `.container` centers content, caps width at 900px
3. `.star`: superscript styling for the * in FCI* (line-height:0 prevents line-height blowout)
4. Header: navy blue (#1a3a6b) with gold accent line
5. Chart sections: spacing, h2 styling, description text, explicit 380px height for Plotly
6. Footer: light gray background, muted text

**Key concept:** CSS rules are "selector { property: value }". Browsers apply their own
defaults before your CSS — the reset section overrides the worst of these for consistency.

## Step 3 — JavaScript: data loading and charts (2026-03-24)

Created `src/main.js`:
- `parseCSV()`: manual CSV parser (no library needed for our simple format)
- `formatDate()`: converts "2026-03-24" → "March 2026" for the header
- `zeroLine()`: Plotly shape for the dashed zero reference line on each chart
- `baseLayout()`: shared Plotly layout defaults (colors, font, axes, hover mode)
- `drawCharts()`: builds all three Plotly charts from parsed data arrays
- `Promise.all([fetch(CSV), fetch(JSON)])`: loads both files in parallel, then renders

**Key concept:** `fetch()` is blocked on file:// protocol. To develop locally, run:
  `python3 -m http.server 8000` from repo root, then open
  `http://localhost:8000/frontend/src/`

**Data path:** `../../backend/data/output/` (relative to frontend/src/)
