/* =============================================================
   FCI* Website — main.js
   Loads fcistar.csv (quarterly), fci_nowcast.csv (recent monthly
   official + daily nowcast extension), and metadata.json, then
   draws three interactive Plotly charts.

   Figure 1 (FCI & FCI*) and Figure 2 (FCI gap) extend FCI past the
   last estimated quarter: a SOLID line through the last official
   monthly FCI-G release, then a DOTTED line (same color) for the
   daily nowcast of the not-yet-released dates. The gap extension
   holds FCI* at its latest estimate (metadata.fcistar_last).
   ============================================================= */

const DATA_URL    = '../../backend/data/output/fcistar.csv';
const NOWCAST_URL = '../../backend/data/output/fci_nowcast.csv';
const META_URL    = '../../backend/data/output/metadata.json';

const CHART_IDS = ['chart-fci-and-star', 'chart-fci-gap', 'chart-ygap'];

// Color palette
const COLORS = {
  fci:     '#7a9cc4',  // steel blue — input series
  fcistar: '#1a3a6b',  // navy      — our main contribution
  gap:     '#c0392b',  // red       — FCI gap
  ygap:    '#1a3a6b',  // navy      — output gap
};

// Per-chart series (x/y arrays of every plotted point), used to rescale
// the y-axis to the data within the visible x-window. Set in drawCharts.
let CHART_SERIES = {};
let END_STR = null;   // latest date across all series (drives the x-range end)


/* ---------------- CSV PARSER ---------------- */
function parseCSV(text) {
  const lines = text.trim().split('\n');
  const headers = lines[0].split(',').map(h => h.trim());
  return lines.slice(1).map(line => {
    const vals = line.split(',');
    const row = {};
    headers.forEach((h, i) => {
      row[h] = (h === 'date' || h === 'kind') ? vals[i].trim() : parseFloat(vals[i]);
    });
    return row;
  });
}


/* ---------------- DATE FORMATTERS ---------------- */
function formatQuarter(iso) {
  const month = parseInt(iso.slice(5, 7), 10);
  return `${iso.slice(0, 4)} Q${Math.ceil(month / 3)}`;
}
const MONTHS = ['January','February','March','April','May','June',
                'July','August','September','October','November','December'];
function formatMonthYear(iso) {
  return `${MONTHS[parseInt(iso.slice(5, 7), 10) - 1]} ${iso.slice(0, 4)}`;
}
function formatLongDate(iso) {
  return `${MONTHS[parseInt(iso.slice(5, 7), 10) - 1]} ${parseInt(iso.slice(8, 10), 10)}, ${iso.slice(0, 4)}`;
}


/* ---------------- LAYOUT HELPERS ---------------- */
function zeroLine() {
  return {
    type: 'line',
    xref: 'paper', x0: 0, x1: 1,
    yref: 'y',     y0: 0, y1: 0,
    line: { color: '#aaa', width: 1, dash: 'dash' }
  };
}

function baseLayout(yAxisTitle) {
  return {
    margin: { t: 24, r: 24, b: 48, l: 60 },
    paper_bgcolor: '#fff',
    plot_bgcolor:  '#fff',
    font: { family: 'Helvetica Neue, Arial, sans-serif', size: 12, color: '#444' },
    xaxis: {
      type: 'date',
      gridcolor: '#eee',
      linecolor: '#ccc',
      tickformatstops: [
        { dtickrange: [null, 'M6'], value: '%b %Y' },
        { dtickrange: ['M6', null], value: '%Y' },
      ],
      hoverformat: '%B %-d, %Y',
    },
    yaxis: {
      title: { text: yAxisTitle, standoff: 8 },
      zeroline: false,
      gridcolor: '#eee',
      linecolor: '#ccc',
    },
    legend: { orientation: 'h', x: 0, y: 1.08, font: { size: 12 }, bgcolor: 'rgba(0,0,0,0)' },
    hovermode: 'x unified',
  };
}


/* ---------------- DRAW CHARTS ---------------- */
function drawCharts(q, nc) {
  const col = (arr, c) => arr.map(d => d[c]);
  const official = nc.filter(d => d.kind === 'official');   // recent monthly official points
  const nowcast  = nc.filter(d => d.kind === 'nowcast');    // daily not-yet-released nowcast

  // last point of the SOLID line = last official month if any, else the last quarter
  const lastSolid = official.length ? official[official.length - 1] : q[q.length - 1];

  // FCI: solid = quarterly history + monthly official; dotted = lastSolid + daily nowcast
  const fciSolidX = col(q, 'date').concat(col(official, 'date'));
  const fciSolidY = col(q, 'fci').concat(col(official, 'fci'));
  const fciDotX   = [lastSolid.date].concat(col(nowcast, 'date'));
  const fciDotY   = [lastSolid.fci].concat(col(nowcast, 'fci'));

  // Gap: same split. Historical gap uses actual FCI*; the extension holds FCI* at its last value.
  const lastSolidGap = lastSolid.fci_gap;
  const gapSolidX = col(q, 'date').concat(col(official, 'date'));
  const gapSolidY = col(q, 'fci_gap').concat(col(official, 'fci_gap'));
  const gapDotX   = [lastSolid.date].concat(col(nowcast, 'date'));
  const gapDotY   = [lastSolidGap].concat(col(nowcast, 'fci_gap'));

  const config = { responsive: true, displayModeBar: 'hover' };
  const dot = (color) => ({ color, width: 2, dash: 'dot' });

  // --- Chart 1: FCI and FCI* ---
  Plotly.newPlot('chart-fci-and-star', [
    { x: fciSolidX, y: fciSolidY, name: 'FCI', type: 'scatter', mode: 'lines',
      line: { color: COLORS.fci, width: 2 }, hovertemplate: '%{y:.2f}<extra>FCI</extra>' },
    { x: fciDotX, y: fciDotY, name: 'FCI (nowcast)', type: 'scatter', mode: 'lines', showlegend: false,
      line: dot(COLORS.fci), hovertemplate: '%{y:.2f}<extra>FCI · nowcast</extra>' },
    { x: col(q, 'date'), y: col(q, 'fcistar'), name: 'FCI*', type: 'scatter', mode: 'lines',
      line: { color: COLORS.fcistar, width: 2.5 }, hovertemplate: '%{y:.2f}<extra>FCI*</extra>' },
  ], { ...baseLayout('Pct. pts. of next-year GDP growth'), shapes: [zeroLine()] }, config);

  // --- Chart 2: FCI gap ---
  Plotly.newPlot('chart-fci-gap', [
    { x: gapSolidX, y: gapSolidY, name: 'FCI gap', type: 'scatter', mode: 'lines',
      line: { color: COLORS.gap, width: 2 }, hovertemplate: '%{y:.2f}<extra>FCI gap</extra>' },
    { x: gapDotX, y: gapDotY, name: 'FCI gap (nowcast)', type: 'scatter', mode: 'lines', showlegend: false,
      line: dot(COLORS.gap), hovertemplate: '%{y:.2f}<extra>FCI gap · nowcast</extra>' },
  ], { ...baseLayout('Pct. pts. of next-year GDP growth'), shapes: [zeroLine()] }, config);

  // --- Chart 3: Output gap (unchanged, quarterly) ---
  Plotly.newPlot('chart-ygap', [
    { x: col(q, 'date'), y: col(q, 'y_gap'), name: 'Output gap', type: 'scatter', mode: 'lines',
      line: { color: COLORS.ygap, width: 2 }, hovertemplate: '%{y:.2f}%<extra>Output gap</extra>' },
  ], { ...baseLayout('Percent'), shapes: [zeroLine()] }, config);

  // record series for y-rescaling, and the global x-range end
  CHART_SERIES = {
    'chart-fci-and-star': [ { x: fciSolidX, y: fciSolidY }, { x: fciDotX, y: fciDotY }, { x: col(q, 'date'), y: col(q, 'fcistar') } ],
    'chart-fci-gap':      [ { x: gapSolidX, y: gapSolidY }, { x: gapDotX, y: gapDotY } ],
    'chart-ygap':         [ { x: col(q, 'date'), y: col(q, 'y_gap') } ],
  };
  END_STR = (nowcast.length ? nowcast : official.length ? official : q).slice(-1)[0].date;
}


/* ---------------- RANGE BUTTONS ---------------- */
// [min,max] of all series' values within [xStart,xEnd], 5% padded; null if empty.
function visibleYRange(seriesList, xStart, xEnd) {
  let lo = Infinity, hi = -Infinity;
  for (const s of seriesList) {
    for (let i = 0; i < s.x.length; i++) {
      const d = new Date(s.x[i]);
      if (d < xStart || d > xEnd) continue;
      const v = s.y[i];
      if (v == null || isNaN(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  if (!isFinite(lo)) return null;
  const pad = (hi - lo) * 0.05 || 0.1;
  return [lo - pad, hi + pad];
}

function applyRange(years) {
  const end = new Date(END_STR);
  let start, startStr;
  if (years === 'all') {
    start = new Date('1990-01-01');
  } else {
    start = new Date(end);
    start.setUTCFullYear(start.getUTCFullYear() - years);
    startStr = start.toISOString().slice(0, 10);
  }
  CHART_IDS.forEach(id => {
    const yrange = visibleYRange(CHART_SERIES[id], start, end);
    const update = { 'yaxis.autorange': false };
    if (years === 'all') {
      update['xaxis.autorange'] = true;
    } else {
      update['xaxis.range'] = [startStr, END_STR];
    }
    if (yrange) update['yaxis.range'] = yrange;
    Plotly.relayout(id, update);
  });
}

function initRangeButtons() {
  document.querySelectorAll('.range-buttons button').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('.range-buttons button').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      const y = btn.dataset.years;
      applyRange(y === 'all' ? 'all' : parseInt(y, 10));
    });
  });
}


/* ---------------- MAIN ENTRY ---------------- */
Promise.all([
  fetch(DATA_URL).then(r => r.text()),
  fetch(NOWCAST_URL).then(r => r.text()),
  fetch(META_URL).then(r => r.json()),
])
.then(([csvText, ncText, meta]) => {
  document.getElementById('last-updated-date').textContent   = formatQuarter(meta.sample_end);
  document.getElementById('last-refreshed-date').textContent = formatMonthYear(meta.last_updated);
  const nowcastEl = document.getElementById('nowcast-through-date');
  if (nowcastEl && meta.nowcast_through) nowcastEl.textContent = formatLongDate(meta.nowcast_through);

  const q  = parseCSV(csvText);
  const nc = parseCSV(ncText);
  drawCharts(q, nc);
  initRangeButtons();
})
.catch(err => {
  console.error('Failed to load data:', err);
  document.getElementById('last-updated-date').textContent   = 'unavailable';
  document.getElementById('last-refreshed-date').textContent = 'unavailable';
});
