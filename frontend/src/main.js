/* =============================================================
   FCI* Website — main.js
   Loads fcistar.csv (quarterly), fci_nowcast.csv (recent monthly
   official + daily nowcast), fci_components.csv (7-factor
   decomposition), and metadata.json, then draws three charts.

   Figure 1 has two views, toggled by the Line/Decomposition control
   or by clicking the FCI line:
     - Line:          FCI (solid official -> dotted nowcast) + FCI*.
     - Decomposition: the 7 contributions as signed stacked bars
                      (nowcast = hatched) + black Total line.
   Figures 2 (FCI gap) and 3 (output gap) are unchanged; the gap is
   extended with the nowcast holding FCI* at metadata.fcistar_last.
   ============================================================= */

const DATA_URL    = '../../backend/data/output/fcistar.csv';
const NOWCAST_URL = '../../backend/data/output/fci_nowcast.csv';
const COMP_URL    = '../../backend/data/output/fci_components.csv';
const META_URL    = '../../backend/data/output/metadata.json';

const CHART_IDS = ['chart-fci-and-star', 'chart-fci-gap', 'chart-ygap'];
const FCI_CHART = 'chart-fci-and-star';
const DAY = 86400000;

const COLORS = { fci: '#7a9cc4', fcistar: '#1a3a6b', gap: '#c0392b', ygap: '#1a3a6b' };

// Decomposition components (order + colors match the paper figure / MATLAB scheme)
const COMP = [
  { key: 'ffr',      name: 'FFR',           color: 'rgb(0,114,189)'  },
  { key: 't10yr',    name: '10Y Rate',      color: 'rgb(217,83,25)'  },
  { key: 'mortgage', name: 'Mortgage Rate', color: 'rgb(237,177,32)' },
  { key: 'bbb',      name: 'BBB',           color: 'rgb(126,47,142)' },
  { key: 'stock',    name: 'Equity',        color: 'rgb(119,172,48)' },
  { key: 'house',    name: 'Housing',       color: 'rgb(77,190,238)' },
  { key: 'dollar',   name: 'Dollar',        color: 'rgb(162,20,47)'  },
];

const PLOT_CONFIG = { responsive: true, displayModeBar: 'hover' };

let Q, NC, COMP_SOLID, COMP_NOW;   // parsed data
let CHART_SERIES = {};             // line-mode y-range series per chart
let END_STR = null;                // latest date across series
let RANGE = 'all';                 // current x-range selection
let decompMode = false;            // Figure 1 view


/* ---------------- CSV PARSER ---------------- */
function parseCSV(text) {
  const lines = text.trim().split('\n');
  const headers = lines[0].split(',').map(h => h.trim());
  return lines.slice(1).map(line => {
    const vals = line.split(',');
    const row = {};
    headers.forEach((h, i) => { row[h] = (h === 'date' || h === 'kind') ? vals[i].trim() : parseFloat(vals[i]); });
    return row;
  });
}


/* ---------------- DATE FORMATTERS ---------------- */
function formatQuarter(iso) {
  return `${iso.slice(0, 4)} Q${Math.ceil(parseInt(iso.slice(5, 7), 10) / 3)}`;
}
const MONTHS = ['January','February','March','April','May','June','July','August','September','October','November','December'];
function formatMonthYear(iso) { return `${MONTHS[parseInt(iso.slice(5, 7), 10) - 1]} ${iso.slice(0, 4)}`; }
function formatLongDate(iso) { return `${MONTHS[parseInt(iso.slice(5, 7), 10) - 1]} ${parseInt(iso.slice(8, 10), 10)}, ${iso.slice(0, 4)}`; }


/* ---------------- LAYOUT HELPERS ---------------- */
function zeroLine() {
  return { type: 'line', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: 0, y1: 0,
           line: { color: '#aaa', width: 1, dash: 'dash' } };
}
function baseLayout(yAxisTitle) {
  return {
    margin: { t: 24, r: 24, b: 48, l: 60 },
    paper_bgcolor: '#fff', plot_bgcolor: '#fff',
    font: { family: 'Helvetica Neue, Arial, sans-serif', size: 12, color: '#444' },
    xaxis: { type: 'date', gridcolor: '#eee', linecolor: '#ccc',
      tickformatstops: [ { dtickrange: [null, 'M6'], value: '%b %Y' }, { dtickrange: ['M6', null], value: '%Y' } ],
      hoverformat: '%B %-d, %Y' },
    yaxis: { title: { text: yAxisTitle, standoff: 8 }, zeroline: false, gridcolor: '#eee', linecolor: '#ccc' },
    legend: { orientation: 'h', x: 0, y: 1.08, font: { size: 12 }, bgcolor: 'rgba(0,0,0,0)' },
    hovermode: 'x unified',
  };
}


/* ---------------- FCI LINE ARRAYS (shared by line view + decomposition total) ---------------- */
function fciLineArrays() {
  const col = (a, c) => a.map(d => d[c]);
  const official = NC.filter(d => d.kind === 'official');
  const nowc = NC.filter(d => d.kind === 'nowcast');
  const lastSolid = official.length ? official[official.length - 1] : Q[Q.length - 1];
  return {
    solidX: col(Q, 'date').concat(col(official, 'date')),
    fciSolid: col(Q, 'fci').concat(col(official, 'fci')),
    gapSolid: col(Q, 'fci_gap').concat(col(official, 'fci_gap')),
    dotX: [lastSolid.date].concat(col(nowc, 'date')),
    fciDot: [lastSolid.fci].concat(col(nowc, 'fci')),
    gapDot: [lastSolid.fci_gap].concat(col(nowc, 'fci_gap')),
  };
}


/* ---------------- FIGURE 1 TRACES ---------------- */
function chart1LineTraces() {
  const L = fciLineArrays();
  return [
    { x: L.solidX, y: L.fciSolid, name: 'FCI', type: 'scatter', mode: 'lines',
      line: { color: COLORS.fci, width: 2 }, hovertemplate: '%{y:.2f}<extra>FCI</extra>' },
    { x: L.dotX, y: L.fciDot, name: 'FCI (nowcast)', type: 'scatter', mode: 'lines', showlegend: false,
      line: { color: COLORS.fci, width: 2, dash: 'dot' }, hovertemplate: '%{y:.2f}<extra>FCI · nowcast</extra>' },
    { x: Q.map(d => d.date), y: Q.map(d => d.fcistar), name: 'FCI*', type: 'scatter', mode: 'lines',
      line: { color: COLORS.fcistar, width: 2.5 }, hovertemplate: '%{y:.2f}<extra>FCI*</extra>' },
  ];
}
function chart1DecompTraces() {
  const sd = COMP_SOLID.map(r => r.date), nd = COMP_NOW.map(r => r.date);
  // Each bar spans the period ENDING at its date (right-aligned: offset = -width), so a quarter
  // bar covers its 3 months, the monthly official bar covers its month, and daily nowcast bars
  // their day — none straddles into the next period (which would overlap the hatched nowcast).
  const wSolid = COMP_SOLID.map(r => (r.kind === 'quarter' ? 88 : 28) * DAY);
  const oSolid = wSolid.map(w => -w);
  const wNow = 1 * DAY;
  const traces = [];
  COMP.forEach(c => {
    traces.push({ x: sd, y: COMP_SOLID.map(r => r[c.key]), name: c.name, type: 'bar',
      marker: { color: c.color }, width: wSolid, offset: oSolid, legendgroup: c.key,
      hovertemplate: '%{y:.2f}<extra>' + c.name + '</extra>' });
    traces.push({ x: nd, y: COMP_NOW.map(r => r[c.key]), type: 'bar', showlegend: false, legendgroup: c.key,
      marker: { color: c.color, pattern: { shape: '/', size: 4, solidity: 0.35 } }, width: wNow, offset: -wNow,
      hovertemplate: '%{y:.2f}<extra>' + c.name + ' · nowcast</extra>' });
  });
  const L = fciLineArrays();
  traces.push({ x: L.solidX, y: L.fciSolid, name: 'Total', type: 'scatter', mode: 'lines',
    line: { color: 'black', width: 2.5 }, hovertemplate: '%{y:.2f}<extra>FCI</extra>' });
  traces.push({ x: L.dotX, y: L.fciDot, type: 'scatter', mode: 'lines', showlegend: false,
    line: { color: 'black', width: 2.5, dash: 'dot' }, hovertemplate: '%{y:.2f}<extra>FCI · nowcast</extra>' });
  return traces;
}
function drawChart1() {
  const layout = { ...baseLayout('Pct. pts. of next-year GDP growth'), shapes: [zeroLine()] };
  if (decompMode) layout.barmode = 'relative';
  Plotly.react(FCI_CHART, decompMode ? chart1DecompTraces() : chart1LineTraces(), layout, PLOT_CONFIG);
}


/* ---------------- DRAW ALL CHARTS ---------------- */
function drawCharts() {
  const L = fciLineArrays();
  drawChart1();

  Plotly.newPlot('chart-fci-gap', [
    { x: L.solidX, y: L.gapSolid, name: 'FCI gap', type: 'scatter', mode: 'lines',
      line: { color: COLORS.gap, width: 2 }, hovertemplate: '%{y:.2f}<extra>FCI gap</extra>' },
    { x: L.dotX, y: L.gapDot, name: 'FCI gap (nowcast)', type: 'scatter', mode: 'lines', showlegend: false,
      line: { color: COLORS.gap, width: 2, dash: 'dot' }, hovertemplate: '%{y:.2f}<extra>FCI gap · nowcast</extra>' },
  ], { ...baseLayout('Pct. pts. of next-year GDP growth'), shapes: [zeroLine()] }, PLOT_CONFIG);

  Plotly.newPlot('chart-ygap', [
    { x: Q.map(d => d.date), y: Q.map(d => d.y_gap), name: 'Output gap', type: 'scatter', mode: 'lines',
      line: { color: COLORS.ygap, width: 2 }, hovertemplate: '%{y:.2f}%<extra>Output gap</extra>' },
  ], { ...baseLayout('Percent'), shapes: [zeroLine()] }, PLOT_CONFIG);

  CHART_SERIES = {
    [FCI_CHART]: [ { x: L.solidX, y: L.fciSolid }, { x: L.dotX, y: L.fciDot }, { x: Q.map(d => d.date), y: Q.map(d => d.fcistar) } ],
    'chart-fci-gap': [ { x: L.solidX, y: L.gapSolid }, { x: L.dotX, y: L.gapDot } ],
    'chart-ygap': [ { x: Q.map(d => d.date), y: Q.map(d => d.y_gap) } ],
  };
  END_STR = L.dotX[L.dotX.length - 1];

  // clicking the FCI line (or any point in Figure 1) toggles the decomposition view
  document.getElementById(FCI_CHART).on('plotly_click', () => setView(decompMode ? 'line' : 'decomp'));
}


/* ---------------- RANGE / Y-AXIS ---------------- */
function lineYRange(seriesList, xStart, xEnd) {
  let lo = Infinity, hi = -Infinity;
  for (const s of seriesList) for (let i = 0; i < s.x.length; i++) {
    const d = new Date(s.x[i]); if (d < xStart || d > xEnd) continue;
    const v = s.y[i]; if (v == null || isNaN(v)) continue;
    if (v < lo) lo = v; if (v > hi) hi = v;
  }
  if (!isFinite(lo)) return null;
  const pad = (hi - lo) * 0.05 || 0.1;
  return [lo - pad, hi + pad];
}
function decompYRange(xStart, xEnd) {
  let lo = 0, hi = 0;
  for (const row of COMP_SOLID.concat(COMP_NOW)) {
    const d = new Date(row.date); if (d < xStart || d > xEnd) continue;
    let p = 0, n = 0;
    for (const c of COMP) { const v = row[c.key]; if (v > 0) p += v; else n += v; }
    if (p > hi) hi = p; if (n < lo) lo = n;
  }
  const pad = (hi - lo) * 0.05 || 0.1;
  return [lo - pad, hi + pad];
}
function applyRange(years) {
  RANGE = years;
  const end = new Date(END_STR);
  let start, startStr;
  if (years === 'all') { start = new Date('1990-01-01'); }
  else { start = new Date(end); start.setUTCFullYear(start.getUTCFullYear() - years); startStr = start.toISOString().slice(0, 10); }

  CHART_IDS.forEach(id => {
    const yrange = (id === FCI_CHART && decompMode) ? decompYRange(start, end) : lineYRange(CHART_SERIES[id], start, end);
    const update = { 'yaxis.autorange': false };
    if (years === 'all') update['xaxis.autorange'] = true; else update['xaxis.range'] = [startStr, END_STR];
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


/* ---------------- VIEW TOGGLE (line / decomposition) ---------------- */
function setView(mode) {
  const want = (mode === 'decomp');
  if (want === decompMode) return;
  decompMode = want;
  document.body.classList.toggle('decomp-on', decompMode);   // CSS swaps the caption note + its button
  drawChart1();
  applyRange(RANGE);
}
function initViewToggle() {
  document.querySelectorAll('.inline-toggle').forEach(btn => {
    btn.addEventListener('click', () => setView(btn.dataset.view));
  });
}


/* ---------------- MAIN ENTRY ---------------- */
Promise.all([
  fetch(DATA_URL).then(r => r.text()),
  fetch(NOWCAST_URL).then(r => r.text()),
  fetch(COMP_URL).then(r => r.text()),
  fetch(META_URL).then(r => r.json()),
])
.then(([csvText, ncText, compText, meta]) => {
  document.getElementById('last-updated-date').textContent   = formatQuarter(meta.sample_end);
  document.getElementById('last-refreshed-date').textContent = formatMonthYear(meta.last_updated);
  const nowcastEl = document.getElementById('nowcast-through-date');
  if (nowcastEl && meta.nowcast_through) nowcastEl.textContent = formatLongDate(meta.nowcast_through);

  Q  = parseCSV(csvText);
  NC = parseCSV(ncText);
  const comp = parseCSV(compText);
  COMP_SOLID = comp.filter(d => d.kind !== 'nowcast');
  COMP_NOW   = comp.filter(d => d.kind === 'nowcast');

  drawCharts();
  initRangeButtons();
  initViewToggle();
})
.catch(err => {
  console.error('Failed to load data:', err);
  document.getElementById('last-updated-date').textContent   = 'unavailable';
  document.getElementById('last-refreshed-date').textContent = 'unavailable';
});
