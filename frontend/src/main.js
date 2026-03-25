/* =============================================================
   FCI* Website — main.js
   Loads fcistar.csv and metadata.json, then draws three
   interactive Plotly charts.
   ============================================================= */

const DATA_URL = '../../backend/data/output/fcistar.csv';
const META_URL = '../../backend/data/output/metadata.json';

// Color palette
const COLORS = {
  fci:     '#7a9cc4',  // steel blue — input series
  fcistar: '#1a3a6b',  // navy      — our main contribution
  gap:     '#c0392b',  // red       — FCI gap
  ygap:    '#1a3a6b',  // navy      — output gap
};


/* -------------------------------------------------------------
   CSV PARSER
   Our CSV is simple (no quoted fields, no commas in values),
   so we parse it manually rather than importing a library.
   ------------------------------------------------------------- */
function parseCSV(text) {
  const lines = text.trim().split('\n');
  const headers = lines[0].split(',').map(h => h.trim());
  return lines.slice(1).map(line => {
    const vals = line.split(',');
    const row = {};
    headers.forEach((h, i) => {
      row[h] = h === 'date' ? vals[i].trim() : parseFloat(vals[i]);
    });
    return row;
  });
}


/* -------------------------------------------------------------
   FORMAT QUARTER
   Converts "2024-12-30" → "2024 Q4" for the header display.
   End-of-quarter dates: Mar=Q1, Jun=Q2, Sep=Q3, Dec=Q4.
   ------------------------------------------------------------- */
function formatQuarter(isoString) {
  const month = parseInt(isoString.slice(5, 7), 10);
  const year  = isoString.slice(0, 4);
  const q     = Math.ceil(month / 3);
  return `${year} Q${q}`;
}


/* -------------------------------------------------------------
   ZERO LINE SHAPE
   Plotly "shapes" are annotations drawn on top of charts.
   xref:'paper' means x0=0,x1=1 spans the full chart width
   regardless of the data range.
   ------------------------------------------------------------- */
function zeroLine() {
  return {
    type: 'line',
    xref: 'paper', x0: 0, x1: 1,
    yref: 'y',     y0: 0, y1: 0,
    line: { color: '#aaa', width: 1, dash: 'dash' }
  };
}


/* -------------------------------------------------------------
   BASE LAYOUT
   Shared Plotly layout properties applied to every chart.
   The spread operator (...baseLayout()) in each chart call
   merges these defaults with any chart-specific overrides.
   ------------------------------------------------------------- */
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
    },
    yaxis: {
      title: { text: yAxisTitle, standoff: 8 },
      zeroline: false,
      gridcolor: '#eee',
      linecolor: '#ccc',
    },
    legend: {
      orientation: 'h',
      x: 0, y: 1.08,
      font: { size: 12 },
      bgcolor: 'rgba(0,0,0,0)',
    },
    hovermode: 'x unified',
  };
}


/* -------------------------------------------------------------
   DRAW CHARTS
   Called once the CSV is parsed. Builds three Plotly charts.
   ------------------------------------------------------------- */
function drawCharts(data) {
  const dates   = data.map(d => d.date);
  const fci     = data.map(d => d.fci);
  const fcistar = data.map(d => d.fcistar);
  const fciGap  = data.map(d => d.fci_gap);
  const yGap    = data.map(d => d.y_gap);

  // responsive:true  → chart resizes when the window is resized
  // displayModeBar:'hover' → Plotly toolbar (zoom, download) appears on hover
  const config = { responsive: true, displayModeBar: 'hover' };

  // --- Chart 1: FCI and FCI* ---
  Plotly.newPlot('chart-fci-and-star', [
    {
      x: dates, y: fci,
      name: 'FCI',
      type: 'scatter', mode: 'lines',
      line: { color: COLORS.fci, width: 2 },
      hovertemplate: '%{y:.2f}<extra>FCI</extra>'
    },
    {
      x: dates, y: fcistar,
      name: 'FCI*',
      type: 'scatter', mode: 'lines',
      line: { color: COLORS.fcistar, width: 2.5 },
      hovertemplate: '%{y:.2f}<extra>FCI*</extra>'
    }
  ], {
    ...baseLayout('Pct. pts. of next-year GDP growth'),
    shapes: [zeroLine()]
  }, config);

  // --- Chart 2: FCI gap ---
  Plotly.newPlot('chart-fci-gap', [
    {
      x: dates, y: fciGap,
      name: 'FCI gap',
      type: 'scatter', mode: 'lines',
      line: { color: COLORS.gap, width: 2 },
      hovertemplate: '%{y:.2f}<extra>FCI gap</extra>'
    }
  ], {
    ...baseLayout('Pct. pts. of next-year GDP growth'),
    shapes: [zeroLine()]
  }, config);

  // --- Chart 3: Output gap ---
  Plotly.newPlot('chart-ygap', [
    {
      x: dates, y: yGap,
      name: 'Output gap',
      type: 'scatter', mode: 'lines',
      line: { color: COLORS.ygap, width: 2 },
      hovertemplate: '%{y:.2f}%<extra>Output gap</extra>'
    }
  ], {
    ...baseLayout('Percent'),
    shapes: [zeroLine()]
  }, config);
}


/* -------------------------------------------------------------
   MAIN ENTRY POINT
   Promise.all fires two fetch requests simultaneously and waits
   for both to finish before proceeding. If either fails, the
   .catch() at the end handles the error gracefully.
   ------------------------------------------------------------- */
Promise.all([
  fetch(DATA_URL).then(r => r.text()),
  fetch(META_URL).then(r => r.json())
])
.then(([csvText, meta]) => {
  document.getElementById('last-updated-date').textContent = formatQuarter(meta.sample_end);
  const data = parseCSV(csvText);
  drawCharts(data);
})
.catch(err => {
  console.error('Failed to load data:', err);
  document.getElementById('last-updated-date').textContent = 'unavailable';
});
