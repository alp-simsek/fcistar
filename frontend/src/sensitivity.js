/* =============================================================
   FCI* sensitivity page — sensitivity.js
   Loads fcistar_sensitivity.json (3x3 grid of FCI* over the SPF
   25/50/75 percentiles for GDP growth x core PCE) and fci_nowcast.csv
   (latest FCI nowcast), then renders the current-estimate cards and
   the percentile table (rows = inflation, cols = GDP growth).
   ============================================================= */

const SENS_URL    = '../../backend/data/output/fcistar_sensitivity.json';
const NOWCAST_URL = '../../backend/data/output/fci_nowcast.csv';

const PCTS = ['p25', 'p50', 'p75'];
const PLBL = { p25: '25th', p50: '50th', p75: '75th' };

function fmt(v) { return (v < 0 ? '−' : '') + Math.abs(v).toFixed(2); }
function pct(v) { return v.toFixed(1) + '%'; }
function spaceQuarter(l) { return l ? l.replace('Q', ' Q') : l; }

// Latest FCI nowcast value from fci_nowcast.csv (last row with kind = nowcast).
function lastNowcastFci(csv) {
  const lines = csv.trim().split('\n');
  const hdr = lines[0].split(',');
  const fi = hdr.indexOf('fci'), ki = hdr.indexOf('kind');
  let val = null;
  for (let i = 1; i < lines.length; i++) {
    const v = lines[i].split(',');
    if (v[ki] && v[ki].trim() === 'nowcast') val = parseFloat(v[fi]);
  }
  return val;
}

// Pale heatmap background: low FCI* (looser) -> blue, high (tighter) -> red.
function cellBg(v, lo, hi) {
  const t = hi > lo ? (v - lo) / (hi - lo) : 0.5;
  const mix = (a, b) => Math.round(a + (b - a) * t);
  return `rgb(${mix(0xd6, 0xf5)},${mix(0xe4, 0xd9)},${mix(0xf0, 0xd4)})`;
}

function buildCards(target, fci, fcistar) {
  document.getElementById('sens-quarter').textContent = spaceQuarter(target);
  const gap = fci - fcistar;
  document.getElementById('sens-cards').innerHTML = `
    <div class="card"><div class="card-label">FCI</div>
      <div class="card-value">${fmt(fci)}</div><div class="card-sub">nowcast → quarter-end</div></div>
    <div class="card"><div class="card-label">FCI*</div>
      <div class="card-value">${fmt(fcistar)}</div><div class="card-sub">forecast (median)</div></div>
    <div class="card"><div class="card-label">FCI gap</div>
      <div class="card-value ${gap < 0 ? 'gap-loose' : 'gap-tight'}">${fmt(gap)}</div>
      <div class="card-sub">${gap < 0 ? 'looser than neutral' : 'tighter than neutral'}</div></div>`;
}

function buildTable(sens) {
  const g = sens.fcistar;
  let lo = Infinity, hi = -Infinity;
  for (const ip of PCTS) for (const gp of PCTS) { const v = g[ip][gp]; if (v < lo) lo = v; if (v > hi) hi = v; }

  let h = '<table class="sens-table"><thead>';
  h += '<tr><th class="corner" rowspan="2" colspan="2"></th>'
     + '<th class="grouphdr" colspan="3">GDP growth (SPF percentile)</th></tr>';
  h += '<tr>' + PCTS.map(gp => `<th>${PLBL[gp]}<span class="pctval">${pct(sens.gdp[gp])}</span></th>`).join('') + '</tr>';
  h += '</thead><tbody>';
  PCTS.forEach((ip, ri) => {
    h += '<tr>';
    if (ri === 0) h += '<th class="grouphdr rot" rowspan="3">Core PCE inflation</th>';
    h += `<th class="rowhdr">${PLBL[ip]}<span class="pctval">${pct(sens.corepce[ip])}</span></th>`;
    PCTS.forEach(gp => {
      const v = g[ip][gp];
      const center = (ip === 'p50' && gp === 'p50') ? ' center' : '';
      h += `<td class="cell${center}" style="background:${cellBg(v, lo, hi)}">${fmt(v)}</td>`;
    });
    h += '</tr>';
  });
  h += '</tbody></table>';
  document.getElementById('sens-table').innerHTML = h;
}

Promise.all([
  fetch(SENS_URL).then(r => r.json()),
  fetch(NOWCAST_URL).then(r => r.text()),
])
.then(([sens, ncText]) => {
  buildCards(sens.target_quarter, lastNowcastFci(ncText), sens.fcistar.p50.p50);
  buildTable(sens);
})
.catch(err => {
  console.error('Failed to load sensitivity data:', err);
  document.getElementById('sens-table').textContent = 'Sensitivity data unavailable.';
});
