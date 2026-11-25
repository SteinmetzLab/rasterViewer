// Multiscale population-raster viewer (WebGL2).
//
// Every displayed cell (block x block screen pixels) shows the mean firing rate over every
// neuron and instant it covers. Exact at every zoom: all spikes live on the GPU, each frame
// draws the visible time range as points, additively counted into a (cell column x neuron
// row) texture (pass 1), then each cell averages its rows with overlap weights, restricted
// to the shank segment under the cell's center (pass 2), then a log color map (pass 3).
// reference.py implements the same semantics on the CPU; ?test=1 compares the two.
// sonify.js plays the view as sound (the Sonify button); ?sonifytest=1 checks it.
"use strict";

const DATA = "data/";
const L = 150, RM = 118, TOP = 52, BOT = 50;  // stage margins (CSS px) around the raster
const OV_Y = 8, OV_H = 30;                      // overview strip (CSS px from stage top)
const IS_MAC = /Mac/.test(navigator.platform || navigator.userAgent);

const DEPTH_H = 15;   // um: height of a unit's band when plotted by depth (NP2.0 site pitch)
// y of the deepest and shallowest recorded site centers on every shank: NP2.0 bank 0, 192 rows
// at 15 um. Plotted by depth, every shank spans this range (widened if a unit lies outside it),
// so all shanks are drawn identically and empty stretches stay visible.
const DEPTH_SITES_UM = [0, 2865];

let meta, ticks, rows, N, FS, TS, TE, rastermap = null, orderings = null;
let gl, cvs, ov, ctx, prog = {}, res = {};
// n0, n1: the visible span of the vertical coordinate Y of the current y-axis layout (rows for
// "rows" and "rastermap", micrometers for "depth"; see buildLayouts)
let state = {
  t0: 0, t1: 1, n0: 0, n1: 1, ymode: "rows", mode: "abs", eps: 0.1, block: 2,
  spikes: "all", burst: { S: 100, M: 4, K: 2 },
  abs: { vmin: 0.3, vmax: 50, cmap: "magma" },
  norm: { vmin: 0.1, vmax: 10, cmap: "viridis" },
};
let layout = {}, lastDraw = { i0: 0, i1: 0, ms: 0 }, dirty = true, lastMouse = null;
let LAY = {};           // y-axis layouts by name
let SETS = {};          // spike sets on the GPU: all, and burst / onset subsets
let burst = { flags: null, params: {} };   // current burst classification
let pinned = null;      // {x, y, key}: tooltip pinned by a click, until the view changes
let hitboxes = [];      // clickable probe / shank labels on the neuron axis
let altHeld = false;    // Alt / Option held: the tooltip follows the pointer

// ---------------------------------------------------------------- loading
async function fetchBinary(url, label) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`${url}: HTTP ${resp.status}`);
  const total = +resp.headers.get("Content-Length");
  const buf = new Uint8Array(total);
  const reader = resp.body.getReader();
  let got = 0, lastShow = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    buf.set(value, got); got += value.length;
    if (got - lastShow > 8e6) {
      lastShow = got;
      document.getElementById("loadMsg").textContent =
        `Loading ${label}: ${(got / 1e6).toFixed(0)} of ${(total / 1e6).toFixed(0)} MB`;
    }
  }
  return buf.buffer;
}

async function load() {
  const raw = RawUnderlay.load();        // optional: raw voltage (raw_underlay.js), found meanwhile
  meta = await (await fetch(DATA + "meta.json")).json();
  N = meta.n_rows; FS = meta.fs; TS = meta.t_start;
  const mr = await fetch(DATA + "manifest.json");
  if (mr.ok) await loadPacked(await mr.json());       // compressed (pack_data.py)
  else {                                               // raw files (prep_data.py)
    ticks = new Uint32Array(await fetchBinary(DATA + "ticks.u32", "spike times"));
    rows = new Uint16Array(await fetchBinary(DATA + "rows.u16", "spike neurons"));
  }
  TE = Math.max(meta.t_end, (ticks[ticks.length - 1] + 1) / FS);  // view end is exclusive
  if (ticks.length !== meta.n_spikes || rows.length !== meta.n_spikes)
    throw new Error("spike file lengths do not match meta.json");
  const rr = await fetch(DATA + "rastermap.json");       // optional (prep_rastermap.py)
  if (rr.ok) rastermap = await rr.json();
  const ro = await fetch(DATA + "orderings.json");       // optional: further orderings of the neurons
  if (ro.ok) orderings = await ro.json();
  await raw;
}

// Spikes packed by pack_data.py: gzip streams split into parts. Decoded while downloading,
// then checked against the checksums pack_data.py computed from the original arrays.
async function loadPacked(man) {
  if (man.format !== "rasterViewer-spikes-1") throw new Error("unknown data format " + man.format);
  if (typeof DecompressionStream === "undefined")
    throw new Error("this browser cannot decompress gzip streams (needs Chrome 80+, Safari 16.4+ or Firefox 113+)");
  const n = man.n_spikes;
  ticks = new Uint32Array(n); rows = new Uint16Array(n);
  const total = Object.values(man.streams).reduce((a, s) => a + s.parts.reduce((b, p) => b + p.bytes, 0), 0);
  let got = 0, shown = 0;
  const progress = () => {
    if (got - shown < 4e6 && got < total) return;
    shown = got;
    document.getElementById("loadMsg").textContent =
      `Loading spikes: ${(got / 1e6).toFixed(0)} of ${(total / 1e6).toFixed(0)} MB`;
  };
  async function readStream(stream, onBytes) {
    const ds = new DecompressionStream("gzip"), writer = ds.writable.getWriter();
    const pump = (async () => {
      try {
        for (const p of stream.parts) {
          const resp = await fetch(DATA + p.file);
          if (!resp.ok) throw new Error(`${p.file}: HTTP ${resp.status}`);
          const rd = resp.body.getReader();
          for (;;) {
            const { done, value } = await rd.read();
            if (done) break;
            got += value.length; progress();
            await writer.write(value);
          }
        }
        await writer.close();
      } catch (err) { writer.abort(err).catch(() => {}); throw err; }
    })();
    const rd = ds.readable.getReader();
    let nb = 0;
    const consume = (async () => {
      for (;;) { const { done, value } = await rd.read(); if (done) break; onBytes(value, nb); nb += value.length; }
    })();
    await Promise.all([pump, consume]);
    if (nb !== stream.raw_bytes) throw new Error(`decoded ${nb} bytes, expected ${stream.raw_bytes}`);
  }
  // spike times: first tick, then uint8 deltas (255 = escape, true value in tick_escapes)
  let cur = man.tick0, i = 0, e = 0;
  const escapes = man.tick_escapes;
  ticks[0] = cur;
  await readStream(man.streams.ticks, (b) => {
    for (let q = 0; q < b.length; q++) {
      const d = b[q] === 255 ? escapes[e++] : b[q];
      cur += d; ticks[++i] = cur;
    }
  });
  // rows: all high bytes, then all low bytes
  await readStream(man.streams.rows, (b, off) => {
    for (let q = 0; q < b.length; q++) {
      const k = off + q;
      if (k < n) rows[k] = b[q] << 8; else rows[k - n] |= b[q];
    }
  });
  let st = 0, sr = 0, sx = 0;
  for (let k = 0; k < n; k++) { st += ticks[k]; sr += rows[k]; sx += rows[k] * (k & 1023); }
  const c = man.checksums;
  if (st !== c.sum_ticks || sr !== c.sum_rows || sx !== c.sum_rows_x_index_mod_1024 ||
      ticks[0] !== c.first_tick || ticks[n - 1] !== c.last_tick)
    throw new Error("spike data failed its checksum after decoding");
}

// ---------------------------------------------------------------- y-axis layouts
// Each layout places every neuron on a vertical coordinate Y: position p shows row order[p]
// and occupies [lo[p], lo[p] + unit). Y is cut into segments; a cell only averages neurons of
// the segment under its center. Same definitions as reference.py's Data.layout.
//   rows:      row i at [i, i+1); segments = (probe, shank)
//   depth:     micrometers, to scale; every (probe, shank) segment spans the same full recorded
//              range of y (all sites, see DEPTH_SITES_UM) plus h/2 at each end; a unit is a band
//              of height h centered on its y, shallow at top; units may overlap
//   rastermap: row isort[k] at [k, k+1); one segment (cells aggregate across probes)
//   Each entry of data/orderings.json is one more layout of the rastermap kind (name
//   "rastermap", its own order), stored under its key and listed in the Neurons menu.
function buildLayouts() {
  const ident = Int32Array.from({ length: N }, (_, i) => i);
  const idx = Float64Array.from(ident);
  const segRows = meta.segments.map((s) => ({ y0: s.start, y1: s.end, probe: s.probe,
                                               shank: s.shank, r0: s.start, r1: s.end }));
  LAY.rows = { name: "rows", Y: N, unit: 1, order: ident, lo: idx, segs: segRows, frac: true };
  if (rastermap) {
    LAY.rastermap = { name: "rastermap", Y: N, unit: 1, order: Int32Array.from(rastermap.isort),
                      lo: idx, segs: [{ y0: 0, y1: N, probe: null, shank: null, r0: 0, r1: N }],
                      frac: true };
  }
  for (const o of (orderings ? orderings.orderings : [])) {
    if (!/^[a-z][a-z0-9_]*$/.test(o.key) || Object.hasOwn(LAY, o.key) || o.key === "depth")
      throw new Error("orderings.json: bad or duplicate key " + o.key);
    LAY[o.key] = { name: "rastermap", key: o.key, label: o.label, axis: o.axis, Y: N, unit: 1,
                   order: Int32Array.from(o.isort), lo: idx,
                   segs: [{ y0: 0, y1: N, probe: null, shank: null, r0: 0, r1: N }], frac: true };
  }
  // depth, in integer units of 0.1 um (y_um is stored to 0.1 um), so positions are exact
  const h10 = Math.round(DEPTH_H * 10), y10 = meta.rows.y_um.map((y) => Math.round(y * 10));
  const lo = new Float64Array(N), local10 = new Float64Array(N), segs = [];
  const mn = Math.min(Math.round(DEPTH_SITES_UM[0] * 10), ...y10);
  const mx = Math.max(Math.round(DEPTH_SITES_UM[1] * 10), ...y10);
  const top10 = mx + h10 / 2, span10 = mx - mn + h10;      // the same for every shank
  let base = 0;
  for (const s of meta.segments) {
    for (let i = s.start; i < s.end; i++) {
      local10[i] = top10 - y10[i] - h10 / 2;            // 0 for a unit at the top site
      lo[i] = (base + local10[i]) / 10;
    }
    segs.push({ y0: base / 10, y1: (base + span10) / 10, probe: s.probe, shank: s.shank,
                r0: s.start, r1: s.end, top10 });
    base += span10;
  }
  LAY.depth = { name: "depth", Y: base / 10, unit: DEPTH_H, order: ident, lo, segs, local10,
                frac: false };
  // assumptions the cell computations rely on
  const check = (ok, msg) => { if (!ok) throw new Error("layout check failed: " + msg); };
  meta.rows.y_um.forEach((y, i) => check(Math.abs(y * 10 - y10[i]) < 1e-6, `y_um[${i}] not a multiple of 0.1`));
  for (const g of segs) for (let i = g.r0 + 1; i < g.r1; i++) check(lo[i] >= lo[i - 1], `depth bands unsorted at row ${i}`);
  meta.segments.forEach((g, i) => check(g.end > g.start && (i === 0 || g.start === meta.segments[i - 1].end), "segments not contiguous"));
  for (const [key, L] of Object.entries(LAY)) if (L.name === "rastermap") {
    const seen = new Uint8Array(N);
    for (const r of L.order) { check(r >= 0 && r < N && !seen[r], key + " order is not a permutation"); seen[r] = 1; }
    check(L.order.length === N, key + " order has the wrong length");
  }
  const sel = document.getElementById("ymode");      // list the further orderings in the menu
  for (const L of Object.values(LAY)) if (L.key) {
    const op = document.createElement("option");
    op.value = L.key; op.textContent = L.label; sel.appendChild(op);
  }
}
const lay = () => LAY[state.ymode];
function segAt(L, y) {               // index of the segment containing y (clamped)
  let a = 0, b = L.segs.length - 1;
  while (a < b) { const m = (a + b + 1) >> 1; if (L.segs[m].y0 <= y) a = m; else b = m - 1; }
  return a;
}

// ---------------------------------------------------------------- WebGL setup
const VS_POINTS = `#version 300 es
layout(location=0) in uint a_tick;
layout(location=1) in uint a_row;
uniform uint u_base;   // ceil(T0), ticks
uniform float u_frac;  // u_base - T0, in [0, 1)
uniform float u_dT;    // ticks per cell column
uniform int u_W, u_rlo, u_R;
// u_bounds[j] = ceil(T0 + j * dT), computed in float64 on the CPU: the first integer tick of
// column j. The float estimate below can be off by one column near an edge (float32 has
// ~1e-7 relative precision); the integer comparisons then make the assignment exact.
uniform highp usampler2D u_bounds;
void main() {
  gl_PointSize = 1.0;
  int r = int(a_row) - u_rlo;
  if (r < 0 || r >= u_R) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  float f = (float(a_tick - u_base) + u_frac) / u_dT;
  int col = clamp(int(floor(f)), 0, u_W - 1);
  for (int it = 0; it < 4; it++) {
    if (col > 0 && a_tick < texelFetch(u_bounds, ivec2(col, 0), 0).r) col--;
    else if (col < u_W - 1 && a_tick >= texelFetch(u_bounds, ivec2(col + 1, 0), 0).r) col++;
    else break;
  }
  gl_Position = vec4((float(col) + 0.5) / float(u_W) * 2.0 - 1.0,
                     (float(r) + 0.5) / float(u_R) * 2.0 - 1.0, 0.0, 1.0);
}`;
const FS_ONE = `#version 300 es
precision highp float;
out vec4 o;
void main() { o = vec4(1.0, 0.0, 0.0, 1.0); }`;
// Pass 1, alternative ("search" counting): one fragment per (column j, row r). Each row's spike
// ticks are stored sorted and contiguous (u_rowticks, row r at [off[r], off[r+1])), so the
// row's count in column j is lowerBound(bounds[j+1]) - lowerBound(bounds[j]), found by binary
// search with integer comparisons. Gives the same exact counts as drawing points, but its cost
// scales with the grid (W x R) instead of the number of spikes in the window. That matters on
// tile-based GPUs (Apple silicon), where rasterizing tens of millions of points is slow.
// TW (texture width, a power of 2) is substituted before compiling.
const FS_SEARCH = `#version 300 es
precision highp float; precision highp int; precision highp usampler2D;
uniform highp usampler2D u_bounds;    // R32UI, (W+1) x 1: first tick of each column
uniform highp usampler2D u_rowticks;  // R32UI, TW x TH: spike ticks sorted by (row, tick)
uniform highp usampler2D u_rowoff;    // R32UI, 1024 wide: first index of each row, N+1 entries
uniform int u_rlo;
out vec4 o;
uint tickAt(uint i) { return texelFetch(u_rowticks, ivec2(int(i & (TWu - 1u)), int(i >> TWLOG2u)), 0).r; }
uint offAt(int i) { return texelFetch(u_rowoff, ivec2(i & 1023, i >> 10), 0).r; }
uint lowerBound(uint lo, uint hi, uint x) {   // first index in [lo, hi) with tick >= x
  for (int it = 0; it < 32; it++) {           // 32 halvings cover any uint range
    if (lo >= hi) break;
    uint m = lo + ((hi - lo) >> 1);
    if (tickAt(m) < x) lo = m + 1u; else hi = m;
  }
  return lo;
}
void main() {
  int j = int(gl_FragCoord.x), row = u_rlo + int(gl_FragCoord.y);
  uint s = offAt(row), e = offAt(row + 1);
  uint a = texelFetch(u_bounds, ivec2(j, 0), 0).r, b = texelFetch(u_bounds, ivec2(j + 1, 0), 0).r;
  uint ia = lowerBound(s, e, a);
  uint ib = lowerBound(ia, e, b);
  o = vec4(float(ib - ia), 0.0, 0.0, 1.0);    // exact: a row has far fewer than 2^24 spikes
}`;
const VS_FULL = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;
// Pass 2: one fragment per cell. Which neurons each cell row averages over, and with what
// weights, depends on the y-axis layout. The CPU (cellSpans, float64) supplies per cell row:
//  - rows / rastermap (u_depth == 0): positions ia..ceil(b)-1 and the fractional overlaps fa, fb
//    at the two ends; ov(p) = (p == ib ? fb : 1) - (p == ia ? fa : 0). Position p shows row
//    u_order[p].
//  - depth (u_depth == 1): candidate rows ia..iend-1 and the cell's span [a, b) in um relative
//    to row ia's band start. Each row's band start within its segment is an integer number of
//    0.1 um (rowinfo.w), so lo(i) = (w(i) - w(ia)) / 10 is exact, and the overlap
//    ov(i) = clamp(min(b, lo + h) - max(a, lo)) / h is accurate even when cells are tiny.
// cr.x < 0 marks a cell with no neurons (value -2).
// Two spike counts per texel come in (u_counts = set A, u_counts2 = set B, e.g. all spikes and
// burst spikes). Two values go out: r from count u_c1 . (A, B), b from u_c2 . (A, B), so that
// e.g. tonic = all - burst is formed from exact integer counts.
const FS_AGG = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D u_counts, u_counts2;  // R32F, W x R: counts of rows u_rlo .. u_rlo + R - 1
uniform sampler2D u_rowinfo;   // RGBA32F, 1024 wide: seg_start, seg_end, r0, depth band start
uniform sampler2D u_cellrows;  // RGBA32F, H x 1 (see above)
uniform sampler2D u_order;     // R32F, 1024 wide: row shown at each position
uniform highp usampler2D u_bounds;  // R32UI, (W+1) x 1: first tick of each column
uniform float u_fs, u_eps, u_h;
uniform vec2 u_c1, u_c2;       // weights of (A, B) for the two outputs
uniform int u_H, u_rlo, u_norm, u_depth;
out vec4 o;
vec2 cnt(int j, int r) {       // (output 1, output 2) counts of texel (j, r)
  vec2 ab = vec2(texelFetch(u_counts, ivec2(j, r), 0).r, texelFetch(u_counts2, ivec2(j, r), 0).r);
  return vec2(dot(u_c1, ab), dot(u_c2, ab));
}
vec4 rowinfo(int i) { return texelFetch(u_rowinfo, ivec2(i % 1024, i / 1024), 0); }
int orderAt(int p) { return int(texelFetch(u_order, ivec2(p % 1024, p / 1024), 0).r); }
float weight(int row) { return (u_norm == 1) ? 1.0 / (rowinfo(row).z + u_eps) : 1.0; }
void main() {
  int j = int(gl_FragCoord.x);
  int k = u_H - 1 - int(gl_FragCoord.y);          // k = 0 is the top cell row
  // the column's duration: its whole number of samples (a tick stands for [s-0.5, s+0.5))
  uint nt = texelFetch(u_bounds, ivec2(j + 1, 0), 0).r - texelFetch(u_bounds, ivec2(j, 0), 0).r;
  if (nt == 0u) { o = vec4(-1.0, 0.0, -1.0, 1.0); return; }
  float dtSec = float(nt) / u_fs;
  vec4 cr = texelFetch(u_cellrows, ivec2(k, 0), 0);
  if (cr.x < 0.0) { o = vec4(-2.0, 0.0, -2.0, 1.0); return; }
  vec2 num = vec2(0.0);
  float den = 0.0;
  if (u_depth == 0) {
    int ia = int(cr.x), ib = int(cr.z);
    float fa = cr.y, fb = cr.w;
    int iend = fb > 0.0 ? ib + 1 : ib;
    for (int p = ia; p < iend; p++) {
      float ov = (p == ib ? fb : 1.0) - (p == ia ? fa : 0.0);
      if (ov <= 0.0) continue;
      int row = orderAt(p);
      num += ov * weight(row) * cnt(j, row - u_rlo);
      den += ov;
    }
  } else {
    int ia = int(cr.x), iend = int(cr.z);
    float a = cr.y, b = cr.w, w0 = rowinfo(ia).w;
    vec2 numEq = vec2(0.0);
    for (int i = ia; i < iend; i++) {
      float lo = (rowinfo(i).w - w0) * 0.1;
      vec2 wc = weight(i) * cnt(j, i - u_rlo);
      numEq += wc;
      float ov = (min(b, lo + u_h) - max(a, lo)) / u_h;
      if (ov <= 0.0) continue;
      num += ov * wc;
      den += ov;
    }
    // Every candidate overlaps the cell (decided in float64 on the CPU). If all the overlaps
    // are slivers that round to 0 in float32, weight them equally (their limit) rather than
    // calling the cell empty.
    if (den <= 0.0) { num = numEq / float(iend - ia); den = 1.0; }
  }
  vec2 v = den > 0.0 ? num / (den * dtSec) : vec2(-2.0);
  o = vec4(v.x, den, v.y, 1.0);
}`;
const FS_CMAP = `#version 300 es
precision highp float; precision highp sampler2D;
uniform sampler2D u_val, u_lut, u_lut2;
uniform float u_lmin, u_lmax;
uniform int u_overlay;         // 1: r through u_lut plus b through u_lut2, added
uniform vec4 u_under;          // fragment-coordinate box (x0, y0, x1, y1) with raw voltage underneath
uniform float u_ualpha;        // opacity of the cells with spikes in that box; < 0: no such box
out vec4 o;
float tOf(float v) { return v > 0.0 ? clamp((log(v) - u_lmin) / (u_lmax - u_lmin), 0.0, 1.0) : 0.0; }
vec4 lut(sampler2D s, float t) { return texelFetch(s, ivec2(int(t * 255.0 + 0.5), 0), 0); }
void main() {
  vec4 val = texelFetch(u_val, ivec2(gl_FragCoord.xy), 0);
  float v = val.r;
  // No neuron at this depth (-2): drawn as a rate of zero, like any cell without spikes.
  // The value stays -2, so the readout and ?test=1 still tell the two cases apart.
  if (v < -1.5) v = 0.0;
  if (v < 0.0) { o = vec4(0.25, 0.25, 0.25, 1.0); return; }    // column holds no samples
  float alpha = 1.0;
  if (u_ualpha >= 0.0 && all(greaterThanEqual(gl_FragCoord.xy, u_under.xy)) &&
      all(lessThan(gl_FragCoord.xy, u_under.zw))) {
    bool spikes = v > 0.0 || (u_overlay == 1 && val.b > 0.0);
    if (!spikes) { o = vec4(0.0); return; }                    // the voltage shows through
    alpha = u_ualpha;
  }
  if (u_overlay == 1) {
    o = vec4(min(lut(u_lut, tOf(v)).rgb + lut(u_lut2, tOf(val.b)).rgb, 1.0), alpha);
    return;
  }
  o = vec4(lut(u_lut, tOf(v)).rgb, alpha);
}`;

function compile(vs, fs) {
  const p = gl.createProgram();
  for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    gl.attachShader(p, s);
  }
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  const u = {};
  const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < n; i++) {
    const name = gl.getActiveUniform(p, i).name;
    u[name] = gl.getUniformLocation(p, name);
  }
  return { p, u };
}

function initGL() {
  cvs = document.getElementById("gl");
  // alpha: where raw voltage is drawn on the canvas underneath (raw_underlay.js), cells
  // without spikes are transparent; everywhere else the output is opaque
  gl = cvs.getContext("webgl2", { antialias: false, preserveDrawingBuffer: true,
                                  premultipliedAlpha: false, alpha: true });
  if (!gl) throw new Error("WebGL2 is not available in this browser");
  if (!gl.getExtension("EXT_color_buffer_float"))
    throw new Error("EXT_color_buffer_float is not available");
  // Counting needs blending into a 32-bit float target. A 16-bit fallback would silently
  // cap counts at 2048, so refuse instead.
  res.floatBlend = !!gl.getExtension("EXT_float_blend");
  if (!res.floatBlend)
    throw new Error("this browser/GPU lacks EXT_float_blend, which exact counting needs");
  prog.points = compile(VS_POINTS, FS_ONE);
  prog.agg = compile(VS_FULL, FS_AGG);
  prog.cmap = compile(VS_FULL, FS_CMAP);

  SETS.all = { name: "all", ticks, rows, search: null, ...makeVao(ticks, rows) };
  res.emptyVao = gl.createVertexArray();

  // per-row info texture
  const H = Math.ceil(N / 1024), info = new Float32Array(1024 * H * 4);
  for (let i = 0; i < N; i++) {
    info[4 * i] = meta.rows.seg_start[i]; info[4 * i + 1] = meta.rows.seg_end[i];
    info[4 * i + 2] = meta.rows.r0[i]; info[4 * i + 3] = LAY.depth.local10[i];
  }
  res.rowinfo = makeTex(1024, H, gl.RGBA32F, gl.RGBA, gl.FLOAT, info);
  for (const L of Object.values(LAY)) {        // position -> row, exact in float32 (N < 2^24)
    const o = new Float32Array(1024 * H); o.set(L.order);
    L.orderTex = makeTex(1024, H, gl.R32F, gl.RED, gl.FLOAT, o);
  }
  res.lut = makeTex(256, 1, gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, null);
  const ramp = (r, g, b) => {   // black to (r, g, b), linear
    const px = new Uint8Array(256 * 4);
    for (let i = 0; i < 256; i++) px.set([r * i, g * i, b * i, 255], 4 * i);
    return makeTex(256, 1, gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, px);
  };
  res.lutBurst = ramp(1, 0, 1); res.lutTonic = ramp(0, 1, 0);   // magenta + green = white
  res.counts = { tex: null, fbo: gl.createFramebuffer(), w: 0, h: 0 };
  res.counts2 = { tex: null, fbo: gl.createFramebuffer(), w: 0, h: 0 };
  res.vals = { tex: null, fbo: gl.createFramebuffer(), w: 0, h: 0 };
  // Search counting is an optimization: if anything about it fails, use points only.
  res.search = null; res.searchNote = "off (?count=points)";
  if (COUNT_MODE !== "points") {
    try { initSearch(); verifySearch(); }
    catch (err) {
      console.warn("search counting disabled:", err);
      if (SETS.all.search) {
        gl.deleteTexture(SETS.all.search.rowticks); gl.deleteTexture(SETS.all.search.rowoff);
        SETS.all.search = null;
      }
      res.search = null; res.searchNote = "off: " + err.message;
    }
  }
}

// Data for "search" counting: spike ticks regrouped by row (a stable counting sort, so each
// row stays in time order), in a 2D integer texture, plus each row's start index. Needs a
// texture big enough for all spikes; without one, only point counting is used.
function initSearch() {
  const maxT = gl.getParameter(gl.MAX_TEXTURE_SIZE);
  const TW = 1 << Math.floor(Math.log2(Math.min(16384, maxT)));
  const src = FS_SEARCH.replaceAll("TWLOG2u", `${Math.log2(TW)}u`).replaceAll("TWu", `${TW}u`);
  if (/\bTW(u|LOG2u)\b/.test(src)) throw new Error("search shader template not fully substituted");
  res.search = { prog: compile(VS_FULL, src), TW, maxT, ratio: SEARCH_RATIO_DEFAULT };
  SETS.all.search = makeSearchTex(ticks, rows);
}

// Row-grouped copy of a set's spikes for search counting (texture of ticks, and row offsets).
function makeSearchTex(tk, rw) {
  const { TW, maxT } = res.search;
  const n = tk.length, TH = Math.max(1, Math.ceil(n / TW));
  if (TH > maxT) throw new Error(`${n} spikes need a ${TW}x${TH} texture, max ${maxT}`);
  while (gl.getError() !== gl.NO_ERROR) {}   // clear stale errors so the check below is ours
  const off = new Uint32Array(1024 * Math.ceil((N + 1) / 1024));   // padded to whole texture rows
  for (let i = 0; i < n; i++) off[rw[i] + 1]++;
  for (let r = 0; r < N; r++) off[r + 1] += off[r];
  const cur = off.slice(0, N), rt = new Uint32Array(TW * TH);
  for (let i = 0; i < n; i++) rt[cur[rw[i]]++] = tk[i];
  const rowticks = makeTex(TW, TH, gl.R32UI, gl.RED_INTEGER, gl.UNSIGNED_INT, rt);
  const rowoff = makeTex(1024, off.length / 1024, gl.R32UI, gl.RED_INTEGER, gl.UNSIGNED_INT, off);
  const e = gl.getError();
  if (e !== gl.NO_ERROR) {        // e.g. out of GPU memory for the 244 MB texture
    gl.deleteTexture(rowticks); gl.deleteTexture(rowoff);
    throw new Error(`texture upload failed (GL error 0x${e.toString(16)})`);
  }
  return { rowticks, rowoff };
}

function makeVao(tk, rw) {
  const vao = gl.createVertexArray(), bufs = [gl.createBuffer(), gl.createBuffer()];
  gl.bindVertexArray(vao);
  gl.bindBuffer(gl.ARRAY_BUFFER, bufs[0]); gl.bufferData(gl.ARRAY_BUFFER, tk, gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0); gl.vertexAttribIPointer(0, 1, gl.UNSIGNED_INT, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, bufs[1]); gl.bufferData(gl.ARRAY_BUFFER, rw, gl.STATIC_DRAW);
  gl.enableVertexAttribArray(1); gl.vertexAttribIPointer(1, 1, gl.UNSIGNED_SHORT, 0, 0);
  gl.bindVertexArray(null);
  return { vao, bufs };
}

function freeSet(S) {
  if (!S) return;
  gl.deleteVertexArray(S.vao); S.bufs.forEach((b) => gl.deleteBuffer(b));
  if (S.search) { gl.deleteTexture(S.search.rowticks); gl.deleteTexture(S.search.rowoff); }
}

// ---------------------------------------------------------------- bursts
// A burst: a run of >= K consecutive spikes of one unit, every ISI <= M ms, whose first spike
// follows >= S ms of silence from that unit (S = 0: no silence needed). Every spike of the run
// is a burst spike; its first spike is the burst onset. Same definition and arithmetic as
// bursts/detect_bursts.py (the self-test compares the two). One pass in time order, keeping per
// unit its last spike and the current run (its first K - 1 spikes are held until the run
// reaches K spikes). Default: the thalamic LTS convention (Lu, Guido & Sherman 1992).
const BURST_DEFAULT = { S: 100, M: 4, K: 2 };
function detectBursts(S_ms, M_ms, K) {
  K = Math.max(2, Math.round(K));    // the run buffer below assumes K >= 2 (the UI enforces it)
  const n = ticks.length, flags = new Uint8Array(n);       // bit 0 burst spike, bit 1 onset
  const M = M_ms * FS / 1000, S = S_ms * FS / 1000, KB = K - 1;
  const last = new Float64Array(N).fill(-Infinity), len = new Int32Array(N);
  const ok = new Uint8Array(N), buf = new Int32Array(N * KB);
  for (let k = 0; k < n; k++) {
    const r = rows[k], t = ticks[k], pre = t - last[r];
    last[r] = t;
    if (pre <= M) {
      const L = ++len[r];
      if (!ok[r]) continue;
      if (L < K) buf[r * KB + L - 1] = k;
      else if (L === K) {
        for (let q = 0; q < KB; q++) flags[buf[r * KB + q]] |= 1;
        flags[buf[r * KB]] |= 2; flags[k] |= 1;
      } else flags[k] |= 1;
    } else {
      len[r] = 1; ok[r] = pre >= S ? 1 : 0; buf[r * KB] = k;
    }
  }
  return flags;
}

// Run detection and build the GPU sets of burst spikes and of burst onsets.
function buildBurstSets(S_ms, M_ms, K) {
  const t0 = performance.now(), flags = detectBursts(S_ms, M_ms, K);
  const sets = {};
  for (const [name, bit] of [["burst", 1], ["onset", 2]]) {
    let m = 0;
    for (let k = 0; k < flags.length; k++) if (flags[k] & bit) m++;
    const tk = new Uint32Array(m), rw = new Uint16Array(m);
    let q = 0, sumIdx = 0;
    for (let k = 0; k < flags.length; k++) if (flags[k] & bit) { tk[q] = ticks[k]; rw[q++] = rows[k]; sumIdx += k; }
    freeSet(SETS[name]);
    SETS[name] = { name, ticks: tk, rows: rw, n: m, sumIdx, ...makeVao(tk, rw), search: null };
    if (res.search) {
      try { SETS[name].search = makeSearchTex(tk, rw); }
      catch (err) { console.warn(`search counting off for ${name}:`, err); }
    }
  }
  burst.flags = flags; burst.params = { S: S_ms, M: M_ms, K };
  burst.ms = performance.now() - t0;
  console.log(`bursts (${S_ms} ms, ${M_ms} ms, ${K}): ${SETS.onset.n} bursts, ` +
              `${SETS.burst.n} burst spikes, ${burst.ms.toFixed(0)} ms`);
}
function ensureBursts() {
  const p = state.burst, q = burst.params;
  if (!burst.flags || q.S !== p.S || q.M !== p.M || q.K !== p.K) buildBurstSets(p.S, p.M, p.K);
}

// Which spike sets a display mode counts, and how the two outputs combine them:
// value 1 = c1 . (A, B), value 2 = c2 . (A, B).
const SPIKE_MODES = {
  all:     { A: "all",   B: null,    c1: [1, 0],  c2: [0, 0],  label: "All spikes" },
  tonic:   { A: "all",   B: "burst", c1: [1, -1], c2: [0, 0],  label: "Tonic spikes" },
  burst:   { A: "burst", B: null,    c1: [1, 0],  c2: [0, 0],  label: "Burst spikes" },
  onset:   { A: "onset", B: null,    c1: [1, 0],  c2: [0, 0],  label: "Burst onsets" },
  overlay: { A: "all",   B: "burst", c1: [0, 1],  c2: [1, -1], label: "Burst + tonic" },
};

// Startup check of search counting on this GPU and browser, and calibration of when to use it.
// Counts the session's last 300 s (largest tick values) for all rows on a 512-column grid with
// both methods. The counts must be bitwise identical and sum to the number of spikes in the
// window; otherwise search is disabled. (A shorter window keeps the startup cost to ~0.1 s.)
// The second timing of each method (the first includes one-time setup) gives this machine's
// cost per point and per search texel, hence the spikes-per-texel break-even for "auto".
function verifySearch() {
  const t0 = performance.now(), span = Math.min(300, TE - TS);
  const W = 512, v = { T0: (TE - span) * FS, dT: span * FS / W, n0: 0, n1: N };
  const sync = () => {
    gl.bindFramebuffer(gl.FRAMEBUFFER, res.counts.fbo);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, new Float32Array(4));
  };
  const run = (m) => { countSpikes(v, W, m, 0, N); sync(); const t = performance.now();
                       const r = countSpikes(v, W, m, 0, N); sync(); return { r, ms: performance.now() - t }; };
  const p = run("points"), a = readCounts(W, N);
  const s = run("search"), b = readCounts(W, N);
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) throw new Error(`startup check: counts differ at texel ${i} (${a[i]} vs ${b[i]})`);
    sum += a[i];
  }
  const nsp = p.r.i1 - p.r.i0;
  if (sum !== nsp) throw new Error(`startup check: counted ${sum} of ${nsp} spikes`);
  if (p.r.path !== "points" || s.r.path !== "search") throw new Error("startup check: wrong method ran");
  // break-even spikes per texel = (search ms per texel) / (points ms per spike), kept in a sane range
  const ratio = (s.ms / (W * N)) / (p.ms / nsp);
  if (Number.isFinite(ratio) && ratio > 0) res.search.ratio = Math.min(3, Math.max(0.05, ratio));
  res.search.calib = { pointsMs: p.ms, searchMs: s.ms, W, rows: N, spikes: nsp, rawRatio: ratio,
                       checkMs: performance.now() - t0 };
  res.searchNote = `on, verified at startup; auto threshold ${res.search.ratio.toPrecision(2)} spikes/texel`;
  console.log("search counting:", res.searchNote, res.search.calib);
}

function readCounts(W, R, target = res.counts) {   // -> Float32Array R x W of pass-1 counts
  gl.bindFramebuffer(gl.FRAMEBUFFER, target.fbo);
  const px = new Float32Array(W * R * 4);
  gl.readPixels(0, 0, W, R, gl.RGBA, gl.FLOAT, px);
  const out = new Float32Array(W * R);
  for (let i = 0; i < W * R; i++) out[i] = px[4 * i];
  return out;
}

// Which pass-1 method to use: "auto" (default) or forced with ?count=points / ?count=search.
const COUNT_MODE = (() => {
  const c = new URLSearchParams(location.search).get("count");
  return c === "points" || c === "search" ? c : "auto";
})();
// In "auto", search when the window holds more than res.search.ratio spikes per (column x row)
// texel of the counts grid. verifySearch() measures it on each machine at startup; this is the
// fallback. On the M3 Max (Chrome/ANGLE Metal): about 2.6 ns per point and 0.75 ns per search
// texel, so about 0.3. Both methods are fast near the break-even.
const SEARCH_RATIO_DEFAULT = 0.3;

function makeTex(w, h, ifmt, fmt, type, data) {
  const t = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texImage2D(gl.TEXTURE_2D, 0, ifmt, w, h, 0, fmt, type, data);
  for (const [k, v] of [[gl.TEXTURE_MIN_FILTER, gl.NEAREST], [gl.TEXTURE_MAG_FILTER, gl.NEAREST],
                        [gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE]])
    gl.texParameteri(gl.TEXTURE_2D, k, v);
  return t;
}

function ensureTarget(t, w, h, ifmt, fmt, type) {
  if (t.w === w && t.h === h && t.tex) return;
  if (t.tex) gl.deleteTexture(t.tex);
  t.tex = makeTex(w, h, ifmt, fmt, type, null); t.w = w; t.h = h;
  gl.bindFramebuffer(gl.FRAMEBUFFER, t.fbo);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, t.tex, 0);
  const st = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
  if (st !== gl.FRAMEBUFFER_COMPLETE) throw new Error("framebuffer incomplete: " + st);
}

function setLUT(name) {
  const src = meta.colormaps[name], px = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) {
    px[4 * i] = src[3 * i]; px[4 * i + 1] = src[3 * i + 1]; px[4 * i + 2] = src[3 * i + 2];
    px[4 * i + 3] = 255;
  }
  gl.bindTexture(gl.TEXTURE_2D, res.lut);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, 256, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, px);
}

// ---------------------------------------------------------------- core computation
function lowerBound(arr, x) {      // first index with arr[i] >= x
  let lo = 0, hi = arr.length;
  while (lo < hi) { const m = (lo + hi) >>> 1; if (arr[m] < x) lo = m + 1; else hi = m; }
  return lo;
}

// Pass 1. v = {T0, dT (ticks)}, W columns, rows rlo..rhi-1. Leaves per-(column, row) spike
// counts in res.counts (W x R) and the column boundaries in res.bounds.
// method: "points", "search" or "auto" (both give identical counts).
function countSpikes(v, W, method, rlo, rhi, set = SETS.all, target = res.counts) {
  const R = rhi - rlo;
  const j0 = v.j0 || 0, T0 = v.T0 + j0 * v.dT;         // j0: whole-column shift of the grid (viewGrid)
  const base = Math.ceil(T0), Tend = v.T0 + (j0 + W) * v.dT;
  const i0 = lowerBound(set.ticks, base), i1 = lowerBound(set.ticks, Math.ceil(Tend));

  // integer column boundaries (float64 here), used by pass 1 and for durations in pass 2
  const bnd = new Uint32Array(W + 1);
  for (let j = 0; j <= W; j++) bnd[j] = Math.ceil(v.T0 + (j0 + j) * v.dT);
  if (!res.bounds) res.bounds = gl.createTexture();
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, res.bounds);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32UI, W + 1, 1, 0, gl.RED_INTEGER, gl.UNSIGNED_INT, bnd);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

  // pass 1: counts per (column, row)
  ensureTarget(target, W, R, gl.R32F, gl.RED, gl.FLOAT);
  gl.bindFramebuffer(gl.FRAMEBUFFER, target.fbo);
  gl.viewport(0, 0, W, R);
  gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
  const useSearch = !!res.search && !!set.search && i1 > i0 &&
    (method === "search" || (method === "auto" && i1 - i0 > res.search.ratio * W * R));
  if (useSearch) {
    const S = set.search, P = res.search.prog; gl.useProgram(P.p);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, res.bounds);
    gl.uniform1i(P.u.u_bounds, 0);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, S.rowticks);
    gl.uniform1i(P.u.u_rowticks, 1);
    gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, S.rowoff);
    gl.uniform1i(P.u.u_rowoff, 2);
    gl.uniform1i(P.u.u_rlo, rlo);
    gl.bindVertexArray(res.emptyVao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  } else if (i1 > i0) {
    gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE); gl.blendEquation(gl.FUNC_ADD);
    const P = prog.points; gl.useProgram(P.p);
    gl.uniform1ui(P.u.u_base, base); gl.uniform1f(P.u.u_frac, base - T0);
    gl.uniform1f(P.u.u_dT, v.dT);
    gl.uniform1i(P.u.u_W, W); gl.uniform1i(P.u.u_rlo, rlo); gl.uniform1i(P.u.u_R, R);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, res.bounds);
    gl.uniform1i(P.u.u_bounds, 0);
    gl.bindVertexArray(set.vao);
    gl.drawArrays(gl.POINTS, i0, i1 - i0);
    gl.bindVertexArray(null);
    gl.disable(gl.BLEND);
  }
  return { i0, i1, rlo, R, path: useSearch ? "search" : "points" };
}

// v = {T0, dT (ticks), n0, n1 (Y of layout L)}, grid W x H, norm (bool), eps. Leaves values in
// res.vals.
function computeValues(v, W, H, norm, eps, method = COUNT_MODE, L = lay(), spikes = state.spikes) {
  const SM = SPIKE_MODES[spikes];
  if (SM.A !== "all" || SM.B) ensureBursts();
  // per cell row: which positions it averages and with what end overlaps (float64 here)
  const cr = new Float32Array(4 * H);
  let pmin = Infinity, pmax = -Infinity;
  for (let k = 0; k < H; k++) {
    const sp = cellSpan(v, H, k, L), [pa, pb] = candidates(L, sp);
    if (pb <= pa) { cr[4 * k] = -1; continue; }
    pmin = Math.min(pmin, pa); pmax = Math.max(pmax, pb);
    if (L.frac) {
      const fa = Math.floor(sp.a), fb = Math.floor(sp.b);
      cr[4 * k] = fa; cr[4 * k + 1] = sp.a - fa; cr[4 * k + 2] = fb; cr[4 * k + 3] = sp.b - fb;
    } else {
      cr[4 * k] = pa; cr[4 * k + 1] = sp.a - L.lo[pa];
      cr[4 * k + 2] = pb; cr[4 * k + 3] = sp.b - L.lo[pa];
    }
  }
  // rows to count in pass 1: those shown at the positions any cell uses
  let rlo = 0, rhi = 1;
  if (pmax > pmin) {
    rlo = Infinity; rhi = -Infinity;
    for (let p = pmin; p < pmax; p++) { const r = L.order[p]; if (r < rlo) rlo = r; if (r >= rhi) rhi = r + 1; }
  }
  const c = countSpikes(v, W, method, rlo, rhi, SETS[SM.A], res.counts);
  const c2 = SM.B ? countSpikes(v, W, method, rlo, rhi, SETS[SM.B], res.counts2) : null;
  // pass 2: aggregate rows within the cell's segment
  ensureTarget(res.vals, W, H, gl.RGBA32F, gl.RGBA, gl.FLOAT);
  gl.bindFramebuffer(gl.FRAMEBUFFER, res.vals.fbo);
  gl.viewport(0, 0, W, H);
  const A = prog.agg; gl.useProgram(A.p);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, res.counts.tex);
  gl.uniform1i(A.u.u_counts, 0);
  gl.activeTexture(gl.TEXTURE5); gl.bindTexture(gl.TEXTURE_2D, (c2 ? res.counts2 : res.counts).tex);
  gl.uniform1i(A.u.u_counts2, 5);
  gl.uniform2f(A.u.u_c1, SM.c1[0], SM.c1[1]); gl.uniform2f(A.u.u_c2, SM.c2[0], SM.c2[1]);
  gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, res.rowinfo);
  gl.uniform1i(A.u.u_rowinfo, 1);
  if (!res.cellrows) res.cellrows = gl.createTexture();
  gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, res.cellrows);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, H, 1, 0, gl.RGBA, gl.FLOAT, cr);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.uniform1i(A.u.u_cellrows, 2);
  gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, res.bounds);
  gl.uniform1i(A.u.u_bounds, 3);
  gl.uniform1f(A.u.u_fs, FS); gl.uniform1f(A.u.u_eps, eps);
  gl.activeTexture(gl.TEXTURE4); gl.bindTexture(gl.TEXTURE_2D, L.orderTex);
  gl.uniform1i(A.u.u_order, 4);
  gl.uniform1f(A.u.u_h, L.unit); gl.uniform1i(A.u.u_depth, L.frac ? 0 : 1);
  gl.uniform1i(A.u.u_H, H); gl.uniform1i(A.u.u_rlo, rlo); gl.uniform1i(A.u.u_norm, norm ? 1 : 0);
  gl.bindVertexArray(res.emptyVao);
  gl.drawArrays(gl.TRIANGLES, 0, 3);
  return { ...c, path: c2 && c2.path !== c.path ? `${c.path}+${c2.path}` : c.path,
           spikesInWindow: c.i1 - c.i0 };
}

function drawColormapped(W, H, s, under = null) {
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.viewport(0, 0, W, H);
  gl.disable(gl.BLEND);
  const C = prog.cmap; gl.useProgram(C.p);
  gl.uniform4fv(C.u.u_under, under ? under.rect : [0, 0, 0, 0]);
  gl.uniform1f(C.u.u_ualpha, under ? under.alpha : -1);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, res.vals.tex);
  gl.uniform1i(C.u.u_val, 0);
  const over = state.spikes === "overlay";
  gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, over ? res.lutBurst : res.lut);
  gl.uniform1i(C.u.u_lut, 1);
  gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, res.lutTonic);
  gl.uniform1i(C.u.u_lut2, 2); gl.uniform1i(C.u.u_overlay, over ? 1 : 0);
  gl.uniform1f(C.u.u_lmin, Math.log(s.vmin)); gl.uniform1f(C.u.u_lmax, Math.log(s.vmax));
  gl.bindVertexArray(res.emptyVao);
  gl.drawArrays(gl.TRIANGLES, 0, 3);
}

function readValues(W, H, ch = 0) {   // -> Float32Array H x W of output ch (0 or 2), row 0 = top
  gl.bindFramebuffer(gl.FRAMEBUFFER, res.vals.fbo);
  const px = new Float32Array(W * H * 4);
  gl.readPixels(0, 0, W, H, gl.RGBA, gl.FLOAT, px);
  const out = new Float32Array(W * H);
  for (let y = 0; y < H; y++)
    for (let x = 0; x < W; x++) out[(H - 1 - y) * W + x] = px[4 * (y * W + x) + ch];
  return out;
}

// Span [a, b) of Y that cell row k averages over (float64, same expressions as reference.py):
// the cell's own extent, cut to the segment containing the cell's center.
function cellSpan(v, H, k, L) {
  const dn = (v.n1 - v.n0) / H;
  const y0 = v.n0 + k * dn, y1 = y0 + dn;
  const seg = L.segs[segAt(L, 0.5 * (y0 + y1))];
  return { a: Math.max(y0, seg.y0), b: Math.min(y1, seg.y1, v.n1), seg };
}

// Positions [pa, pb) whose intervals overlap [a, b) by more than zero.
function candidates(L, sp) {
  if (!(sp.b > sp.a)) return [0, 0];
  if (L.frac) return [Math.max(0, Math.floor(sp.a)), Math.min(L.order.length, Math.ceil(sp.b))];
  // depth: bands start in ascending order within a segment; overlap iff lo + h > a and lo < b
  const lo = L.lo, first = (pred) => {      // first i in [r0, r1) with pred(i), else r1
    let x = sp.seg.r0, y = sp.seg.r1;
    while (x < y) { const m = (x + y) >> 1; if (pred(m)) y = m; else x = m + 1; }
    return x;
  };
  return [first((i) => lo[i] + L.unit > sp.a), first((i) => lo[i] >= sp.b)];
}

// CPU value of one cell (float64), the same semantics; used for the hover readout and tests.
function cellCPU(v, W, H, j, k, norm, eps, L = lay(), spikes = state.spikes) {
  const SM = SPIKE_MODES[spikes];
  if (SM.A !== "all" || SM.B) ensureBursts();
  const sp = cellSpan(v, H, k, L), [pa, pb] = candidates(L, sp);
  const n = Math.max(0, pb - pa), ovs = new Float64Array(n);
  const cnt = { all: new Float64Array(n), burst: new Float64Array(n), onset: new Float64Array(n) };
  const fl = burst.flags;
  const slot = new Map();
  let den = 0;
  for (let p = pa; p < pb; p++) {
    const lo = L.lo[p];
    ovs[p - pa] = Math.max(0, Math.min(sp.b, lo + L.unit) - Math.max(sp.a, lo)) / L.unit;
    den += ovs[p - pa]; slot.set(L.order[p], p - pa);
  }
  const jg = (v.j0 || 0) + j, c0 = v.T0 + jg * v.dT, c1 = v.T0 + (jg + 1) * v.dT;
  const lo = lowerBound(ticks, Math.ceil(c0)), hi = lowerBound(ticks, Math.ceil(c1));
  for (let s = lo; s < hi; s++) {
    const q = slot.get(rows[s]);
    if (q === undefined) continue;
    cnt.all[q]++;
    if (fl) { if (fl[s] & 1) cnt.burst[q]++; if (fl[s] & 2) cnt.onset[q]++; }
  }
  const num = [0, 0], tot = { all: 0, burst: 0, onset: 0 };
  for (let p = pa; p < pb; p++) {
    const q = p - pa, w = norm ? 1 / (meta.rows.r0[L.order[p]] + eps) : 1;
    const a = cnt[SM.A][q], b = SM.B ? cnt[SM.B][q] : 0;
    num[0] += ovs[q] * w * (SM.c1[0] * a + SM.c1[1] * b);
    num[1] += ovs[q] * w * (SM.c2[0] * a + SM.c2[1] * b);
    if (ovs[q] > 0) for (const key in tot) tot[key] += cnt[key][q];
  }
  const nt = Math.ceil(c1) - Math.ceil(c0);
  const val = (x) => nt <= 0 ? -1 : den > 0 ? x / (den * nt / FS) : -2;
  return { value: val(num[0]), value2: val(num[1]), neurons: den, spikes: tot.all,
           burstSpikes: fl ? tot.burst : null, onsets: fl ? tot.onset : null,
           a: sp.a, b: sp.b, seg: sp.seg, dtSec: nt / FS };
}

// ---------------------------------------------------------------- layout and view
function computeLayout() {
  const stage = document.getElementById("stage");
  const dpr = window.devicePixelRatio || 1;
  const sw = stage.clientWidth, sh = stage.clientHeight;
  const pw = Math.max(50, sw - L - RM), ph = Math.max(50, sh - TOP - BOT);
  const b = state.block;
  const W = Math.max(4, Math.floor(pw * dpr / b)), H = Math.max(4, Math.floor(ph * dpr / b));
  const x0 = Math.round(L * dpr) / dpr, y0 = Math.round(TOP * dpr) / dpr;  // whole device px
  layout = { dpr, sw, sh, W, H, x: x0, y: y0, w: W * b / dpr, h: H * b / dpr };
  cvs.width = W; cvs.height = H;
  Object.assign(cvs.style, { left: layout.x + "px", top: layout.y + "px",
                             width: layout.w + "px", height: layout.h + "px" });
  ov = document.getElementById("overlay");
  ov.width = Math.round(sw * dpr); ov.height = Math.round(sh * dpr);
  ov.style.width = sw + "px"; ov.style.height = sh + "px";
  ctx = ov.getContext("2d");
}

// The grid of the view. j0 (optional, an integer): the columns are those of the grid starting at
// T0, shifted by j0 columns: column j covers ticks T0 + (j0 + j) * dT onward. While sonify plays
// it moves the view by whole columns this way, so each column keeps exactly the same ticks.
function viewGrid() {
  const g = Sonify.grid();
  if (g) return { ...g, n0: state.n0, n1: state.n1 };
  return { T0: state.t0 * FS, dT: (state.t1 - state.t0) * FS / layout.W,
           n0: state.n0, n1: state.n1 };
}

// A column can't be narrower than one 30 kHz sample, so the shortest window is W samples.
// Vertically, at least 2 neurons' worth (2 rows, or 2 band heights when plotted by depth).
const MIN_N = 2;
function clampView() {
  const Lc = lay(), minT = Math.max(0.001, ((layout.W || 1000) + 0.5) / FS);  // +0.5 guards rounding
  [state.t0, state.t1] = clampAxis(state.t0, state.t1, minT, TS, TE);
  [state.n0, state.n1] = clampAxis(state.n0, state.n1, MIN_N * Lc.unit, 0, Lc.Y);
}
// Keep [a, b) inside [lo, hi) with a span of at least minS. Leaves a valid range untouched, so
// clamping twice changes nothing (the pinned readout relies on that).
function clampAxis(a, b, minS, lo, hi) {
  const tol = 1e-12 * Math.max(1, Math.abs(hi));
  if (b - a >= minS - tol && a >= lo && b <= hi) return [a, b];
  const s = Math.min(Math.max(b - a, minS), hi - lo), c = 0.5 * (a + b);
  const na = Math.max(lo, Math.min(Math.max(c - s / 2, lo), hi - s));
  return [na, Math.min(na + s, hi)];
}
const viewKey = () => [state.t0, state.t1, state.n0, state.n1, state.ymode, state.block,
                       layout.W, layout.H].join();

function render() {
  dirty = false;
  if (cvs.width !== layout.W || cvs.height !== layout.H) computeLayout();
  clampView();     // the shortest allowed window depends on the current grid width
  const v = viewGrid(), s = state[state.mode];
  const t = performance.now();
  const r = computeValues(v, layout.W, layout.H, state.mode === "norm", state.eps);
  drawColormapped(layout.W, layout.H, s, RawUnderlay.draw());
  lastDraw = { ...r, ms: performance.now() - t };
  Sonify.afterRender();    // while playing: schedule the sound of the columns just drawn
  if (pinned && pinned.key !== viewKey()) pinned = null;   // a pin lasts until the view changes
  drawOverlay();
  updateStatus();
  showHover(altHeld && lastMouse ? lastMouse : pinned);
  scheduleHash();
}
function requestRender() { if (!dirty) { dirty = true; requestAnimationFrame(render); } }

// ---------------------------------------------------------------- overlay: axes etc.
function niceStep(span, target) {
  const raw = span / target, p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p;
  return 10 * p;
}
function fmtNum(x, step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return x.toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d });
}
function fmtDur(s) {
  if (s >= 1) return `${+s.toPrecision(4)} s`;
  if (s >= 1e-3) return `${+(s * 1e3).toPrecision(3)} ms`;
  return `${+(s * 1e6).toPrecision(3)} µs`;
}
function fmtCbar(x) { return x >= 1 ? x.toLocaleString("en-US") : String(+x.toPrecision(2)); }

function drawOverlay() {
  const { dpr, x: px, y: py, w: pw, h: ph } = layout;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, layout.sw, layout.sh);
  const ink = "#e6e9ed", ink2 = "#a9b0b8", ink3 = "#6c747d";
  ctx.font = "12px Arial"; ctx.fillStyle = ink2; ctx.strokeStyle = ink3; ctx.lineWidth = 1;

  // overview strip: population rate over the whole session, current window marked
  const pop = meta.pop_rate_1s, pmax = Math.max(...pop.slice(2, -1));
  const ox = (t) => px + (t - TS) / (TE - TS) * pw;
  ctx.fillStyle = "#2f3842";
  ctx.beginPath(); ctx.moveTo(ox(TS), OV_Y + OV_H);
  for (let i = Math.ceil(TS); i < Math.floor(TE); i++)
    ctx.lineTo(ox(i + 0.5), OV_Y + OV_H - Math.min(1, pop[i] / pmax) * OV_H);
  ctx.lineTo(ox(TE), OV_Y + OV_H); ctx.closePath(); ctx.fill();
  const wx0 = ox(state.t0), wx1 = Math.max(ox(state.t1), wx0 + 2);
  ctx.fillStyle = "rgba(78,168,255,0.22)"; ctx.fillRect(wx0, OV_Y, wx1 - wx0, OV_H);
  ctx.strokeStyle = "#4ea8ff"; ctx.strokeRect(wx0 + 0.5, OV_Y + 0.5, wx1 - wx0 - 1, OV_H - 1);
  ctx.fillStyle = ink2; ctx.textAlign = "right"; ctx.textBaseline = "middle";
  ctx.fillText("Population rate,", px - 8, OV_Y + 8);
  ctx.fillText("whole session", px - 8, OV_Y + 22);
  RawUnderlay.drawLegend(ctx);
  ctx.font = "12px Arial"; ctx.lineWidth = 1;

  // time axis
  const span = state.t1 - state.t0, step = niceStep(span, Math.max(2, pw / 110));
  ctx.strokeStyle = ink3; ctx.beginPath();
  ctx.moveTo(px, py + ph + 0.5); ctx.lineTo(px + pw, py + ph + 0.5); ctx.stroke();
  ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillStyle = ink2;
  for (let t = Math.ceil(state.t0 / step) * step; t <= state.t1 + 1e-12; t += step) {
    const xx = px + (t - state.t0) / span * pw - Sonify.shiftCss;   // moves with the raster while playing
    if (xx < px - 0.5) continue;
    ctx.beginPath(); ctx.moveTo(xx, py + ph); ctx.lineTo(xx, py + ph + 5); ctx.stroke();
    ctx.fillText(fmtNum(t, step), xx, py + ph + 8);
  }
  ctx.fillStyle = ink; ctx.fillText("Time (s)", px + pw / 2, py + ph + 28);

  // neuron axis
  const L = lay(), nspan = state.n1 - state.n0;
  const yOf = (n) => py + (n - state.n0) / nspan * ph;
  hitboxes = [];
  ctx.save(); ctx.beginPath(); ctx.rect(0, py, px, ph); ctx.clip();
  if (L.name === "rastermap") {
    // probe of each neuron as a color strip, plus position ticks
    const p0 = Math.max(0, Math.floor(state.n0)), p1 = Math.min(N, Math.ceil(state.n1));
    for (let p = p0; p < p1; p++) {
      const ya = yOf(p), yb = yOf(p + 1);
      ctx.fillStyle = meta.probe_hues[meta.rows.probe[L.order[p]]];
      ctx.fillRect(px - 12, ya, 8, Math.max(0.6, yb - ya - (yb - ya > 3 ? 1 : 0)));
    }
    const st = niceStep(nspan, Math.max(2, ph / 60));
    ctx.fillStyle = ink2; ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (let q = Math.ceil(state.n0 / st) * st; q <= state.n1; q += st)
      ctx.fillText(q.toLocaleString("en-US"), px - 16, yOf(q));
  } else {
    // probe and shank bands with labels; clicking a label (or band) zooms to it
    const probeExt = {};
    for (const sg of L.segs) {
      if (sg.y1 <= state.n0 || sg.y0 >= state.n1) continue;
      const ya = Math.max(py, yOf(sg.y0)), yb = Math.min(py + ph, yOf(sg.y1));
      ctx.globalAlpha = sg.shank % 2 ? 0.55 : 1.0;
      ctx.fillStyle = meta.probe_hues[sg.probe];
      ctx.fillRect(px - 12, ya, 8, Math.max(1, yb - ya - 1));
      ctx.globalAlpha = 1;
      hitboxes.push({ x0: px - 30, x1: px, y0: ya, y1: yb, n0: sg.y0, n1: sg.y1 });
      if (yb - ya >= 13) {
        ctx.fillStyle = ink2; ctx.textAlign = "right"; ctx.textBaseline = "middle";
        ctx.fillText(`S${sg.shank}`, px - 16, (ya + yb) / 2);
      }
      const e = probeExt[sg.probe] || (probeExt[sg.probe] = [ya, yb, sg.y0, sg.y1]);
      e[0] = Math.min(e[0], ya); e[1] = Math.max(e[1], yb);
    }
    ctx.font = "13px Arial";
    for (const [p, [ya, yb]] of Object.entries(probeExt)) {
      if (yb - ya < 15) continue;
      ctx.fillStyle = meta.probe_hues[p]; ctx.textAlign = "right"; ctx.textBaseline = "middle";
      const txt = `Probe ${+p + 1}`, tw = ctx.measureText(txt).width, yc = (ya + yb) / 2;
      ctx.fillText(txt, px - 40, yc);
      const segs = L.segs.filter((g) => g.probe === +p);
      hitboxes.push({ x0: px - 44 - tw, x1: px - 36, y0: yc - 9, y1: yc + 9,
                      n0: segs[0].y0, n1: segs[segs.length - 1].y1 });
    }
  }
  ctx.restore();
  ctx.save(); ctx.translate(16, py + ph / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = ink; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.font = "12px Arial";
  ctx.fillText(L.name === "rows" ? "Neuron (probe, shank, depth: shallow at top)"
    : L.name === "depth" ? `Position along shank, to scale (${DEPTH_H} µm per neuron; shallow at top)`
    : "Neuron, " + (L.axis || "Rastermap order"), 0, 0);
  ctx.restore();

  // pinned tooltip: outline the cell it describes
  if (pinned) {
    const { W, H } = layout, j = Math.floor((pinned.x - px) / pw * W), k = Math.floor((pinned.y - py) / ph * H);
    ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 1;
    ctx.strokeRect(px + j * pw / W - 1.5, py + k * ph / H - 1.5, pw / W + 3, ph / H + 3);
  }

  // colorbar
  const s = state[state.mode], cx = px + pw + 18, cw = 14, over = state.spikes === "overlay";
  const lut = meta.colormaps[s.cmap];
  for (let i = 0; i < ph; i++) {
    const t = 1 - (i + 0.5) / ph, q = Math.round(t * 255);
    if (over) {     // burst (magenta) and tonic (green) ramps side by side
      ctx.fillStyle = `rgb(${q},0,${q})`; ctx.fillRect(cx, py + i, cw / 2, 1.5);
      ctx.fillStyle = `rgb(0,${q},0)`; ctx.fillRect(cx + cw / 2, py + i, cw / 2, 1.5);
    } else {
      ctx.fillStyle = `rgb(${lut[3 * q]},${lut[3 * q + 1]},${lut[3 * q + 2]})`;
      ctx.fillRect(cx, py + i, cw, 1.5);
    }
  }
  if (over) {
    ctx.font = "10px Arial"; ctx.textAlign = "center"; ctx.textBaseline = "bottom";
    ctx.fillStyle = "#e05ad0"; ctx.fillText("B", cx + cw / 4, py - 2);
    ctx.fillStyle = "#5ad05a"; ctx.fillText("T", cx + 3 * cw / 4, py - 2);
  }
  ctx.font = "11px Arial"; ctx.fillStyle = ink2; ctx.textAlign = "left"; ctx.textBaseline = "middle";
  const l0 = Math.log10(s.vmin), l1 = Math.log10(s.vmax);
  const mult = (l1 - l0) > 2.5 ? [1] : (l1 - l0) > 1.2 ? [1, 3] : [1, 2, 5];
  for (let d = Math.floor(l0); d <= Math.ceil(l1); d++)
    for (const m of mult) {
      const val = m * Math.pow(10, d), lv = Math.log10(val);
      if (lv < l0 - 1e-9 || lv > l1 + 1e-9) continue;
      const yy = py + ph - (lv - l0) / (l1 - l0) * ph;
      ctx.strokeStyle = ink2; ctx.beginPath(); ctx.moveTo(cx + cw, yy); ctx.lineTo(cx + cw + 4, yy); ctx.stroke();
      ctx.fillText(fmtCbar(val), cx + cw + 7, yy);
    }
  ctx.save(); ctx.translate(cx + cw + 62, py + ph / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = ink; ctx.textAlign = "center"; ctx.font = "12px Arial";
  const what = { all: "firing rate", tonic: "tonic firing rate", burst: "burst-spike rate",
                 onset: "burst rate", overlay: "burst (B) and tonic (T) rates" }[state.spikes];
  const unit = state.spikes === "onset" ? "bursts/s" : "spikes/s";
  ctx.fillText(state.mode === "abs" ? `Mean ${what} (${unit})`
               : `${what[0].toUpperCase() + what.slice(1)} / (session mean + ${state.eps} spikes/s)`, 0, 0);
  ctx.restore();
  Sonify.drawHead(ctx);
}

let testNote = null;      // the result of a self-test (?test=1, ?rawtest=1, ?sonifytest=1): kept in the status line
function updateStatus() {
  if (testNote) { document.getElementById("status").textContent = testNote; return; }
  const { W, H, dpr } = layout, b = state.block;
  const dt = (state.t1 - state.t0) / W, dn = (state.n1 - state.n0) / H;
  const nsp = lastDraw.i1 - lastDraw.i0;
  document.getElementById("status").textContent =
    `${state.t0.toFixed(3)} to ${state.t1.toFixed(3)} s (${fmtDur(state.t1 - state.t0)}) · ` +
    (state.ymode === "depth"
      ? `${state.n0.toFixed(0)} to ${state.n1.toFixed(0)} µm of stacked shanks · ` +
        `1 cell = ${b}×${b} px = ${dn.toPrecision(3)} µm × ${fmtDur(dt)} · `
      : `${lay().name === "rastermap" ? (lay().key ? "positions" : "Rastermap positions") : "rows"} ${state.n0.toFixed(1)} to ` +
        `${state.n1.toFixed(1)} (${(state.n1 - state.n0).toFixed(1)} neurons) · ` +
        `1 cell = ${b}×${b} px = ${dn.toPrecision(3)} neurons × ${fmtDur(dt)} · `) +
    `grid ${W}×${H} · ${nsp >= 1e5 ? (nsp / 1e6).toFixed(2) + " M" : nsp.toLocaleString("en-US")} ` +
    `${{ all: "spikes", tonic: "spikes", overlay: "spikes", burst: "burst spikes", onset: "bursts" }[state.spikes]} in window · ` +
    (state.spikes !== "all" && burst.flags ? `bursts: ${SETS.onset.n.toLocaleString("en-US")} ` +
      `(${(100 * SETS.burst.n / ticks.length).toFixed(1)}% of spikes) at ≥${burst.params.S} ms silence, ` +
      `ISI ≤${burst.params.M} ms, ≥${burst.params.K} spikes · ` : "") +
    `counted by ${lastDraw.path} · ${lastDraw.ms.toFixed(1)} ms to submit` + Sonify.status() + RawUnderlay.status();
}

// ---------------------------------------------------------------- hover
function showHover(m) {
  const el = document.getElementById("hover");
  const { x: px, y: py, w: pw, h: ph, W, H } = layout;
  if (!m) { el.style.display = "none"; return; }
  const lx = m.x - px, ly = m.y - py;
  if (lx < 0 || ly < 0 || lx >= pw || ly >= ph) { el.style.display = "none"; return; }
  const j = Math.min(W - 1, Math.floor(lx / pw * W)), k = Math.min(H - 1, Math.floor(ly / ph * H));
  const v = viewGrid(), norm = state.mode === "norm", L = lay();
  const c = cellCPU(v, W, H, j, k, norm, state.eps, L);
  const t = state.t0 + lx / pw * (state.t1 - state.t0);
  const y = state.n0 + ly / ph * (state.n1 - state.n0);
  const R = meta.rows;
  const desc = (row) => `cluster ${R.cluster[row]}, Probe ${R.probe[row] + 1} shank ${R.shank[row]}, ` +
    `${R.y_um[row].toFixed(0)} µm above deepest site, session mean ${R.r0[row].toFixed(2)} spikes/s`;
  let neuronTxt;
  // the neuron line names a neuron from the segment the cell's value uses
  const yc = Math.min(Math.max(y, c.seg.y0), c.seg.y1 - 1e-9);
  if (L.frac) {
    const p = Math.min(N - 1, Math.floor(yc)), row = L.order[p];
    neuronTxt = (L.name === "rastermap" ? `position ${p}, row ${row}, ` : `row ${row}, `) + desc(row);
  } else {
    // the units whose bands contain the pointer; name the one centered nearest
    const seg = c.seg, here = [];
    for (let i = seg.r0; i < seg.r1; i++) if (L.lo[i] <= yc && yc < L.lo[i] + L.unit) here.push(i);
    const yum = seg.top10 / 10 - (yc - seg.y0);
    if (!here.length) neuronTxt = `none at ${yum.toFixed(0)} µm above deepest site (Probe ${seg.probe + 1} shank ${seg.shank})`;
    else {
      here.sort((a, b) => Math.abs(L.lo[a] + L.unit / 2 - yc) - Math.abs(L.lo[b] + L.unit / 2 - yc));
      neuronTxt = `row ${here[0]}, ${desc(here[0])}` +
        (here.length > 1 ? ` (+${here.length - 1} more overlapping here)` : "");
    }
  }
  const unit = state.spikes === "onset" ? "bursts/s" : "spikes/s";
  const fmtV = (x) => norm ? `${x.toPrecision(3)} (rate / (session mean + ε))` : `${x.toPrecision(4)} ${unit}`;
  const valTxt = c.value === -1 ? "no samples in this column" : c.value === -2 ? "no neurons in this cell"
    : state.spikes === "overlay" ? `burst ${fmtV(c.value)}, tonic ${fmtV(c.value2)}`
    : `${SPIKE_MODES[state.spikes].label.toLowerCase()}: ${fmtV(c.value)}`;
  const burstTxt = c.burstSpikes === null ? "" :
    ` (${c.burstSpikes} in ${c.onsets} burst${c.onsets === 1 ? "" : "s"}, ${c.spikes - c.burstSpikes} tonic)`;
  el.innerHTML =
    `<span class="k">Time</span> ${t.toFixed(4)} s<br>` +
    `<span class="k">Neuron</span> ${neuronTxt}<br>` +
    `<span class="k">Cell</span> ${c.spikes} spikes${burstTxt} over ${c.neurons.toPrecision(3)} neurons × ${fmtDur(c.dtSec)}<br>` +
    `<span class="k">Value</span> ${valTxt}`;
  el.style.display = "block";
  const bx = m.x + 16 + el.offsetWidth > layout.sw ? m.x - el.offsetWidth - 12 : m.x + 16;
  const by = m.y + 16 + el.offsetHeight > layout.sh ? m.y - el.offsetHeight - 12 : m.y + 16;
  el.style.left = bx + "px"; el.style.top = by + "px";
}

// ---------------------------------------------------------------- interaction
function zoomTime(f, tc) {
  state.t0 = tc - (tc - state.t0) * f; state.t1 = tc + (state.t1 - tc) * f;
  clampView(); requestRender();
}
function zoomNeurons(f, nc) {
  state.n0 = nc - (nc - state.n0) * f; state.n1 = nc + (state.n1 - nc) * f;
  clampView(); requestRender();
}

function setupInteraction() {
  const stage = document.getElementById("stage"), box = document.getElementById("box");
  const toLocal = (e) => { const r = stage.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };
  const tAt = (x) => state.t0 + (x - layout.x) / layout.w * (state.t1 - state.t0);
  const nAt = (y) => state.n0 + (y - layout.y) / layout.h * (state.n1 - state.n0);

  const wheelAxis = { t: -1e9, horiz: false };
  let inGesture = false;
  stage.addEventListener("wheel", (e) => {
    e.preventDefault();
    const m = toLocal(e);
    const unit = e.deltaMode === 1 ? 40 : e.deltaMode === 2 ? 400 : 1;
    const dx = e.deltaX * unit, dy = e.deltaY * unit;
    // Horizontal scrolling (two-finger sideways swipe on a trackpad, tilt wheel) pans time.
    // Shift + wheel stays neuron zoom: macOS (and Chrome on Windows) deliver it as deltaX.
    // Trackpad scrolls are often diagonal, so the axis chosen by a gesture's first event sticks
    // until the gesture pauses (150 ms); otherwise one swipe would alternate pan and zoom.
    const now = performance.now();
    if (now - wheelAxis.t > 150) wheelAxis.horiz = !e.shiftKey && Math.abs(dx) > Math.abs(dy);
    wheelAxis.t = now;
    if (wheelAxis.horiz && !e.shiftKey && !e.ctrlKey) {
      const s = state.t1 - state.t0, dt = dx / layout.w * s;
      state.t0 += dt; state.t1 += dt; clampView(); requestRender(); return;
    }
    if (e.ctrlKey && inGesture) return;    // Safari: the gesture handlers already zoom
    let d = dy || dx;
    // A trackpad pinch arrives as Ctrl + wheel with small deltas (about 100 x log of the
    // pinch scale in Chrome), so on macOS it gets a higher gain. Ctrl + mouse wheel on Windows
    // is unchanged.
    if (IS_MAC && e.ctrlKey) d *= 6;
    const f = Math.exp(d * 0.0015);
    const inLeft = m.x < layout.x, inPlotX = m.x >= layout.x && m.x <= layout.x + layout.w;
    if (inLeft || (e.shiftKey && inPlotX))
      zoomNeurons(f, Math.min(Math.max(nAt(m.y), state.n0), state.n1));
    else zoomTime(f, Math.min(Math.max(tAt(m.x), state.t0), state.t1));
  }, { passive: false });

  // Safari sends trackpad pinches as its own gesture events (not Ctrl + wheel) and would
  // otherwise zoom the whole page. Pinching zooms time around the pointer, like the wheel.
  // Listening on the document also stops a pinch over the top bar from zooming the page.
  let gestureScale = 1;
  document.addEventListener("gesturestart", (e) => { e.preventDefault(); gestureScale = 1; inGesture = true; });
  document.addEventListener("gesturechange", (e) => {
    e.preventDefault();
    if (!(e.scale > 0)) return;
    const m = toLocal(e), f = gestureScale / e.scale;
    gestureScale = e.scale;
    if (m.y < 0 || m.y > layout.sh) return;
    zoomTime(f, Math.min(Math.max(tAt(m.x), state.t0), state.t1));
  });
  document.addEventListener("gestureend", (e) => { e.preventDefault(); inGesture = false; });

  let drag = null;
  stage.addEventListener("contextmenu", (e) => e.preventDefault());
  stage.addEventListener("mousedown", (e) => {
    const m = toLocal(e);
    // preventDefault below keeps focus where it was; release it so the keys work
    if (document.activeElement && document.activeElement !== document.body) document.activeElement.blur();
    const hb = e.button === 0 && hitAt(m);
    if (hb) { state.n0 = hb.n0; state.n1 = hb.n1; clampView(); requestRender(); e.preventDefault(); return; }
    if (m.y >= OV_Y && m.y <= OV_Y + OV_H && m.x >= layout.x && m.x <= layout.x + layout.w) {
      drag = { kind: "overview" }; jumpOverview(m.x); return;
    }
    // Ctrl + drag, Cmd + drag (macOS) or right-drag
    const boxMode = e.ctrlKey || (IS_MAC && e.metaKey) || e.button === 2;
    drag = { kind: boxMode ? "box" : "pan", x: m.x, y: m.y,
             t0: state.t0, t1: state.t1, n0: state.n0, n1: state.n1 };
    e.preventDefault();
  });
  const hitAt = (m) => hitboxes.find((h) => m.x >= h.x0 && m.x <= h.x1 && m.y >= h.y0 && m.y <= h.y1);
  window.addEventListener("mousemove", (e) => {
    const m = toLocal(e);
    lastMouse = m; altHeld = e.altKey;
    stage.style.cursor = !drag && hitAt(m) ? "pointer" : "";
    if (drag) drag.maxd = Math.max(drag.maxd || 0, Math.hypot(m.x - drag.x, m.y - drag.y));
    if (drag && drag.kind === "pan" && Sonify.playing && drag.maxd < 4) {
      // a click while the view scrolls: nothing to pan yet
    } else if (drag && drag.kind === "pan") {
      const dt = -(m.x - drag.x) / layout.w * (drag.t1 - drag.t0);
      const dn = -(m.y - drag.y) / layout.h * (drag.n1 - drag.n0);
      state.t0 = drag.t0 + dt; state.t1 = drag.t1 + dt;
      state.n0 = drag.n0 + dn; state.n1 = drag.n1 + dn;
      clampView(); requestRender();
    } else if (drag && drag.kind === "box") {
      Object.assign(box.style, { display: "block", left: Math.min(m.x, drag.x) + "px",
        top: Math.min(m.y, drag.y) + "px", width: Math.abs(m.x - drag.x) + "px",
        height: Math.abs(m.y - drag.y) + "px" });
    } else if (drag && drag.kind === "overview") {
      jumpOverview(m.x);
    } else showHover(altHeld ? m : pinned);
  });
  window.addEventListener("mouseup", (e) => {
    if (drag && drag.kind === "box") {
      box.style.display = "none";
      const m = toLocal(e);
      if (Math.abs(m.x - drag.x) > 4) {
        const a = tAt(Math.min(m.x, drag.x)), b = tAt(Math.max(m.x, drag.x));
        state.t0 = a; state.t1 = b;
      }
      if (Math.abs(m.y - drag.y) > 4) {
        const a = nAt(Math.min(m.y, drag.y)), b = nAt(Math.max(m.y, drag.y));
        state.n0 = a; state.n1 = b;
      }
      clampView(); requestRender();
    } else if (drag && drag.kind === "pan") {
      // a click (no real drag) pins the tooltip to that cell, until the view changes
      const m = toLocal(e);
      if (Sonify.playing) {
        // no pinned readout while the view scrolls (it would be cleared at the next step)
      } else if (Math.max(drag.maxd || 0, Math.hypot(m.x - drag.x, m.y - drag.y)) < 4) {
        Object.assign(state, { t0: drag.t0, t1: drag.t1, n0: drag.n0, n1: drag.n1 });
        clampView();
        const inPlot = m.x >= layout.x && m.x < layout.x + layout.w &&
                       m.y >= layout.y && m.y < layout.y + layout.h;
        pinned = inPlot ? { x: m.x, y: m.y, key: viewKey() } : null;
        requestRender();
      }
    }
    drag = null;
  });
  stage.addEventListener("mouseleave", () => { lastMouse = null; showHover(pinned); });
  // Alt (Option on a Mac) held: the tooltip follows the pointer
  const setAlt = (on) => { altHeld = on; showHover(on && lastMouse ? lastMouse : pinned); };
  window.addEventListener("keydown", (e) => { if (e.key === "Alt") { e.preventDefault(); setAlt(true); } });
  window.addEventListener("keyup", (e) => { if (e.key === "Alt") { e.preventDefault(); setAlt(false); } });
  window.addEventListener("blur", () => setAlt(false));

  function jumpOverview(x) {
    const t = TS + (x - layout.x) / layout.w * (TE - TS), s = state.t1 - state.t0;
    state.t0 = t - s / 2; state.t1 = t + s / 2; clampView(); requestRender();
  }

  window.addEventListener("keydown", (e) => {
    if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;
    if (e.ctrlKey || e.metaKey) return;        // Ctrl / Cmd combinations belong to the browser
    const ts = state.t1 - state.t0, ns = state.n1 - state.n0;
    const k = e.key;
    if (k === "ArrowLeft") { state.t0 -= ts * 0.2; state.t1 -= ts * 0.2; }
    else if (k === "ArrowRight") { state.t0 += ts * 0.2; state.t1 += ts * 0.2; }
    else if (k === "ArrowUp") { state.n0 -= ns * 0.2; state.n1 -= ns * 0.2; }
    else if (k === "ArrowDown") { state.n0 += ns * 0.2; state.n1 += ns * 0.2; }
    else if (k === "+" || k === "=") return zoomTime(0.5, 0.5 * (state.t0 + state.t1));
    else if (k === "-" || k === "_") return zoomTime(2, 0.5 * (state.t0 + state.t1));
    else if (k === "r" || k === "R") resetView();
    else if (k === "a" || k === "A") { state.n0 = 0; state.n1 = lay().Y; }
    else if (k === "Escape") { pinned = null; }
    else if (k === "w" || k === "W") { state.t0 = TS; state.t1 = TE; }
    else if ((k === "s" || k === "S" || k === " ") && !e.altKey) {
      e.preventDefault(); if (!e.repeat) Sonify.toggle(); return;
    }
    else if (RawUnderlay.key(k)) { /* raw voltage underlay: U, O, M, V, [, ] */ }
    else return;
    e.preventDefault(); clampView(); requestRender();
  });

  window.addEventListener("resize", () => { computeLayout(); requestRender(); });
}

function resetView() { state.t0 = TS; state.t1 = TE; state.n0 = 0; state.n1 = lay().Y; }

// ---------------------------------------------------------------- controls and URL state
function syncControls() {
  const s = state[state.mode];
  document.getElementById("mode").value = state.mode;
  document.getElementById("eps").value = state.eps;
  document.getElementById("eps").disabled = state.mode !== "norm";
  document.getElementById("epsWrap").style.opacity = state.mode === "norm" ? 1 : 0.45;
  document.getElementById("block").value = String(state.block);
  document.getElementById("ymode").value = state.ymode;
  document.getElementById("spikes").value = state.spikes;
  document.getElementById("burstWrap").style.display = state.spikes === "all" ? "none" : "flex";
  document.getElementById("bS").value = state.burst.S;
  document.getElementById("bM").value = state.burst.M;
  document.getElementById("bK").value = state.burst.K;
  document.getElementById("vmin").value = s.vmin;
  document.getElementById("vmax").value = s.vmax;
  document.getElementById("cmap").value = s.cmap;
  setLUT(s.cmap);
}

function setupControls() {
  const cm = document.getElementById("cmap");
  for (const name of Object.keys(meta.colormaps)) {
    const o = document.createElement("option"); o.value = o.textContent = name; cm.appendChild(o);
  }
  const on = (id, ev, fn) => document.getElementById(id).addEventListener(ev, fn);
  // a select keeps keyboard focus after a choice, which would swallow the arrow keys
  for (const el of document.querySelectorAll("#bar select"))
    el.addEventListener("change", () => el.blur());
  on("mode", "change", (e) => { state.mode = e.target.value; syncControls(); requestRender(); });
  on("eps", "change", (e) => {
    const x = parseFloat(e.target.value); if (x >= 0) state.eps = x; syncControls(); requestRender();
  });
  on("block", "change", (e) => { state.block = +e.target.value; computeLayout(); requestRender(); });
  for (const id of ["vmin", "vmax"]) on(id, "change", (e) => {
    const x = parseFloat(e.target.value), s = state[state.mode];
    if (x > 0) s[id] = x;
    if (s.vmin >= s.vmax) s.vmax = s.vmin * 10;
    syncControls(); requestRender();
  });
  on("cmap", "change", (e) => { state[state.mode].cmap = e.target.value; syncControls(); requestRender(); });
  on("reset", "click", () => { resetView(); requestRender(); });
  on("allN", "click", () => { state.n0 = 0; state.n1 = lay().Y; requestRender(); });
  if (!LAY.rastermap) document.querySelector('#ymode option[value="rastermap"]').disabled = true;
  // Burst detection takes about 1 s; say so, and let the message paint before it runs.
  const renderAfterDetect = () => {
    const p = state.burst, q = burst.params;
    const needs = state.spikes !== "all" && (!burst.flags || p.S !== q.S || p.M !== q.M || p.K !== q.K);
    if (!needs) return requestRender();
    document.getElementById("status").textContent = "Detecting bursts...";
    requestAnimationFrame(() => setTimeout(() => { ensureBursts(); requestRender(); }, 0));
  };
  on("spikes", "change", (e) => { state.spikes = e.target.value; syncControls(); renderAfterDetect(); });
  for (const [id, key, lo, hi, intg] of [["bS", "S", 0, 10000, false], ["bM", "M", 0.1, 1000, false],
                                          ["bK", "K", 2, 16, true]])
    on(id, "change", (e) => {
      let x = parseFloat(e.target.value);
      if (Number.isFinite(x)) state.burst[key] = Math.min(hi, Math.max(lo, intg ? Math.round(x) : x));
      syncControls(); renderAfterDetect();
    });
  on("ymode", "change", (e) => {     // a new vertical coordinate: show all neurons
    state.ymode = e.target.value; state.n0 = 0; state.n1 = lay().Y; pinned = null;
    clampView(); requestRender();
  });
  on("allT", "click", () => { state.t0 = TS; state.t1 = TE; requestRender(); });
  on("helpBtn", "click", () => {
    const h = document.getElementById("help"); h.style.display = h.style.display === "block" ? "none" : "block";
  });
}

let hashTimer = null;
function scheduleHash() {
  clearTimeout(hashTimer);
  hashTimer = setTimeout(() => {
    const s = state[state.mode];
    const p = new URLSearchParams({ t0: state.t0.toFixed(6), t1: state.t1.toFixed(6),
      n0: state.n0.toFixed(3), n1: state.n1.toFixed(3), ymode: state.ymode, spikes: state.spikes,
      bS: state.burst.S, bM: state.burst.M, bK: state.burst.K, mode: state.mode, eps: state.eps,
      block: state.block, vmin: s.vmin, vmax: s.vmax, cmap: s.cmap });
    RawUnderlay.toHash(p);
    Sonify.toHash(p);
    history.replaceState(null, "", location.pathname + location.search + "#" + p);
  }, 300);
}
function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  RawUnderlay.fromHash(p);
  if (Object.hasOwn(LAY, p.get("ymode") || "")) { state.ymode = p.get("ymode"); state.n0 = 0; state.n1 = lay().Y; }
  if (Object.hasOwn(SPIKE_MODES, p.get("spikes") || "")) state.spikes = p.get("spikes");
  for (const [k, lo, hi] of [["S", 0, 10000], ["M", 0.1, 1000], ["K", 2, 16]]) {
    const x = parseFloat(p.get("b" + k));
    if (Number.isFinite(x)) state.burst[k] = Math.min(hi, Math.max(lo, k === "K" ? Math.round(x) : x));
  }
  for (const k of ["t0", "t1", "n0", "n1"]) {
    const x = parseFloat(p.get(k)); if (Number.isFinite(x)) state[k] = x;
  }
  const e = parseFloat(p.get("eps")); if (Number.isFinite(e) && e >= 0) state.eps = e;
  if (p.has("mode")) state.mode = p.get("mode") === "norm" ? "norm" : "abs";
  if (p.has("block")) state.block = Math.max(1, parseInt(p.get("block")));
  const s = state[state.mode];
  const a = parseFloat(p.get("vmin")), b = parseFloat(p.get("vmax"));
  if (a > 0) s.vmin = a;
  if (b > s.vmin) s.vmax = b; else if (s.vmax <= s.vmin) s.vmax = s.vmin * 10;
  if (p.has("cmap") && meta.colormaps[p.get("cmap")]) s.cmap = p.get("cmap");
  Sonify.fromHash(p);
}

// ---------------------------------------------------------------- self-test (?test=1)
// Compares GPU cell values against reference.py output written by make_refs.py, once per
// pass-1 counting method (points, and search when available).
async function runTests() {
  const views = await (await fetch("refs/views.json")).json();
  const methods = res.search ? ["points", "search"] : ["points"];
  const report = [], savedBurst = { ...state.burst };
  // burst detection: the same spikes as bursts/detect_bursts.py (counts and index sums)
  const bj = await (await fetch("refs/bursts.json")).json();
  state.burst = { S: bj.params[0], M: bj.params[1], K: bj.params[2] };
  ensureBursts();
  const burstCheck = { n_burst: [SETS.burst.n, bj.n_burst], n_onset: [SETS.onset.n, bj.n_onset],
                       sum_index_burst: [SETS.burst.sumIdx, bj.sum_index_burst],
                       sum_index_onset: [SETS.onset.sumIdx, bj.sum_index_onset] };
  burstCheck.match = Object.values(burstCheck).every(([a, b]) => a === b);
  const compare = (g, ref) => {
    let maxRel = 0, nBad = 0, sumG = 0, sumR = 0, maxAbs = 0;
    for (let i = 0; i < g.length; i++) {
      const d = Math.abs(g[i] - ref[i]), rel = d / Math.max(Math.abs(ref[i]), 1e-12);
      sumG += g[i]; sumR += ref[i];
      if (d > maxAbs) maxAbs = d;
      if (!(ref[i] === 0 && g[i] === 0) && rel > maxRel) maxRel = rel;
      if (rel > 1e-3 && d > 1e-6) nBad++;
    }
    return { nBad, maxRel, maxAbs, sumRatio: sumG / sumR };
  };
  for (const tv of views) {
    const ref = new Float32Array(await (await fetch("refs/" + tv.file)).arrayBuffer());
    const ref2 = tv.file2 ? new Float32Array(await (await fetch("refs/" + tv.file2)).arrayBuffer()) : null;
    const spikes = tv.spikes || "all";
    if (tv.burst_params) { state.burst = { S: tv.burst_params[0], M: tv.burst_params[1], K: tv.burst_params[2] }; ensureBursts(); }
    const v = { T0: tv.T0, dT: tv.dT, n0: tv.n0, n1: tv.n1 };
    const L = LAY[tv.layout || "rows"];
    if (!L) { report.push({ name: tv.name, skipped: "layout not available" }); continue; }
    if (tv.h !== undefined && tv.h !== L.unit && !L.frac) throw new Error("reference used another band height");
    let pointCounts = null;
    for (const method of methods) {
      const cs = computeValues(v, tv.W, tv.H, tv.norm, tv.eps, method, L, spikes);
      const g = readValues(tv.W, tv.H), cmp2 = ref2 ? compare(readValues(tv.W, tv.H, 2), ref2) : null;
      // raw pass-1 counts must be bitwise identical between methods (integers in float32)
      const twoSets = !!SPIKE_MODES[spikes].B;
      const counts = readCounts(tv.W, cs.R), counts2 = twoSets ? readCounts(tv.W, cs.R, res.counts2) : null;
      let countsDiffer = 0;
      if (!pointCounts) pointCounts = [counts, counts2];
      else {
        for (let i = 0; i < counts.length; i++) if (counts[i] !== pointCounts[0][i]) countsDiffer++;
        if (twoSets) for (let i = 0; i < counts2.length; i++) if (counts2[i] !== pointCounts[1][i]) countsDiffer++;
      }
      let { nBad, maxRel, maxAbs, sumRatio } = compare(g, ref);
      if (cmp2) { nBad += cmp2.nBad; maxRel = Math.max(maxRel, cmp2.maxRel); maxAbs = Math.max(maxAbs, cmp2.maxAbs); }
      // CPU path (hover) at 40 random cells
      let cpuMax = 0, seed = 12345;
      const rnd = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296);
      for (let q = 0; q < 40; q++) {
        const j = Math.floor(rnd() * tv.W), k = Math.floor(rnd() * tv.H);
        const cc = cellCPU(v, tv.W, tv.H, j, k, tv.norm, tv.eps, L, spikes), r = ref[k * tv.W + j];
        cpuMax = Math.max(cpuMax, Math.abs(cc.value - r) / Math.max(Math.abs(r), 1e-12));
        if (ref2) {
          const r2 = ref2[k * tv.W + j];
          cpuMax = Math.max(cpuMax, Math.abs(cc.value2 - r2) / Math.max(Math.abs(r2), 1e-12));
        }
      }
      report.push({ name: tv.name, layout: L.name, spikes, method, path: cs.path, countsDiffer, W: tv.W, H: tv.H, cells: g.length, badCells: nBad,
                    fracBad: nBad / g.length, maxRel, maxAbs, sumRatio, cpuMaxRel: cpuMax });
    }
  }
  state.burst = savedBurst;
  if (burst.flags) ensureBursts();      // back to the user's burst definition
  // timing of the whole-session, all-neuron view at the current grid size, per method
  const W = layout.W, H = layout.H, fullViewMs = {};
  const v = { T0: TS * FS, dT: (TE - TS) * FS / W, n0: 0, n1: N };
  const px = new Float32Array(4);
  for (const method of methods) {
    const times = (fullViewMs[method] = []);
    for (let r = 0; r < 6; r++) {
      const t = performance.now();
      computeValues(v, W, H, false, 0.1, method, LAY.rows, "all");
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, px);   // forces completion
      times.push(performance.now() - t);
    }
  }
  const out = { report, burstCheck, fullViewMs, grid: [W, H], floatBlend: res.floatBlend,
                search: !!res.search, searchNote: res.searchNote,
                searchRatio: res.search && res.search.ratio, calib: res.search && res.search.calib,
                renderer: gl.getParameter(gl.RENDERER) };
  const dbg = gl.getExtension("WEBGL_debug_renderer_info");
  if (dbg) out.renderer = gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL);
  window.__testResults = out;
  console.log(JSON.stringify(out, null, 1));
  return out;
}

// ---------------------------------------------------------------- main
(async function main() {
  try {
    await load();
    buildLayouts();
    initGL();
    resetView();
    readHash();
    clampView();
    setupControls();
    Sonify.setup();
    syncControls();
    computeLayout();
    setupInteraction();
    if (state.spikes !== "all") {        // a burst mode from the URL: detect before showing
      document.getElementById("loadMsg").textContent = "Detecting bursts";
      await new Promise((r) => requestAnimationFrame(() => setTimeout(r, 0)));
      ensureBursts();
    }
    document.getElementById("loading").style.display = "none";
    render();
    if (new URLSearchParams(location.search).has("rawtest")) {
      const r = window.__rawTest = await RawUnderlay.test();
      console.log("raw underlay self-test", r);
      document.getElementById("status").textContent = testNote = r.error ? "Raw self-test: " + r.error :
        `Raw self-test ${r.ok ? "passed" : "FAILED"} (${r.base}): clock and motion of ` +
        `${r.probes.filter((q) => q.ok).length} of ${r.probes.length} probes load; ` +
        `${r.decode.filter((q) => q.ok).length} of ${r.decode.length} ` +
        `chunks decode to the stored sum; processing of ${r.process.length} tiles differs from Python by at most ` +
        `${Math.max(...r.process.map((q) => q.maxAbsDiff_uV)).toExponential(1)} µV ` +
        `(${Math.max(...r.process.map((q) => q.ms)).toFixed(0)} ms per tile); drawing of ${r.draw.length} views: ` +
        `${Math.max(...r.draw.map((q) => q.error ? 1 : q.fracDiffer)).toExponential(1)} of the pixels differ at most; ` +
        (!r.display ? "no display copy; " : r.displayError ? `display copy: ${r.displayError}; ` :
          `display copy (${r.display}): ${r.ddecode.filter((q) => q.ok).length} of ${r.ddecode.length} tiles decode to the stored sum; ` +
          `level 0 differs from the recorded samples by ${Math.max(...r.dcross.map((q) => q.rms_uV)).toFixed(1)} µV rms at most; ` +
          `${r.dcross.reduce((a, q) => a + q.level1Wrong, 0)} values of level 1 disagree with level 0; drawing of ${r.ddraw.length} views: ` +
          `${Math.max(...r.ddraw.map((q) => q.error ? 1 : q.fracDiffer)).toExponential(1)} of the pixels differ at most; `) +
        `details in window.__rawTest`;
    } else if (new URLSearchParams(location.search).has("sonifytest")) {
      document.getElementById("status").textContent = "Running the sonify self-test";
      const r = await Sonify.test();
      const bad = r.checks.filter((c) => !c.pass).map((c) => c.name);
      document.getElementById("status").textContent = testNote = `Sonify self-test: ${r.passed} of ${r.total} checks passed` +
        (bad.length ? `; FAILED: ${bad.join("; ")}` : "") + "; details in console / window.__sonifyTest";
      requestRender();
    } else if (new URLSearchParams(location.search).has("test")) {
      document.getElementById("status").textContent = "Running self-test";
      const r = await runTests();
      const worst = Math.max(...r.report.map((x) => x.skipped ? 0 : x.fracBad));
      const nDiff = r.report.reduce((a, x) => a + (x.countsDiffer || 0), 0);
      const cpuWorst = Math.max(...r.report.map((x) => x.skipped ? 0 : x.cpuMaxRel));
      const ms = Object.entries(r.fullViewMs)
        .map(([m, t]) => `${m} ${Math.min(...t).toFixed(1)} ms`).join(", ");
      document.getElementById("status").textContent = testNote =
        `Self-test: ${r.report.length} view x method runs, worst fraction of mismatched cells ` +
        `${worst.toExponential(2)}; counts differing between methods: ${nDiff}; ` +
        `readout (CPU) max relative error ${cpuWorst.toExponential(1)}; ` +
        `burst detection ${r.burstCheck.match ? "matches" : "DIFFERS FROM"} Python; full view ${r.grid[0]}x${r.grid[1]}: ${ms}; ` +
        `details in console / window.__testResults`;
      requestRender();
    }
  } catch (err) {
    document.getElementById("loadMsg").textContent = "Error: " + err.message;
    console.error(err);
  }
})();
