// Behavioral traces under the raster (optional).
//
// If data/traces.json exists (written by traces/prep_traces.py), each trace it lists is drawn in
// its own lane below the raster, on the raster's time axis. A trace is a regularly sampled
// series: sample i belongs to time t0 + i * dt on the viewer's clock; NaN marks missing data.
//
// How a lane is drawn:
//   - black background, a gray dashed line at zero, the trace as a white line;
//   - the area between zero and the trace is shaded, green above zero and red below, and the
//     shade deepens with distance from zero. The vertical range of a lane never changes with
//     the view, so in a window too short to show any change over time the height and shade of
//     the band still show the level;
//   - when a screen column holds two or more samples, the white line is their mean, and a
//     translucent band spans their smallest to their largest value, so brief events stay
//     visible at any zoom. The mean, minimum and maximum of a column are exact (no sample is
//     skipped or counted twice): see aggregate().
//
// Right-click a lane to hide it; traces that are not shown are listed at the bottom left, click
// one to show it. A trace listed with "show": false in traces.json starts out not shown.
// Drag the bar between the raster and the traces to divide the height between them.
// ?tracetest=1 runs a self-test.
"use strict";

const Traces = (() => {
  const GAP = 12;            // CSS px between the raster and the trace panel (holds the drag bar)
  const LANE_GAP = 6;        // CSS px between lanes
  const MIN_LANE = 10;       // CSS px: the least height of a lane
  const F = 4;               // each pyramid level summarizes F entries of the level below
  const FRAC_DEFAULT = 0.25, FRAC_MIN = 0.04, FRAC_MAX = 0.8;
  const POS = "40,210,90", NEG = "245,60,50", SHADE_MAX = 0.95;
  const FILE = "traces.json", FORMAT = "rasterViewer-traces-1";

  let list = [];             // the traces, in the order of the file
  let frac = FRAC_DEFAULT;   // share of the plot height given to the trace panel
  let cv = null, g = null, handle = null, chips = null;
  let hoverX = null;         // pointer x over the panel (CSS px in the stage), or null
  let loadError = null;

  // ------------------------------------------------------------ data
  // levels[0] is the samples themselves; entry j of level k summarizes samples
  // [j * F^k, (j + 1) * F^k): smallest, largest, sum and number of the non-missing ones.
  function buildPyramid(v) {
    const levels = [{ n: v.length }];
    let n = v.length, prev = null;
    while (n > F) {
      const m = Math.ceil(n / F);
      const mn = new Float32Array(m).fill(Infinity), mx = new Float32Array(m).fill(-Infinity);
      const sum = new Float64Array(m), cnt = new Uint32Array(m);
      for (let i = 0; i < n; i++) {
        const j = Math.floor(i / F);
        if (!prev) {
          const x = v[i];
          if (x === x) { if (x < mn[j]) mn[j] = x; if (x > mx[j]) mx[j] = x; sum[j] += x; cnt[j]++; }
        } else if (prev.cnt[i]) {
          if (prev.mn[i] < mn[j]) mn[j] = prev.mn[i];
          if (prev.mx[i] > mx[j]) mx[j] = prev.mx[i];
          sum[j] += prev.sum[i]; cnt[j] += prev.cnt[i];
        }
      }
      prev = { n: m, mn, mx, sum, cnt };
      levels.push(prev);
      n = m;
    }
    return levels;
  }

  // Smallest, largest, sum and count of the non-missing samples with index in [ia, ib).
  // Exact: the range is split into whole pyramid entries, at most 2 (F - 1) per level.
  function aggregate(tr, ia, ib, o = {}) {
    const v = tr.v, lv = tr.levels, top = lv.length - 1;
    let lo = Math.max(0, ia), hi = Math.min(tr.n, ib), k = 0;
    let mn = Infinity, mx = -Infinity, sum = 0, cnt = 0;
    const take = (i) => {
      if (k === 0) {
        const x = v[i];
        if (x === x) { if (x < mn) mn = x; if (x > mx) mx = x; sum += x; cnt++; }
      } else {
        const L = lv[k], c = L.cnt[i];
        if (c) { if (L.mn[i] < mn) mn = L.mn[i]; if (L.mx[i] > mx) mx = L.mx[i]; sum += L.sum[i]; cnt += c; }
      }
    };
    while (lo < hi) {
      if (k === top) { while (lo < hi) take(lo++); break; }
      while (lo < hi && lo % F) take(lo++);
      while (lo < hi && hi % F) take(--hi);
      lo /= F; hi /= F; k++;
    }
    o.mn = mn; o.mx = mx; o.sum = sum; o.cnt = cnt;
    return o;
  }

  function makeTrace(d, v) {
    if (!(d.dt > 0) || !Number.isFinite(d.t0)) throw new Error(`trace ${d.name}: bad t0 or dt`);
    if (v.length !== d.n) throw new Error(`trace ${d.name}: ${v.length} samples, the list says ${d.n}`);
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < v.length; i++) { const x = v[i]; if (x === x) { if (x < lo) lo = x; if (x > hi) hi = x; } }
    let ylim = Array.isArray(d.ylim) && d.ylim[1] > d.ylim[0] ? d.ylim.slice(0, 2) : null;
    if (!ylim) ylim = lo < hi ? [Math.min(lo, 0), Math.max(hi, 0)] : [-1, 1];
    return { name: d.name, label: d.label || d.name, unit: d.unit || "", t0: d.t0, dt: d.dt, n: d.n,
             ylim, v, levels: buildPyramid(v), hidden: d.show === false, rect: null, info: d };
  }

  async function load() {
    try {
      const resp = await fetch(DATA + FILE);
      if (!resp.ok) return;                                   // no traces for this data set
      const man = await resp.json();
      if (man.format !== FORMAT) throw new Error("unknown format " + man.format);
      const out = [];
      for (const d of man.traces) {
        const r = await fetch(DATA + d.file);
        if (!r.ok) throw new Error(`${d.file}: HTTP ${r.status}`);
        out.push(makeTrace(d, new Float32Array(await r.arrayBuffer())));
      }
      list = out;
    } catch (err) {
      loadError = err.message; list = [];
      console.error("behavioral traces not loaded:", err);
    }
  }

  // ------------------------------------------------------------ drawing one lane
  // ctx: 2D context in CSS px. R: {x, y, w, h}, the lane. vw: {t0, span, shift, dpr}: time t is
  // drawn at x = R.x + (t - t0) / span * R.w - shift (shift: see Sonify.shiftCss).
  // Returns what it drew, for the readout and the self-test.
  function lanePath(tr, R, vw) {
    const tOfX = (x) => vw.t0 + (x - R.x + vw.shift) / R.w * vw.span;
    const xOfT = (t) => R.x + (t - vw.t0) / vw.span * R.w - vw.shift;
    const nCol = Math.max(1, Math.round(R.w * vw.dpr));
    const perCol = vw.span / nCol / tr.dt;                    // samples per device column
    const segs = [];                                          // runs of consecutive points with data
    let cur = null;
    const push = (x, m, lo, hi) => {
      if (!cur) { cur = { x: [], m: [], lo: [], hi: [] }; segs.push(cur); }
      cur.x.push(x); cur.m.push(m); cur.lo.push(lo); cur.hi.push(hi);
    };
    if (perCol < 2) {                                         // every sample is a point of the line
      const ia = Math.max(0, Math.floor((tOfX(R.x) - tr.t0) / tr.dt) - 1);
      const ib = Math.min(tr.n, Math.ceil((tOfX(R.x + R.w) - tr.t0) / tr.dt) + 2);
      for (let i = ia; i < ib; i++) {
        const x = tr.v[i];
        if (x === x) push(xOfT(tr.t0 + i * tr.dt), x, x, x); else cur = null;
      }
      return { segs, aggregated: false, perCol, nCol };
    }
    // column c holds the samples with time in [time of its left edge, time of its right edge)
    const o = {}, colMean = new Float64Array(nCol).fill(NaN), colCnt = new Uint32Array(nCol), colSum = new Float64Array(nCol);
    let ia = Math.ceil((tOfX(R.x) - tr.t0) / tr.dt);
    for (let c = 0; c < nCol; c++) {
      const ib = Math.ceil((tOfX(R.x + (c + 1) / vw.dpr) - tr.t0) / tr.dt);
      aggregate(tr, ia, ib, o);
      colCnt[c] = o.cnt; colSum[c] = o.sum;
      if (o.cnt) { colMean[c] = o.sum / o.cnt; push(R.x + (c + 0.5) / vw.dpr, colMean[c], o.mn, o.mx); } else cur = null;
      ia = ib;
    }
    return { segs, aggregated: true, perCol, nCol, colMean, colCnt, colSum };
  }

  function shadeGradient(ctx, R, lo, hi) {
    const vmax = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
    // an end of the range at zero takes the color of the other end (transparent, but canvas
    // gradients mix the colors themselves, so a green stop would tint a red lane)
    const col = (v, rgb) => `rgba(${rgb},${(SHADE_MAX * Math.min(1, Math.abs(v) / vmax)).toFixed(4)})`;
    const grd = ctx.createLinearGradient(0, R.y, 0, R.y + R.h);
    grd.addColorStop(0, col(hi, hi > 0 ? POS : NEG));
    if (lo < 0 && hi > 0) {
      const z = hi / (hi - lo);
      grd.addColorStop(z, `rgba(${POS},0)`); grd.addColorStop(z, `rgba(${NEG},0)`);
    }
    grd.addColorStop(1, col(lo, lo < 0 ? NEG : POS));
    return grd;
  }

  function drawLane(ctx, tr, R, vw) {
    const [lo, hi] = tr.ylim;
    const yOf = (v) => R.y + (hi - v) / (hi - lo) * R.h;
    const y0 = Math.min(R.y + R.h, Math.max(R.y, yOf(0)));    // zero, or the nearest edge of the lane
    const P = lanePath(tr, R, vw);
    ctx.save();
    ctx.beginPath(); ctx.rect(R.x, R.y, R.w, R.h); ctx.clip();
    ctx.fillStyle = "#000000"; ctx.fillRect(R.x, R.y, R.w, R.h);
    // a lone point (data on one side only) is drawn half a sample wide
    const half = Math.max(0.5 / vw.dpr, 0.5 * tr.dt / vw.span * R.w);
    const xs = (s) => s.x.length > 1 ? s.x : [s.x[0] - half, s.x[0] + half];
    const ys = (s, a) => s.x.length > 1 ? a : [a[0], a[0]];
    ctx.fillStyle = shadeGradient(ctx, R, lo, hi);
    for (const s of P.segs) {
      const x = xs(s), m = ys(s, s.m);
      ctx.beginPath(); ctx.moveTo(x[0], y0);
      for (let i = 0; i < x.length; i++) ctx.lineTo(x[i], yOf(m[i]));
      ctx.lineTo(x[x.length - 1], y0); ctx.closePath(); ctx.fill();
    }
    if (P.aggregated) {                                       // smallest to largest value of each column
      ctx.fillStyle = "rgba(255,255,255,0.30)";
      for (const s of P.segs) {
        const x = xs(s), a = ys(s, s.hi), b = ys(s, s.lo);
        ctx.beginPath(); ctx.moveTo(x[0], yOf(a[0]));
        for (let i = 1; i < x.length; i++) ctx.lineTo(x[i], yOf(a[i]));
        for (let i = x.length - 1; i >= 0; i--) ctx.lineTo(x[i], yOf(b[i]));
        ctx.closePath(); ctx.fill();
      }
    }
    if (lo <= 0 && hi >= 0) {
      ctx.strokeStyle = "#8a8f96"; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
      const yz = Math.round(yOf(0) * vw.dpr) / vw.dpr + 0.5 / vw.dpr;
      ctx.beginPath(); ctx.moveTo(R.x, yz); ctx.lineTo(R.x + R.w, yz); ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 1.25; ctx.lineJoin = "round";
    for (const s of P.segs) {
      const x = xs(s), m = ys(s, s.m);
      ctx.beginPath(); ctx.moveTo(x[0], yOf(m[0]));
      for (let i = 1; i < x.length; i++) ctx.lineTo(x[i], yOf(m[i]));
      ctx.stroke();
    }
    ctx.restore();
    return P;
  }

  // The value shown for the pointer at x (CSS px), given what drawLane drew (P): the mean of
  // the device column under the pointer when the line is a mean, otherwise the nearest sample.
  // Returns {v, x}: the value (null where there is no data) and where to draw the cursor line.
  function valueAt(tr, P, x, R, vw) {
    if (P.aggregated) {
      const c = Math.min(P.nCol - 1, Math.max(0, Math.floor((x - R.x) * vw.dpr)));
      return { v: P.colCnt[c] ? P.colMean[c] : null, x: R.x + (c + 0.5) / vw.dpr };
    }
    const t = vw.t0 + (x - R.x + vw.shift) / R.w * vw.span;
    const i = Math.round((t - tr.t0) / tr.dt);
    return { v: i >= 0 && i < tr.n && tr.v[i] === tr.v[i] ? tr.v[i] : null, x };
  }

  // ------------------------------------------------------------ the panel
  const visible = () => list.filter((t) => !t.hidden);
  const fmtVal = (x) => (Math.abs(x) >= 1000 ? x.toFixed(0) : Math.abs(x) < 1e-3 ? (x === 0 ? "0" : x.toExponential(1))
    : String(+x.toPrecision(3))).replace(/-/g, "−");

  // Height taken from the plot area (CSS px): {gap, h}, both 0 when no trace is shown.
  function panel(plotH) {
    const n = visible().length;
    if (!n) return { gap: 0, h: 0 };
    const hMin = n * MIN_LANE + (n - 1) * LANE_GAP;          // every lane keeps some height
    const h = Math.round(Math.min(FRAC_MAX, Math.max(FRAC_MIN, frac)) * plotH);
    return { gap: GAP, h: Math.max(hMin, Math.min(h, plotH - GAP - 50)) };
  }

  function setup() {
    if (cv || (!list.length && !loadError)) return;
    const stage = document.getElementById("stage");
    cv = document.createElement("canvas");
    Object.assign(cv.style, { position: "absolute", left: "0", cursor: "crosshair", display: "none" });
    handle = document.createElement("div");
    handle.title = "Drag to divide the height between the raster and the traces (double-click: default)";
    Object.assign(handle.style, { position: "absolute", height: GAP + "px", cursor: "ns-resize", display: "none" });
    handle.innerHTML = `<div style="position:absolute;left:50%;top:${GAP / 2 - 2}px;width:56px;height:4px;` +
      `margin-left:-28px;border-radius:2px;background:#6c747d"></div>`;
    chips = document.createElement("div");
    Object.assign(chips.style, { position: "absolute", left: "8px", bottom: "3px", font: "11px Arial",
      color: "#6c747d", whiteSpace: "nowrap", display: "none" });
    stage.insertBefore(cv, document.getElementById("box"));
    stage.appendChild(handle); stage.appendChild(chips);
    g = cv.getContext("2d");

    const grip = handle.firstChild;
    handle.addEventListener("mouseenter", () => { grip.style.background = "#4ea8ff"; });
    handle.addEventListener("mouseleave", () => { if (!dragging) grip.style.background = "#6c747d"; });
    let dragging = false;
    handle.addEventListener("mousedown", (e) => {
      if (e.button !== 0) return;
      dragging = true; e.preventDefault(); e.stopPropagation();
    });
    handle.addEventListener("dblclick", (e) => { e.preventDefault(); frac = FRAC_DEFAULT; relayout(); });
    window.addEventListener("mousemove", (e) => {
      if (!dragging) return;
      const r = stage.getBoundingClientRect(), plotH = Math.max(50, layout.sh - TOP - BOT);
      const y = e.clientY - r.top + GAP / 2;                  // where the panel's top edge should be
      frac = Math.min(FRAC_MAX, Math.max(FRAC_MIN, (TOP + plotH - y) / plotH));
      relayout();
    });
    window.addEventListener("mouseup", () => { if (dragging) { dragging = false; grip.style.background = "#6c747d"; } });

    cv.addEventListener("mousemove", (e) => {
      hoverX = e.clientX - stage.getBoundingClientRect().left; draw();
    });
    cv.addEventListener("mouseleave", () => { hoverX = null; draw(); });
    chips.addEventListener("mousedown", (e) => e.stopPropagation());   // not a click on the plot
    chips.addEventListener("click", (e) => {
      const name = e.target && e.target.dataset ? e.target.dataset.trace : null;
      const tr = list.find((t) => t.name === name);
      if (tr) { tr.hidden = false; relayout(); }
    });
  }

  function relayout() { computeLayout(); requestRender(); }

  // Called at the end of the viewer's computeLayout: put the panel where the layout says.
  function place() {
    if (!cv) return;
    const { dpr, sw, x: px, w: pw } = layout, H = layout.traceH || 0;
    cv.style.display = handle.style.display = H ? "block" : "none";
    if (H) {
      cv.width = Math.round(sw * dpr); cv.height = Math.round(H * dpr);
      Object.assign(cv.style, { top: layout.traceY + "px", width: sw + "px", height: H + "px" });
      Object.assign(handle.style, { left: px + "px", width: pw + "px", top: (layout.traceY - GAP) + "px" });
    }
    const hid = list.filter((t) => t.hidden);
    chips.style.display = list.length || loadError ? "block" : "none";
    chips.textContent = loadError ? `Behavioral traces could not be loaded: ${loadError}`
      : !hid.length ? "Right-click a trace to hide it" : "Behavioral traces not shown (click one to show it): ";
    if (!loadError) for (const t of hid) {
      const b = document.createElement("span");
      b.dataset.trace = t.name; b.title = "Show this trace"; b.textContent = t.label;
      b.style.cssText = "cursor:pointer;color:#a9b0b8;border:1px solid #3a424b;border-radius:3px;padding:1px 5px;margin-right:4px";
      chips.appendChild(b);
    }
  }

  function lanes() {
    const vis = visible(), n = vis.length, H = layout.traceH || 0;
    const lh = (H - LANE_GAP * (n - 1)) / n;
    return vis.map((tr, i) => ({ tr, R: { x: layout.x, y: i * (lh + LANE_GAP), w: layout.w, h: lh } }));
  }

  function draw() {
    if (!cv || !layout.traceH) return;
    const { dpr, sw, x: px, w: pw } = layout, H = layout.traceH;
    const ink = "#e6e9ed", ink2 = "#a9b0b8", ink3 = "#6c747d";
    const vw = { t0: state.t0, span: state.t1 - state.t0, shift: Sonify.shiftCss, dpr };
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, sw, H);
    const hx = hoverX !== null && hoverX >= px && hoverX < px + pw ? hoverX : null;
    let lineX = hx;
    for (const { tr, R } of lanes()) {
      tr.rect = R;
      const P = drawLane(g, tr, R, vw);
      const [lo, hi] = tr.ylim, yOf = (v) => R.y + (hi - v) / (hi - lo) * R.h;
      // name and unit, and the values at the top, at zero and at the bottom of the lane
      g.textAlign = "right"; g.textBaseline = "middle";
      const two = R.h >= 30 && tr.unit;
      g.font = "12px Arial"; g.fillStyle = ink;
      g.fillText(tr.label, px - 40, R.y + R.h / 2 - (two ? 7 : 0));
      if (two) { g.font = "11px Arial"; g.fillStyle = ink2; g.fillText(`(${tr.unit})`, px - 40, R.y + R.h / 2 + 8); }
      g.font = "10px Arial"; g.fillStyle = ink3; g.strokeStyle = ink3; g.lineWidth = 1;
      const tick = (v, base) => {
        const y = yOf(v);
        g.beginPath(); g.moveTo(px - 4, y); g.lineTo(px, y); g.stroke();
        g.textBaseline = base; g.fillText(fmtVal(v), px - 7, y);
      };
      if (R.h >= 24) { tick(hi, "top"); tick(lo, "bottom"); }
      if (lo < 0 && hi > 0 && yOf(0) - R.y >= 14 && R.y + R.h - yOf(0) >= 14) tick(0, "middle");
      if (hx !== null) {
        const q = valueAt(tr, P, hx, R, vw);
        if (P.aggregated) lineX = q.x;                       // the middle of the column that is read
        g.textAlign = "left"; g.textBaseline = "middle"; g.font = "12px Arial"; g.fillStyle = ink;
        g.fillText(q.v === null ? "no data" : `${fmtVal(q.v)} ${tr.unit}`, px + pw + 8, R.y + R.h / 2);
      }
    }
    const vline = (x, style) => {
      g.strokeStyle = style; g.lineWidth = 1;
      g.beginPath(); g.moveTo(x, 0); g.lineTo(x, H); g.stroke();
    };
    if (hx !== null) vline(lineX, "rgba(255,255,255,0.45)");
    if (Sonify.headX !== null) vline(Sonify.headX, "rgba(255,255,255,0.9)");
  }

  // ------------------------------------------------------------ pointer (called by the viewer)
  // m: pointer position in the stage (CSS px). True when it is over the trace panel.
  const hit = (m) => !!layout.traceH && m.y >= layout.traceY - GAP && m.y <= layout.traceY + layout.traceH;
  function hideAt(m) {
    const y = m.y - layout.traceY;
    const L = lanes().find(({ R }) => y >= R.y - LANE_GAP / 2 && y <= R.y + R.h + LANE_GAP / 2);
    if (!L || m.x < layout.x - 140 || m.x > layout.x + layout.w) return false;
    L.tr.hidden = true; hoverX = null; relayout();
    return true;
  }

  // ------------------------------------------------------------ URL
  // tshow lists the traces that are shown (empty: none), so a link shows what its maker saw
  // whatever the data set's defaults are; without tshow the defaults of traces.json apply.
  function toHash(p) {
    if (!list.length) return;
    p.set("tfrac", frac.toFixed(3));
    p.set("tshow", visible().map((t) => t.name).join(","));
  }
  function fromHash(p) {
    const f = parseFloat(p.get("tfrac"));
    if (Number.isFinite(f)) frac = Math.min(FRAC_MAX, Math.max(FRAC_MIN, f));
    const shown = p.has("tshow") ? new Set(p.get("tshow").split(",")) : null;
    for (const t of list) t.hidden = shown ? !shown.has(t.name) : t.info.show === false;
  }

  // ------------------------------------------------------------ self-test (?tracetest=1)
  function test() {
    const checks = [];
    const check = (name, pass, detail) => checks.push({ name, pass: !!pass, detail });
    let seed = 12345;
    const rnd = () => (seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296;
    const brute = (v, ia, ib) => {
      let mn = Infinity, mx = -Infinity, sum = 0, cnt = 0;
      for (let i = Math.max(0, ia); i < Math.min(v.length, ib); i++) {
        const x = v[i]; if (x === x) { if (x < mn) mn = x; if (x > mx) mx = x; sum += x; cnt++; }
      }
      return { mn, mx, sum, cnt };
    };
    const same = (a, b) => a.mn === b.mn && a.mx === b.mx && a.cnt === b.cnt &&
      Math.abs(a.sum - b.sum) <= 1e-9 * (1 + Math.abs(b.sum));

    // 1. the files: counts and sums recorded by prep_traces.py
    for (const tr of list) {
      const b = brute(tr.v, 0, tr.n), d = tr.info;
      check(`${tr.name}: file matches the list (${tr.n} samples, ${b.cnt} with data)`,
            b.cnt === d.n_valid && Math.abs(b.sum - d.sum) <= 1e-6 * (1 + Math.abs(d.sum)),
            { n_valid: [b.cnt, d.n_valid], sum: [b.sum, d.sum] });
    }
    // 2. aggregate() against a plain loop: synthetic series with gaps, and every loaded trace
    const synth = [1, 3, 4, 5, 16, 17, 1000, 4097].map((n) => {
      const v = new Float32Array(n);
      for (let i = 0; i < n; i++) v[i] = rnd() < 0.15 ? NaN : (rnd() - 0.5) * 20;
      if (n > 600) for (let i = 300; i < 600; i++) v[i] = NaN;        // a long gap
      return { name: `synthetic series of ${n}`, v, n, levels: buildPyramid(v) };
    });
    for (const tr of [...synth, ...list]) {
      let bad = 0, tried = 0;
      for (let q = 0; q < 2000; q++) {
        const a = Math.floor(rnd() * (tr.n + 20)) - 10;
        const len = q % 3 === 0 ? Math.floor(rnd() * 12) : Math.floor(rnd() * (tr.n + 10));
        tried++; if (!same(aggregate(tr, a, a + len), brute(tr.v, a, a + len))) bad++;
      }
      if (!same(aggregate(tr, 0, tr.n), brute(tr.v, 0, tr.n))) bad++;
      check(`${tr.name}: range summaries equal a plain loop (${tried + 1} ranges)`, bad === 0, { bad });
    }
    // 3. columns partition the samples: every sample in view lands in exactly one column
    {
      const n = 5000, v = new Float32Array(n);
      for (let i = 0; i < n; i++) v[i] = Math.sin(i / 37) * 3 + (i % 97 === 0 ? 8 : 0);
      const tr = { name: "t", v, n, t0: 10.003, dt: 0.01, ylim: [-4, 8], levels: buildPyramid(v) };
      const R = { x: 20, y: 0, w: 333, h: 60 };
      let ok = true; const detail = [];
      for (const [t0, span, shift, dpr] of [[10, 50, 0, 1], [12.5, 31.7, 0.5, 2], [5, 70, 0, 1.25], [20.2, 9, 0.8, 1]]) {
        const vw = { t0, span, shift, dpr }, P = lanePath(tr, R, vw);
        const tl = t0 + shift / R.w * span, te = tl + span;
        const want = brute(v, Math.ceil((tl - tr.t0) / tr.dt), Math.ceil((te - tr.t0) / tr.dt));
        let mn = Infinity, mx = -Infinity, pts = 0, cnt = 0, sum = 0;
        for (const s of P.segs) for (let i = 0; i < s.x.length; i++) { mn = Math.min(mn, s.lo[i]); mx = Math.max(mx, s.hi[i]); pts++; }
        for (let c = 0; c < P.nCol; c++) { cnt += P.colCnt[c]; sum += P.colSum[c]; }
        // the readout at the middle of every column is that column's mean
        let readBad = 0;
        for (let c = 0; c < P.nCol; c++) {
          const q = valueAt(tr, P, R.x + (c + 0.5) / dpr, R, vw);
          if (!(q.v === (P.colCnt[c] ? P.colMean[c] : null))) readBad++;
        }
        const good = P.aggregated && mn === want.mn && mx === want.mx && cnt === want.cnt &&
                     Math.abs(sum - want.sum) < 1e-9 * (1 + Math.abs(want.sum)) && readBad === 0;
        ok = ok && good; detail.push({ t0, span, good, pts, cnt, wantCnt: want.cnt, sum, wantSum: want.sum, readBad });
      }
      check("screen columns hold every sample in view exactly once (count, sum, smallest, largest), and the readout gives each column's mean", ok, detail);
    }
    // 4. drawing: a step from +2 to -1 in a lane from -2 to +4
    {
      const n = 400, v = new Float32Array(n);
      for (let i = 0; i < n; i++) v[i] = i < 200 ? 2 : i < 300 ? -1 : NaN;
      const tr = { name: "step", v, n, t0: 0, dt: 0.01, ylim: [-2, 4], levels: buildPyramid(v) };
      const c = document.createElement("canvas"); c.width = 400; c.height = 122;
      const x = c.getContext("2d", { willReadFrequently: true });
      {   // a lane wholly below zero is shaded red only, one wholly above zero green only
        const lane = { x: 0, y: 0, w: 4, h: 100 }, v1 = new Float32Array(4).fill(-2), v2 = new Float32Array(4).fill(2);
        const mk = (v, ylim) => ({ name: "g", v, n: 4, t0: 0, dt: 1, ylim, levels: buildPyramid(v) });
        drawLane(x, mk(v1, [-3, 0]), lane, { t0: -0.5, span: 4, shift: 0, dpr: 1 });
        const pn = [20, 40, 55].map((yy) => Array.from(x.getImageData(2, yy, 1, 1).data));
        drawLane(x, mk(v2, [0, 3]), lane, { t0: -0.5, span: 4, shift: 0, dpr: 1 });
        const pp = [45, 60, 80].map((yy) => Array.from(x.getImageData(2, yy, 1, 1).data));
        check("drawing: a lane that does not cross zero is shaded in one color",
              pn.every((q) => q[1] <= q[2] + 6 && q[0] > q[1]) && pp.every((q) => q[0] <= q[2] + 6 && q[1] > q[0]), { pn, pp });
        x.clearRect(0, 0, 400, 122);
      }
      const R = { x: 0, y: 0.5, w: 400, h: 120 };            // values at half-integer y: lines fill a pixel row
      drawLane(x, tr, R, { t0: 0, span: 4, shift: 0, dpr: 1 });   // one sample per pixel
      const px = (xx, yy) => Array.from(x.getImageData(xx, yy, 1, 1).data);
      const yOf = (val) => Math.floor(0.5 + (4 - val) / 6 * 120);
      const white = (p) => p[0] > 200 && p[1] > 200 && p[2] > 200;
      const green = (p) => p[1] > p[0] + 15 && p[1] > p[2] + 15;
      const red = (p) => p[0] > p[1] + 25 && p[0] > p[2] + 25;
      const black = (p) => p[0] < 12 && p[1] < 12 && p[2] < 12;
      const a = px(100, yOf(2)), b = px(100, yOf(1.8)), b2 = px(100, yOf(1)), c0 = px(100, yOf(3)), d = px(100, yOf(-1));
      const e = px(250, yOf(-1)), f = px(250, yOf(-0.9)), h = px(250, yOf(1)), k = px(350, yOf(1)), k2 = px(350, yOf(-1));
      check("drawing: the line is white at the trace's value", white(a) && white(e), { a, e });
      check("drawing: green between zero and a positive value, deeper further from zero",
            green(b) && green(b2) && b[1] > b2[1] + 20, { b, b2 });
      check("drawing: red between zero and a negative value", red(f), { f });
      check("drawing: nothing beyond the trace, on the other side of zero, or where data are missing",
            black(c0) && black(d) && black(h) && black(k) && black(k2), { c0, d, h, k, k2 });
      // the same series seen from far away, 8 samples per column: column c holds samples
      // 8c - 3 to 8c + 4, so column 25 holds three samples of +2 and five of -1
      x.clearRect(0, 0, 400, 122);
      const P = drawLane(x, tr, { x: 0, y: 0.5, w: 50, h: 120 }, { t0: -0.035, span: 4, shift: 0, dpr: 1 });
      const s0 = P.segs[0], col = (i) => [s0.m[i], s0.lo[i], s0.hi[i]];
      const eq = (u, w) => u.every((q, i) => Math.abs(q - w[i]) < 1e-12);
      check("far view: a column that spans the step shows its mean, smallest and largest value",
            P.aggregated && P.segs.length === 1 && s0.x.length === 38 && eq(col(0), [2, 2, 2]) && eq(col(24), [2, 2, 2]) &&
            eq(col(25), [0.125, -1, 2]) && eq(col(26), [-1, -1, -1]),
            { n: s0.x.length, c0: col(0), c24: col(24), c25: col(25), c26: col(26) });
    }
    // 5. which traces are shown: the defaults of the file, and the URL
    if (list.length) {
      const was = list.map((t) => t.hidden), wasFrac = frac;
      const shownNow = () => visible().map((t) => t.name).join(",");
      const names = list.map((t) => t.name), dflt = list.filter((t) => t.info.show !== false).map((t) => t.name).join(",");
      fromHash(new URLSearchParams(""));
      const a = shownNow();
      fromHash(new URLSearchParams("tshow=" + names[names.length - 1]));
      const b = shownNow();
      fromHash(new URLSearchParams("tshow="));
      const c = shownNow();
      fromHash(new URLSearchParams("tshow=" + names.join(",") + "&tfrac=0.4"));
      const d = shownNow(), p = new URLSearchParams(); toHash(p);
      check("shown traces: the file's defaults without tshow in the URL; otherwise exactly those tshow lists; the URL written back says the same",
            a === dflt && b === names[names.length - 1] && c === "" && d === names.join(",") &&
            p.get("tshow") === d && p.get("tfrac") === "0.400", { dflt, a, b, c, d, written: p.toString() });
      // 6. the panel on screen: lanes inside the panel, time lined up with the raster's axis
      computeLayout();                                       // all traces shown
      const L = lanes(), last = L[L.length - 1].R;
      check("panel: lanes fill the panel and share the raster's left and right edges",
            L[0].R.y === 0 && Math.abs(last.y + last.h - layout.traceH) < 1e-6 && L.every(({ R }) => R.x === layout.x && R.w === layout.w),
            { traceH: layout.traceH, n: L.length });
      check("panel: sits below the raster, above the time axis",
            layout.traceY >= layout.y + layout.h && Math.abs(layout.traceY + layout.traceH - layout.axisY) < 1e-6,
            { rasterBottom: layout.y + layout.h, traceY: layout.traceY, axisY: layout.axisY });
      const hAll = layout.h;
      for (const t of list) t.hidden = true;
      computeLayout();
      check("panel: with no trace shown it takes no height (the raster reaches the time axis)",
            layout.traceH === 0 && layout.axisY === layout.y + layout.h && layout.h > hAll && cv.style.display === "none",
            { traceH: layout.traceH, rasterH: [hAll, layout.h] });
      list.forEach((t, i) => { t.hidden = was[i]; }); frac = wasFrac;
      relayout();
    }
    const passed = checks.filter((c) => c.pass).length;
    return { passed, total: checks.length, checks, traces: list.map((t) => t.name) };
  }

  return { load, setup, panel, place, draw, hit, hideAt, toHash, fromHash, test,
           get count() { return list.length; },
           get debug() { return { list, frac, lanes: lanes(), aggregate, lanePath, drawLane, valueAt }; } };
})();
