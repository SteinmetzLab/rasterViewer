"use strict";
// Raw-voltage underlay (prototype). When the neurons are plotted by depth and the window is
// short, the recorded voltage of each site is drawn under the spikes, on its own canvas below
// the main one. The data are optional: rawdata/manifest.json, written by
// raw_underlay/make_raw_proto.py and served by `python serve.py --raw <folder>`. Without it
// (e.g. on the public site) none of this is active.
//
// Drawing. A tile is 0.2 s of one shank (384 sites), processed in a worker (raw_worker.js:
// high-pass, median reference) and kept as a texture in microvolts. Every device pixel shows
// the sample of largest magnitude among the samples and sites it covers.
//   Time: the viewer's clock -> probe sample through the pairs in clock_p<probe>.f64 (linear
//         between pairs), the same map the spike times went through. Sample s is drawn over
//         [s - 0.5, s + 0.5).
//   Depth: a site at height y occupies the same 15 um band as a neuron at that height. Its two
//         sites (x = 27 and 59 um) share the band, upper and lower half.
// Three stored versions of the data can be compared (V), and so can the difference between a
// lossy version and the lossless one (D).

const RawUnderlay = (() => {
  const MAX_SPAN_S = 0.5;        // longest window with the underlay
  const MAX_TEX = 40;            // tiles kept on the GPU (6000 x 384 half floats, 4.6 MB each)
  const RETRY_MS = 3000;         // a tile that failed to load is asked for again after this
  const ALPHAS = [1, 0.6, 0];    // opacity of the spike cells over the voltage (key O)
  let man = null, clock = null, cv = null, g = null, P = null, boundsTex = null, vao = null;
  let covT = [0, 0];             // viewer time covered by the raw data
  const R = { on: true, variant: 0, diff: false, uv: 100, alpha: 0 };
  const tiles = new Map();       // "variant|shank|tile" -> { state, tex, used }
  let workers = [], jobs = new Map(), queue = [], nextId = 1, frame = 0;
  let info = { active: false, hint: "", loading: 0, failed: 0 };

  const VS = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;
  const FRAG = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D u_a, u_b;     // R16F, tile x 384: microvolts (u_b: the version to subtract)
uniform sampler2D u_bounds;     // R32F, (PW + 1) x 1: sample at each pixel-column edge, minus base
uniform float u_tileOff;        // first sample of this tile, minus base
uniform float u_rho0, u_drho;   // half-site rows from the top of the shank: at pixel row 0, per row
uniform float u_range;          // microvolts at black and white
uniform int u_PH, u_n, u_diff;
out vec4 o;
void main() {
  int x = int(gl_FragCoord.x);
  float p = float(u_PH - 1 - int(gl_FragCoord.y));
  float sa = texelFetch(u_bounds, ivec2(x, 0), 0).r - u_tileOff;
  float sb = texelFetch(u_bounds, ivec2(x + 1, 0), 0).r - u_tileOff;
  int i0 = max(0, int(floor(sa + 0.5))), i1 = min(u_n, int(ceil(sb + 0.5)));
  i1 = min(i1, i0 + 64);     // a cap the views offered do not reach (0.5 s on 235 px)
  float ra = u_rho0 + p * u_drho;
  int r0 = max(0, int(floor(ra))), r1 = min(384, int(ceil(ra + u_drho)));
  r1 = min(r1, r0 + 24);     // all 28 shanks on 450 px is 24 half-rows per pixel
  if (i1 <= i0 || r1 <= r0) discard;
  float best = 0.0;
  for (int r = r0; r < r1; r++) {
    int c = (191 - r / 2) * 2 + (r & 1);      // site row from the top -> channel
    for (int i = i0; i < i1; i++) {
      float v = texelFetch(u_a, ivec2(i, c), 0).r;
      if (u_diff == 1) v -= texelFetch(u_b, ivec2(i, c), 0).r;
      if (abs(v) > abs(best)) best = v;
    }
  }
  float gr = clamp(0.5 + 0.5 * best / u_range, 0.0, 1.0);
  o = vec4(gr, gr, gr, 1.0);
}`;

  function compileRaw(vs, fs) {
    const p = g.createProgram();
    for (const [type, src] of [[g.VERTEX_SHADER, vs], [g.FRAGMENT_SHADER, fs]]) {
      const s = g.createShader(type);
      g.shaderSource(s, src); g.compileShader(s);
      if (!g.getShaderParameter(s, g.COMPILE_STATUS)) throw new Error(g.getShaderInfoLog(s));
      g.attachShader(p, s);
    }
    g.linkProgram(p);
    if (!g.getProgramParameter(p, g.LINK_STATUS)) throw new Error(g.getProgramInfoLog(p));
    const u = {};
    for (let i = 0; i < g.getProgramParameter(p, g.ACTIVE_UNIFORMS); i++) {
      const name = g.getActiveUniform(p, i).name;
      u[name] = g.getUniformLocation(p, name);
    }
    return { p, u };
  }

  // ------------------------------------------------------------ loading
  async function load(base = "rawdata/") {
    try {
      const r = await fetch(base + "manifest.json");
      if (!r.ok) return;
      const m = await r.json();
      const cb = await (await fetch(base + m.clock)).arrayBuffer();
      const pairs = new Float64Array(cb), n = pairs.length / 2;
      clock = { s: new Float64Array(n), t: new Float64Array(n) };
      for (let i = 0; i < n; i++) { clock.s[i] = pairs[2 * i]; clock.t[i] = pairs[2 * i + 1]; }
      for (let i = 1; i < n; i++)
        if (!(clock.s[i] > clock.s[i - 1]) || !(clock.t[i] > clock.t[i - 1])) throw new Error("clock pairs are not increasing");
      m.base = base;
      // the drawing assumes NP2.0 geometry: two sites per 15 um row, channel c at row c >> 1
      for (const s of m.shanks) {
        s.ymax = Math.max(...s.y_um);
        if (m.n_channels !== 384 || s.y_um.some((y, c) => y !== 15 * (c >> 1)) ||
            s.x_um.some((x, c) => x !== s.x_um[c & 1]) || !(s.x_um[0] < s.x_um[1]))
          throw new Error("unexpected site geometry on shank " + s.shank);
      }
      man = m;
      covT = [timeAt(m.s0), timeAt(m.s0 + m.n_samples)];
      const nw = Math.max(1, Math.min(4, (navigator.hardwareConcurrency || 4) - 1));
      for (let i = 0; i < nw; i++) {
        const w = new Worker("raw_worker.js");
        w.busy = false;
        w.onmessage = (e) => finish(w, e.data);
        w.onerror = (e) => {                  // the worker itself failed: fail its job, free it
          console.error("raw worker", e.message);
          for (const [id, job] of jobs) if (job.worker === w) finish(w, { id, ok: false, error: "worker failed: " + e.message });
        };
        workers.push(w);
      }
    } catch (err) {
      console.warn("raw underlay not available:", err);
      man = null;
    }
  }
  function interp(xs, ys, x) {          // linear between pairs, extended linearly at the ends
    let a = 0, b = xs.length - 1;
    while (b - a > 1) { const m = (a + b) >> 1; if (xs[m] <= x) a = m; else b = m; }
    return ys[a] + (x - xs[a]) * (ys[b] - ys[a]) / (xs[b] - xs[a]);
  }
  const sampleAt = (t) => interp(clock.t, clock.s, t);
  const timeAt = (s) => interp(clock.s, clock.t, s);

  function initGLRaw() {
    cv = document.getElementById("rawcv");
    g = cv.getContext("webgl2", { antialias: false, alpha: true, premultipliedAlpha: false });
    if (!g) throw new Error("no second WebGL2 context for the raw underlay");
    P = compileRaw(VS, FRAG);
    vao = g.createVertexArray();
    boundsTex = g.createTexture();
    g.pixelStorei(g.UNPACK_ALIGNMENT, 1);
  }
  function texParams() {
    for (const [k, v] of [[g.TEXTURE_MIN_FILTER, g.NEAREST], [g.TEXTURE_MAG_FILTER, g.NEAREST],
                          [g.TEXTURE_WRAP_S, g.CLAMP_TO_EDGE], [g.TEXTURE_WRAP_T, g.CLAMP_TO_EDGE]])
      g.texParameteri(g.TEXTURE_2D, k, v);
  }

  // ------------------------------------------------------------ tiles
  const vkey = (i) => man.variants[i].key;
  const tileKey = (vi, sh, T) => `${vkey(vi)}|${sh}|${T}`;
  const chunkUrl = (vi, sh, T) => man.base + man.path.replace("{probe}", man.probe)
    .replace("{shank}", sh).replace("{variant}", vkey(vi)).replace("{tile}", T);
  function tileParams(vi, shank) {
    const v = man.variants[vi], s = man.shanks.find((q) => q.shank === shank);
    return { nch: man.n_channels, tile: man.tile, margin: man.margin, sos: man.highpass.sos,
             dead: s.dead, uv: s.uv_per_count, scale: v.scale, offset: v.offset };
  }
  function post(job) {
    const w = workers.find((q) => !q.busy);
    if (!w) return false;
    w.busy = true; job.worker = w; jobs.set(job.id, job);
    w.postMessage(job.msg);
    return true;
  }
  function pump() { while (queue.length && post(queue[queue.length - 1])) queue.pop(); }
  function call(msg) {                 // -> Promise of the worker's reply
    return new Promise((resolve, reject) => {
      const id = nextId++;
      queue.push({ id, msg: { ...msg, id }, resolve, reject });     // newest first
      pump();
    });
  }
  function finish(w, d) {
    w.busy = false;
    const job = jobs.get(d.id); jobs.delete(d.id);
    if (job) { if (d.ok) job.resolve(d); else job.reject(new Error(d.error)); }
    pump();
  }
  function tileData(vi, shank, T) {    // processed tile as Float32Array (channels x samples, uV)
    const nT = man.n_samples / man.tile;
    const urls = [T > 0 ? chunkUrl(vi, shank, T - 1) : null, chunkUrl(vi, shank, T),
                  T < nT - 1 ? chunkUrl(vi, shank, T + 1) : null];
    return call({ type: "tile", urls, p: tileParams(vi, shank) });
  }
  function wantTile(vi, shank, T) {    // -> texture, or null while it is being made
    const key = tileKey(vi, shank, T);
    let t = tiles.get(key);
    if (!t) {
      t = { state: "pending", tex: null, used: frame };
      tiles.set(key, t);
      tileData(vi, shank, T).then((d) => {
        if (!tiles.has(key)) return;                 // dropped meanwhile
        t.tex = g.createTexture();
        g.bindTexture(g.TEXTURE_2D, t.tex);
        g.texImage2D(g.TEXTURE_2D, 0, g.R16F, man.tile, man.n_channels, 0, g.RED, g.FLOAT, d.data);
        texParams();
        t.state = "ready"; t.ms = d.ms;
        requestRender();
      }).catch((err) => {
        if (err.message === "dropped") { tiles.delete(key); return; }   // asked for again if needed
        t.state = "error"; t.failedAt = performance.now();
        console.error("raw tile", key, err); requestRender();
      });
    } else if (t.state === "error" && performance.now() - t.failedAt > RETRY_MS) {
      tiles.delete(key);                             // try again on the next frame
      requestRender();
    }
    t.used = frame;
    return t.state === "ready" ? t.tex : null;
  }
  function trim() {
    // drop the oldest queued work when panning has outrun the workers, and the oldest textures
    while (queue.length > 48) queue.shift().reject(new Error("dropped"));
    if (tiles.size <= MAX_TEX) return;
    const old = [...tiles].filter(([, t]) => t.used !== frame && t.state !== "pending")
      .sort((a, b) => a[1].used - b[1].used);
    for (const [k, t] of old) {
      if (tiles.size <= MAX_TEX) break;
      if (t.tex) g.deleteTexture(t.tex);
      tiles.delete(k);
    }
  }

  // ------------------------------------------------------------ drawing
  // Called by render() before the spikes are drawn. Returns null, or the region (in cells of
  // the main canvas, fragment coordinates) under which voltage is drawn, with the opacity the
  // spike cells should have there.
  function draw() {
    info = { active: false, hint: "", loading: 0, failed: 0 };
    if (!man) return null;
    if (!cv) {
      try { initGLRaw(); }
      catch (err) { console.warn("raw underlay switched off:", err); man = null; cv = null; return null; }
    }
    const hide = () => { cv.style.display = "none"; return null; };
    if (!R.on || state.ymode !== "depth") return hide();
    const span = state.t1 - state.t0, L = LAY.depth, nspan = state.n1 - state.n0;
    const segs = L.segs.filter((s) => s.probe === man.probe && s.y1 > state.n0 && s.y0 < state.n1 &&
                                      man.shanks.some((q) => q.shank === s.shank));
    const ta = Math.max(state.t0, covT[0]), tb = Math.min(state.t1, covT[1]);
    if (!segs.length || tb <= ta) return hide();
    if (span > MAX_SPAN_S) { info.hint = `raw voltage is available here: zoom in to ${MAX_SPAN_S} s or less`; return hide(); }
    frame++;
    const b = state.block, PW = layout.W * b, PH = layout.H * b;
    if (cv.width !== PW || cv.height !== PH) { cv.width = PW; cv.height = PH; }
    Object.assign(cv.style, { display: "block", left: layout.x + "px", top: layout.y + "px",
                              width: layout.w + "px", height: layout.h + "px" });
    // sample at every pixel-column edge, relative to the first tile needed
    const sEdge = new Float64Array(PW + 1);
    for (let x = 0; x <= PW; x++) sEdge[x] = sampleAt(state.t0 + span * x / PW) - man.s0;
    const T0 = Math.max(0, Math.floor((sEdge[0] - 0.5) / man.tile));
    const T1 = Math.min(man.n_samples / man.tile - 1, Math.floor((sEdge[PW] + 0.5) / man.tile));
    const base = T0 * man.tile, rel = new Float32Array(PW + 1);
    for (let x = 0; x <= PW; x++) rel[x] = sEdge[x] - base;

    g.viewport(0, 0, PW, PH);
    g.disable(g.SCISSOR_TEST);
    g.clearColor(0.5, 0.5, 0.5, 1); g.clear(g.COLOR_BUFFER_BIT);
    g.useProgram(P.p); g.bindVertexArray(vao);
    g.activeTexture(g.TEXTURE2); g.bindTexture(g.TEXTURE_2D, boundsTex);
    g.texImage2D(g.TEXTURE_2D, 0, g.R32F, PW + 1, 1, 0, g.RED, g.FLOAT, rel);
    texParams();
    g.uniform1i(P.u.u_a, 0); g.uniform1i(P.u.u_b, 1); g.uniform1i(P.u.u_bounds, 2);
    g.uniform1i(P.u.u_PH, PH); g.uniform1i(P.u.u_n, man.tile);
    g.uniform1f(P.u.u_range, R.uv);
    const diff = R.diff && R.variant !== 0;
    g.uniform1i(P.u.u_diff, diff ? 1 : 0);
    g.enable(g.SCISSOR_TEST);
    // pixel columns of each tile: those whose central sample lies in it
    const colOf = (s) => {             // first column with center sample >= s
      let a = 0, c = PW;
      while (a < c) { const m = (a + c) >> 1; if (0.5 * (sEdge[m] + sEdge[m + 1]) < s) a = m + 1; else c = m; }
      return a;
    };
    let rowLo = PH, rowHi = 0;
    for (const sg of segs) {
      // pixel rows whose center lies in this shank's segment (row 0 at the top)
      const p0 = Math.max(0, Math.ceil((sg.y0 - state.n0) / nspan * PH - 0.5));
      const p1 = Math.min(PH, Math.ceil((sg.y1 - state.n0) / nspan * PH - 0.5));
      if (p1 <= p0) continue;
      rowLo = Math.min(rowLo, p0); rowHi = Math.max(rowHi, p1);
      // half-site rows: a site at height y starts (top10 - 75) / 10 - y um into the segment,
      // like a neuron at that height (buildLayouts); the top site is half-row 0
      const off = sg.y0 + (sg.top10 - 75) / 10 - man.shanks.find((q) => q.shank === sg.shank).ymax;
      g.uniform1f(P.u.u_rho0, (state.n0 - off) / 7.5);
      g.uniform1f(P.u.u_drho, nspan / PH / 7.5);
      for (let T = T0; T <= T1; T++) {
        const ta_ = wantTile(R.variant, sg.shank, T);
        const tb_ = diff ? wantTile(0, sg.shank, T) : null;
        if (!ta_ || (diff && !tb_)) {
          const failed = [tileKey(R.variant, sg.shank, T), tileKey(0, sg.shank, T)]
            .some((k) => tiles.get(k) && tiles.get(k).state === "error");
          if (failed) info.failed++; else info.loading++;
          continue;
        }
        const x0 = colOf(T * man.tile - 0.5), x1 = colOf((T + 1) * man.tile - 0.5);
        if (x1 <= x0) continue;
        g.scissor(x0, PH - p1, x1 - x0, p1 - p0);
        g.activeTexture(g.TEXTURE0); g.bindTexture(g.TEXTURE_2D, ta_);
        g.activeTexture(g.TEXTURE1); g.bindTexture(g.TEXTURE_2D, tb_ || ta_);
        g.uniform1f(P.u.u_tileOff, T * man.tile - base);
        g.drawArrays(g.TRIANGLES, 0, 3);
      }
    }
    g.disable(g.SCISSOR_TEST);
    trim();
    info.active = true;
    // the region, in cells of the main canvas: columns and rows whose centers are covered
    const W = layout.W, H = layout.H;
    const cx0 = Math.max(0, Math.ceil((ta - state.t0) / span * W - 0.5));
    const cx1 = Math.min(W, Math.ceil((tb - state.t0) / span * W - 0.5));
    const cy0 = Math.max(0, Math.ceil((segs[0].y0 - state.n0) / nspan * H - 0.5));
    const cy1 = Math.min(H, Math.ceil((segs[segs.length - 1].y1 - state.n0) / nspan * H - 0.5));
    return { rect: [cx0, H - cy1, cx1, H - cy0], alpha: ALPHAS[R.alpha] };
  }

  // ------------------------------------------------------------ keys, legend, status, URL
  function key(k) {
    if (!man || state.ymode !== "depth") return false;
    if (k === "v" || k === "V") { R.variant = (R.variant + 1) % man.variants.length; if (R.variant === 0) R.diff = false; }
    else if (k === "d" || k === "D") { R.diff = !R.diff; if (R.diff && R.variant === 0) R.variant = 1; }
    else if (k === "u" || k === "U") R.on = !R.on;
    else if (k === "o" || k === "O") R.alpha = (R.alpha + 1) % ALPHAS.length;
    else if (k === "[") R.uv = Math.max(5, R.uv / 2);
    else if (k === "]") R.uv = Math.min(3200, R.uv * 2);
    else return false;
    return true;
  }
  function variantText() {
    const v = man.variants[R.variant];
    const size = `${v.ratio.toFixed(1)}× smaller than the recorded data`;
    return R.diff && R.variant !== 0 ? `${v.label} minus lossless` : `${v.label} (${size})`;
  }
  function drawLegend(ctx, ox, OV_Y, OV_H) {
    if (!man) return;
    // where the raw data are, under the overview strip
    ctx.fillStyle = "#f2c14e";
    ctx.fillRect(ox(covT[0]), OV_Y + OV_H + 1, Math.max(2, ox(covT[1]) - ox(covT[0])), 2);
    if (!info.active) return;
    const { x: px, y: py } = layout;
    const lines = [`Raw voltage, probe ${man.probe + 1}: ${variantText()}`,
                   `Black −${R.uv} µV, white +${R.uv} µV · spikes drawn at ${Math.round(100 * ALPHAS[R.alpha])}% opacity` +
                   (info.loading ? ` · loading ${info.loading} tile${info.loading === 1 ? "" : "s"}` : "") +
                   (info.failed ? ` · ${info.failed} tile${info.failed === 1 ? "" : "s"} failed to load` : ""),
                   "Keys: V version · D difference from lossless · [ ] range · O spikes · U hide"];
    ctx.font = "12px Arial"; ctx.textAlign = "left"; ctx.textBaseline = "top";
    const w = Math.max(...lines.map((s) => ctx.measureText(s).width)) + 16;
    ctx.fillStyle = "rgba(14,17,21,0.82)"; ctx.fillRect(px + 8, py + 8, w, 16 * lines.length + 10);
    ctx.fillStyle = "#e6e9ed";
    lines.forEach((s, i) => ctx.fillText(s, px + 16, py + 14 + 16 * i));
  }
  const status = () => !man ? "" : info.active ? ` · raw voltage: ${variantText()}, ±${R.uv} µV`
    : info.hint ? ` · ${info.hint}` : "";
  function toHash(p) {
    if (!man) return;
    p.set("raw", R.on ? vkey(R.variant) : "off");
    p.set("rawdiff", R.diff ? 1 : 0); p.set("rawuv", R.uv); p.set("rawsp", R.alpha);
  }
  function fromHash(p) {
    if (!man) return;
    const v = p.get("raw");
    if (v === "off") R.on = false;
    else { const i = man.variants.findIndex((q) => q.key === v); if (i >= 0) R.variant = i; }
    R.diff = p.get("rawdiff") === "1" && R.variant !== 0;
    const u = parseFloat(p.get("rawuv")); if (u >= 5 && u <= 3200) R.uv = u;
    const a = parseInt(p.get("rawsp")); if (a >= 0 && a < ALPHAS.length) R.alpha = a;
  }

  // ------------------------------------------------------------ self-test (?rawtest=1)
  // Decoding: the sum of chunk 0 of every array against the manifest. Processing: one tile of
  // each version against make_raw_proto.py's process_tile.
  async function test() {
    if (!man) return { error: "no raw data (rawdata/manifest.json)" };
    const out = { decode: [], process: [] };
    for (const s of man.shanks) for (let vi = 0; vi < man.variants.length; vi++) {
      const d = await call({ type: "sum", url: chunkUrl(vi, s.shank, 0), nch: man.n_channels, tile: man.tile });
      out.decode.push({ shank: s.shank, variant: vkey(vi), ok: d.sum === s.chunk0_sum[vkey(vi)] });
    }
    const tt = man.test;
    for (let vi = 0; vi < man.variants.length; vi++) {
      const d = await tileData(vi, tt.shank, tt.tile);
      const ref = new Float32Array(await (await fetch(man.base + tt.files[vkey(vi)])).arrayBuffer());
      let maxAbs = 0, maxRef = 0;
      for (let i = 0; i < ref.length; i++) {
        maxAbs = Math.max(maxAbs, Math.abs(ref[i] - d.data[i])); maxRef = Math.max(maxRef, Math.abs(ref[i]));
      }
      out.process.push({ variant: vkey(vi), n: ref.length, sameLength: ref.length === d.data.length,
                         maxAbsDiff_uV: maxAbs, maxRef_uV: maxRef, ms: d.ms });
    }
    out.ok = out.decode.every((q) => q.ok) && out.process.every((q) => q.sameLength && q.maxAbsDiff_uV < 0.01);
    return out;
  }

  return { load, draw, key, drawLegend, status, toHash, fromHash, test,
           get available() { return !!man; }, get coverage() { return covT; }, get info() { return info; },
           sampleAt: (t) => sampleAt(t) };
})();
