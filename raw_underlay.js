"use strict";
// Raw-voltage underlay. When the neurons are plotted by depth and the window is short, the
// recorded voltage of each site is drawn under the spikes, on its own canvas below the main
// one. The data (raw_underlay/make_raw.py) are hosted in a bucket, one file per 0.2 s of one
// shank; a folder served locally (`python serve.py --raw <folder>`) or ?rawbase=<address>
// takes their place. Without the data the viewer works as before.
//
// Drawing. A tile is 0.2 s of one shank (384 sites), processed in a worker (raw_worker.js:
// high-pass, median reference) and kept as a texture in microvolts. Every device pixel shows
// the sample of largest magnitude among the samples and sites it covers.
//   Time:  the viewer's clock -> probe sample through the knots in clock_p<probe>.f64 (linear
//          between knots), the same map the spike times went through; every probe has its own.
//          Sample s is drawn over [s - 0.5, s + 0.5).
//   Depth: a site at height y occupies the same 15 um band as a neuron at that height. Its two
//          sites (x = 27 and 59 um) share the band, upper and lower half.
//   Motion: the sort ran on motion-corrected data, so a unit's depth is in corrected
//          coordinates: corrected(y) = recorded(y + d), with d the motion estimate of that
//          second at that depth (motion_p<probe>.f32; linear between its depth windows,
//          constant beyond the outer ones). The voltage is drawn accordingly: a pixel at depth
//          y shows the sites at y + d. Key M switches the shift off and on.

const RawUnderlay = (() => {
  const HOSTED = "https://pub-6f7885d91fac46b4a0aad6e22005454a.r2.dev/v1/";
  const MAX_SPAN_S = 0.5;        // longest window with the underlay
  const MAX_VIEW_TILES = 24;     // most tiles (shanks x 0.2 s) one view may need (about 1 MB each)
  const MAX_TEX = 48;            // tiles kept on the GPU (6000 x 384 half floats, 4.6 MB each)
  const MAX_FILES = 64;          // downloaded chunk files kept in memory (about 1 MB each)
  const RETRY_MS = 3000;         // a tile that failed to load is asked for again after this
  const ALPHAS = [1, 0.6, 0];    // opacity of the spike cells over the voltage (key O)
  const DEF = { on: true, uv: 100, alpha: 0, motion: true };
  let man = null, cv = null, g = null, P = null, boundsTex = null, vao = null;
  const R = { ...DEF };
  const probes = new Map();      // probe -> { state, s, t (clock knots), motion (Float32Array) }
  const tiles = new Map();       // "probe|shank|tile" -> { state, tex, used }
  const files = new Map();       // url -> Promise<ArrayBuffer>, most recently used last
  let workers = [], jobs = new Map(), queue = [], nextId = 1, frame = 0;
  let info = { active: false, hint: "", loading: 0, failed: 0 };

  const VS = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;
  const frag = (nMot) => `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D u_a;          // R16F, tile x 384: microvolts
uniform sampler2D u_bounds;     // R32F, (PW + 1) x 1: sample at each pixel-column edge, minus base
uniform float u_tileOff;        // first sample of this tile, minus base
uniform float u_rho0, u_drho;   // half-site rows from the top of the shank: at pixel row 0, per row
uniform float u_ytop;           // height (um) of the top edge of the top site's band
uniform float u_mot[${nMot}];   // motion estimate (um) at the depth windows, for this tile's second
uniform float u_motY0, u_motDy; // height of the first window, spacing (um)
uniform float u_range;          // microvolts at black and white
uniform int u_PH, u_n;
out vec4 o;
void main() {
  int x = int(gl_FragCoord.x);
  float p = float(u_PH - 1 - int(gl_FragCoord.y));
  float sa = texelFetch(u_bounds, ivec2(x, 0), 0).r - u_tileOff;
  float sb = texelFetch(u_bounds, ivec2(x + 1, 0), 0).r - u_tileOff;
  int i0 = max(0, int(floor(sa + 0.5))), i1 = min(u_n, int(ceil(sb + 0.5)));
  i1 = min(i1, i0 + 64);     // a cap the views offered do not reach (0.5 s on 235 px)
  float ra = u_rho0 + p * u_drho;
  // the pixel is at corrected height yc; it shows the sites at yc + d, which are d / 7.5
  // half-rows nearer the top
  float yc = u_ytop - (ra + 0.5 * u_drho) * 7.5;
  float u = clamp((yc - u_motY0) / u_motDy, 0.0, float(${nMot - 1}));
  int k = min(int(floor(u)), ${nMot - 2});
  ra -= mix(u_mot[k], u_mot[k + 1], u - float(k)) / 7.5;
  int r0 = max(0, int(floor(ra))), r1 = min(384, int(ceil(ra + u_drho)));
  r1 = min(r1, r0 + 24);     // all 28 shanks on 450 px is 24 half-rows per pixel
  if (i1 <= i0 || r1 <= r0) discard;
  float best = 0.0;
  for (int r = r0; r < r1; r++) {
    int c = (191 - r / 2) * 2 + (r & 1);      // site row from the top -> channel
    for (int i = i0; i < i1; i++) {
      float v = texelFetch(u_a, ivec2(i, c), 0).r;
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
      u[name.replace("[0]", "")] = g.getUniformLocation(p, name);
    }
    return { p, u };
  }

  // ------------------------------------------------------------ loading
  // Where the data are: ?rawbase=<address>; a folder served by serve.py --raw (only looked for
  // on a local server); otherwise the hosted copy. ?rawbase=off switches the underlay off.
  async function load() {
    const q = new URLSearchParams(location.search).get("rawbase");
    if (q === "off") return;
    const local = /^(localhost|127\.0\.0\.1)$/.test(location.hostname);
    const bases = q ? [q.replace(/\/?$/, "/")] : local ? ["rawdata/", HOSTED] : [HOSTED];
    for (const base of bases) {
      try {
        const r = await fetch(base + "manifest.json", { signal: AbortSignal.timeout(8000) });
        if (!r.ok) continue;
        const m = await r.json();
        check(m);
        m.base = base;
        man = m;
        break;
      } catch (err) {
        console.warn("raw underlay: no data at " + base + ":", err.message || err);
      }
    }
    if (!man) return;
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
  }
  // What the drawing relies on. NP2.0 geometry (two sites per 15 um row, channel c at row
  // c >> 1, left site first) is checked by make_raw.py for every shank.
  function check(m) {
    const bad = (msg) => { throw new Error("manifest: " + msg); };
    if (m.version !== 2) bad("version " + m.version);
    if (m.n_channels !== 384 || m.row_pitch_um !== 15) bad("unexpected site geometry");
    for (const pr of m.probes) {
      const mo = pr.motion, per = mo.dt_s * m.fs, first = (mo.t0_s - mo.dt_s / 2) * m.fs;
      // a tile must lie inside one time bin of the motion estimate
      if (per % m.tile !== 0 || first % m.tile !== 0) bad("motion bins do not line up with tiles");
      if (mo.depths_um.length < 2 || mo.depths_um.length > 16) bad("motion depth windows");
      if (mo.depths_um.length !== m.probes[0].motion.depths_um.length) bad("motion windows differ between probes");
      const dy = mo.depths_um[1] - mo.depths_um[0];
      if (mo.depths_um.some((y, i) => Math.abs(y - mo.depths_um[0] - i * dy) > 1e-6)) bad("motion windows unevenly spaced");
      if (!(pr.x_um[0] < pr.x_um[1]) || pr.shanks.some((s, i) => s.shank !== i)) bad("shanks");
    }
    if (m.probes.some((pr, i) => pr.probe !== i)) bad("probes");
  }
  // Clock knots and motion estimate of a probe, fetched when it is first needed.
  function probe(p) {
    let q = probes.get(p);
    if (q) {
      if (q.state === "error" && performance.now() - q.failedAt > RETRY_MS) probes.delete(p);
      else return q;
    }
    const pr = man.probes[p];
    q = { state: "pending", pr };
    probes.set(p, q);
    (async () => {
      const get = async (f) => {
        const r = await fetch(man.base + f, { signal: AbortSignal.timeout(30000) });
        if (!r.ok) throw new Error(`${f}: HTTP ${r.status}`);
        return r.arrayBuffer();
      };
      const [cb, mb] = await Promise.all([get(pr.clock), get(pr.motion.file)]);
      const pairs = new Float64Array(cb), n = pairs.length / 2;
      q.s = new Float64Array(n); q.t = new Float64Array(n);
      for (let i = 0; i < n; i++) { q.s[i] = pairs[2 * i]; q.t[i] = pairs[2 * i + 1]; }
      for (let i = 1; i < n; i++)
        if (!(q.s[i] > q.s[i - 1]) || !(q.t[i] > q.t[i - 1])) throw new Error("clock knots are not increasing");
      q.motion = new Float32Array(mb);
      const [a, b, c] = pr.motion.shape;
      if (q.motion.length !== a * b * c || a !== pr.shanks.length || c !== pr.motion.depths_um.length)
        throw new Error("motion file has the wrong size");
      q.state = "ready";
    })().catch((err) => {
      q.state = "error"; q.failedAt = performance.now();
      console.error("raw underlay, probe " + p, err);
    }).finally(() => requestRender());
    return q;
  }
  function interp(xs, ys, x) {          // linear between knots, extended linearly at the ends
    let a = 0, b = xs.length - 1;
    while (b - a > 1) { const m = (a + b) >> 1; if (xs[m] <= x) a = m; else b = m; }
    return ys[a] + (x - xs[a]) * (ys[b] - ys[a]) / (xs[b] - xs[a]);
  }
  // Motion estimate (um) at the depth windows for tile T of a shank: the time bin that holds it.
  function motionOf(q, shank, T) {
    const mo = q.pr.motion, [, nt, nd] = mo.shape;
    const b = Math.min(nt - 1, Math.max(0, Math.floor((T * man.tile / man.fs - mo.t0_s) / mo.dt_s + 0.5)));
    return q.motion.subarray((shank * nt + b) * nd, (shank * nt + b + 1) * nd);
  }

  function initGLRaw() {
    cv = document.getElementById("rawcv");
    g = cv.getContext("webgl2", { antialias: false, alpha: true, premultipliedAlpha: false });
    if (!g) throw new Error("no second WebGL2 context for the raw underlay");
    P = compileRaw(VS, frag(man.probes[0].motion.depths_um.length));
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
  const tileKey = (p, sh, T) => `${p}|${sh}|${T}`;
  const chunkUrl = (p, sh, T) => man.base + man.path.replace("{probe}", p).replace("{shank}", sh).replace("{tile}", T);
  function tileParams(p, shank) {
    const s = man.probes[p].shanks[shank], v = man.variant;
    return { nch: man.n_channels, tile: man.tile, margin: man.margin, sos: man.highpass.sos,
             dead: s.dead, uv: s.uv_per_count, scale: v.scale, offset: v.offset };
  }
  // A chunk file as downloaded (still compressed). Fetched here and not in the workers, so
  // that a file needed by several tiles (each tile also reads the edges of its two neighbors)
  // is downloaded once.
  function file(url) {
    let p = files.get(url);
    if (p) { files.delete(url); files.set(url, p); return p; }
    p = (async () => {
      const r = await fetch(url, { signal: AbortSignal.timeout(30000) });
      if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
      return r.arrayBuffer();
    })();
    files.set(url, p);
    p.catch(() => { if (files.get(url) === p) files.delete(url); });
    while (files.size > MAX_FILES) files.delete(files.keys().next().value);
    return p;
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
  async function tileData(p, shank, T) {    // processed tile as Float32Array (channels x samples, uV)
    const pr = man.probes[p];
    const urls = [T > pr.tile_first ? chunkUrl(p, shank, T - 1) : null, chunkUrl(p, shank, T),
                  T < pr.tile_last ? chunkUrl(p, shank, T + 1) : null];
    const bufs = await Promise.all(urls.map((u) => u ? file(u) : null));
    return call({ type: "tile", urls, bufs, p: tileParams(p, shank) });
  }
  function wantTile(p, shank, T) {     // -> texture, or null while it is being made
    const key = tileKey(p, shank, T);
    let t = tiles.get(key);
    if (!t) {
      t = { state: "pending", tex: null, used: frame };
      tiles.set(key, t);
      tileData(p, shank, T).then((d) => {
        if (tiles.get(key) !== t) return;            // dropped meanwhile
        t.tex = g.createTexture();
        g.bindTexture(g.TEXTURE_2D, t.tex);
        g.texImage2D(g.TEXTURE_2D, 0, g.R16F, man.tile, man.n_channels, 0, g.RED, g.FLOAT, d.data);
        texParams();
        t.state = "ready"; t.ms = d.ms;
        requestRender();
      }).catch((err) => {
        if (tiles.get(key) !== t) return;
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
    const segs = L.segs.filter((s) => s.y1 > state.n0 && s.y0 < state.n1 && man.probes[s.probe] &&
                                      man.probes[s.probe].shanks[s.shank]);
    if (!segs.length) return hide();
    if (span > MAX_SPAN_S) { info.hint = `raw voltage: zoom in to ${MAX_SPAN_S} s or less`; return hide(); }
    // tiles a view of this length can need per shank, whatever the probe's clock
    const perShank = Math.ceil(span * man.fs / man.tile) + 1;
    if (segs.length * perShank > MAX_VIEW_TILES) {
      const ns = Math.floor(MAX_VIEW_TILES / perShank);
      info.hint = `raw voltage: zoom in to ${ns} shanks or fewer` + (perShank > 2 ? ", or to a shorter window" : "");
      return hide();
    }
    frame++;
    const b = state.block, PW = layout.W * b, PH = layout.H * b;
    if (cv.width !== PW || cv.height !== PH) { cv.width = PW; cv.height = PH; }
    Object.assign(cv.style, { display: "block", left: layout.x + "px", top: layout.y + "px",
                              width: layout.w + "px", height: layout.h + "px" });
    g.viewport(0, 0, PW, PH);
    g.disable(g.SCISSOR_TEST);
    g.clearColor(0.5, 0.5, 0.5, 1); g.clear(g.COLOR_BUFFER_BIT);
    g.useProgram(P.p); g.bindVertexArray(vao);
    g.uniform1i(P.u.u_a, 0); g.uniform1i(P.u.u_bounds, 2);
    g.uniform1i(P.u.u_PH, PH); g.uniform1i(P.u.u_n, man.tile);
    g.uniform1f(P.u.u_range, R.uv);
    const ymax = (man.n_channels / 2 - 1) * man.row_pitch_um;      // height of the top site
    g.uniform1f(P.u.u_ytop, ymax + man.row_pitch_um / 2);
    g.enable(g.SCISSOR_TEST);
    const sEdge = new Float64Array(PW + 1), rel = new Float32Array(PW + 1);
    let lastProbe = -1, T0 = 0, T1 = -1, base = 0;
    const colOf = (s) => {             // first column with center sample >= s
      let a = 0, c = PW;
      while (a < c) { const m = (a + c) >> 1; if (0.5 * (sEdge[m] + sEdge[m + 1]) < s) a = m + 1; else c = m; }
      return a;
    };
    for (const sg of segs) {
      const q = probe(sg.probe), pr = man.probes[sg.probe];
      if (q.state !== "ready") { if (q.state === "error") info.failed++; else info.loading++; continue; }
      if (sg.probe !== lastProbe) {
        // sample at every pixel-column edge on this probe's clock, relative to the first tile needed
        lastProbe = sg.probe;
        for (let x = 0; x <= PW; x++) sEdge[x] = interp(q.t, q.s, state.t0 + span * x / PW);
        T0 = Math.max(pr.tile_first, Math.floor((sEdge[0] - 0.5) / man.tile));
        T1 = Math.min(pr.tile_last, Math.floor((sEdge[PW] + 0.5) / man.tile));
        base = T0 * man.tile;
        for (let x = 0; x <= PW; x++) rel[x] = sEdge[x] - base;
        g.activeTexture(g.TEXTURE2); g.bindTexture(g.TEXTURE_2D, boundsTex);
        g.texImage2D(g.TEXTURE_2D, 0, g.R32F, PW + 1, 1, 0, g.RED, g.FLOAT, rel);
        texParams();
        const mo = pr.motion;
        g.uniform1f(P.u.u_motY0, mo.depths_um[0]);
        g.uniform1f(P.u.u_motDy, mo.depths_um[1] - mo.depths_um[0]);
      }
      // pixel rows whose center lies in this shank's segment (row 0 at the top)
      const p0 = Math.max(0, Math.ceil((sg.y0 - state.n0) / nspan * PH - 0.5));
      const p1 = Math.min(PH, Math.ceil((sg.y1 - state.n0) / nspan * PH - 0.5));
      if (p1 <= p0) continue;
      // half-site rows: a site at height y starts (top10 - 75) / 10 - y um into the segment,
      // like a neuron at that height (buildLayouts); the top site is half-row 0
      const off = sg.y0 + (sg.top10 - 75) / 10 - ymax;
      g.uniform1f(P.u.u_rho0, (state.n0 - off) / 7.5);
      g.uniform1f(P.u.u_drho, nspan / PH / 7.5);
      for (let T = T0; T <= T1; T++) {
        const tex = wantTile(sg.probe, sg.shank, T);
        if (!tex) {
          const t = tiles.get(tileKey(sg.probe, sg.shank, T));
          if (t && t.state === "error") info.failed++; else info.loading++;
          continue;
        }
        const x0 = colOf(T * man.tile - 0.5), x1 = colOf((T + 1) * man.tile - 0.5);
        if (x1 <= x0) continue;
        g.scissor(x0, PH - p1, x1 - x0, p1 - p0);
        g.activeTexture(g.TEXTURE0); g.bindTexture(g.TEXTURE_2D, tex);
        g.uniform1f(P.u.u_tileOff, T * man.tile - base);
        const m = motionOf(q, sg.shank, T);
        g.uniform1fv(P.u.u_mot, R.motion ? m : new Float32Array(m.length));
        g.drawArrays(g.TRIANGLES, 0, 3);
      }
    }
    g.disable(g.SCISSOR_TEST);
    trim();
    info.active = true;
    // the region, in cells of the main canvas: all columns, and the rows whose centers are
    // in the shanks shown
    const W = layout.W, H = layout.H;
    const cy0 = Math.max(0, Math.ceil((segs[0].y0 - state.n0) / nspan * H - 0.5));
    const cy1 = Math.min(H, Math.ceil((segs[segs.length - 1].y1 - state.n0) / nspan * H - 0.5));
    return { rect: [0, H - cy1, W, H - cy0], alpha: ALPHAS[R.alpha] };
  }

  // ------------------------------------------------------------ keys, legend, status, URL
  function key(k) {
    if (!man || state.ymode !== "depth") return false;
    if (k === "u" || k === "U") R.on = !R.on;
    else if (k === "o" || k === "O") R.alpha = (R.alpha + 1) % ALPHAS.length;
    else if (k === "m" || k === "M") R.motion = !R.motion;
    else if (k === "[") R.uv = Math.max(5, R.uv / 2);
    else if (k === "]") R.uv = Math.min(3200, R.uv * 2);
    else return false;
    return true;
  }
  function variantText() {
    const v = man.variant, step = v.scale * man.probes[0].shanks[0].uv_per_count;
    return `stored in steps of ${step.toFixed(0)} µV` + (v.ratio ? `, ${v.ratio.toFixed(1)}× smaller than recorded` : "");
  }
  function drawLegend(ctx) {
    if (!man || !info.active) return;
    const { x: px, y: py } = layout;
    const lines = [`Raw voltage (${variantText()})` + (R.motion ? "" : " · motion shift off"),
                   `Black −${R.uv} µV, white +${R.uv} µV · spikes drawn at ${Math.round(100 * ALPHAS[R.alpha])}% opacity` +
                   (info.loading ? ` · loading ${info.loading}` : "") +
                   (info.failed ? ` · ${info.failed} failed to load` : ""),
                   "Keys: [ ] range · O spikes · M motion shift · U hide"];
    ctx.font = "12px Arial"; ctx.textAlign = "left"; ctx.textBaseline = "top";
    const w = Math.max(...lines.map((s) => ctx.measureText(s).width)) + 16;
    ctx.fillStyle = "rgba(14,17,21,0.82)"; ctx.fillRect(px + 8, py + 8, w, 16 * lines.length + 10);
    ctx.fillStyle = "#e6e9ed";
    lines.forEach((s, i) => ctx.fillText(s, px + 16, py + 14 + 16 * i));
  }
  const status = () => !man ? "" : info.active ? ` · raw voltage ±${R.uv} µV` : info.hint ? ` · ${info.hint}` : "";
  // Only what differs from the defaults goes into the URL.
  function toHash(p) {
    if (!R.on) p.set("raw", "off");
    if (R.uv !== DEF.uv) p.set("rawuv", R.uv);
    if (R.alpha !== DEF.alpha) p.set("rawsp", R.alpha);
    if (!R.motion) p.set("rawmotion", 0);
  }
  function fromHash(p) {
    R.on = p.get("raw") !== "off";
    const u = parseFloat(p.get("rawuv")); R.uv = u >= 5 && u <= 3200 ? u : DEF.uv;
    const a = parseInt(p.get("rawsp")); R.alpha = a >= 0 && a < ALPHAS.length ? a : DEF.alpha;
    R.motion = p.get("rawmotion") !== "0";
  }

  // ------------------------------------------------------------ self-test (?rawtest=1)
  // Clock and motion files of every probe load and pass their checks. Decoding: the sum of one
  // chunk of every shank against the manifest. Processing: the tiles listed in the manifest
  // against make_raw_proto.py's process_tile.
  async function test() {
    if (!man) return { error: "no raw data" };
    const out = { base: man.base, decode: [], process: [], probes: [] };
    for (const pr of man.probes) {
      probe(pr.probe);
      for (let i = 0; i < 300 && probes.get(pr.probe).state === "pending"; i++) await new Promise((r) => setTimeout(r, 100));
      out.probes.push({ probe: pr.probe, ok: probes.get(pr.probe).state === "ready" });
      for (const s of pr.shanks) {
        if (!s.check) { out.decode.push({ probe: pr.probe, shank: s.shank, ok: false, missing: true }); continue; }
        const url = chunkUrl(pr.probe, s.shank, s.check.tile);
        try {
          const d = await call({ type: "sum", url, buf: await file(url), nch: man.n_channels, tile: man.tile });
          out.decode.push({ probe: pr.probe, shank: s.shank, ok: d.sum === s.check.sum });
        } catch (err) { out.decode.push({ probe: pr.probe, shank: s.shank, ok: false, error: String(err.message || err) }); }
      }
    }
    for (const tt of man.test || []) {
      const d = await tileData(tt.probe, tt.shank, tt.tile);
      const ref = new Float32Array(await (await fetch(man.base + tt.file)).arrayBuffer());
      let maxAbs = 0, maxRef = 0;
      for (let i = 0; i < ref.length; i++) {
        maxAbs = Math.max(maxAbs, Math.abs(ref[i] - d.data[i])); maxRef = Math.max(maxRef, Math.abs(ref[i]));
      }
      out.process.push({ probe: tt.probe, shank: tt.shank, tile: tt.tile, n: ref.length,
                         sameLength: ref.length === d.data.length, maxAbsDiff_uV: maxAbs, maxRef_uV: maxRef, ms: d.ms });
    }
    out.ok = out.probes.every((q) => q.ok) && out.decode.every((q) => q.ok) && out.process.length > 0 &&
      out.process.every((q) => q.sameLength && q.maxAbsDiff_uV < 0.01);
    return out;
  }

  return { load, draw, key, drawLegend, status, toHash, fromHash, test,
           get available() { return !!man; }, get info() { return info; }, get settings() { return { ...R }; },
           get manifest() { return man; },
           // for checks: probe sample at a viewer time, and the motion estimate used for a tile
           sampleAt: (p, t) => { const q = probes.get(p); return q && q.state === "ready" ? interp(q.t, q.s, t) : NaN; },
           motionAt: (p, shank, T) => { const q = probes.get(p); return q && q.state === "ready" ? Array.from(motionOf(q, shank, T)) : null; } };
})();
