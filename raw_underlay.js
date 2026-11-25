"use strict";
// Raw-voltage underlay. When the neurons are plotted by depth and the window is short, the
// recorded voltage of each site is drawn under the spikes, on its own canvas below the main
// one. The data are hosted in a bucket; a folder served locally (`python serve.py --raw
// <folder> --disp <folder>`) or ?rawbase=<address> takes their place. Without the data the
// viewer works as before.
//
// Two sources, switched with key V:
//   display copy (raw_underlay/make_display.py; the default when it is there): already
//     high-passed and median-referenced, one byte per value in steps of a few microvolts.
//     Level 0 has every sample in tiles of 50 ms; level 1 has, for every block of 8 samples,
//     the sample of largest magnitude, in tiles of 0.5 s, and is used for longer windows.
//     A tile is one small file and needs no processing.
//   recorded samples (raw_underlay/make_raw.py): tiles of 0.2 s, processed in a worker
//     (raw_worker.js: high-pass, median reference); a tile also needs its two neighbors in
//     time for the filter's edges.
//
// Drawing. A tile of one shank (384 sites) is kept as a texture. Every device pixel shows the
// value of largest magnitude among the samples (or blocks) and sites it covers.
//   Time:  the viewer's clock -> probe sample through the knots in clock_p<probe>.f64 (linear
//          between knots), the same map the spike times went through; every probe has its own.
//          Sample s is drawn over [s - 0.5, s + 0.5). At level 1 a pixel shows every block that
//          holds one of its samples. A pixel column is drawn with the tile that holds its first
//          sample, and reads on into the next tile where it reaches across the edge.
//   Depth: a site at height y occupies the same 15 um band as a neuron at that height. Its two
//          sites (x = 27 and 59 um) share the band, upper and lower half.
//   Motion: the sort ran on motion-corrected data, so a unit's depth is in corrected
//          coordinates: corrected(y) = recorded(y + d), with d the motion estimate of that
//          second at that depth (motion_p<probe>.f32; linear between its depth windows,
//          constant beyond the outer ones). The voltage is drawn accordingly: a pixel at depth
//          y shows the sites at y + d. Key M switches the shift off and on.
// Tiles next to the view in time are fetched ahead once the view is complete, so that panning
// does not wait.

const RawUnderlay = (() => {
  const HOSTED = "https://pub-6f7885d91fac46b4a0aad6e22005454a.r2.dev/v1/";
  const TEX_BYTES = 230e6;       // tiles kept on the GPU (a 0.2 s tile of recorded samples is 4.6 MB,
                                 // a tile of the display copy 0.6 to 0.7 MB)
  const MAX_FILES = 96;          // downloaded files kept in memory (0.3 to 1.2 MB each)
  const MAX_WAIT = 160;          // downloads waiting to start; older ones are given up
  const FETCH_MS = 30000;        // a download may take this long from the moment it starts
  const RETRY_MS = 3000;         // a tile that failed to load is asked for again after this
  const MAN_RETRY_MS = 15000;    // the data's index is looked for again after this
  const HINT_SPAN_S = 0.5;       // in other layouts, say how to get the voltage from this zoom on
  const ALPHAS = [1, 0.6, 0];    // opacity of the spike cells over the voltage (key O)
  const DEF = { on: true, uv: 100, alpha: 0, motion: true, src: "disp" };
  let man = null, dman = null, cv = null, g = null, P = null, boundsTex = null, vao = null;
  let off = false, manTry = { busy: false, at: -1e9 };
  let dispBase = null, dispTry = { busy: false, at: -1e9, wait: MAN_RETRY_MS };   // looked for less and less often
  let RAWLV = null, DLV = [];    // what can be drawn: the recorded samples, the display copy's levels
  const R = { ...DEF };
  const probes = new Map();      // probe -> { state, s, t (clock knots), motion (Float32Array) }
  const tiles = new Map();       // "level|probe|shank|tile" -> { state, tex, used, keep, bytes }
  const files = new Map();       // url -> Promise<ArrayBuffer>, most recently used last
  let workers = [], jobs = new Map(), queue = [], nextId = 1, frame = 0, texBytes = 0;
  let info = { active: false, hint: "", loading: 0, failed: 0 };
  let shown = null;              // the level of the last drawing
  let loadNow = null, loadLast = null;   // downloads of the view being completed, and of the last one completed
  let retryTimer = null;

  const VS = `#version 300 es
void main() {
  vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;
  const frag = (nMot) => `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D u_a;          // values x 384: R16F microvolts, or R8_SNORM counts / 127
uniform sampler2D u_b;          // the next tile in time, for a column that reaches across the edge
uniform int u_hasB;
uniform float u_scale;          // microvolts per unit of u_a
uniform int u_k;                // samples per value (1, or the block of level 1)
uniform sampler2D u_bounds;     // R32F, (PW + 1) x 1: sample at each pixel-column edge, minus base
uniform float u_tileOff;        // first sample of this tile, minus base
uniform float u_rho0, u_drho;   // half-site rows from the top of the shank: at pixel row 0, per row
uniform float u_ytop;           // height (um) of the top edge of the top site's band
uniform float u_mot[${nMot}];   // motion estimate (um) at the depth windows, for this tile's second
uniform float u_motY0, u_motDy; // height of the first window, spacing (um)
uniform float u_range;          // microvolts at black and white
uniform int u_PH, u_n;          // u_n: values per tile
out vec4 o;
void main() {
  int x = int(gl_FragCoord.x);
  float p = float(u_PH - 1 - int(gl_FragCoord.y));
  float sa = texelFetch(u_bounds, ivec2(x, 0), 0).r - u_tileOff;
  float sb = texelFetch(u_bounds, ivec2(x + 1, 0), 0).r - u_tileOff;
  int nv = u_hasB == 1 ? 2 * u_n : u_n;                     // values that can be read
  int i0 = max(0, int(floor(sa + 0.5))), i1 = min(nv * u_k, int(ceil(sb + 0.5)));
  int j0 = i0 / u_k, j1 = min(nv, (i1 + u_k - 1) / u_k);    // the values that hold these samples
  j1 = min(j1, j0 + 64);     // a cap the views offered do not reach (chooseLevel: 60 per column)
  float ra = u_rho0 + p * u_drho;
  // the pixel is at corrected height yc; it shows the sites at yc + d, which are d / 7.5
  // half-rows nearer the top
  float yc = u_ytop - (ra + 0.5 * u_drho) * 7.5;
  float u = clamp((yc - u_motY0) / u_motDy, 0.0, float(${nMot - 1}));
  int k = min(int(floor(u)), ${nMot - 2});
  ra -= mix(u_mot[k], u_mot[k + 1], u - float(k)) / 7.5;
  int r0 = max(0, int(floor(ra))), r1 = min(384, int(ceil(ra + u_drho)));
  r1 = min(r1, r0 + 24);     // all 28 shanks on 450 px is 24 half-rows per pixel
  if (i1 <= i0 || j1 <= j0 || r1 <= r0) discard;
  float best = 0.0;
  for (int r = r0; r < r1; r++) {
    int c = (191 - r / 2) * 2 + (r & 1);      // site row from the top -> channel
    for (int i = j0; i < j1; i++) {
      float v = i < u_n ? texelFetch(u_a, ivec2(i, c), 0).r : texelFetch(u_b, ivec2(i - u_n, c), 0).r;
      if (abs(v) > abs(best)) best = v;
    }
  }
  float gr = clamp(0.5 + 0.5 * best * u_scale / u_range, 0.0, 1.0);
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
  // The display copy is looked for next to it (d1/ beside v1/, dispdata/ beside rawdata/) or at
  // ?dispbase=<address>; ?dispbase=off does without it. If the index cannot be fetched (a slow
  // connection while the spikes load, say), draw() looks again later.
  async function load() {
    try { await find(); }
    catch (err) { console.warn("raw underlay: not available yet:", err); man = null; }
    manTry.at = performance.now();
  }
  async function getJson(url) {
    const r = await fetch(url, { signal: AbortSignal.timeout(20000) });
    if (!r.ok) throw new Error(`HTTP ${r.status}`);
    return r.json();
  }
  async function find() {
    const par = new URLSearchParams(location.search), q = par.get("rawbase");
    if (q === "off") { off = true; return; }
    const slash = (u) => u.replace(/\/?$/, "/");
    const local = /^(localhost|127\.0\.0\.1)$/.test(location.hostname);
    const bases = q ? [slash(q)] : local ? ["rawdata/", HOSTED] : [HOSTED];
    let m = null;
    for (const base of bases) {
      try {
        m = await getJson(base + "manifest.json");
        check(m);
        m.base = base;
        break;
      } catch (err) {
        m = null;
        console.warn("raw underlay: no data at " + base + ":", err.message || err);
      }
    }
    if (!m) return;
    const dq = par.get("dispbase");
    dispBase = dq === "off" ? null : dq ? slash(dq) : m.base === "rawdata/" ? "dispdata/" :
      /v1\/$/.test(m.base) ? m.base.replace(/v1\/$/, "d1/") : null;
    const nw = Math.max(1, Math.min(4, (navigator.hardwareConcurrency || 4) - 1));
    for (let i = workers.length; i < nw; i++) {
      const w = new Worker("raw_worker.js");
      w.busy = false;
      w.onmessage = (e) => finish(w, e.data);
      w.onerror = (e) => {                  // the worker itself failed: fail its job and stop using it
        console.error("raw worker", e.message);
        workers = workers.filter((x) => x !== w);
        for (const [id, job] of jobs) if (job.worker === w) finish(w, { id, ok: false, error: "worker failed: " + e.message });
        if (!workers.length) while (queue.length) queue.pop().reject(new Error("no worker is left"));
      };
      workers.push(w);
    }
    man = m;
    RAWLV = { id: "r", src: "raw", block: 1, n: m.tile, spt: m.tile, scale: 1, bytes: m.tile * m.n_channels * 2,
              maxSpan: 0.5, maxTiles: 24, maxFetch: 6, range: (p) => m.probes[p] };
    await findDisplay();
  }
  // The display copy's index. Without it the recorded samples are drawn, and draw() looks again.
  async function findDisplay() {
    if (!dispBase || dman || dispTry.busy) return;
    dispTry.busy = true;
    try {
      const d = await getJson(dispBase + "manifest.json");
      checkDisplay(d, man);
      d.base = dispBase;
      DLV = d.levels.map((l, i) => ({
        id: "d" + l.level, src: "disp", block: l.block, n: l.tile, spt: l.tile * l.block,
        scale: d.clip * d.step_uv, step: d.step_uv, bytes: l.tile * man.n_channels, path: l.path,
        maxSpan: i === 0 ? 0.5 : 4, maxTiles: i === 0 ? 48 : 24, maxFetch: 12, range: (p) => d.probes[p].levels[i] }));
      dman = d;
    } catch (err) {
      if (dispTry.wait === MAN_RETRY_MS)
        console.warn("raw underlay: no display copy at " + dispBase + " (the recorded samples are used):", err.message || err);
      dispTry.wait = Math.min(2 * dispTry.wait, 300000);
    }
    dispTry.at = performance.now(); dispTry.busy = false;
  }
  // The display copy must describe the same probes, and a tile of either level must lie inside
  // one time bin of the motion estimate.
  function checkDisplay(d, m) {
    const bad = (msg) => { throw new Error("display manifest: " + msg); };
    if (d.version !== 1 || d.kind !== "display") bad("version " + d.version);
    if (d.n_channels !== m.n_channels || d.fs !== m.fs) bad("does not match the recorded samples");
    if (d.levels.length !== 2 || d.levels[0].block !== 1 || !(d.levels[1].block > 1)) bad("levels");
    if (!(d.step_uv > 0) || d.clip !== 127) bad("scaling");
    if (d.probes.length !== m.probes.length) bad("probes");
    d.probes.forEach((pr, i) => {
      const mo = m.probes[i].motion, per = mo.dt_s * m.fs, first = (mo.t0_s - mo.dt_s / 2) * m.fs;
      if (pr.probe !== i || pr.shanks.length !== m.probes[i].shanks.length) bad("probes");
      for (const l of d.levels)
        if (per % (l.tile * l.block) !== 0 || first % (l.tile * l.block) !== 0) bad("motion bins do not line up with tiles");
    });
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
      setTimeout(requestRender, RETRY_MS + 50);          // a render asks for it again
    }).finally(() => requestRender());
    return q;
  }
  function interp(xs, ys, x) {          // linear between knots, extended linearly at the ends
    let a = 0, b = xs.length - 1;
    while (b - a > 1) { const m = (a + b) >> 1; if (xs[m] <= x) a = m; else b = m; }
    return ys[a] + (x - xs[a]) * (ys[b] - ys[a]) / (xs[b] - xs[a]);
  }
  // Motion estimate (um) at the depth windows for the tile of a shank that starts at a sample:
  // the time bin that holds it.
  function motionOf(q, shank, sample) {
    const mo = q.pr.motion, [, nt, nd] = mo.shape;
    const b = Math.min(nt - 1, Math.max(0, Math.floor((sample / man.fs - mo.t0_s) / mo.dt_s + 0.5)));
    return q.motion.subarray((shank * nt + b) * nd, (shank * nt + b + 1) * nd);
  }

  function initGLRaw() {
    cv = document.getElementById("rawcv");
    g = cv.getContext("webgl2", { antialias: false, alpha: true, premultipliedAlpha: false });
    if (!g) throw new Error("no second WebGL2 context for the raw underlay");
    if (g.getParameter(g.MAX_TEXTURE_SIZE) < man.tile) throw new Error("the textures of this GPU are too small for a tile");
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
  const tileKey = (lv, p, sh, T) => `${lv.id}|${p}|${sh}|${T}`;
  const chunkUrl = (p, sh, T) => man.base + man.path.replace("{probe}", p).replace("{shank}", sh).replace("{tile}", T);
  const dispUrl = (lv, p, sh, T) => dman.base + lv.path.replace("{probe}", p).replace("{shank}", sh).replace("{tile}", T);
  function tileParams(p, shank) {
    const s = man.probes[p].shanks[shank], v = man.variant;
    return { nch: man.n_channels, tile: man.tile, margin: man.margin, sos: man.highpass.sos,
             dead: s.dead, uv: s.uv_per_count, scale: v.scale, offset: v.offset };
  }
  // A file as downloaded (still compressed). Fetched here and not in the workers, so that a
  // file needed by several tiles (a tile of recorded samples also reads the edges of its two
  // neighbors) is downloaded once. A limited number of downloads run at a time, the one asked
  // for last first, and each has FETCH_MS from the moment it starts.
  const waiting = [];            // downloads not started yet, the next one last
  let fetching = 0, maxFetch = 6;
  function pumpFiles() {
    while (fetching < maxFetch && waiting.length) {
      const w = waiting.pop();
      fetching++;
      (async () => {
        const r = await fetch(w.url, { signal: AbortSignal.timeout(FETCH_MS) });
        if (!r.ok) throw new Error(`${w.url}: HTTP ${r.status}`);
        const buf = await r.arrayBuffer();
        if (loadNow) { loadNow.files++; loadNow.bytes += buf.byteLength; }
        return buf;
      })().then(w.resolve, w.reject).finally(() => { fetching--; pumpFiles(); });
    }
  }
  // forView: asked for by a tile of the view (or one fetched ahead), so it may be given up
  // when the view has moved on.
  function file(url, forView = false) {
    let p = files.get(url);
    if (p) {
      files.delete(url); files.set(url, p);
      const i = waiting.findIndex((w) => w.url === url);        // asked for again: next in line
      if (i >= 0) { const w = waiting.splice(i, 1)[0]; w.forView = w.forView && forView; waiting.push(w); }
      return p;
    }
    p = new Promise((resolve, reject) => waiting.push({ url, resolve, reject, forView }));
    pumpFiles();
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
  // -> Promise of the worker's reply. tile: the key of the view's tile this is for, if any; the
  // job is then given up if it is still waiting when the tile has left the view.
  function call(msg, tile = null) {
    if (!workers.length) return Promise.reject(new Error("no worker is left"));
    return new Promise((resolve, reject) => {
      const id = nextId++;
      queue.push({ id, msg: { ...msg, id }, resolve, reject, tile });     // newest first
      pump();
    });
  }
  function finish(w, d) {
    w.busy = false;
    const job = jobs.get(d.id); jobs.delete(d.id);
    if (job) { if (d.ok) job.resolve(d); else job.reject(new Error(d.error)); }
    pump();
  }
  // The files of a tile. Recorded samples: its own and, for the filter's margins, its
  // neighbors in time. Display copy: its own.
  function tileUrls(lv, p, shank, T) {
    if (lv.src === "disp") return [dispUrl(lv, p, shank, T)];
    const pr = man.probes[p];
    return [T > pr.tile_first ? chunkUrl(p, shank, T - 1) : null, chunkUrl(p, shank, T),
            T < pr.tile_last ? chunkUrl(p, shank, T + 1) : null];
  }
  // -> { data }: a tile's values, channels x values. Recorded samples: Float32Array of
  // microvolts, processed in a worker. Display copy: Int8Array of counts (microvolts / step).
  async function tileData(lv, p, shank, T, forView = false) {
    const urls = tileUrls(lv, p, shank, T);
    const bufs = await Promise.all(urls.map((u) => u ? file(u, forView) : null));
    const key = forView ? tileKey(lv, p, shank, T) : null;
    if (key && !inView(key)) throw new Error("dropped");
    if (lv.src === "disp") return call({ type: "dtile", url: urls[0], buf: bufs[0], nch: man.n_channels, n: lv.n }, key);
    return call({ type: "tile", urls, bufs, p: tileParams(p, shank) }, key);
  }
  // -> texture, or null while it is being made. ahead: fetched ahead of the view, not drawn now.
  function wantTile(lv, p, shank, T, ahead = false) {
    const key = tileKey(lv, p, shank, T);
    let t = tiles.get(key);
    if (t && t.state === "error" && !ahead && performance.now() - t.failedAt > RETRY_MS) {
      tiles.delete(key); t = null;                   // try again
    }
    if (!t) {
      t = { state: "pending", tex: null, used: -1, keep: frame, bytes: 0, urls: tileUrls(lv, p, shank, T) };
      tiles.set(key, t);
      tileData(lv, p, shank, T, true).then((d) => {
        if (tiles.get(key) !== t) return;            // dropped meanwhile
        if (t.keep !== frame) { tiles.delete(key); return; }   // not wanted any more: only its files are kept
        t.tex = g.createTexture();
        g.bindTexture(g.TEXTURE_2D, t.tex);
        if (lv.src === "disp") g.texImage2D(g.TEXTURE_2D, 0, g.R8_SNORM, lv.n, man.n_channels, 0, g.RED, g.BYTE, d.data);
        else g.texImage2D(g.TEXTURE_2D, 0, g.R16F, lv.n, man.n_channels, 0, g.RED, g.FLOAT, d.data);
        texParams();
        t.bytes = lv.bytes; texBytes += lv.bytes;
        t.state = "ready"; t.ms = d.ms;
        requestRender();
      }).catch((err) => {
        if (tiles.get(key) !== t) return;
        if (err.message === "dropped") { tiles.delete(key); requestRender(); return; }   // asked for again if needed
        t.state = "error"; t.failedAt = performance.now();
        if (t.used === frame) {                          // a tile of the view: say so, and ask again
          console.error("raw tile", key, err); requestRender();
          setTimeout(requestRender, RETRY_MS + 50);
        }
      });
    }
    t.keep = frame;
    if (!ahead) t.used = frame;
    return t.state === "ready" ? t.tex : null;
  }
  const inView = (key) => { const t = tiles.get(key); return !!t && t.keep === frame; };
  function trim() {
    // give up downloads that have not started and that no tile wanted now needs any more;
    // drop the oldest queued work when panning has outrun the downloads or the workers, and
    // the oldest textures
    const wanted = new Set();
    for (const t of tiles.values()) if (t.keep === frame) for (const u of t.urls) wanted.add(u);
    for (let i = waiting.length - 1; i >= 0; i--)
      if (waiting[i].forView && !wanted.has(waiting[i].url)) waiting.splice(i, 1)[0].reject(new Error("dropped"));
    while (waiting.length > MAX_WAIT) waiting.shift().reject(new Error("dropped"));
    for (let i = queue.length - 1; i >= 0; i--)
      if (queue[i].tile && !inView(queue[i].tile)) queue.splice(i, 1)[0].reject(new Error("dropped"));
    while (queue.length > 96) queue.shift().reject(new Error("dropped"));
    for (const [k, t] of tiles) if (t.state === "error" && t.keep !== frame) tiles.delete(k);
    if (texBytes <= TEX_BYTES) return;
    const old = [...tiles].filter(([, t]) => t.keep !== frame && t.state === "ready")
      .sort((a, b) => a[1].keep - b[1].keep);
    for (const [k, t] of old) {
      if (texBytes <= TEX_BYTES) break;
      g.deleteTexture(t.tex); texBytes -= t.bytes;
      tiles.delete(k);
    }
  }

  // Which level a view is drawn from, or null with the reason in info.hint. Display copy: every
  // sample (level 0) while a block of level 1 would be wider than 2 device pixels and the tiles
  // are few enough; otherwise level 1.
  const source = () => (dman && R.src === "disp" ? "disp" : "raw");
  function chooseLevel(span, PW, nSeg) {
    const per = (lv) => Math.ceil(span * man.fs / lv.spt - 1e-6) + 1;   // tiles a view can need per shank
    const narrow = (lv) => span * man.fs / (PW * lv.block) > 60;        // values per pixel column: the shader reads 64 at most
    const fits = (lv) => span <= lv.maxSpan + 1e-9 && nSeg * per(lv) <= lv.maxTiles && !narrow(lv);
    const lvs = source() === "disp" ? DLV.slice() : [RAWLV];
    if (lvs.length === 2 && span * man.fs / PW >= lvs[1].block / 2) lvs.shift();
    for (const lv of lvs) if (fits(lv)) return lv;
    const last = lvs[lvs.length - 1];
    if (span > last.maxSpan + 1e-9) info.hint = `raw voltage: zoom in to ${last.maxSpan} s or less`;
    else if (narrow(last)) info.hint = "raw voltage: zoom in to a shorter window (the plot is narrow)";
    else {
      const ns = Math.floor(last.maxTiles / per(last));
      info.hint = `raw voltage: zoom in to ${ns} shank${ns === 1 ? "" : "s"} or fewer` + (per(last) > 2 ? ", or to a shorter window" : "");
    }
    if (source() === "raw" && DLV.some(fits)) info.hint += ", or press V for the filtered copy";
    return null;
  }

  // ------------------------------------------------------------ drawing
  // Called by render() before the spikes are drawn. Returns null, or the region (in cells of
  // the main canvas, fragment coordinates) under which voltage is drawn, with the opacity the
  // spike cells should have there.
  function draw() {
    info = { active: false, hint: "", loading: 0, failed: 0 };
    if (off) return null;
    const wide = state.t1 - state.t0 > HINT_SPAN_S;
    if (!man) {                          // the index did not load: say so when it matters, look again
      if (state.ymode === "depth" && !wide) info.hint = "raw voltage: the data could not be reached, trying again";
      if (!manTry.busy && performance.now() - manTry.at > MAN_RETRY_MS) {
        manTry.busy = true;
        load().finally(() => { manTry.busy = false; requestRender(); });
      }
      return null;
    }
    if (!dman && dispBase && !dispTry.busy && performance.now() - dispTry.at > dispTry.wait)
      findDisplay().finally(requestRender);
    if (!cv) {
      try { initGLRaw(); }
      catch (err) { console.warn("raw underlay switched off:", err); off = true; cv = null; return null; }
    }
    frame++;                           // tiles not asked for in this call are no longer wanted
    const hide = () => { cv.style.display = "none"; loadNow = null; trim(); return null; };
    if (!R.on) return hide();
    if (state.ymode !== "depth") {
      if (!wide) info.hint = 'raw voltage: set Neurons to "Plotted by depth"';
      return hide();
    }
    const span = state.t1 - state.t0, L = LAY.depth, nspan = state.n1 - state.n0;
    const segs = L.segs.filter((s) => s.y1 > state.n0 && s.y0 < state.n1 && man.probes[s.probe] &&
                                      man.probes[s.probe].shanks[s.shank]);
    if (!segs.length) return hide();
    const b = state.block, PW = layout.W * b, PH = layout.H * b;
    const lv = chooseLevel(span, PW, segs.length);
    if (!lv) return hide();
    shown = lv;
    maxFetch = lv.maxFetch;
    if (cv.width !== PW || cv.height !== PH) { cv.width = PW; cv.height = PH; }
    Object.assign(cv.style, { display: "block", left: layout.x + "px", top: layout.y + "px",
                              width: layout.w + "px", height: layout.h + "px" });
    g.viewport(0, 0, PW, PH);
    g.disable(g.SCISSOR_TEST);
    g.clearColor(0.5, 0.5, 0.5, 1); g.clear(g.COLOR_BUFFER_BIT);
    g.useProgram(P.p); g.bindVertexArray(vao);
    g.uniform1i(P.u.u_a, 0); g.uniform1i(P.u.u_b, 1); g.uniform1i(P.u.u_bounds, 2);
    g.uniform1i(P.u.u_PH, PH); g.uniform1i(P.u.u_n, lv.n); g.uniform1i(P.u.u_k, lv.block);
    g.uniform1f(P.u.u_scale, lv.scale);
    g.uniform1f(P.u.u_range, R.uv);
    const ymax = (man.n_channels / 2 - 1) * man.row_pitch_um;      // height of the top site
    g.uniform1f(P.u.u_ytop, ymax + man.row_pitch_um / 2);
    g.enable(g.SCISSOR_TEST);
    const sEdge = new Float64Array(PW + 1), rel = new Float32Array(PW + 1);
    let lastProbe = -1, T0 = 0, T1 = -1, base = 0;
    const drawnSegs = [];
    const colOf = (s) => {             // first column whose left edge is at sample s or later
      let a = 0, c = PW;
      while (a < c) { const m = (a + c) >> 1; if (sEdge[m] < s) a = m + 1; else c = m; }
      return a;
    };
    const count = (lvl, p, sh, T) => {   // a tile that is not there yet: loading, or failed
      const t = tiles.get(tileKey(lvl, p, sh, T));
      if (t && t.state === "error") info.failed++; else info.loading++;
    };
    for (const sg of segs) {
      const q = probe(sg.probe), pr = man.probes[sg.probe];
      if (q.state !== "ready") { if (q.state === "error") info.failed++; else info.loading++; continue; }
      if (sg.probe !== lastProbe) {
        // sample at every pixel-column edge on this probe's clock, relative to the first tile needed
        lastProbe = sg.probe;
        for (let x = 0; x <= PW; x++) sEdge[x] = interp(q.t, q.s, state.t0 + span * x / PW);
        // the tiles that hold the first sample of some column (a column's first sample is
        // floor(left edge + 0.5))
        T0 = Math.max(lv.range(sg.probe).tile_first, Math.floor(Math.floor(sEdge[0] + 0.5) / lv.spt));
        T1 = Math.min(lv.range(sg.probe).tile_last, Math.floor(Math.floor(sEdge[PW - 1] + 0.5) / lv.spt));
        base = T0 * lv.spt;
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
      drawnSegs.push([sg, T0, T1]);
      // half-site rows: a site at height y starts (top10 - 75) / 10 - y um into the segment,
      // like a neuron at that height (buildLayouts); the top site is half-row 0
      const off = sg.y0 + (sg.top10 - 75) / 10 - ymax;
      g.uniform1f(P.u.u_rho0, (state.n0 - off) / 7.5);
      g.uniform1f(P.u.u_drho, nspan / PH / 7.5);
      for (let T = T0; T <= T1; T++) {
        // the columns whose first sample is in this tile; the last of them may reach into the next
        const x0 = colOf(T * lv.spt - 0.5), x1 = colOf((T + 1) * lv.spt - 0.5);
        if (x1 <= x0) continue;
        const tex = wantTile(lv, sg.probe, sg.shank, T);
        const across = T < lv.range(sg.probe).tile_last && Math.ceil(sEdge[x1] + 0.5) > (T + 1) * lv.spt;
        const texB = across ? wantTile(lv, sg.probe, sg.shank, T + 1) : null;
        if (across && !texB && T === T1) count(lv, sg.probe, sg.shank, T + 1);   // otherwise counted in its own turn
        if (!tex) { count(lv, sg.probe, sg.shank, T); continue; }
        g.scissor(x0, PH - p1, x1 - x0, p1 - p0);
        g.activeTexture(g.TEXTURE0); g.bindTexture(g.TEXTURE_2D, tex);
        g.activeTexture(g.TEXTURE1); g.bindTexture(g.TEXTURE_2D, texB || tex);
        g.uniform1i(P.u.u_hasB, texB ? 1 : 0);
        g.uniform1f(P.u.u_tileOff, T * lv.spt - base);
        const m = motionOf(q, sg.shank, T * lv.spt);
        g.uniform1fv(P.u.u_mot, R.motion ? m : new Float32Array(m.length));
        g.drawArrays(g.TRIANGLES, 0, 3);
      }
    }
    g.disable(g.SCISSOR_TEST);
    // how long the view took to complete, and what it downloaded
    const now = performance.now();
    if (info.failed) loadNow = null;
    else if (info.loading && !loadNow) loadNow = { t0: now, files: 0, bytes: 0 };
    if (!info.loading && !info.failed && loadNow) { loadLast = { ms: now - loadNow.t0, files: loadNow.files, bytes: loadNow.bytes }; loadNow = null; }
    // the view is complete: fetch the tiles next to it in time, as many on each side as it shows
    if (!info.loading && !info.failed) {
      const most = Math.max(1, Math.floor(lv.maxTiles / (2 * drawnSegs.length)));
      for (const [sg, a, c] of drawnSegs) {
        const n = Math.min(c - a + 1, most), rg = lv.range(sg.probe);
        for (let i = 1; i <= n; i++) {
          if (c + i <= rg.tile_last) wantTile(lv, sg.probe, sg.shank, c + i, true);
          if (a - i >= rg.tile_first) wantTile(lv, sg.probe, sg.shank, a - i, true);
        }
      }
    }
    if (info.failed && !retryTimer)      // what failed is asked for again by a later drawing
      retryTimer = setTimeout(() => { retryTimer = null; requestRender(); }, RETRY_MS + 50);
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
    if (!man || off || state.ymode !== "depth") return false;
    // V whenever the layout allows the voltage (the other source may fit a view this one does not)
    if ((k === "v" || k === "V") && dman && R.on) { R.src = R.src === "disp" ? "raw" : "disp"; return true; }
    // the others only while the voltage is shown; U also brings it back after it was hidden with U
    if (!info.active && !((k === "u" || k === "U") && !R.on)) return false;
    if (k === "u" || k === "U") R.on = !R.on;
    else if (k === "o" || k === "O") R.alpha = (R.alpha + 1) % ALPHAS.length;
    else if (k === "m" || k === "M") R.motion = !R.motion;
    else if (k === "[") R.uv = Math.max(6.25, R.uv / 2);
    else if (k === "]") R.uv = Math.min(3200, R.uv * 2);
    else return false;
    return true;
  }
  function variantText() {
    if (shown && shown.src === "disp")
      return `filtered copy in steps of ${shown.step} µV` + (shown.block > 1 ? `, largest of every ${shown.block} samples` : "");
    const v = man.variant, step = v.scale * man.probes[0].shanks[0].uv_per_count;
    return `recorded samples in steps of ${step.toFixed(0)} µV, filtered in the browser`;
  }
  function loadText() {
    if (info.loading) return ` · loading ${info.loading}`;
    if (!loadLast) return "";
    const mb = loadLast.bytes / 1e6;
    return ` · this view: ${loadLast.files} file${loadLast.files === 1 ? "" : "s"}, ` +
      `${mb < 10 ? mb.toFixed(1) : mb.toFixed(0)} MB, ${(loadLast.ms / 1000).toFixed(1)} s`;
  }
  function drawLegend(ctx) {
    if (!man || !info.active) return;
    const { x: px, y: py } = layout;
    const lines = [`Raw voltage (${variantText()})` + (R.motion ? "" : " · motion shift off"),
                   `Black −${R.uv} µV, white +${R.uv} µV · spikes drawn at ${Math.round(100 * ALPHAS[R.alpha])}% opacity` +
                   loadText() + (info.failed ? ` · ${info.failed} failed to load` : ""),
                   "Keys: [ ] range · O spikes · M motion shift · " + (dman ? "V source · " : "") + "U hide"];
    ctx.font = "12px Arial"; ctx.textAlign = "left"; ctx.textBaseline = "top";
    const w = Math.max(...lines.map((s) => ctx.measureText(s).width)) + 16;
    ctx.fillStyle = "rgba(14,17,21,0.82)"; ctx.fillRect(px + 8, py + 8, w, 16 * lines.length + 10);
    ctx.fillStyle = "#e6e9ed";
    lines.forEach((s, i) => ctx.fillText(s, px + 16, py + 14 + 16 * i));
  }
  const status = () => info.active ? ` · raw voltage ±${R.uv} µV` : info.hint ? ` · ${info.hint}` : "";
  // Only what differs from the defaults goes into the URL.
  function toHash(p) {
    if (!R.on) p.set("raw", "off");
    if (R.uv !== DEF.uv) p.set("rawuv", R.uv);
    if (R.alpha !== DEF.alpha) p.set("rawsp", R.alpha);
    if (!R.motion) p.set("rawmotion", 0);
    if (R.src !== DEF.src) p.set("rawsrc", R.src);
  }
  function fromHash(p) {
    R.on = p.get("raw") !== "off";
    const u = parseFloat(p.get("rawuv")); R.uv = u >= 5 && u <= 3200 ? u : DEF.uv;
    const a = parseInt(p.get("rawsp")); R.alpha = a >= 0 && a < ALPHAS.length ? a : DEF.alpha;
    R.motion = p.get("rawmotion") !== "0";
    R.src = p.get("rawsrc") === "raw" ? "raw" : DEF.src;
  }

  // ------------------------------------------------------------ self-test (?rawtest=1)
  // Clock and motion files of every probe load and pass their checks. Recorded samples: the sum
  // of one chunk of every shank against the manifest (decoding), and the tiles listed in the
  // manifest against make_raw_proto.py's process_tile (processing). Display copy: the sum of
  // one tile of every shank and level against its manifest; level 1 against level 0; level 0
  // against the recorded samples processed here (the two differ by the rounding of both
  // copies, so this catches a shift in time or site, not small errors). Then the drawing of a
  // few views from each source.
  async function test() {
    if (!man) return { error: "no raw data" };
    const out = { base: man.base, display: dman ? dman.base : null, decode: [], process: [], probes: [] };
    const keep = { ...R };
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
      const d = await tileData(RAWLV, tt.probe, tt.shank, tt.tile);
      const ref = new Float32Array(await (await fetch(man.base + tt.file)).arrayBuffer());
      let maxAbs = 0, maxRef = 0;
      for (let i = 0; i < ref.length; i++) {
        maxAbs = Math.max(maxAbs, Math.abs(ref[i] - d.data[i])); maxRef = Math.max(maxRef, Math.abs(ref[i]));
      }
      out.process.push({ probe: tt.probe, shank: tt.shank, tile: tt.tile, n: ref.length,
                         sameLength: ref.length === d.data.length, maxAbsDiff_uV: maxAbs, maxRef_uV: maxRef, ms: d.ms });
    }
    const drawOk = (q) => !q.error && q.fracDiffer < 0.002 && q.drawn > 0.5 * q.pixels && q.cellsWrong === 0 && q.aligned;
    R.src = "raw";
    out.draw = await testViews(false);
    out.ok = out.probes.every((q) => q.ok) && out.decode.every((q) => q.ok) && out.process.length > 0 &&
      out.process.every((q) => q.sameLength && q.maxAbsDiff_uV < 0.01) && out.draw.every(drawOk);
    if (dman) {
      try { await testDisplay(out); }
      catch (err) { out.displayError = String(err.message || err); }
      R.src = "disp";
      out.ddraw = await testViews(true);
      out.ok = out.ok && !out.displayError && out.ddecode.every((q) => q.ok) && out.dcross.length > 0 &&
        out.dcross.every((q) => q.rms_uV < 6 && q.corr > 0.9 && q.level1Wrong === 0) &&
        out.ddraw.every(drawOk) && out.ddraw.some((q) => q.level === "d0") && out.ddraw.some((q) => q.level === "d1");
    }
    Object.assign(R, keep);
    render();
    return out;
  }
  async function testDisplay(out) {
    out.ddecode = []; out.dcross = [];
    for (const pr of dman.probes) for (const s of pr.shanks) for (const lv of DLV) {
      const ck = (s.check || []).find((c) => "d" + c.level === lv.id), r = { probe: pr.probe, shank: s.shank, level: lv.id };
      if (!ck) { out.ddecode.push({ ...r, ok: false, missing: true }); continue; }
      try {
        const d = (await tileData(lv, pr.probe, s.shank, ck.tile)).data;
        let sum = 0;
        for (let i = 0; i < d.length; i++) sum += d[i];
        out.ddecode.push({ ...r, ok: sum === ck.sum && d.length === lv.n * man.n_channels });
      } catch (err) { out.ddecode.push({ ...r, ok: false, error: String(err.message || err) }); }
    }
    const [l0, l1] = DLV, nch = man.n_channels;
    if (man.tile % l0.spt || l1.spt % l0.spt) throw new Error("tile lengths do not nest");
    for (const tt of man.test || []) {
      // level 0 against the recorded samples processed in the browser
      const ref = (await tileData(RAWLV, tt.probe, tt.shank, tt.tile)).data, per = man.tile / l0.spt, d0 = [];
      for (let i = 0; i < per; i++) d0.push((await tileData(l0, tt.probe, tt.shank, tt.tile * per + i)).data);
      let see = 0, maxAbs = 0, sxy = 0, sxx = 0, syy = 0;
      for (let c = 0; c < nch; c++) for (let i = 0; i < man.tile; i++) {
        const v = d0[Math.floor(i / l0.n)][c * l0.n + i % l0.n] * l0.step, r = ref[c * man.tile + i], e = v - r;
        see += e * e; maxAbs = Math.max(maxAbs, Math.abs(e)); sxy += v * r; sxx += v * v; syy += r * r;
      }
      // level 1 against level 0: each value is a sample of its block, of the largest magnitude there
      const T1 = Math.floor(tt.tile * man.tile / l1.spt), n0 = l1.spt / l0.spt, a0 = [];
      const a1 = (await tileData(l1, tt.probe, tt.shank, T1)).data;
      for (let i = 0; i < n0; i++) a0.push((await tileData(l0, tt.probe, tt.shank, T1 * n0 + i)).data);
      let wrong = 0;
      for (let c = 0; c < nch; c++) for (let j = 0; j < l1.n; j++) {
        const got = a1[c * l1.n + j];
        let big = 0, found = false;
        for (let k = 0; k < l1.block; k++) {
          const smp = j * l1.block + k, v = a0[Math.floor(smp / l0.n)][c * l0.n + smp % l0.n];
          big = Math.max(big, Math.abs(v)); found = found || v === got;
        }
        if (!found || Math.abs(got) !== big) wrong++;
      }
      out.dcross.push({ probe: tt.probe, shank: tt.shank, tile: tt.tile, rms_uV: Math.sqrt(see / (nch * man.tile)),
                        maxAbsDiff_uV: maxAbs, corr: sxy / Math.sqrt(sxx * syy), level1Tile: T1, level1Wrong: wrong });
    }
  }

  // The drawing of the current view, checked against the rule computed here in double
  // precision: every device pixel of the raw canvas (see the head of this file for the rule).
  // Samples are looked up by their number in the recording, whatever tile holds them, so this
  // does not repeat the drawing's bookkeeping of tiles. The GPU works in single precision (and
  // keeps recorded samples as half floats), so a gray level may differ by 1, and a few pixels
  // whose edge falls on a sample or site boundary may differ by more.
  async function checkDraw() {
    render();
    if (!info.active || info.loading || info.failed) return { error: "the view is not drawn completely" };
    const b = state.block, PW = layout.W * b, PH = layout.H * b;
    const lv = shown, spt = lv.spt, n = lv.n, blk = lv.block, mult = lv.src === "disp" ? lv.step : 1;
    const px = new Uint8Array(PW * PH * 4);
    render();                                           // draw and read in the same task
    g.readPixels(0, 0, PW, PH, g.RGBA, g.UNSIGNED_BYTE, px);
    const span = state.t1 - state.t0, nspan = state.n1 - state.n0;
    const ymax = (man.n_channels / 2 - 1) * man.row_pitch_um;
    let pixels = 0, drawn = 0, differ = 0, maxDiff = 0;
    for (const sg of LAY.depth.segs) {
      if (!(sg.y1 > state.n0 && sg.y0 < state.n1)) continue;
      const q = probes.get(sg.probe), pr = man.probes[sg.probe], mo = pr.motion, nd = mo.depths_um.length;
      const rg = lv.range(sg.probe);
      const sE = new Float64Array(PW + 1);
      for (let x = 0; x <= PW; x++) sE[x] = interp(q.t, q.s, state.t0 + span * x / PW);
      const data = new Map();
      for (let T = Math.floor((sE[0] - 0.5) / spt); T <= Math.floor((sE[PW] + 0.5) / spt); T++)
        if (T >= rg.tile_first && T <= rg.tile_last) data.set(T, (await tileData(lv, sg.probe, sg.shank, T)).data);
      const value = (j) => {             // value number j of the recording at this level, or null
        const d = data.get(Math.floor(j / n));
        return d ? d : null;
      };
      const off = sg.y0 + (sg.top10 - 75) / 10 - ymax, drho = nspan / PH / 7.5;
      for (let p = 0; p < PH; p++) {
        const nc = state.n0 + (p + 0.5) * nspan / PH;
        if (!(nc >= sg.y0 && nc < sg.y1)) continue;
        const ra0 = (state.n0 + p * nspan / PH - off) / 7.5;
        const yc = ymax + 7.5 - (ra0 + 0.5 * drho) * 7.5;
        const u = Math.min(nd - 1, Math.max(0, (yc - mo.depths_um[0]) / (mo.depths_um[1] - mo.depths_um[0])));
        const k = Math.min(Math.floor(u), nd - 2);
        for (let x = 0; x < PW; x++) {
          // the samples whose 1-sample band overlaps the column, and the values that hold them
          const i0 = Math.floor(sE[x] + 0.5), i1 = Math.ceil(sE[x + 1] + 0.5);
          const j0 = Math.floor(i0 / blk), j1 = Math.min(Math.ceil(i1 / blk), j0 + 64);
          let want = 128;                               // nothing drawn: the canvas's mid gray
          if (value(j0)) {                              // no data where the column starts: not drawn
            const m = motionOf(q, sg.shank, i0), shift = R.motion ? m[k] + (m[k + 1] - m[k]) * (u - k) : 0;
            const ra = ra0 - shift / 7.5;
            const r0 = Math.max(0, Math.floor(ra)), r1 = Math.min(384, Math.ceil(ra + drho), r0 + 24);
            if (j1 > j0 && r1 > r0) {
              let best = 0;
              for (let r = r0; r < r1; r++) {
                const c = (191 - (r >> 1)) * 2 + (r & 1);
                for (let j = j0; j < j1; j++) {
                  const d = value(j);
                  if (!d) break;
                  const v = d[c * n + j % n];
                  if (Math.abs(v) > Math.abs(best)) best = v;
                }
              }
              want = Math.round(255 * Math.min(1, Math.max(0, 0.5 + 0.5 * best * mult / R.uv)));
              drawn++;
            }
          }
          const got = px[((PH - 1 - p) * PW + x) * 4], e = Math.abs(got - want);
          pixels++;
          if (e > 1) { differ++; maxDiff = Math.max(maxDiff, e); }
        }
      }
    }
    // The spike cells on the main canvas: opaque where a cell holds spikes, transparent elsewhere
    // in the shanks shown, untouched outside them. And the two canvases cover the same rectangle.
    const W = layout.W, H = layout.H, v0 = readValues(W, H, 0), v2 = readValues(W, H, 2);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    const mp = new Uint8Array(W * H * 4);
    gl.readPixels(0, 0, W, H, gl.RGBA, gl.UNSIGNED_BYTE, mp);
    const segsShown = LAY.depth.segs.filter((sg) => sg.y1 > state.n0 && sg.y0 < state.n1);
    let cellsWrong = 0;
    for (let y = 0; y < H; y++) {
      const nc = state.n0 + (y + 0.5) * nspan / H, inside = segsShown.some((sg) => nc >= sg.y0 && nc < sg.y1);
      for (let x = 0; x < W; x++) {
        const spikes = v0[y * W + x] > 0 || (state.spikes === "overlay" && v2[y * W + x] > 0);
        const want = !inside ? 255 : spikes ? Math.round(255 * ALPHAS[R.alpha]) : 0;
        if (Math.abs(mp[((H - 1 - y) * W + x) * 4 + 3] - want) > 1) cellsWrong++;
      }
    }
    const ra = cvs.getBoundingClientRect(), rb = cv.getBoundingClientRect();
    const aligned = ["left", "top", "width", "height"].every((q) => Math.abs(ra[q] - rb[q]) < 0.01);
    return { level: lv.id, pixels, drawn, differ, fracDiffer: differ / Math.max(1, pixels), maxDiff, cells: W * H, cellsWrong, aligned };
  }
  // checkDraw on a few views: one shank; a shank that had moved by tens of micrometers; four
  // shanks across a tile edge and across a change of the motion estimate (a whole second of
  // probe time); two probes, each on its own clock; and, from the display copy, two longer
  // windows, which are drawn from level 1 (in the longest a pixel column covers about 60
  // samples and 7 tile edges are crossed). Times are found on the probe's clock.
  async function testViews(display) {
    const seg = (p, sh) => LAY.depth.segs.find((q) => q.probe === p && q.shank === sh);
    const timeAt = (p, smp) => { const q = probes.get(p); return interp(q.s, q.t, smp); };
    const views = [
      { name: "one shank, 33 ms", t0: 1125.037, t1: 1125.070, n0: seg(0, 0).y0 + 520, n1: seg(0, 0).y0 + 920 },
      { name: "moved shank, 80 ms", t0: 273.50, t1: 273.58, n0: seg(1, 0).y1 - 420, n1: seg(1, 0).y1 },
      { name: "four shanks across a second", t0: timeAt(4, 1749 * 30000) - 0.15, t1: timeAt(4, 1749 * 30000) + 0.15,
        n0: seg(4, 0).y0, n1: seg(4, 3).y1 },
      { name: "motion shift off", t0: 273.50, t1: 273.58, n0: seg(1, 0).y1 - 420, n1: seg(1, 0).y1, motion: false },
      { name: "two probes (two clocks)", t0: 1125.04, t1: 1125.10, n0: seg(0, 3).y1 - 300, n1: seg(1, 0).y0 + 300 },
    ];
    if (display) views.push({ name: "one probe, 1.8 s", t0: 1124.6, t1: 1126.4, n0: seg(0, 0).y0, n1: seg(0, 3).y1 },
                            { name: "one shank, 3.6 s", t0: 1120.2, t1: 1123.8, n0: seg(0, 1).y0, n1: seg(0, 1).y1 });
    const out = [], keep = { ...R };
    for (const v of views) {
      Object.assign(R, DEF, { motion: v.motion !== false, src: keep.src });
      state.ymode = "depth"; state.t0 = v.t0; state.t1 = v.t1; state.n0 = v.n0; state.n1 = v.n1;
      render();
      for (let i = 0; i < 600 && (!info.active || info.loading); i++) { await new Promise((r) => setTimeout(r, 100)); render(); }
      out.push({ view: v.name, ...(await checkDraw()) });
    }
    Object.assign(R, keep);
    render();
    return out;
  }

  return { load, draw, key, drawLegend, status, toHash, fromHash, test, checkDraw,
           get available() { return !!man; }, get info() { return info; }, get settings() { return { ...R }; },
           get manifest() { return man; }, get displayManifest() { return dman; },
           get level() { return info.active && shown ? shown.id : null; }, get lastLoad() { return loadLast; },
           get stats() {      // for checks: what is held and what is under way
             return { tiles: tiles.size, textures: [...tiles.values()].filter((t) => t.tex).length, files: files.size,
                      textureMB: texBytes / 1e6, waiting: waiting.length, fetching, queued: queue.length, jobs: jobs.size, frame,
                      inView: [...tiles.values()].filter((t) => t.used === frame).length,
                      ahead: [...tiles.values()].filter((t) => t.keep === frame && t.used !== frame).length };
           },
           // for checks: probe sample at a viewer time, and the motion estimate used for a tile
           sampleAt: (p, t) => { const q = probes.get(p); return q && q.state === "ready" ? interp(q.t, q.s, t) : NaN; },
           motionAt: (p, shank, T) => { const q = probes.get(p); return q && q.state === "ready" ? Array.from(motionOf(q, shank, T * man.tile)) : null; } };
})();
