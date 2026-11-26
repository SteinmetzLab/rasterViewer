// Sonification for the raster viewer: the view scrolls under a playhead line and what crosses
// the line is played as sound. Uses viewer.js's globals (state, layout, gl, res, ...).
//
// What you hear is what you see:
//  - Time. The raster moves left at a fixed speed on the screen: one plot width in `cross`
//    seconds (default 30). So the playback rate is (window duration) / cross seconds of
//    recording per second, and it changes with the zoom. The view advances by whole columns of
//    cells on a fixed grid, so a column holds the same spikes at every step; while a column
//    crosses the line, its cells set the sound.
//  - Pitch. The vertical axis of the current view maps to pitch, high at the top, over a range
//    of `octaves` centered on A4. A cell row sounds the note of the selected scale nearest to
//    its pitch ("No rounding": up to 64 equal bands). When a neuron is taller than a cell, all of
//    its cells sound the one note of the neuron's middle, so a neuron is one note.
//  - Loudness. A cell's level follows its color: with t its position on the color scale (0 to 1,
//    the same formula as the color map), level = 10^(-36 dB x (1 - t) / 20), and silence where
//    t = 0 (no spikes, or a rate at or below the bottom of the color range). A note's level is
//    the RMS of the levels of the cells that sound it.
//  - One mechanism at every zoom. Each note is an oscillator whose loudness follows the note's
//    level with a fast attack and a slower release, set by the instrument. Zoomed in, a spike is
//    one bright cell among dark ones: the level jumps up for one column and falls back, which is
//    a struck note that rings out. Zoomed out, every column has some rate: the level changes a
//    little from column to column and the note is a sustained tone whose loudness follows the
//    firing rate. In between, spikes arriving faster than the release fuse into a tremolo.
//
// Sound and picture share one clock (the audio clock): column g of the scrolling grid sounds
// from audio time aA + (g - PA) / rate on, and the picture is drawn for the time now heard
// (audio time minus the output latency the browser reports).
"use strict";

const Sonify = (() => {
  const CENTER_MIDI = 69;    // A4 (440 Hz): the middle of the pitch range
  const DR_DB = 36;          // the bottom of the color scale is this far below the top, in dB
  const LOOK_S = 0.3;        // sound is scheduled this far ahead of the playhead
  const LEAD_S = 0.08;       // from pressing play to the first sound
  const MAX_BANDS = 64;      // voices when pitch is not rounded to a scale
  const NOTE_NAMES = ["C", "C♯ / D♭", "D", "D♯ / E♭", "E", "F", "F♯ / G♭", "G", "G♯ / A♭", "A", "A♯ / B♭", "B"];
  const SCALES = {           // pitch classes relative to the key
    pentatonic: { label: "Major pentatonic", pcs: [0, 2, 4, 7, 9] },
    minor_pentatonic: { label: "Minor pentatonic", pcs: [0, 3, 5, 7, 10] },
    major: { label: "Major", pcs: [0, 2, 4, 5, 7, 9, 11] },
    minor: { label: "Natural minor", pcs: [0, 2, 3, 5, 7, 8, 10] },
    whole: { label: "Whole tone", pcs: [0, 2, 4, 6, 8, 10] },
    chromatic: { label: "12-tone (chromatic)", pcs: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] },
    none: { label: "No rounding (continuous)", pcs: null },
  };
  // harm: amplitudes of the harmonics 1, 2, 3, ...; attack and release: time constants (s) with
  // which a note's loudness rises to and falls back from its level
  const INSTRUMENTS = {
    marimba: { label: "Marimba", harm: [1, 0, 0, 0.3], attack: 0.004, release: 0.16 },
    harp: { label: "Harp", harm: [1, 0.5, 0.33, 0.2, 0.12, 0.08], attack: 0.004, release: 0.45 },
    organ: { label: "Organ", harm: [1, 0.6, 0.3, 0.25, 0, 0.1], attack: 0.03, release: 0.09 },
    sine: { label: "Pure tone", harm: [1], attack: 0.006, release: 0.12 },
  };
  const cfg = { key: 0, scale: "pentatonic", inst: "marimba", octaves: 4, cross: 30, vol: 0.7 };

  let G = null;              // audio graph: { ctx, master, voices, attack, release, key }
  let playing = false, A = null, raf = 0, suspendTimer = 0;
  let M = null;              // cell row -> voice map of the current view
  let shiftCss = 0;          // the canvas is drawn this far left (CSS px): the part of a column already passed
  let head = null;           // { x (CSS px in the stage), t (s) }: the playhead
  let drawn = { x: -1, shift: -1 };   // what the overlay last showed
  let rendered = "";         // viewKey() of the view the value texture holds
  let lat = -1;              // output latency (s), smoothed

  // ------------------------------------------------------------ pitch
  const midiHz = (m) => 440 * Math.pow(2, (m - 69) / 12);
  const pitchRange = (octaves) => [CENTER_MIDI - 6 * octaves, CENTER_MIDI + 6 * octaves];

  // The voices' pitches (MIDI numbers, ascending) for settings c and a grid of H cell rows.
  function voicePitches(c, H) {
    const [lo, hi] = pitchRange(c.octaves);
    if (c.scale === "none") {
      const n = Math.max(1, Math.min(H, MAX_BANDS));
      return Array.from({ length: n }, (_, i) => lo + (i + 0.5) / n * (hi - lo));
    }
    const pcs = SCALES[c.scale].pcs, out = [];
    for (let m = lo; m <= hi; m++) if (pcs.includes((((m - c.key) % 12) + 12) % 12)) out.push(m);
    return out;
  }

  function nearest(pitches, p) {      // index of the pitch nearest to p (a tie goes to the lower)
    let a = 0, b = pitches.length - 1;
    while (a < b) { const m = (a + b) >> 1; if (pitches[m] < p) a = m + 1; else b = m; }
    return a > 0 && p - pitches[a - 1] <= pitches[a] - p ? a - 1 : a;
  }

  // Where on the vertical axis cell row k sounds: the cell's middle, or, when neurons are taller
  // than cells, the middle of the neuron under the cell's middle.
  function soundY(k, H, n0, n1, L) {
    const dn = (n1 - n0) / H;
    let y = n0 + (k + 0.5) * dn;
    if (L.frac) { if (dn < 1) y = Math.floor(y) + 0.5; return y; }
    if (dn >= L.unit) return y;
    // plotted by depth: of the units whose band holds y (same segment), the one centered nearest
    const seg = L.segs[segAt(L, y)];
    let a = seg.r0, b = seg.r1;
    while (a < b) { const m = (a + b) >> 1; if (L.lo[m] + L.unit > y) b = m; else a = m + 1; }
    let best = null;
    for (let i = a; i < seg.r1 && L.lo[i] <= y; i++) {
      const c = L.lo[i] + L.unit / 2;
      if (best === null || Math.abs(c - y) < Math.abs(best - y)) best = c;
    }
    return best === null ? y : best;
  }

  // voice[k]: the voice cell row k sounds (k = 0 at the top, highest pitch); count[v]: its rows.
  function rowVoices(pitches, c, H, n0, n1, L) {
    const [lo, hi] = pitchRange(c.octaves), voice = new Int16Array(H), count = new Int32Array(pitches.length);
    for (let k = 0; k < H; k++) {
      const f = Math.min(1, Math.max(0, (soundY(k, H, n0, n1, L) - n0) / (n1 - n0)));
      voice[k] = nearest(pitches, hi - f * (hi - lo));
      count[voice[k]]++;
    }
    return { voice, count };
  }

  // ------------------------------------------------------------ levels
  // position on the color scale, as in the color-map shader
  const tOf = (v, lmin, lmax) => v > 0 ? Math.min(1, Math.max(0, (Math.log(v) - lmin) / (lmax - lmin))) : 0;
  const levelOf = (t) => t > 0 ? Math.pow(10, -DR_DB * (1 - t) / 20) : 0;

  // Level of every voice for column c of px: RGBA float values read from the value texture, n
  // columns by H rows, bottom row first (channel 0 = value, 2 = second value of the overlay).
  function columnLevels(px, n, H, c, map, lmin, lmax, overlay, out) {
    out.fill(0);
    for (let y = 0; y < H; y++) {
      const o = 4 * (y * n + c);
      let t = tOf(px[o], lmin, lmax);
      if (overlay) t = Math.max(t, tOf(px[o + 2], lmin, lmax));
      if (t > 0) { const a = levelOf(t); out[map.voice[H - 1 - y]] += a * a; }
    }
    for (let v = 0; v < out.length; v++) out[v] = map.count[v] ? Math.sqrt(out[v] / map.count[v]) : 0;
    return out;
  }

  function readColumns(j, n, H) {     // values of grid columns j .. j + n - 1 as drawn
    const px = new Float32Array(n * H * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, res.vals.fbo);
    gl.readPixels(j, 0, n, H, gl.RGBA, gl.FLOAT, px);
    return px;
  }

  // ------------------------------------------------------------ audio graph
  // ctx: an AudioContext, or an OfflineAudioContext in the self-test (then without the limiter)
  function buildGraph(ctx, limiter = true) {
    const master = ctx.createGain();
    master.gain.value = 0;
    if (limiter) {           // catches peaks when many notes are loud at once (browsers add a fixed make-up gain to it)
      const c = ctx.createDynamicsCompressor();
      c.threshold.value = -6; c.knee.value = 6; c.ratio.value = 20; c.attack.value = 0.002; c.release.value = 0.1;
      master.connect(c); c.connect(ctx.destination);
    } else master.connect(ctx.destination);
    return { ctx, master, voices: [], attack: 0.005, release: 0.15, key: "" };
  }

  function setVoices(g, pitches, inst) {
    const ctx = g.ctx, now = ctx.currentTime;
    for (const v of g.voices) {      // let the old notes die away
      v.gain.gain.cancelScheduledValues(now); v.gain.gain.setTargetAtTime(0, now, 0.03);
      v.osc.stop(now + 0.4);
    }
    const real = new Float32Array(inst.harm.length + 1), imag = new Float32Array(inst.harm.length + 1);
    imag.set(inst.harm, 1);
    const wave = ctx.createPeriodicWave(real, imag);     // normalized to a peak of 1
    g.attack = inst.attack; g.release = inst.release;
    g.voices = pitches.map((m, i) => {
      const osc = ctx.createOscillator(), gain = ctx.createGain(), hz = midiHz(m);
      osc.setPeriodicWave(wave); osc.frequency.value = hz;
      gain.gain.value = 0;
      osc.connect(gain); gain.connect(g.master);
      // each starts at a different point of its cycle, so the notes' peaks do not all coincide
      osc.start(now + ((i * 0.618034) % 1) / hz);
      // low notes a little louder, high ones a little softer: equal levels then sound more equal
      const tilt = Math.min(1.6, Math.max(0.6, Math.pow(440 / hz, 0.3)));
      return { osc, gain, hz, midi: m, tilt, last: 0 };
    });
  }

  // One column: each voice's loudness heads for its level from audio time T on.
  function scheduleColumn(g, levels, T) {
    for (let v = 0; v < g.voices.length; v++) {
      const V = g.voices[v], a = levels[v] * V.tilt;
      if (Math.abs(a - V.last) < 1e-4) continue;         // unchanged; last is NaN after a restart
      const rising = V.last === V.last ? a > V.last : a > 0;
      V.gain.gain.setTargetAtTime(a, T, rising ? g.attack : g.release);
      V.last = a;
    }
  }

  const masterLevel = (g) => cfg.vol * 0.9 / Math.sqrt(Math.max(1, g.voices.length));

  // ------------------------------------------------------------ playing
  function ensureVoices() {
    const H = layout.H, L = lay();
    const vkey = [cfg.key, cfg.scale, cfg.inst, cfg.octaves, cfg.scale === "none" ? H : 0].join();
    if (G.key !== vkey) {
      setVoices(G, voicePitches(cfg, H), INSTRUMENTS[cfg.inst]); G.key = vkey; M = null;
      if (A) A.gNext = Math.floor(A.PA);      // the new notes take over from the column now sounding (see due)
    }
    const mkey = [vkey, H, state.n0, state.n1, state.ymode].join();
    if (!M || M.key !== mkey) {
      M = rowVoices(G.voices.map((v) => v.midi), cfg, H, state.n0, state.n1, L);
      M.key = mkey; M.levels = new Float64Array(G.voices.length);
    }
  }

  // ---- transport (pure functions of the anchor and a time; the self-test runs them on a fake clock)
  // An anchor fixes the grid of columns: column g covers ticks T0a + g * dT to T0a + (g + 1) * dT,
  // and the view shows columns k .. k + W - 1 (step k). The playhead is at column position PA at
  // audio time aA and moves at `rate` columns per second. endTick: the session's end.
  function makeAnchor(t0, t1, W, tHead, now, cross, endTick) {
    const s = t1 - t0, dT = s * FS / W, T0a = t0 * FS, PA = (tHead - t0) * FS / dT;
    return { W, s, dT, T0a, PA, aA: now + LEAD_S, rate: W / cross, cross, gNext: Math.floor(PA), k: 0,
             fresh: true, own: { t0, t1 },
             // the last step whose window still ends inside the session
             kmax: Math.max(0, Math.floor((endTick - T0a) / dT - W + 1e-6)) };
  }
  // Where the view and the line are at audio time `heard`. The view scrolls with the line at its
  // middle; where it cannot (the window touches the start or the end of the session), the view
  // stays and the line moves instead. col: the line's position in columns from the view's left
  // edge; frac: the part of a column the scrolling view has moved past its step.
  function position(A, heard) {
    const P = A.PA + Math.max(0, heard - A.aA) * A.rate;
    let k = Math.floor(P - A.W / 2), frac = P - A.W / 2 - k;
    if (k < 0) { k = 0; frac = 0; }
    if (k >= A.kmax) { k = A.kmax; frac = 0; }
    return { P, k, frac, col: P - k, done: P - k >= A.W };
  }
  const columnTime = (A, g) => A.aA + (g - A.PA) / A.rate;     // audio time at which column g starts to sound
  // Columns [g0, g1) to schedule at audio time `now`: those that start within LOOK_S, not yet
  // scheduled, inside the view at step A.k, and not already over (after a stall).
  function due(A, now) {
    return [Math.max(A.gNext, A.k, Math.floor(A.PA + (now - A.aA) * A.rate)),
            Math.min(A.k + A.W, Math.ceil(A.PA + (now + LOOK_S - A.aA) * A.rate))];
  }

  // (Re)start from playhead time tHead (s) with the window now shown.
  function anchor(tHead) {
    const now = G.ctx.currentTime;
    for (const v of G.voices) { v.gain.gain.cancelScheduledValues(now); v.last = NaN; }
    A = makeAnchor(state.t0, state.t1, layout.W, tHead, now, cfg.cross, TE * FS);
  }
  const mine = () => A && state.t0 === A.own.t0 && state.t1 === A.own.t1 && layout.W === A.W;
  // While playing, the cells are counted on the anchor's grid shifted by a whole number of
  // columns, so a column holds exactly the same ticks at every step (viewer.js viewGrid).
  const grid = () => (playing && mine() && !A.fresh ? { T0: A.T0a, dT: A.dT, j0: A.k } : null);

  // The raster and the raw-voltage canvas under it (raw_underlay.js) are drawn shiftCss to the left.
  function shiftCanvases() {
    for (const c of [cvs, document.getElementById("rawcv")]) {
      if (!c) continue;
      c.style.transform = shiftCss ? `translateX(${-shiftCss}px)` : "";
      c.style.clipPath = shiftCss ? `inset(0 0 0 ${shiftCss}px)` : "";
    }
  }

  function tick() {
    if (!playing) return;
    raf = requestAnimationFrame(tick);
    // the user moved or zoomed the view, or the grid changed: go on from the middle of the new view
    if (!mine() || cfg.cross !== A.cross) anchor(0.5 * (state.t0 + state.t1));
    // the picture shows what is heard now: audio time less the output latency (smoothed, so a
    // changing estimate does not make the picture jump)
    const l = G.ctx.outputLatency || G.ctx.baseLatency || 0;
    lat = lat < 0 ? l : lat + 0.05 * (l - lat);
    const p = position(A, G.ctx.currentTime - lat);
    if (p.done) { stop(); return; }                       // the line reached the end of the session
    const cellCss = state.block / layout.dpr;
    shiftCss = Math.floor(p.frac * state.block) / layout.dpr;   // whole device pixels
    head = { x: layout.x + (p.col - p.frac) * cellCss, t: (A.T0a + p.P * A.dT) / FS };
    shiftCanvases();
    if (p.k !== A.k || A.fresh) {
      A.k = p.k; A.fresh = false;
      state.t0 = (A.T0a + p.k * A.dT) / FS; state.t1 = state.t0 + A.s;
      clampView();
      A.own = { t0: state.t0, t1: state.t1 };
      render();
    } else if (rendered !== viewKey()) render();      // something else about the view changed
    else if (shiftCss !== drawn.shift || Math.abs(head.x - drawn.x) > 0.3)
      { drawOverlay(); Traces.draw(); }   // the line, the time ticks and the traces move between steps
    schedule();
  }

  // Schedule the sound of the columns that come under the line within the next LOOK_S, from the
  // cell values as last drawn. Runs every frame: also while the view stands still.
  function schedule() {
    if (!playing || !mine() || A.fresh || rendered !== viewKey()) return;
    ensureVoices();
    const now = G.ctx.currentTime, H = layout.H, s = state[state.mode];
    const [g0, g1] = due(A, now);
    if (g1 > g0) {
      const n = g1 - g0, px = readColumns(g0 - A.k, n, H);
      const lmin = Math.log(s.vmin), lmax = Math.log(s.vmax), over = state.spikes === "overlay";
      for (let c = 0; c < n; c++)
        scheduleColumn(G, columnLevels(px, n, H, c, M, lmin, lmax, over, M.levels), columnTime(A, g0 + c));
      A.gNext = g1;
    }
    // If scheduling stops (a stalled page), the sound fades out instead of holding its last notes.
    const end = Math.max(now, columnTime(A, A.gNext)), m = G.master.gain;
    m.cancelScheduledValues(now);
    m.setTargetAtTime(masterLevel(G), now, 0.03);
    m.setTargetAtTime(0, end + 0.15, 0.08);
  }

  // Called at the end of every render: note which view the value texture now holds.
  function afterRender() { rendered = viewKey(); }

  function start() {
    if (playing) return;
    try {
      if (!G) G = buildGraph(new (window.AudioContext || window.webkitAudioContext)({ latencyHint: "interactive" }));
    } catch (err) {
      document.getElementById("status").textContent = "Sonify: this browser cannot play sound (" + err.message + ")";
      return;
    }
    clearTimeout(suspendTimer);
    G.ctx.resume().catch((err) => {
      stop(); document.getElementById("status").textContent = "Sonify: the browser did not start the sound (" + err.message + ")";
    });
    playing = true; lat = -1;
    // from the middle of the view; from its left edge when the view starts at the session's start
    anchor(state.t0 <= TS ? state.t0 : 0.5 * (state.t0 + state.t1));
    syncUI();
    raf = requestAnimationFrame(tick);
  }

  function stop() {
    if (!playing) return;
    playing = false; cancelAnimationFrame(raf);
    const now = G.ctx.currentTime;
    G.master.gain.cancelScheduledValues(now); G.master.gain.setTargetAtTime(0, now, 0.03);
    for (const v of G.voices) { v.gain.gain.cancelScheduledValues(now); v.gain.gain.setTargetAtTime(0, now, 0.05); v.last = 0; }
    suspendTimer = setTimeout(() => { if (!playing) G.ctx.suspend(); }, 600);
    shiftCss = 0; head = null; A = null;
    shiftCanvases();
    syncUI(); requestRender();
  }
  const toggle = () => (playing ? stop() : start());

  // ------------------------------------------------------------ drawing, status, controls
  function drawHead(ctx) {              // the playhead line, on the plot and on the overview strip
    if (!playing || !head) return;
    drawn = { x: head.x, shift: shiftCss };
    const { y: py, h: ph, x: px, w: pw } = layout;
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.9)"; ctx.fillStyle = "#ffffff"; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(head.x, py - 5); ctx.lineTo(head.x, py + ph + 5); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(head.x - 5, py - 9); ctx.lineTo(head.x + 5, py - 9); ctx.lineTo(head.x, py - 2); ctx.closePath(); ctx.fill();
    const ox = px + (head.t - TS) / (TE - TS) * pw;
    ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(ox, OV_Y - 2); ctx.lineTo(ox, OV_Y + OV_H + 2); ctx.stroke();
    ctx.restore();
  }

  function status() {
    if (!playing) return "";
    const sp = (state.t1 - state.t0) / cfg.cross;
    const rate = sp >= 1 ? `${+sp.toPrecision(3)}× real time` : `1/${+(1 / sp).toPrecision(3)} of real time`;
    return ` · sonify: ${rate}, ${voicePitches(cfg, layout.H).length} notes`;
  }

  function syncUI() {
    const b = document.getElementById("sonify");
    if (!b) return;
    b.textContent = playing ? "Stop sound" : "Sonify";
    document.getElementById("sonWrap").style.display = playing ? "flex" : "none";
    document.getElementById("sonKey").value = String(cfg.key);
    document.getElementById("sonKey").disabled = cfg.scale === "chromatic" || cfg.scale === "none";
    document.getElementById("sonScale").value = cfg.scale;
    document.getElementById("sonInst").value = cfg.inst;
    document.getElementById("sonOct").value = String(cfg.octaves);
    document.getElementById("sonCross").value = String(cfg.cross);
    document.getElementById("sonVol").value = String(cfg.vol);
  }

  function setup() {
    const fill = (id, entries) => {
      const el = document.getElementById(id);
      for (const [value, label] of entries) {
        const o = document.createElement("option"); o.value = value; o.textContent = label; el.appendChild(o);
      }
      return el;
    };
    fill("sonKey", NOTE_NAMES.map((n, i) => [i, n]));
    fill("sonScale", Object.entries(SCALES).map(([k, s]) => [k, s.label]));
    fill("sonInst", Object.entries(INSTRUMENTS).map(([k, s]) => [k, s.label]));
    fill("sonOct", [2, 3, 4, 5].map((n) => [n, n]));
    fill("sonCross", [10, 20, 30, 60].map((n) => [n, n]));
    const on = (id, fn) => document.getElementById(id).addEventListener("change", (e) => {
      fn(e.target.value); e.target.blur(); syncUI(); scheduleHash();
    });
    on("sonKey", (v) => { cfg.key = +v; });
    on("sonScale", (v) => { cfg.scale = v; });
    on("sonInst", (v) => { cfg.inst = v; });
    on("sonOct", (v) => { cfg.octaves = +v; });
    on("sonCross", (v) => { cfg.cross = +v; });
    document.getElementById("sonVol").addEventListener("input", (e) => { cfg.vol = +e.target.value; });
    on("sonVol", (v) => { cfg.vol = +v; });
    document.getElementById("sonify").addEventListener("click", toggle);
    // a button that keeps the focus would be pressed again by the space bar
    for (const b of document.querySelectorAll("#bar button")) b.addEventListener("click", () => b.blur());
    document.addEventListener("visibilitychange", () => { if (document.hidden) stop(); });
    syncUI();
  }

  function toHash(p) {
    p.set("skey", cfg.key); p.set("sscale", cfg.scale); p.set("sinst", cfg.inst);
    p.set("soct", cfg.octaves); p.set("scross", cfg.cross); p.set("svol", cfg.vol);
  }
  function fromHash(p) {
    const k = parseInt(p.get("skey")); if (k >= 0 && k < 12) cfg.key = k;
    if (Object.hasOwn(SCALES, p.get("sscale") || "")) cfg.scale = p.get("sscale");
    if (Object.hasOwn(INSTRUMENTS, p.get("sinst") || "")) cfg.inst = p.get("sinst");
    const o = parseInt(p.get("soct")); if ([2, 3, 4, 5].includes(o)) cfg.octaves = o;
    const c = parseInt(p.get("scross")); if ([10, 20, 30, 60].includes(c)) cfg.cross = c;
    const v = parseFloat(p.get("svol")); if (v >= 0 && v <= 1.5) cfg.vol = v;
  }

  // ------------------------------------------------------------ self-test (?sonifytest=1)
  async function test() {
    const checks = [], ok = (name, pass, detail = "") => checks.push({ name, pass: !!pass, detail: String(detail) });
    const same = (a, b) => a.length === b.length && a.every((x, i) => x === b[i]);

    // 1. notes of a scale, written out by hand
    ok("C major pentatonic, 4 octaves", same(voicePitches({ key: 0, scale: "pentatonic", octaves: 4 }, 100),
      [45, 48, 50, 52, 55, 57, 60, 62, 64, 67, 69, 72, 74, 76, 79, 81, 84, 86, 88, 91, 93]));
    ok("D natural minor, 2 octaves", same(voicePitches({ key: 2, scale: "minor", octaves: 2 }, 100),
      [57, 58, 60, 62, 64, 65, 67, 69, 70, 72, 74, 76, 77, 79, 81]));
    ok("chromatic, 3 octaves: 37 notes", voicePitches({ key: 5, scale: "chromatic", octaves: 3 }, 100).length === 37);
    ok("no rounding: one band per row up to 64",
       voicePitches({ scale: "none", octaves: 4 }, 30).length === 30 && voicePitches({ scale: "none", octaves: 4 }, 500).length === 64);

    // 2. rows to notes: high at the top, in order, and a neuron taller than a cell is one note
    for (const [name, L] of Object.entries(LAY)) {
      const c = { key: 0, scale: "pentatonic", octaves: 4 }, P = voicePitches(c, 300);
      const m = rowVoices(P, c, 300, 0, L.Y, L);
      let mono = true;
      for (let k = 1; k < 300; k++) if (m.voice[k] > m.voice[k - 1]) mono = false;
      ok(`${name}: top row highest note, bottom row lowest, never rising downward`,
         mono && m.voice[0] === P.length - 1 && m.voice[299] === 0, `${m.voice[0]}, ${m.voice[299]}`);
    }
    {
      const c = { key: 0, scale: "chromatic", octaves: 4 }, P = voicePitches(c, 300), L = LAY.rows;
      const m = rowVoices(P, c, 300, 100.25, 103.25, L);      // 3 neurons' worth on 300 rows
      const notes = new Map();
      for (let k = 0; k < 300; k++) {
        const p = Math.floor(100.25 + (k + 0.5) * 0.01);
        if (!notes.has(p)) notes.set(p, new Set());
        notes.get(p).add(m.voice[k]);
      }
      ok("zoomed in: every cell of a neuron sounds the same note", [...notes.values()].every((s) => s.size === 1),
         [...notes].map(([p, s]) => `${p}: ${[...s]}`).join("; "));
    }

    // 2b. where a cell row sounds, by brute force: of the units of the segment under the row's
    //     middle whose band holds it, the one centered nearest (plotted by depth); the row of
    //     the position under it (other layouts)
    const soundYRef = (k, H, n0, n1, L) => {
      const dn = (n1 - n0) / H, y = n0 + (k + 0.5) * dn;
      if (L.frac) return dn < 1 ? Math.floor(y) + 0.5 : y;
      if (dn >= L.unit) return y;
      const seg = L.segs.find((g, i) => y < g.y1 || i === L.segs.length - 1);
      let best = null;
      for (let i = 0; i < N; i++) {
        if (i < seg.r0 || i >= seg.r1 || !(L.lo[i] <= y && y < L.lo[i] + L.unit)) continue;
        const c = L.lo[i] + L.unit / 2;
        if (best === null || Math.abs(c - y) < Math.abs(best - y)) best = c;
      }
      return best === null ? y : best;
    };
    {
      let nBad = 0, nUnit = 0, n = 0, seed = 7;
      const rnd = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296);
      for (let q = 0; q < 300; q++) {
        const L = LAY.depth, span = 30 + rnd() * 400, n0 = rnd() * (L.Y - span), H = 40 + Math.floor(rnd() * 200);
        for (let k = 0; k < H; k += 7) {
          const a = soundY(k, H, n0, n0 + span, L), b = soundYRef(k, H, n0, n0 + span, L);
          n++; if (a !== b) nBad++;
          if (b !== n0 + (k + 0.5) * span / H) nUnit++;
        }
      }
      ok("plotted by depth, zoomed in: a row sounds the middle of the nearest unit under it (brute force)",
         nBad === 0 && nUnit > 100, `${nBad} of ${n} rows differ; ${nUnit} rows on a unit`);
    }

    // 2c. the transport on a fake clock: every column is scheduled exactly once, before it
    //     sounds, at its own time, from inside the view; the line is over the column heard
    const sim = (name, t0, t1, W, tHead, cross, stalls) => {
      const An = makeAnchor(t0, t1, W, tHead, 0, cross, TE * FS), latency = 0.04, seen = new Map(), bad = [];
      const first = Math.floor(An.PA);
      let t = 0, steps = 0, scrolled = 0, lastCol = -1, lastK = 0, maxOff = 0, done = false;
      for (; steps < 4e6; steps++) {
        t += stalls && steps % 97 === 50 ? 0.2 : 1 / 60;     // now and then a stall shorter than the lookahead
        const p = position(An, t - latency);
        if (p.done) { done = true; break; }
        if (p.k < lastK) bad.push(`view stepped back at ${t}`);
        if (p.k > 0 && p.k < An.kmax) { scrolled++; maxOff = Math.max(maxOff, Math.abs(p.col - p.frac - An.W / 2)); }
        if (p.col < 0 || p.col >= An.W) bad.push(`line outside the view at ${t}`);
        lastK = p.k; lastCol = p.col; An.k = p.k;
        const [g0, g1] = due(An, t);
        for (let g = g0; g < g1; g++) {
          if (seen.has(g)) bad.push(`column ${g} twice`);
          if (g < An.k || g >= An.k + An.W) bad.push(`column ${g} outside the view`);
          if (g !== first && columnTime(An, g) < t - 1e-9) bad.push(`column ${g} late`);
          seen.set(g, t);
        }
        if (g1 > g0) An.gNext = g1;
        const gh = Math.floor(p.P);                        // the column heard now, under the line
        if (!seen.has(gh)) bad.push(`column ${gh} heard before it was scheduled`);
        if (Math.floor(p.col) + p.k !== gh) bad.push(`line not over column ${gh}`);
      }
      const last = An.kmax + An.W - 1;
      let missing = 0;
      for (let g = first; g <= last; g++) if (!seen.has(g)) missing++;
      const wantScroll = An.kmax > 1;
      ok(`transport: ${name}`, done && !bad.length && !missing && seen.size === last - first + 1 &&
         (wantScroll ? scrolled > 0 && maxOff < 1e-6 : scrolled === 0) && lastCol > An.W - 1 - 0.25 * An.rate,
         `${seen.size} columns, ${missing} missing, ${scrolled} scrolling frames, line at column ${lastCol.toFixed(2)} of ${An.W} at the end` +
         (bad.length ? "; " + bad.slice(0, 3).join("; ") : ""));
    };
    sim("60 s window in mid-session, to the end", TE - 400, TE - 340, 738, TE - 370, 30, false);
    sim("the same with stalled frames", TE - 400, TE - 340, 738, TE - 370, 30, true);
    sim("window at the session start: the line starts at the left edge", TS, TS + 2, 300, TS, 10, false);
    sim("whole session in view: the line crosses a standing view", TS, TE, 737, TS, 30, true);
    sim("window at the session end: the line runs to the right edge", TE - 3, TE, 500, TE - 1.5, 30, false);
    sim("shortest window (one sample per column)", TE - 0.5, TE - 0.5 + 500.5 / FS, 500, TE - 0.5 + 250.25 / FS, 60, false);

    // 3. the synthesis, rendered offline: known levels in, measured tones out
    {
      const sr = 48000, oc = new OfflineAudioContext(1, 3 * sr, sr), g = buildGraph(oc, false);
      setVoices(g, [57, 69, 81], INSTRUMENTS.sine);
      g.master.gain.value = 1;
      const col = 0.05, lev = new Float64Array(3);
      for (let c = 0; c < 60; c++) {
        const t = c * col + 1e-9;
        lev[0] = t >= 0.5 && t < 0.55 ? 1 : 0;        // one column: a struck note
        lev[1] = t >= 1.0 && t < 2.0 ? 0.5 : 0;       // 20 columns: a held note
        lev[2] = t >= 1.5 && t < 1.6 ? 0.25 : 0;
        scheduleColumn(g, lev, c * col);
      }
      const x = (await oc.startRendering()).getChannelData(0);
      const amp = (hz, t) => {       // amplitude of the tone at hz around time t (50 ms Hann window)
        const n = Math.round(0.05 * sr), i0 = Math.round(t * sr) - (n >> 1);
        let re = 0, im = 0, sw = 0;
        for (let i = 0; i < n; i++) {
          const w = 0.5 - 0.5 * Math.cos(2 * Math.PI * (i + 0.5) / n), ph = 2 * Math.PI * hz * (i0 + i) / sr;
          re += w * x[i0 + i] * Math.cos(ph); im += w * x[i0 + i] * Math.sin(ph); sw += w;
        }
        return 2 * Math.hypot(re, im) / sw;
      };
      const [v0, v1, v2] = g.voices, { attack, release } = g;
      const held = amp(v1.hz, 1.9) / (0.5 * v1.tilt);
      ok("held note: amplitude = level × tilt", Math.abs(held - 1) < 0.02, held.toFixed(4));
      ok("held note: silent before its first column", amp(v1.hz, 0.95) < 1e-3, amp(v1.hz, 0.95).toExponential(1));
      const rise = amp(v1.hz, 1.0 + 0.025 + 5 * attack) / (0.5 * v1.tilt);
      ok("held note: reached within the attack", rise > 0.9, rise.toFixed(3));
      const fall = amp(v1.hz, 2.0 + release) / (0.5 * v1.tilt);
      ok("held note: falls to 1/e one release time after its last column", Math.abs(fall - Math.exp(-1)) < 0.04, fall.toFixed(3));
      const struck = amp(v0.hz, 0.55) / v0.tilt, ring = amp(v0.hz, 0.55 + release) / v0.tilt;
      // around the end of its column the window holds the last 25 ms of the rise and the first 25 ms of the fall
      ok("struck note: sounds during its column", struck > 0.8, struck.toFixed(3));
      ok("struck note: rings out with the release time", Math.abs(ring / Math.exp(-1) - 1) < 0.15, ring.toFixed(3));
      ok("struck note: silent before", amp(v0.hz, 0.44) < 1e-3, amp(v0.hz, 0.44).toExponential(1));
      const third = amp(v2.hz, 1.575) / (0.25 * v2.tilt);
      ok("third note: level 0.25 during its two columns", third > 0.9 && third < 1.02, third.toFixed(3));
      ok("notes do not leak into each other", amp(v0.hz, 1.9) < 2e-3 && amp(v2.hz, 1.2) < 2e-3,
         `${amp(v0.hz, 1.9).toExponential(1)}, ${amp(v2.hz, 1.2).toExponential(1)}`);
      ok("frequencies: A3, A4, A5", same(g.voices.map((v) => Math.round(v.hz)), [220, 440, 880]));
    }

    // 4. from the picture to the levels: columns of real views, against the float64 readout
    //    (cellCPU, itself checked against reference.py by ?test=1) pooled with separate code
    const views = [
      { name: "whole session, rows", L: LAY.rows, t0: TS, t1: TE, n0: 0, n1: N, W: 240, H: 150, norm: false, spikes: "all" },
      { name: "20 s, 300 rows, normalized", L: LAY.rows, t0: TS + 600, t1: TS + 620, n0: 900.3, n1: 1200.8, W: 200, H: 180, norm: true, spikes: "all" },
      { name: "1 s, 12 neurons (one note per neuron)", L: LAY.rows, t0: TS + 900, t1: TS + 901, n0: 2000.4, n1: 2012.4, W: 300, H: 240, norm: false, spikes: "all" },
      { name: "5 s, plotted by depth, 600 µm", L: LAY.depth, t0: TS + 300, t1: TS + 305, n0: 1000, n1: 1600, W: 200, H: 200, norm: false, spikes: "all" },
      { name: "0.5 s, plotted by depth, 150 µm (one note per neuron)", L: LAY.depth, t0: TS + 300, t1: TS + 300.5, n0: 1300, n1: 1450, W: 200, H: 200, norm: false, spikes: "all" },
    ];
    if (LAY.rastermap) views.push({ name: "60 s, Rastermap order, burst + tonic", L: LAY.rastermap, t0: TS + 1200, t1: TS + 1260, n0: 0, n1: N, W: 200, H: 160, norm: false, spikes: "overlay" });
    const s = { vmin: 0.3, vmax: 50 }, lmin = Math.log(s.vmin), lmax = Math.log(s.vmax);
    for (const tv of views) {
      const v = { T0: tv.t0 * FS, dT: (tv.t1 - tv.t0) * FS / tv.W, n0: tv.n0, n1: tv.n1 };
      const c = { key: 0, scale: "pentatonic", octaves: 4 }, P = voicePitches(c, tv.H);
      const map = rowVoices(P, c, tv.H, tv.n0, tv.n1, tv.L);
      computeValues(v, tv.W, tv.H, tv.norm, 0.1, COUNT_MODE, tv.L, tv.spikes);
      // columns to compare: the middle one, and the first, the last and the loudest with any sound
      const all = readColumns(0, tv.W, tv.H), tmp = new Float64Array(P.length), tot = new Float64Array(tv.W);
      for (let j = 0; j < tv.W; j++)
        tot[j] = columnLevels(all, tv.W, tv.H, j, map, lmin, lmax, tv.spikes === "overlay", tmp).reduce((a, b) => a + b, 0);
      const cols = new Set([Math.floor(tv.W / 2), tot.findIndex((x) => x > 0), tot.findLastIndex((x) => x > 0), tot.indexOf(Math.max(...tot))]);
      cols.delete(-1);
      let worst = 0, loud = 0;
      for (const j of cols) {
        const j0 = Math.min(Math.max(j - 1, 0), tv.W - 3);   // read through the path playing uses: a few columns, not from 0
        const got = columnLevels(readColumns(j0, 3, tv.H), 3, tv.H, j - j0, map, lmin, lmax, tv.spikes === "overlay", new Float64Array(P.length));
        const sum = new Float64Array(P.length), cnt = new Float64Array(P.length);
        for (let k = 0; k < tv.H; k++) {
          const cc = cellCPU(v, tv.W, tv.H, j, k, tv.norm, 0.1, tv.L, tv.spikes);
          let t = cc.value > 0 ? Math.min(1, Math.max(0, Math.log(cc.value / s.vmin) / Math.log(s.vmax / s.vmin))) : 0;
          if (tv.spikes === "overlay" && cc.value2 > 0)
            t = Math.max(t, Math.min(1, Math.max(0, Math.log(cc.value2 / s.vmin) / Math.log(s.vmax / s.vmin))));
          const y = soundYRef(k, tv.H, tv.n0, tv.n1, tv.L);
          const p = 93 - Math.min(1, Math.max(0, (y - tv.n0) / (tv.n1 - tv.n0))) * 48;
          let best = 0;
          for (let q = 1; q < P.length; q++) if (Math.abs(P[q] - p) < Math.abs(P[best] - p)) best = q;
          cnt[best]++;
          if (t > 0) sum[best] += Math.pow(10, -1.8 * (1 - t)) ** 2;
        }
        for (let q = 0; q < P.length; q++) {
          const ref = cnt[q] ? Math.sqrt(sum[q] / cnt[q]) : 0;
          worst = Math.max(worst, Math.abs(got[q] - ref)); loud = Math.max(loud, ref);
        }
      }
      ok(`levels from the picture: ${tv.name}`, worst < 2e-3 && loud > 0.01, `max difference ${worst.toExponential(1)}, loudest ${loud.toFixed(3)}`);
    }
    const out = { passed: checks.filter((c) => c.pass).length, total: checks.length, checks };
    window.__sonifyTest = out;
    console.log(JSON.stringify(out, null, 1));
    return out;
  }

  return { setup, toggle, start, stop, afterRender, grid, drawHead, status, toHash, fromHash, test,
           get playing() { return playing; }, get shiftCss() { return shiftCss; },
           get headX() { return playing && head ? head.x : null; },    // x of the playhead (CSS px), or null
           get debug() { return { A, head, cfg, voices: G ? G.voices.map((v) => ({ midi: v.midi, last: v.last })) : [], ctx: G && G.ctx, master: G && G.master }; } };
})();
