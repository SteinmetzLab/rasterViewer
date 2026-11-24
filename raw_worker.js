"use strict";
// Worker for the raw-voltage underlay (raw_underlay.js): decodes chunks of raw data (downloaded
// by the page) and applies the display processing. raw_underlay/ holds the Python reference of
// both steps (make_raw.py decode_tile, make_raw_proto.py process_tile); ?rawtest=1 compares.
//
// A chunk is one file of a Zarr v2 array (channels x samples, int16, C order) stored with a
// delta filter and zlib: inflate, then a running sum over the flattened array (int16 wraps).

const chunks = new Map();          // url -> Promise<Int16Array>, most recently used last
const MAX_CHUNKS = 16;

function getChunk(url, buf, nch, n) {
  let p = chunks.get(url);
  if (p) { chunks.delete(url); chunks.set(url, p); return p; }
  p = (async () => {
    const out = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream("deflate"))).arrayBuffer();
    if (out.byteLength !== nch * n * 2) throw new Error(`${url}: ${out.byteLength} bytes, expected ${nch * n * 2}`);
    const a = new Int16Array(out);
    for (let i = 1; i < a.length; i++) a[i] += a[i - 1];
    return a;
  })();
  chunks.set(url, p);
  p.catch(() => chunks.delete(url));
  while (chunks.size > MAX_CHUNKS) chunks.delete(chunks.keys().next().value);
  return p;
}

// One second-order section, direct form II transposed (as scipy.signal.sosfilt), from rest,
// over v[off + k * step] for k = 0 .. len-1.
function section(v, off, len, step, s) {
  const b0 = s[0], b1 = s[1], b2 = s[2], a1 = s[4], a2 = s[5];
  let z1 = 0, z2 = 0;
  for (let k = 0, i = off; k < len; k++, i += step) {
    const x = v[i], y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    v[i] = y;
  }
}

// k-th smallest of a[0 .. n-1] (quickselect; reorders a so that a[0 .. k-1] <= a[k]).
function select(a, n, k) {
  let lo = 0, hi = n - 1;
  while (lo < hi) {
    const p = a[(lo + hi) >> 1];
    let i = lo, j = hi;
    while (i <= j) {
      while (a[i] < p) i++;
      while (a[j] > p) j--;
      if (i <= j) { const t = a[i]; a[i] = a[j]; a[j] = t; i++; j--; }
    }
    if (k <= j) hi = j; else if (k >= i) lo = i; else break;
  }
  return a[k];
}

// Display processing of one tile. prev / next: the neighboring chunks, or null at the ends.
//   1. the tile plus up to `margin` samples on each side; counts = raw * scale + offset
//   2. per channel: subtract the mean of that stretch; high-pass forward, then backward
//   3. dead channels = 0; per sample, subtract the median over the other channels
//   4. cut the margins off; microvolts
function processTile(prev, cur, next, p) {
  const nch = p.nch, n = p.tile, a = prev ? p.margin : 0, b = next ? p.margin : 0, len = a + n + b;
  const buf = new Float64Array(nch * len);
  const isDead = new Uint8Array(nch);
  for (const c of p.dead) isDead[c] = 1;
  for (let c = 0; c < nch; c++) {
    const o = c * len, src = c * n;
    let sum = 0;
    for (let i = 0; i < a; i++) buf[o + i] = prev[src + n - a + i] * p.scale + p.offset;
    for (let i = 0; i < n; i++) buf[o + a + i] = cur[src + i] * p.scale + p.offset;
    for (let i = 0; i < b; i++) buf[o + a + n + i] = next[src + i] * p.scale + p.offset;
    for (let i = 0; i < len; i++) sum += buf[o + i];
    const mean = sum / len;
    for (let i = 0; i < len; i++) buf[o + i] -= mean;
    for (const s of p.sos) section(buf, o, len, 1, s);
    for (const s of p.sos) section(buf, o + len - 1, len, -1, s);
  }
  const good = [];
  for (let c = 0; c < nch; c++) if (!isDead[c]) good.push(c);
  const ng = good.length, tmp = new Float64Array(ng), out = new Float32Array(nch * n);
  const h = ng >> 1, even = ng % 2 === 0;
  // The medians of BLK samples at a time, so that memory is read along each channel.
  const BLK = 64, blk = new Float64Array(ng * BLK), meds = new Float64Array(BLK);
  for (let i0 = 0; i0 < n; i0 += BLK) {
    const m = Math.min(BLK, n - i0);
    for (let g = 0; g < ng; g++) {
      const src = good[g] * len + a + i0;
      for (let j = 0; j < m; j++) blk[j * ng + g] = buf[src + j];
    }
    for (let j = 0; j < m; j++) {
      const col = blk.subarray(j * ng, (j + 1) * ng);
      tmp.set(col);
      let med = select(tmp, ng, h);
      if (even) {                      // mean of the two middle values, as numpy.median
        let lower = tmp[0];
        for (let g = 1; g < h; g++) if (tmp[g] > lower) lower = tmp[g];
        med = 0.5 * (lower + med);
      }
      meds[j] = med;
    }
    for (let g = 0; g < ng; g++) {
      const c = good[g], src = c * len + a + i0, dst = c * n + i0;
      for (let j = 0; j < m; j++) out[dst + j] = (buf[src + j] - meds[j]) * p.uv;
    }
  }
  return out;
}

onmessage = async (e) => {
  const m = e.data;
  try {
    if (m.type === "tile") {
      const [prev, cur, next] = await Promise.all(m.urls.map((u, i) => u ? getChunk(u, m.bufs[i], m.p.nch, m.p.tile) : null));
      const t = performance.now();
      const data = processTile(prev, cur, next, m.p);
      postMessage({ id: m.id, ok: true, data, ms: performance.now() - t }, [data.buffer]);
    } else if (m.type === "sum") {       // checksum of a decoded chunk (for ?rawtest=1)
      const a = await getChunk(m.url, m.buf, m.nch, m.tile);
      let s = 0;
      for (let i = 0; i < a.length; i++) s += a[i];
      postMessage({ id: m.id, ok: true, sum: s });
    }
  } catch (err) {
    postMessage({ id: m.id, ok: false, error: String(err && err.message || err) });
  }
};
