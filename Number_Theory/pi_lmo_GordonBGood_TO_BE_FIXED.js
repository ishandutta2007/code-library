"use strict";

const MAXVALUE = 9007199254740991; // 2**53 - 1

// creates a function returning a lazily memoized value from a thunk...
function lazy(thunk) {
  let value = undefined;
  return function() {
if (value === undefined) { value = thunk(); thunk = null; }
return value;
  }
}

// a page-segmented odds-only bit-packed Sieve of Eratosthenes;

const PGSZBITS = 262144; // about CPU l1 cache size in bits (power of two)

const CLUT = function () { // fast "pop count" Counting Look Up Table...
  const arr = new Uint8Array(65536);
  for (let i = 0; i < 65536; ++i) {
let nmbts = 0 | 0; let v = i;
while (v > 0) { ++nmbts; v &= (v - 1) | 0; }
arr[i] = nmbts | 0; }
  return arr;
}();

function countPageFromTo(bitstrt, bitlmt, sb) {
  const fst = bitstrt >> 5; const lst = bitlmt >> 5;
  const pg = new Uint32Array(sb.buffer);
  let v0 = (pg[fst] | ((0xFFFFFFFF >>> 0) << (bitstrt & 31))) >>> 0;
  let cnt = ((lst - fst) << 5) + CLUT[v0 & 0xFFFF]; cnt += CLUT[v0 >>> 16];
  for (let i = fst | 0; i < lst; ++i) {
let v = pg[i] >>> 0;
cnt -= CLUT[v & 0xFFFF]; cnt -= CLUT[v >>> 16];
  }
  let v1 = (pg[lst] | ((0xFFFFFFFE >>> 0) << (bitlmt & 31))) >>> 0;
  cnt -= CLUT[v1 & 0xFFFF]; cnt -= CLUT[v1 >>> 16]; return cnt | 0;
}

function partialSievePage(lwi, bp, sb) {
  const btsz = sb.length << 3;
  let s = Math.trunc((bp * bp - 3) / 2); // compute the start index...
  if (s >= lwi) s -= lwi; // adjust start index based on page lower limit...   
  else { // for the case where this isn't the first prime squared instance
let r = ((lwi - s) % bp) >>> 0;
s = (r != (0 >>> 0) ? bp - r : 0) >>> 0; }
  if (bp <= 32) {
for (let slmt = Math.min(btsz, s + (bp << 3)); s < slmt; s += bp) {
  const shft = s & 7; const msk = ((1 >>> 0) << shft) >>> 0;
  for (let c = s >> 3, clmt = sb.length; c < clmt | 0; c += bp)
    sb[c] |= msk; } }
  else
for (let slmt = sb.length << 3; s < slmt; s += bp)
  sb[s >> 3] |= ((1 >>> 0) << (s & 7)) >>> 0;
}

function partialSieveCountPage(lwi, bp, cntarr, sb) {
  const btsz = sb.length << 3; let cullcnt = 0;
  let s = Math.trunc((bp * bp - 3) / 2); // compute the start index...
  if (s >= lwi) // adjust start index based on page lower limit...
s -= lwi;
  else { // for the case where this isn't the first prime squared instance
let r = ((lwi - s) % bp) >>> 0;
s = (r != (0 >>> 0) ? bp - r : 0) >>> 0; }
  if (bp <= 32) {
for (let slmt = Math.min(btsz, s + (bp << 3)); s < slmt; s += bp) {
  const shft = s & 7; const msk = ((1 >>> 0) << shft) >>> 0;
  for (let c = s >>> 3, clmt = sb.length; c < clmt | 0; c += bp) {
    const isbit = ((sb[c] >>> shft) ^ 1) & 1;
    cntarr[c >> 6] -= isbit; cullcnt += isbit; sb[c] |= msk; }
}
  }
  else
for (let slmt = sb.length << 3; s < slmt; s += bp) {
  const sba = s >>> 3; const shft = s & 7;
  const isbit = ((sb[sba] >>> shft) ^ 1) & 1;
  cntarr[s >> 9] -= isbit; cullcnt += isbit;
  sb[sba] |= ((1 >>> 0) << shft) >>> 0; }
  return cullcnt;
}

// pre-culled pattern of small wheel primes...
const WHLPRMS = [ 2, 3, 5, 7, 11, 13, 17 ];
const WHLPTRNLEN = WHLPRMS.reduce((s, v) => s * v, 1) >>> 1; // odds only!
const WHLPTRN = function() { // larger than WHLPTRN by one buffer for overflow
  const len = (WHLPTRNLEN + (PGSZBITS >>> 3) + 3) & (-4); // up 2 even 32 bits!
  const arr = new Uint8Array(len);
  for (let bp of WHLPRMS.slice(1)) partialSievePage(0, bp, arr);
  arr[0] |= ~(-2 << ((WHLPRMS[WHLPRMS.length - 1] - 3) >> 1)) >>> 0; return arr;
}();

function fillPage(lwi, sb) {
  const mod = (lwi / 8) % WHLPTRNLEN;
  sb.set(new Uint8Array(WHLPTRN.buffer, mod, sb.length));
}

function cullPage(lwi, bpras, sb) {
  const btsz = sb.length << 3; let bp = 3;
  const nxti = lwi + btsz; // just beyond the current page 
  for (let bpra of bpras()) {
for (let bpri = 0; bpri < bpra.length; ++bpri) {
  const bpr = bpra[bpri]; bp += bpr + bpr;
  let s = (bp * bp - 3) / 2; // compute start index of prime squared
  if (s >= nxti) return; // enough bp's
  partialSievePage(lwi, bp, sb);
}
  }
}

function soePages(bitsz, bpras) {
  const buf =  new Uint8Array(bitsz >> 3); let lowi = 0;
  const gen = bpras === undefined ? makeBasePrimeRepArrs() : bpras;
  return function*() {
while (true) {
  fillPage(lowi, buf); cullPage(lowi, gen, buf);
  yield { lwi: lowi, sb: buf }; lowi += bitsz; }
  };
}

function makeBasePrimeRepArrs() {
  const buf = new Uint8Array(128); let gen = undefined; // avoid data race!
  fillPage(0, buf);
  for (let i = 8, bp = 19, sqr = bp * bp; sqr < 2048+3;
                                      ++i, bp += 2, sqr = bp * bp)
if (((buf[i >> 3] >>> 0) & ((1 << (i & 7)) >>> 0)) === 0 >>> 0)
  for (let c = (sqr - 3) >> 1; c < 1024; c += bp)
    buf[c >> 3] |= (1 << (c & 7)) >>> 0; // init zeroth buf
  function sb2bprs(sb) {
const btsz = sb.length << 3; let oi = 0;
const arr = new Uint8Array(countPageFromTo(0, btsz - 1, sb));
for (let i = 0, j = 0; i < btsz; ++i)
  if (((sb[i >> 3] >>> 0) & ((1 << (i & 7)) >>> 0)) === 0 >>> 0) {
    arr[j++] = (i - oi) >>> 0; oi = i; }
return { bpra: arr, lastgap: (btsz - oi) | 0 };
  }
  let { bpra, lastgap } = sb2bprs(buf);
  function next() {
const nxtpg = sb2bprs(gen.next().value.sb);
nxtpg.bpra[0] += lastgap; lastgap = nxtpg.lastgap;
return { head: nxtpg.bpra, tail: lazy(next) };
  }
  const lazylist = { head: bpra, tail: lazy(function() {
if (gen === undefined) {
  gen = soePages(1024)(); gen.next() } // past first page
return next();
  }) };
  return function*() { // return a generator of rep pages...
let ll = lazylist; while (true) {  yield ll.head; ll = ll.tail(); }
  };
}

function *revPrimesFrom(top, bpras) {
  const topndx = (top - 3) >>> 1;
  const buf = new Uint8Array(PGSZBITS >>> 3);
  let lwi = (((topndx / PGSZBITS) >>> 0) * PGSZBITS) >>> 0;
  let si = (topndx - lwi) >>> 0;
  for (; lwi >= 0; lwi -= PGSZBITS) { // usually external limit!
const base = 3 + lwi + lwi;
fillPage(lwi, buf); cullPage(lwi, bpras, buf);
for (; si >= 0 >>> 0; --si)
  if (((buf[si >> 3] >>> 0) & ((1 << (si & 7)) >>> 0)) === (0 >>> 0))
    yield base + si + si;
si = PGSZBITS - 1;
  }
};

const TinyPrimes = [ 2, 3, 5, 7, 11, 13, 17, 19 ]; // degree eight
const TinyProduct = TinyPrimes.reduce((s, v) => s * v) >>> 0;
const TinyHalfProduct = TinyProduct >>> 1;
const TinyTotient = TinyPrimes.reduce((s, v) => s * (v - 1), 1) >>> 0;
const TinyLength = (TinyProduct + 8) >>> 2; // include zero and half point!
const TinyTotients = function() {
  const arr = new Uint32Array(TinyLength | 0);
  arr[TinyLength - 1] = 1; // mark mid point value as not prime - never is
  let spn = 3 * 5 * 7; arr[0] = 1; // mark zeroth value as not prime!
  for (let bp of [ 3, 5, 7 ]) // cull small base prime values...
for (let c = (bp + 1) >>> 1; c <= spn; c += bp) arr[c] |= 1;
  for (let bp of [ 11, 13, 17, 19 ]) {
for (let i = 1 + spn; i < TinyLength; i += spn) {
  const rng = i + spn > TinyLength ? spn >> 1 : spn;
  arr.set(new  Uint32Array(arr.buffer, 4, rng), i); }
spn *= bp;
for (let c = (bp + 1) >>> 1; // eliminate prime in pattern!
       c < (spn > TinyLength ? TinyLength : spn + 1); c += bp)
  arr[c] |= 1;
  }
  arr.reduce((s, v, i) => { // accumulate sums...
const ns = s + (v ^ 1); arr[i] = ns; return ns; }, 0);
  return arr;
}();  

function tinyPhi(m) {
  const d = Math.trunc(m / TinyProduct);
  const ti = (m - d * TinyProduct + 1) >>> 1;
  const t = ti < TinyLength
          ? TinyTotients[ti]
          : TinyTotient - TinyTotients[TinyHalfProduct - ti];
  return d * TinyTotient + t;
}

function *countPrimesTo(limit) {
  if (limit <= WHLPRMS[WHLPRMS.length - 1]) {
let cnt = 0; for (let p of WHLPRMS) { if (p > limit) break; else ++cnt; }
return cnt; }

  const bpras = makeBasePrimeRepArrs();
  if (limit < 1024**2 + 3) { // for limit < about a million, just sieve...
let p = 3; let cnt = WHLPRMS.length;
for (let bpra of bpras())
  for (let bpr of bpra) { // just count base prime values to limit
    p += bpr + bpr; if (p > limit) return cnt; ++cnt; }
  }

  if (limit <= 32 * 2 * PGSZBITS + 3) { // count sieve to about 32 million...
const lmti = (limit - 3) / 2;
let cnt = WHLPRMS.length; // just use page counting to limit as per usual...
for (let pg of soePages(PGSZBITS, bpras)()) {
  const nxti = pg.lwi + (pg.sb.length << 3);
  if (nxti > lmti) { cnt += countPageFromTo(0, lmti - pg.lwi, pg.sb); break; }
  cnt += countPageFromTo(0, PGSZBITS - 1, pg.sb);
}
return cnt;
  }

  // Actual LMO prime counting code starts here...
  const sqrt = Math.trunc(Math.sqrt(limit)) >>> 0;
  const cbrt = Math.trunc(Math.cbrt(limit)) >>> 0;
  const sqrtcbrt = Math.trunc(Math.sqrt(cbrt)) >>> 0;
  const top = Math.trunc(limit / cbrt); // sized for maximum required!
  const bsprms = function() {
let bp = 3; let cnt = WHLPRMS.length + 1; for (let bpra of bpras())
  for (let bpr of bpra) {
    bp += bpr + bpr; if (bp > cbrt) return new Uint32Array(cnt); ++cnt; }
  }();
  bsprms.set(WHLPRMS, 1); // index zero not used == 0!
  const pisqrtcbrt = function() {
let cnt = WHLPRMS.length; let i = cnt + 1; let bp = 3;
stop: for (let bpra of bpras())
  for (let bpr of bpra) {
    bp += bpr + bpr; if (bp > cbrt) break stop;
    if (bp <= sqrtcbrt) ++cnt; bsprms[i++] = bp >>> 0; }
return cnt;
  }();
  const pis = function() { // starts with index 0!
const arr = new Uint32Array(cbrt + 2); let j = 0;
for (let i = 1; i < bsprms.length; ++i)
  for (; j < bsprms[i]; ) arr[j++] = (i - 1) >>> 0;
for (; j < arr.length; ) arr[j++] = (bsprms.length - 1) >>> 0;
return arr;
  }();
  const phis = function() { // index below TinyPhi degree never used...
const arr = (new Array(bsprms.length)).fill(1);
arr[0] = 0; arr[1] = 3; arr[2] = 2; // unused
for (let i = WHLPRMS.length + 2; i < arr.length; ++i) {
  arr[i] -= i - WHLPRMS.length - 1; } // account for non phi primes!
return arr;
  }();
  // indexed by `m`, contains lpf and Moebius value bit; one is negative...
  const specialroots = new Uint16Array(cbrt + 1); // filled in with S1 below...
  const S1 = function() { // it is very easy to calculate S1 recursively...
let s1acc = tinyPhi(limit);
function level(lpfni, lmtlpfni, mfv, m) {
  while (lpfni < lmtlpfni) {
    const pn = bsprms[lpfni]; const nm = m * pn;
    if (nm > cbrt) { // don't split, found S2 root leaf...
      specialroots[m] = (lmtlpfni << 1) | (mfv < 0 ? 1 : 0); return; }
    else { // recurse for S1; never more than 11 levels deep...
      s1acc += mfv * tinyPhi(Math.trunc(limit / nm)); // down level...
      level(9, lpfni, -mfv, nm); // Moebius sign change on down level!
      ++lpfni; } // up prime value, same level!
  }
}
level(9, bsprms.length, -1, 1); return s1acc;
  }();

  // at last, calculate the more complex parts of the final answer:
  function *complex() {
let s2acc = 0; let p2acc = 0; let p2cnt = 0; // for "P2" calculation
const buf = new Uint8Array(PGSZBITS >>> 3); let ttlcnt = 0;
const cnts = new Uint8Array(PGSZBITS >>> 9);
const cntaccs = new Uint32Array(cnts.length);
const revgen = revPrimesFrom(sqrt, bpras);
let rp = revgen.next().value; let p2v = Math.trunc(limit / rp);
const lwilmt = Math.trunc((top - 3) / 2);

const updtmsk = ((PGSZBITS << 3) - 1) >>> 0;
let nxttm = Date.now() + 1000;
for (let lwi = 0; lwi <= lwilmt; lwi += PGSZBITS) { // for all pages
  if ((lwi & updtmsk) == updtmsk && Date.now() >= nxttm) {
    nxttm = Date.now() + 1000; yield lwi / lwilmt * 100.0 };
  let pgcnt = 0; const low = 3 + lwi + lwi;
  const high = Math.min(low + (PGSZBITS << 1) - 2, top);
  let cntstrti = 0 >>> 0;
  function countTo(stop) {
    const cntwrd = stop >>> 9; const bsndx = stop & (-512);
    const xtr = countPageFromTo(bsndx, stop, buf);
    while (cntstrti < cntwrd) {
      const ncnt = cntaccs[cntstrti] + cnts[cntstrti];
      cntaccs[++cntstrti] = ncnt; }
    return cntaccs[cntwrd] + xtr;
  }
  const bpilmt = pis[Math.trunc(Math.sqrt(high)) >>> 0] >>> 0;
  const maxbpi = pis[Math.min(cbrt, Math.trunc(Math.sqrt(limit/low)))]>>>0;
  const tminbpi = pis[Math.min(Math.trunc(top / (high + 2)),
                               bsprms[maxbpi]) >>> 0];
  const minbpi = Math.max(TinyPrimes.length, tminbpi) + 1;
  fillPage(lwi, buf); let bpi = (WHLPRMS.length + 1) >>> 0;
 
  if (minbpi <= maxbpi) { // jump to doing P2 if not range

    // for bpi < minbpi there are no special leaves...
    for (; bpi < minbpi; ++bpi) { // eliminate all Tiny Phi primes...
      const bp = bsprms[bpi]; const i = (bp - 3) >>> 1; // cull base primes!
      phis[bpi] += countPageFromTo(0, PGSZBITS - 1, buf);
      partialSievePage(lwi, bp, buf); }
    for (let i = 0; i < cnts.length; ++i) { // init cnts arr...
      const s = i << 9; const c = countPageFromTo(s, s + 511, buf);
      cnts[i] = c; pgcnt += c; }

    // for all base prime values up to limit**(1/6) in the page,
    // add all special leaves composed of this base prime value and
    // any number of other higher base primes, all different,
    // that qualify as special leaves...
    let brkchkr = false;
    for (; bpi <= Math.min(pisqrtcbrt, maxbpi) >>> 0; ++bpi) {
      const bp = bsprms[bpi];
      const minm = Math.max(Math.trunc(limit / (bp * (high + 2))),
                            Math.trunc(cbrt / bp)) >>> 0;
      const maxm = Math.min(Math.trunc(limit / (bp * low)), cbrt) >>> 0;
      if (bp >= maxm) { brkchkr = true; break; }
      for (let m = maxm; m > minm; --m) {
        const rt = specialroots[m];
        if (rt != 0 && bpi < rt >>> 1) {
          const stop = Math.trunc(limit / (bp * m) - low) >>> 1;
          const mu = ((rt & 1) << 1) - 1; // one bit means negative!
          s2acc -= mu * (phis[bpi] + countTo(stop));
        } }
      phis[bpi] += pgcnt; // update intermediate base prime counters
      pgcnt -= partialSieveCountPage(lwi, bp, cnts, buf);
      cntstrti = 0; cntaccs[0] = 0;
    }
    // for all base prime values > limit**(1/6) in the page,
    // add results of all special levaes compoaed using only two primes...
    if (!brkchkr)
    for (; bpi <= maxbpi; ++bpi) {
      const bp = bsprms[bpi];
      let l = pis[Math.min(Math.trunc(limit / (bp * low)), cbrt)>>>0]>>>0;
      if (bp >= bsprms[l]) break;
      const piminm = pis[Math.max(Math.trunc(limit / (bp * (high + 2))),
                                  bp) >>> 0] >>> 0;
      for (; l > piminm; --l) {          
        const stop = Math.trunc(limit / (bp * bsprms[l]) - low) >>> 1;
        s2acc += phis[bpi] + countTo(stop);
      }
      phis[bpi] += pgcnt; // update intermediate base prime counters
      if (bpi <= bpilmt) {
        pgcnt -= partialSieveCountPage(lwi, bp, cnts, buf);
        cntstrti = 0; cntaccs[0] = 0; }
    }
  }

  // complete cull page segment, then count up "P2" terms in range...
  for (; bpi <= bpilmt; ++bpi) partialSievePage(lwi, bsprms[bpi], buf);
  let ndx = 0 >>> 0;
  while (p2v >= low && p2v < high + 2) {
    const nndx = (p2v - low) >>> 1;
    ++p2cnt; ttlcnt += countPageFromTo(ndx, nndx, buf);
    p2acc += ttlcnt;
    ndx = (nndx + 1) >>> 0; p2v = Math.trunc(limit / revgen.next().value);
  }
  if (ndx < PGSZBITS) ttlcnt += countPageFromTo(ndx, PGSZBITS - 1, buf);
}
// adjust for now known delta picbrt to pisqrt!
p2acc -= p2cnt * (p2cnt - 1) / 2 +
           p2cnt * (bsprms.length - 1 - WHLPRMS.length);
const Piy = bsprms.length - 1;
return S1 + s2acc - p2acc + Piy - 1;
  }
  const gen = complex(); let lst = gen.next();
  for (; !lst.done; lst = gen.next()) yield lst.value;
  return lst.value;
}

let cancelled = false;

function doit() {
  const limit =  Math.floor(parseFloat(document.getElementById('limit').value));
  const start = Date.now();
  const pggen = countPrimesTo(limit);
  function pgfnc() {
if (cancelled) {            
  document.getElementById('output').innerText = "Cancelled!!!";
}
else {
  const pgrslt = pggen.next();
  if (!pgrslt.done) {
    // value is a percent done...            
    document.getElementById('output').innerText = "Sieved " + (pgrslt.value.toFixed(3)) + "%";
    setTimeout(pgfnc, 7); return;
  }
  // value is the count result...
  const elpsd = Date.now() - start;
  document.getElementById('output').innerText = "Found " + pgrslt.value 
+ " primes up to " + limit + " in " + elpsd + " milliseconds.";      
}
cancelled = false;    
document.getElementById('go').onclick = strtclick;
document.getElementById('go').value = "Start Sieve...";            
document.getElementById('go').disabled = false;
  }
  pgfnc();
}

const cancelclick = function () {
  cancelled = true;
  document.getElementById('go').disabled = true;
  document.getElementById('go').value = "Cancelled!!!";
  document.getElementById('go').onclick = strtclick;
}

const strtclick = function () {
  const limit =  Math.floor(parseFloat(document.getElementById('limit').value));
  if (!Number.isInteger(limit) || (limit < 0) || (limit > MAXVALUE)) {
    
  document.getElementById('output').innerText = "Top limit must be an integer between 0 and " + MAXVALUE + "!";
return;
  }
  document.getElementById('output').innerText = "Sieved 0%";
  document.getElementById('go').onclick = cancelclick;
  document.getElementById('go').value = "Running, click to cancel...";
  cancelled = false;
  setTimeout(doit, 7);
};

document.getElementById('go').onclick = strtclick;

